import os
import shutil
from datetime import datetime
import gzip
import tempfile
from pathlib import Path
from typhon.command_utils import run_command


def decompress_if_gzipped(input_path, log):
    """
    If the file is gzipped, decompress to a temporary file and return its path.
    Otherwise, return the original path.
    The temp file will be deleted by the caller.
    """
    # Convert PosixPath to string if needed
    input_path_str = str(input_path)
    if input_path_str.endswith('.gz'):
        # Use only the file extension as the suffix (e.g., .fastq)
        suffix = os.path.splitext(os.path.splitext(input_path_str)[0])[1]
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        log(f"Decompressing {input_path_str} to temp file {tmp.name}")
        with gzip.open(input_path_str, 'rb') as f_in:
            with open(tmp.name, 'wb') as f_out:
                # Read and write in chunks to handle large files
                while True:
                    chunk = f_in.read(8192)
                    if not chunk:
                        break
                    # Ensure we're writing bytes
                    if isinstance(chunk, str):
                        chunk = chunk.encode('utf-8')
                    f_out.write(chunk)
        return tmp.name
    return input_path_str


def get_genion_bin():
    """
    Resolve the Genion binary path relative to the Typhon package root.
    """
    module_dir = os.path.dirname(os.path.abspath(__file__))
    typhon_root = os.path.abspath(os.path.join(module_dir, '..', '..'))
    bin_path = os.path.join(typhon_root, 'bin', 'genion')
    if not os.path.isfile(bin_path):
        raise FileNotFoundError(f"Genion binary not found at {bin_path}. Please run setup_genion.py or check your installation.")
    return bin_path


def run_genion(
    input_fastq,
    input_sam,
    gtf_for_genion,
    selfalign_paf,
    selfalign_tsv,
    output_dir,
    threads=1,
    genion_bin=None,
    genomic_superdups=None,
    keep_intermediate=False,
    log_path=None,
    min_support=1
):
    """
    Run the Genion pipeline step for Typhon, for a single sample.
    Assumes reference files (GTF, self-align PAF, self-align TSV) are already prepared.
    Handles compressed (.gz) FASTQ and SAM files.
    All output files go directly to output_dir.
    
    Args:
        input_fastq: Path to input FASTQ file
        input_sam: Path to input SAM file (from LongGF/minimap2)
        gtf_for_genion: Path to GTF file for Genion
        selfalign_paf: Path to self-alignment PAF file
        selfalign_tsv: Path to self-alignment TSV file
        output_dir: Output directory for results
        threads: Number of threads (not used by Genion itself)
        genion_bin: Path to custom Genion binary
        genomic_superdups: Path to genomic segmental duplications file
        keep_intermediate: Whether to keep intermediate files
        log_path: Path to log file
        min_support: Minimum supporting reads for fusion calls (default: 1)
    """
    # Set up logging
    if log_path is None:
        log_path = os.path.join(output_dir, 'run_genion.log')
    def log(msg):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(log_path, 'a') as f:
            f.write(f"[{timestamp}] {msg}\n")
        print(f"[{timestamp}] {msg}")

    os.makedirs(output_dir, exist_ok=True)
    temp_files = []
    try:
        # Use consistent sample naming for output files
        sample_name = Path(input_fastq).stem.replace('.fastq', '').replace('.fq', '').replace('.gz', '')

        # Decompress FASTQ and SAM if needed (for tools that require files)
        fastq_for_use = decompress_if_gzipped(input_fastq, log)
        if fastq_for_use != input_fastq:
            temp_files.append(fastq_for_use)
        sam_for_use = decompress_if_gzipped(input_sam, log)
        if sam_for_use != input_sam:
            temp_files.append(sam_for_use)
        # Decompress GTF and self-align reference files if needed
        gtf_for_use = decompress_if_gzipped(gtf_for_genion, log)
        if gtf_for_use != gtf_for_genion:
            temp_files.append(gtf_for_use)
        selfalign_tsv_for_use = decompress_if_gzipped(selfalign_tsv, log)
        if selfalign_tsv_for_use != selfalign_tsv:
            temp_files.append(selfalign_tsv_for_use)

        # 1. Convert SAM to PAF using paftools.js
        log('Converting SAM to PAF with paftools.js...')
        paf_file = os.path.join(output_dir, f'{sample_name}.paf')
        # TODO: Set paftools.js path as needed
        paftools_cmd = f"paftools.js sam2paf {sam_for_use} > {paf_file}"
        run_command(paftools_cmd)
        log(f'PAF file generated: {paf_file}')

        # 2. Prepare Genion arguments
        if genion_bin is None:
            genion_bin = get_genion_bin()
        if not genomic_superdups:
            # Create an empty file if not provided
            genomic_superdups = os.path.join(output_dir, 'genomicSuperDups.txt')
            open(genomic_superdups, 'a').close()
        genion_out = os.path.join(output_dir, f'{sample_name}_genion.tsv')
        genion_cmd = [
            genion_bin,
            '-i', fastq_for_use,
            '--gtf', gtf_for_use,
            '-g', paf_file,
            '-s', selfalign_tsv_for_use,
            '-d', genomic_superdups,
            '-o', genion_out,
            '--min-support', str(min_support)
        ]
        log(f"Running Genion: {' '.join(str(x) for x in genion_cmd)}")
        run_command(genion_cmd, shell=False)
        log(f'Genion output: {genion_out}')

        # 3. Cleanup intermediate files if requested
        if not keep_intermediate:
            for f in [paf_file]:
                try:
                    os.remove(f)
                    log(f'Removed intermediate file: {f}')
                except Exception as e:
                    log(f'WARNING: Could not remove {f}: {e}')
            # Remove decompressed temp files
            for f in temp_files:
                try:
                    os.remove(f)
                    log(f'Removed decompressed temp file: {f}')
                except Exception as e:
                    log(f'WARNING: Could not remove temp file {f}: {e}')
        
        # 4. Handle debug output files - FIXED: Keep .fail files for debug mode
        # Only remove .log files (Genion's standard log, not our pipeline log)
        # The .fail files contain important debug information and should be preserved
        genion_log_file = os.path.join(output_dir, f'{sample_name}_genion.tsv.log')
        if os.path.exists(genion_log_file):
            try:
                os.remove(genion_log_file)
                log(f'Removed Genion log file: {genion_log_file}')
            except Exception as e:
                log(f'WARNING: Could not remove {genion_log_file}: {e}')
        
        # Check if .fail file exists and log its presence (don't remove it!)
        genion_fail_file = os.path.join(output_dir, f'{sample_name}_genion.tsv.fail')
        if os.path.exists(genion_fail_file):
            # Check file size to report if it contains debug data
            fail_size = os.path.getsize(genion_fail_file)
            if fail_size > 0:
                log(f'Debug file preserved: {genion_fail_file} ({fail_size} bytes)')
            else:
                # Remove empty .fail files as per our custom Genion patch
                try:
                    os.remove(genion_fail_file)
                    log(f'Removed empty debug file: {genion_fail_file}')
                except Exception as e:
                    log(f'WARNING: Could not remove empty {genion_fail_file}: {e}')
        
        # 5. Placeholder for post-processing (R scripts)
        log('Post-processing step not yet implemented.')

    except Exception as e:
        log(f'ERROR: {e}')
        raise 