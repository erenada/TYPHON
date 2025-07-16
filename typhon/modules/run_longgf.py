import os
import glob
import gzip
import shutil
import tempfile
from pathlib import Path
from typhon.command_utils import run_command


def decompress_if_gzipped(input_path):
    """
    If the file is gzipped, decompress to a temporary file and return its path.
    Otherwise, return the original path.
    """
    if input_path.endswith('.gz'):
        # Use only the file extension as the suffix (e.g., .fa, .gtf)
        suffix = os.path.splitext(os.path.splitext(input_path)[0])[1]
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        tmp.close()  # Close the file handle so we can write to it
        
        with gzip.open(input_path, 'rb') as f_in:
            with open(tmp.name, 'wb') as f_out:
                f_out.write(f_in.read())
        return tmp.name
    return input_path


def postprocess(output_dir):
    """
    Import and call postprocess function to avoid circular import issues.
    """
    try:
        from .postprocess import postprocess as _postprocess
        _postprocess(output_dir)
    except ImportError:
        # Fallback: call R script directly
        import subprocess
        module_dir = Path(__file__).parent
        r_script_path = module_dir / "LongGF_process_results.R"
        subprocess.run(["Rscript", str(r_script_path), output_dir], check=True)


def run_longgf(fastq_dir, genome, gtf, output_dir, threads=1, keep_intermediate=False, log_path=None):
    """
    Run LongGF pipeline step:
    - Aligns each FASTQ with minimap2
    - Converts SAM to BAM and sorts
    - Runs LongGF
    - Post-processes LongGF logs
    - Further processes results and writes Excel files
    
    Args:
        fastq_dir: Directory containing FASTQ files
        genome: Path to reference genome FASTA file
        gtf: Path to GTF annotation file
        output_dir: Output directory for results
        threads: Number of threads to use (default: 1)
        keep_intermediate: Whether to keep intermediate files (default: False)
        log_path: Path to log file (optional)
    
    Returns:
        dict: Results summary with file paths and statistics
    """
    def log(msg):
        if log_path:
            with open(log_path, 'a') as f:
                f.write(f"{msg}\n")
        print(msg)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Support gzipped fastq
    fastq_files = glob.glob(os.path.join(fastq_dir, '*.fastq')) + glob.glob(os.path.join(fastq_dir, '*.fastq.gz'))
    
    if not fastq_files:
        raise ValueError(f"No FASTQ files found in {fastq_dir}")
    
    log(f"Found {len(fastq_files)} FASTQ files for processing")
    
    # Decompress reference files if needed
    genome_path = decompress_if_gzipped(genome)
    gtf_path = decompress_if_gzipped(gtf)
    temp_files = []
    
    if genome_path != genome:
        temp_files.append(genome_path)
    if gtf_path != gtf:
        temp_files.append(gtf_path)
    
    try:
        processed_samples = []
        
        for fastq_file in fastq_files:
            filename = Path(fastq_file).stem.replace('.fastq', '').replace('.gz', '')
            log(f"Processing sample: {filename}")
            
            input_fastq = fastq_file
            output_unsorted_sam = os.path.join(output_dir, f"{filename}.sam")
            output_unsorted_bam = os.path.join(output_dir, f"{filename}.bam")
            output_sorted_bam = os.path.join(output_dir, f"{filename}_sorted.bam")
            output_longgf_total = os.path.join(output_dir, f"{filename}.log")
            
            # 1. Minimap2 alignment (support gzipped fastq)
            log("Running Minimap2 alignment...")
            if fastq_file.endswith('.gz'):
                cmd = f"gzip -dc {input_fastq} | minimap2 -ax splice -uf -k14 --secondary=no -G 50k -t {threads} {genome_path} - -o {output_unsorted_sam}"
            else:
                cmd = f"minimap2 -ax splice -uf -k14 --secondary=no -G 50k -t {threads} {genome_path} {input_fastq} -o {output_unsorted_sam}"
            run_command(cmd)
            
            # 2. Convert SAM to BAM
            log("Converting SAM to BAM...")
            run_command(f"samtools view -S -b -@ {threads} {output_unsorted_sam} -o {output_unsorted_bam}")
            
            # 3. Sort BAM by name
            log("Sorting BAM by name...")
            run_command(f"samtools sort -n -@ {threads} {output_unsorted_bam} -o {output_sorted_bam}")
            
            # 4. Run LongGF
            log("Running LongGF...")
            run_command(f"LongGF {output_sorted_bam} {gtf_path} 100 50 100 2 0 1 0 > {output_longgf_total}")
            
            processed_samples.append({
                'sample': filename,
                'sam_file': output_unsorted_sam,
                'bam_file': output_unsorted_bam,
                'sorted_bam': output_sorted_bam,
                'log_file': output_longgf_total
            })
            
            log(f"Completed processing for sample: {filename}")
        
        # Post-process all LongGF logs and generate Excel files using the R script
        log("Post-processing LongGF results...")
        postprocess(output_dir)
        
        # Cleanup intermediate files if requested
        if not keep_intermediate:
            log("Cleaning up intermediate files...")
            
            # Remove _results.txt files
            results_files = glob.glob(os.path.join(output_dir, '*_results.txt'))
            for file_path in results_files:
                try:
                    os.remove(file_path)
                    log(f"Removed: {file_path}")
                except Exception as e:
                    log(f"Warning: Could not remove {file_path}: {e}")
            
            # Remove _sorted.bam files
            sorted_bam_files = glob.glob(os.path.join(output_dir, '*_sorted.bam'))
            for file_path in sorted_bam_files:
                try:
                    os.remove(file_path)
                    log(f"Removed: {file_path}")
                except Exception as e:
                    log(f"Warning: Could not remove {file_path}: {e}")
        
        # Cleanup temporary files
        for temp_file in temp_files:
            try:
                os.remove(temp_file)
                log(f"Removed temp file: {temp_file}")
            except Exception as e:
                log(f"Warning: Could not remove temp file {temp_file}: {e}")
        
        log("LongGF processing completed successfully")
        
        return {
            'samples_processed': len(processed_samples),
            'output_directory': output_dir,
            'excel_outputs': [
                os.path.join(output_dir, 'Combined_LongGF_chimera_results_total.xlsx'),
                os.path.join(output_dir, 'Combined_LongGF_chimera_results_with_sample_info.xlsx')
            ],
            'sam_files': [s['sam_file'] for s in processed_samples]
        }
        
    except Exception as e:
        # Cleanup temp files on error
        for temp_file in temp_files:
            try:
                os.remove(temp_file)
            except:
                pass
        log(f"Error in LongGF processing: {e}")
        raise 