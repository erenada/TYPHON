import os
import subprocess
from pathlib import Path

def prepare_genion_reference_files(gtf, transcriptome_fasta, output_dir, threads=1, reference_type='gencode', log=None):
    """
    Prepare reference files for Genion (GTF conversion, self-align PAF, self-align TSV).
    Returns (gtf_for_genion, selfalign_paf, selfalign_tsv)
    """
    os.makedirs(output_dir, exist_ok=True)
    def _log(msg):
        if log:
            log(msg)
        else:
            print(msg)
    # 1. GTF conversion if needed
    if reference_type == 'gencode':
        gtf_mod1 = os.path.join(output_dir, 'Genion_modified_gtf_file_1.gtf')
        gtf_mod2 = os.path.join(output_dir, 'Genion_modified_gtf_file_2.gtf')
        gtf_final = os.path.join(output_dir, 'Genion_modified_gtf_final.gtf')
        if not os.path.exists(gtf_mod2):
            _log('Converting Gencode GTF to Genion-compatible format (utility)...')
            sed1 = f"sed -E 's/gene_id \"(ENS(MUS)?[GTE][0-9]+)\\.([0-9]+)\";/gene_id \"\\1\"; gene_version \"\\3\";/' {gtf} > {gtf_mod1}"
            sed2 = f"sed -E 's/transcript_id \"(ENS(MUS)?[GTE][0-9]+)\\.([0-9]+)\";/transcript_id \"\\1\"; transcript_version \"\\3\";/' {gtf_mod1} > {gtf_mod2}"
            sed3 = f"sed -E 's/exon_id \"(ENS(MUS)?[GTE][0-9]+)\\.([0-9]+)\";/exon_id \"\\1\"; exon_version \"\\3\";/' {gtf_mod2} > {gtf_final}"
            sed4 = f"sed 's/^chrM/MT/;s/^chrX/X/;s/^chrY/Y/;s/^chr//' {gtf_final} > {gtf_mod1}"
            gtftk_cmd = f"gtftk convert_ensembl -i {gtf_mod1} -o {gtf_mod2}"
            subprocess.run(sed1, shell=True, check=True)
            subprocess.run(sed2, shell=True, check=True)
            subprocess.run(sed3, shell=True, check=True)
            subprocess.run(sed4, shell=True, check=True)
            subprocess.run(gtftk_cmd, shell=True, check=True)
        gtf_for_genion = gtf_mod2
    elif reference_type == 'ensembl':
        _log('Using Ensembl GTF as-is (utility).')
        gtf_for_genion = gtf
    else:
        raise ValueError("reference_type must be 'gencode' or 'ensembl'")
    # 2. Generate self-align PAF with minimap2
    selfalign_paf = os.path.join(output_dir, 'selfalign.paf')
    if not os.path.exists(selfalign_paf):
        _log('Generating self-align PAF with minimap2 (utility)...')
        minimap2_cmd = f"minimap2 {transcriptome_fasta} {transcriptome_fasta} -t {threads} -X -c -o {selfalign_paf}"
        subprocess.run(minimap2_cmd, shell=True, check=True)
    # 3. Run Genion_selfalign_paf_generation.R to produce TSV
    selfalign_tsv = os.path.join(output_dir, 'selfalign.tsv')
    if not os.path.exists(selfalign_tsv):
        _log('Running Genion_selfalign_paf_generation.R (utility)...')
        rscript_path = os.path.join(os.path.dirname(__file__), '../modules/Genion_selfalign_paf_generation.R')
        rscript_cmd = f"Rscript {rscript_path} {selfalign_paf} {selfalign_tsv}"
        subprocess.run(rscript_cmd, shell=True, check=True)
    return gtf_for_genion, selfalign_paf, selfalign_tsv 