#!/usr/bin/env python3
"""
JaffaL Module for TYPHON Pipeline

Converts FASTQ files to FASTA, runs JaffaL via bpipe, and integrates results
with LongGF and Genion outputs using the original R scripts.

This module converts the bash workflow from the original TYPHON Main.sh script
(lines 120-170) to Python while preserving the exact same functionality.

Authors: Harry Kane, PhD; Eren Ada, PhD
"""

import os
import sys
import logging
import subprocess
import shutil
import glob
from pathlib import Path


def convert_fastq_to_fasta(fastq_file, output_dir):
    """
    Convert FASTQ to FASTA using sed equivalent.
    
    Original bash command:
    sed -n '1~4s:^@:>:p;2~4p' "$file" >> "$file.fasta"
    
    Args:
        fastq_file: Path to input FASTQ file
        output_dir: Directory to write FASTA file
    
    Returns:
        Path to generated FASTA file
    """
    fastq_path = Path(fastq_file)
    fasta_file = os.path.join(output_dir, f"{fastq_path.stem}.fasta")
    
    logging.info(f"Converting {fastq_path.name} to FASTA format...")
    
    try:
        with open(fastq_file, 'r') as infile, open(fasta_file, 'w') as outfile:
            lines = infile.readlines()
            
            # Process every 4th line starting from line 1 (header) and line 2 (sequence)
            for i in range(0, len(lines), 4):
                if i < len(lines):
                    # Header line: replace @ with >
                    header = lines[i].strip()
                    if header.startswith('@'):
                        header = '>' + header[1:]
                    outfile.write(header + '\n')
                    
                    # Sequence line
                    if i + 1 < len(lines):
                        sequence = lines[i + 1].strip()
                        outfile.write(sequence + '\n')
        
        logging.info(f"Created FASTA file: {fasta_file}")
        return fasta_file
        
    except Exception as e:
        logging.error(f"Failed to convert {fastq_file} to FASTA: {e}")
        raise


def run_bpipe_jaffal(fasta_file, jaffal_dir, threads=16, max_memory="24G"):
    """
    Run bpipe with JAFFAL.groovy for a single FASTA file.
    
    Original bash command:
    $8/tools/bin/bpipe run $8/JAFFAL.groovy $input_fastq_fasta
    
    Args:
        fasta_file: Path to input FASTA file
        jaffal_dir: Path to JaffaL installation directory
        threads: Number of threads to use
        max_memory: Maximum memory for Java process (e.g., "28G")
    
    Returns:
        Path to results directory for this sample
    """
    fasta_path = Path(fasta_file)
    sample_name = fasta_path.stem
    
    # Change to results directory (bpipe creates subdirectories here)
    results_dir = os.path.join(jaffal_dir, 'results')
    jaffal_groovy = os.path.abspath(os.path.join(jaffal_dir, 'JAFFAL.groovy'))
    bpipe_cmd = os.path.abspath(os.path.join(jaffal_dir, 'tools', 'bin', 'bpipe'))
    
    logging.info(f"Running JaffaL for sample {sample_name}...")
    logging.debug(f"bpipe command: {bpipe_cmd}")
    logging.debug(f"JAFFAL groovy: {jaffal_groovy}")
    logging.debug(f"Results directory: {results_dir}")
    
    try:
        # Run bpipe command with threading parameter - use absolute path for fasta_file
        abs_fasta_file = os.path.abspath(fasta_file)
        cmd = [bpipe_cmd, 'run', '-p', f'threads={threads}', jaffal_groovy, abs_fasta_file]
        logging.info(f"Running command: {' '.join(cmd)}")
        logging.info(f"Working directory: {results_dir}")
        logging.info(f"Setting MAX_JAVA_MEM to: {max_memory}")
        
        # Set environment variable to override bpipe's default memory setting
        env = os.environ.copy()
        env['MAX_JAVA_MEM'] = max_memory
        
        result = subprocess.run(
            cmd,
            cwd=results_dir,
            capture_output=True,
            text=True,
            check=True,
            env=env
        )
        
        logging.info(f"JaffaL completed for {sample_name}")
        logging.debug(f"bpipe stdout: {result.stdout}")
        
        # Return path to sample-specific results directory
        sample_results_dir = os.path.join(results_dir, sample_name)
        return sample_results_dir
        
    except subprocess.CalledProcessError as e:
        logging.error(f"bpipe failed for {sample_name}: {e}")
        logging.error(f"stdout: {e.stdout}")
        logging.error(f"stderr: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error running bpipe for {sample_name}: {e}")
        raise


def aggregate_jaffal_results(jaffal_dir, output_dir):
    """
    Aggregate JaffaL results from all samples into a single file.
    
    Original bash command:
    find $8/results -mindepth 2 -type f -name "*.fastq.summary" -exec cat {} \; > JaffaL_combined_results.txt
    
    Args:
        jaffal_dir: Path to JaffaL installation directory
        output_dir: Directory to write combined results
    
    Returns:
        Path to combined results file
    """
    results_dir = os.path.join(jaffal_dir, 'results')
    combined_file = os.path.join(output_dir, 'JaffaL_combined_results.txt')
    
    logging.info("Aggregating JaffaL results from all samples...")
    
    try:
        # Find all .summary files in subdirectories
        summary_files = []
        for root, dirs, files in os.walk(results_dir):
            # Only look in subdirectories (mindepth 2 equivalent)
            if root != results_dir:
                for file in files:
                    if file.endswith('.summary'):
                        summary_files.append(os.path.join(root, file))
        
        if not summary_files:
            logging.warning("No .summary files found in JaffaL results")
            # Create empty file
            with open(combined_file, 'w') as f:
                f.write("")
            return combined_file
        
        logging.info(f"Found {len(summary_files)} summary files to combine")
        
        # Concatenate all summary files
        with open(combined_file, 'w') as outfile:
            for summary_file in summary_files:
                logging.debug(f"Adding {summary_file} to combined results")
                with open(summary_file, 'r') as infile:
                    outfile.write(infile.read())
        
        logging.info(f"Created combined JaffaL results: {combined_file}")
        return combined_file
        
    except Exception as e:
        logging.error(f"Failed to aggregate JaffaL results: {e}")
        raise





def run_jaffal(fastq_dir, jaffal_dir, output_dir, threads=1, keep_intermediate=False, config=None):
    """
    Main JaffaL execution function.
    
    Converts the bash workflow from Main.sh (Step 3. Running JaffaL) to Python.
    This module only handles JaffaL execution - integration is handled separately.
    
    Args:
        fastq_dir: Directory containing input FASTQ files
        jaffal_dir: Path to JaffaL installation directory  
        output_dir: Main pipeline output directory
        threads: Number of threads (for future use)
        keep_intermediate: Whether to preserve temporary files
        config: Configuration dictionary (optional, for memory management settings)
    
    Returns:
        Dictionary with paths to generated files
    """
    logging.info("=" * 50)
    logging.info("Starting JaffaL Analysis")
    logging.info("=" * 50)
    
    # Extract JaffaL-specific configuration
    jaffal_config = {}
    if config:
        jaffal_config = config.get('modules', {}).get('jaffal', {})
    
    # Memory management settings
    process_sequentially = jaffal_config.get('process_samples_sequentially', True)
    max_memory = jaffal_config.get('max_memory', '20G')  # Reduced to 20G
    bpipe_memory = jaffal_config.get('bpipe_memory', '20G')  # Use bpipe_memory for Java heap
    
    # Thread configuration - use config value, fallback to parameter, then default
    config_threads = jaffal_config.get('threads', threads)
    if config_threads != threads:
        logging.info(f"Using threads from config: {config_threads} (parameter was: {threads})")
        threads = config_threads
    
    logging.info(f"Processing mode: {'Sequential' if process_sequentially else 'Parallel'}")
    logging.info(f"Memory limits: max={max_memory}, bpipe={bpipe_memory}")
    logging.info(f"Thread configuration: {threads}")
    
    # Validate inputs
    if not os.path.exists(fastq_dir):
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")
    
    if not os.path.exists(jaffal_dir):
        raise FileNotFoundError(f"JaffaL directory not found: {jaffal_dir}")
    
    # Create JaffaL output directories
    jaffal_output_dir = os.path.join(output_dir, 'jaffal_results')
    os.makedirs(jaffal_output_dir, exist_ok=True)
    
    # Step 1: Create temporary fasta_files directory
    fasta_files_dir = os.path.join(jaffal_dir, 'fasta_files')
    if os.path.exists(fasta_files_dir):
        shutil.rmtree(fasta_files_dir)
    os.makedirs(fasta_files_dir)
    
    logging.info(f"Created temporary directory: {fasta_files_dir}")
    
    try:
        # Step 2: Copy FASTQ files to fasta_files directory
        fastq_files = glob.glob(os.path.join(fastq_dir, '*.fastq')) + \
                     glob.glob(os.path.join(fastq_dir, '*.fq'))
        
        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
        
        logging.info(f"Found {len(fastq_files)} FASTQ files to process")
        
        copied_fastq_files = []
        for fastq_file in fastq_files:
            dest_file = os.path.join(fasta_files_dir, os.path.basename(fastq_file))
            shutil.copy2(fastq_file, dest_file)
            copied_fastq_files.append(dest_file)
        
        # Step 3: Convert FASTQ files to FASTA
        fasta_files = []
        for fastq_file in copied_fastq_files:
            fasta_file = convert_fastq_to_fasta(fastq_file, fasta_files_dir)
            fasta_files.append(fasta_file)
        
        # Remove FASTQ files from temporary directory
        for fastq_file in copied_fastq_files:
            os.remove(fastq_file)
        
        # Step 4: Run bpipe for each FASTA file
        sample_results = []
        
        if process_sequentially:
            # Sequential processing with memory management
            logging.info("Processing samples sequentially for memory management")
            for fasta_file in fasta_files:
                try:
                    # Clean up any previous results
                    results_dir = os.path.join(jaffal_dir, 'results')
                    if os.path.exists(results_dir):
                        shutil.rmtree(results_dir)
                    os.makedirs(results_dir)
                    
                    # Process single sample
                    sample_result_dir = run_bpipe_jaffal(fasta_file, jaffal_dir, threads, bpipe_memory)
                    sample_results.append(sample_result_dir)
                    
                    # Clean up FASTA file after processing
                    os.remove(fasta_file)
                    
                    # Force garbage collection and cleanup
                    import gc
                    import time
                    gc.collect()
                    
                    # Kill any lingering Java processes from bpipe
                    subprocess.run(['pkill', '-f', 'bpipe'], capture_output=True)
                    subprocess.run(['pkill', '-f', 'java.*bpipe'], capture_output=True)
                    
                    # Brief pause to ensure cleanup
                    time.sleep(2)
                    
                    # Check memory status after cleanup
                    try:
                        import psutil
                        memory_info = psutil.virtual_memory()
                        logging.info(f"Completed sample {Path(fasta_file).stem} - Available memory: {memory_info.available / (1024**3):.1f}GB")
                    except ImportError:
                        logging.info(f"Completed sample {Path(fasta_file).stem}, memory and processes cleaned up")
                    
                except Exception as e:
                    logging.error(f"Failed to process sample {fasta_file}: {e}")
                    raise
        else:
            # Parallel processing (original behavior)
            logging.info("Processing samples in parallel")
            for fasta_file in fasta_files:
                sample_result_dir = run_bpipe_jaffal(fasta_file, jaffal_dir, threads, bpipe_memory)
                sample_results.append(sample_result_dir)
        
        # Step 5: Aggregate results from all samples
        combined_results_file = aggregate_jaffal_results(jaffal_dir, jaffal_output_dir)
        
        logging.info("=" * 50)
        logging.info("JaffaL Analysis Completed Successfully")
        logging.info(f"Combined results: {combined_results_file}")
        logging.info("=" * 50)
        
        return {
            'combined_results': combined_results_file,
            'sample_results': sample_results
        }
        
    finally:
        # Cleanup temporary directory
        if not keep_intermediate and os.path.exists(fasta_files_dir):
            shutil.rmtree(fasta_files_dir)
            logging.info(f"Cleaned up temporary directory: {fasta_files_dir}")


if __name__ == "__main__":
    # Test function - can be used for debugging
    import argparse
    
    parser = argparse.ArgumentParser(description="Run JaffaL analysis")
    parser.add_argument('--fastq-dir', required=True, help='Directory with FASTQ files')
    parser.add_argument('--jaffal-dir', required=True, help='JaffaL installation directory')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--longgf-excel', help='LongGF Excel results file (not used in this module)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    run_jaffal(
        fastq_dir=args.fastq_dir,
        jaffal_dir=args.jaffal_dir,
        output_dir=args.output_dir,
        keep_intermediate=args.debug
    ) 