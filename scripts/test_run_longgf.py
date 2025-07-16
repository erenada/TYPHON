#!/usr/bin/env python3
"""
Test script for run_longgf.py module functionality.
Tests the LongGF integration with real data from typhon_old.
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

# Add the parent directory to the path to import typhon modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from typhon.modules.run_longgf import run_longgf


def test_longgf_basic():
    """Test that LongGF dependencies are available."""
    print("Testing LongGF dependencies...")
    try:
        import subprocess
        
        # Test minimap2
        result = subprocess.run(['minimap2', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ minimap2 version: {result.stdout.strip()}")
        else:
            print(f"✗ minimap2 not working properly")
            return False
            
        # Test samtools
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            print(f"✓ {version_line}")
        else:
            print(f"✗ samtools not working properly")
            return False
            
        # Test LongGF
        result = subprocess.run(['LongGF'], capture_output=True, text=True)
        if 'Usage:' in result.stderr or 'Usage:' in result.stdout:
            print("✓ LongGF is available")
        else:
            print(f"⚠ LongGF may not be working properly")
            
        # Test R and required packages
        r_test_script = """
        required_packages <- c("data.table", "magrittr", "openxlsx", "stringr", "tidyr", "dplyr")
        missing_packages <- required_packages[!sapply(required_packages, require, character.only=TRUE, quietly=TRUE)]
        if (length(missing_packages) > 0) {
            cat("Missing R packages:", paste(missing_packages, collapse=", "), "\\n")
            quit(status=1)
        } else {
            cat("All required R packages are available\\n")
        }
        """
        
        result = subprocess.run(['Rscript', '-e', r_test_script], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ R and required packages are available")
        else:
            print(f"✗ R packages missing: {result.stdout.strip()}")
            return False
            
        return True
        
    except Exception as e:
        print(f"✗ Dependency test failed: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Test the LongGF step of Typhon with real data.")
    parser.add_argument('--fastq_dir', required=True, help='Directory with input FASTQ files')
    parser.add_argument('--genome', required=True, help='Reference genome fasta')
    parser.add_argument('--gtf', required=True, help='Reference GTF file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--keep_intermediate', action='store_true', help='Keep intermediate files')
    parser.add_argument('--test_deps_only', action='store_true', help='Only test dependencies, do not run pipeline')
    args = parser.parse_args()

    # Test dependencies first
    if not test_longgf_basic():
        print("❌ Dependency tests failed. Please check your environment.")
        return 1
    
    if args.test_deps_only:
        print("✅ All dependency tests passed!")
        return 0

    # Validate input files
    if not os.path.isdir(args.fastq_dir):
        print(f"❌ FASTQ directory not found: {args.fastq_dir}")
        return 1
        
    if not os.path.isfile(args.genome):
        print(f"❌ Genome file not found: {args.genome}")
        return 1
        
    if not os.path.isfile(args.gtf):
        print(f"❌ GTF file not found: {args.gtf}")
        return 1

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    log_file = os.path.join(args.output_dir, 'test_longgf.log')

    print(f"🚀 Starting LongGF integration test...")
    print(f"   FASTQ dir: {args.fastq_dir}")
    print(f"   Genome: {args.genome}")
    print(f"   GTF: {args.gtf}")
    print(f"   Output: {args.output_dir}")
    print(f"   Threads: {args.threads}")
    print(f"   Log: {log_file}")

    try:
        # Run LongGF pipeline
        results = run_longgf(
            fastq_dir=args.fastq_dir,
            genome=args.genome,
            gtf=args.gtf,
            output_dir=args.output_dir,
            threads=args.threads,
            keep_intermediate=args.keep_intermediate,
            log_path=log_file
        )
        
        print("✅ LongGF pipeline completed successfully!")
        print(f"   Samples processed: {results['samples_processed']}")
        print(f"   SAM files created: {len(results['sam_files'])}")
        
        # List SAM files for Genion integration
        print("\n📁 SAM files available for Genion integration:")
        for sam_file in results['sam_files']:
            if os.path.exists(sam_file):
                size_mb = os.path.getsize(sam_file) / (1024*1024)
                print(f"   ✓ {sam_file} ({size_mb:.1f} MB)")
            else:
                print(f"   ✗ {sam_file} (missing)")
        
        # Check Excel outputs
        print("\n📊 Excel outputs:")
        for excel_file in results['excel_outputs']:
            if os.path.exists(excel_file):
                size_kb = os.path.getsize(excel_file) / 1024
                print(f"   ✓ {excel_file} ({size_kb:.1f} KB)")
            else:
                print(f"   ✗ {excel_file} (missing)")
        
        return 0
        
    except Exception as e:
        print(f"❌ LongGF pipeline failed: {e}")
        print(f"   Check log file: {log_file}")
        return 1


if __name__ == '__main__':
    exit(main()) 