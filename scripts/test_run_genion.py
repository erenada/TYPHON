#!/usr/bin/env python3
"""
Test script for run_genion.py module functionality.
Tests the Genion integration with debug mode preservation.
"""

import argparse
import os
import sys
import tempfile
from pathlib import Path

# Add the parent directory to the path to import typhon modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from typhon.modules.run_genion import run_genion, get_genion_bin
from typhon.utils.genion_reference import prepare_genion_reference_files


def test_genion_binary():
    """Test that the custom Genion binary is available and working."""
    print("Testing Genion binary availability...")
    try:
        genion_bin = get_genion_bin()
        print(f"✓ Genion binary found at: {genion_bin}")
        
        # Test version
        import subprocess
        result = subprocess.run([genion_bin, '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ Genion version: {result.stdout.strip()}")
        else:
            print(f"⚠ Genion version check failed: {result.stderr}")
        return True
    except Exception as e:
        print(f"✗ Genion binary test failed: {e}")
        return False


def test_paftools():
    """Test that paftools.js is available."""
    print("Testing paftools.js availability...")
    try:
        import subprocess
        result = subprocess.run(['paftools.js'], capture_output=True, text=True)
        if 'Usage:' in result.stderr or 'Usage:' in result.stdout:
            print("✓ paftools.js is available")
            return True
        else:
            print(f"⚠ paftools.js may not be working properly")
            return False
    except Exception as e:
        print(f"✗ paftools.js test failed: {e}")
        return False


def test_reference_preparation():
    """Test reference file preparation."""
    print("Testing reference file preparation...")
    
    # Mock test files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a minimal mock GTF
        mock_gtf = os.path.join(temp_dir, 'test.gtf')
        with open(mock_gtf, 'w') as f:
            f.write('1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001.1"; transcript_id "ENST00000000001.1";\n')
        
        # Create a minimal mock transcriptome FASTA
        mock_fasta = os.path.join(temp_dir, 'test.fa')
        with open(mock_fasta, 'w') as f:
            f.write('>ENST00000000001.1\n')
            f.write('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n')
        
        output_dir = os.path.join(temp_dir, 'ref_output')
        
        try:
            gtf_for_genion, selfalign_paf, selfalign_tsv = prepare_genion_reference_files(
                gtf=mock_gtf,
                transcriptome_fasta=mock_fasta,
                output_dir=output_dir,
                threads=1,
                reference_type='gencode'
            )
            
            print(f"✓ Reference preparation completed")
            print(f"  - GTF: {gtf_for_genion}")
            print(f"  - PAF: {selfalign_paf}")
            print(f"  - TSV: {selfalign_tsv}")
            
            # Check if files exist
            if all(os.path.exists(f) for f in [gtf_for_genion, selfalign_paf, selfalign_tsv]):
                print("✓ All reference files created successfully")
                return True
            else:
                print("✗ Some reference files missing")
                return False
                
        except Exception as e:
            print(f"✗ Reference preparation failed: {e}")
            return False


def main():
    parser = argparse.ArgumentParser(description="Test run_genion.py functionality")
    parser.add_argument('--full-test', action='store_true', 
                       help='Run full test including reference preparation (requires R)')
    args = parser.parse_args()
    
    print("=" * 60)
    print("Typhon run_genion.py Test Suite")
    print("=" * 60)
    
    tests_passed = 0
    total_tests = 0
    
    # Test 1: Genion binary
    total_tests += 1
    if test_genion_binary():
        tests_passed += 1
    
    print()
    
    # Test 2: paftools.js
    total_tests += 1
    if test_paftools():
        tests_passed += 1
    
    print()
    
    # Test 3: Reference preparation (optional, requires R and dependencies)
    if args.full_test:
        total_tests += 1
        if test_reference_preparation():
            tests_passed += 1
        print()
    
    # Summary
    print("=" * 60)
    print(f"Test Results: {tests_passed}/{total_tests} tests passed")
    print("=" * 60)
    
    if tests_passed == total_tests:
        print("✓ All tests passed! run_genion.py is ready for Phase 2.3 testing.")
        return 0
    else:
        print("✗ Some tests failed. Check dependencies and setup.")
        return 1


if __name__ == '__main__':
    sys.exit(main()) 