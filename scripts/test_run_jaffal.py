#!/usr/bin/env python3
"""
Test script for JaffaL module

Tests the refactored run_jaffal.py module with test data.
"""

import os
import sys
import logging
import tempfile
import shutil
from pathlib import Path

# Add typhon module to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from typhon.modules.run_jaffal import run_jaffal


def setup_logging():
    """Setup logging for the test."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def test_jaffal_basic():
    """Test basic JaffaL functionality."""
    setup_logging()
    
    # Define paths
    project_root = Path(__file__).parent.parent.resolve()
    fastq_dir = project_root / "test_data" / "FASTQ"
    output_dir = project_root / "test_output_pipeline"
    jaffal_dir = project_root / "jaffal" / "JAFFAL-version-2.3"
    
    logging.info("Testing JaffaL module with refactored code...")
    logging.info(f"FASTQ directory: {fastq_dir}")
    logging.info(f"JaffaL directory: {jaffal_dir}")
    logging.info(f"Output directory: {output_dir}")
    
    # Check if required directories exist
    if not fastq_dir.exists():
        logging.error(f"FASTQ directory not found: {fastq_dir}")
        return False
    
    # Debug JaffaL directory path
    logging.info(f"Current working directory: {os.getcwd()}")
    logging.info(f"Script file location: {__file__}")
    logging.info(f"Project root: {project_root}")
    logging.info(f"Checking JaffaL directory: {jaffal_dir}")
    logging.info(f"JaffaL directory resolved: {jaffal_dir.resolve()}")
    logging.info(f"JaffaL directory exists: {jaffal_dir.exists()}")
    logging.info(f"OS path exists: {os.path.exists(str(jaffal_dir))}")
    logging.info(f"Absolute path exists: {os.path.exists(str(jaffal_dir.resolve()))}")
    
    # Force the path to the known working absolute path for testing
    jaffal_dir = Path("/home/eren/Desktop/Typhon_pipeline/jaffal/JAFFA-version-2.3")
    
    if not os.path.exists(str(jaffal_dir)):
        logging.error(f"JaffaL directory not found: {jaffal_dir}")
        logging.info("Please run 'python setup_jaffal.py' first")
        return False
    else:
        logging.info(f"JaffaL directory confirmed to exist: {jaffal_dir}")
    
    # List FASTQ files
    fastq_files = list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq"))
    if not fastq_files:
        logging.error(f"No FASTQ files found in {fastq_dir}")
        return False
    
    logging.info(f"Found {len(fastq_files)} FASTQ files:")
    for f in fastq_files:
        logging.info(f"  - {f.name}")
    
    try:
        # Debug paths before calling run_jaffal
        fastq_path = str(fastq_dir.resolve())
        jaffal_path = str(jaffal_dir.resolve())
        output_path = str(output_dir.resolve())
        
        logging.info(f"About to call run_jaffal with:")
        logging.info(f"  fastq_dir: {fastq_path} (exists: {os.path.exists(fastq_path)})")
        logging.info(f"  jaffal_dir: {jaffal_path} (exists: {os.path.exists(jaffal_path)})")
        logging.info(f"  output_dir: {output_path} (exists: {os.path.exists(output_path)})")
        
        # Run JaffaL (using same threading as old TYPHON)
        results = run_jaffal(
            fastq_dir=fastq_path,
            jaffal_dir=jaffal_path,
            output_dir=output_path,
            threads=30,  # Match old TYPHON threading
            keep_intermediate=True  # Keep files for debugging
        )
        
        logging.info("JaffaL test completed successfully!")
        logging.info("Results:")
        for key, value in results.items():
            logging.info(f"  {key}: {value}")
        
        # Verify output files exist
        combined_results = results.get('combined_results')
        if combined_results and os.path.exists(combined_results):
            logging.info(f"Combined results file created: {combined_results}")
            # Show file size
            size = os.path.getsize(combined_results)
            logging.info(f"File size: {size} bytes")
        else:
            logging.warning("Combined results file not found or empty")
        
        return True
        
    except Exception as e:
        logging.error(f"JaffaL test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_jaffal_basic()
    if success:
        print("\n✓ JaffaL test completed successfully")
        sys.exit(0)
    else:
        print("\n✗ JaffaL test failed")
        sys.exit(1) 