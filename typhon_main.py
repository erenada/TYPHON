#!/usr/bin/env python3
"""
TYPHON - Chimeric RNA Detection Pipeline

Main entry point for the TYPHON pipeline that integrates LongGF, Genion, and JaffaL
for comprehensive chimeric RNA detection from direct RNA sequencing data.

Author: Eren Ada, PhD
"""

import argparse
import sys
import os
import yaml
import logging
from pathlib import Path
from datetime import datetime

# Add the typhon module to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from typhon.modules.run_longgf import run_longgf
from typhon.modules.run_genion import run_genion


def setup_logging(log_file=None, debug=False):
    """Setup logging configuration."""
    level = logging.DEBUG if debug else logging.INFO
    format_str = '%(asctime)s - %(levelname)s - %(message)s'
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=handlers
    )


def load_config(config_file):
    """Load and validate YAML configuration file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logging.info(f"Loaded configuration from {config_file}")
        return config
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {config_file}")
        sys.exit(1)
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML configuration: {e}")
        sys.exit(1)


def validate_config(config):
    """Validate configuration parameters."""
    required_sections = ['project', 'input', 'references']
    
    for section in required_sections:
        if section not in config:
            logging.error(f"Missing required configuration section: {section}")
            return False
    
    # Validate required input files
    required_files = {
        'fastq_dir': config['input']['fastq_dir'],
        'genome': config['references']['genome'],
        'gtf': config['references']['gtf']
    }
    
    for name, path in required_files.items():
        if not os.path.exists(path):
            logging.error(f"Required {name} path does not exist: {path}")
            return False
    
    logging.info("Configuration validation passed")
    return True



def run_longgf_step(config):
    """Execute LongGF module."""
    logging.info("Starting LongGF analysis...")
    
    longgf_config = config.get('modules', {}).get('longgf', {})
    
    try:
        sam_files = run_longgf(
            fastq_dir=config['input']['fastq_dir'],
            genome=config['references']['genome'],
            gtf=config['references']['gtf'],
            output_dir=config['project']['output_dir'],
            threads=config['project'].get('threads', 4),
            keep_intermediate=longgf_config.get('keep_intermediate', False)
        )
        
        logging.info(f"LongGF completed successfully. Generated {len(sam_files)} SAM files")
        return sam_files
        
    except Exception as e:
        logging.error(f"LongGF failed: {e}")
        raise


def run_genion_step(config, sam_files):
    """Execute Genion module."""
    logging.info("Starting Genion analysis...")
    
    genion_config = config.get('modules', {}).get('genion', {})
    
    try:
        # TODO: Need to prepare Genion reference files (selfalign_paf, selfalign_tsv)
        # This will be implemented when we integrate the genion_reference.py utilities
        
        logging.warning("Genion integration not fully implemented yet")
        logging.info("Required for full Genion integration:")
        logging.info("  - Self-alignment PAF generation")
        logging.info("  - GTF conversion for Genion format")
        logging.info("  - Genomic superduplicates file")
        
        # Placeholder implementation
        results = []
        for sam_file in sam_files:
            sample_name = Path(sam_file).stem
            logging.info(f"Would process {sample_name} with Genion")
            # TODO: Implement actual Genion execution
            results.append(f"genion_result_{sample_name}.tsv")
        
        logging.info(f"Genion step completed (placeholder) for {len(results)} samples")
        return results
        
    except Exception as e:
        logging.error(f"Genion failed: {e}")
        raise


def run_jaffal_step(config):
    """Execute JaffaL module (placeholder)."""
    logging.info("JaffaL analysis not yet implemented")
    # TODO: Implement JaffaL integration
    pass


def create_output_directory(output_dir):
    """Create output directory structure."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    subdirs = ['longgf_results', 'genion_results', 'jaffal_results', 'logs']
    for subdir in subdirs:
        Path(output_dir, subdir).mkdir(exist_ok=True)
    
    logging.info(f"Created output directory: {output_dir}")


def main():
    """Main pipeline execution function."""
    parser = argparse.ArgumentParser(
        description="TYPHON - Chimeric RNA Detection Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default configuration file (config.yaml)
  python typhon_main.py

  # Run with custom configuration file
  python typhon_main.py --config my_config.yaml

  # Run with config and override threads
  python typhon_main.py --threads 30

  # Run only LongGF module
  python typhon_main.py --modules longgf
        """
    )
    
    parser.add_argument('--config', '-c', default='config.yaml',
                        help='YAML configuration file (default: config.yaml)')
    parser.add_argument('--threads', '-t', type=int,
                        help='Number of threads (overrides config)')
    parser.add_argument('--output', '-o',
                        help='Output directory (overrides config)')
    parser.add_argument('--modules', nargs='+', 
                        choices=['longgf', 'genion', 'jaffal'],
                        help='Run specific modules only')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be executed without running')
    
    args = parser.parse_args()
    
    # Check if config file exists
    if not os.path.exists(args.config):
        parser.error(f"Configuration file '{args.config}' not found. Please ensure it exists or specify a different path with --config.")
    
    # Load configuration
    config = load_config(args.config)
    
    # Apply command line overrides
    if args.threads:
        config['project']['threads'] = args.threads
    if args.output:
        config['project']['output_dir'] = args.output
    if args.debug:
        config['options']['debug'] = True
    
    # Setup basic logging first
    setup_logging(debug=config['options'].get('debug', False))
    
    # Setup file logging after output directory is created
    os.makedirs(config['project']['output_dir'], exist_ok=True)
    log_file = os.path.join(config['project']['output_dir'], 'logs', 'typhon.log')
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    setup_logging(log_file, config['options'].get('debug', False))
    
    # Log pipeline start
    logging.info("=" * 60)
    logging.info("TYPHON Pipeline Starting")
    logging.info(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Project: {config['project']['name']}")
    logging.info(f"Output: {config['project']['output_dir']}")
    logging.info(f"Threads: {config['project']['threads']}")
    logging.info("=" * 60)
    
    # Validate configuration
    if not validate_config(config):
        sys.exit(1)
    
    # Create output directory
    create_output_directory(config['project']['output_dir'])
    
    if args.dry_run:
        logging.info("DRY RUN - Would execute the following steps:")
        enabled_modules = []
        if args.modules:
            enabled_modules = args.modules
        else:
            for module, settings in config.get('modules', {}).items():
                if settings.get('enabled', False):
                    enabled_modules.append(module)
        
        for module in enabled_modules:
            logging.info(f"  - {module.upper()} analysis")
        
        logging.info("Use --no-dry-run to actually run the pipeline")
        sys.exit(0)
    
    try:
        sam_files = []
        
        # Determine which modules to run
        modules_to_run = args.modules or []
        if not modules_to_run:
            # Run all enabled modules from config
            for module, settings in config.get('modules', {}).items():
                if settings.get('enabled', False):
                    modules_to_run.append(module)
        
        # Execute modules in order
        if 'longgf' in modules_to_run:
            sam_files = run_longgf_step(config)
        
        if 'genion' in modules_to_run:
            if not sam_files:
                # Look for existing SAM files if LongGF wasn't run
                output_dir = config['project']['output_dir']
                sam_files = list(Path(output_dir).glob("*.sam"))
                if not sam_files:
                    logging.error("No SAM files found for Genion. Run LongGF first or provide SAM files.")
                    sys.exit(1)
            
            run_genion_step(config, sam_files)
        
        if 'jaffal' in modules_to_run:
            run_jaffal_step(config)
        
        # Pipeline completion
        logging.info("=" * 60)
        logging.info("TYPHON Pipeline completed successfully!")
        logging.info(f"Results saved to: {config['project']['output_dir']}")
        logging.info(f"Log file: {log_file}")
        logging.info("=" * 60)
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        if config['options'].get('debug', False):
            import traceback
            logging.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main() 