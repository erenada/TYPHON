#!/usr/bin/env python3
"""
TYPHON - Chimeric RNA Detection Pipeline

Main entry point for the TYPHON pipeline that integrates LongGF, Genion, and JaffaL
for comprehensive chimeric RNA detection from direct RNA sequencing data.

Authors: Harry Kane, PhD; Eren Ada, PhD
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
from typhon.modules.exon_repair import run_exon_repair
from typhon.command_utils import setup_module_logger


def parse_args():
    """Parse command line arguments."""
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
    parser.add_argument('--debug', action='store_true', default=True,
                        help='Enable debug logging (enabled by default)')
    parser.add_argument('--no-debug', action='store_true',
                        help='Disable debug logging')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be executed without running')
    
    return parser.parse_args()


def setup_logging(log_file=None, debug=False):
    """Setup logging configuration."""
    level = logging.DEBUG if debug else logging.INFO
    format_str = '%(asctime)s - %(levelname)s - %(message)s'
    
    # Get root logger
    root_logger = logging.getLogger()
    
    # Clear existing handlers to avoid duplicates
    root_logger.handlers.clear()
    
    # Create handlers
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    # Set up formatter and add handlers
    formatter = logging.Formatter(format_str)
    for handler in handlers:
        handler.setFormatter(formatter)
        root_logger.addHandler(handler)
    
    root_logger.setLevel(level)


def load_config(config_file):
    """Load and validate YAML configuration file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        logging.info(f"Loaded configuration from {config_file}")
        
        # Convert relative paths to absolute paths for robustness
        config = resolve_config_paths(config)
        
        return config
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {config_file}")
        sys.exit(1)
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML configuration: {e}")
        sys.exit(1)


def resolve_config_paths(config):
    """Convert relative paths in config to absolute paths for robustness."""
    # Get the directory containing the config file (project root)
    project_root = os.path.dirname(os.path.abspath('config.yaml'))
    
    # Paths that should be resolved to absolute paths
    path_fields = [
        ('project', 'output_dir'),
        ('input', 'fastq_dir'),
        ('references', 'genome'),
        ('references', 'gtf'),
        ('references', 'transcriptome'),
        ('modules', 'jaffal', 'jaffal_dir'),
        ('modules', 'jaffal', 'reference_files', 'genome_fasta_gz'),
        ('modules', 'jaffal', 'reference_files', 'transcriptome_fasta'),
        ('modules', 'jaffal', 'reference_files', 'annotation_bed'),
        ('modules', 'jaffal', 'reference_files', 'annotation_tab'),
    ]
    
    for path_spec in path_fields:
        # Navigate to the nested config value
        current = config
        for key in path_spec[:-1]:
            if key in current:
                current = current[key]
            else:
                break
        else:
            # If we got here, all intermediate keys exist
            final_key = path_spec[-1]
            if final_key in current and current[final_key]:
                path_value = current[final_key]
                # Convert to absolute path if it's relative
                if not os.path.isabs(path_value):
                    current[final_key] = os.path.abspath(os.path.join(project_root, path_value))
                    logging.debug(f"Resolved path {'.'.join(path_spec)}: {path_value} -> {current[final_key]}")
    
    return config


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
    # Set up module-specific logging
    module_logger = setup_module_logger('longgf', config['project']['output_dir'])
    module_logger.info("Starting LongGF analysis...")
    
    longgf_config = config.get('modules', {}).get('longgf', {})
    
    try:
        # Use longgf_results subdirectory
        longgf_output_dir = os.path.join(config['project']['output_dir'], 'longgf_results')
        results = run_longgf(
            fastq_dir=config['input']['fastq_dir'],
            genome=config['references']['genome'],
            gtf=config['references']['gtf'],
            output_dir=longgf_output_dir,
            threads=config['project'].get('threads', 4),
            keep_intermediate=longgf_config.get('keep_intermediate', False),
            min_overlap_len=longgf_config.get('min_overlap_len', 100),
            bin_size=longgf_config.get('bin_size', 50),
            min_map_len=longgf_config.get('min_map_len', 100),
            pseudogene=longgf_config.get('pseudogene', 2),
            secondary_alignment=longgf_config.get('secondary_alignment', 0),
            min_sup_read=longgf_config.get('min_sup_read', 1),
            output_flag=longgf_config.get('output_flag', 0)
        )
        
        sam_files = results['sam_files']
        module_logger.info(f"LongGF completed successfully. Generated {len(sam_files)} SAM files")
        return sam_files
        
    except Exception as e:
        module_logger.error(f"LongGF failed: {e}")
        raise


def run_genion_step(config, sam_files):
    """Execute Genion module."""
    # Set up module-specific logging
    module_logger = setup_module_logger('genion', config['project']['output_dir'])
    module_logger.info("Starting Genion analysis...")
    
    genion_config = config.get('modules', {}).get('genion', {})
    references = config.get('references', {})
    project = config.get('project', {})
    
    try:
        from typhon.utils.genion_reference import prepare_genion_reference_files
        from typhon.modules.run_genion import run_genion
        
        # Prepare Genion reference files
        genion_ref_dir = os.path.join(project['output_dir'], 'genion_references')
        module_logger.info("Preparing Genion reference files...")
        
        gtf_for_genion, selfalign_paf, selfalign_tsv = prepare_genion_reference_files(
            gtf=references['gtf'],
            transcriptome_fasta=references['transcriptome'],
            output_dir=genion_ref_dir,
            threads=project.get('threads', 1),
            reference_type='gencode'
        )
        
        # Process each SAM file with Genion
        genion_output_dir = os.path.join(project['output_dir'], 'genion_results')
        results = []
        
        for sam_file in sam_files:
            sample_name = Path(sam_file).stem
            
            # Determine corresponding FASTQ file
            fastq_dir = config['input']['fastq_dir']
            fastq_file = None
            for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
                potential_fastq = os.path.join(fastq_dir, f"{sample_name}{ext}")
                if os.path.exists(potential_fastq):
                    fastq_file = potential_fastq
                    break
            
            if not fastq_file:
                module_logger.error(f"No FASTQ file found for sample {sample_name}")
                continue
                
            module_logger.info(f"Processing {sample_name} with Genion...")
            
            run_genion(
                input_fastq=fastq_file,
                input_sam=sam_file,
                gtf_for_genion=gtf_for_genion,
                selfalign_paf=selfalign_paf,
                selfalign_tsv=selfalign_tsv,
                output_dir=genion_output_dir,
                threads=genion_config.get('threads', 1),  # Use threads from config
                keep_intermediate=genion_config.get('keep_debug', True),
                log_path=os.path.join(genion_output_dir, 'run_genion.log'),
                min_support=genion_config.get('min_support', 1)
            )
            
            result_file = os.path.join(genion_output_dir, f'{sample_name}_genion.tsv')
            if os.path.exists(result_file):
                results.append(result_file)
                
        module_logger.info(f"Genion analysis completed for {len(results)} samples")
        return results
        
    except Exception as e:
        module_logger.error(f"Genion failed: {e}")
        raise


def run_jaffal_step(config, longgf_excel_file=None):
    """Execute JaffaL module."""
    # Set up module-specific logging
    module_logger = setup_module_logger('jaffal', config['project']['output_dir'])
    module_logger.info("Starting JaffaL analysis...")
    
    jaffal_config = config.get('modules', {}).get('jaffal', {})
    
    # Validate JaffaL is enabled
    if not jaffal_config.get('enabled', False):
        module_logger.warning("JaffaL is disabled in configuration")
        return []
    
    try:
        from typhon.modules.run_jaffal import run_jaffal
        
        # Get configuration parameters
        jaffal_dir = jaffal_config.get('jaffal_dir', './jaffal')
        fastq_dir = config['input']['fastq_dir']
        output_dir = config['project']['output_dir']
        
        # Find LongGF results file if not provided (prefer CSV, fallback to Excel)
        if not longgf_excel_file:
            longgf_results_dir = os.path.join(output_dir, 'longgf_results')
            if os.path.exists(longgf_results_dir):
                # Try CSV first (faster)
                csv_files = list(Path(longgf_results_dir).glob('Combined_LongGF_chimera_results_total.csv'))
                if csv_files:
                    longgf_excel_file = str(csv_files[0])
                else:
                    # Fallback to Excel
                    excel_files = list(Path(longgf_results_dir).glob('Combined_LongGF_chimera_results_total.xlsx'))
                    if excel_files:
                        longgf_excel_file = str(excel_files[0])
                    else:
                        # Try alternative locations (CSV first, then Excel)
                        for basename in ['Combined_LongGF_chimera_results_total.csv', 'Combined_LongGF_chimera_results_total.xlsx']:
                            for alt_dir in ['', 'test_output_longgf/', f"{output_dir}/"]:
                                alt_path = alt_dir + basename
                                if os.path.exists(alt_path):
                                    longgf_excel_file = alt_path
                                    break
                            if longgf_excel_file:
                                break
        
        if not longgf_excel_file or not os.path.exists(longgf_excel_file):
            module_logger.error("LongGF results file not found. Run LongGF first.")
            return []
        
        # Validate JaffaL installation
        if not os.path.exists(jaffal_dir):
            module_logger.error(f"JaffaL directory not found: {jaffal_dir}")
            module_logger.error("Run 'python setup_jaffal.py' first to install JaffaL")
            return []
        
        jaffal_executable = os.path.join(jaffal_dir, 'JAFFAL.groovy')
        if not os.path.exists(jaffal_executable):
            module_logger.error(f"JAFFAL.groovy not found: {jaffal_executable}")
            module_logger.error("JaffaL installation may be incomplete")
            return []
        
        module_logger.info(f"Using JaffaL installation: {jaffal_dir}")
        module_logger.info(f"Using LongGF results: {longgf_excel_file}")
        
        # Run JaffaL
        results = run_jaffal(
            fastq_dir=fastq_dir,
            jaffal_dir=jaffal_dir,
            output_dir=output_dir,
            threads=jaffal_config.get('threads', config['project'].get('threads', 4)),
            keep_intermediate=jaffal_config.get('keep_intermediate', False),
            config=config
        )
        
        module_logger.info("JaffaL analysis completed successfully")
        module_logger.info(f"Combined results: {results.get('combined_results')}")
        module_logger.info(f"Overlap analysis: {results.get('overlap_excel')}")
        
        return results
        
    except Exception as e:
        module_logger.error(f"JaffaL analysis failed: {e}")
        raise


def create_output_directory(output_dir):
    """Create output directory structure."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    subdirs = ['longgf_results', 'genion_results', 'jaffal_results', 'logs']
    for subdir in subdirs:
        Path(output_dir, subdir).mkdir(exist_ok=True)
    
    logging.info(f"Created output directory: {output_dir}")


def main():
    try:
        # Parse command line arguments
        args = parse_args()
        
        # Check if config file exists
        if not os.path.exists(args.config):
            print(f"Configuration file '{args.config}' not found. Please ensure it exists or specify a different path with --config.")
            sys.exit(1)
        
        # Load configuration
        config = load_config(args.config)
        
        # Apply command line overrides
        if args.threads:
            config['project']['threads'] = args.threads
        if args.output:
            config['project']['output_dir'] = args.output
        
        # Handle debug flags - --no-debug overrides --debug
        if args.no_debug:
            config.setdefault('options', {})['debug'] = False
        elif args.debug:
            config.setdefault('options', {})['debug'] = True
        
        # Validate configuration
        if not validate_config(config):
            sys.exit(1)
        
        # Create output directory structure
        output_dir = config['project']['output_dir']
        create_output_directory(output_dir)
        
        # Set up logging
        log_file = os.path.join(output_dir, 'logs', 'typhon_main.log')
        setup_logging(log_file=log_file, debug=config.get('options', {}).get('debug', False))
        logger = logging.getLogger(__name__)
        
        # Log pipeline start
        logger.info("=" * 60)
        logger.info("TYPHON Pipeline Starting")
        logger.info(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Project: {config['project'].get('name', 'Unnamed')}")
        logger.info(f"Output: {config['project']['output_dir']}")
        logger.info(f"Threads: {config['project']['threads']}")
        logger.info("=" * 60)
        
        # Handle dry run
        if args.dry_run:
            logger.info("DRY RUN - Would execute the following steps:")
            enabled_modules = []
            if args.modules:
                enabled_modules = args.modules
            else:
                for module, settings in config.get('modules', {}).items():
                    if settings.get('enabled', False):
                        enabled_modules.append(module)
            
            for module in enabled_modules:
                logger.info(f"  - {module.upper()} analysis")
            
            logger.info("Use without --dry-run to actually run the pipeline")
            return 0
        
        # Determine which modules to run
        modules_to_run = args.modules or []
        if not modules_to_run:
            # Run all enabled modules from config
            for module, settings in config.get('modules', {}).items():
                if settings.get('enabled', False):
                    modules_to_run.append(module)
        
        # Execute modules in order
        sam_files = []
        longgf_excel_file = None
        
        if 'longgf' in modules_to_run:
            logger.info("Running LongGF module")
            sam_files = run_longgf_step(config)
            # Find the generated LongGF results file (prefer CSV)
            # Search in multiple locations in order of preference
            search_locations = [
                os.path.join(output_dir, 'longgf_results'),  # Primary location
                output_dir,  # Fallback for backward compatibility
                'test_output_longgf'  # Legacy location
            ]
            
            for location in search_locations:
                if os.path.exists(location):
                    # Try CSV first (faster)
                    csv_files = list(Path(location).glob('Combined_LongGF_chimera_results_total.csv'))
                    if csv_files:
                        longgf_excel_file = str(csv_files[0])
                        break
                    
                    # Fallback to Excel
                    excel_files = list(Path(location).glob('Combined_LongGF_chimera_results_total.xlsx'))
                    if excel_files:
                        longgf_excel_file = str(excel_files[0])
                        break
        
        if 'genion' in modules_to_run:
            if not sam_files:
                # Look for existing SAM files if LongGF wasn't run
                sam_files = list(Path(output_dir).glob("*.sam"))
                
                # Also check common LongGF output locations
                if not sam_files:
                    for location in ['test_output_longgf', 'longgf_results', f"{output_dir}/longgf_results"]:
                        if os.path.exists(location):
                            sam_files = list(Path(location).glob("*.sam"))
                            if sam_files:
                                logger.info(f"Found existing SAM files in {location}")
                                break
                
                if not sam_files:
                    logger.error("No SAM files found for Genion. Run LongGF first or provide SAM files.")
                    sys.exit(1)
            
            logger.info("Running Genion module")
            run_genion_step(config, sam_files)
            
        if 'jaffal' in modules_to_run:
            logger.info("Running JaffaL module")
            run_jaffal_step(config, longgf_excel_file)

        # Integration and exon repair (if enabled)
        if config['options'].get('enable_integration', True):
            integration_method = config['options'].get('overlap_analysis_method', 'exon_repair')
            
            if integration_method == 'exon_repair' and config['options'].get('exon_repair', {}).get('enabled', True):
                logger.info("=" * 60)
                logger.info("Starting Exon Repair Protocol")
                logger.info("=" * 60)
                
                try:
                    exon_repair_results = run_exon_repair(
                        config=config,
                        output_dir=config['project']['output_dir']
                    )
                    
                    logger.info("Exon repair completed successfully")
                    logger.info(f"High-confidence chimeras: {exon_repair_results['statistics']['total_chimeras']}")
                    
                except Exception as e:
                    logger.error(f"Exon repair failed: {e}")
                    if config['options'].get('debug', False):
                        import traceback
                        logger.error(traceback.format_exc())
                    # Don't exit - allow pipeline to complete with warning
                    logger.warning("Continuing pipeline without exon repair")

        # Pipeline completion
        logger.info("=" * 60)
        logger.info("TYPHON Pipeline completed successfully!")
        logger.info(f"Results saved to: {config['project']['output_dir']}")
        logger.info(f"Log file: {log_file}")
        logger.info("=" * 60)
        
        return 0
        
    except Exception as e:
        # Use basic logging if logger not initialized yet
        try:
            logger.error(f"Pipeline failed: {e}")
            if config.get('options', {}).get('debug', False):
                import traceback
                logger.error(traceback.format_exc())
        except NameError:
            logging.error(f"Pipeline failed: {e}")
        return 1

if __name__ == '__main__':
    sys.exit(main()) 