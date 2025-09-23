#!/usr/bin/env python3
"""
Automated setup script for custom Genion build for TYPHON Pipeline.

This script implements the complete TYPHON Genion customization process:
- Clones official Genion repo
- Replaces annotate.cpp with TYPHON custom version (read ID tracking)
- Applies additional cleanup patch
- Compiles with configurable debug flags
- Installs binary to configurable location

The TYPHON modifications enable:
1. Read-level output (one line per supporting read ID)
2. Full debug output enabled by default
3. Enhanced tracking for fusion validation

Authors: Harry Kane, PhD; Eren Ada, PhD
Based on AI_GENION_CUSTOMIZATION_GUIDE.md and detailed analysis
"""
import os
import sys
import subprocess
import shutil
import tempfile
import yaml
from datetime import datetime
import argparse


def load_config(config_path="config.yaml"):
    """Load configuration from config.yaml."""
    if not os.path.exists(config_path):
        return {}
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def setup_logging(log_path, debug=False):
    """Setup logging with configurable path."""
    import logging
    
    level = logging.DEBUG if debug else logging.INFO
    
    # Create log directory if it doesn't exist
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger()


def log_to_file(msg, log_path):
    """Log message to file."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(log_path, 'a') as f:
        f.write(f"[{timestamp}] {msg}\n")
    print(f"[{timestamp}] {msg}")


def run_cmd(cmd, cwd=None, check=True, log_path=None):
    """Run command with logging."""
    cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
    if log_path:
        log_to_file(f"Running: {cmd_str} (in {cwd or os.getcwd()})", log_path)
    
    try:
        result = subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)
        if log_path and result.stdout:
            log_to_file(result.stdout, log_path)
        if log_path and result.stderr:
            log_to_file(f"STDERR: {result.stderr}", log_path)
        return result
    except subprocess.CalledProcessError as e:
        if log_path:
            log_to_file(f"ERROR: Command failed: {cmd_str}", log_path)
            log_to_file(f"STDOUT: {e.stdout}", log_path)
            log_to_file(f"STDERR: {e.stderr}", log_path)
        raise


def verify_typhon_files(script_dir, log_path):
    """Verify that required TYPHON customization files exist."""
    required_files = {
        'custom_annotate': os.path.join(script_dir, 'Genion_files', 'annotate.cpp'),
        'patch_file': os.path.join(script_dir, 'Genion_files', 'genion_custom.patch')
    }
    
    missing_files = []
    for name, path in required_files.items():
        if not os.path.exists(path):
            missing_files.append(f"{name}: {path}")
    
    if missing_files:
        error_msg = f"Missing required TYPHON customization files:\n" + "\n".join(missing_files)
        log_to_file(f"ERROR: {error_msg}", log_path)
        raise FileNotFoundError(error_msg)
    
    log_to_file("All required TYPHON customization files found", log_path)
    return required_files


def apply_typhon_modifications(temp_dir, custom_annotate_path, patch_path, log_path):
    """
    Apply TYPHON modifications to Genion source code.
    
    This implements the complete customization process:
    1. Backup original annotate.cpp
    2. Replace with TYPHON version (includes read ID tracking)
    3. Apply additional cleanup patch
    """
    src_dir = os.path.join(temp_dir, "src")
    original_annotate = os.path.join(src_dir, "annotate.cpp")
    backup_annotate = os.path.join(src_dir, "annotate.cpp.original")
    
    # Step 1: Backup original file
    if os.path.exists(original_annotate):
        shutil.copy2(original_annotate, backup_annotate)
        log_to_file("Backed up original annotate.cpp", log_path)
    
    # Step 2: Replace with TYPHON custom version
    shutil.copy2(custom_annotate_path, original_annotate)
    log_to_file("Replaced annotate.cpp with TYPHON custom version", log_path)
    log_to_file("  - Enables full debug output by default", log_path)
    log_to_file("  - Adds read ID tracking to fusion output", log_path)
    log_to_file("  - Outputs one line per supporting read", log_path)
    
    # Step 3: Apply additional cleanup patch
    try:
        with open(patch_path, 'r') as patch_file:
            result = subprocess.run(
                ["patch", "-p0"], 
                input=patch_file.read(),
                text=True,
                cwd=temp_dir, 
                capture_output=True
            )
            if result.returncode == 0:
                log_to_file("Applied additional cleanup patch", log_path)
            else:
                log_to_file(f"WARNING: Patch application had issues: {result.stderr}", log_path)
    except subprocess.CalledProcessError as e:
        log_to_file(f"WARNING: Patch application failed (may be already applied): {e}", log_path)
        # Continue - the main customization (file replacement) is already done


def cleanup_files(patch_path, custom_annotate_path, log_path):
    """Remove customization files if requested."""
    files_to_clean = [patch_path, custom_annotate_path]
    for f in files_to_clean:
        if os.path.exists(f):
            try:
                os.remove(f)
                log_to_file(f"Removed file: {f}", log_path)
            except Exception as e:
                log_to_file(f"WARNING: Failed to remove {f}: {e}", log_path)


def main():
    parser = argparse.ArgumentParser(
        description="Setup custom TYPHON Genion build with read ID tracking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
TYPHON Genion Customizations:
  - Enables full debug output by default
  - Adds read ID tracking for each fusion candidate
  - Outputs one line per supporting read (not per fusion)
  - Enables detailed downstream analysis and validation

Examples:
  # Basic setup with config.yaml (recommended)
  python setup_genion.py --debug
  
  # Setup with custom configuration
  python setup_genion.py --output-dir ./custom_bin --debug-compilation --threads 8
        """
    )
    parser.add_argument('--config', default='config.yaml',
                        help='YAML configuration file (default: config.yaml)')
    parser.add_argument('--cleanup', action='store_true', 
                        help='Remove patch and source files after build')
    parser.add_argument('--debug', action='store_true', 
                        help='Enable debug logging')
    parser.add_argument('--output-dir', 
                        help='Output directory for binary (overrides config.yaml)')
    parser.add_argument('--debug-compilation', action='store_true',
                        help='Compile with debug flags -g -O0 (overrides config.yaml)')
    parser.add_argument('--threads', type=int,
                        help='Number of threads for compilation (overrides config.yaml)')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Get configuration sections
    genion_config = config.get('modules', {}).get('genion', {})
    project_config = config.get('project', {})
    
    # Determine configuration values (command line overrides config.yaml)
    output_bin_dir = args.output_dir or genion_config.get('output_bin_dir') or os.path.join(os.path.dirname(__file__), 'bin')
    debug_compilation = args.debug_compilation or genion_config.get('debug_compilation', True)
    threads = args.threads or genion_config.get('threads') or project_config.get('threads', 1)
    
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Setup logging
    log_dir = os.path.join(project_config.get('output_dir', script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)  # Create log directory if it doesn't exist
    log_path = os.path.join(log_dir, 'setup_genion.log')
    
    # Initialize logging
    if args.debug:
        logger = setup_logging(log_path, debug=True)
    
    # Repository URL
    genion_repo = "https://github.com/vpc-ccg/genion"
    
    log_to_file("=" * 70, log_path)
    log_to_file("TYPHON Custom Genion Setup", log_path)
    log_to_file(f"Timestamp: {datetime.now()}", log_path)
    log_to_file("=" * 70, log_path)
    log_to_file("Configuration:", log_path)
    log_to_file(f"  Output directory: {output_bin_dir}", log_path)
    log_to_file(f"  Debug compilation: {debug_compilation}", log_path)
    log_to_file(f"  Threads: {threads}", log_path)
    log_to_file("", log_path)
    log_to_file("TYPHON Customizations:", log_path)
    log_to_file("  - Read ID tracking per fusion candidate", log_path)
    log_to_file("  - Full debug output enabled by default", log_path)
    log_to_file("  - One line per supporting read (not per fusion)", log_path)
    
    temp_dir = tempfile.mkdtemp(prefix="typhon_genion_build_")
    try:
        # Verify TYPHON customization files exist
        required_files = verify_typhon_files(script_dir, log_path)
        
        # Step 1: Clone Genion repo
        log_to_file("", log_path)
        log_to_file("Step 1: Cloning official Genion repository...", log_path)
        run_cmd(["git", "clone", genion_repo, temp_dir], log_path=log_path)

        # Step 2: Apply TYPHON modifications
        log_to_file("", log_path)
        log_to_file("Step 2: Applying TYPHON customizations...", log_path)
        apply_typhon_modifications(
            temp_dir, 
            required_files['custom_annotate'],
            required_files['patch_file'],
            log_path
        )

        # Step 3: Compile with optional debug flags
        log_to_file("", log_path)
        log_to_file("Step 3: Compiling custom Genion...", log_path)
        
        make_cmd = ["make"]
        if threads > 1:
            make_cmd.extend(["-j", str(threads)])
        
        # Set compilation flags if debug mode requested
        env = os.environ.copy()
        if debug_compilation:
            env["CXXFLAGS"] = "-g -O0 -DDEBUG -Wall"
            env["CFLAGS"] = "-g -O0 -DDEBUG -Wall"
            log_to_file(f"Using debug compilation flags: {env['CXXFLAGS']}", log_path)
        else:
            log_to_file("Using default optimized compilation flags", log_path)
            
        result = subprocess.run(
            make_cmd, 
            cwd=temp_dir, 
            env=env,
            capture_output=True, 
            text=True
        )
        if result.returncode != 0:
            log_to_file(f"ERROR: Compilation failed: {result.stderr}", log_path)
            raise subprocess.CalledProcessError(result.returncode, make_cmd)
        log_to_file("Compilation completed successfully", log_path)

        # Step 4: Install binary
        log_to_file("", log_path)
        log_to_file("Step 4: Installing custom Genion binary...", log_path)
        
        bin_path = os.path.join(temp_dir, "genion")
        if not os.path.exists(bin_path):
            raise FileNotFoundError(f"Compiled genion binary not found at {bin_path}")
        
        os.makedirs(output_bin_dir, exist_ok=True)
        dest_bin = os.path.join(output_bin_dir, "genion")
        shutil.move(bin_path, dest_bin)
        os.chmod(dest_bin, 0o755)
        log_to_file(f"Installed to: {dest_bin}", log_path)
        
        # Step 5: Verify installation
        try:
            result = run_cmd([dest_bin, "--help"], log_path=log_path, check=False)
            if result.returncode == 0 or "usage" in result.stderr.lower():
                log_to_file("Binary installation verified", log_path)
            else:
                log_to_file("WARNING: Binary verification inconclusive", log_path)
        except Exception:
            log_to_file("WARNING: Could not verify binary (may still be functional)", log_path)

        log_to_file("=" * 70, log_path)
        log_to_file("TYPHON Custom Genion Setup Completed Successfully!", log_path)
        log_to_file("", log_path)
        log_to_file("Genion setup completed successfully", log_path)
        log_to_file(f"Binary location: {dest_bin}", log_path)
        log_to_file(f"Compilation mode: {'Debug' if debug_compilation else 'Optimized'}", log_path)
        log_to_file("=" * 70, log_path)
        
        # Optional cleanup
        if args.cleanup:
            log_to_file("", log_path)
            log_to_file("Performing cleanup of source files...", log_path)
            cleanup_files(
                required_files['patch_file'], 
                required_files['custom_annotate'], 
                log_path
            )
            
    except Exception as e:
        log_to_file("", log_path)
        log_to_file("=" * 70, log_path)
        log_to_file(f"SETUP FAILED: {e}", log_path)
        if args.debug:
            import traceback
            log_to_file("Full traceback:", log_path)
            log_to_file(traceback.format_exc(), log_path)
        log_to_file("=" * 70, log_path)
        sys.exit(1)
    finally:
        # Clean up temp directory
        try:
            shutil.rmtree(temp_dir)
            log_to_file(f"Cleaned up temporary directory: {temp_dir}", log_path)
        except Exception as cleanup_err:
            log_to_file(f"WARNING: Failed to clean up temp dir: {cleanup_err}", log_path)


if __name__ == "__main__":
    main() 