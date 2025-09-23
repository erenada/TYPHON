import subprocess
import logging
import time
from datetime import datetime


def setup_module_logger(module_name, output_dir):
    """
    Set up a module-specific logger that writes to both the main log and a module-specific log.
    
    Args:
        module_name: Name of the module (e.g., 'longgf', 'genion', 'jaffal')
        output_dir: Base output directory for logs
    
    Returns:
        logger: Module-specific logger instance
    """
    import os
    
    # Create module-specific logger
    logger = logging.getLogger(f"typhon.{module_name}")
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Create module-specific log file
    module_log_dir = os.path.join(output_dir, 'logs')
    os.makedirs(module_log_dir, exist_ok=True)
    module_log_file = os.path.join(module_log_dir, f'{module_name}.log')
    
    # Create file handler for module-specific log
    file_handler = logging.FileHandler(module_log_file)
    file_handler.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add handler to logger
    logger.addHandler(file_handler)
    
    # Also add to root logger so it appears in main typhon.log
    root_logger = logging.getLogger()
    if root_logger.handlers:
        # Copy the root logger's handlers to ensure messages go to main log
        for handler in root_logger.handlers:
            if handler not in logger.handlers:
                logger.addHandler(handler)
    
    logger.info(f"Starting {module_name} module")
    
    return logger


def run_command(cmd, shell=True, check=True, capture_output=False, cwd=None, log_output=True):
    """
    Run a shell command with enhanced logging and error checking.
    
    Args:
        cmd: Command string or list
        shell: Whether to use shell
        check: Whether to check return code
        capture_output: Whether to capture output for return
        cwd: Working directory
        log_output: Whether to log command output (default: True)
    
    Returns:
        subprocess.CompletedProcess or stdout string if capture_output=True
    """
    if isinstance(cmd, list):
        cmd_str = ' '.join(str(x) for x in cmd)
    else:
        cmd_str = str(cmd)
    
    # Log command start with timestamp
    start_time = time.time()
    logging.info(f"[CMD START] {cmd_str}")
    
    try:
        if isinstance(cmd, list) and shell:
            # If cmd is a list but shell=True, join it
            cmd = cmd_str
        elif isinstance(cmd, str) and not shell:
            # If cmd is a string but shell=False, split it
            cmd = cmd.split()
        
        # Always capture output for logging, even if not requested by caller
        result = subprocess.run(
            cmd, 
            shell=shell, 
            check=check, 
            capture_output=True,  # Always capture for logging
            text=True, 
            cwd=cwd
        )
        
        # Calculate execution time
        end_time = time.time()
        duration = end_time - start_time
        
        # Log command completion and timing
        logging.info(f"[CMD DONE] Completed in {duration:.2f}s (exit code: {result.returncode})")
        
        # Log stdout if present and logging enabled
        if log_output and result.stdout.strip():
            # Split long output for better readability
            stdout_lines = result.stdout.strip().split('\n')
            logging.info(f"[CMD STDOUT] {len(stdout_lines)} lines of output:")
            for i, line in enumerate(stdout_lines):
                # Log first 50 and last 10 lines for very long output
                if len(stdout_lines) > 60 and 50 <= i < len(stdout_lines) - 10:
                    if i == 50:
                        logging.info(f"[CMD STDOUT] ... ({len(stdout_lines) - 60} lines omitted) ...")
                    continue
                logging.info(f"[CMD STDOUT] {line}")
        
        # Log stderr if present (always log errors)
        if result.stderr.strip():
            stderr_lines = result.stderr.strip().split('\n')
            logging.warning(f"[CMD STDERR] {len(stderr_lines)} lines of stderr:")
            for line in stderr_lines:
                logging.warning(f"[CMD STDERR] {line}")
        
        # Return based on caller's request
        if capture_output:
            return result.stdout
        return result
        
    except subprocess.CalledProcessError as e:
        end_time = time.time()
        duration = end_time - start_time
        
        # Enhanced error logging
        logging.error(f"[CMD FAILED] Command failed after {duration:.2f}s: {cmd_str}")
        logging.error(f"[CMD FAILED] Exit code: {e.returncode}")
        
        if e.stdout:
            logging.error(f"[CMD FAILED] Last stdout: {e.stdout.strip()}")
        if e.stderr:
            logging.error(f"[CMD FAILED] Error output: {e.stderr.strip()}")
            
        raise
    except Exception as e:
        end_time = time.time()
        duration = end_time - start_time
        logging.error(f"[CMD ERROR] Unexpected error after {duration:.2f}s: {cmd_str}")
        logging.error(f"[CMD ERROR] Exception: {e}")
        raise


def run_command_with_realtime_output(cmd, shell=True, check=True, cwd=None):
    """
    Run a command with real-time output streaming to both console and log.
    Useful for long-running commands like minimap2, samtools, etc.
    """
    if isinstance(cmd, list):
        cmd_str = ' '.join(str(x) for x in cmd)
    else:
        cmd_str = str(cmd)
    
    logging.info(f"[CMD REALTIME] Starting: {cmd_str}")
    start_time = time.time()
    
    try:
        if isinstance(cmd, list) and shell:
            cmd = cmd_str
        elif isinstance(cmd, str) and not shell:
            cmd = cmd.split()
        
        # Start process with pipes
        process = subprocess.Popen(
            cmd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Merge stderr with stdout
            text=True,
            cwd=cwd,
            bufsize=1,  # Line buffered
            universal_newlines=True
        )
        
        # Read and log output in real-time
        output_lines = []
        if process.stdout:
            for line in process.stdout:
                line = line.rstrip()
                if line:  # Only log non-empty lines
                    print(line)  # Console output
                    logging.info(f"[CMD OUTPUT] {line}")  # Log output
                    output_lines.append(line)
        
        # Wait for process to complete
        return_code = process.wait()
        
        end_time = time.time()
        duration = end_time - start_time
        
        if return_code == 0:
            logging.info(f"[CMD REALTIME] Completed successfully in {duration:.2f}s")
        else:
            logging.error(f"[CMD REALTIME] Failed with exit code {return_code} after {duration:.2f}s")
            if check:
                raise subprocess.CalledProcessError(return_code, cmd_str)
        
        return return_code
        
    except Exception as e:
        end_time = time.time()
        duration = end_time - start_time
        logging.error(f"[CMD REALTIME] Error after {duration:.2f}s: {e}")
        raise 