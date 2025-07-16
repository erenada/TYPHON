import subprocess
import logging

def run_command(cmd, shell=True, check=True, capture_output=False, cwd=None):
    """Run a shell command with error checking and optional output capture."""
    if isinstance(cmd, list):
        cmd_str = ' '.join(str(x) for x in cmd)
    else:
        cmd_str = str(cmd)
    
    logging.info(f"Running command: {cmd_str}")
    try:
        if isinstance(cmd, list) and shell:
            # If cmd is a list but shell=True, join it
            cmd = cmd_str
        elif isinstance(cmd, str) and not shell:
            # If cmd is a string but shell=False, split it
            cmd = cmd.split()
            
        result = subprocess.run(cmd, shell=shell, check=check, capture_output=capture_output, text=True, cwd=cwd)
        if capture_output:
            return result.stdout
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd_str}\nError: {e}")
        raise 