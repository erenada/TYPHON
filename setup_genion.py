#!/usr/bin/env python3
"""
Automated setup script for custom Genion build for Typhon.
- Clones official Genion repo
- Applies custom patch to annotate.cpp
- Compiles Genion
- Moves binary to Typhon/bin
- Cleans up temp files (optional)
- Logs all steps
"""
import os
import sys
import subprocess
import shutil
import tempfile
from datetime import datetime
import argparse

# Paths (edit as needed)
GENION_REPO = "https://github.com/vpc-ccg/genion"
PATCH_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Genion_files/genion_custom.patch'))
EXTRA_ANNOTATE = os.path.abspath(os.path.join(os.path.dirname(__file__), 'Genion_files/annotate.cpp'))
OUTPUT_BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'bin'))
LOG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), 'setup_genion.log'))


def log(msg):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(LOG_PATH, 'a') as f:
        f.write(f"[{timestamp}] {msg}\n")
    print(f"[{timestamp}] {msg}")


def run_cmd(cmd, cwd=None, check=True):
    log(f"Running: {' '.join(cmd)} (in {cwd or os.getcwd()})")
    try:
        result = subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)
        log(result.stdout)
        if result.stderr:
            log(f"STDERR: {result.stderr}")
        return result
    except subprocess.CalledProcessError as e:
        log(f"ERROR: Command failed: {' '.join(cmd)}")
        log(f"STDOUT: {e.stdout}")
        log(f"STDERR: {e.stderr}")
        raise


def cleanup_files():
    # Remove patch file and extra annotate.cpp if they exist
    for f in [PATCH_PATH, EXTRA_ANNOTATE]:
        if os.path.exists(f):
            try:
                os.remove(f)
                log(f"Removed file: {f}")
            except Exception as e:
                log(f"WARNING: Failed to remove {f}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Setup custom Genion build for Typhon.")
    parser.add_argument('--cleanup', action='store_true', help='Remove patch and extra files after build')
    args = parser.parse_args()

    log("--- Genion setup started ---")
    temp_dir = tempfile.mkdtemp(prefix="genion_build_")
    try:
        # 1. Clone Genion repo
        log(f"Cloning Genion repo to {temp_dir}")
        run_cmd(["git", "clone", GENION_REPO, temp_dir])

        # 2. Apply patch
        src_annotate = os.path.join(temp_dir, "src", "annotate.cpp")
        log(f"Applying patch {PATCH_PATH} to {src_annotate}")
        run_cmd(["patch", src_annotate, PATCH_PATH], cwd=os.path.join(temp_dir, "src"))

        # 3. Compile Genion
        log("Compiling Genion with make")
        run_cmd(["make"], cwd=temp_dir)

        # 4. Move binary to OUTPUT_BIN_DIR
        bin_path = os.path.join(temp_dir, "genion")
        if not os.path.exists(bin_path):
            log("ERROR: Compiled genion binary not found!")
            raise FileNotFoundError("Compiled genion binary not found!")
        os.makedirs(OUTPUT_BIN_DIR, exist_ok=True)
        dest_bin = os.path.join(OUTPUT_BIN_DIR, "genion")
        shutil.move(bin_path, dest_bin)
        log(f"Moved genion binary to {dest_bin}")

        log("Genion setup completed successfully.")
        # Optional cleanup
        if args.cleanup:
            log("Performing cleanup of patch and extra files.")
            cleanup_files()
    except Exception as e:
        log(f"SETUP FAILED: {e}")
        sys.exit(1)
    finally:
        # 5. Clean up temp dir
        try:
            shutil.rmtree(temp_dir)
            log(f"Cleaned up temp directory {temp_dir}")
        except Exception as cleanup_err:
            log(f"WARNING: Failed to clean up temp dir: {cleanup_err}")

if __name__ == "__main__":
    main() 