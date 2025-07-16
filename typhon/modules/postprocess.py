import os
import glob
import re
import subprocess
from pathlib import Path


def process_log_files(log_dir):
    """
    Processes LongGF .log files to generate _results.txt files,
    mimicking the provided Bash script's logic.
    """
    log_files = glob.glob(os.path.join(log_dir, "*.log"))

    for log_file_path in log_files:
        filename_base = os.path.basename(log_file_path).rsplit('.', 1)[0]
        results_file_path = os.path.join(log_dir, f"{filename_base}_results.txt")
        intermediate_lines = []
        
        # Step 1 & 2: Extract lines between "GF" and "SumGF", excluding lines containing "SumGF"
        in_gf_block = False
        with open(log_file_path, "r") as f:
            for line in f:
                # DO NOT strip here, preserve original whitespace for later processing
                if "SumGF" in line:
                    in_gf_block = False
                    continue
                if line.startswith("GF"):
                    in_gf_block = True
                
                if in_gf_block:
                    intermediate_lines.append(line) # Append raw line

        # Step 3: Awk logic - join lines not starting with "GF"
        # This will now concatenate raw lines with their original internal spaces/tabs.
        processed_awk_lines = []
        current_record_parts = []
        for line in intermediate_lines:
            if line.startswith("GF"):
                if current_record_parts:
                    # Join previous parts; remove trailing newline from last part if present
                    processed_awk_lines.append("".join(current_record_parts).rstrip('\n'))
                current_record_parts = [line.rstrip('\n') + " "] # Append raw GF line, remove its newline, add one space
            else:
                current_record_parts.append(line.rstrip('\n')) # Append raw subsequent line, remove its newline
        if current_record_parts:
            processed_awk_lines.append("".join(current_record_parts).rstrip('\n'))

        # Step 4: Remove first 3 characters
        processed_awk_lines = [line[3:] for line in processed_awk_lines]

        # Step 5: Replace spaces with commas
        processed_awk_lines = [line.replace(' ', ',') for line in processed_awk_lines]

        # Step 6: Replace tabs with commas
        processed_awk_lines = [line.replace('\t', ',') for line in processed_awk_lines]

        # Step 7: Extract specific fields (1,2,13,14-)
        final_lines = []
        for line in processed_awk_lines:
            fields = line.split(',')
            if len(fields) >= 14:
                selected_fields = [fields[0], fields[1]] + fields[12:]  # 1,2,13,14- (0-indexed: 0,1,12,13...)
                final_lines.append(','.join(selected_fields))

        # Step 8: Remove trailing comma from last line if present
        if final_lines and final_lines[-1].endswith(','):
            final_lines[-1] = final_lines[-1][:-1]

        # Write the results
        with open(results_file_path, "w") as results_file:
            for line in final_lines:
                results_file.write(line + '\n')


def postprocess(output_dir):
    """
    Post-process LongGF results:
    1. Process log files to create _results.txt files
    2. Call R script to generate Excel files
    """
    # Process log files first
    process_log_files(output_dir)
    
    # Find the R script in the modules directory
    module_dir = Path(__file__).parent
    r_script_path = module_dir / "LongGF_process_results.R"
    
    if not r_script_path.exists():
        raise FileNotFoundError(f"LongGF_process_results.R not found at {r_script_path}")
    
    # Call R script to generate Excel files
    try:
        result = subprocess.run(
            ["Rscript", str(r_script_path), output_dir],
            capture_output=True,
            text=True,
            check=True
        )
        print("R script executed successfully")
        if result.stdout:
            print(f"R script output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error running R script: {e}")
        if e.stderr:
            print(f"R script error: {e.stderr}")
        raise
    except FileNotFoundError:
        print("Error: Rscript not found. Make sure R is installed and in PATH.")
        raise 