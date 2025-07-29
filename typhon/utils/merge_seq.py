#!/usr/bin/env python3
"""
Sequence Merging Utility

Merges FASTA sequences with the same header/name by concatenating their sequences.
Based on original script from TYPHON wrapper script.

Authors: Harry Kane, PhD; Eren Ada, PhD
Original credits preserved below:

Credit to user "phngs" for their great public script on bioinformatics.stackexchange.com
See: https://bioinformatics.stackexchange.com/questions/4714/how-to-merge-transcript-sequence-with-same-name-in-a-fasta-file
"""

import sys
import argparse
import logging
from collections import defaultdict


def merge_sequences_from_file(fasta_file, output_file=None):
    """
    Merge sequences with the same header in a FASTA file.
    
    Args:
        fasta_file (str): Path to input FASTA file
        output_file (str, optional): Path to output file. If None, prints to stdout.
        
    Returns:
        dict: Dictionary of merged sequences {header: sequence}
    """
    data = defaultdict(str)
    current_header = ""
    
    try:
        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_header = line
                elif current_header:  # Only process sequence lines if we have a header
                    data[current_header] += line
                    
        # Output results
        if output_file:
            with open(output_file, 'w') as out_f:
                for header, sequence in data.items():
                    out_f.write(f"{header}\n{sequence}\n")
        else:
            # Print to stdout (original behavior)
            for header, sequence in data.items():
                print(f"{header}{sequence}")
                
        return dict(data)
        
    except Exception as e:
        logging.error(f"Error processing FASTA file {fasta_file}: {e}")
        return {}


def merge_sequences_from_data(fasta_data):
    """
    Merge sequences with the same header from FASTA data in memory.
    
    Args:
        fasta_data (str): FASTA format string
        
    Returns:
        dict: Dictionary of merged sequences {header: sequence}
    """
    data = defaultdict(str)
    current_header = ""
    
    for line in fasta_data.split('\n'):
        line = line.strip()
        if line.startswith(">"):
            current_header = line
        elif current_header and line:  # Skip empty lines
            data[current_header] += line
            
    return dict(data)


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description="Merge FASTA sequences with the same header"
    )
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')
    
    # Process file
    merged_seqs = merge_sequences_from_file(args.fasta_file, args.output)
    
    if args.verbose:
        logging.info(f"Merged {len(merged_seqs)} unique sequences")


if __name__ == "__main__":
    main() 