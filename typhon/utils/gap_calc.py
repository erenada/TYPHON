#!/usr/bin/env python3
"""
Gap Calculation Utility

Calculates maximum internal gap lengths in aligned sequences for QC analysis.
Based on original script from TYPHON wrapper script.

Authors: Harry Kane, PhD; Eren Ada, PhD
Original credits preserved below:

Credit to user "Laurent H." for their great public script on stackoverflow.com
See: https://stackoverflow.com/questions/49303791/finding-the-positions-and-lengths-of-gaps-indels-in-a-sequence-alignment-with

Credit to user "pippo1980" for their great public script on bioinformatics.stackexchange.com  
See: https://bioinformatics.stackexchange.com/questions/20442/find-open-reading-frames-in-a-dna-sequence
"""

from Bio import SeqIO
import re
import sys
import argparse
import logging


def calculate_max_gap_length(sequence_str):
    """
    Calculate the maximum internal gap length in a sequence.
    
    Args:
        sequence_str (str): DNA/RNA sequence string
        
    Returns:
        int: Maximum gap length found, or 0 if no gaps
    """
    # Pattern: letter(s) + gap(s) + letter(s)
    pattern_A = '[A-Z]'      # Letters before gap
    pattern_B = '[-]+'       # Gap(s)
    pattern_C = '[A-Z]'      # Letters after gap
    
    full_pattern = pattern_A + pattern_B + pattern_C
    
    # Find all matches and get the longest gap
    matches = re.findall(full_pattern, str(sequence_str))
    if matches:
        max_gap = max(matches, key=len, default="22")  # Default "22" for empty case
        return len(max_gap) - 2  # Subtract 2 for the flanking letters
    else:
        return 0


def process_fasta_file(fasta_file):
    """
    Process FASTA file and calculate max gap lengths for all sequences.
    
    Args:
        fasta_file (str): Path to FASTA file
        
    Returns:
        list: Maximum gap lengths for each sequence
    """
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        gap_lengths = []
        
        for record in records:
            gap_length = calculate_max_gap_length(record.seq)
            gap_lengths.append(gap_length)
            
        return gap_lengths
        
    except Exception as e:
        logging.error(f"Error processing FASTA file {fasta_file}: {e}")
        return []


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(
        description="Calculate maximum internal gap lengths in aligned sequences"
    )
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')
    
    # Process FASTA file
    gap_lengths = process_fasta_file(args.fasta_file)
    
    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            for gap_length in gap_lengths:
                f.write(f"{gap_length}\n")
    else:
        # Print to stdout (original behavior)
        for gap_length in gap_lengths:
            print(gap_length)


if __name__ == "__main__":
    main() 