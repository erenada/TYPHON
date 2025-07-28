#!/usr/bin/env python3
"""
Sequence Utilities for TYPHON Exon Repair Protocol

Contains utility functions for FASTA sequence manipulation, particularly
for merging sequences with identical headers as required by the exon repair protocol.

Author: Eren Ada, PhD

Credit: merge_sequences_by_header function based on public script by user "phngs" 
from bioinformatics.stackexchange.com - thank you for the elegant solution!
"""

import os
import logging
from typing import Dict, List


def merge_sequences_by_header(input_fasta: str, output_fasta: str) -> None:
    """
    Merge FASTA sequences with identical headers into single entries.
    
    This function reads a FASTA file and concatenates all sequences that share
    the same header name. This is essential for the exon repair protocol where
    multiple exon segments per Read_ID need to be merged into single sequences.
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output FASTA file with merged sequences
        
    Raises:
        FileNotFoundError: If input FASTA file doesn't exist
        IOError: If there are issues reading/writing files
    """
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    
    logger.debug(f"Merging sequences in {input_fasta} -> {output_fasta}")
    
    # Dictionary to store merged sequences by header
    data = {}
    
    try:
        with open(input_fasta, "r") as fasta:
            current_header = ""
            for line in fasta:
                if line.startswith(">"):
                    current_header = line.rstrip()  # Remove newline but keep header
                    if current_header not in data:
                        data[current_header] = ""
                else:
                    # Strip newlines and concatenate sequence
                    data[current_header] += line.strip()
        
        # Write merged sequences to output file
        with open(output_fasta, "w") as out:
            for fasta_header, fasta_sequence in data.items():
                out.write(f"{fasta_header}\n{fasta_sequence}\n")
                
        logger.debug(f"Successfully merged {len(data)} unique sequences")
        
    except IOError as e:
        logger.error(f"Error processing FASTA files: {e}")
        raise


def validate_fasta_file(fasta_path: str) -> bool:
    """
    Validate that a FASTA file exists and has valid format.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        bool: True if valid, False otherwise
    """
    if not os.path.exists(fasta_path):
        return False
        
    try:
        with open(fasta_path, 'r') as f:
            first_line = f.readline().strip()
            return first_line.startswith('>')
    except Exception:
        return False


def count_sequences_in_fasta(fasta_path: str) -> int:
    """
    Count the number of sequences in a FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        int: Number of sequences (headers starting with >)
    """
    if not os.path.exists(fasta_path):
        return 0
        
    count = 0
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        return 0
        
    return count


def get_sequence_headers(fasta_path: str) -> List[str]:
    """
    Extract all sequence headers from a FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        List[str]: List of headers (without > symbol)
    """
    headers = []
    
    if not os.path.exists(fasta_path):
        return headers
        
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    headers.append(line[1:].strip())  # Remove > and whitespace
    except Exception:
        pass
        
    return headers 