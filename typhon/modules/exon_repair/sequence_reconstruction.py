#!/usr/bin/env python3
"""
Phase 5: Sequence Reconstruction for TYPHON Exon Repair Protocol

This module implements the sequence reconstruction pipeline from the R workflow,
including bedtools getfasta extraction, multi-step sequence merging, and final
chimeric sequence assembly.


Authors: Harry Kane, PhD; Eren Ada, PhD
"""

import os
import sys
import logging
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# Add typhon utils to path for sequence utilities
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from utils.sequence_utils import merge_sequences_by_header, validate_fasta_file, count_sequences_in_fasta


class SequenceReconstructor:
    """Handles Phase 5: Sequence Reconstruction for exon repair protocol."""
    
    def __init__(self, config: dict, output_dir: str):
        """
        Initialize sequence reconstructor.
        
        Args:
            config: TYPHON configuration dictionary
            output_dir: Main pipeline output directory
        """
        self.config = config
        self.output_dir = output_dir
        self.logger = logging.getLogger(__name__)
        
        # Create phase 5 working directory
        self.work_dir = os.path.join(output_dir, 'exon_repair')
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Reference files from config
        self.genome_fasta = config['references']['genome']
        
        # Processing parameters
        self.threads = config['project'].get('threads', 4)
        self.keep_intermediate = config['options'].get('exon_repair', {}).get('keep_intermediate', True)


    def reconstruct_sequences(self, phase4_results: dict, chimera_library: pd.DataFrame) -> dict:
        """
        Main sequence reconstruction pipeline implementing R workflow lines 311-337.
        
        Args:
            phase4_results: Output from Phase 4 with BED files and summary data
            chimera_library: Original chimera library for final filtering
            
        Returns:
            dict: Complete results including reconstructed sequences and filtered library
        """
        self.logger.info("=" * 50)
        self.logger.info("Phase 5: Sequence Reconstruction")
        self.logger.info("=" * 50)
        
        try:
            # 5.1: Extract sequences using bedtools getfasta
            geneA_fasta, geneB_fasta = self._extract_sequences_with_bedtools(
                phase4_results['bed_file_A'], 
                phase4_results['bed_file_B']
            )
            
            # 5.2: Clean strand indicators from FASTA headers
            self._clean_strand_indicators(geneA_fasta, geneB_fasta)
            
            # 5.3: Multi-step sequence merging process
            final_fasta = self._merge_sequences_multistep(geneA_fasta, geneB_fasta)
            
            # 5.4: Final filtering and status assignment
            filtered_results = self._filter_and_assign_status(
                phase4_results['summary_data'], 
                chimera_library
            )
            
            # Generate results
            results = {
                'reconstructed_sequences': final_fasta,
                'filtered_chimeras': filtered_results['filtered_chimeras'],
                'final_chimeras_path': filtered_results['output_path'],
                'statistics': self._generate_statistics(filtered_results, final_fasta),
                'intermediate_files': {
                    'geneA_fasta': geneA_fasta,
                    'geneB_fasta': geneB_fasta,
                    'geneA_collapsed': os.path.join(self.work_dir, 'Fasta_geneA_collapse.fa'),
                    'geneB_collapsed': os.path.join(self.work_dir, 'Fasta_geneB_collapse.fa'),
                    'combined_fasta': os.path.join(self.work_dir, 'Fasta_cat_not_merged_seqs_exon_repair.fa')
                }
            }
            
            self.logger.info("Phase 5 completed successfully:")
            self.logger.info(f"  - Reconstructed sequences: {results['statistics']['total_reconstructed_sequences']}")
            self.logger.info(f"  - Filtered chimeras: {results['statistics']['total_filtered_chimeras']}")
            self.logger.info(f"  - Intrachromosomal: {results['statistics']['intrachromosomal_count']}")
            self.logger.info(f"  - Interchromosomal: {results['statistics']['interchromosomal_count']}")
            self.logger.info(f"  - Final sequences: {final_fasta}")
            
            return results
            
        except Exception as e:
            self.logger.error(f"Phase 5 sequence reconstruction failed: {str(e)}")
            raise


    def _extract_sequences_with_bedtools(self, bed_file_A: str, bed_file_B: str) -> Tuple[str, str]:
        """
        Extract sequences using bedtools getfasta with strand information.
        
        Implements R workflow lines 311-317:
        bedtools getfasta -s -fi genome.fa -bed bed_file_A.bed -fo geneA.fa -nameOnly
        bedtools getfasta -s -fi genome.fa -bed bed_file_B.bed -fo geneB.fa -nameOnly
        
        Args:
            bed_file_A: BED file for GeneA segments
            bed_file_B: BED file for GeneB segments
            
        Returns:
            Tuple of (geneA_fasta_path, geneB_fasta_path)
        """
        self.logger.info("Step 5.1: Extracting sequences with bedtools getfasta...")
        
        # Output file paths
        geneA_fasta = os.path.join(self.work_dir, 'Fasta_geneA_not_collapse.fa')
        geneB_fasta = os.path.join(self.work_dir, 'Fasta_geneB_not_collapse.fa')
        
        try:
            # Extract GeneA sequences
            cmd_geneA = [
                'bedtools', 'getfasta', '-s', '-fi', str(self.genome_fasta),
                '-bed', str(bed_file_A), '-fo', str(geneA_fasta), '-nameOnly'
            ]
            self.logger.info(f"Running: {' '.join(cmd_geneA)}")
            result_A = subprocess.run(cmd_geneA, check=True, capture_output=True, text=True)
            
            # Extract GeneB sequences  
            cmd_geneB = [
                'bedtools', 'getfasta', '-s', '-fi', str(self.genome_fasta),
                '-bed', str(bed_file_B), '-fo', str(geneB_fasta), '-nameOnly'
            ]
            self.logger.info(f"Running: {' '.join(cmd_geneB)}")
            result_B = subprocess.run(cmd_geneB, check=True, capture_output=True, text=True)
            
            # Validate output files
            if not validate_fasta_file(geneA_fasta):
                raise FileNotFoundError(f"GeneA FASTA extraction failed: {geneA_fasta}")
            if not validate_fasta_file(geneB_fasta):
                raise FileNotFoundError(f"GeneB FASTA extraction failed: {geneB_fasta}")
            
            # Log extraction statistics
            geneA_count = count_sequences_in_fasta(geneA_fasta)
            geneB_count = count_sequences_in_fasta(geneB_fasta)
            self.logger.info(f"Extracted {geneA_count} GeneA segments and {geneB_count} GeneB segments")
            
            return geneA_fasta, geneB_fasta
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"bedtools getfasta failed: {e}")
            self.logger.error(f"stderr: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Sequence extraction failed: {e}")
            raise


    def _clean_strand_indicators(self, geneA_fasta: str, geneB_fasta: str) -> None:
        """
        Clean strand indicators from FASTA headers using sed.
        
        Implements R workflow lines 314-315:
        sed -i 's/(+)//;s/(-)//' Fasta_geneA_not_collapse.fa
        sed -i 's/(+)//;s/(-)//' Fasta_geneB_not_collapse.fa
        
        Args:
            geneA_fasta: Path to GeneA FASTA file
            geneB_fasta: Path to GeneB FASTA file
        """
        self.logger.info("Step 5.2: Cleaning strand indicators from FASTA headers...")
        
        try:
            # Clean GeneA FASTA
            cmd_clean_A = ['sed', '-i', 's/(+)//;s/(-)//', str(geneA_fasta)]
            self.logger.debug(f"Running: {' '.join(cmd_clean_A)}")
            subprocess.run(cmd_clean_A, check=True)
            
            # Clean GeneB FASTA
            cmd_clean_B = ['sed', '-i', 's/(+)//;s/(-)//', str(geneB_fasta)]
            self.logger.debug(f"Running: {' '.join(cmd_clean_B)}")
            subprocess.run(cmd_clean_B, check=True)
            
            self.logger.info("Successfully cleaned strand indicators")
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to clean strand indicators: {e}")
            raise


    def _merge_sequences_multistep(self, geneA_fasta: str, geneB_fasta: str) -> str:
        """
        Multi-step sequence merging process.
        
        Implements R workflow lines 320-327:
        1. Merge GeneA sequences by Read_ID
        2. Merge GeneB sequences by Read_ID  
        3. Concatenate GeneA and GeneB files
        4. Final merge by Read_ID to create chimeric sequences
        5. Rename sequences for final output
        
        Args:
            geneA_fasta: Path to GeneA FASTA file
            geneB_fasta: Path to GeneB FASTA file
            
        Returns:
            str: Path to final reconstructed FASTA file
        """
        self.logger.info("Step 5.3: Multi-step sequence merging process...")
        
        # Define intermediate file paths
        geneA_collapsed = os.path.join(self.work_dir, 'Fasta_geneA_collapse.fa')
        geneB_collapsed = os.path.join(self.work_dir, 'Fasta_geneB_collapse.fa')
        combined_fasta = os.path.join(self.work_dir, 'Fasta_cat_not_merged_seqs_exon_repair.fa')
        merged_fasta = os.path.join(self.work_dir, 'Merged_seqs_exon_repair.fa')
        final_fasta = os.path.join(self.work_dir, 'Merged_seqs_exon_repair_renamed.fa')
        
        try:
            # Step 1: Merge GeneA sequences by Read_ID
            self.logger.info("  1. Merging GeneA sequences by Read_ID...")
            merge_sequences_by_header(str(geneA_fasta), str(geneA_collapsed))
            
            # Step 2: Merge GeneB sequences by Read_ID
            self.logger.info("  2. Merging GeneB sequences by Read_ID...")
            merge_sequences_by_header(str(geneB_fasta), str(geneB_collapsed))
            
            # Step 3: Concatenate GeneA and GeneB files
            self.logger.info("  3. Concatenating GeneA and GeneB segments...")
            cmd_cat = ['cat', str(geneA_collapsed), str(geneB_collapsed)]
            with open(str(combined_fasta), 'w') as f:
                subprocess.run(cmd_cat, stdout=f, check=True)
            
            # Step 4: Final merge by Read_ID to create chimeric sequences
            self.logger.info("  4. Final merge to create complete chimeric sequences...")
            merge_sequences_by_header(str(combined_fasta), str(merged_fasta))
            
            # Step 5: Rename sequences for final output
            self.logger.info("  5. Renaming sequences for final output...")
            cmd_rename = ['seqkit', 'rename', str(merged_fasta), '-o', str(final_fasta)]
            subprocess.run(cmd_rename, check=True, capture_output=True, text=True)
            
            # Validate final output
            if not validate_fasta_file(final_fasta):
                raise FileNotFoundError(f"Final FASTA file validation failed: {final_fasta}")
            
            final_count = count_sequences_in_fasta(final_fasta)
            self.logger.info(f"Successfully reconstructed {final_count} chimeric sequences")
            
            return final_fasta
            
        except Exception as e:
            self.logger.error(f"Multi-step sequence merging failed: {e}")
            raise


    def _filter_and_assign_status(self, summary_data: pd.DataFrame, 
                                 chimera_library: pd.DataFrame) -> dict:
        """
        Filter original chimera library by passing reads and assign chromosomal status.
        
        Implements R workflow lines 298-308:
        ChRNAs <- subset(Chimera_library, Chimera_library$Read_ID %in% Summary_total$Read_ID)
        data$Chromosomal_Status <- ifelse(data$Chromosome_Gene_A == data$Chromosome_Gene_B, 
                                          "Intrachromosomal", "Interchromosomal")
        
        Args:
            summary_data: Summary data from Phase 4 with passing Read_IDs
            chimera_library: Original chimera library from Phase 1
            
        Returns:
            dict: Filtered results and output path
        """
        self.logger.info("Step 5.4: Final filtering and chromosomal status assignment...")
        
        try:
            # Get set of reads that passed all processing steps
            passing_reads = set(summary_data['Read_ID'].unique())
            self.logger.info(f"Filtering {len(chimera_library)} total chimeras by {len(passing_reads)} passing reads")
            
            # Filter chimera library to only passing reads
            filtered_chimeras = chimera_library[chimera_library['Read_ID'].isin(passing_reads)].copy()
            self.logger.info(f"Filtered to {len(filtered_chimeras)} passing chimeras")
            
            if len(filtered_chimeras) == 0:
                self.logger.warning("No chimeras passed all filtering steps")
                return {
                    'filtered_chimeras': pd.DataFrame(),
                    'output_path': None
                }
            
            # Add chromosomal status classification
            filtered_chimeras['Chromosomal_Status'] = np.where(
                filtered_chimeras['Chromosome_Gene_A'] == filtered_chimeras['Chromosome_Gene_B'],
                'Intrachromosomal',
                'Interchromosomal'
            )
            
            # Export results as CSV and Excel
            # R line 305: write.xlsx(ChRNAs, file = "All_chRNAs_passing_blast_exon_repair.xlsx")
            output_path = os.path.join(self.work_dir, 'All_chRNAs_passing_blast_exon_repair.csv')
            filtered_chimeras.to_csv(output_path, index=False)
            
            # Also save as Excel format
            excel_output_path = os.path.join(self.work_dir, 'All_chRNAs_passing_blast_exon_repair.xlsx')
            filtered_chimeras.to_excel(excel_output_path, index=False)
            self.logger.info(f"Results exported to: {output_path} and {excel_output_path}")
            
            # Remove duplicate Read_ID/Chimera_ID combinations for final statistics
            # R line 303: data_subset <- data %>% distinct(Read_ID, Chimera_ID, .keep_all=TRUE)
            final_chimeras = filtered_chimeras.drop_duplicates(subset=['Read_ID', 'Chimera_ID'], keep='first')
            self.logger.info(f"Final dataset: {len(final_chimeras)} unique chimeras")
            
            # Generate summary statistics
            intra_count = len(final_chimeras[final_chimeras['Chromosomal_Status'] == 'Intrachromosomal'])
            inter_count = len(final_chimeras[final_chimeras['Chromosomal_Status'] == 'Interchromosomal'])
            self.logger.info(f"Chromosomal status: {intra_count} intrachromosomal, {inter_count} interchromosomal")
            
            return {
                'filtered_chimeras': final_chimeras,  # Use deduplicated for statistics
                'exported_chimeras': filtered_chimeras,  # What was actually exported
                'output_path': output_path
            }
            
        except Exception as e:
            self.logger.error(f"Final filtering and status assignment failed: {e}")
            raise


    def _generate_statistics(self, filtered_results: dict, final_fasta: str) -> dict:
        """
        Generate comprehensive statistics for Phase 5 results.
        
        Args:
            filtered_results: Results from filtering step
            final_fasta: Path to final reconstructed FASTA file
            
        Returns:
            dict: Comprehensive statistics
        """
        try:
            filtered_chimeras = filtered_results['filtered_chimeras']
            
            if len(filtered_chimeras) == 0:
                return {
                    'total_reconstructed_sequences': 0,
                    'total_filtered_chimeras': 0,
                    'intrachromosomal_count': 0,
                    'interchromosomal_count': 0,
                    'reconstruction_success_rate': 0.0
                }
            
            stats = {
                'total_reconstructed_sequences': count_sequences_in_fasta(final_fasta),
                'total_filtered_chimeras': len(filtered_chimeras),
                'intrachromosomal_count': len(filtered_chimeras[filtered_chimeras['Chromosomal_Status'] == 'Intrachromosomal']),
                'interchromosomal_count': len(filtered_chimeras[filtered_chimeras['Chromosomal_Status'] == 'Interchromosomal']),
                'reconstruction_success_rate': (count_sequences_in_fasta(final_fasta) / len(filtered_chimeras)) * 100 if len(filtered_chimeras) > 0 else 0.0
            }
            
            # Add tool-specific statistics if available
            if 'Origin' in filtered_chimeras.columns:
                for tool in ['LongGF', 'Genion', 'JaffaL']:
                    tool_count = len(filtered_chimeras[filtered_chimeras['Origin'] == tool])
                    stats[f'{tool.lower()}_final_count'] = tool_count
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Failed to generate statistics: {e}")
            return {
                'total_reconstructed_sequences': 0,
                'total_filtered_chimeras': 0,
                'intrachromosomal_count': 0,
                'interchromosomal_count': 0,
                'reconstruction_success_rate': 0.0
            }


def run_phase5_sequence_reconstruction(phase4_results: dict, chimera_library: pd.DataFrame,
                                     work_dir: str, config: dict) -> dict:
    """
    Standalone function to run Phase 5 sequence reconstruction.
    
    Args:
        phase4_results: Output from Phase 4 processing
        chimera_library: Original chimera library from Phase 1
        work_dir: Working directory for outputs
        config: TYPHON configuration dictionary
        
    Returns:
        dict: Phase 5 results including reconstructed sequences
    """
    reconstructor = SequenceReconstructor(config, work_dir)
    return reconstructor.reconstruct_sequences(phase4_results, chimera_library) 