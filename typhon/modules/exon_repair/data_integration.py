#!/usr/bin/env python3
"""
Phase 1: Data Integration and Preparation

Implements R code lines 48-96: Load and integrate results from all three tools
with origin tracking and GTF annotation mapping.

Authors: Harry Kane, PhD; Eren Ada, PhD
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"Required packages not available: {e}")
    print("Install with: pip install pandas numpy")
    exit(1)


class DataIntegrator:
    """Handles Phase 1: Data Integration and Preparation for exon repair."""
    
    def __init__(self, config: dict, output_dir: str):
        """
        Initialize data integrator.
        
        Args:
            config: TYPHON configuration dictionary
            output_dir: Main pipeline output directory
        """
        self.config = config
        self.output_dir = output_dir
        self.logger = logging.getLogger(__name__)
        
        # Create exon repair working directory
        self.work_dir = os.path.join(output_dir, 'exon_repair')
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Reference files from config
        self.gtf_file = config['references']['gtf']


    def integrate_tool_data(self, longgf_file: str = None, genion_dir: str = None, 
                           jaffal_file: str = None) -> pd.DataFrame:
        """
        Phase 1: Data Integration and Preparation
        
        Implements R code lines 48-96: Load and integrate results from all three tools
        with origin tracking and GTF annotation mapping.
        
        Returns:
            pd.DataFrame: Integrated chimera library with metadata
        """
        self.logger.info("=" * 50)
        self.logger.info("Phase 1: Data Integration and Preparation")
        self.logger.info("=" * 50)
        
        # Auto-detect input files if not provided
        longgf_file = longgf_file or self._find_longgf_file()
        genion_dir = genion_dir or os.path.join(self.output_dir, 'genion_results')
        jaffal_file = jaffal_file or self._find_jaffal_file()
        
        # Validate input files
        self._validate_input_files(longgf_file, genion_dir, jaffal_file)
        
        # 1.1 Tool Data Integration (R code lines 48-58)
        self.logger.info("Step 1.1: Loading tool results...")
        
        # Load LongGF results (columns 2,3: Read_ID, Chimera_ID)
        longgf_df = self._load_longgf_results(longgf_file)
        
        # Load Genion results (columns 28,8: Read_ID, Chimera_ID - but we'll be flexible)
        genion_df = self._load_genion_results(genion_dir)
        
        # Load JaffaL results (columns 1,2: Read_ID, Chimera_ID)
        jaffal_df = self._load_jaffal_results(jaffal_file)
        
        # 1.2 Tool Origin Tracking (R code lines 70-81)
        self.logger.info("Step 1.2: Adding origin tracking...")
        
        # Add origin column to each dataset
        longgf_df['Origin'] = 'LongGF'
        genion_df['Origin'] = 'Genion'
        jaffal_df['Origin'] = 'JaffaL'
        
        # Merge all tools (R equivalent: rbind)
        total_df = pd.concat([longgf_df, jaffal_df, genion_df], ignore_index=True)
        
        # Data cleaning and deduplication (R code lines 59-62)
        self.logger.info("Step 1.3: Data cleaning and deduplication...")
        
        # Handle NaN values and ensure all Chimera_IDs are strings
        # Note: Genion :: formatting is now standardized during data loading
        total_df = total_df.dropna(subset=['Chimera_ID']).copy()
        total_df['Chimera_ID'] = total_df['Chimera_ID'].astype(str)
        
        # Remove duplicates by Read_ID (keep first occurrence)
        # R equivalent: Total <- Total %>% distinct(Read_ID, .keep_all=TRUE)
        total_df = total_df.drop_duplicates(subset=['Read_ID'], keep='first')
        
        # Filter out multi-gene fusions (3+ genes) before GTF mapping and downstream processing.
        # Keep only 2-gene fusions (exactly one colon in Chimera_ID). Save filtered multi-gene
        # entries to a separate CSV for transparency and downstream analysis if needed.
        colon_counts = total_df['Chimera_ID'].str.count(':')
        multi_gene_mask = colon_counts > 1
        multi_gene_count = int(multi_gene_mask.sum())
        if multi_gene_count > 0:
            try:
                multi_gene_df = total_df.loc[multi_gene_mask].copy()
                multi_gene_file = os.path.join(self.work_dir, 'multi_gene_fusions.csv')
                multi_gene_df.to_csv(multi_gene_file, index=False)
                self.logger.info(f"Saved {multi_gene_count} multi-gene fusions (>=3 genes) to: {multi_gene_file}")
            except Exception as e:
                # Do not fail the pipeline if saving the side file encounters an issue
                self.logger.warning(f"Could not save multi-gene fusions file: {e}")
            
            # Proceed with only 2-gene fusions in the main pipeline
            total_df = total_df.loc[~multi_gene_mask].copy()
            self.logger.info(f"Continuing with {len(total_df)} two-gene fusions (exactly 2 genes)")
        
        # Create context dataframe for provenance tracking (matches filtered set)
        context_df = total_df[['Read_ID', 'Chimera_ID', 'Origin']].copy()
        
        # 1.3 GTF Annotation Mapping (R code lines 84-96)
        self.logger.info("Step 1.4: GTF annotation mapping...")
        
        # Split Chimera_ID into GeneA and GeneB
        # R equivalent: Data[c('GeneA', 'GeneB')] <- str_split_fixed(Data$Chimera_ID, ':', 2)
        # Handle cases where Chimera_ID may not have exactly 2 parts
        split_chimera = total_df['Chimera_ID'].str.split(':', expand=True)
        if split_chimera.shape[1] >= 2:
            total_df['GeneA'] = split_chimera.iloc[:, 0]
            total_df['GeneB'] = split_chimera.iloc[:, 1]
        else:
            # Handle malformed Chimera_IDs by setting them to empty strings
            total_df['GeneA'] = ''
            total_df['GeneB'] = ''
            self.logger.warning(f"Some Chimera_IDs don't contain ':' separator. Setting GeneA/GeneB to empty.")
        
        # Load and process GTF for gene metadata
        gtf_genes = self._load_gtf_gene_metadata()
        
        # Map chromosome and strand information
        # R equivalent: Data$Chromosome_Gene_A <- gtf3[match(Data$GeneA, gtf3$gene_name), "seqnames"]
        # Handle duplicate gene names by keeping first occurrence
        gtf_genes_unique = gtf_genes.drop_duplicates(subset=['gene_name'], keep='first')
        
        total_df['Chromosome_Gene_A'] = total_df['GeneA'].map(
            gtf_genes_unique.set_index('gene_name')['seqnames']
        )
        total_df['Strand_Gene_A'] = total_df['GeneA'].map(
            gtf_genes_unique.set_index('gene_name')['strand']
        )
        total_df['Chromosome_Gene_B'] = total_df['GeneB'].map(
            gtf_genes_unique.set_index('gene_name')['seqnames']
        )
        total_df['Strand_Gene_B'] = total_df['GeneB'].map(
            gtf_genes_unique.set_index('gene_name')['strand']
        )
        
        # Save intermediate files
        self._save_phase1_outputs(total_df, context_df)
        
        self.logger.info(f"Phase 1 completed: {len(total_df)} chimeras from all tools")
        self.logger.info(f"LongGF: {len(longgf_df)}, Genion: {len(genion_df)}, JaffaL: {len(jaffal_df)}")
        
        return total_df


    def _find_longgf_file(self) -> str:
        """Auto-detect LongGF results file (prefer CSV, fallback to Excel)."""
        search_paths = [
            # Try CSV first (faster)
            os.path.join(self.output_dir, 'longgf_results', 'Combined_LongGF_chimera_results_total.csv'),
            os.path.join(self.output_dir, 'Combined_LongGF_chimera_results_total.csv'),
            'Combined_LongGF_chimera_results_total.csv',
            # Fallback to Excel
            os.path.join(self.output_dir, 'longgf_results', 'Combined_LongGF_chimera_results_total.xlsx'),
            os.path.join(self.output_dir, 'Combined_LongGF_chimera_results_total.xlsx'),
            'Combined_LongGF_chimera_results_total.xlsx'
        ]
        
        for path in search_paths:
            if os.path.exists(path):
                return path
        
        raise FileNotFoundError("LongGF results file not found. Expected: Combined_LongGF_chimera_results_total.csv or .xlsx")
    
    
    def _find_jaffal_file(self) -> str:
        """Auto-detect JaffaL results file."""
        search_paths = [
            os.path.join(self.output_dir, 'jaffal_results', 'JaffaL_combined_results.txt'),
            os.path.join(self.output_dir, 'JaffaL_combined_results.txt')
        ]
        
        for path in search_paths:
            if os.path.exists(path):
                return path
        
        raise FileNotFoundError("JaffaL results file not found. Expected: JaffaL_combined_results.txt")
    
    
    def _validate_input_files(self, longgf_file: str, genion_dir: str, jaffal_file: str):
        """Validate that required input files exist."""
        if not os.path.exists(longgf_file):
            raise FileNotFoundError(f"LongGF file not found: {longgf_file}")
        if not os.path.exists(genion_dir):
            raise FileNotFoundError(f"Genion directory not found: {genion_dir}")
        if not os.path.exists(jaffal_file):
            raise FileNotFoundError(f"JaffaL file not found: {jaffal_file}")
        
        self.logger.info(f"Input validation passed")
        self.logger.info(f"  LongGF: {longgf_file}")
        self.logger.info(f"  Genion: {genion_dir}")
        self.logger.info(f"  JaffaL: {jaffal_file}")
    
    
    def _load_longgf_results(self, file_path: str) -> pd.DataFrame:
        """Load LongGF results from CSV or Excel file."""
        # R equivalent: LongGF <- read.xlsx("file.xlsx")
        # LongGF_subset <- LongGF[, c(2, 3)]
        try:
            # Try CSV first (faster), then Excel
            if file_path.endswith('.xlsx') or file_path.endswith('.xls'):
                # Check if CSV version exists
                csv_path = file_path.rsplit('.', 1)[0] + '.csv'
                if os.path.exists(csv_path):
                    self.logger.debug(f"Using CSV version: {csv_path}")
                    df = pd.read_csv(csv_path)
                else:
                    self.logger.debug(f"Using Excel version: {file_path}")
                    df = pd.read_excel(file_path)
            else:
                # Assume CSV
                df = pd.read_csv(file_path)
            
            # R code uses LongGF[, c(2, 3)] - columns 2,3 (1-indexed) = columns 1,2 (0-indexed)
            # But our file has Read_ID, Chimera_ID in columns 1,2 - let's try both approaches
            if 'Read_ID' in df.columns and 'Chimera_ID' in df.columns:
                longgf_subset = df[['Read_ID', 'Chimera_ID']].copy()
                self.logger.info("Using named columns for LongGF (Read_ID, Chimera_ID)")
            else:
                # Use R code exact positions: columns 2,3 (1-indexed) = columns 1,2 (0-indexed)
                if len(df.columns) < 3:
                    raise ValueError("LongGF Excel file must have at least 3 columns for R code compatibility")
                self.logger.warning("Using positional columns for LongGF (R code: columns 2,3)")
                longgf_subset = df.iloc[:, [1, 2]].copy()  # R code columns 2,3 
                longgf_subset.columns = ['Read_ID', 'Chimera_ID']
            
            # Clean data
            longgf_subset = longgf_subset.dropna().copy()
            
            self.logger.info(f"Loaded {len(longgf_subset)} LongGF results")
            return longgf_subset
            
        except Exception as e:
            self.logger.error(f"Failed to load LongGF results: {e}")
            raise
    
    
    def _load_genion_results(self, genion_dir: str) -> pd.DataFrame:
        """
        Load Genion results from .tsv and .fail files.
        
        Standardizes Chimera_ID format by converting :: to : for consistency with other tools.
        This ensures Genion results are properly integrated instead of being filtered out.
        """
        # R equivalent: Genion <- read.xlsx("file.xlsx")
        # Genion_subset <- Genion[, c(28, 8)]
        # Based on actual data: Read_ID is last column, Chimera_ID with :: is column 8 (0-indexed: 7)
        
        genion_chimeras = []
        files_found = []
        
        # Find .tsv files (passed Genion results) and .fail files
        tsv_files = list(Path(genion_dir).glob('*_genion.tsv'))
        fail_files = list(Path(genion_dir).glob('*_genion.tsv.fail'))
        
        # Process TSV files
        for tsv_file in tsv_files:
            try:
                df = pd.read_csv(tsv_file, sep='\t', header=None)
                if len(df.columns) >= 28:
                    # R code: Genion[, c(28, 8)] - column 28 for Read_ID, column 8 for Chimera_ID
                    read_ids = df.iloc[:, 27].dropna().tolist()  # Column 28 (0-indexed: 27)
                    chimera_ids = df.iloc[:, 7].dropna().tolist()  # Column 8 (0-indexed: 7)
                    
                    for read_id, chimera_id in zip(read_ids, chimera_ids):
                        # Standardize :: to : format for consistency with other tools
                        standardized_chimera_id = str(chimera_id).replace('::', ':')
                        genion_chimeras.append({'Read_ID': str(read_id), 'Chimera_ID': standardized_chimera_id})
                    
                    files_found.append(str(tsv_file))
                else:
                    self.logger.warning(f"Genion file {tsv_file} has only {len(df.columns)} columns, expected at least 28")
            except Exception as e:
                self.logger.warning(f"Could not read {tsv_file}: {e}")
        
        # Process .fail files
        for fail_file in fail_files:
            try:
                df = pd.read_csv(fail_file, sep='\t', header=None)
                if len(df.columns) >= 28:
                    # R code: Genion[, c(28, 8)] - column 28 for Read_ID, column 8 for Chimera_ID  
                    read_ids = df.iloc[:, 27].dropna().tolist()  # Column 28 (0-indexed: 27)
                    chimera_ids = df.iloc[:, 7].dropna().tolist()  # Column 8 (0-indexed: 7)
                    
                    for read_id, chimera_id in zip(read_ids, chimera_ids):
                        # Standardize :: to : format for consistency with other tools
                        standardized_chimera_id = str(chimera_id).replace('::', ':')
                        genion_chimeras.append({'Read_ID': str(read_id), 'Chimera_ID': standardized_chimera_id})
                    
                    files_found.append(str(fail_file))
                else:
                    self.logger.warning(f"Genion file {fail_file} has only {len(df.columns)} columns, expected at least 28")
            except Exception as e:
                self.logger.warning(f"Could not read {fail_file}: {e}")
        
        if not genion_chimeras:
            self.logger.warning("No Genion results found")
            return pd.DataFrame(columns=['Read_ID', 'Chimera_ID'])
        
        genion_df = pd.DataFrame(genion_chimeras)
        self.logger.info(f"Loaded {len(genion_df)} Genion results from {len(files_found)} files")
        return genion_df
    
    
    def _load_jaffal_results(self, jaffal_file: str) -> pd.DataFrame:
        """Load JaffaL results from combined results file."""
        # R equivalent: JaffaL <- read.xlsx("file.xlsx")
        # JaffaL_subset <- JaffaL[, c(1, 2)]
        try:
            df = pd.read_csv(jaffal_file, sep='\t')
            
            # JaffaL has headers: transcript (read_id) and fusion_genes (chimera_id)
            if 'transcript' in df.columns and 'fusion_genes' in df.columns:
                jaffal_subset = df[['transcript', 'fusion_genes']].copy()
                jaffal_subset.columns = ['Read_ID', 'Chimera_ID']
                self.logger.info("Using named columns for JaffaL (transcript, fusion_genes)")
            else:
                # Fall back to positional columns (R code uses columns 1,2)
                self.logger.warning("Using positional columns for JaffaL (columns 1,2)")
                if len(df.columns) < 2:
                    raise ValueError("JaffaL file must have at least 2 columns")
                jaffal_subset = df.iloc[:, [0, 1]].copy()
                jaffal_subset.columns = ['Read_ID', 'Chimera_ID']
            
            # Clean data
            jaffal_subset = jaffal_subset.dropna().copy()
            
            self.logger.info(f"Loaded {len(jaffal_subset)} JaffaL results")
            return jaffal_subset
            
        except Exception as e:
            self.logger.error(f"Failed to load JaffaL results: {e}")
            raise
    
    
    def _load_gtf_gene_metadata(self) -> pd.DataFrame:
        """Load GTF file and extract gene metadata."""
        # R equivalent: gtf <- rtracklayer::import("gencode.vM28.annotation.gtf")
        # gtf2 <- as.data.frame(gtf)
        # gtf3 <- subset(gtf2, type == "gene")
        
        self.logger.info(f"Loading GTF gene metadata from {self.gtf_file}")
        
        try:
            # Parse GTF file for gene features
            genes = []
            
            with open(self.gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == 'gene':
                        seqname = fields[0]
                        strand = fields[6]
                        attributes = fields[8]
                        
                        # Parse gene_name from attributes
                        gene_name = None
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if attr.startswith('gene_name'):
                                gene_name = attr.split('"')[1]
                                break
                        
                        if gene_name:
                            genes.append({
                                'gene_name': gene_name,
                                'seqnames': seqname,
                                'strand': strand
                            })
            
            gtf_genes = pd.DataFrame(genes)
            self.logger.info(f"Loaded {len(gtf_genes)} gene annotations")
            return gtf_genes
            
        except Exception as e:
            self.logger.error(f"Failed to load GTF file: {e}")
            raise
    
    
    def _save_phase1_outputs(self, total_df: pd.DataFrame, context_df: pd.DataFrame):
        """Save Phase 1 outputs for downstream processing."""
        # R equivalent file outputs from lines 77-81 and 96
        
        # Save read ID list for sequence extraction
        read_id_file = os.path.join(self.work_dir, 'all_chimera_read_ids.txt')
        with open(read_id_file, 'w') as f:
            for read_id in total_df['Read_ID']:
                f.write(f"{read_id}\n")
        
        # Save read ID and chimera ID pairs
        read_chimera_file = os.path.join(self.work_dir, 'read_chimera_pairs.txt')
        total_df[['Read_ID', 'Chimera_ID']].to_csv(read_chimera_file, sep='\t', index=False, header=False)
        
        # Save context file with origin tracking
        context_file = os.path.join(self.work_dir, 'chimera_context.txt')
        context_df.to_csv(context_file, sep='\t', index=False)
        
        # Save complete chimera library as CSV
        library_file = os.path.join(self.work_dir, 'chimera_library.csv')
        total_df.to_csv(library_file, index=False)
        
        self.logger.info(f"Phase 1 outputs saved:")
        self.logger.info(f"  Read IDs: {read_id_file}")
        self.logger.info(f"  Library: {library_file}") 