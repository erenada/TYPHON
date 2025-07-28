#!/usr/bin/env python3
"""
Phase 3: BLAST Analysis and Transcript Selection

Implements R code lines 190-225: Parse BLAST results, filter transcripts,
and determine chimera gene order.

Author: Eren Ada, PhD
"""

import os
import logging
from typing import Dict, List, Tuple, Optional

try:
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"Required packages not available: {e}")
    print("Install with: pip install pandas numpy")
    exit(1)


class TranscriptSelector:
    """Handles Phase 3: BLAST Analysis and Transcript Selection for exon repair."""
    
    def __init__(self, config: dict, output_dir: str):
        """
        Initialize transcript selector.
        
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


    def analyze_and_select_transcripts(self, chimera_library: pd.DataFrame, 
                                     phase2_results: dict) -> dict:
        """
        Phase 3: BLAST Analysis and Transcript Selection
        
        Implements R code lines 190-225: Parse BLAST results, filter transcripts,
        and determine chimera gene order.
        
        Args:
            chimera_library: Output from Phase 1 with Read_IDs and chimera metadata
            phase2_results: Output from Phase 2 with BLAST results and transcript metadata
            
        Returns:
            dict: Filtered BLAST results with gene order assignments
        """
        self.logger.info("=" * 50)
        self.logger.info("Phase 3: BLAST Analysis and Transcript Selection")
        self.logger.info("=" * 50)
        
        # 3.1 BLAST Results Processing (R code lines 190-198)
        blast_results = self._process_blast_results(phase2_results['blast_results'], chimera_library)
        
        # 3.2 Transcript Filtering (R code lines 207-220)
        filtered_results = self._filter_and_annotate_transcripts(blast_results, phase2_results['transcript_metadata'])
        
        # 3.3 Gene Order Determination (R code lines 221-225)
        ordered_results = self._determine_gene_order_and_selection(filtered_results)
        
        # Save outputs before returning results
        self._save_phase3_outputs(ordered_results)
        
        results = {
            'blast_processed': blast_results,
            'filtered_transcripts': filtered_results,
            'ordered_results': ordered_results,
            'selected_transcripts_with_order': ordered_results,  # Full data with BLAST coordinates
            'ordered_results_csv': getattr(self, 'ordered_results_csv', None),  # CSV path for future use
            'phase3_complete': True
        }
        
        self.logger.info("Phase 3 completed successfully")
        return results


    def _process_blast_results(self, blast_results_file: str, chimera_library: pd.DataFrame) -> pd.DataFrame:
        """
        Process BLAST results and cross-reference with chimera library.
        
        Implements R code lines 190-198:
        Data <- read.table("blast_result.txt", sep = "\t")
        colnames(Data) <- c("Read_ID", "Transcript_ID", "%_identity", ...)
        Data[c('Gene', 'Junk')] <- str_split_fixed(Data$Transcript_ID, '-2+', 2)
        Data$Chimera_ID <- Chimera_library[match(Data$Read_ID, Chimera_library$Read_ID), "Chimera_ID"]
        Data[c('GeneA', 'GeneB')] <- str_split_fixed(Data$Chimera_ID, ':', 2)
        """
        self.logger.info("Step 3.1: Processing BLAST results...")
        
        try:
            # Load BLAST results with standard column names (R code line 191)
            blast_df = pd.read_csv(blast_results_file, sep='\t', header=None)
            blast_df.columns = [
                "Read_ID", "Transcript_ID", "%_identity", "alignment_length", 
                "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", 
                "e.value", "bit.score"
            ]
            
            # Extract gene names from transcript IDs (R code lines 196-197)
            # R equivalent: Data[c('Gene', 'Junk')] <- str_split_fixed(Data$Transcript_ID, '-2+', 2)
            # Split on last dash followed by digits (e.g., "Cd68-201" -> "Cd68", "201")
            split_result = blast_df['Transcript_ID'].str.rsplit('-', n=1, expand=True)
            if split_result.shape[1] == 2:
                blast_df['Gene'] = split_result[0]
                blast_df['Transcript_Number'] = split_result[1]
            else:
                # Fallback: use whole transcript ID as gene name
                blast_df['Gene'] = blast_df['Transcript_ID']
                blast_df['Transcript_Number'] = ''
            
            # Select relevant columns (R code line 198)
            blast_df = blast_df[["Read_ID", "Transcript_ID", "Gene", "%_identity", "alignment_length", 
                               "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", 
                               "e.value", "bit.score"]]
            
            # Cross-reference with chimera library (R code line 200)
            # R equivalent: Data$Chimera_ID <- Chimera_library[match(Data$Read_ID, Chimera_library$Read_ID), "Chimera_ID"]
            blast_df['Chimera_ID'] = blast_df['Read_ID'].map(
                chimera_library.set_index('Read_ID')['Chimera_ID']
            )
            
            # Split Chimera_ID into GeneA and GeneB (R code line 201)
            # R equivalent: Data[c('GeneA', 'GeneB')] <- str_split_fixed(Data$Chimera_ID, ':', 2)
            blast_df[['GeneA', 'GeneB']] = blast_df['Chimera_ID'].str.split(':', expand=True)
            
            # Remove rows without chimera information
            blast_df = blast_df.dropna(subset=['Chimera_ID'])
            
            self.logger.info(f"Processed {len(blast_df)} BLAST hits")
            return blast_df
            
        except Exception as e:
            self.logger.error(f"Failed to process BLAST results: {e}")
            raise


    def _filter_and_annotate_transcripts(self, blast_df: pd.DataFrame, transcript_metadata: dict) -> pd.DataFrame:
        """
        Filter transcripts and add metadata annotations.
        
        Implements R code lines 207-220:
        Data$Type <- Transcripts[match(Data$Transcript_ID, Transcripts$V5), "V8"]
        Data$T_length <- Transcripts[match(Data$Transcript_ID, Transcripts$V5), "V7"]
        Data$Tag <- gtf2[match(Data$Transcript_ID, gtf2$transcript_name), "tag"]
        Data <- Data[!Data$Type=="retained_intron", ]
        Data <- filter(Data, Gene==GeneA | Gene==GeneB)
        """
        self.logger.info("Step 3.2: Filtering and annotating transcripts...")
        
        try:
            # Load transcript metadata (R code lines 203-206)
            transcripts_file = transcript_metadata['transcripts_metadata']
            transcripts_df = pd.read_csv(transcripts_file, sep='|', header=None)
            
            # R equivalent: Transcripts$V9 <- NULL; Transcripts$V1<-gsub(">","",as.character(Transcripts$V1))
            if transcripts_df.shape[1] > 8:  # Only drop if V9 exists
                transcripts_df = transcripts_df.drop(columns=[8])  # V9 equivalent (0-indexed: 8)
            transcripts_df.iloc[:, 0] = transcripts_df.iloc[:, 0].str.replace('>', '')
            
            # Add transcript type and length annotations (R code lines 207-209)
            # R equivalent: Data$Type <- Transcripts[match(Data$Transcript_ID, Transcripts$V5), "V8"]
            # V5 is column 4 (0-indexed), V8 is column 7 (0-indexed), V7 is column 6 (0-indexed)
            transcript_lookup_type = transcripts_df.set_index(transcripts_df.columns[4])[transcripts_df.columns[7]]
            transcript_lookup_length = transcripts_df.set_index(transcripts_df.columns[4])[transcripts_df.columns[6]]
            
            blast_df['Type'] = blast_df['Transcript_ID'].map(transcript_lookup_type)
            blast_df['T_length'] = blast_df['Transcript_ID'].map(transcript_lookup_length)
            
            # Add GTF tag information (R code line 210)
            # R equivalent: Data$Tag <- gtf2[match(Data$Transcript_ID, gtf2$transcript_name), "tag"]
            gtf_tags = self._load_gtf_transcript_tags()
            blast_df['Tag'] = blast_df['Transcript_ID'].map(gtf_tags)
            
            # Filter out retained_intron transcripts (R code line 211)
            # R equivalent: Data <- Data[!Data$Type=="retained_intron", ]
            blast_df = blast_df[blast_df['Type'] != 'retained_intron']
            
            # Keep only transcripts matching GeneA or GeneB (R code line 212)
            # R equivalent: Data <- filter(Data, Gene==GeneA | Gene==GeneB)
            blast_df = blast_df[
                (blast_df['Gene'] == blast_df['GeneA']) | 
                (blast_df['Gene'] == blast_df['GeneB'])
            ]
            
            self.logger.info(f"After filtering: {len(blast_df)} transcript matches")
            return blast_df
            
        except Exception as e:
            self.logger.error(f"Failed to filter and annotate transcripts: {e}")
            raise


    def _load_gtf_transcript_tags(self) -> dict:
        """Load GTF transcript tags for annotation."""
        self.logger.info("Loading GTF transcript tags...")
        
        try:
            transcript_tags = {}
            
            with open(self.gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == 'transcript':
                        attributes = fields[8]
                        
                        # Parse transcript_name and tag from attributes
                        transcript_name = None
                        tag = None
                        
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if attr.startswith('transcript_name'):
                                transcript_name = attr.split('"')[1]
                            elif attr.startswith('tag'):
                                tag = attr.split('"')[1]
                        
                        if transcript_name:
                            transcript_tags[transcript_name] = tag
            
            self.logger.info(f"Loaded {len(transcript_tags)} transcript tags")
            return transcript_tags
            
        except Exception as e:
            self.logger.error(f"Failed to load GTF transcript tags: {e}")
            raise


    def _determine_gene_order_and_selection(self, blast_df: pd.DataFrame) -> pd.DataFrame:
        """
        Determine gene order and select best transcripts.
        
        Implements R code lines 213-225:
        Data$LongGF_Order <- ifelse(Data$Gene==Data$GeneA, "A", "B")
        Data$Pick <- paste(Data$Read_ID, Data$LongGF_Order, sep="_")
        Data$Prefer <- ifelse(Data$Tag=="GENCODE_Primary", "A", "B")
        Data <- arrange(Data, Read_ID, -bit.score, Prefer, desc(T_length))
        Data_subset <- Data %>% group_by(Pick) %>% slice_head(n = 1) %>% ungroup()
        Remove_singles <- Data_subset[Data_subset$Read_ID %in% Data_subset$Read_ID[duplicated(Data_subset$Read_ID)],]
        Remove_singles <- arrange(Remove_singles, desc(Read_ID), q.start)
        Remove_singles$Actual_order <- rep(c("A", "B"), length.out = nrow(Remove_singles))
        """
        self.logger.info("Step 3.3: Determining gene order and selecting transcripts...")
        
        try:
            # Create order assignment (R code line 213)
            # R equivalent: Data$LongGF_Order <- ifelse(Data$Gene==Data$GeneA, "A", "B")
            blast_df['LongGF_Order'] = np.where(
                blast_df['Gene'] == blast_df['GeneA'], 'A', 'B'
            )
            
            # Create Pick identifier for grouping (R code line 214)
            # R equivalent: Data$Pick <- paste(Data$Read_ID, Data$LongGF_Order, sep="_")
            blast_df['Pick'] = blast_df['Read_ID'] + '_' + blast_df['LongGF_Order']
            
            # Create preference for GENCODE_Primary transcripts (R code line 215)
            # R equivalent: Data$Prefer <- ifelse(Data$Tag=="GENCODE_Primary", "A", "B")
            blast_df['Prefer'] = np.where(
                blast_df['Tag'] == 'GENCODE_Primary', 'A', 'B'
            )
            
            # Sort by criteria (R code line 216)
            # R equivalent: Data <- arrange(Data, Read_ID, -bit.score, Prefer, desc(T_length))
            blast_df = blast_df.sort_values([
                'Read_ID', 'bit.score', 'Prefer', 'T_length'
            ], ascending=[True, False, True, False])
            
            # Select best transcript per gene per read (R code lines 217-218)
            # R equivalent: Data_subset <- Data %>% group_by(Pick) %>% slice_head(n = 1) %>% ungroup()
            data_subset = blast_df.groupby('Pick').first().reset_index()
            
            # Keep only reads with both GeneA and GeneB (R code line 219)
            # R equivalent: Remove_singles <- Data_subset[Data_subset$Read_ID %in% Data_subset$Read_ID[duplicated(Data_subset$Read_ID)],]
            duplicate_read_ids = data_subset['Read_ID'][data_subset['Read_ID'].duplicated()].unique()
            remove_singles = data_subset[data_subset['Read_ID'].isin(duplicate_read_ids)]
            
            # Sort and assign actual order (R code lines 220-221)
            # R equivalent: Remove_singles <- arrange(Remove_singles, desc(Read_ID), q.start)
            remove_singles = remove_singles.sort_values(['Read_ID', 'q.start'], ascending=[False, True])
            
            # R equivalent: Remove_singles$Actual_order <- rep(c("A", "B"), length.out = nrow(Remove_singles))
            remove_singles['Actual_order'] = ['A', 'B'] * (len(remove_singles) // 2) + ['A'] * (len(remove_singles) % 2)
            
            self.logger.info(f"Selected transcripts for {len(remove_singles)} reads with both genes")
            
            return remove_singles
            
        except Exception as e:
            self.logger.error(f"Failed to determine gene order and selection: {e}")
            raise


    def _save_phase3_outputs(self, ordered_results: pd.DataFrame):
        """Save Phase 3 outputs for downstream processing."""
        
        # Save ordered results as CSV
        ordered_file = os.path.join(self.work_dir, 'ordered_blast_results.csv')
        ordered_results.to_csv(ordered_file, index=False)
        
        # Save simplified version for reference
        simple_file = os.path.join(self.work_dir, 'phase3_transcript_selection.txt')
        ordered_results[['Read_ID', 'Transcript_ID', 'Gene', 'Chimera_ID', 'Actual_order']].to_csv(
            simple_file, sep='\t', index=False
        )
        
        self.logger.info(f"Phase 3 outputs saved:")
        self.logger.info(f"  Ordered results: {ordered_file}")
        self.logger.info(f"  Simple selection: {simple_file}")
        
        # Store CSV path for Phase 4
        self.ordered_results_csv = ordered_file 