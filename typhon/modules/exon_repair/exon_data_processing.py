"""
Phase 4: Exon Data Processing and Breakpoint Detection
=====================================================

This module implements the complex exon data processing pipeline from the R workflow,
including GTF attribute parsing, breakpoint detection, and BED file generation.

"""

import pandas as pd
import numpy as np
import logging
import re
from pathlib import Path
from typing import Tuple, Dict, Any
import subprocess


def parse_gtf_attributes_to_bed(exon_bed_path: Path) -> pd.DataFrame:
    """
    Parse GTF attributes from BED file using a reliable approach.
    
    Instead of trying to reproduce complex R data.table operations,
    directly parse GTF attributes to extract transcript_name and gene_name.
    
    Args:
        exon_bed_path: Path to exon BED file from convert2bed
        
    Returns:
        DataFrame with standardized exon annotation columns
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Parsing GTF attributes from BED file: {exon_bed_path}")
    
    # Load exon BED file
    exon_df = pd.read_csv(exon_bed_path, sep='\t', header=None, low_memory=False)
    logger.info(f"Loaded {len(exon_df)} total BED entries")
    
    # Filter for exon features (V8 == "exon")
    exon_df = exon_df[exon_df.iloc[:, 7] == "exon"].copy()
    logger.info(f"Filtered to {len(exon_df)} exon entries")
    
    if len(exon_df) == 0:
        raise ValueError("No exon entries found in BED file")
    
    # Extract GTF attributes from column 9 (V10 in R, 0-indexed)
    attributes_col = exon_df.iloc[:, 9].astype(str)
    
    def extract_gtf_attribute(attr_string, attr_name):
        """Extract a specific attribute from GTF attribute string."""
        import re
        pattern = f'{attr_name} "([^"]+)"'
        match = re.search(pattern, attr_string)
        return match.group(1) if match else None
    
    # Extract the attributes we need
    logger.info("Extracting GTF attributes...")
    
    results = []
    for idx, row in exon_df.iterrows():
        attr_string = row.iloc[9]  # GTF attributes column
        
        # Extract key attributes
        transcript_name = extract_gtf_attribute(attr_string, 'transcript_name')
        gene_name = extract_gtf_attribute(attr_string, 'gene_name') 
        transcript_id = extract_gtf_attribute(attr_string, 'transcript_id')
        transcript_type = extract_gtf_attribute(attr_string, 'transcript_type')
        
        # Extract exon_number (can be quoted or unquoted)
        exon_number_match = re.search(r'exon_number (\d+)', attr_string)
        exon_number = exon_number_match.group(1) if exon_number_match else None
        
        # Use transcript_name if available, fall back to transcript_id
        final_transcript_id = transcript_name if transcript_name else transcript_id
        final_gene = gene_name if gene_name else 'Unknown'
        
        results.append({
            'Exon_chromosome': row.iloc[0],
            'Exon_chromosome_start': int(row.iloc[1]),
            'Exon_chromosome_end': int(row.iloc[2]),
            'Exon_ID': row.iloc[3],
            'Exon_strand': row.iloc[5],
            'transcript_type': transcript_type or 'Unknown',
            'Gene': final_gene,
            'Transcript_ID': final_transcript_id,  # This should match BLAST results
            'exon_number': exon_number
        })
    
    result_df = pd.DataFrame(results)
    logger.info(f"Parsed {len(result_df)} exon annotations")
    logger.info(f"Sample transcript IDs: {result_df['Transcript_ID'].head().tolist()}")
    
    return result_df


def calculate_cumulative_exon_lengths(transcript_exon_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate individual and cumulative exon lengths per transcript.
    
    Args:
        transcript_exon_df: DataFrame with exon information joined to transcripts
        
    Returns:
        DataFrame with exon length calculations
    """
    logger = logging.getLogger(__name__)
    
    # Calculate individual exon lengths
    # R equivalent: New$Exon_length <- abs(New$Exon_chromosome_start - New$Exon_chromosome_end)
    transcript_exon_df['Exon_length'] = abs(
        transcript_exon_df['Exon_chromosome_start'] - transcript_exon_df['Exon_chromosome_end']
    )
    
    # Calculate cumulative exon lengths per transcript
    # R equivalent: New2 <- New %>% group_by(Read_ID, Transcript_ID) %>% mutate(Cumulative_exon_length = cumsum(Exon_length))
    transcript_exon_df['Cumulative_exon_length'] = transcript_exon_df.groupby(
        ['Read_ID', 'Transcript_ID']
    )['Exon_length'].cumsum()
    
    # Initialize distance column
    transcript_exon_df['Abs_exon_distance'] = np.nan
    
    logger.info(f"Calculated cumulative exon lengths for {len(transcript_exon_df)} exon entries")
    
    return transcript_exon_df


def calculate_transcript_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate transcript-relative exon coordinates and breakpoint distances.
    
    Implements the complex data.table operations from R code lines 256-267:
    - GeneB logic: Exon_start_in_transcript = cumulative offset + 1
    - GeneA logic: Exon_end_in_transcript = cumulative length
    - Map BLAST coordinates to cumulative exon positions
    
    Args:
        df: DataFrame with cumulative exon length calculations
        
    Returns:
        DataFrame with transcript coordinates and breakpoint distances
    """
    logger = logging.getLogger(__name__)
    
    def calculate_positions_for_group(group):
        """Calculate transcript positions for a single Read_ID/Transcript_ID group."""
        group = group.sort_values('exon_number')
        
        if group['Actual_order'].iloc[0] == 'B':
            # For GeneB: start positions are cumulative offsets + 1
            # R equivalent: c(0, Cumulative_exon_length[-.N]) + 1
            cum_lens = group['Cumulative_exon_length'].shift(1).fillna(0)
            group['Exon_start_in_transcript'] = cum_lens + 1
        else:  # GeneA
            # For GeneA: end positions are cumulative lengths
            group['Exon_end_in_transcript'] = group['Cumulative_exon_length']
        
        return group
    
    # Initialize the coordinate columns before groupby
    df['Exon_start_in_transcript'] = np.nan
    df['Exon_end_in_transcript'] = np.nan
    df['Abs_exon_distance'] = np.nan
    
    # Apply position calculations by group
    # R equivalent: setDT(New2); setkey(New2, Read_ID, Transcript_ID, exon_number)
    df = df.groupby(['Read_ID', 'Transcript_ID']).apply(calculate_positions_for_group).reset_index(drop=True)
    
    # Calculate absolute distances to exon boundaries
    # R equivalent: New2$Abs_exon_distance[cond1]<- abs(New2$Exon_end_in_transcript-New2$s.end)[cond1]
    cond_A = (df['Actual_order'] == 'A') & df['Exon_end_in_transcript'].notna()
    if cond_A.any():
        df.loc[cond_A, 'Abs_exon_distance'] = abs(
            df.loc[cond_A, 'Exon_end_in_transcript'] - df.loc[cond_A, 's.end']
        )
    
    # R equivalent: New2$Abs_exon_distance[cond2]<- abs(New2$Exon_start_in_transcript-New2$s.start)[cond2]
    cond_B = (df['Actual_order'] == 'B') & df['Exon_start_in_transcript'].notna()
    if cond_B.any():
        df.loc[cond_B, 'Abs_exon_distance'] = abs(
            df.loc[cond_B, 'Exon_start_in_transcript'] - df.loc[cond_B, 's.start']
        )
    
    logger.info(f"Calculated transcript coordinates and breakpoint distances for {len(df)} entries")
    
    return df


def select_breakpoint_exons(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Select breakpoint exons with minimum distance for each read.
    
    Args:
        df: DataFrame with breakpoint distance calculations
        
    Returns:
        Tuple of (GeneA breakpoint info, GeneB breakpoint info)
    """
    logger = logging.getLogger(__name__)
    
    # Split data into GeneA and GeneB segments
    # R equivalent: DataA <- as.data.frame(New2[New2$Actual_order=="A", ])
    data_A = df[df['Actual_order'] == 'A'].copy()
    data_B = df[df['Actual_order'] == 'B'].copy()
    
    logger.info(f"Split into {len(data_A)} GeneA and {len(data_B)} GeneB entries")
    
    # Select breakpoint exons (minimum distance)
    # R equivalent: data_newA <- DataA %>% arrange(Abs_exon_distance) %>% group_by(Read_ID) %>% dplyr::slice(1)
    breakpoint_A = data_A.loc[data_A.groupby('Read_ID')['Abs_exon_distance'].idxmin()]
    breakpoint_B = data_B.loc[data_B.groupby('Read_ID')['Abs_exon_distance'].idxmin()]
    
    logger.info(f"Selected breakpoint exons for {len(breakpoint_A)} GeneA and {len(breakpoint_B)} GeneB reads")
    
    return breakpoint_A, breakpoint_B


def determine_exon_ranges(data_A: pd.DataFrame, data_B: pd.DataFrame, 
                         breakpoint_A: pd.DataFrame, breakpoint_B: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Determine exon ranges for sequence reconstruction.
    
    Args:
        data_A: All GeneA exon data
        data_B: All GeneB exon data  
        breakpoint_A: GeneA breakpoint exons
        breakpoint_B: GeneB breakpoint exons
        
    Returns:
        Tuple of (processed GeneA ranges, processed GeneB ranges)
    """
    logger = logging.getLogger(__name__)
    
    # Include all exons from start to breakpoint (GeneA)
    # R equivalent: NewA <- inner_join(DataA, subset_datanew_A, by = "Read_ID")
    subset_A = breakpoint_A[['Read_ID', 'exon_number']].rename(columns={'exon_number': 'breakpoint_exon'})
    merged_A = data_A.merge(subset_A, on='Read_ID', how='inner')
    
    # R equivalent: New_processedA <- subset(NewA, exon_number.x <= exon_number.y)
    processed_A = merged_A[merged_A['exon_number'] <= merged_A['breakpoint_exon']].copy()
    
    # Include all exons from breakpoint to end (GeneB)
    subset_B = breakpoint_B[['Read_ID', 'exon_number']].rename(columns={'exon_number': 'breakpoint_exon'})
    merged_B = data_B.merge(subset_B, on='Read_ID', how='inner')
    
    # R equivalent: New_processedB <- subset(NewB, exon_number.x >= exon_number.y)
    processed_B = merged_B[merged_B['exon_number'] >= merged_B['breakpoint_exon']].copy()
    
    # Sort by Read_ID and exon_number
    # R equivalent: New_processedA <- arrange(New_processedA, Read_ID, exon_number.x)
    processed_A = processed_A.sort_values(['Read_ID', 'exon_number'])
    processed_B = processed_B.sort_values(['Read_ID', 'exon_number'])
    
    # Add quality score column
    processed_A['Quality'] = '.'
    processed_B['Quality'] = '.'
    
    logger.info(f"Determined exon ranges: {len(processed_A)} GeneA, {len(processed_B)} GeneB entries")
    
    return processed_A, processed_B


def create_bed_files(processed_A: pd.DataFrame, processed_B: pd.DataFrame, 
                    output_dir: Path) -> Tuple[Path, Path]:
    """
    Create BED format files for sequence extraction.
    
    Args:
        processed_A: Processed GeneA exon ranges
        processed_B: Processed GeneB exon ranges
        output_dir: Output directory for BED files
        
    Returns:
        Tuple of (GeneA BED path, GeneB BED path)
    """
    logger = logging.getLogger(__name__)
    
    # Create BED format files
    # R equivalent: bed_file_A <- select(New_processedA, "Exon_chromosome", "Exon_chromosome_start", ...)
    bed_A = processed_A[['Exon_chromosome', 'Exon_chromosome_start', 'Exon_chromosome_end', 
                        'Read_ID', 'Quality', 'Exon_strand']].copy()
    bed_B = processed_B[['Exon_chromosome', 'Exon_chromosome_start', 'Exon_chromosome_end', 
                        'Read_ID', 'Quality', 'Exon_strand']].copy()
    
    # Write BED files
    bed_a_path = Path(output_dir) / "bed_file_A.bed"
    bed_b_path = Path(output_dir) / "bed_file_B.bed" 
    
    bed_A.to_csv(bed_a_path, sep='\t', header=False, index=False)
    bed_B.to_csv(bed_b_path, sep='\t', header=False, index=False)
    
    logger.info(f"Created BED files: {bed_a_path} ({len(bed_A)} entries), {bed_b_path} ({len(bed_B)} entries)")
    
    return bed_a_path, bed_b_path


def run_phase4_exon_processing(selected_transcripts: pd.DataFrame, exon_bed_path: Path, 
                              output_dir: Path, config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run Phase 4: Exon Data Processing and Breakpoint Detection.
    
    Args:
        selected_transcripts: Output from Phase 3 transcript selection
        exon_bed_path: Path to exon BED file from convert2bed  
        output_dir: Output directory for results
        config: Configuration dictionary
        
    Returns:
        Dictionary with results including BED file paths and summary data
    """
    logger = logging.getLogger(__name__)
    logger.info("Starting Phase 4: Exon Data Processing and Breakpoint Detection")
    
    try:
        # 4.1: Parse GTF attributes from BED file
        logger.info("Step 4.1: Parsing GTF attributes from BED file")
        exon_annotations = parse_gtf_attributes_to_bed(exon_bed_path)
        
        # 4.2: Join BLAST results with exon data
        logger.info("Step 4.2: Joining transcript data with exon annotations")
        # R equivalent: New <- inner_join(Remove_singles, Data_subset8, by = "Transcript_ID")
        transcript_exon_df = selected_transcripts.merge(exon_annotations, on='Transcript_ID', how='inner')
        
        if len(transcript_exon_df) == 0:
            raise ValueError("No matching transcripts found in exon annotations")
        
        # Convert exon numbers to numeric and sort
        transcript_exon_df['exon_number'] = pd.to_numeric(transcript_exon_df['exon_number'], errors='coerce')
        transcript_exon_df = transcript_exon_df.dropna(subset=['exon_number'])
        transcript_exon_df = transcript_exon_df.sort_values(['Transcript_ID', 'exon_number'])
        
        logger.info(f"Joined data: {len(transcript_exon_df)} transcript-exon associations")
        
        # 4.3: Calculate cumulative exon lengths
        logger.info("Step 4.3: Calculating cumulative exon lengths")
        transcript_exon_df = calculate_cumulative_exon_lengths(transcript_exon_df)
        
        # 4.4: Calculate transcript coordinates and breakpoint distances
        logger.info("Step 4.4: Calculating transcript coordinates and breakpoint detection")
        transcript_exon_df = calculate_transcript_coordinates(transcript_exon_df)
        
        # 4.5: Select breakpoint exons
        logger.info("Step 4.5: Selecting breakpoint exons")
        breakpoint_A, breakpoint_B = select_breakpoint_exons(transcript_exon_df)
        
        # Determine exon ranges for reconstruction
        data_A = transcript_exon_df[transcript_exon_df['Actual_order'] == 'A'].copy()
        data_B = transcript_exon_df[transcript_exon_df['Actual_order'] == 'B'].copy()
        
        processed_A, processed_B = determine_exon_ranges(data_A, data_B, breakpoint_A, breakpoint_B)
        
        # Create summary file
        summary_total = pd.concat([processed_A, processed_B], ignore_index=True)
        
        # Clean up summary file
        columns_to_remove = ['Quality', 'Pick', 'Prefer']
        for col in columns_to_remove:
            if col in summary_total.columns:
                summary_total = summary_total.drop(columns=[col])
        
        # Remove duplicate Gene column if present
        gene_cols = [col for col in summary_total.columns if col.startswith('Gene') and col.endswith('.y')]
        if gene_cols:
            summary_total = summary_total.drop(columns=gene_cols)
        
        # 4.6: Create BED files
        logger.info("Step 4.6: Creating BED files for sequence extraction")
        bed_a_path, bed_b_path = create_bed_files(processed_A, processed_B, output_dir)
        
        # Save summary file as CSV
        summary_path = Path(output_dir) / "phase4_summary_exon_data.csv"
        summary_total.to_csv(summary_path, index=False)
        
        # Export debug data as CSV
        debug_path = Path(output_dir) / "phase4_transcript_exon_debug.csv"
        transcript_exon_df.to_csv(debug_path, index=False)
        
        results = {
            'bed_file_A': bed_a_path,
            'bed_file_B': bed_b_path, 
            'summary_data': summary_total,
            'summary_path': summary_path,
            'transcript_exon_data': transcript_exon_df,
            'debug_path': debug_path,
            'breakpoint_A': breakpoint_A,
            'breakpoint_B': breakpoint_B,
            'processed_A': processed_A,
            'processed_B': processed_B,
            'n_reads_with_breakpoints': len(summary_total['Read_ID'].unique()),
            'n_total_exon_segments': len(summary_total)
        }
        
        logger.info(f"Phase 4 completed successfully:")
        logger.info(f"  - {results['n_reads_with_breakpoints']} reads with identified breakpoints")
        logger.info(f"  - {results['n_total_exon_segments']} total exon segments for reconstruction")
        logger.info(f"  - BED files: {bed_a_path}, {bed_b_path}")
        logger.info(f"  - Summary: {summary_path}")
        
        return results
        
    except Exception as e:
        logger.error(f"Phase 4 failed: {str(e)}")
        raise 