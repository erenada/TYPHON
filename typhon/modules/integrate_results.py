#!/usr/bin/env python3
"""
Results Integration Module for TYPHON Pipeline

Handles 3-way overlap analysis between LongGF, Genion, and JaffaL results.
This module is separated from run_jaffal.py to allow for independent improvements
and enhancements to the integration logic.

Features:
- Standardized chimera ID formatting across tools
- Flexible input file detection
- Comprehensive overlap statistics
- Extensible architecture for additional analysis

Author: Eren Ada, PhD
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import pandas as pd
except ImportError:
    pd = None
import argparse


class ChimeraResult:
    """Container for chimera detection results from a single tool."""
    
    def __init__(self, tool_name: str, chimera_ids: List[str], source_file: Optional[str] = None):
        self.tool_name = tool_name
        self.chimera_ids = set(self._standardize_chimera_ids(chimera_ids))
        self.source_file = source_file
        self.count = len(self.chimera_ids)
    
    def _standardize_chimera_ids(self, chimera_ids: List[str]) -> List[str]:
        """Standardize chimera ID format across tools."""
        standardized = []
        for chimera_id in chimera_ids:
            # Convert :: to : (Genion format standardization)
            standardized_id = str(chimera_id).replace('::', ':')
            standardized.append(standardized_id)
        return standardized
    
    def __str__(self):
        return f"{self.tool_name}: {self.count} chimeras"


class OverlapAnalyzer:
    """Performs overlap analysis between multiple chimera detection tools."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.results = {}
        self.overlap_stats = {}
    
    def add_result(self, result: ChimeraResult):
        """Add results from a detection tool."""
        self.results[result.tool_name] = result
        logging.info(f"Added {result}")
    
    def calculate_overlaps(self) -> Dict:
        """Calculate all possible overlaps between tools."""
        tools = list(self.results.keys())
        overlaps = {}
        
        # Pairwise overlaps
        for i, tool1 in enumerate(tools):
            for j, tool2 in enumerate(tools[i+1:], i+1):
                overlap_key = f"{tool1}_∩_{tool2}"
                overlap = self.results[tool1].chimera_ids & self.results[tool2].chimera_ids
                overlaps[overlap_key] = {
                    'chimeras': overlap,
                    'count': len(overlap),
                    'tools': [tool1, tool2]
                }
        
        # Three-way overlap (if we have 3 tools)
        if len(tools) >= 3:
            three_way = set.intersection(*[self.results[tool].chimera_ids for tool in tools])
            overlaps['LongGF_∩_Genion_∩_JaffaL'] = {
                'chimeras': three_way,
                'count': len(three_way),
                'tools': tools
            }
        
        self.overlap_stats = overlaps
        return overlaps
    
    def generate_statistics_report(self) -> str:
        """Generate a comprehensive statistics report."""
        report = []
        report.append("=" * 60)
        report.append("TYPHON Results Integration - Statistics Report")
        report.append("=" * 60)
        
        # Individual tool results
        report.append("\nIndividual Tool Results:")
        report.append("-" * 30)
        for tool_name, result in self.results.items():
            report.append(f"{tool_name:>12}: {result.count:>6} chimeras")
        
        # Overlap analysis
        if self.overlap_stats:
            report.append("\nOverlap Analysis:")
            report.append("-" * 30)
            for overlap_name, data in self.overlap_stats.items():
                report.append(f"{overlap_name:>25}: {data['count']:>6} chimeras")
        
        # Calculate precision (overlap ratios)
        if 'LongGF_∩_Genion_∩_JaffaL' in self.overlap_stats:
            three_way_count = self.overlap_stats['LongGF_∩_Genion_∩_JaffaL']['count']
            report.append("\nConsensus Analysis:")
            report.append("-" * 30)
            for tool_name, result in self.results.items():
                if result.count > 0:
                    precision = (three_way_count / result.count) * 100
                    report.append(f"{tool_name} precision: {precision:.1f}% ({three_way_count}/{result.count})")
        
        report.append("=" * 60)
        return "\n".join(report)


def load_longgf_results(excel_file: str) -> ChimeraResult:
    """Load LongGF results from Excel file."""
    logging.info(f"Loading LongGF results from {excel_file}")
    
    if pd is None:
        raise ImportError("pandas is required for Excel file processing. Install with: pip install pandas")
    
    try:
        df = pd.read_excel(excel_file)
        
        if 'Chimera_ID' not in df.columns:
            raise ValueError("LongGF Excel file missing 'Chimera_ID' column")
        
        chimera_ids = df['Chimera_ID'].dropna().tolist()
        logging.info(f"Found {len(chimera_ids)} LongGF chimeras")
        
        return ChimeraResult("LongGF", chimera_ids, excel_file)
        
    except Exception as e:
        logging.error(f"Failed to load LongGF results: {e}")
        raise


def load_genion_results(output_dir: str) -> ChimeraResult:
    """Load Genion results from .tsv and .fail files."""
    logging.info(f"Loading Genion results from {output_dir}")
    
    # Look for Genion results in multiple locations
    search_dirs = [
        os.path.join(output_dir, 'genion_results'),
        output_dir,
        os.path.join(output_dir, '..', 'test_output_pipeline', 'genion_results')
    ]
    
    genion_chimeras = []
    files_found = []
    
    for search_dir in search_dirs:
        if not os.path.exists(search_dir):
            continue
            
        # Find .tsv files (passed Genion results)
        tsv_files = list(Path(search_dir).glob('*_genion.tsv'))
        fail_files = list(Path(search_dir).glob('*_genion.tsv.fail'))
        
        # Also check for combined files
        total_passed = os.path.join(search_dir, 'total_passed_genion_reads.txt')
        total_failed = os.path.join(search_dir, 'total_failed_genion_reads.txt')
        
        if os.path.exists(total_passed):
            tsv_files.append(Path(total_passed))
        if os.path.exists(total_failed):
            fail_files.append(Path(total_failed))
        
        # Process TSV files
        for tsv_file in tsv_files:
            try:
                df = pd.read_csv(tsv_file, sep='\t', header=None)
                if len(df.columns) >= 2:
                    chimeras = df.iloc[:, 1].dropna().tolist()  # V2 column
                    genion_chimeras.extend(chimeras)
                    files_found.append(str(tsv_file))
            except Exception as e:
                logging.warning(f"Could not read {tsv_file}: {e}")
        
        # Process .fail files
        for fail_file in fail_files:
            try:
                df = pd.read_csv(fail_file, sep='\t', header=None)
                if len(df.columns) >= 2:
                    chimeras = df.iloc[:, 1].dropna().tolist()  # V2 column
                    genion_chimeras.extend(chimeras)
                    files_found.append(str(fail_file))
            except Exception as e:
                logging.warning(f"Could not read {fail_file}: {e}")
        
        if files_found:
            break
    
    if not files_found:
        raise FileNotFoundError(f"No Genion result files found in any search directory")
    
    logging.info(f"Found {len(genion_chimeras)} Genion chimeras from {len(files_found)} files")
    
    return ChimeraResult("Genion", genion_chimeras, f"{len(files_found)} files")


def load_jaffal_results(jaffal_results_file: str) -> ChimeraResult:
    """Load JaffaL results from combined results file."""
    logging.info(f"Loading JaffaL results from {jaffal_results_file}")
    
    try:
        df = pd.read_csv(jaffal_results_file, sep='\t')
        
        # Look for fusion_genes column or similar
        fusion_col = None
        for col in df.columns:
            if 'fusion' in col.lower() or 'gene' in col.lower():
                fusion_col = col
                break
        
        if fusion_col is None:
            # Try first column if no fusion column found
            fusion_col = df.columns[0]
            logging.warning(f"No 'fusion_genes' column found, using {fusion_col}")
        
        chimera_ids = df[fusion_col].dropna().tolist()
        logging.info(f"Found {len(chimera_ids)} JaffaL chimeras")
        
        return ChimeraResult("JaffaL", chimera_ids, jaffal_results_file)
        
    except Exception as e:
        logging.error(f"Failed to load JaffaL results: {e}")
        raise


def call_r_script_integration(output_dir: str, longgf_excel: str, jaffal_results: str) -> Dict[str, str]:
    """
    Call the original R script for backward compatibility.
    
    This maintains compatibility with the original pipeline while allowing
    for enhanced Python-based analysis.
    """
    logging.info("Running original R script for backward compatibility...")
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(script_dir, 'Combine_chimera_results_and_create_overlap_file.R')
    
    if not os.path.exists(r_script):
        raise FileNotFoundError(f"R script not found: {r_script}")
    
    try:
        cmd = ['Rscript', r_script, output_dir, longgf_excel, jaffal_results]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        logging.info("R script integration completed")
        logging.debug(f"R script output: {result.stdout}")
        
        # Return paths to generated files
        return {
            'overlap_excel': os.path.join(output_dir, 'Overlapping_mRNA_chimeras.xlsx'),
            'read_id_file': os.path.join(output_dir, 'Read_ID_Overlap.txt'),
            'read_chimera_convert': os.path.join(output_dir, 'Read_and_Chimera_ID_convert.txt')
        }
        
    except subprocess.CalledProcessError as e:
        logging.error(f"R script failed: {e}")
        logging.error(f"stderr: {e.stderr}")
        raise


def integrate_results(output_dir: str, longgf_excel: str = None, jaffal_results: str = None, 
                     use_r_script: bool = True, generate_report: bool = True) -> Dict:
    """
    Main integration function that performs comprehensive overlap analysis.
    
    Args:
        output_dir: Main pipeline output directory
        longgf_excel: Path to LongGF Excel results (auto-detected if None)
        jaffal_results: Path to JaffaL results (auto-detected if None)
        use_r_script: Whether to also run the original R script
        generate_report: Whether to generate statistics report
    
    Returns:
        Dictionary with integration results and statistics
    """
    logging.info("=" * 50)
    logging.info("Starting Results Integration Analysis")
    logging.info("=" * 50)
    
    # Auto-detect input files if not provided
    if not longgf_excel:
        longgf_search_paths = [
            os.path.join(output_dir, 'longgf_results', 'Combined_LongGF_chimera_results_total.xlsx'),
            os.path.join(output_dir, '..', 'test_output_longgf', 'Combined_LongGF_chimera_results_total.xlsx'),
            'Combined_LongGF_chimera_results_total.xlsx'
        ]
        
        for path in longgf_search_paths:
            if os.path.exists(path):
                longgf_excel = path
                break
    
    if not jaffal_results:
        jaffal_search_paths = [
            os.path.join(output_dir, 'jaffal_results', 'JaffaL_combined_results.txt'),
            os.path.join(output_dir, 'JaffaL_combined_results.txt')
        ]
        
        for path in jaffal_search_paths:
            if os.path.exists(path):
                jaffal_results = path
                break
    
    # Validate required files
    if not longgf_excel or not os.path.exists(longgf_excel):
        raise FileNotFoundError(f"LongGF Excel file not found: {longgf_excel}")
    
    if not jaffal_results or not os.path.exists(jaffal_results):
        raise FileNotFoundError(f"JaffaL results file not found: {jaffal_results}")
    
    # Load results from all tools
    analyzer = OverlapAnalyzer(output_dir)
    
    try:
        longgf_result = load_longgf_results(longgf_excel)
        analyzer.add_result(longgf_result)
        
        genion_result = load_genion_results(output_dir)
        analyzer.add_result(genion_result)
        
        jaffal_result = load_jaffal_results(jaffal_results)
        analyzer.add_result(jaffal_result)
        
    except Exception as e:
        logging.error(f"Failed to load tool results: {e}")
        raise
    
    # Perform overlap analysis
    overlaps = analyzer.calculate_overlaps()
    
    # Generate and save statistics report
    if generate_report:
        report = analyzer.generate_statistics_report()
        logging.info(f"\n{report}")
        
        report_file = os.path.join(output_dir, 'integration_statistics.txt')
        with open(report_file, 'w') as f:
            f.write(report)
        logging.info(f"Statistics report saved to: {report_file}")
    
    # Run original R script for backward compatibility
    r_script_files = {}
    if use_r_script:
        try:
            r_script_files = call_r_script_integration(output_dir, longgf_excel, jaffal_results)
        except Exception as e:
            logging.warning(f"R script integration failed: {e}")
    
    logging.info("=" * 50)
    logging.info("Results Integration Completed")
    
    if '3-way' in overlaps or 'LongGF_∩_Genion_∩_JaffaL' in overlaps:
        three_way_key = 'LongGF_∩_Genion_∩_JaffaL' if 'LongGF_∩_Genion_∩_JaffaL' in overlaps else '3-way'
        logging.info(f"High-confidence chimeras: {overlaps[three_way_key]['count']}")
    
    logging.info("=" * 50)
    
    return {
        'overlaps': overlaps,
        'analyzer': analyzer,
        'r_script_files': r_script_files,
        'statistics_report': report if generate_report else None
    }


def main():
    """Standalone integration tool."""
    parser = argparse.ArgumentParser(
        description="TYPHON Results Integration Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect all files
  python integrate_results.py --output-dir ./test_output_pipeline
  
  # Specify specific files
  python integrate_results.py --output-dir ./test_output_pipeline \
    --longgf-excel ./longgf_results/Combined_LongGF_chimera_results_total.xlsx \
    --jaffal-results ./jaffal_results/JaffaL_combined_results.txt
  
  # Python-only analysis (skip R script)
  python integrate_results.py --output-dir ./test_output_pipeline --no-r-script
        """
    )
    
    parser.add_argument('--output-dir', required=True,
                        help='Main pipeline output directory')
    parser.add_argument('--longgf-excel',
                        help='Path to LongGF Excel results file')
    parser.add_argument('--jaffal-results',
                        help='Path to JaffaL combined results file')
    parser.add_argument('--no-r-script', action='store_true',
                        help='Skip running the original R script')
    parser.add_argument('--no-report', action='store_true',
                        help='Skip generating statistics report')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        results = integrate_results(
            output_dir=args.output_dir,
            longgf_excel=args.longgf_excel,
            jaffal_results=args.jaffal_results,
            use_r_script=not args.no_r_script,
            generate_report=not args.no_report
        )
        
        # Print summary to console
        if '3-way' in results['overlaps'] or 'LongGF_∩_Genion_∩_JaffaL' in results['overlaps']:
            three_way_key = 'LongGF_∩_Genion_∩_JaffaL' if 'LongGF_∩_Genion_∩_JaffaL' in results['overlaps'] else '3-way'
            overlap_count = results['overlaps'][three_way_key]['count']
            print(f"\nSUCCESS: Found {overlap_count} high-confidence chimeric RNAs")
        
    except Exception as e:
        logging.error(f"Integration failed: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main() 