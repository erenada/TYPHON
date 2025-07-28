"""
TYPHON Exon Repair Module

Molecular-level sequence reconstruction for chimeric RNA detection.
Replaces simple statistical integration with comprehensive exon repair protocol.

Author: Eren Ada, PhD
"""

from .data_integration import DataIntegrator
from .blast_setup import BlastSetupProcessor
from .transcript_selection import TranscriptSelector
from .sequence_reconstruction import SequenceReconstructor

__all__ = ['DataIntegrator', 'BlastSetupProcessor', 'TranscriptSelector', 'SequenceReconstructor', 'run_exon_repair']


def run_exon_repair(config: dict, output_dir: str, **kwargs) -> dict:
    """
    Main entry point for exon repair protocol using refactored modules.
    
    Args:
        config: TYPHON configuration dictionary
        output_dir: Main pipeline output directory
        **kwargs: Additional parameters (longgf_file, genion_dir, jaffal_file)
    
    Returns:
        dict: Results including reconstructed sequences and statistics
    """
    import logging
    
    logger = logging.getLogger(__name__)
    logger.info("Starting exon repair with refactored modules")
    
    try:
        # Phase 1: Data Integration and Preparation
        logger.info("Phase 1: Data Integration")
        integrator = DataIntegrator(config, output_dir)
        chimera_library = integrator.integrate_tool_data(
            longgf_file=kwargs.get('longgf_file'),
            genion_dir=kwargs.get('genion_dir'),
            jaffal_file=kwargs.get('jaffal_file')
        )
        
        # Phase 2: Sequence Extraction and BLAST Setup
        logger.info("Phase 2: BLAST Setup")
        blast_processor = BlastSetupProcessor(config, output_dir)
        phase2_results = blast_processor.setup_blast_analysis(chimera_library)
        
        # Phase 3: BLAST Analysis and Transcript Selection
        logger.info("Phase 3: Transcript Selection")
        transcript_selector = TranscriptSelector(config, output_dir)
        phase3_results = transcript_selector.analyze_and_select_transcripts(chimera_library, phase2_results)
        
        # Phase 4: Exon Data Processing and Breakpoint Detection
        from .exon_data_processing import run_phase4_exon_processing
        logger.info("Phase 4: Exon Data Processing and Breakpoint Detection")
        phase4_results = run_phase4_exon_processing(
            phase3_results['selected_transcripts_with_order'],
            phase2_results['transcript_metadata']['exons_bed'],
            integrator.work_dir,
            config
        )
        
        # Phase 5: Sequence Reconstruction
        logger.info("Phase 5: Sequence Reconstruction")
        reconstructor = SequenceReconstructor(config, output_dir)
        phase5_results = reconstructor.reconstruct_sequences(phase4_results, chimera_library)
        
        # Return comprehensive results
        results = {
            'chimera_library': chimera_library,
            'phase2_results': phase2_results,
            'phase3_results': phase3_results,
            'phase4_results': phase4_results,
            'phase5_results': phase5_results,
            'work_dir': integrator.work_dir,
            'status': 'complete',
            'final_sequences': phase5_results['reconstructed_sequences'],
            'final_chimeras': phase5_results['filtered_chimeras'],
            'final_chimeras_path': phase5_results['final_chimeras_path'],
            'statistics': {
                'total_chimeras': len(chimera_library),
                'longgf_count': len(chimera_library[chimera_library['Origin'] == 'LongGF']) if 'Origin' in chimera_library.columns else 0,
                'genion_count': len(chimera_library[chimera_library['Origin'] == 'Genion']) if 'Origin' in chimera_library.columns else 0,
                'jaffal_count': len(chimera_library[chimera_library['Origin'] == 'JaffaL']) if 'Origin' in chimera_library.columns else 0,
                'reads_with_breakpoints': phase4_results['n_reads_with_breakpoints'],
                'total_exon_segments': phase4_results['n_total_exon_segments'],
                'final_reconstructed_sequences': phase5_results['statistics']['total_reconstructed_sequences'],
                'final_filtered_chimeras': phase5_results['statistics']['total_filtered_chimeras'],
                'intrachromosomal_count': phase5_results['statistics']['intrachromosomal_count'],
                'interchromosomal_count': phase5_results['statistics']['interchromosomal_count'],
                'reconstruction_success_rate': phase5_results['statistics']['reconstruction_success_rate']
            }
        }
        
        logger.info("Exon repair completed successfully - molecular sequences reconstructed")
        logger.info(f"Final results: {results['statistics']['final_reconstructed_sequences']} sequences, "
                   f"{results['statistics']['final_filtered_chimeras']} filtered chimeras")
        return results
        
    except Exception as e:
        logger.error(f"Exon repair failed: {e}")
        raise 