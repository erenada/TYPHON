#!/usr/bin/env python3
"""
Phase 2: Sequence Extraction and BLAST Setup

Implements R code lines 99-140: Extract sequences from BAM files and prepare BLAST database.

Authors: Harry Kane, PhD; Eren Ada, PhD
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"Required packages not available: {e}")
    print("Install with: pip install pandas numpy")
    sys.exit(1)


class BlastSetupProcessor:
    """Handles Phase 2: Sequence Extraction and BLAST Setup for exon repair."""
    
    def __init__(self, config: dict, output_dir: str):
        """
        Initialize BLAST setup processor.
        
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
        self.genome_fasta = config['references']['genome']
        self.gtf_file = config['references']['gtf']
        self.transcriptome_fasta = config['references']['transcriptome']
        
        # Processing parameters
        self.threads = config['project'].get('threads', 4)


    def setup_blast_analysis(self, chimera_library: pd.DataFrame) -> dict:
        """
        Phase 2: Sequence Extraction and BLAST Setup
        
        Implements R code lines 99-140: Extract sequences from BAM files and prepare BLAST database.
        
        Args:
            chimera_library: Output from Phase 1 with Read_IDs
            
        Returns:
            dict: Paths to extracted sequences and BLAST database
        """
        self.logger.info("=" * 50)
        self.logger.info("Phase 2: Sequence Extraction and BLAST Setup")
        self.logger.info("=" * 50)
        
        # 2.1 BAM to FASTA Extraction (R code lines 107-113)
        fasta_file = self._extract_sequences_from_bam(chimera_library)
        
        # 2.2 BLAST Database Preparation (R code lines 116-127)
        blast_db_path = self._prepare_blast_database()
        
        # 2.3 Transcript Metadata Processing (R code lines 135-140, 201-206)
        transcript_metadata = self._prepare_transcript_metadata()
        
        # 2.4 BLAST Execution (R code lines 133-134)
        blast_results = self._run_blast_analysis(fasta_file, blast_db_path)
        
        results = {
            'sequences_fasta': fasta_file,
            'blast_database': blast_db_path,
            'transcript_metadata': transcript_metadata,
            'blast_results': blast_results
        }
        
        self.logger.info("Phase 2 completed successfully")
        return results


    def _extract_sequences_from_bam(self, chimera_library: pd.DataFrame) -> str:
        """
        Extract sequences from BAM files using Read_IDs.
        
        Implements R code lines 107-113:
        samtools view -@ 28 -h -N read_ids.txt merged.bam > output.sam
        samtools collate output.sam -o collated.sam  
        samtools fasta collated.sam > sequences.fa
        seqkit sort sequences.fa > sorted.fa
        seqkit rmdup -n sorted.fa > final.fa
        
        If multiple BAM files exist, they are merged first to match original workflow.
        """
        self.logger.info("Step 2.1: BAM to FASTA extraction...")
        
        # Find BAM files
        bam_files = self._find_bam_files()
        if not bam_files:
            raise FileNotFoundError("No BAM files found for sequence extraction")
        
        # Merge BAM files if multiple exist (matches original workflow)
        if len(bam_files) > 1:
            bam_file = self._merge_bam_files(bam_files)
            self.logger.info(f"Merged {len(bam_files)} BAM files into: {bam_file}")
        else:
            bam_file = bam_files[0]
            self.logger.info(f"Using single BAM file: {bam_file}")
        
        # Read IDs file from Phase 1
        read_ids_file = os.path.join(self.work_dir, 'all_chimera_read_ids.txt')
        if not os.path.exists(read_ids_file):
            raise FileNotFoundError(f"Read IDs file not found: {read_ids_file}")
        
        # Output files
        sam_file = os.path.join(self.work_dir, 'extracted_reads.sam')
        collated_sam = os.path.join(self.work_dir, 'collated_reads.sam')
        fasta_unsorted = os.path.join(self.work_dir, 'sequences_unsorted.fa')
        fasta_sorted = os.path.join(self.work_dir, 'sequences_sorted.fa')
        final_fasta = os.path.join(self.work_dir, 'sequences_final.fa')
        
        try:
            # Step 1: samtools view -@ threads -h -N read_ids.txt bam_file > sam_file
            cmd1 = [
                'samtools', 'view', f'-@{self.threads}', '-h', '-N', 
                read_ids_file, bam_file
            ]
            self.logger.info(f"Running: {' '.join(cmd1)} > {sam_file}")
            with open(sam_file, 'w') as f:
                result = subprocess.run(cmd1, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            
            # Step 2: samtools collate sam_file -o collated_sam
            cmd2 = ['samtools', 'collate', sam_file, '-o', collated_sam]
            self.logger.info(f"Running: {' '.join(cmd2)}")
            subprocess.run(cmd2, check=True, capture_output=True, text=True)
            
            # Step 3: samtools fasta collated_sam > fasta_unsorted
            cmd3 = ['samtools', 'fasta', collated_sam]
            self.logger.info(f"Running: {' '.join(cmd3)} > {fasta_unsorted}")
            with open(fasta_unsorted, 'w') as f:
                subprocess.run(cmd3, stdout=f, check=True, text=True)
            
            # Step 4: seqkit sort fasta_unsorted > fasta_sorted
            cmd4 = ['seqkit', 'sort', fasta_unsorted]
            self.logger.info(f"Running: {' '.join(cmd4)} > {fasta_sorted}")
            with open(fasta_sorted, 'w') as f:
                subprocess.run(cmd4, stdout=f, check=True, text=True)
            
            # Step 5: seqkit rmdup -n fasta_sorted > final_fasta
            cmd5 = ['seqkit', 'rmdup', '-n', fasta_sorted]
            self.logger.info(f"Running: {' '.join(cmd5)} > {final_fasta}")
            with open(final_fasta, 'w') as f:
                subprocess.run(cmd5, stdout=f, check=True, text=True)
            
            self.logger.info(f"Sequence extraction completed: {final_fasta}")
            
            # Create Original_ONT file with Chimera_ID headers for QC compatibility
            self._create_original_ont_file(final_fasta)
            
            return final_fasta
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Sequence extraction failed: {e}")
            if e.stderr:
                self.logger.error(f"Error details: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in sequence extraction: {e}")
            raise


    def _find_bam_files(self) -> List[str]:
        """Find BAM files for sequence extraction."""
        # Check config first
        exon_repair_config = self.config.get('options', {}).get('exon_repair', {})
        bam_file = exon_repair_config.get('bam_file')
        
        if bam_file and os.path.exists(bam_file):
            return [bam_file]
        
        # Auto-detect BAM files in output directory
        bam_files = []
        
        # Look for BAM files in main output directory
        for bam_file in Path(self.output_dir).glob('*.bam'):
            bam_files.append(str(bam_file))
        
        # Also check longgf_results subdirectory
        longgf_dir = os.path.join(self.output_dir, 'longgf_results')
        if os.path.exists(longgf_dir):
            for bam_file in Path(longgf_dir).glob('*.bam'):
                bam_files.append(str(bam_file))
        
        self.logger.info(f"Found {len(bam_files)} BAM files: {bam_files}")
        return bam_files


    def _merge_bam_files(self, bam_files: list) -> str:
        """
        Merge multiple BAM files into a single merged BAM file.
        
        Implements the original workflow's use of a merged BAM file containing
        all samples, which is essential for multi-sample exon repair.
        
        Args:
            bam_files: List of BAM file paths to merge
            
        Returns:
            str: Path to merged BAM file
        """
        merged_bam_unsorted = os.path.join(self.work_dir, 'merged_samples_unsorted.bam')
        merged_bam = os.path.join(self.work_dir, 'merged_samples.bam')
        
        try:
            # Step 1: Use samtools merge to combine all BAM files
            cmd_merge = ['samtools', 'merge', '-@', str(self.threads), merged_bam_unsorted] + bam_files
            self.logger.info(f"Merging BAM files: {' '.join(cmd_merge)}")
            subprocess.run(cmd_merge, check=True, capture_output=True, text=True)
            
            # Step 2: Sort the merged BAM file (required for indexing)
            cmd_sort = ['samtools', 'sort', '-@', str(self.threads), '-o', merged_bam, merged_bam_unsorted]
            self.logger.info(f"Sorting merged BAM: {' '.join(cmd_sort)}")
            subprocess.run(cmd_sort, check=True, capture_output=True, text=True)
            
            # Step 3: Index the sorted merged BAM file
            cmd_index = ['samtools', 'index', merged_bam]
            self.logger.info(f"Indexing sorted merged BAM: {' '.join(cmd_index)}")
            subprocess.run(cmd_index, check=True, capture_output=True, text=True)
            
            # Clean up intermediate unsorted file
            if os.path.exists(merged_bam_unsorted):
                os.remove(merged_bam_unsorted)
            
            self.logger.info(f"Successfully merged, sorted, and indexed {len(bam_files)} BAM files into: {merged_bam}")
            return merged_bam
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BAM file merging/sorting failed: {e}")
            if e.stderr:
                self.logger.error(f"Error details: {e.stderr}")
            raise


    def _prepare_blast_database(self) -> str:
        """
        Prepare BLAST database from transcriptome.
        
        Implements R code lines 116-127: Complex transcriptome processing pipeline
        cut -d"|" -f5- transcripts.fa > blast_ref_reduced_first_line.fa
        sed '/|/ s/^/>/' blast_ref_reduced_first_line.fa > blast_ref_reduced_first_line2.fa
        sed 's/|.*//' blast_ref_reduced_first_line2.fa > blast_ref_reduced_first_line3.fa
        seqkit rename blast_ref_reduced_first_line3.fa > blast_ref_reduced_first_line4.fa
        # R script processing (we'll implement this in Python)
        makeblastdb -in processed.fa -parse_seqids -blastdb_version 5 -title "All_chimera_db" -dbtype nucl -out All_chimera_db
        """
        self.logger.info("Step 2.2: BLAST database preparation...")
        
        # Create blast directory
        blast_dir = os.path.join(self.work_dir, 'blast_reference')
        os.makedirs(blast_dir, exist_ok=True)
        
        # Input transcriptome FASTA
        transcriptome_fasta = self.transcriptome_fasta
        if not os.path.exists(transcriptome_fasta):
            raise FileNotFoundError(f"Transcriptome FASTA not found: {transcriptome_fasta}")
        
        # Intermediate files
        blast_ref1 = os.path.join(blast_dir, 'blast_ref_reduced_first_line.fa')
        blast_ref2 = os.path.join(blast_dir, 'blast_ref_reduced_first_line2.fa')
        blast_ref3 = os.path.join(blast_dir, 'blast_ref_reduced_first_line3.fa')
        blast_ref4 = os.path.join(blast_dir, 'blast_ref_reduced_first_line4.fa')
        final_blast_ref = os.path.join(blast_dir, 'my_blast_fasta_reference.fa')
        blast_db_prefix = os.path.join(blast_dir, 'All_chimera_db')
        
        try:
            # Step 1: cut -d"|" -f5- transcripts.fa > blast_ref_reduced_first_line.fa
            cmd1 = ['cut', '-d|', '-f5-', transcriptome_fasta]
            self.logger.info(f"Running: {' '.join(cmd1)} > {blast_ref1}")
            with open(blast_ref1, 'w') as f:
                subprocess.run(cmd1, stdout=f, check=True, text=True)
            
            # Step 2: sed '/|/ s/^/>/' blast_ref_reduced_first_line.fa > blast_ref_reduced_first_line2.fa
            cmd2 = ['sed', '/|/ s/^/>/', blast_ref1]
            self.logger.info(f"Running: {' '.join(cmd2)} > {blast_ref2}")
            with open(blast_ref2, 'w') as f:
                subprocess.run(cmd2, stdout=f, check=True, text=True)
            
            # Step 3: sed 's/|.*//' blast_ref_reduced_first_line2.fa > blast_ref_reduced_first_line3.fa
            cmd3 = ['sed', 's/|.*//', blast_ref2]
            self.logger.info(f"Running: {' '.join(cmd3)} > {blast_ref3}")
            with open(blast_ref3, 'w') as f:
                subprocess.run(cmd3, stdout=f, check=True, text=True)
            
            # Step 4: seqkit rename blast_ref_reduced_first_line3.fa > blast_ref_reduced_first_line4.fa
            cmd4 = ['seqkit', 'rename', blast_ref3]
            self.logger.info(f"Running: {' '.join(cmd4)} > {blast_ref4}")
            with open(blast_ref4, 'w') as f:
                subprocess.run(cmd4, stdout=f, check=True, text=True)
            
            # Step 5: R script processing - implement the header standardization in Python
            # Original R script: Blast_reference_setup.R 
            # We'll implement this logic directly in Python
            self._process_blast_reference_headers(blast_ref4, final_blast_ref)
            
            # Step 6: makeblastdb
            cmd6 = [
                'makeblastdb', 
                '-in', final_blast_ref,
                '-parse_seqids',
                '-blastdb_version', '5',
                '-title', 'All_chimera_db',
                '-dbtype', 'nucl',
                '-out', blast_db_prefix
            ]
            self.logger.info(f"Running: {' '.join(cmd6)}")
            subprocess.run(cmd6, check=True, capture_output=True, text=True)
            
            # Cleanup intermediate files
            for temp_file in [blast_ref1, blast_ref2, blast_ref3, blast_ref4]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            
            self.logger.info(f"BLAST database created: {blast_db_prefix}")
            return blast_db_prefix
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST database preparation failed: {e}")
            if e.stderr:
                self.logger.error(f"Error details: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in BLAST database preparation: {e}")
            raise


    def _process_blast_reference_headers(self, input_fasta: str, output_fasta: str):
        """
        Process BLAST reference headers to standardize format and remove duplicates.
        
        Replaces the R script: Blast_reference_setup.R
        This removes duplicate sequences based on lowercase sequence names,
        keeping only the first occurrence of each unique sequence name.
        """
        self.logger.info("Processing BLAST reference headers and removing duplicates...")
        
        try:
            sequences = {}
            seen_lowercase = set()
            duplicate_count = 0
            
            with open(input_fasta, 'r') as infile:
                current_header = None
                current_sequence = []
                
                for line in infile:
                    line = line.strip()
                    if line.startswith('>'):
                        # Process previous sequence if exists
                        if current_header is not None:
                            seq_name = current_header[1:]  # Remove >
                            seq_name_lower = seq_name.lower()
                            
                            if seq_name_lower not in seen_lowercase:
                                # Keep this sequence (first occurrence)
                                sequences[current_header] = ''.join(current_sequence)
                                seen_lowercase.add(seq_name_lower)
                            else:
                                # Skip duplicate
                                duplicate_count += 1
                        
                        # Start new sequence
                        current_header = line
                        current_sequence = []
                    else:
                        # Add sequence line
                        current_sequence.append(line)
                
                # Process final sequence
                if current_header is not None:
                    seq_name = current_header[1:]  # Remove >
                    seq_name_lower = seq_name.lower()
                    
                    if seq_name_lower not in seen_lowercase:
                        sequences[current_header] = ''.join(current_sequence)
                        seen_lowercase.add(seq_name_lower)
                    else:
                        duplicate_count += 1
            
            # Write deduplicated sequences
            with open(output_fasta, 'w') as outfile:
                for header, sequence in sequences.items():
                    outfile.write(f"{header}\n")
                    # Write sequence in 60-character lines (matching R script nbchar=60)
                    for i in range(0, len(sequence), 60):
                        outfile.write(f"{sequence[i:i+60]}\n")
            
            self.logger.info(f"BLAST reference processed: {len(sequences)} unique sequences kept, {duplicate_count} duplicates removed")
            
        except Exception as e:
            self.logger.error(f"Failed to process BLAST reference headers: {e}")
            raise


    def _prepare_transcript_metadata(self) -> dict:
        """
        Prepare transcript metadata for filtering.
        
        Implements R code lines 135-140, 201-206:
        awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' transcripts.fa > transcripts_for_exon_repair_prep.txt
        sed '/^[0-9]/d' transcripts_for_exon_repair_prep.txt > transcripts_for_exon_repair.txt
        convert2bed --input=gtf --attribute-key=exon_id --max-mem=32G < annotation.gtf > all_exons.bed
        """
        self.logger.info("Step 2.3: Transcript metadata preparation...")
        
        # Create subdirectory for exon repair files
        exon_repair_dir = os.path.join(self.work_dir, 'modified_exon_repair')
        os.makedirs(exon_repair_dir, exist_ok=True)
        
        # Output files
        transcripts_prep = os.path.join(exon_repair_dir, 'transcripts_for_exon_repair_prep.txt')
        transcripts_final = os.path.join(exon_repair_dir, 'transcripts_for_exon_repair.txt')
        exons_bed = os.path.join(exon_repair_dir, 'all_exons.bed')
        
        try:
            # Step 1: Extract transcript lengths using awk
            # awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' transcripts.fa
            self._extract_transcript_lengths(self.transcriptome_fasta, transcripts_prep)
            
            # Step 2: Remove numeric lines (lengths) keeping only headers
            # sed '/^[0-9]/d' transcripts_for_exon_repair_prep.txt > transcripts_for_exon_repair.txt
            cmd2 = ['sed', '/^[0-9]/d', transcripts_prep]
            self.logger.info(f"Running: {' '.join(cmd2)} > {transcripts_final}")
            with open(transcripts_final, 'w') as f:
                subprocess.run(cmd2, stdout=f, check=True, text=True)
            
            # Step 3: Convert GTF to BED format for exons
            # convert2bed --input=gtf --attribute-key=exon_id --max-mem=32G < annotation.gtf > all_exons.bed
            cmd3 = [
                'convert2bed', 
                '--input=gtf', 
                '--attribute-key=exon_id',
                f'--max-mem={min(32, self.config.get("options", {}).get("max_memory_gb", 32))}G'
            ]
            self.logger.info(f"Running: {' '.join(cmd3)} < {self.gtf_file} > {exons_bed}")
            with open(self.gtf_file, 'r') as infile, open(exons_bed, 'w') as outfile:
                subprocess.run(cmd3, stdin=infile, stdout=outfile, check=True, text=True)
            
            # Cleanup prep file
            if os.path.exists(transcripts_prep):
                os.remove(transcripts_prep)
            
            results = {
                'transcripts_metadata': transcripts_final,
                'exons_bed': exons_bed,
                'exon_repair_dir': exon_repair_dir
            }
            
            self.logger.info("Transcript metadata preparation completed")
            return results
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Transcript metadata preparation failed: {e}")
            if e.stderr:
                self.logger.error(f"Error details: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in transcript metadata preparation: {e}")
            raise


    def _extract_transcript_lengths(self, fasta_file: str, output_file: str):
        """
        Extract transcript lengths from FASTA file.
        
        Implements the awk command:
        awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' file.fa
        """
        self.logger.info("Extracting transcript lengths...")
        
        try:
            with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
                seqlen = 0
                first_header = True
                
                for line in infile:
                    line = line.strip()
                    if line.startswith('>'):
                        # Print previous sequence length (except for first header)
                        if not first_header and seqlen > 0:
                            outfile.write(f"{seqlen}\n")
                        # Print header
                        outfile.write(f"{line}\n")
                        seqlen = 0
                        first_header = False
                    else:
                        # Add sequence length
                        seqlen += len(line)
                
                # Print final sequence length
                if seqlen > 0:
                    outfile.write(f"{seqlen}\n")
            
            self.logger.info(f"Transcript lengths extracted: {output_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to extract transcript lengths: {e}")
            raise


    def _run_blast_analysis(self, query_fasta: str, blast_db: str) -> str:
        """
        Run BLAST analysis against the prepared database.
        
        Implements R code lines 133-134:
        blastn -query sequences.fa -db blast_db -outfmt 6 -num_threads 30 -out blast_results.txt
        """
        self.logger.info("Step 2.4: Running BLAST analysis...")
        
        # Create blast results directory
        blast_results_dir = os.path.join(self.work_dir, 'blast_result')
        os.makedirs(blast_results_dir, exist_ok=True)
        
        # Output file
        blast_output = os.path.join(blast_results_dir, 'chimera_blast_result.txt')
        
        # Get blast threads from config
        blast_threads = self.config.get('options', {}).get('exon_repair', {}).get('blast_threads', self.threads)
        
        try:
            # Run BLAST: blastn -query fasta -db blast_db -outfmt 6 -num_threads threads -out results
            cmd = [
                'blastn',
                '-query', query_fasta,
                '-db', blast_db,
                '-outfmt', '6',
                '-num_threads', str(blast_threads),
                '-out', blast_output
            ]
            
            self.logger.info(f"Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Check if results file was created and has content
            if os.path.exists(blast_output) and os.path.getsize(blast_output) > 0:
                # Count number of BLAST hits
                with open(blast_output, 'r') as f:
                    hit_count = sum(1 for line in f)
                self.logger.info(f"BLAST analysis completed: {hit_count} hits found")
            else:
                self.logger.warning("BLAST analysis completed but no hits found")
            
            return blast_output
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST analysis failed: {e}")
            if e.stderr:
                self.logger.error(f"BLAST error details: {e.stderr}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in BLAST analysis: {e}")
            raise 
    
    
    def _create_original_ont_file(self, sequences_file: str):
        """
        Create Original_ONT_overlapping_mRNA_chimeras_fasta.fa with Chimera_ID headers.
        
        Converts Read_IDs to Chimera_IDs using the read_chimera_pairs.txt mapping,
        matching the original bash script functionality.
        """
        try:
            import pandas as pd
            from Bio import SeqIO
            
            # Load Read_ID → Chimera_ID mapping
            mapping_file = os.path.join(self.work_dir, 'read_chimera_pairs.txt')
            if not os.path.exists(mapping_file):
                self.logger.warning("Read-chimera mapping file not found, skipping Original_ONT creation")
                return
            
            # Read the mapping (tab-separated: Read_ID → Chimera_ID)
            mapping_df = pd.read_csv(mapping_file, sep='\t', header=None, names=['Read_ID', 'Chimera_ID'])
            id_mapping = dict(zip(mapping_df['Read_ID'], mapping_df['Chimera_ID']))
            
            # Create Original_ONT file with converted headers
            original_ont_file = os.path.join(self.work_dir, 'Original_ONT_overlapping_mRNA_chimeras_fasta.fa')
            
            sequences_converted = 0
            with open(original_ont_file, 'w') as out_f:
                for record in SeqIO.parse(sequences_file, 'fasta'):
                    read_id = record.id
                    if read_id in id_mapping:
                        chimera_id = id_mapping[read_id]
                        # Write with Chimera_ID header
                        out_f.write(f">{chimera_id}\n{record.seq}\n")
                        sequences_converted += 1
                    else:
                        # Keep original Read_ID if no mapping found
                        out_f.write(f">{read_id}\n{record.seq}\n")
                        sequences_converted += 1
            
            self.logger.info(f"Created Original_ONT file with {sequences_converted} sequences: {original_ont_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to create Original_ONT file: {e}")
            # Continue without failing the whole process