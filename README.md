# TYPHON: Chimeric RNA Detection Pipeline

**Authors:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 1.0.0  

Modular bioinformatics pipeline for chimeric RNA detection from long-read RNA sequencing data, integrating LongGF, custom Genion, JaffaL, and exon repair protocols.

## Installation

### Prerequisites
- Linux (Ubuntu 18.04+)
- Conda/Mamba package manager
- Java 11+ (for JaffaL)
- 16+ GB RAM (32+ GB recommended)

### Setup

```bash
# Clone repository
git clone https://github.com/erenada/TYPHON.git
cd TYPHON

# Create conda environment
conda env create -f environment.yml
conda activate typhon_env

# Setup custom Genion
python setup_genion.py

# Setup JaffaL (optional)
python setup_jaffal.py

# Verify installation
./bin/genion --version  # Should show: 1.2.3-dirty
```

## Configuration

Copy and edit the configuration template:
```bash
cp config_template.yaml config.yaml
# Edit config.yaml with your data paths
```

### Key Configuration Sections

**Input/Output:**
```yaml
project:
  name: Your_Analysis_Name
  output_dir: ./results
  threads: 20

input:
  fastq_dir: ./data/fastq_files

references:
  genome: ./references/genome.fa
  gtf: ./references/annotation.gtf
  transcriptome: ./references/transcripts.fa
```

**Pipeline Modules:**
```yaml
modules:
  longgf:
    enabled: true
    min_overlap_len: 100
    
  genion:
    enabled: true
    min_support: 1
    keep_debug: true
    
  jaffal:
    enabled: true
    jaffal_dir: ./jaffal/JAFFA-version-2.3
    
    # Memory management
    max_memory: "28G"
    bpipe_memory: "24G"
    process_samples_sequentially: true
    
    # JaffaL reference files
    reference_files:
      genome_fasta_gz: ./test_data/FILES_FOR_JAFFAL/mm39.fa.gz
      transcriptome_fasta: ./test_data/FILES_FOR_JAFFAL/mm39_gencode_M28.fasta
      annotation_bed: ./test_data/FILES_FOR_JAFFAL/mm39_gencode_M28.bed
      annotation_tab: ./test_data/FILES_FOR_JAFFAL/mm39_gencode_M28.tab
```

**Exon Repair Protocol:**
```yaml
options:
  enable_integration: true
  cleanup_intermediate: true
  debug: true
  
  exon_repair:
    enabled: true
    keep_intermediate: true
    blast_threads: 20
    min_blast_identity: 80
```

## Usage

### Basic Execution
```bash
# Run complete pipeline
python typhon_main.py

# Run with custom config
python typhon_main.py --config my_config.yaml

# Run specific modules
python typhon_main.py --modules longgf genion

# Override configuration
python typhon_main.py --threads 30 --output ./custom_results

# Dry run (test configuration)
python typhon_main.py --dry-run
```

### Command-Line Options
- `--config`, `-c`: Configuration file (default: config.yaml)
- `--threads`, `-t`: Override thread count
- `--output`, `-o`: Override output directory
- `--modules`: Run specific modules (longgf, genion, jaffal)
- `--debug`: Enable debug logging (default)
- `--no-debug`: Disable debug logging
- `--dry-run`: Show execution plan without running

## Pipeline Modules

### LongGF
Direct RNA-seq fusion detection using long-read alignments with minimap2.
- **Input:** FASTQ files, genome FASTA, GTF annotation
- **Process:** Aligns reads to genome, identifies fusion candidates through split alignments
- **Parameters:** Configurable overlap length (100bp default), pseudogene filtering, minimum support reads
- **Output:** SAM alignments, Excel/CSV fusion results with read-level evidence
- **Post-processing:** R-based result aggregation and chimera classification

### Custom Genion
Graph-based fusion detection with TYPHON-specific enhancements and debug output.
- **Input:** FASTQ files, SAM alignments (from LongGF), processed reference transcriptome
- **Process:** Builds splice graphs, detects fusion events using custom binary with enhanced logging
- **Features:** Self-alignment PAF generation, TSV-formatted detailed output, failure analysis (.fail files)
- **Output:** TSV fusion results with read-level detail, comprehensive debug information
- **Integration:** Custom compilation with debug flags for detailed fusion characterization

### JaffaL (JAFFA-Long)
JAFFA pipeline optimized for Nanopore/PacBio long-read data using bpipe workflow.
- **Input:** FASTQ files, reference genome/transcriptome, annotation databases
- **Process:** Assembly-based fusion detection with transcript reconstruction
- **Tools:** Integrates Velvet/Oases assembly, Bowtie2/Minimap2 alignment, custom fusion calling
- **Memory Management:** Configurable sequential processing to prevent memory overload on large datasets
- **Features:** Bpipe memory allocation control, per-sample cleanup, aggressive garbage collection
- **Output:** Combined fusion results with confidence scoring and breakpoint resolution
- **Validation:** Cross-references with known fusion databases and genomic repeat regions

### Exon Repair Protocol
Five-phase molecular-level sequence reconstruction for chimeric RNA validation.
- **Phase 1:** Data integration from LongGF, Genion, and JaffaL results
- **Phase 2:** BLAST database setup and sequence extraction from fusion candidates
- **Phase 3:** Transcript selection using BLAST analysis and confidence scoring
- **Phase 4:** Exon boundary detection and breakpoint-aware coordinate calculation
- **Phase 5:** Sequence reconstruction using bedtools getfasta and multi-step merging
- **Output:** Validated chimeric sequences with precise breakpoint coordinates, filtered fusion library
- **Features:** Handles complex splice variants, validates fusion feasibility, generates high-confidence chimeric sequences

## Output Structure

```
results/
├── longgf_results/                    # LongGF outputs
│   ├── *.sam                         # Alignment files for each sample
│   ├── *.log                         # Detailed alignment logs
│   ├── *_results.txt                 # Processed fusion candidates
│   └── Combined_LongGF_chimera_results_total.xlsx  # Aggregated results
├── genion_results/                   # Genion outputs  
│   ├── *_genion.tsv                 # Main fusion results per sample
│   ├── *_genion.tsv.fail            # Debug output for failed candidates
│   └── genion_references/           # Processed reference files
├── jaffal_results/                  # JaffaL outputs
│   ├── jaffa_results.csv           # Primary fusion calls
│   ├── *.fastq/                    # Per-sample bpipe output directories
│   └── overlap_analysis.xlsx       # Cross-tool comparison results
├── exon_repair/                     # Exon repair outputs
│   ├── blast_results/              # BLAST analysis files
│   ├── bed_files/                  # Breakpoint coordinate files
│   ├── reconstructed_sequences/    # Final chimeric sequences
│   │   ├── *_geneA.fa             # 5' partner sequences
│   │   ├── *_geneB.fa             # 3' partner sequences  
│   │   └── *_merged.fa            # Complete chimeric sequences
│   └── phase[1-5]_results/         # Intermediate processing outputs
└── logs/                           # Pipeline logs
    ├── typhon.log                  # Main pipeline log
    ├── longgf.log                  # Module-specific logs
    ├── genion.log
    └── jaffal.log
```

## Final Results

The pipeline produces two primary final outputs in the `exon_repair/` directory:

### Validated Chimeric Sequences
- **File:** `Merged_seqs_exon_repair_renamed.fa`
- **Format:** FASTA
- **Content:** High-confidence reconstructed chimeric sequences with precise breakpoint coordinates
- **Features:** 5' and 3' gene partners merged into complete chimeric transcripts

### Comprehensive Chimera Metadata
- **Files:** 
  - `All_chRNAs_passing_blast_exon_repair.xlsx` (Excel format)
  - `All_chRNAs_passing_blast_exon_repair.csv` (CSV format)
- **Content:** Validated chimeras that passed all filtering and BLAST analysis steps
- **Includes:** 
  - Chromosomal classification (Intrachromosomal/Interchromosomal)
  - BLAST validation metrics
  - Breakpoint coordinates and exon boundaries
  - Read-level evidence and support statistics
  - Tool origin tracking (LongGF, Genion, JaffaL)

These files represent the final, publication-ready results with molecular-level validation and can be used directly for downstream analysis, visualization, or experimental validation.

## Troubleshooting

**Setup Issues:**
- Ensure conda environment is activated: `conda activate typhon_env`
- Run setup scripts in order: `setup_genion.py` then `setup_jaffal.py`
- Verify Genion binary: `./bin/genion --version`

**Runtime Issues:**
- Check logs in `{output_dir}/logs/`
- Use `--debug` for detailed logging
- Validate paths in `config.yaml`
- Use `--dry-run` to test configuration

**Common Errors:**
- **Path not found:** Use absolute paths in configuration
- **Missing SAM files:** Run LongGF before Genion
- **JaffaL setup failure:** Ensure Java 11+ is installed
- **Out of memory errors:** Enable `process_samples_sequentially: true` for JaffaL or reduce `max_memory` settings

## Requirements

**System:**
- Linux OS (Ubuntu 18.04+)
- 16+ GB RAM (32+ GB recommended for large datasets)
- Storage requirements:
  - 50+ GB free disk space (minimum)
  - Additional space for FASTQ files (typically 5-50+ GB per sample for long-read data)
  - Temporary processing space: 2-3x the size of input FASTQ files
  - Consider that pipeline generates multiple intermediate files during processing

**Software:**
- Python 3.9+
- Conda/Mamba
- Java 11+ (OpenJDK recommended)
- All bioinformatics tools installed via conda environment

**Performance Notes:**
- Long-read FASTQ files are typically large (several GB to 50+ GB per sample)
- Ensure sufficient RAM for genome indexing and alignment steps
- SSD storage recommended for faster I/O during intensive processing steps
- Monitor disk space during pipeline execution as intermediate files can be substantial

## License

Creative Commons Attribution-NonCommercial 4.0 International License  
**Academic and Research Use Only** 