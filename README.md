# TYPHON: Chimeric RNA Detection Pipeline

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux-green.svg)](https://www.linux.org/)

**Authors:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 1.0.0  

## Overview

TYPHON is a comprehensive modular bioinformatics pipeline designed for robust chimeric RNA detection from long-read RNA sequencing data (Nanopore/PacBio). The pipeline integrates multiple complementary fusion detection tools—LongGF, custom Genion, and JaffaL—followed by a sophisticated five-phase exon repair protocol for molecular-level sequence validation.

**Key Features:**
- **Multi-tool integration:** Combines three fusion detection algorithms for comprehensive coverage
- **Long-read optimized:** Specifically designed for Nanopore and PacBio sequencing technologies
- **Molecular validation:** Advanced exon repair protocol reconstructs complete chimeric sequences with precise breakpoints
- **Memory management:** Configurable processing modes to handle large datasets efficiently
- **High-confidence results:** Cross-validation between tools and BLAST-based filtering ensures reliable fusion calls

## About the Name

TYPHON is named after [Typhon](https://en.wikipedia.org/wiki/Typhon), the monstrous father of the Chimera in Greek mythology. According to Hesiod, Typhon and Echidna were the parents of the [Chimera](https://en.wikipedia.org/wiki/Chimera_(mythology))—a fire-breathing hybrid creature composed of different animal parts. The name reflects this pipeline's purpose: detecting and analyzing chimeric RNA molecules, which are hybrid transcripts formed by the fusion of different genes, much like the mythological chimera combines parts from different creatures.

## Installation

### Prerequisites
- Linux (Ubuntu 18.04+)
- Conda/Mamba package manager
- 16+ GB RAM (32+ GB recommended for large datasets)
- 50+ GB free disk space (minimum)
- Java 11+ and Perl rename utility (installation instructions below)

### System Dependencies Installation

**Java 11:**
```bash
# Ubuntu/Debian
sudo apt install openjdk-11-jre

# Verify installation
java -version
```

**Perl Rename Utility:**
The rename utility is included in the conda environment and installed automatically. Alternative system installation:
```bash
# Ubuntu/Debian (if needed)
sudo apt install rename

# Verify it's Perl-based
rename --version
```

### Setup

```bash
# Clone repository
git clone https://github.com/erenada/TYPHON.git
cd TYPHON

# Create conda environment
conda env create -f environment.yml
conda activate typhon_env
```

## Configuration

**IMPORTANT:** Configure the pipeline before running setup scripts, as they read paths from the configuration file.

Copy and edit the configuration template:
```bash
cp config_template.yaml config.yaml
# Edit config.yaml with your data paths and settings
```

**Required configuration steps:**
1. Set correct paths for `input.fastq_dir`, `references.genome`, `references.gtf`, `references.transcriptome`
2. Configure `jaffal.jaffal_dir` and `jaffal.reference_files` paths (see [JaffaL Reference Files Setup Guide](docs/jaffal_reference_setup.md) for download instructions)
3. Set `genion.output_bin_dir` path
4. Adjust thread counts and memory settings for your system

## JaffaL Reference Files

JaffaL requires four specific reference files that must be downloaded separately from UCSC databases:

1. **Genome FASTA** (`.fa.gz`) - Genomic reference sequences
2. **Transcriptome FASTA** (`.fasta`) - Transcript sequences  
3. **Annotation BED** (`.bed`) - Exon coordinates
4. **Annotation TAB** (`.tab`) - Gene/transcript metadata

**CRITICAL:** All files must use the **same genome build** and **annotation version** as your other reference files.

### Quick Setup

```bash
# Create directory
mkdir -p ./references/jaffal

# Download files following the detailed guide
# Update config.yaml with correct paths
```

**For complete download instructions, see: [JaffaL Reference Files Setup Guide](docs/jaffal_reference_setup.md)**

This guide provides step-by-step instructions for downloading from UCSC Genome Browser and Table Browser, with specific settings for each file type.

### Complete Setup

After configuring `config.yaml`, run the setup scripts:
```bash
# Setup custom Genion (reads paths from config.yaml)
python setup_genion.py

# Setup JaffaL (reads paths from config.yaml)
python setup_jaffal.py

# Verify installation
./bin/genion --version  # Should show: 1.2.3-dirty
java -version           # Should show Java 11+
conda run -n typhon_env which minimap2 longgf samtools
conda run -n typhon_env rename --man | head -5  # Check rename utility
```

## Configuration Reference

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
    
    # JaffaL reference files (see docs/jaffal_reference_setup.md)
    reference_files:
      genome_fasta_gz: ./references/jaffal/mm39.fa.gz
      transcriptome_fasta: ./references/jaffal/mm39_gencode_M28.fasta
      annotation_bed: ./references/jaffal/mm39_gencode_M28.bed
      annotation_tab: ./references/jaffal/mm39_gencode_M28.tab
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
- Check system dependencies: `java -version` and `rename --version`

**Runtime Issues:**
- Check logs in `{output_dir}/logs/`
- Use `--debug` for detailed logging
- Validate paths in `config.yaml`
- Use `--dry-run` to test configuration

**Common Errors:**
- **Path not found:** Use absolute paths in configuration
- **Missing SAM files:** Run LongGF before Genion
- **JaffaL setup failure:** Check Java 11+ installation and JaffaL reference files
- **Out of memory errors:** Enable `process_samples_sequentially: true` for JaffaL

## Performance Notes

**Storage Considerations:**
- FASTQ files: 5-50+ GB per sample for long-read data
- Temporary processing space: 2-3x the size of input FASTQ files
- Monitor disk space during execution as intermediate files can be substantial

**Optimization:**
- SSD storage recommended for faster I/O during intensive processing steps
- Ensure sufficient RAM for genome indexing and alignment steps
- Use `process_samples_sequentially: true` for JaffaL on memory-constrained systems

## Citations and References

TYPHON integrates several published bioinformatics tools. Please cite the original publications when using this pipeline:

### TYPHON Pipeline

**Citation for TYPHON will be provided once the manuscript is published.**

### Core Tools

**LongGF** (Long-read Gene Fusion detection):
> Liu Q, Hu Y, Stucky A, Fang L, Zhong JF, Wang K. LongGF: computational algorithm and software tool for fast and accurate detection of gene fusions by long-read transcriptome sequencing. BMC Genomics. 2020;21:793. doi:10.1186/s12864-020-07207-4

**Genion** (Gene fusion detection for long reads):
> Karaoglanoglu F, Chauve C, Hach F. Genion, an accurate tool to detect gene fusion from long transcriptomics reads. BMC Genomics. 2022;23:144. doi:10.1186/s12864-022-08339-5

**JAFFA/JaffaL** (Fusion gene detection):
> Davidson NM, Majewski IJ, Oshlack A. JAFFA: High sensitivity transcriptome-focused fusion gene detection. Genome Med. 2015;7:43. doi:10.1186/s13073-015-0167-x
> 
> Davidson NM, Chen Y, Sadras T, et al. JAFFAL: detecting fusion genes with long-read transcriptome sequencing. Genome Biol. 2022;23:10. doi:10.1186/s13059-021-02588-5

### Supporting Tools

**Minimap2** (Sequence alignment - used by LongGF and JaffaL):
> Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191

## License

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)

Creative Commons Attribution-NonCommercial 4.0 International License  
**Academic and Research Use Only**