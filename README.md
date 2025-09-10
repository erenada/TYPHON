# TYPHON: Chimeric RNA Detection Pipeline

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux-green.svg)](https://www.linux.org/)

**Authors:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 1.0.0  

## Overview

TYPHON is a pipeline for chimeric RNA detection from long-read sequencing data (Nanopore/PacBio). It integrates three fusion detection tools (LongGF, Genion, JaffaL) with a five-phase exon repair protocol for sequence-based validation.

**Key Features:**
- **Multi-tool integration** - Combines three fusion detection algorithms
- **Long-read optimized** - Designed for Nanopore and PacBio technologies  
- **Sequence-based validation** - Reconstructs chimeric sequences with breakpoint coordinates
- **Filtering** - Cross-validation and BLAST-based filtering

## About the Name

TYPHON is named after the mythological father of the Chimera. Like the mythological chimera that combines parts from different creatures, this pipeline detects chimeric RNA molecules formed by the fusion of different genes.

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

```bash
# Copy template and edit with your specific paths
cp config_template.yaml config.yaml
# Edit config.yaml with your data paths and settings
```

**Essential setup:**
1. Update file paths in `references:` section (genome, GTF, transcriptome)
2. Set input FASTQ directory path
3. Configure JaffaL reference files (see [JaffaL Reference Setup Guide](docs/jaffal_reference_setup.md))
4. Adjust thread counts and memory for your system

**For detailed configuration options, see: [Configuration Guide](docs/configuration.md)**

## JaffaL Reference Files

JaffaL requires four specific reference files from UCSC databases: genome FASTA (`.fa.gz`), transcriptome FASTA (`.fasta`), annotation BED (`.bed`), and annotation TAB (`.tab`).

**CRITICAL:** All files must use the **same genome build** and **annotation version** as your other reference files.

**For complete download instructions, see: [JaffaL Reference Files Setup Guide](docs/jaffal_reference_setup.md)**

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
Direct RNA-seq fusion detection using minimap2 alignments. Identifies fusion candidates through split alignments with configurable overlap thresholds and pseudogene filtering.

### Custom Genion
Graph-based fusion detection with TYPHON-specific enhancements. Builds splice graphs from SAM alignments with enhanced debug output and failure analysis.

### JaffaL (JAFFA-Long)
Assembly-based fusion detection optimized for long-read data. Uses bpipe workflow with Velvet/Oases assembly and configurable memory management for large datasets.

### Exon Repair Protocol
Five-phase sequence reconstruction: (1) data integration, (2) BLAST setup, (3) transcript selection, (4) exon boundary detection, (5) sequence reconstruction. Produces chimeric sequences with breakpoint coordinates.

## Output Structure

```
results/
├── longgf_results/                    # LongGF outputs
│   ├── *.bam                         # BAM alignment files per sample
│   ├── *.sam                         # SAM alignment files per sample
│   ├── *.log                         # Detailed alignment logs per sample
│   ├── Combined_LongGF_chimera_results_total.xlsx  # Aggregated results (Excel)
│   ├── Combined_LongGF_chimera_results_total.csv   # Aggregated results (CSV)
│   └── Combined_LongGF_chimera_results_with_sample_info.xlsx  # Results with sample tracking
├── genion_results/                   # Genion outputs  
│   ├── *_genion.tsv                 # Main fusion results per sample
│   ├── *_genion.tsv.fail            # Debug output for failed candidates
│   ├── *.paf                        # PAF alignment files per sample
│   └── run_genion.log               # Genion execution log
├── genion_references/               # Processed reference files
│   ├── Genion_modified_gtf_final.gtf  # Modified GTF for Genion
│   ├── selfalign.paf                # Self-alignment reference
│   └── selfalign.tsv                # Self-alignment data
├── jaffal_results/                  # JaffaL outputs
│   └── JaffaL_combined_results.txt # Combined fusion results from all samples
├── exon_repair/                     # Exon repair outputs
│   ├── blast_reference/            # BLAST database files
│   │   ├── All_chimera_db.*        # BLAST database components
│   │   └── my_blast_fasta_reference.fa  # Processed transcriptome reference
│   ├── blast_result/               # BLAST analysis results
│   │   └── chimera_blast_result.txt  # Raw BLAST output
│   ├── modified_exon_repair/       # Exon data processing
│   │   ├── all_exons.bed           # Exon coordinate data
│   │   └── transcripts_for_exon_repair.txt  # Transcript metadata
│   ├── All_chRNAs_passing_blast_exon_repair.csv   # Final validated chimeras (CSV)
│   ├── All_chRNAs_passing_blast_exon_repair.xlsx  # Final validated chimeras (Excel)
│   ├── Merged_seqs_exon_repair_renamed.fa         # Final reconstructed sequences
│   ├── chimera_library.csv         # Integrated chimera data from all tools
│   ├── bed_file_A.bed              # Gene A breakpoint coordinates
│   ├── bed_file_B.bed              # Gene B breakpoint coordinates
│   ├── Fasta_geneA_collapse.fa     # Gene A sequences (collapsed by read ID)
│   ├── Fasta_geneB_collapse.fa     # Gene B sequences (collapsed by read ID)
│   ├── merged_samples.bam          # Multi-sample merged BAM file
│   └── [intermediate processing files]  # Phase-specific outputs and temporary files
└── logs/                           # Pipeline logs
    ├── typhon_main.log             # Main pipeline log
    ├── longgf.log                  # LongGF module log
    ├── genion.log                  # Genion module log
    └── jaffal.log                  # JaffaL module log
```

## Key Output Files

**Primary Results (in `exon_repair/` directory):**
- `All_chRNAs_passing_blast_exon_repair.csv/.xlsx` - Validated chimeras with chromosomal classification, breakpoint coordinates, and tool origin tracking
- `Merged_seqs_exon_repair_renamed.fa` - Reconstructed chimeric sequences with breakpoint coordinates

These are the primary outputs of the exon repair module and can be used for downstream analysis.

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

- **Storage:** Requires 2-3x input FASTQ size for temporary processing space
- **Memory:** Enable `process_samples_sequentially: true` for JaffaL on memory-constrained systems  
- **Optimization:** SSD storage recommended for faster I/O performance

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