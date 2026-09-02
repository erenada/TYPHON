# TYPHON: Chimeric RNA Detection Pipeline

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux-green.svg)](https://www.linux.org/)

**Authors:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 1.0.0  

## Overview

TYPHON is a pipeline for chimeric RNA detection from long-read direct RNA-sequencing data. It integrates three fusion detection tools (LongGF, Genion, JAFFAL) with a multi-phase exon repair protocol for identification of candidate chimeric RNA molecules.

**Key Features:**
- **Multi-tool integration** - Combines three fusion detection algorithms
- **Long-read optimized** - Designed for Nanopore direct RNA sequencing data 
- **Exon repair** - Reconstructs chimeric sequences with breakpoint coordinates
- **Filtering** - BLAST-based filtering to confirm that chimeric RNA reads align to transcripts from designated parent genes

## What's with the name?

TYPHON is named after the mythological father of the Chimera from greek myth. Like the mythological chimera that combines parts from different creatures, this pipeline detects chimeric RNA molecules formed by fusing sequences derived from different genes.

## Installation

### Prerequisites
- Linux (Ubuntu 18.04+)
- Conda/Mamba package manager
- 64+ GB RAM (128+ GB recommended for large datasets or if other processes running etc.)
- Free disk space of roughly 8-10x your gzipped input size
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
1. Download your GENCODE reference files of interest (see https://www.gencodegenes.org/). Currently, TYPHON only supports GENCODE references for the primary reference files, support for ENSEMBL and other formats may be added in the future.
2. Update file paths in the `references:` section of the config.yaml file (genome, GTF, transcriptome). The config.yaml is used to specify all of the information that TYPHON requires.
3. Set input FASTQ directory path.
4. Configure JAFFAL reference files (see [JAFFAL Reference Setup Guide](docs/jaffal_reference_setup.md)).
5. Adjust the thread counts and memory for your system.

## JAFFAL Reference Files

JAFFAL requires four specific reference files from UCSC databases: genome FASTA (`.fa.gz`), transcriptome FASTA (`.fasta`), annotation BED (`.bed`), and annotation TAB (`.tab`).

**CRITICAL:** All files must use the **same genome build** and **annotation version** as your other reference files.

**For complete download instructions, see: [JAFFAL Reference Files Setup Guide](docs/jaffal_reference_setup.md)**

### Complete Setup

After configuring `config.yaml`, run the setup scripts:
```bash
# Setup custom Genion (reads paths from config.yaml)
python setup_genion.py

# Setup JAFFAL (reads paths from config.yaml)
python setup_jaffal.py

# Verify installation
./bin/genion --version  # Should show: 1.2.3-dirty
java -version           # Should show Java 11+
python typhon_main.py --help # Should print the help message
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
Chimeric RNA detection using LongGF. LongGF carries out fusion detection on alignments and was published by Liu *et al*., 2020; doi:10.1186/s12864-020-07207-4. Provides configurable overlap thresholds and pseudogene filtering options.

### Genion (slightly modified for usage with TYPHON)
Chimeric RNA detection using Genion. Genion carries out fusion detection utilizing a dynamic programming exon chaining algorithm and was published by Karaoglanoglu *et al*., 2022; doi:10.1186/s12864-022-08339-5. Please note that Genion is slightly modified during TYPHON setup to enable the inclusion of candidate chimeric RNA reads which normally fail to pass Genion's stringent filtering protocols, as TYPHON opts to be permissive by default with respect to prospective chimeric RNA sequences. As such, this results in a much larger number of potential chimeras being called which is not indicative of typical Genion outputs.

### JAFFAL (Long read fusion detection protocol of JAFFA)
Chimeric RNA detection using JAFFAL. JAFFAL carries out fusion detection using both transcriptomic and genomic alignment and was published by Davidson *et al*., 2022; doi:10.1186/s13073-015-0167-x. Uses bpipe workflow and runs as a complete end-to-end pipeline.

### Exon Repair Protocol
Five-phase sequence reconstruction: (1) data integration, (2) BLAST setup, (3) transcript selection, (4) exon boundary detection, (5) sequence reconstruction. Produces exon-repaired chimeric sequences with breakpoint coordinates.

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
- `All_chRNAs_passing_blast_exon_repair.csv.xlsx` - Validated chimeras with chromosomal classification, breakpoint coordinates, and tool origin tracking
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
- **JAFFAL setup failure:** Check Java 11+ installation and JAFFAL reference files

## Performance Notes

- **Storage:** Requires 2-3x input FASTQ size for temporary processing space
- **Optimization:** SSD storage recommended for faster I/O performance

## Citations and References

TYPHON integrates several published bioinformatics tools. Please cite the original publications when using this pipeline:

### TYPHON Pipeline

**Citation for TYPHON will be provided once the manuscript is published.**

### Core Tools

**LongGF** (Long-read gene fusion detection):
> Liu Q, Hu Y, Stucky A, Fang L, Zhong JF, Wang K. LongGF: computational algorithm and software tool for fast and accurate detection of gene fusions by long-read transcriptome sequencing. BMC Genomics. 2020;21:793. doi:10.1186/s12864-020-07207-4.

**Genion** (Gene fusion detection for long reads):
> Karaoglanoglu F, Chauve C, Hach F. Genion, an accurate tool to detect gene fusion from long transcriptomics reads. BMC Genomics. 2022;23:144. doi:10.1186/s12864-022-08339-5.

**JAFFA/JAFFAL** (Fusion gene detection):
> Davidson NM, Majewski IJ, Oshlack A. JAFFA: High sensitivity transcriptome-focused fusion gene detection. Genome Med. 2015;7:43. doi:10.1186/s13073-015-0167-x.
> 
> Davidson NM, Chen Y, Sadras T, et al. JAFFAL: detecting fusion genes with long-read transcriptome sequencing. Genome Biol. 2022;23:10. doi:10.1186/s13059-021-02588-5.
>
**BLAST+** (Alignment of chimeric transcripts to the reference transcriptome for exon repair):
> Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990;215:403-10. doi:10.1016/S0022-2836(05)80360-2.
>
> Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics. 2009;10:421. doi:10.1186/1471-2105-10-421.
> 
**Minimap2** (Sequence alignment - used for input to/by LongGF, Genion, and JAFFAL):
> Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191.

### Supporting Tools

**SAMtools** (File format conversiona and data processing):
> Danecek P, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi: 10.1093/gigascience/giab008.
>
**BEDTools** (Pulldown of sequence information from reference files):
> Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-2. doi:10.1093/bioinformatics/btq033.
>
**SeqKit** (FASTA/Q file manipulation and data processing):
> Shen W, et al. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE. 2016;11(10):e0163962. doi:10.1371/journal.pone.0163962.
>
**BEDOPS** (Data processing):
> Neph S, et al. BEDOPS: high-performance genomic feature operations. Bioinformatics. 2012;28(14):1919-20. doi: 10.1093/bioinformatics/bts277.
>
**Biopython** (Python operations):
> Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-3. doi: 10.1093/bioinformatics/btp163.
>
**Bowtie2** (Used by JAFFAL):
> Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012;9(4):357-9. doi: 10.1038/nmeth.1923. 
>
**Trimmomatic** (Used by JAFFAL):
> Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114-20. doi: 10.1093/bioinformatics/btu170. 
>
**Velvet** (Used by JAFFAL):
> Zerbino DR, Birney E. Velvet: algorithms for de novo short read assembly using de Bruijn graphs. Genome Res. 2008;18(5):821-9. doi: 10.1101/gr.074492.107. 
>
**Oases** (Used by JAFFAL):
> Schulz MH, et al. Oases: robust de novo RNA-seq assembly across the dynamic range of expression levels. Bioinformatics. 2012;28(8):1086-92. doi: 10.1093/bioinformatics/bts094.
>
**BLAT** (Used by JAFFAL):
> Kent WJ. BLAT--the BLAST-like alignment tool. Genome Res. 2002;12(4):656-64. doi: 10.1101/gr.229202.

## License

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)

Creative Commons Attribution-NonCommercial 4.0 International License  
**Academic and Research Use Only**
