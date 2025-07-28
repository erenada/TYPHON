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
    
  # JaffaL reference files
  reference_files:
    genome_fasta_gz: ./references/genome.fa.gz
    transcriptome_fasta: ./references/transcripts.fa
    annotation_bed: ./references/annotation.bed
    annotation_tab: ./references/annotation.tab
```

**Exon Repair Protocol:**
```yaml
options:
  enable_integration: true
  overlap_analysis_method: "exon_repair"
  
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
Direct RNA-seq fusion detection using long-read alignments.
- **Input:** FASTQ files, genome FASTA, GTF annotation
- **Output:** SAM alignments, Excel/CSV fusion results

### Custom Genion
Graph-based fusion detection with TYPHON enhancements.
- **Input:** FASTQ files, SAM alignments, processed references
- **Output:** TSV fusion results with read-level detail

### JaffaL
JAFFA-Long pipeline for Nanopore/PacBio data.
- **Input:** FASTQ files, reference genome/transcriptome
- **Output:** Combined fusion results

### Exon Repair
Molecular-level sequence reconstruction protocol.
- **Input:** Results from LongGF, Genion, JaffaL
- **Output:** Validated chimeric sequences, filtered results

## Output Structure

```
results/
├── longgf_results/           # LongGF outputs
├── genion_results/           # Genion outputs
├── jaffal_results/           # JaffaL outputs
├── exon_repair/              # Exon repair outputs
│   ├── blast_results/
│   ├── bed_files/
│   └── reconstructed_sequences/
└── logs/                     # Pipeline logs
```

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

## Requirements

**System:**
- Linux OS (Ubuntu 18.04+)
- 16+ GB RAM (32+ GB recommended)
- 50+ GB free disk space

**Software:**
- Python 3.9+
- Conda/Mamba
- Java 11+ (OpenJDK recommended)
- All bioinformatics tools installed via conda environment

## License

Creative Commons Attribution-NonCommercial 4.0 International License  
**Academic and Research Use Only** 