# Typhon: Modular Pipeline for Chimeric RNA Detection

**Developers:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 0.1.0  

A modular bioinformatics pipeline designed to provide robust, comprehensive chimeric RNA detection for long-read RNA sequencing data using tools such as LongGF, Genion (custom build), and JaffaL.

## Project Status

🚧 **Under Development** - This project is currently being refactored and improved.

## Directory Structure

```
Typhon_pipeline/
├── typhon/                    # Main pipeline modules
│   ├── modules/              # Pipeline step modules (longgf, genion, jaffal)
│   ├── utils/                # Utility functions
│   └── scripts/              # R scripts and other supporting scripts
├── tests/                    # Test suite
│   ├── unit/                 # Unit tests
│   ├── integration/          # Integration tests
│   ├── e2e/                  # End-to-end tests
│   └── data/                 # Test data
├── examples/                 # Example configurations and data
│   ├── configs/              # Example YAML configurations
│   └── data/                 # Small example datasets
├── docs/                     # Documentation
├── scripts/                  # Development and setup scripts
├── bin/                      # Compiled binaries (custom Genion)
├── Genion_files/             # Custom Genion patch and source files
└── .archive/                 # Archive for old files during refactoring
```

## Development Roadmap

See `pipeline_dev_roadmap.md` for detailed development progress tracking.

## Installation

### Prerequisites
- Java 11 (required for JaffaL)
- Perl rename utility (automatically included in conda environment)

See [System Dependencies](docs/system_dependencies.md) for detailed installation instructions.

### Quick Start
```bash
# Clone the repository
git clone https://github.com/erenada/TYPHON.git
cd TYPHON

# Install Java 11 (if not already installed)
sudo apt install openjdk-11-jre  # Ubuntu/Debian

# Create conda environment (includes all bioinformatics tools and dependencies)
conda env create -f environment.yml
conda activate typhon_env

# Build custom Genion with debug mode enabled
python setup_genion.py

# Setup JaffaL (optional - only if you plan to use JaffaL module)
# IMPORTANT: Configure the 'jaffal' section in config.yaml first with your reference file paths
# See examples/mouse_config.yaml for a template
python setup_jaffal.py --debug  # --debug provides detailed logging (optional)

# Verify installation
./bin/genion --version  # Should show: 1.2.3-dirty
```

## Usage

### Main Pipeline Script

The pipeline now includes a unified main script with YAML configuration support:

```bash
# Edit the provided configuration file with your paths and settings
# examples/mouse_config.yaml is provided as a template

# Run the pipeline with dry-run to see what would be executed
python typhon_main.py --dry-run

# Run the full pipeline
python typhon_main.py

# Run specific modules only
python typhon_main.py --modules longgf

# Use a custom configuration file
python typhon_main.py --config my_custom_config.yaml

# Override configuration settings
python typhon_main.py --threads 30 --output /custom/path
```

### Configuration File Format

The pipeline uses YAML configuration files for easy parameter management:

```yaml
project:
  name: "My_RNA_Fusion_Analysis"
  output_dir: "/path/to/output"
  threads: 20

input:
  fastq_dir: "/path/to/fastq/files"

references:
  genome: "/path/to/genome.fa"
  gtf: "/path/to/annotation.gtf"
  transcriptome: "/path/to/transcripts.fa"

modules:
  longgf:
    enabled: true
    keep_intermediate: false
  genion:
    enabled: true
    min_support: 1
  jaffal:
    enabled: false

options:
  cleanup_intermediate: true
  debug: false
```

### JaffaL Setup (Optional)

If you plan to use the JaffaL module for fusion detection:

1. **Configure reference files** in your config.yaml:
   ```yaml
   modules:
     jaffal:
       enabled: true
       jaffal_dir: ./jaffal/JAFFA-version-2.3
       genome_build: mm39  # or hg38
       annotation: gencode_M28  # or gencode43
       reference_files:
         genome_fasta_gz: /path/to/genome.fa.gz
         transcriptome_fasta: /path/to/transcripts.fasta
         annotation_bed: /path/to/annotation.bed
         annotation_tab: /path/to/annotation.tab
   ```

2. **Run the JaffaL setup script**:
   ```bash
   python setup_jaffal.py --debug  # --debug for detailed logging
   ```

The setup script will:
- Download and install JaffaL (JAFFA-version-2.3)
- Configure it for your specified genome build (mouse/human)
- Process reference files and build Bowtie2 indices
- Apply TYPHON-specific modifications

### Current Status

- **LongGF Integration**: ✅ Fully implemented and tested with production data
- **Genion Integration**: ✅ Fully implemented and tested - produces identical results to original pipeline
- **JaffaL Integration**: 🚧 Setup script implemented, execution module in progress
- **Configuration System**: ✅ YAML-based configuration with validation and example templates

## Contributing

*Coming soon - contribution guidelines will be added*

## License

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License - see the [LICENSE](LICENSE) file for details.

**Academic and Research Use Only** - Commercial use is prohibited.

## Citation

When using Typhon, please cite the original tools and the upcoming Typhon paper:

### Typhon Pipeline
- **Typhon**: Ada, E., & Kane, H. (2025). Typhon: Modular Pipeline for Chimeric RNA Detection. *[Citation will be provided upon publication]*

### Integrated Tools
- **LongGF**: Liu, Q., et al. (2020). LongGF: computational algorithm and software tool for fast and accurate detection of gene fusions by long-read transcriptome sequencing. *BMC Genomics* 21:793. https://doi.org/10.1186/s12864-020-07207-4

- **Genion**: Karaoglanoglu, F., et al. (2022). Genion, an accurate tool to detect gene fusion from long transcriptomics reads. *BMC Genomics* 23:129. https://doi.org/10.1186/s12864-022-08339-5

- **JAFFA**: Davidson, N.M., et al. (2015). JAFFA: High sensitivity transcriptome-focused fusion gene detection. *Genome Medicine* 7:43. https://doi.org/10.1186/s13073-015-0167-x

- **JAFFAL**: Davidson, N.M., et al. (2022). JAFFAL: detecting fusion genes with long-read transcriptome sequencing. *Genome Biology* 23:10. https://doi.org/10.1186/s13059-021-02588-5 