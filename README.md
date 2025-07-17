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

### Current Status

- **LongGF Integration**: Fully implemented and tested with production data
- **Genion Integration**: Fully implemented and tested - produces identical results to original pipeline
- **JaffaL Integration**: Planned for future implementation
- **Configuration System**: YAML-based configuration with validation and example templates

## Contributing

*Coming soon - contribution guidelines will be added*

## License

*To be determined* 