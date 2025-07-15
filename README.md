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

*Coming soon - pipeline is under development*

The pipeline will support:
- Long-read RNA sequencing data (FASTQ format)
- Multiple reference types (Gencode, Ensembl)
- Comprehensive chimeric RNA detection using three complementary tools
- Automated reference preparation and processing

## Contributing

*Coming soon - contribution guidelines will be added*

## License

*To be determined* 