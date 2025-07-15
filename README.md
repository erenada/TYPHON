# Typhon: Modular Pipeline for Chimeric RNA Detection

**Developers:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 0.1.0  

A Python package designed to provide a robust, modular, and user-friendly pipeline for chimeric RNA detection using tools such as LongGF, Genion, and JaffaL.

## Project Status

🚧 **Under Development** - This project is currently being refactored and improved.

## Directory Structure

```
Typhon_pipeline/
├── typhon/                    # Main Python package
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
├── packaging/                # Packaging resources (conda, docker)
├── bin/                      # Compiled binaries (custom Genion)
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

# Create conda environment (includes rename utility)
conda env create -f environment.yml
conda activate typhon_env

# Install the package in development mode (coming soon)
# pip install -e .
```

## Usage

*Coming soon - package is under development*

## Contributing

*Coming soon - contribution guidelines will be added*

## License

*To be determined* 