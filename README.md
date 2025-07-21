# TYPHON: Research-Grade Chimeric RNA Detection Pipeline

**Developers:** Harry Kane, PhD; Eren Ada, PhD  
**Version:** 1.0.0  
**Last Updated:** 07/21/2025

A comprehensive, research-grade bioinformatics pipeline for chimeric RNA detection from long-read RNA sequencing data, integrating LongGF, custom Genion, and JaffaL with enhanced logging and debugging capabilities.

## Directory Structure

```
Typhon_pipeline/
├── typhon/                    # Main pipeline modules
│   ├── modules/              # Pipeline step modules (longgf, genion, jaffal)
│   ├── utils/                # Utility functions
│   └── scripts/              # Supporting scripts and data processing
├── tests/                    # Test suite (unit, integration, e2e)
├── examples/                 # Example configurations and data
│   ├── configs/              # Example YAML configurations
│   └── data/                 # Small example datasets
├── docs/                     # Documentation
├── scripts/                  # Development and setup scripts
├── bin/                      # Compiled binaries (custom Genion)
├── jaffal/                   # JaffaL installation directory
├── Genion_files/             # Custom Genion patch and source files
└── test_data/                # Test datasets and references
```

## Quick Start

### Prerequisites
- Linux (Ubuntu 18.04+)
- Conda/Mamba package manager  
- Java 11+ (for JaffaL)
- 16+ GB RAM (32+ GB recommended)

### Installation

```bash
# 1. Clone and setup environment
git clone https://github.com/erenada/TYPHON.git
cd TYPHON
conda env create -f environment.yml
conda activate typhon_env

# 2. Configure your data paths
cp config_template.yaml config.yaml
# Edit config.yaml with your FASTQ files and reference paths

# 3. Setup custom Genion (required - must run first)
python setup_genion.py

# 4. Setup JaffaL (required for complete pipeline)
python setup_jaffal.py

# 5. Verify installation
./bin/genion --version  # Should show: 1.2.3-dirty

# 6. Run TYPHON pipeline
python typhon_main.py
```

### Setup Workflow

The TYPHON pipeline requires a specific setup order:

1. **Environment Setup**: Create conda environment and activate it
2. **Configuration**: Copy and edit config template with your data paths
3. **Setup Genion**: Install custom Genion binary with TYPHON modifications
4. **Setup JaffaL**: Download, compile, and configure JaffaL with reference processing
5. **Run Pipeline**: Execute the integrated TYPHON pipeline

#### Setup Script Details

**setup_genion.py**:
- Clones official Genion repository
- Applies TYPHON custom modifications (read ID tracking, debug output)
- Compiles with research-grade enhancements
- Installs binary to `./bin/genion`

**setup_jaffal.py**:
- Downloads JAFFA v2.3 and extracts to `./jaffal/`
- Compiles custom C++ tools for fusion detection
- Processes reference files and builds Bowtie2 indices
- Applies TYPHON-specific configuration modifications

**typhon_main.py**:
- Integrates LongGF, custom Genion, and JaffaL modules
- Manages data flow between pipeline stages
- Provides comprehensive logging and error handling
- Supports modular execution and configuration overrides

### Basic Usage

```bash
# After completing setup steps above, run the pipeline:

# Run complete pipeline (debug enabled by default)
python typhon_main.py

# Run specific modules
python typhon_main.py --modules longgf genion

# Run with custom settings and disable debug logging
python typhon_main.py --config my_config.yaml --threads 30 --no-debug

# Test run without execution
python typhon_main.py --dry-run

# Override output directory
python typhon_main.py --output ./custom_results --threads 16
```

## Pipeline Modules

### LongGF
**Direct RNA-seq fusion detection**
- Input: FASTQ files, genome, GTF annotation
- Output: SAM alignments, Excel/CSV fusion results
- Status: Production ready

### Custom Genion  
**Enhanced graph-based fusion detection**
- Input: FASTQ files, SAM alignments, processed references
- Output: Read-level fusion results with individual read IDs
- Enhancement: Research-mode comprehensive reporting (141x more candidates)
- Status: Production ready

### JaffaL
**JAFFA-Long pipeline integration**
- Input: FASTQ files, reference genome/transcriptome
- Output: JaffaL fusion results with overlap analysis
- Status: Production ready

## Enhanced Features

### Comprehensive Logging
- **Command Output Capture**: All external tool output logged with timing
- **Module-Specific Logs**: Separate log files (longgf.log, genion.log, jaffal.log)
- **Error Debugging**: Complete stderr capture for troubleshooting
- **Reproducibility**: Full audit trail for scientific reproducibility

### Robust Path Handling
- **Absolute Path Resolution**: Automatic conversion eliminates common failures
- **Portable Configuration**: Works with both relative and absolute paths
- **Enhanced Validation**: Comprehensive input file checking

### Research-Grade Output
- **Custom Genion**: Debug mode provides all potential fusion candidates
- **Read-Level Detail**: One line per supporting read with unique IDs
- **Complete Transparency**: PASS/FAIL reasoning for all candidates
- **No Hidden Filtering**: Researcher controls final filtering criteria

## Configuration

TYPHON uses YAML configuration files. See `config_template.yaml` for a complete example.

For JaffaL integration, ensure you configure the reference files paths that will be processed during setup:

```yaml
project:
  name: TYPHON_Analysis
  output_dir: ./results
  threads: 20

options:
  debug: true                 # Enable debug logging (can override with --no-debug)

input:
  fastq_dir: ./data/FASTQ

references:
  genome: ./references/genome.fa
  gtf: ./references/annotation.gtf
  transcriptome: ./references/transcripts.fa

modules:
  longgf:
    enabled: true
  genion:
    enabled: true
  jaffal:
    enabled: true
    jaffal_dir: ./jaffal/JAFFA-version-2.3
    genome_build: mm39                    # or hg38 for human
    annotation: gencode_M28               # or gencode43 for human
    min_low_spanning_reads: 1             # Minimum spanning reads threshold
    reference_files:
      genome_fasta_gz: /path/to/genome.fa.gz
      transcriptome_fasta: /path/to/transcripts.fasta
      annotation_bed: /path/to/annotation.bed
      annotation_tab: /path/to/annotation.tab
```

## Documentation

- **Module Documentation**: Detailed docs for each module available separately
- **Genion Customization**: See `Genion_customization.md` for implementation details
- **Development Roadmap**: See `pipeline_dev_roadmap.md` for development status
- **System Dependencies**: See `docs/system_dependencies.md` for detailed installation requirements
- **Integration Testing**: See `docs/genion_integration_improvements.md` for comprehensive testing documentation

## Current Development Status

- **LongGF Integration**: Fully implemented and tested with production data
- **Genion Integration**: Fully implemented and tested - produces identical results to original pipeline  
- **JaffaL Integration**: Setup script implemented, execution module fully functional
- **Configuration System**: YAML-based configuration with validation and example templates
- **Testing Suite**: Comprehensive integration testing completed

## Scientific Impact

TYPHON transforms standard fusion detection from clinical-grade conservative filtering to research-grade comprehensive analysis:

- **LongGF**: Validated identical results (scientific accuracy maintained)
- **Genion**: 141x more fusion candidates (4,816 vs 34) with complete reasoning
- **JaffaL**: 99.6% consistency with enhanced robustness
- **Integration**: All three modules working harmoniously with full logging

## Troubleshooting

### Common Issues
- **Setup order**: Ensure you run `setup_genion.py` before `setup_jaffal.py`
- **Module errors**: Check module-specific log files in `logs/` directory
- **Path issues**: Use absolute paths in configuration
- **Missing tools**: Ensure conda environment is activated
- **Setup failures**: Run setup scripts with `--debug` flag for detailed logging
- **Binary not found**: Verify `./bin/genion --version` shows `1.2.3-dirty` after setup

### Getting Help
1. Debug logging is enabled by default (disable with `--no-debug` if needed)
2. Use dry-run mode to test configuration: `python typhon_main.py --dry-run`
3. Check log files in your output directory's `logs/` folder
4. Verify all paths in config.yaml exist and are accessible
5. Ensure conda environment activation: `conda activate typhon_env`

### Command-Line Options
- `--config`, `-c`: Specify configuration file (default: config.yaml)
- `--threads`, `-t`: Override thread count from configuration
- `--output`, `-o`: Override output directory from configuration  
- `--modules`: Run specific modules only (longgf, genion, jaffal)
- `--debug`: Enable debug logging (enabled by default)
- `--no-debug`: Disable debug logging
- `--dry-run`: Show execution plan without running



## Citation

When using TYPHON, please cite:

**TYPHON Pipeline**: Citation will be updated upon publication

**Integrated Tools**:
- **LongGF**: Liu, Q., et al. (2020). LongGF: computational algorithm and software tool for fast and accurate detection of gene fusions by long-read transcriptome sequencing. *BMC Genomics* 21:793. https://doi.org/10.1186/s12864-020-07207-4

- **Genion**: Karaoglanoglu, F., et al. (2022). Genion, an accurate tool to detect gene fusion from long transcriptomics reads. *BMC Genomics* 23:129. https://doi.org/10.1186/s12864-022-08339-5

- **JAFFA**: Davidson, N.M., et al. (2015). JAFFA: High sensitivity transcriptome-focused fusion gene detection. *Genome Medicine* 7:43. https://doi.org/10.1186/s13073-015-0167-x

- **JAFFAL**: Davidson, N.M., et al. (2022). JAFFAL: detecting fusion genes with long-read transcriptome sequencing. *Genome Biology* 23:10. https://doi.org/10.1186/s13059-021-02588-5

## License

Creative Commons Attribution-NonCommercial 4.0 International License  
**Academic and Research Use Only** - Commercial use prohibited 