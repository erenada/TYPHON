# System Dependencies for Typhon Pipeline

This document outlines the system dependencies required for the Typhon pipeline and provides installation instructions for different scenarios.

## Required System Dependencies

### 1. Java 11 (Required for JaffaL)

**Installation:**
```bash
# Ubuntu/Debian
sudo apt install openjdk-11-jre

# Check installation
java -version
```

**Important:** Java 11 is specifically required for JaffaL compatibility. Other versions may cause issues.

### 2. Perl Rename Utility (Required for file processing)

The Typhon pipeline requires the Perl-based `rename` utility for advanced file renaming operations using regular expressions.

#### Option A: Conda Installation (Recommended)
The rename utility is included in the conda environment and will be installed automatically:

```bash
conda env create -f environment.yml
# OR update existing environment
conda env update -n typhon_env -f environment.yml
```

This installs the `rename` package from Bioconda (version 1.601).

#### Option B: System Installation (Alternative)
If you prefer system-wide installation or the conda version doesn't meet your needs:

```bash
# Ubuntu/Debian
sudo apt install rename

# Verify it's the Perl version
rename --version
# Should show: "/usr/bin/rename using File::Rename version X.XX"
```

**Note:** Make sure you get the Perl-based rename utility, not the util-linux version.

#### Differences Between Versions:
- **System version (`apt install rename`)**: Usually File::Rename, Perl-based with `--version` flag
- **Conda version (`bioconda::rename`)**: Version 1.601, different command syntax, uses `--man` for help

Both versions provide Perl regex functionality required by the pipeline.

## Verification Commands

After installation, verify all dependencies:

```bash
# Check Java
java -version

# Check rename utility (conda environment)
conda run -n typhon_env rename --man | head -5

# OR check system rename
which rename && rename --version

# Verify conda environment tools
conda run -n typhon_env which minimap2 longgf samtools
```

## Troubleshooting

### Java Issues
- **Wrong Java version:** Ensure Java 11 is installed and set as default
- **PATH issues:** Make sure Java is accessible from conda environment

### Rename Utility Issues  
- **Command not found:** Install via conda (recommended) or system package manager
- **Wrong version:** Verify you have Perl-based rename, not util-linux rename
- **Permission issues:** System installation may require sudo privileges

### Platform-Specific Notes

#### Linux (Ubuntu/Debian)
- Use `apt install openjdk-11-jre rename` for system dependencies
- Conda environment handles most bioinformatics tools

#### macOS
- Use `brew install openjdk@11` for Java
- Rename utility available via conda or `brew install rename`

#### Windows (WSL recommended)
- Use Windows Subsystem for Linux with Ubuntu
- Follow Linux installation instructions
- Keep all files within WSL filesystem for best performance

## Docker Alternative

For complex environments, consider using the provided Docker container which includes all system dependencies pre-configured.

## Need Help?

If you encounter issues with system dependencies:
1. Check the troubleshooting section above
2. Verify your operating system and package manager
3. Consult the main README.md for additional setup guidance
4. Open an issue on the GitHub repository with your system details 