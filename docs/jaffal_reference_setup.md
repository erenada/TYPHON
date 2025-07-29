# JaffaL Reference Files Setup Guide

JaffaL requires specific reference files that must be downloaded separately from UCSC databases. This guide provides step-by-step instructions for downloading and organizing these files.

## Overview

You need **four files** that must match your genome build and annotation version:

1. **Genome FASTA** (`.fa.gz`) - Genomic reference sequences
2. **Transcriptome FASTA** (`.fasta`) - Transcript sequences  
3. **Annotation BED** (`.bed`) - Exon coordinates
4. **Annotation TAB** (`.tab`) - Gene/transcript metadata

**CRITICAL:** All files must use the **same genome build** (e.g., mm39, hg38) and **same annotation version** (e.g., GENCODE M28, M31) as your other reference files used by LongGF and Genion.

## Step 1: Download Genome FASTA

Navigate to [UCSC Genome Downloads](https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/):

**Note:** For older genome builds, use the [UCSC Archive](https://hgdownload.soe.ucsc.edu/goldenPath/archive/) instead.

### For Mouse (mm39)
1. Select `Mus_musculus` → `bigZips/`
2. Download `mm39.fa.gz`

### For Human (hg38)
1. Select `Homo_sapiens` → `bigZips/` 
2. Choose the patch version (e.g., `p13/`, `p14/`) that matches your GENCODE version
3. Download `hg38.pXX.fa.gz` (where XX matches your patch version)

**Finding your patch version:**
- Visit [GENCODE releases](https://www.gencodegenes.org/human/releases.html)
- Click "Show all releases"
- Find your GENCODE version in the "GENCODE release" column
- Check the adjacent "Genome assembly version" column (e.g., "GRCh38.p13")
- The "p.XX" number corresponds to the UCSC patch directories

## Step 2: Download Annotation Files

Navigate to [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)

### Initial Setup (same for all three files)
- **Organism:** Human or Mouse
- **Assembly:** Your genome build (e.g., mm39, hg38)  
- **Group:** "Genes and Gene Predictions"
- **Track:** Your GENCODE version (e.g., "ALL GENCODE VM28")
- **Region:** "genome"

### A) Download TAB file
1. **Output format:** "all fields from selected table"
2. **Output filename:** `<genome>_<annotation>.tab`
   - Example: `mm39_gencode_M28.tab`
3. **Output field separator:** "tsv (tab-separated)"
4. **File type:** "plain text"
5. Click **"get output"**

### B) Download FASTA file
1. **Output format:** "sequence"  
2. **Output filename:** `<genome>_<annotation>.fasta`
   - Example: `mm39_gencode_M28.fasta`
3. **File type:** "plain text"
4. Click **"get output"**
5. On next page: Select **"genomic"** → **"submit"**
6. **IMPORTANT:** Uncheck **"Introns"** box
7. Click **"get sequence"**

### C) Download BED file
1. **Output format:** "BED – browser extensible data"
2. **Output filename:** `<genome>_<annotation>.bed`
   - Example: `mm39_gencode_M28.bed`
3. **File type:** "plain text"
4. Click **"get output"**
5. Under "Create one BED record per:" select **"Exons plus 0 bases at each end"**
6. Click **"get BED"**

## Step 3: File Organization

Create a dedicated directory for JaffaL reference files:

```bash
mkdir -p ./references/jaffal
```

**IMPORTANT:** Place all four downloaded files in this directory together. All files must be in the same location for the setup script to work properly.

```
./references/jaffal/
├── mm39.fa.gz                    # Genome FASTA (compressed)
├── mm39_gencode_M28.fasta        # Transcriptome FASTA
├── mm39_gencode_M28.bed          # Annotation BED
└── mm39_gencode_M28.tab          # Annotation TAB
```

## Step 4: Update Configuration

Update your `config.yaml` with the correct paths:

```yaml
jaffal:
  jaffal_dir: ./jaffal/JAFFA-version-2.3
  reference_files:
    genome_fasta_gz: ./references/jaffal/mm39.fa.gz
    transcriptome_fasta: ./references/jaffal/mm39_gencode_M28.fasta
    annotation_bed: ./references/jaffal/mm39_gencode_M28.bed
    annotation_tab: ./references/jaffal/mm39_gencode_M28.tab
```

## Naming Convention Examples

### Mouse (mm39 + GENCODE M28)
- `mm39.fa.gz`
- `mm39_gencode_M28.fasta`
- `mm39_gencode_M28.bed`
- `mm39_gencode_M28.tab`

### Human (hg38.p14 + GENCODE 43)
- `hg38.p14.fa.gz`
- `hg38.p14_gencode43.fasta`
- `hg38.p14_gencode43.bed`
- `hg38.p14_gencode43.tab`

## Troubleshooting

### Version Mismatch Errors
- Ensure all reference files (genome, GTF, transcriptome, JaffaL files) use identical genome builds and annotation versions
- Check that your main references and JaffaL references match exactly

### File Not Found Errors
- Verify all paths in `config.yaml` point to existing files
- Use absolute paths if relative paths cause issues
- Check file permissions

### JaffaL Setup Failures
- Ensure the genome FASTA file remains compressed (`.gz`)
- Ensure annotation files (.fasta, .bed, .tab) are uncompressed
- Verify file integrity after download

### Download Issues
- Large files may take time to download from UCSC
- If downloads fail, try again during off-peak hours
- Verify internet connection stability for large files

## Additional Resources

- [JaffaL Official Wiki](https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-run-jaffa-with-hg19-or-mm10)
- [UCSC Genome Browser Help](https://genome.ucsc.edu/FAQ/)
- [GENCODE Release Information](https://www.gencodegenes.org/)

## Notes

- This process needs to be done once per genome build/annotation combination
- Downloaded files can be reused across multiple TYPHON analyses
- Consider organizing reference files by genome build for easy management
- The setup process is automated by `setup_jaffal.py` once files are downloaded 