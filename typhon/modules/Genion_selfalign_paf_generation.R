#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pafr)
})

args <- commandArgs(trailingOnly = TRUE)

# Expect two arguments: input PAF path, output TSV path
if (length(args) != 2) {
  stop("Usage: Genion_selfalign_paf_generation.R <selfalign.paf> <selfalign.tsv>")
}

paf_path <- args[1]
out_path <- args[2]

if (!file.exists(paf_path)) {
  stop(paste("PAF file not found:", paf_path))
}

# Read PAF using pafr; PAF is tab-delimited by specification
ali <- read_paf(paf_path)

# If empty, write an empty TSV and exit cleanly
if (is.null(ali) || nrow(ali) == 0) {
  write.table(data.frame(), file = out_path, quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = FALSE)
  quit(save = "no")
}

# Extract base IDs from qname and tname by trimming everything after the first pipe `|` if present
# If no pipe exists, the full name is kept intact.
# Trim after first pipe, then remove trailing Ensembl version suffix like .2
qbase <- sub("\\|.*$", "", ali$qname)
qbase <- sub("\\.[0-9]+$", "", qbase)
tbase <- sub("\\|.*$", "", ali$tname)
tbase <- sub("\\.[0-9]+$", "", tbase)

out <- data.frame(qbase = as.character(qbase),
                  tbase = as.character(tbase),
                  stringsAsFactors = FALSE)

# Drop duplicate pairs while preserving order
out <- unique(out)

# Write two-column TSV without header
write.table(out, file = out_path, quote = FALSE, sep = "\t",
            col.names = FALSE, row.names = FALSE)

quit(save = "no")
