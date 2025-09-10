#!/usr/bin/env Rscript

# LongGF_process_results.R
# Purpose: Aggregate per-sample LongGF *_results.txt files into combined summary tables.
# Inputs:  Single argument <output_dir> that contains the LongGF outputs.
# Outputs: Two Excel files and CSV counterparts written into <output_dir>:
#          - Combined_LongGF_chimera_results_with_sample_info.(xlsx|csv)
#          - Combined_LongGF_chimera_results_total.(xlsx|csv)
# Note:    Functionality and output schema remain unchanged. Edits add
#          comments, argument checks, and minor defensive improvements.

library(data.table)
library(magrittr)
library(openxlsx)
library(stringr)
library(tidyr)
library(dplyr)

# Ensure data frames are not implicitly converted to factors
options(stringsAsFactors = FALSE)

# Basic usage/argument validation
args <- commandArgs(trailingOnly = TRUE) 
if (length(args) != 1) {
  cat("Usage: LongGF_process_results.R <output_dir>\n")
  quit(status = 1)
}
path <- (args[1])

# Output directory preparation
folderOutput <- paste0(path, "/")
excel_file_output <- paste0(folderOutput, "Combined_LongGF_chimera_results_with_sample_info.xlsx")
excel_file_output2 <- paste0(folderOutput, "Combined_LongGF_chimera_results_total.xlsx")

# Discover per-sample results files produced earlier in the pipeline
file_names <- list.files(path=path, pattern="*_results.txt", full.names=TRUE, recursive=FALSE)

if (length(file_names) == 0) {
  cat("No *_results.txt files found in", path, "\n")
  quit(status = 1)
}

# Read and stack all results into one table
df <- lapply(file_names, function(f) fread(f, header=FALSE, fill=TRUE, showProgress = FALSE)) %>% 
  set_names(basename(file_names)) %>% 
  rbindlist(idcol = "file_name", fill=TRUE)

# Expand comma-separated entries in the first column into multiple rows
df <- as.data.frame(df)
df = separate_rows(df,1,sep = ",")

# First three columns are: Sample, Chimera_ID, Num_reads_total_chimera
df2 <- df[, -(1:3)]
df3 <- df[, (1:3)]

# Label alternating x./y. columns to match the expected melt pattern
col_odd <- seq_len(ncol(df2)) %% 2
suppressWarnings(colnames(df2)[col_odd == 1] <- paste("x", 1:1000000, sep = "."))
suppressWarnings(colnames(df2)[col_odd == 0] <- paste("y", 1:1000000, sep = "."))
df4 <- as.data.table(cbind(df3, df2))

# Melt paired x./y. columns into two value columns
# (value1/value2 correspond to the A/B sides used below)
df <- as.data.frame(melt.data.table(df4, id.vars = NULL, measure.vars = patterns("x.", "y.")))
# Remove empty rows before parsing
df <- df[!(is.na(df$value1) | df$value1==""), ]
# Restore original schema
colnames(df) <- c("Sample", "Chimera_ID", "Num_reads_total_chimera", "Num_reads_this_read_ID", "To_be_split_A", "To_be_split_B")

# Parse A-side composite string to extract genomic fields
# Expected pattern resembles: (...:.../...:...)-(positions)(strand)
df[c('split1', 'split2')] <- str_split_fixed(df$To_be_split_A, '\\(', 2)
df <- separate(df, split2, into = c("split3", "split4"), sep = 1)
df[c('split5', 'split6')] <- str_split_fixed(df$split4, ':', 2)
df[c('split7', 'split8')] <- str_split_fixed(df$split6, '/', 2)
df[c('split9', 'split10')] <- str_split_fixed(df$split8, ':', 2)
df[c('split11', 'split12')] <- str_split_fixed(df$split7, '\\-', 2)

# Select and rename columns to the original schema; carry B-side string forward
df <- select(df, c("Sample", 'split9', 'Chimera_ID','Num_reads_total_chimera', 'Num_reads_this_read_ID','split5', 'split11', 'split12', 'split3', 'To_be_split_B'))
colnames(df) <- c("Sample", "Read_ID", "Chimera_ID", "Num_reads_total_chimera", "Num_reads_this_read_ID", "Chromosome_Gene_A", "Position_1_Gene_A", "Position_2_Gene_A", "Strand_Gene_A", "To_be_split_B")

# Parse B-side composite string to extract genomic fields
df[c('split1', 'split2')] <- str_split_fixed(df$To_be_split_B, '\\(', 2)
df <- separate(df, split2, into = c("split3", "split4"), sep = 1)
df[c('split5', 'split6')] <- str_split_fixed(df$split4, ':', 2)
df[c('split7', 'split8')] <- str_split_fixed(df$split6, '/', 2)
df[c('split9', 'split10')] <- str_split_fixed(df$split7, '-', 2)

# Construct outputs with required columns and ordering
df_complicated <- select(df, c("Sample", "Read_ID", "Chimera_ID", "Num_reads_total_chimera", "Chromosome_Gene_A", "Position_1_Gene_A", "Position_2_Gene_A", "Strand_Gene_A", "split5", "split9", "split10", "split3"))
colnames(df_complicated) <- c("Sample", "Read_ID", "Chimera_ID", "Num_reads_total_chimera", "Chromosome_Gene_A", "Position_1_Gene_A", "Position_2_Gene_A", "Strand_Gene_A", "Chromosome_Gene_B", "Position_1_Gene_B", "Position_2_Gene_B", "Strand_Gene_B")
df_complicated <- df_complicated[order(df_complicated$Sample, df_complicated$Chimera_ID, df_complicated$Read_ID),]

df_simple <- select(df, c("Read_ID", "Chimera_ID", "Num_reads_total_chimera", "Chromosome_Gene_A", "Position_1_Gene_A", "Position_2_Gene_A", "Strand_Gene_A", "split5", "split9", "split10", "split3"))
colnames(df_simple) <- c("Read_ID", "Chimera_ID", "Num_reads_total_chimera", "Chromosome_Gene_A", "Position_1_Gene_A", "Position_2_Gene_A", "Strand_Gene_A", "Chromosome_Gene_B", "Position_1_Gene_B", "Position_2_Gene_B", "Strand_Gene_B")
df_simple <- df_simple[order(df_simple$Chimera_ID, df_simple$Read_ID),]

# Write outputs (Excel and CSV)
write.xlsx(df_complicated, file = excel_file_output)
write.xlsx(df_simple, file = excel_file_output2)

# Also create CSV versions for better tool compatibility
csv_file_output <- gsub(".xlsx", ".csv", excel_file_output)
csv_file_output2 <- gsub(".xlsx", ".csv", excel_file_output2)
write.csv(df_complicated, file = csv_file_output, row.names = FALSE)
write.csv(df_simple, file = csv_file_output2, row.names = FALSE)

# Completion messages
cat("Excel files created successfully:\n")
cat(" -", excel_file_output, "\n")
cat(" -", excel_file_output2, "\n")
cat("CSV files created successfully:\n")
cat(" -", csv_file_output, "\n")
cat(" -", csv_file_output2, "\n") 