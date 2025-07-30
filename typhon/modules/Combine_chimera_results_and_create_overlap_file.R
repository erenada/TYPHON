#!/usr/bin/env Rscript

library(data.table)
library(magrittr)
library(openxlsx)
library(stringr)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE) 
path_genion_dir <- (args[1])
LongGF_excel_file_path <- (args[2])
JaffaL_results_text_file_path <- (args[3])
# Load in the Genion files
print("Loading in the Genion files")
genion_fail_path <- paste0(path_genion_dir, "/total_failed_genion_reads.txt")
genion_pass_path <- paste0(path_genion_dir, "/total_passed_genion_reads.txt")
Genion_fail <- fread(genion_fail_path, sep = "\t")
Genion_pass <- fread(genion_pass_path, sep = "\t")
Genion_fail_chimeras <- Genion_fail$V2
Genion_pass_chimeras <- Genion_pass$V2
Genion_total <- as.data.frame(c(Genion_fail_chimeras, Genion_pass_chimeras))
colnames(Genion_total) <- "Chimera_ID"
Genion_total$Chimera_ID <- str_replace(Genion_total$Chimera_ID, "::", ":")
# Load in the JaffaL file
print("Loading in the JaffaL file")
JaffaL_results <- fread(JaffaL_results_text_file_path, sep = "\t")
JaffaL_total_chimeras <- as.data.frame(JaffaL_results$fusion_genes)
colnames(JaffaL_total_chimeras) <- "Chimera_ID"
# Load in the LongGF file
print("Loading in the LongGF file")
# Try CSV first (faster), then Excel as fallback
csv_path <- gsub("\\.xlsx$", ".csv", LongGF_excel_file_path)
if (file.exists(csv_path)) {
  print(paste("Using CSV version:", csv_path))
  LongGF_result <- read.csv(csv_path, stringsAsFactors = FALSE)
} else {
  print(paste("Using Excel version:", LongGF_excel_file_path))
  LongGF_result <- read.xlsx(LongGF_excel_file_path)
}
# Identify overlapping chimeras
print("Checking for overlap between the the three sets of chimeras")
Subset_LongGF_A <- subset(LongGF_result, LongGF_result$Chimera_ID %in% Genion_total$Chimera_ID)
Subset_LongGF <- subset(Subset_LongGF_A, Subset_LongGF_A$Chimera_ID %in% JaffaL_total_chimeras$Chimera_ID)
Read_ID_file <- as.data.frame(Subset_LongGF$Read_ID)
Read_Chimera_ID_file <- Subset_LongGF[,c("Read_ID","Chimera_ID")]
print("Producing an excel file of overlapping mRNA chimeras")
excel_file_output <- paste0(path_genion_dir, "/Overlapping_mRNA_chimeras.xlsx")
read_id_text_file_output <- paste0(path_genion_dir, "/Read_ID_Overlap.txt")
read_chimera_id_convert_file_output <- paste0(path_genion_dir, "/Read_and_Chimera_ID_convert.txt")
write.xlsx(Subset_LongGF, file = excel_file_output)
write.table(Read_Chimera_ID_file,  file = read_chimera_id_convert_file_output, quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE)
write.table(Read_ID_file,  file = read_id_text_file_output, quote = FALSE, sep ="\t", row.names = FALSE, col.names = FALSE)
q()
