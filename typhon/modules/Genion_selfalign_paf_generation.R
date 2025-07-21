#!/usr/bin/env Rscript

library(pafr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE) 
path <- (args[1])
ali <- read_paf(path)
Export_path <- (args[2])
datas <- as.data.frame(ali)
suppressWarnings(df1 <- datas %>% separate(qname, c('col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8')))
suppressWarnings(df2 <- datas %>% separate(tname, c('col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8')))
data_for_export <- data.frame(df1$col1,df2$col1)
write.table(data_for_export, file=Export_path, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
q()
