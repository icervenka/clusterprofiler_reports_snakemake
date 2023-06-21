#!/usr/bin/Rscript

input_files <- snakemake@input[["files"]]

all_data <- lapply(input_files, readRDS)
all_data <- do.call(c, all_data)

saveRDS(all_data, snakemake@output)