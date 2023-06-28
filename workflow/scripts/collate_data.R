#!/usr/bin/Rscript

# Snakemake puts in the input list both unnamed and named inputs, effectively
# duplicating them, I select only the named ones 
input_list <- snakemake@input[nchar(names(snakemake@input)) > 0]

all_data <- setNames(
    lapply(input_list, readRDS),
    names(input_list)
)

saveRDS(all_data, snakemake@output[[1]])
