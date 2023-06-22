#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/load_config.R")

input_file <- snakemake@input[[1]]
sp_info <- get_species_info(species)
contrast <- snakemake@params[["contrast"]]

preprocessed_data <- input_file %>%
  load_data(
    sp_info,
    id_column_name = id_column_name,
    id_type = id_type
  ) %>%
  add_human_ids(sp_info)

saveRDS(preprocessed_data, snakemake@output[[1]])