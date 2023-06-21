#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/load_config.R")

input_file <- snakemake@input[["file"]]
sp_info <- get_species_info(species)
contrast <- snakemake@params[["input_contrast"]]

preprocessed_data <- input_file %>%
  load_data(
    sp_info,
    id_column_name = id_column_name,
    id_type = id_type
  ) %>%
  add_human_ids(sp_info)

saveRDS(list(contrast = preprocessed_data), snakemake@output)