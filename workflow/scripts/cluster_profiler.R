#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/wrappers.R")

# snakemake parameters ---------------------------------------------------------
source("workflow/scripts/load_config.R")
source("workflow/scripts/load_extradb.R")

# snakemake inputs and params --------------------------------------------------
input_file <- snakemake@input[["file"]]
input_type <- snakemake@input[["type"]]

input_contrast <- snakemake@params[["contrast"]]
template <- snakemake@params[["template"]]
rds_outdir <- snakemake@params[["rds_outdir"]]

# load data and metadata -------------------------------------------------------
diffexp_data <- readRDS(input_file)
meta <- read.table(metadata, header = TRUE, stringsAsFactors = FALSE)

payloads <- merge(
  read.csv(input_type, stringsAsFactors = FALSE),
  template_to_df(get_template(template))
)

# run analysis -----------------------------------------------------------------
if (template == "all") {
  run_cp(
    diffexp_data,
    get_species_info(species),
    payloads
  ) %>%
    export_cp(
      snakemake@output[[1]]
    )
} else if (template == "unique") {
  all_data <- readRDS(snakemake@input[["all_data"]])

  unique_data <- get_unique_expr_data(
    all_data,
    input_contrast,
    min_set_size = min_set_size
  )

  run_cp(
    unique_data,
    get_species_info(species),
    payloads
  ) %>%
    export_cp(
      snakemake@output[[1]]
    )
}