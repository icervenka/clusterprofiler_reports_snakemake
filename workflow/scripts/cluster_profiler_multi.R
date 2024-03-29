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

#input_contrast <- snakemake@params[["contrast"]]
rds_outdir <- snakemake@params[["rds_outdir"]]
template <- snakemake@params[["template"]]

# load data and metadata -------------------------------------------------------
diffexp_data <- readRDS(input_file)
all_data <- readRDS(all_data_file)
meta <- read.table(metadata, header = TRUE, stringsAsFactors = FALSE)

upset_sets <- get_upset_sets(all_data, min_set_size = min_set_size)

unique_multi_payloads <- merge(
  read.csv(input_type, stringsAsFactors = FALSE),
  template_to_df(get_template(template))
)

cp_data_multi <- map(upset_sets$set_name, function(input_contrast) {
  multi_data <- get_unique_expr_data(
    all_data,
    input_contrast,
    min_set_size = min_set_size
    )

  if (nrow(multi_data) > 0) {
    # dir.create(paste0("reports/", input_contrast), showWarnings = FALSE)
    run_cp(
      multi_data,
      get_species_info(species),
      unique_multi_payloads
    )
  }
})
