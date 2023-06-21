#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/wrappers.R")

# snakemake parameters ---------------------------------------------------------
source("workflow/scripts/load_config.R")
source("workflow/scripts/load_extradb.R")

# snakemake inputs -------------------------------------------------------------
input_file <- snakemake@input[["file"]]
input_type <- snakemake@input[["type"]]

input_contrast <- snakemake@params[["contrast"]]
rds_outdir <- snakemake@params[["rds_outdir"]]

# load data and metadata -------------------------------------------------------
diffexp_data <- readRDS(input_file)[[1]]
meta <- read.table(metadata, header = TRUE, stringsAsFactors = FALSE)

# knitr_output_options <- list(
#   mathjax = NULL,
#   self_contained = self_contained,
#   lib_dir = paste0("../../", report_outdir, input_contrast, "/_libs")
# )

# run analysis -----------------------------------------------------------------
template <- "contrast"
payloads <- merge(
  read.csv(input_type, stringsAsFactors = FALSE),
  template_to_df(get_template(template))
)

cp_data_contrast <- run_cp(
  diffexp_data,
  get_species_info(species),
  payloads
)

export_cp(
  cp_data_contrast,
  paste0(rds_outdir, input_contrast, "/"),
  template
)

# run analysis upset unique set ------------------------------------------------
# only run if there is more than one contrast
if (length(meta$contrast %>% unique()) > 1) {
  template <- "unique"
  payloads <- merge(
    read.csv(input_type, stringsAsFactors = FALSE),
    template_to_df(get_template(template))
  )

  unique_data <- get_unique_expr_data(
    all_data,
    input_contrast,
    min_set_size = min_set_size
  )

  cp_data_unique <- run_cp(
    unique_data,
    get_species_info(species),
    payloads
  )

  export_cp(
    cp_data_contrast,
    paste0(rds_outdir, input_contrast, "/"),
    template
  )
}
