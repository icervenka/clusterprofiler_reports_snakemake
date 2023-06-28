#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/load_config.R")

cp_data_filename <- snakemake@input[["cp_data"]]
cp_data <- readRDS(cp_data_filename)
# output_filename <- paste0(
#   tools::file_path_sans_ext(basename(cp_data_filename)),
#   ".html")
output_filename <- snakemake@output[[1]]

diffexp_data <- readRDS(snakemake@input[["diffexp_data"]])
all_data <- readRDS(snakemake@input[["all_data"]])
contrast <- snakemake@params[["contrast"]]
template <- snakemake@params[["template"]]
report_outdir <- snakemake@params[["report_outdir"]]
sp_info <- get_species_info(species)

knitr_output_options <- list(
  mathjax = NULL,
  self_contained = self_contained,
  lib_dir = paste0("../../", report_outdir, contrast, "/_libs")
)

render_reports(
  cp_data,
  diffexp_data,
  all_data,
  sp_info,
  contrast,
  template,
  report_outdir,
  output_filename,
  cp_script_path,
  knitr_output_options
)
