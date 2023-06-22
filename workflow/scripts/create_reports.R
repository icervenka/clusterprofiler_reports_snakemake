#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/load_config.R")

cp_data <- readRDS(snakemake@input[["cp_data"]])
diffexp_data <- readRDS(snakemake@input[["diffexp_data"]])
all_data <- readRDS(snakemake@input[["all_data"]])
contrast <- snakemake@params[["contrast"]]
report_outdir <- snakemake@params[["report_outdir"]]
sp_info <- get_species_info(species)

render_reports <- function(
    cp_data, diffexp_data, all_data, sp_info, contrast,
    cp_script_path, output_opts = list()) {
  template <- ""

  render(cp_script_path,
    output_file = paste0(x, "_", template, ".html"),
    output_format = "all",
    output_dir = paste0(report_outdir, contrast),
    output_options = output_opts
  )
}

render_reports(
  cp_data,
  diffexp_data,
  all_data,
  sp_info,
  contrast,
  cp_script_path
)
