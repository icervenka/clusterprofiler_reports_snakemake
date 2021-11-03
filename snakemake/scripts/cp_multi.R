#!/usr/bin/Rscript

# imports ----------------------------------------------------------------------------------
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(flexdashboard))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))

source("snakemake/scripts/functions.R")
source("snakemake/scripts/wrappers.R")

# snakemake parameters ---------------------------------------------------------------------
source("snakemake/scripts/load_params.R")

# snakemake inputs --------------------------------------------------------------------------
input_type = snakemake@input[['type']]

# load metadata -----------------------------------------------------------------------------
meta = read.table(metadata, header = T, stringsAsFactors = F)

# run analysis upset unique set -------------------------------------------------------------
template = "multi"

all_data = paste0(input_dir, meta$file) %>%
  map(~ load_data(.x, get_species_info(species), id_column_name = id_column_name, id_type = id_type)) %>%
  setNames(meta$contrast)

upset_sets = get_upset_sets(all_data, min_set_size = min_set_size)

unique_multi_payloads = merge(read.csv(input_type, stringsAsFactors = F),
                               template_to_df(get_template(template)))

map(upset_sets$set_name, function(input_contrast) {
  multi_data = get_unique_expr_data(all_data, input_contrast, min_set_size = min_set_size)

  knitr_output_options = list(
    mathjax = NULL,
    highlight = NULL,
    self_contained = self_contained,
    lib_dir = paste0("../../", report_outdir, input_contrast, "/_libs")
  )

  if(nrow(multi_data) > 0) {
    dir.create(paste0("reports/", input_contrast), showWarnings = F)
    run_cp(multi_data, get_species_info(species), unique_multi_payloads, template, input_contrast, output_opts = knitr_output_options)
  }
})
