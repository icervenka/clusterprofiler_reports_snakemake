#!/usr/bin/Rscript
# imports ----------------------------------------------------------------------
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(flexdashboard))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyverse))

source("snakemake/scripts/functions.R")
source("snakemake/scripts/wrappers.R")

# snakemake parameters ---------------------------------------------------------
source("snakemake/scripts/load_params.R")

# snakemake inputs -------------------------------------------------------------
input_file = snakemake@input[['file']]
input_type = snakemake@input[['type']]

# load data and metadata -------------------------------------------------------
meta = read.table(metadata, header = T, stringsAsFactors = F)

data = input_file %>%
  load_data(get_species_info(species), 
            id_column_name = id_column_name, 
            id_type = id_type)

all_data = paste0(snakemake@params[['input_dir']], meta$file) %>%
  map( ~ load_data(
    .x,
    get_species_info(species),
    id_column_name = id_column_name,
    id_type = id_type
  )) %>%
  setNames(meta$contrast)

knitr_output_options = list(
  mathjax = NULL,
  self_contained = self_contained,
  lib_dir = paste0("../../", report_outdir, input_contrast, "/_libs")
)

# run analysis -----------------------------------------------------------------
template = "contrast"
payloads = merge(read.csv(input_type, stringsAsFactors = F),
                 template_to_df(get_template(template)))

run_cp(data,
       get_species_info(species),
       payloads,
       template,
       input_contrast,
       output_opts = knitr_output_options)

# run analysis upset unique set ------------------------------------------------
# only run if there is more than one contrast
if(length(meta$contrast %>% unique) > 1) {
  template = "unique"
  payloads = merge(read.csv(input_type, stringsAsFactors = F),
                   template_to_df(get_template(template)))
  
  unique_data = get_unique_expr_data(all_data, input_contrast, min_set_size = min_set_size)
  run_cp(
    unique_data,
    get_species_info(species),
    payloads,
    template,
    input_contrast,
    output_opts = knitr_output_options
  )
}
