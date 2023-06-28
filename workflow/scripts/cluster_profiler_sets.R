#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")
source("workflow/scripts/wrappers.R")

# snakemake parameters ---------------------------------------------------------
source("workflow/scripts/load_config.R")
source("workflow/scripts/load_extradb.R")

# snakemake inputs and params --------------------------------------------------

input_type <- snakemake@input[["type"]]
all_data <- readRDS(snakemake@input[["all_data"]])

output_file <- snakemake@output[[1]]

template <- snakemake@params[["template"]]
param_type <- snakemake@params[["type"]]
rds_outdir <- snakemake@params[["rds_outdir"]]
log_dir <- snakemake@params[["log_dir"]]

# load data and metadata -------------------------------------------------------
meta <- read.table(metadata, header = TRUE, stringsAsFactors = FALSE)

payloads <- merge(
  read.csv(input_type, stringsAsFactors = FALSE),
  template_to_df(get_template(template))
)

# run analysis -----------------------------------------------------------------
upset_sets <- get_upset_sets(all_data, min_set_size = 1)

selected_data <- map(upset_sets$set_name, function(x) {
  result <- get_unique_expr_data(
    all_data,
    x,
    min_set_size = min_set_size
  )
}) %>%
  setNames(upset_sets$set_name)
print(selected_data)

walk(names(selected_data), function(x) {
  print(x)
  dir.create(paste0(rds_outdir, x), showWarnings = FALSE)
  if (nrow(selected_data[[x]]) > 0) {
  run_cp(
    x,
    get_species_info(species),
    payloads
  ) %>%
    export_cp(
      paste0(rds_outdir, x, "/", input_type)
    )
  }
  write.table(
    data.frame(),
    file = paste0(log_dir, param_type, ".done"),
    col.names = FALSE
  )
})
