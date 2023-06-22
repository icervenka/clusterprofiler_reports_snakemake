#!/usr/bin/Rscript

source("workflow/scripts/load_packages.R")
source("workflow/scripts/functions.R")

cp_data_filename <- snakemake@input[["cp_data"]]
cp_data <- readRDS(cp_data_filename)

dfs_to_export <- map_if(
  cp_data[[1]],
  is.not.empty.list,
  ~ process_df_for_csv_export(.x),
  .else = NULL
)

# export ORA to text file
discard(dfs_to_export, is.null) %>%
  map_dfr(function(x) {
    if (x$analysis %>% unique() == "ORA") {
      return(x)
    }
  }) %>%
  write.table(
      snakemake@output[["ora"]],
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
  )

# export GSEA to text file
discard(dfs_to_export, is.null) %>%
  map_dfr(function(x) {
    if (x$analysis %>% unique() == "GSEA") {
      return(x)
    }
  }) %>%
  write.table(
      snakemake@output[["gsea"]],
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
  )
