preserve_order_as_factor <- function(data, cols) {
  gc <- enquos(cols)
  data <- data %>%
    mutate(across(!!!gc, function(x) {
      factor(x, levels = unique(x))
    }))
  return(data)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
testdir = gsub("\\.R","", basename(rstudioapi::getActiveDocumentContext()$path))

config_file = "config.yaml"
test_metadata = "ko_ga_vs_wt_ga_all_metadata.tsv"
input_contrast = "ko_wt_ga"
input_type = "MeSH_g2p"
template = "contrast"

source("setup_test_env.R", local = TRUE)

# call individual analyses
w1 = get_wrapper("MSigDB", "ORA")
r1 = map(c("all", "up", "down"), function(x) {
  df = w1(filter_gene_list(geneList, x), "C2@CP", sp_info)
})

# # run a single iteration of run_cp
# run_cp(
#   data,
#   get_species_info(species),
#   payloads,
#   template,
#   input_contrast,
#   output_opts = knitr_output_options
# )

# knit a report
render(
  paste0(testdir, ".Rmd"), 
  output_dir = paste0(testdir, "/output/reports/"))
