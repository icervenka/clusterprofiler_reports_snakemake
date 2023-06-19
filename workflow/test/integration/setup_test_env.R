library(yaml)
# load packages ----------------------------------------------------------------
print(getwd())
source("../../scripts/load_packages.R")
source("../../scripts/functions.R")
source("../../scripts/wrappers.R")

# set parameters ---------------------------------------------------------------
test_config <- read_yaml(paste0("config_files/", config_file))
mapply(assign, names(test_config), test_config, MoreArgs = list(envir = .GlobalEnv))

meta = read.table(
  paste0("input_files/", test_metadata), 
  header = T, 
  stringsAsFactors = F)

input_file <- paste0("input_files/", meta$file[meta$contrast == input_contrast])
sp_info = get_species_info(species)
cp_script_path <- "../../scripts/cp.Rmd"
report_outdir <- paste0(testdir, "/", "reports")
csv_outdir <- paste0(testdir, "/", "csvs")

data <- input_file %>%
  load_data(get_species_info(species),
    id_column_name = id_column_name,
    id_type = id_type
  )

all_data <- all_data <- paste0("input_files/", meta$file) %>%
  map(~ load_data(
    .x,
    get_species_info(species),
    id_column_name = id_column_name,
    id_type = id_type
  )) %>%
  setNames(meta$contrast)

payloads <- merge(
  read.csv(
    paste0("input_types/", input_type, ".csv"),
    stringsAsFactors = F
  ),
  template_to_df(get_template(template))
)

knitr_output_options <- list(
  mathjax = NULL,
  self_contained = self_contained,
  lib_dir = paste0("../../", report_outdir, input_contrast, "/_libs")
)

geneList <- structure(
  data[, fc_column_name],
  names = as.character(data[, "ENTREZID"])
) %>%
  sort(decreasing = T)
