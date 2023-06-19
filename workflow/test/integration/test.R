setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(yaml)
# load packages ----------------------------------------------------------------
source("../scripts/load_packages.R", local = TRUE)
source("../scripts/functions.R", local = TRUE)
source("../scripts/wrappers.R", local = TRUE)

# set parameters ---------------------------------------------------------------
source("param_files/param_file_1.R", local = TRUE)
test_config = read_yaml("config_files/config.yaml")
mapply(assign, names(test_config), test_config, MoreArgs = list(envir = .GlobalEnv))

report_outdir = "output/reports"
csv_outdir = "output/csvs"

metadata = "input_files/ko_ga_vs_wt_ga_all_metadata.tsv"
meta = read.table(metadata, header = T, stringsAsFactors = F)

input_file = "input_files/ko_ga_vs_wt_ga_all.txt"
input_type = "input_types/MeSH_g2p.csv"




data = input_file %>%
  load_data(get_species_info(species), 
            id_column_name = id_column_name, 
            id_type = id_type)


all_data <- paste0(input_dir, meta$file) %>%
  map(~ load_data(
    .x,
    get_species_info(species),
    id_column_name = id_column_name,
    id_type = id_type
  )) %>%
  setNames(meta$contrast)

template = "contrast"
payloads = merge(read.csv(input_type, stringsAsFactors = F),
                 template_to_df(get_template(template)))



setwd("~/pathway_analysis/shamim_gpr35ko_gastroc_pathways/")
report_outdir = "."
input_contrast = "ko_wt_ga"
output_dir = "."

knitr_output_options = list(
  mathjax = NULL,
  self_contained = F,
  lib_dir = paste0("../../", report_outdir, input_contrast, "/_libs")
)

cp_script_path = "~/programming/clusterprofiler_reports_snakemake/workflow/scripts/cp.Rmd"

run_cp(data,
       get_species_info(species),
       payloads,
       template,
       input_contrast,
       output_opts = knitr_output_options)


select_entrez_id = "ENTREZID"
geneList = structure(data[, fc_column_name], names = as.character(data[, select_entrez_id])) %>%
  sort(decreasing = T)


# Generate enrichResult or gseaResult dataframe to pass as data to subpage
map(c("all", "up", "down"), function(x) {
  df = get_wrapper("MeSH_g2p", "ORA")(filter_gene_list(geneList, x),
                                                       "G@gendoo",
                                                       get_species_info("mouse"))
  print(x)
  print(df)
})

sp_info = get_species_info(species)
#payloads = payloads[8,]
all_cp = map(
  split(payloads, 1:nrow(payloads)) %>% unname %>% map( ~ c(.x)),
  ~ do.call(
    create_subpage_data,
    list(
      expr_data = data,
      sp_info = sp_info,
      params = .x
    )
  )
)

preserve_order_as_factor <- function(data, cols) {
  gc <- enquos(cols)
  data <- data %>%
    mutate(across(!!!gc, function(x) {
      factor(x, levels = unique(x))
    }))
  return(data)
}


df = get_wrapper("GeneOntology", "ORA")(filter_gene_list(geneList, "down"),
                                    "CC",
                                    get_species_info("mouse"))


df_f = df@result %>%
  mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>%
  arrange(GeneRatio) %>% 
  slice_tail(n = 10) %>%
  preserve_order_as_factor(Description)

(ggplot(df_f, aes(x = GeneRatio, y = Description)) +
  geom_segment(aes(xend=0, yend = Description), color = "gray50") +
  geom_point(aes(color=-log10(p.adjust), size = Count)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_color_viridis() +
  xlab("Gene Ratio") +
  ylab("")) %>%
  suppressWarnings(print(ggplotly()))

df2 = get_wrapper("GeneOntology", "GSEA")(geneList,
                                        "CC",
                                        get_species_info("mouse"))

(df2_f = df2@result %>%
  arrange(abs(NES)) %>% 
  slice_tail(n = 10) %>%
  arrange(NES) %>%
  preserve_order_as_factor(Description) %>%

ggplot(aes(NES, Description, fill=qvalues)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_viridis() +
  xlab("Normalized Enrichment Score") +
  ylab(NULL)) %>%
  suppressWarnings(print(ggplotly()))


gseaplot2(df2)

gseaplot2(df2, "GO:0000151", title = "GO:0000151", subplots = 1:2) +
  theme(plot.title=element_text(hjust=0.5))
      
       