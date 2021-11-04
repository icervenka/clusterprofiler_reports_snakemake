library(readr)
library(pheatmap)
library(purrr)
library(dplyr)

metadata_file = "~/rnaseq/igor_hectd1_mko/metadata.tsv"
sample_expression_file = "~/rnaseq/igor_hectd1_mko/diffexp/deseq_default/degfiles/sample_expression.csv"
pathways_folder = "~/programming/clusterprofiler_reports_snakemake/output/csv/ko_wt/"
file_pattern = ""
pathway_ids_file = "~/temp/pathway_ids.txt"
output_path = "~/temp/"

metadata = read_delim(metadata_file,
                      delim = "\t", 
                      escape_double = FALSE,
                      trim_ws = TRUE) %>%
  dplyr::select(sample, group) %>%
  unique

sample_expression = read_delim(sample_expression_file,
                               delim = "\t", 
                               escape_double = FALSE,
                               trim_ws = TRUE) %>%
  dplyr::mutate(SYMBOL = toupper(SYMBOL))

pathway_files = list.files(pathways_folder, file_pattern)
pathways = purrr::map_dfr(pathway_files, function(x) {
  full_path = paste0(pathways_folder, x)
  readr::read_delim(full_path, 
                    delim = "\t", 
                    escape_double = FALSE, 
                    trim_ws = TRUE) %>%
    dplyr::mutate(ID = as.character(ID)) %>%
    dplyr::mutate(SYMBOL = toupper(SYMBOL),
                  analysis = gsub(".txt", "", x))
    
})

pathway_ids = read_delim(pathway_ids_file,
                         delim = "\t", 
                         col_names = c("id"),
                         escape_double = FALSE,
                         trim_ws = TRUE)

# pathways = read_delim("~/programming/clusterprofiler_reports_snakemake/output/csv/ko_wt/ConsensusPathDB_HumanCyc_GSEA_all_contrast.txt",
#                       delim = "\t", escape_double = FALSE,
#                       trim_ws = TRUE) %>%
#   dplyr::mutate(SYMBOL = toupper(SYMBOL))
# 
# sample_expression = read_delim("~/rnaseq/igor_hectd1_mko/diffexp/deseq_default/degfiles/sample_expression.csv",
#                                delim = "\t", escape_double = FALSE,
#                                trim_ws = TRUE) %>%
#   dplyr::mutate(SYMBOL = toupper(SYMBOL))

# TODO add annotation
walk(pathway_ids$id, function(x, pathways, sample_expression) {
  a = pathways %>%
    filter(ID == x) %>%
    left_join(sample_expression %>% 
                dplyr::select(SYMBOL, matches(metadata$sample)), 
              by = "SYMBOL") %>%
    dplyr::select(-c(NES, GeneRatio)) %>%
    drop_na()
  
  pathway_name = a$Description[1]
  
  walk(a$analysis %>% unique, function(y) {
    b = a %>%
      dplyr::filter(analysis == y) %>% 
      group_by(SYMBOL) %>%
      arrange(-abs(log2FoldChange)) %>%
      slice_head() %>%
      dplyr::select(SYMBOL, matches(metadata$sample)) %>%
      column_to_rownames(var = "SYMBOL")
    
    num_genes = nrow(b)
    num_samples = length(metadata$sample)
    
    p = pheatmap(
      b,
      color = colorRampPalette(rev(brewer.pal(
        n = 7, name = "RdBu")))(100),
      cluster_cols = F,
      scale = "row",
      cellwidth = 10,
      cellheight = 10,
      border_color = "white",
      treeheight_row = 10
    )
    ggsave(paste0(output_path, 
                  gsub("\\s", "_", pathway_name),
                  "_",
                  y,
                  ".png"), 
           p, 
           units = "mm", 
           width = num_samples * 3.5 + 40, 
           height = num_genes * 3.5 + 20)
    
  })
}, pathways = pathways, sample_expression = sample_expression)

