setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
testdir = gsub("\\.R","", basename(rstudioapi::getActiveDocumentContext()$path))

config_file = "config.yaml"
test_metadata = "ko_ga_vs_wt_ga_all_metadata.tsv"
input_contrast = "ko_wt_ga"
input_type = "ConsensusPathDB"
template = "contrast"

source("setup_test_env.R", local = TRUE)

cpdb <- read.table(
  "../../../resources/pathway_databases/CPDB_pathways_genes_entrez.tab",
  sep = "\t",
  header = T,
  stringsAsFactors = F
)


payloads_sub = payloads_sub <- payloads %>% dplyr::filter(type == "ConsensusPathDB")

create_subpage_data <- function(expr_data, sp_info, params) {
  print("mapping ids")
  if (params$type %in% c("Disease", "ConsensusPathDB") &&
      sp_info["common"] != "human") {
    expr_data <- expr_data %>%
      add_human_ids(sp_info)
    select_entrez_id <- "ENTREZID_HUMAN"
  } else {
    select_entrez_id <- "ENTREZID"
  }
  
  if (params$analysis %in% c("ORA")) {
    expr_data <- filter_data(
      expr_data,
      padj_column = padj_column_name,
      padj_threshold = gene_padj_threshold,
      fc_column = fc_column_name,
      fc_threshold = gene_fc_threshold
    )
  } else if (params$analysis %in% c("GSEA")) {
    # TODO implement custom ordering
  }
  
  geneList <- structure(
    expr_data[, fc_column_name],
    names = as.character(expr_data[, select_entrez_id])
  ) %>%
    sort(decreasing = T)
  geneList <- filter_gene_list(geneList, params$gene_set)
  
  # if (length(filter_gene_list(geneList, params$gene_set)) < 1) {
  #   return(NULL)
  # }
  
  if (length(geneList) < 1) {
    print("no gene list")
    return(NULL)
  }
  
  print(paste0("running analysis", " ", params$type, " ", params$analysis, " ", params$category))
  # Generate enrichResult or gseaResult dataframe to pass as data to subpage
  tryCatch(
    {
      df <- get_wrapper(params$type, params$analysis)(
        geneList,
        params$category,
        sp_info)
    },
    error = function(e) {
      df <- NULL
    }
  )
  
  print("setting readable")
  if (nrow(data.frame(df)) == 0 | is.null(df)) {
    return(NULL)
  }
  if (df@readable == FALSE) {
    df <- setReadable(df, sp_info["orgdb"], keyType = "ENTREZID")
  }
  
  return(list(
    params = params,
    geneList = geneList,
    df = df,
    expr_data = expr_data
  ))
}

all_cp <- map(
  split(payloads_sub, 1:nrow(payloads_sub)) %>% unname() %>% map(~ c(.x)),
  ~ do.call(
    create_subpage_data,
    list(
      expr_data = data,
      sp_info = sp_info,
      params = .x
    )
  )
)
 
w1 = get_wrapper("MSigDB", "ORA")
r1 = map(c("all", "up", "down"), function(x) {
  df = w1(filter_gene_list(geneList, x), "M2@CP", sp_info)
})

df2 = w1(filter_gene_list(geneList, "all"), "HumanCyc")

subpage_df = all_cp[[29]]$df
dotplot_categories = max_categories$dotplot

ddd = subpage_df@result

if (class(subpage_df)[1] == "enrichResult") {
  (subpage_df@result %>%
     mutate(GeneRatio = DOSE::parse_ratio(GeneRatio)) %>%
     arrange(GeneRatio) %>%
     slice_tail(n = dotplot_categories) %>%
     preserve_order_as_factor(Description) %>%
     ggplot(aes(x = GeneRatio, y = Description)) +
     geom_segment(aes(xend = 0, yend = Description), color = "gray50") +
     geom_point(aes(color = p.adjust, size = Count)) +
     scale_size_continuous(range = c(2, 10)) +
     theme_bw(base_size = 12) +
     theme(legend.position = "none") +
     scale_color_viridis() +
     scale_y_discrete(label = function(x) str_trunc(x, 43)) +
     xlab("Gene Ratio") +
     ylab("")
  ) %>%
    ggplotly()
} else if (class(subpage_df)[1] == "gseaResult") {
  (subpage_df@result %>%
     arrange(abs(NES)) %>%
     slice_tail(n = dotplot_categories) %>%
     arrange(NES) %>%
     preserve_order_as_factor(Description) %>%
     ggplot(aes(NES, Description, fill = qvalues)) +
     geom_bar(stat = "identity") +
     theme_bw(base_size = 12) +
     theme(legend.position = "none") +
     scale_fill_viridis() +
     scale_y_discrete(label = function(x) str_trunc(x, 43)) +
     xlab("Normalized Enrichment Score") +
     ylab(NULL)
  ) %>%
    ggplotly()
}

