data_basepath <- "resources/pathway_databases/"

cpdb <- read.table(
  paste0(data_basepath, "CPDB_pathways_genes_entrez.tab"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

uniprot <- read.table(
  paste0(data_basepath, "uniprot_keywords_db.tab"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)