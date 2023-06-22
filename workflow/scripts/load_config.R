cp_script_path <- "workflow/scripts/create_reports.Rmd"

# snakemake common parameters --------------------------------------------------
rmarkdown::find_pandoc(dir = snakemake@config[["pandoc_path"]])

metadata <- snakemake@config[["metadata"]]
species <- snakemake@config[["species"]]
self_contained <- snakemake@config[["self_contained"]]
id_type <- snakemake@config[["id_type"]]
id_column_name <- snakemake@config[["id_column_name"]]
fc_column_name <- snakemake@config[["fc_column_name"]]
pval_column_name <- snakemake@config[["pval_column_name"]]
padj_column_name <- snakemake@config[["padj_column_name"]]
gene_fc_threshold <- snakemake@config[["gene_fc_threshold"]]
gene_padj_threshold <- snakemake@config[["gene_padj_threshold"]]

enrich_minGS <- snakemake@config[["enrich_minGS"]]
enrich_maxGS <- snakemake@config[["enrich_maxGS"]]
enrich_pval_cutoff <- snakemake@config[["enrich_pval_cutoff"]]
go_simplify_cutoff <- snakemake@config[["go_simplify_cutoff"]]
gse_nperm <- snakemake@config[["gse_nperm"]]
gse_minGS <- snakemake@config[["gse_minGS"]]
gse_maxGS <- snakemake@config[["gse_maxGS"]]
gse_pval_cutoff <- snakemake@config[["gse_pval_cutoff"]]

min_set_size <- snakemake@config[["min_set_size"]]

dotplot_categories <- snakemake@config[["max_categories"]][["dotplot"]]
networkplot_categories <- snakemake@config[["max_categories"]][["networkplot"]]
heatplot_categories <- snakemake@config[["max_categories"]][["heatplot"]]
emapplot_categories <- snakemake@config[["max_categories"]][["emapplot"]]

hide_columns <- snakemake@config[["hide_columns"]]
