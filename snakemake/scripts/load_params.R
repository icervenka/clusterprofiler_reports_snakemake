# global variables -------------------------------------------------------------------------
# scripts_basepath = "snakemaka/scripts"
data_basepath = "snakemake/data/"
cpdb = read.table(paste0(data_basepath, "CPDB_pathways_genes.tab"),
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F)

# snakemake common parameters ------------------------------------------
input_contrast = snakemake@params[['contrast']]
input_dir = snakemake@params[['input_dir']]
report_outdir = snakemake@params['report_outdir']
csv_outdir = snakemake@params['csv_outdir']
metadata = snakemake@params[['metadata']]
species = snakemake@params[['species']]
self_contained = snakemake@params[["self_contained"]]
id_type = snakemake@params[['id_type']]
id_column_name = snakemake@params[['id_column_name']]
fc_column_name = snakemake@params[['fc_column_name']]
pval_column_name = snakemake@params[['pval_column_name']]
padj_column_name = snakemake@params[['padj_column_name']]
gene_fc_threshold = snakemake@params[['gene_fc_threshold']]
gene_padj_threshold = snakemake@params[['gene_padj_threshold']]

enrich_minGS = snakemake@params[['enrich_minGS']]
enrich_maxGS = snakemake@params[['enrich_maxGS']]
enrich_pval_cutoff = snakemake@params[['enrich_pval_cutoff']]
gse_nperm = snakemake@params[['gse_nperm']]
gse_minGS = snakemake@params[['gse_minGS']]
gse_maxGS = snakemake@params[['gse_maxGS']]
gse_pval_cutoff = snakemake@params[['gse_pval_cutoff']]

templates = snakemake@params[['templates']]
min_set_size = snakemake@params[['min_set_size']]

dotplot_categories = snakemake@params[['dotplot_categories']]
networkplot_categories = snakemake@params[['networkplot_categories']]
heatplot_categories = snakemake@params[['heatplot_categories']]
emapplot_categories = snakemake@params[['emapplot_categories']]

hide_columns = snakemake@params[['hide_columns']]
