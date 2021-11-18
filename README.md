# ClusterProfiler report generator (Snakemake pipeline)

This program serves as batch processing Snakemake pipeline for pathway
enrichment of bulk RNASeq data. It  makes use of clusterProfiler and flexdashboard
markdown with widgets to generate reports on pathway analysis for gene expression data.
Pathways analysis is split into several standalone html reports. Two types of pathway
analysis are offered under clusterProfiler package: Over-representation analysis
(ORA) and Gene Set Enrichment Analysis (GSEA).

## Usage
You will need a `metadata.tsv` file with two columns `contrast` and `file` containing
the name of the contrast from RNASeq data to analyse the file corresponding file from gene
differential expression analysis containing data for this particular contrast. The
required structure of the file is described in the next part.

Next the `config.yaml` file needs to be edited to reflect your desired output and
analysis parameters. Please see the parameter section of the readme file for detailed
description.

Run the snakemake pipeline by navigating to the folder with `Snakefile` and executing
it by command `snakemake`. One can optionally perform a test dry-run to see if pipeline
was correcly configured by running `snakemake -n`

Please see the documentation to [Snakemake pipeline language](https://snakemake.github.io/)
for detailed introduction how to run the program.

## Input
A text based file containing the output of a gene differential expression analysis must contain the following columns:
  - gene identifier, must be one of supported by org.db
  - log2 fold-change expression values
  - corresponding p-values or adjusted p-values 

## Output
Program creates an output directory with three subdirectories.
- `reports` contains the contrast subdirectories with html report from individual
specified pathway collections
- `csv` contains the contrast subdirectories with a csv file for each combination of
pathway collection, type of analysis (ORA or GSEA) and gene list (all differentially
  expressed gene, only upregulated or only downregulated genes)
- `analysis_params` contains copies of `metadata.tsv` and `config.yaml` files

Finally the program creates a `pathway_archive.tar.gz` file containg the output directory

## Configuration and Paramters

#### species
Species of the samples the RNASeq data comes from. Currently supported species are
'human', 'mouse' and partially 'rat'. Some of the pathway collections only accept human gene IDs,
so they need to by translated using homology databases.

#### create_archive
TRUE/FALSE, If `tar.gz` archive should be created at the end.

#### self_contained
TRUE/FALSE, should the html output be self contained or have a separate lib directory
separate lib directory saves space, but individual html files can't be shared with other
people.

#### types
Pathway collections to be used in pathway analysis, currently supported are subsets of:
- Disease
- GeneOntology
- MSigDB
- MSigDB_cancer
- MeSH_g2p
- MeSH_gendoo
- Pathways
- ConsensusPathDB

For detailed subcategories included in the default pathway analysis see either
`snakemake/data/csv` files or `snakemake/utilities/export_structures.R`


### Input description

#### id_column_name
Name of the gene ID column in the differential expression file

#### id_type
Type of ID column in the deg file. Can be one of supported by org.db, usually
ENTREZID, ENSEMBL or SYMBOL

#### fc_column_name
Name of the log2-fold-change column in the differential expression file

#### padj_column_name
Name of the p-value column in the differential expression file. Used for
volcano plots.

#### pval_column_name
Name of the adjusted p-value column in the differential expression file


### Analysis parameters

#### gene_fc_threshold
Gene log2 Fold-Change differential expression threshold for ORA.
Program will only consider genes above the specified threshold. Usually set to
0.58, which roughly corresponds to 1.5 times increase in expression

#### gene_padj_threshold
Threshold for adjusted p-value of differentially expressed genes for ORA. Only
genes with value lower than threshold will be considered for analysis. Usually
set to 0.05.

#### enrich_minGS
Minimal number of genes in pathway for ORA. Pathways with lower amount of genes
will be omitted from analysis. Usually set to 5-10.

#### enrich_maxGS
Maximum number of genes in pathway for ORA. Pathways with higher amount of genes
will be omitted from analysis. Usually set to around 500.

#### enrich_pval_cutoff
Minimal p-value threshold for over-represented pathways analyzed by ORA.

#### go_simplify_cutoff
Cutoff value to simplify redundant terms from Gene Ontology analysis based on
semantic similarity. See [clusterProfiler vignette](https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html?#reduce-redundancy-of-enriched-go-terms) for more detailed explanation.

#### gse_minGS
Minimal number of genes in pathway for GSEA. Pathways with lower amount of genes
will be omitted from analysis. Usually set to 5-10.

#### gse_maxGS
Maximum number of genes in pathway for GSEA. Pathways with higher amount of genes
will be omitted from analysis. Usually set to around 500.

#### gse_pval_cutoff
Minimal p-value threshold for over-represented pathways analyzed by GSEA.

#### min_set_size
In case of multiple contrasts, clusterProfiler will analyse them seperately as
well as their unions, differences and intersections. Min_set_size deteremines the
minimum amount of genes such subset must have in order to appear in the summary
Upset plot.

### Graphing parameters
clusterProfiler has several visualisations availabe for enrichment results, please
see the [vignette](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)
for detailed information.


#### dotplot: categories
Maximum number of pathways/rows shown in the dotplot.

#### networkplot: categories
Maximum number of pathway nodes shown in the network plot. All genes in the pathway
will be shown.

#### heatplot: categories
Maximum number of pathways/rows shown in the heatmap plot.

#### enrichmap: categories
Maximum number of pathway nodes shown in the plot depicting pathway enrichment and
similarities through common genes.

#### table: hide_columns
Hide columns in the summary table.


## Utilities
There are several utilities included with the program in the `snakemake/utilities` folder
to aid either pre- or post-processing.
- `export_structures.R` - the program makes use of csv files containing the descriptions
of pathway collections to process. All supported csv files have been included, but if you
want to create your own subsets of pathway collections you can create a structured R list
and export it into a set of csv files.
- `uniprot_keywords_db.R` - I was unable to locate a ready-to-use database containing uniprot
keywords, this scripts builds the simple database using organism constrained gene list with
entrez IDs and keyword IDs matched with controlled vocabulary of uniprot keywords. This mapping
does not preserve hierarchical structure, but it's sufficient for basic analysis. Please see the files
for links how to obtain the source data.
- `pathways2heatmaps.R` - After the pathway analysis, you can plot expression heatmaps
for selected pathways if you provide an expression file with 'symbol' column, metadata
file with 'sample' and 'group' column, folder with pathway csv files and text file
with pathway IDs separated by a newline. There are some additional parameters that can be specified,
please see the `pathways2heatmaps.R --help` for description.

## Issues
This program is still under development with new features and pathway collections being added.
Clustercompare feature comparing pathway analysis between different contrasts coming soon!

## Dependencies
