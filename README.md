# ClusterProfiler report generator (Snakemake pipeline)

This program serves as batch processing Snakemake pipelne for pathway
enrichment of bulk RNASeq data. It uses ClusterProfiler package to iterate through
pathways and samples and outputs html reports with interactive graphs.

## Usage
You will need a `metadata.tsv` file with two columns `contrast` and `file` containing the name of the contrast from RNASeq data to analyse the file corresponding text file containing
differential expression data for this particular contrast. Required columns in
this file are a 'gene identifier' (ensembl, refseq or gene symbol), 'log2-fold-change'
and 'p-value'.

Next the `config.yaml` file needs to be edited to reflect your desired output and
analysis parameters. Please see the parameter section of the readme file for detailed
description.

Run the snakemake pipeline by navigating to the folder with `Snakefile` and executing
it by command `snakemake`. One can optionally perform a test dry-run to see if pipeline
was correcly configured by running `snakemake -n`

Please see the documentation to [Snakemake pipeline language](https://snakemake.github.io/)
for detailed introduction how to run the program.

## Input


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

#### gene_fc_threshold


#### gene_padj_threshold


#### enrich_minGS
#### enrich_maxGS
#### enrich_pval_cutoff
#### go_simplify_cutoff
#### gse_nperm
#### gse_minGS
#### gse_maxGS
#### gse_pval_cutoff
#### min_set_size
#### dotplot: categories
#### networkplot: categories
#### heatplot: categories
#### enrichmap: categories
#### table: hide_columns

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
