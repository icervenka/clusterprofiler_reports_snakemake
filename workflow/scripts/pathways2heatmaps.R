#!/usr/bin/Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

option_list <- list(
  make_option(c("-m", "--metadata"),
    type = "character", default = "metadata.tsv",
    help = "metadata file describing sample names and groups.
              Must contain 'sample' and 'group' column. [default=%default]",
    metavar = "character"
  ),
  make_option(c("--grouporder"),
    type = "character", default = "wt,ko",
    help = "comma separated list of groups in the order they shoul be
              arranged on the heatmap. [default=%default]", metavar = "character"
  ),
  make_option(c("-e", "--expression"),
    type = "character", default = "sample_expression.csv",
    help = "file contating normalized gene expression data for samples in metadata.
              Must contain 'symbol' column with gene names and columns corresponding
              to sample names. [default=%default]",
    metavar = "character"
  ),
  make_option(c("-i", "--input"),
    type = "character", default = ".",
    help = "directory where output csvs from clusterprofiler pipeline are stored.
               [default=%default]",
    metavar = "character"
  ),
  make_option(c("--pattern"),
    type = "character", default = "",
    help = "regular expression used to match filenames in input directory. Only
    matching files will be processed. [default=%default]",
    metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "_",
    help = "output directory where pathway heatmaps will be stored.
    [default=%default]",
    metavar = "character"
  )
)
# pathway_id_file is a text file containing IDs of reported by clusterprofiler pipeline,
# one ID per line
opt_parser <- OptionParser(
  usage = "%prog [options] pathway_id_file",
  option_list = option_list
)
opt <- parse_args(opt_parser, positional_arguments = 1)
# metadata_file = "~/rnaseq/igor_hectd1_mko/metadata.tsv"
# sample_expression_file = "~/rnaseq/igor_hectd1_mko/diffexp/deseq_default/degfiles/sample_expression.csv"
# pathways_folder = "~/programming/clusterprofiler_reports_snakemake/output/csv/ko_wt/"
# file_pattern = ""
# pathway_ids_file = "~/temp/pathway_ids.txt"
# output_path = "~/temp/"

metadata_file <- opt$options$metadata
grouporder <- data.frame(
  group = stringr::str_split(opt$options$grouporder, ",")[[1]],
  stringsAsFactors = F
)
sample_expression_file <- opt$options$expression
pathways_folder <- opt$options$input
file_pattern <- opt$options$pattern
output_path <- opt$options$output
pathway_ids_file <- opt$args

# read metadata file, extract sample and group annotation only
metadata <- readr::read_delim(metadata_file,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
  dplyr::select(sample, group) %>%
  unique()
# ordering the metadata dataframe based on the grouporder specified on commandline
metadata <- grouporder %>%
  left_join(metadata, by = "group")

# read sample expression file, unify the case of 'symbol' column
sample_expression <- readr::read_delim(sample_expression_file,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
) %>%
  dplyr::mutate(SYMBOL = toupper(SYMBOL))

# list the pathway files containing the pattern
# files are read without full paths, we need to add the file names to the
# large data frame later, it's easier to only remove the extensions
pathway_files <- list.files(pathways_folder, file_pattern)

# read in the individual files and bind them together by rows into a large dataframe
# files from enrichr and gsea functions have non-corresponding columns 'GeneRatio' and NES'
# which are filled with NA. This is ignored since we don't use them anywhere, but they
# have to be removed before any drop_na() call
cat("\nMerging pathway files.\n\n")
pathways <- purrr::map_dfr(pathway_files, function(x) {
  # add the full paths in order to read the files
  full_path <- paste0(pathways_folder, x)
  readr::read_delim(full_path,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  ) %>%
    # cast all the IDs to character, since some of the IDs are numeric
    dplyr::mutate(ID = as.character(ID)) %>%
    # unify the case of 'symbol' column
    # add the filenames the data came from to the dataframe (withtout ext)
    # to remember which analysis the pathway comes from (same pathway can appear
    # more than once)
    dplyr::mutate(
      SYMBOL = toupper(SYMBOL),
      analysis = gsub(".txt", "", x)
    )
})

# read the IDs of pathways to be processed
pathway_ids <- readr::read_delim(pathway_ids_file,
  delim = "\t",
  col_names = c("id"),
  escape_double = FALSE,
  trim_ws = TRUE,
  show_col_types = FALSE
)

# TODO add annotation
# walk through all the pathway IDs
cat("Exporting pathway:\n")
purrr::walk(pathway_ids$id, function(x, pathways, sample_expression, metadata) {
  # select one ID and join with the sample expression table that
  # with only symbol and sample columns retained
  # has to be merged on symbol, because some clusterprofiler databases only
  # work in human gene IDs which in case of mouse data have to be translated
  # by homology
  a <- pathways %>%
    dplyr::filter(ID == x) %>%
    dplyr::left_join(
      sample_expression %>%
        dplyr::select(SYMBOL, dplyr::matches(metadata$sample)),
      by = "SYMBOL"
    ) %>%
    # remove the NES and GeneRatio columns that will contain NAs
    dplyr::select(-c(NES, GeneRatio)) %>%
    # remove the genes from the table that don't have corresponding
    # symbols on both sided of merge
    tidyr::drop_na()

  # store the pathway description for later use in filename generation
  pathway_name <- a$Description[1]
  pathway_name <- gsub("\\s", "_", pathway_name)
  pathway_name <- gsub("[\\/#:\"*?<>|.,;&]", "", pathway_name)
  pathway_name <- substr(pathway_name, 1, 50)

  if (is.null(pathway_name)) {
    pathway_name <- x
  }

  cat(pathway_name)
  cat("\n")

  # with the selected pathway stored in a, we iterate through all the analyses
  # that might contain the pathway
  walk(a$analysis %>% unique(), function(y) {
    b <- a %>%
      # select one analyses
      dplyr::filter(analysis == y) %>%
      # due to human-mouse homologue mapping, some genes might appear more than
      # once. We select the one with the highest fold-change and discard the rest,
      # since we need the gene names as row names of data frame which doesn't
      # allow duplicates
      dplyr::group_by(SYMBOL) %>%
      dplyr::arrange(-abs(log2FoldChange)) %>%
      dplyr::slice_head() %>%
      # select only columns with symbol and samples
      dplyr::select(SYMBOL, matches(metadata$sample)) %>%
      tibble::column_to_rownames(var = "SYMBOL")

    # store the number of rows and columns of future heatmap for
    # calculating the dimensions of the image
    num_genes <- nrow(b)
    num_samples <- length(metadata$sample)

    # generate heatmap, no column clustering
    p <- pheatmap::pheatmap(
      b,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(
        n = 7, name = "RdBu"
      )))(100),
      cluster_cols = F,
      scale = "row",
      cellwidth = 10,
      cellheight = 10,
      border_color = "white",
      treeheight_row = 10,
      annotation_col = metadata %>% column_to_rownames("sample"),
      annotation_legend = FALSE,
      show_colnames = FALSE
    )
    # save as image, paste filename from pathway name
    # and original analysis name
    # dimensions are empirically determined
    ggplot2::ggsave(
      paste0(
        output_path,
        pathway_name,
        "_",
        y,
        ".png"
      ), # \\/:"*?<>|
      p,
      units = "mm",
      width = num_samples * 3.5 + 40,
      height = num_genes * 3.5 + 20
    )
  })
},
pathways = pathways,
sample_expression = sample_expression,
metadata = metadata
)
