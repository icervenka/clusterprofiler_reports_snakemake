#!/usr/bin/Rscript

# keywlist controlled vocabulary was downloaded from https://www.uniprot.org/docs/keywlist
# gene_mapping is exported from uniprot search constrained by swissprot and organism
# with added 'Cross-reference (GeneID)' and 'Keyword ID' column
# Hierarchical keyword relationships are not preserved in this script, which makes it impossible
# to remove redundancies in the final pathway enrichment

# load packages ----------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# commandline options ----------------------------------------------------------
option_list <- list(
  make_option(c("-k", "--keywlist"),
    type = "character", default = "keywlist.txt",
    help = "Text file with controlled vocabulary for uniprot keywords
              [default=%default]",
    metavar = "character"
  ),
  make_option(c("-s", "--source_names"),
    type = "character", default = "mmusculus",
    help = "Name of the source organism which the swissprot gene-keyword
              mappings come from. Currently uses abbreviated species form
              (hsapiens, mmusculus, ...) [default=%default]",
    metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "uniprot_keywords_db.tab",
    help = "Output file to write the compiled database to [default=%default]",
    metavar = "character"
  )
)

opt_parser <- OptionParser(
  usage = "%prog [options] swissprot_tab_db_files",
  option_list = option_list
)
opt <- parse_args(opt_parser, positional_arguments = T)
# replace opt with shorter variable names
keywlist_path <- opt$options$keywlist
source_names <- str_split(opt$options$source_names, ",")[[1]]
output <- opt$options$output
gene_mapping <- opt$args

# stop if number source names doesn't correspond  to number of loaded database files
if (length(source_names) != length(gene_mapping)) {
  stop("Number of source names is not equal to number of gene mapping files to process.")
}

# main script ------------------------------------------------------------------
# parse the keywlist controlled vocabulary text file
# only extracts ids and names without hierarchical info
# if the controlled vocabulary changes, this might need updating
keywlist <- readLines(keywlist_path)
keywlist <- data.frame(
  X1 = keywlist[grep("^ID |^AC |^IC ", keywlist)],
  stringsAsFactors = F
) %>%
  dplyr::mutate(X1 = stringr::str_trim(gsub("^AC|ID|IC", "", X1)))


# there are alternating lines of ID/IC and AC (names)
# this creates a dataframe in which each IDs are taken from odd lines
# and names from even lines. In case of previous parsing errors this will throw
# an error due different number of rows
uniprot_keyword_mapping <- data.frame(
  id = keywlist$X1[c(FALSE, TRUE)],
  keyword = keywlist$X1[c(TRUE, FALSE)]
) %>%
  dplyr::mutate(keyword = gsub("\\.", "", keyword))

# iterate through all the corresponding source-gene_mapping values
uniprot_gene_mapping <- purrr::map2_dfr(gene_mapping, source_names, function(x, y) {
  gm <- readr::read_delim(x,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  ) %>%
    tidyr::drop_na() %>%
    # this column name is provided when exporting from uniprot
    dplyr::mutate(ENTREZID = gsub("(\\d);.*", "\\1", `Cross-reference (GeneID)`))

  # the keyword column is provided as colon separated list of IDs, it is split
  # into the list.
  genes <- gm$ENTREZID
  ids <- stringr::str_split(gm$`Keyword ID`, ";")

  # the genes and keyword IDs are pasted together into dataframe, where gene ID
  # repeats for each keyword. Source name is added as a column to differentiate
  # individual species
  purrr::map2_dfr(genes, ids, function(x, y) {
    data.frame(id = stringr::str_trim(y), gene = stringr::str_trim(x))
  }) %>%
    dplyr::mutate(source_name = y) %>%
    dplyr::arrange(source_name, id)
})

# keyword names are joined from uniprot_keyword_mapping dataframe
uniprot_db <- uniprot_gene_mapping %>%
  dplyr::left_join(uniprot_keyword_mapping, by = "id") %>%
  dplyr::select(source_name, id, gene, keyword)

# export -----------------------------------------------------------------------
write.table(uniprot_db, output, sep = "\t", quote = FALSE, row.names = FALSE)
