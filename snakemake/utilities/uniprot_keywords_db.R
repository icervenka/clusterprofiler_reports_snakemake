#!/usr/bin/Rscript

# keywlist controlled vocabulary was downloaded from https://www.uniprot.org/docs/keywlist
# gene_mapping is exported from uniprot search constrained by swissprot and organism
# with added 'Cross-reference (GeneID)' and 'Keyword ID' column

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

option_list = list(
  make_option(c("-k", "--kewylist"), type="character", default="kewylist.txt", 
              help=" [default=%default]", 
              metavar="character"),
  make_option(c("-g", "--gene_mapping"), type="character", default="mus_musculus.tab", 
              help=" [default=%default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="uniprot_keywords_db.csv", 
              help=" [default=%default]", 
              metavar="character")
)

opt_parser = OptionParser(usage = "%prog [options]", option_list=option_list);
opt = parse_args(opt_parser);

keywlist = readLines(opt$keywlist)
keywlist = data.frame(X1 = keywlist[grep('^ID |^AC |^IC ', keywlist)], 
             stringsAsFactors = F) %>%
  dplyr::mutate(X1 = stringr::str_trim(gsub("^AC|ID|IC", "", X1)))

uniprot_keyword_mapping = data.frame(
  id = keywlist$X1[c(FALSE, TRUE)],
  keyword = keywlist$X1[c(TRUE, FALSE)]
) %>%
  dplyr::mutate(keyword = gsub("\\.", "", keyword))

rnor_uniprot <- readr::read_delim(opt$gene_mapping,
                           delim = "\t", 
                           escape_double = FALSE,
                           trim_ws = TRUE) %>%
  tidyr::drop_na %>%
  dplyr::mutate(ENTREZID = gsub("(\\d);.*","\\1",`Cross-reference (GeneID)`))

genes = rnor_uniprot$ENTREZID
ids = stringr::str_split(rnor_uniprot$`Keyword ID`, ";")

uniprot_gene_mapping = purrr::map2_dfr(genes, ids, function(x, y) {
  data.frame(id = stringr::str_trim(y), gene = stringr::str_trim(x))
}) %>%
  dplyr::arrange(id)

uniprot_db = uniprot_gene_mapping %>%
  left_join(uniprot_keyword_mapping, by = "id")

write.csv(uniprot_db, opt$output, quote = FALSE, row.names = FALSE)
