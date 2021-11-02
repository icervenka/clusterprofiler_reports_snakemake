#!/usr/bin/Rscript
### List type enumeration of analyses to be performed that can be exported as csv
### files to be read by clusterProfiler report script

export_structures = function(s, path) {
  df = map_dfr(s, function(x) as.data.frame(x, stringsAsFactors = F))
  map(df$type %>% unique, function(x) {
    write.csv(df %>% dplyr::filter(type == x), paste0(path, x,".csv"), quote = F, row.names = F)
  })
}

# next iteration
# TODO add pathfinder
# TODO add pathifier
# TODO pathway topology based analysis (SPIA)
# TODO add WGCNA
# TODO add GAGE
# SeqGSEA
#
# PADOG
#
# GAGE:
#
# GOPlot:
#
# pathfindeR:
#
# Uniprot:
#
# Enrichment Browser:


analysis_structure = list(
  list(
    "type" = "Disease",
    "category" = "DO",
    "heading" = "Disease Ontology",
    "flavour_text" = ""
  ),
  list(
    "type" = "Disease",
    "category" = "NCG",
    "heading" = "Network of Cancer Gene",
    "flavour_text" = ""
  ),
  list(
    "type" = "Disease",
    "category" = "DGN",
    "heading" = "DisGeNET",
    "flavour_text" = ""
  ),
  list(
    "type" = "GeneOntology",
    "category" = "MF",
    "heading" = "Molecular Function",
    "flavour_text" = ""
  ),
  list(
    "type" = "GeneOntology",
    "category" = "CC",
    "heading" = "Cellular Compartment",
    "flavour_text" = ""
  ),
  list(
    "type" = "GeneOntology",
    "category" = "BP",
    "heading" = "Biological Process",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "H",
    "heading" = "Hallmark",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "C2@CGP",
    "heading" = "Chemical and Genetic perturbations",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "C2@CP",
    "heading" = "Canonical pathways",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "C3@TFT:GTRD",
    "heading" = "Transcription Factor Targets",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "C5",
    "heading" = "Gene Ontology",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB",
    "category" = "C7",
    "heading" = "Immunologic",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB_cancer",
    "category" = "C6",
    "heading" = "Oncogenic",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB_cancer",
    "category" = "C4@CGN",
    "heading" = "Cancer Gene Neighborhoods",
    "flavour_text" = ""
  ),
  list(
    "type" = "MSigDB_cancer",
    "category" = "C4@CM",
    "heading" = "Cancer Modules",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_g2p",
    "category" = "C@gene2pubmed",
    "heading" = "Diseases - Gene2Pubmed",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_g2p",
    "category" = "D@gene2pubmed",
    "heading" = "Chemicals - Gene2Pubmed",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_g2p",
    "category" = "G@gene2pubmed",
    "heading" = "Processes - Gene2Pubmed",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_gendoo",
    "category" = "C@gendoo",
    "heading" = "Diseases - Gendoo",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_gendoo",
    "category" = "D@gendoo",
    "heading" = "Chemicals - Gendoo",
    "flavour_text" = ""
  ),
  list(
    "type" = "MeSH_gendoo",
    "category" = "G@gendoo",
    "heading" = "Processes - Gendoo",
    "flavour_text" = ""
  ),
  list(
    "type" = "Pathways",
    "category" = "KEGG",
    "heading" = "KEGG",
    "flavour_text" = ""
  ),
  list(
    "type" = "Pathways",
    "category" = "Reactome",
    "heading" = "Reactome",
    "flavour_text" = ""
  ),
  list(
    "type" = "Pathways",
    "category" = "Wiki",
    "heading" = "Wiki Pathways",
    "flavour_text" = ""
  ),
  # list(
  #   "type" = "Pathways",
  #   "category" = "Uniprot",
  #   "heading" = "Uniprot keywords",
  #   "flavour_text" = ""
  # ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "PharmGKB",
    "heading" = "PharmGKB",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "HumanCyc",
    "heading" = "HumanCyc",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "SMPDB",
    "heading" = "SMPDB",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "Signalink",
    "heading" = "Signalink",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "NetPath",
    "heading" = "NetPath",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "EHMN",
    "heading" = "EHMN",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "INOH",
    "heading" = "INOH",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "BioCarta",
    "heading" = "BioCarta",
    "flavour_text" = ""
  ),
  list(
    "type" = "ConsensusPathDB",
    "category" = "PID",
    "heading" = "PID",
    "flavour_text" = ""
  )
)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("The export path of the structures has to be specified as an argument.", call.=FALSE)
} elif (length(args) > 1) {
  stop("Too many commandline arguments specified.", call.=FALSE)
}

export_structures(analysis_structure, args[1])
