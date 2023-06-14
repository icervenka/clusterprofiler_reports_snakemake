pkgs <- c(
    "clusterProfiler",
    "rmarkdown",
    "htmltools",
    "flexdashboard",
    "DT",
    "plotly",
    "ComplexHeatmap",
    "networkD3",
    "scales",
    "viridis",
    "RColorBrewer",
    "biomaRt",
    "homologene",
    "DOSE",
    "enrichplot",
    "ReactomePA",
    "msigdbr",
    "rWikiPathways",
    "meshes",
    "AnnotationHub",
    "MeSHDbi",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "ggplot2",
    "tibble",
    "tidyr",
    "stringr",
    "purrr",
    "dplyr"
)

lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})