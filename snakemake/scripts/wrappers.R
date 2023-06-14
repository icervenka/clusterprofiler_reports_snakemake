# different clusterProfiler enrichr and gsea functions take different number of
# arguments. The wrappers present the same interface, so they are more easily
# used in loops.

# analysis functions  ----------------------------------------------------------
get_wrapper <- function(type, analysis) {
  function_map <- list(
    "Disease" = list(
      "ORA" = enrichDisease_wrapper,
      "GSEA" = gseaDisease_wrapper
    ),
    "GeneOntology" = list(
      "ORA" = enrichGO_wrapper,
      "GSEA" = gseaGO_wrapper
    ),
    "MSigDB" = list(
      "ORA" = enrichMSIG_wrapper,
      "GSEA" = gseaMSIG_wrapper
    ),
    "MSigDB_cancer" = list(
      "ORA" = enrichMSIG_wrapper,
      "GSEA" = gseaMSIG_wrapper
    ),
    "MeSH_g2p" = list(
      "ORA" = enrichMESH_wrapper,
      "GSEA" = gseaMESH_wrapper
    ),
    "MeSH_gendoo" = list(
      "ORA" = enrichMESH_wrapper,
      "GSEA" = gseaMESH_wrapper
    ),
    "Pathways" = list(
      "ORA" = enrichPathways_wrapper,
      "GSEA" = gseaPathways_wrapper
    ),
    "ConsensusPathDB" = list(
      "ORA" = enrichConsensusdb_wrapper,
      "GSEA" = gseaConsensusdb_wrapper
    )
  )
  return(function_map[[type]][[analysis]])
}

# wrapper functions  -----------------------------------------------------------
enrichDisease_wrapper <- function(list, category, species) {
  enrichDO_wrapper <- function(list, category, species) {
    df <- enrichDO(
      gene = names(list),
      ont = "DO",
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }

  enrichNCG_wrapper <- function(list, category, species) {
    df <- enrichNCG(
      names(list),
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }

  enrichDGN_wrapper <- function(list, category, species) {
    df <- enrichDGN(
      names(list),
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }

  function_map <- list(
    "DO" = enrichDO_wrapper,
    "NCG" = enrichNCG_wrapper,
    "DGN" = enrichDGN_wrapper
  )
  function_map[[category]](list, category, species)
}

gseaDisease_wrapper <- function(list, category, species) {
  gseaDO_wrapper <- function(list, category, species) {
    df <- gseDO(
      list,
      #  nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS,
      verbose = FALSE
    )
    if (!is.null(df)) {
      df <- setReadable(df, "org.Hs.eg.db", keyType = "ENTREZID")
    }
  }

  gseaNCG_wrapper <- function(list, category, species) {
    df <- gseNCG(list,
      # nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS,
      verbose = FALSE
    )
    if (!is.null(df)) {
      df <- setReadable(df, "org.Hs.eg.db", keyType = "ENTREZID")
    }
  }

  gseaDGN_wrapper <- function(list, category, species) {
    df <- gseDGN(list,
      # nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS,
      verbose = FALSE
    )
    if (!is.null(df)) {
      df <- setReadable(df, "org.Hs.eg.db", keyType = "ENTREZID")
    }
  }

  function_map <- list(
    "DO" = gseaDO_wrapper,
    "NCG" = gseaNCG_wrapper,
    "DGN" = gseaDGN_wrapper
  )
  function_map[[category]](list, category, species)
}

enrichGO_wrapper <- function(list, category, species) {
  df <- enrichGO(
    gene = names(list),
    ont = category,
    OrgDb = eval(parse(text = species[["orgdb"]])),
    pvalueCutoff = enrich_pval_cutoff,
    pAdjustMethod = "BH",
    minGSSize = enrich_minGS,
    maxGSSize = enrich_maxGS,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  if (!is.null(df)) {
    df <- clusterProfiler::simplify(df,
      cutoff = go_simplify_cutoff,
      by = "p.adjust",
      select_fun = min
    )
  }
}

gseaGO_wrapper <- function(list, category, species) {
  df <- gseGO(
    geneList = list,
    OrgDb = eval(parse(text = species[["orgdb"]])),
    ont = category,
    # nPerm = gse_nperm,
    pvalueCutoff = gse_pval_cutoff,
    minGSSize = gse_minGS,
    maxGSSize = gse_maxGS,
    verbose = FALSE
  )
  if (!is.null(df)) {
    df <- clusterProfiler::simplify(df,
      cutoff = go_simplify_cutoff,
      by = "p.adjust",
      select_fun = min
    )
  }
}

enrichMSIG_wrapper <- function(list, category, species) {
  category_split <- unlist(strsplit(category, "@", fixed = T))
  if (length(category_split) == 2) {
    subcategory <- category_split[2]
  } else {
    subcategory <- NULL
  }

  db <- msigdbr(
    species = species[["full_name"]],
    category = category_split[1],
    subcategory = subcategory
  ) %>%
    dplyr::select(gs_name, entrez_gene)

  df <- enricher(
    names(list),
    TERM2GENE = db,
    pvalueCutoff = enrich_pval_cutoff,
    pAdjustMethod = "BH",
    minGSSize = enrich_minGS,
    maxGSSize = enrich_maxGS
  )
}

gseaMSIG_wrapper <- function(list, category, species) {
  category_split <- unlist(strsplit(category, "@", fixed = T))
  if (length(category_split) == 2) {
    subcategory <- category_split[2]
  } else {
    subcategory <- NULL
  }

  db <- msigdbr(
    species = species[["full_name"]],
    category = category_split[1],
    subcategory = subcategory
  ) %>%
    dplyr::select(gs_name, entrez_gene)
  df <- GSEA(
    list,
    TERM2GENE = db,
    # nPerm = gse_nperm,
    pvalueCutoff = gse_pval_cutoff,
    minGSSize = gse_minGS,
    maxGSSize = gse_maxGS
  )
}

enrichMESH_wrapper <- function(list, category, species) {
  category_split <- unlist(strsplit(category, "@", fixed = T))
  df <- enrichMeSH(
    names(list),
    MeSHDb = get_mesh_dbi(species),
    database = category_split[2],
    category = category_split[1],
    pvalueCutoff = enrich_pval_cutoff,
    pAdjustMethod = "BH",
    minGSSize = enrich_minGS,
    maxGSSize = enrich_maxGS
  )
}

gseaMESH_wrapper <- function(list, category, species) {
  category_split <- unlist(strsplit(category, "@", fixed = T))
  df <- gseMeSH(
    list,
    MeSHDb = get_mesh_dbi(species),
    database = category_split[2],
    category = category_split[1],
    # nPerm = gse_nperm,
    pvalueCutoff = gse_pval_cutoff,
    minGSSize = gse_minGS,
    maxGSSize = gse_maxGS
  )
}

enrichPathways_wrapper <- function(list, category, species) {
  enrichKEGG_wrapper <- function(list, category, species) {
    df <- enrichKEGG(
      gene = names(list),
      organism = species[["kegg"]],
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS
    )
  }

  enrichReactome_wrapper <- function(list, category, species) {
    df <- enrichPathway(
      gene = names(list),
      organism = species[["common"]],
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS,
      readable = TRUE
    )
  }

  enrichWiki_wrapper <- function(list, category, species) {
    suppressPackageStartupMessages(library(rWikiPathways))
    wiki_archive <- downloadPathwayArchive(
      organism = species[["full_name"]],
      format = "gmt"
    )
    wp2gene <- read.gmt(wiki_archive)
    if (file.exists(wiki_archive)) {
      file.remove(wiki_archive)
    }
    wp2gene <- wp2gene %>%
      tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME

    df <- enricher(
      names(list),
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name,
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS
    )
  }

  enrichUniprot_wrapper <- function(list, category, species) {
    # up <- UniProt.ws(species['tax'])
    # all_keys = keys(up, "ENTREZ_GENE")
    # res <- select(up, all_keys, columns = c("KEYWORDS"), keytype = "ENTREZ_GENE") %>%
    #   na.omit
    # # res <- select(up, keys = c("22627","22629"), columns = c("KEYWORDS"), keytype = "ENTREZ_GENE") %>%
    # #   na.omit
    # uniprot_db = uniprot_db %>%
    #   separate_rows(KEYWORDS, sep = ";") %>%
    #   group_by(KEYWORDS) %>%
    #   mutate(group_id = match(KEYWORDS, unique(KEYWORDS)))
    #
    # term2gene = uniprot_db %>% dplyr::select(group_id, ENTREZ_GENE)
    # term2name = uniprot_db %>% dplyr::select(group_id, KEYWORDS)

    head(uniprot)
    term2gene <- uniprot %>%
      dplyr::filter(source_name == species[["abbreviation"]]) %>%
      dplyr::select(id, gene)

    term2name <- uniprot %>%
      dplyr::filter(source_name == species[["abbreviation"]]) %>%
      dplyr::select(id, keyword)

    df <- enricher(
      names(list),
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = enrich_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = enrich_minGS,
      maxGSSize = enrich_maxGS
    )
  }

  pathways_function_map <- list(
    "KEGG" = enrichKEGG_wrapper,
    "Reactome" = enrichReactome_wrapper,
    "Wiki" = enrichWiki_wrapper,
    "Uniprot" = enrichUniprot_wrapper
  )
  pathways_function_map[[category]](list, category, species)
}

gseaPathways_wrapper <- function(list, category, species) {
  gseaKEGG_wrapper <- function(list, category, species) {
    df <- gseKEGG(
      geneList = list,
      organism = species[["kegg"]],
      # nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS,
      verbose = FALSE
    )
  }

  gseaReactome_wrapper <- function(list, category, species) {
    df <- gsePathway(
      list,
      organism = species["common"],
      #  nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      pAdjustMethod = "BH",
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS,
      verbose = FALSE
    )
  }

  gseaWiki_wrapper <- function(list, category, species) {
    suppressPackageStartupMessages(library(rWikiPathways))
    wiki_archive <- downloadPathwayArchive(
      organism = species["full_name"],
      format = "gmt"
    )
    wp2gene <- read.gmt(wiki_archive)
    if (file.exists(wiki_archive)) {
      file.remove(wiki_archive)
    }
    wp2gene <- wp2gene %>% 
      tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME

    df <- GSEA(
      list,
      TERM2GENE = wpid2gene,
      TERM2NAME = wpid2name,
      # nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS
    )
  }

  gseaUniprot_wrapper <- function(list, category, species) {
    term2gene <- uniprot %>%
      dplyr::filter(source_name == species[["abbreviation"]]) %>%
      dplyr::select(id, gene)

    term2name <- uniprot %>%
      dplyr::filter(source_name == species[["abbreviation"]]) %>%
      dplyr::select(id, keyword)

    df <- GSEA(
      list,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      # nPerm = gse_nperm,
      pvalueCutoff = gse_pval_cutoff,
      minGSSize = gse_minGS,
      maxGSSize = gse_maxGS
    )
  }

  pathways_function_map <- list(
    "KEGG" = gseaKEGG_wrapper,
    "Reactome" = gseaReactome_wrapper,
    "Wiki" = gseaWiki_wrapper,
    "Uniprot" = gseaUniprot_wrapper
  )
  pathways_function_map[[category]](list, category, species)
}

enrichConsensusdb_wrapper <- function(list, category, species) {
  cpdb_cat <- cpdb %>%
    dplyr::filter(source == category) %>%
    dplyr::mutate(pathway = gsub(" - Homo sapiens (human)", "", pathway)) %>%
    dplyr::mutate(external_id = ifelse(external_id == "None",
      dplyr::row_number(),
      external_id
    )) %>%
    separate_rows(entrez_gene_ids, sep = ",")

  term2gene <- cpdb_cat %>% dplyr::select(external_id, entrez_gene_ids)
  term2name <- cpdb_cat %>% dplyr::select(external_id, pathway)

  df <- enricher(
    names(list),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = enrich_pval_cutoff,
    pAdjustMethod = "BH",
    minGSSize = enrich_minGS,
    maxGSSize = enrich_maxGS
  )
  if (!is.null(df)) {
    df <- setReadable(df, "org.Hs.eg.db", keyType = "ENTREZID")
  }
}

gseaConsensusdb_wrapper <- function(list, category, species) {
  cpdb_cat <- cpdb %>%
    dplyr::filter(source == category) %>%
    dplyr::mutate(pathway = gsub(" - Homo sapiens (human)", "", pathway)) %>%
    dplyr::mutate(external_id = ifelse(external_id == "None",
      dplyr::row_number(),
      external_id
    )) %>%
    separate_rows(entrez_gene_ids, sep = ",")

  term2gene <- cpdb_cat %>% dplyr::select(external_id, entrez_gene_ids)
  term2name <- cpdb_cat %>% dplyr::select(external_id, pathway)

  df <- GSEA(
    list,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    # nPerm = gse_nperm,
    pvalueCutoff = gse_pval_cutoff,
    minGSSize = gse_minGS,
    maxGSSize = gse_maxGS
  )
  if (!is.null(df)) {
    df <- setReadable(df, "org.Hs.eg.db", keyType = "ENTREZID")
  }
}
