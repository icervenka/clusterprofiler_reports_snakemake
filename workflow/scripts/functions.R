# DESeq related functions ------------------------------------------------------
modify_tilda <- function(x, add = T) {
  if (add == T) {
    if (!startsWith(x, "~")) {
      x <- c("~ ", x)
    }
  } else {
    x <- gsub("~", "", x) %>% trimws()
  }
  return(x)
}

name_contrast <- function(c1, c2, sep = "_vs_") {
  return(paste0(c1, sep, c2))
}

is.not.empty.list <- function(l) {
  if (length(l) == 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# deg data manipulation function -----------------------------------------------
load_data <- function(filepath, sp_info, ...) {
  data <- read.delim(filepath, sep = "\t", stringsAsFactors = F) %>%
    drop_na() %>%
    add_idcolnames(sp_info, ...)
  return(data)
}

add_idcolnames <- function(
    d,
    sp_info,
    id_column_name = "ENSEMBL",
    id_type = "ENSEMBL",
    required_colnames = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")) {
  d <- d %>% dplyr::rename(!!id_type := all_of(id_column_name))
  data_colnames <- names(d)

  columns_in <- required_colnames %in% data_colnames

  if (!any(columns_in)) {
    stop(paste0(
      "No ID types present in the file. Has to be one of:\n",
      paste(required_colnames, collapse = " ")
    ))
  }

  if (!all(columns_in)) {
    id_column <- required_colnames[columns_in][1]
    ids <- translate_gene_ids(
      d[, id_column], sp_info,
      id_column,
      to_type = c(required_colnames)
    )

    # TODO remove duplicates
    merge_ids <- merge(ids, d, by = required_colnames[columns_in]) %>%
      dplyr::select(all_of(required_colnames), everything())
    return(merge_ids)
  } else {
    return(d %>% dplyr::select(all_of(required_colnames), everything()))
  }
}

translate_gene_ids <- function(
    gene_ids,
    sp_info,
    from_type,
    to_type = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME"),
    method = "bitr",
    drop = TRUE) {
  if (method == "bitr") {
    ids <- suppressMessages(bitr(
      gene_ids,
      fromType = from_type,
      toType = to_type,
      OrgDb = sp_info["orgdb"],
      drop = drop
    ))
    #   mutate(ENTREZID = as.numeric(ENTREZID)) %>%
    #   group_by(ENTREZID) %>%
    #   arrange(ENTREZID) %>%
    #   slice_head(n = 1) %>%
    #   ungroup

    # TODO needs to be changed, depending on the starting id, produces duplicates
    ids <- ids[!duplicated(ids$ENSEMBL), ]
  } else if (method == "biomart") {
    # TODO default to_type
    mart <- useMart("ensembl", paste0(sp_info[["abbreviation"]], "_gene_ensembl"))
    ids <- getBM(
      filters = "ensembl_gene_id",
      attributes = c(
        "entrezgene_id",
        "ensembl_gene_id",
        paste0(sp_info[["symbol"]], "_symbol")
      ),
      values = gene_ids,
      mart = mart
    ) %>%
      dplyr::rename(c("ENTREZID", "ENSEMBL", "SYMBOL"))

    ids <- ids[!duplicated(ids$ENSEMBL), ]

    if (drop == T) {
      ids <- ids %>% na.omit()
    }
  }
  return(ids)
}

# TODO add rat IDs
add_human_ids <- function(ids, species) {
  if (species["common"] == "mouse") {
    genes <- ids[["ENTREZID"]]
    m2h <- mouse2human(genes, db = homologeneData2) %>%
      dplyr::select(-mouseGene) %>%
      dplyr::rename(
        ENTREZID = mouseID,
        ENTREZID_HUMAN = humanID,
        SYMBOL_HUMAN = humanGene
      )
    id_merge <- merge(ids, m2h, by = c("ENTREZID"))
  } else if (species["common"] == "rat") {
    stop("Rat conversion to human homologues is not yet implemented.")
  } else {
    # id_merge = ids %>%
    #   dplyr::mutate(ENTREZID_HUMAN = "ENTREZID")
  }
  return(id_merge)
}

cp_to_df <- function(cp) {
  if (class(cp) == "gseaResult") {
    df <- data.frame(cp)[1:nrow(data.frame(cp)), ]
  } else {
    df <- fortify(cp, showCategory = nrow(data.frame(cp)))
  }
  return(df)
}

get_fcs <- function(cp, geneList) {
  if (length(cp@gene2Symbol) == 0) {
    stop("Entrez Gene IDs are not translated to readable symbols")
  }
  fcs <- merge(cp@gene2Symbol %>% as_tibble(rownames = "ENTREZID"),
    data.frame(
      ENTREZID = names(geneList),
      log2FoldChange = geneList,
      stringsAsFactors = F
    ),
    by = "ENTREZID"
  )
  return(fcs)
}

get_genes <- function(cp) {
  if (class(cp) == "gseaResult") {
    value <- "core_enrichment"
  } else {
    value <- "geneID"
  }
  genes <- cp_to_df(cp) %>%
    tidyr::separate_rows(!!value, sep = "/") %>%
    pull(!!value)
  return(genes)
}

# analysis options -------------------------------------------------------------
get_template <- function(x) {
  templates <- list(
    "contrast" = list(
      "analyses" = c("ORA", "GSEA"),
      "ORA_gene_sets" = c("all", "up", "down"),
      "GSEA_gene_sets" = c("all")
    ),
    "unique" = list(
      "analyses" = c("ORA"),
      "ORA_gene_sets" = c("all", "up", "down")
    ),
    "multi" = list(
      "analyses" = c("ORA"),
      "ORA_gene_sets" = c("all")
    )
  )
  return(templates[[x]])
}

template_to_df <- function(template) {
  map_dfr(template$analyses, function(x) {
    merge(
      tibble(analysis = x),
      tibble(gene_set = template[[paste0(x, "_gene_sets")]])
    )
  })
}

filter_gene_list <- function(geneList, how) {
  filtered <- list(
    "all" = geneList,
    "up" = geneList[geneList > 0],
    "down" = geneList[geneList < 0]
  )
  return(filtered[[how]])
}

# Species ID accessor ----------------------------------------------------------
get_species_info <- function(species_identifier) {
  species_info <- list(
    mouse_info = c(
      "full_name" = "Mus musculus",
      "abbreviation" = "mmusculus",
      "genus" = "Mus",
      "common" = "mouse",
      "kegg" = "mmu",
      "orgdb" = "org.Mm.eg.db",
      "orgdb_short" = "Mm",
      "meshdb" = "MeSH.Mmu.eg.db",
      "symbol" = "mgi",
      "tax" = 10090
    ),
    human_info = c(
      "full_name" = "Homo sapiens",
      "abbreviation" = "hsapiens",
      "genus" = "Homo",
      "common" = "human",
      "kegg" = "hsa",
      "orgdb" = "org.Hs.eg.db",
      "orgdb_short" = "Hs",
      "meshdb" = "MeSH.Hsa.eg.db",
      "symbol" = "hngc",
      "tax" = 9606
    ),
    rat_info = c(
      "full_name" = "Rattus Norvegicus",
      "abbreviation" = "rnorvegicus",
      "genus" = "Rattus",
      "common" = "rat",
      "kegg" = "rno",
      "orgdb" = "org.Rn.eg.db",
      "orgdb_short" = "Rn",
      "meshdb" = "MeSH.Rno.eg.db",
      "symbol" = "rds",
      "tax" = 10116
    )
  )
  for (x in species_info) {
    if (species_identifier %in% x) {
      sp_arr <- x
      break
    }
  }
  return(sp_arr)
}

get_mesh_dbi <- function(species) {
  ah <- AnnotationHub(localHub = TRUE)
  ah_query <- query(ah, c("MeSHDb", species[["full_name"]]))
  db <- MeSHDbi::MeSHDb(ah_query[[1]])
  return(db)
}

# graph functions  -------------------------------------------------------------
palette_to_hex <- function(palette = "RdBu", no_colors = 100, rev = F) {
  if (typeof(palette) == "closure") {
    pal <- palette(no_colors)
  } else if (typeof(palette) == "character") {
    pal <- colorRampPalette(brewer.pal(n = 9, name = palette))(no_colors)
  } else {
    stop("Unknown type of color palette for graph.")
  }

  if (rev) {
    return(rev(pal))
  } else {
    return(pal)
  }
}

generate_js_color_string <- function(gene_node_types,
                                     add_front = "#333333FF",
                                     ...) {
  pal <- palette_to_hex(...)

  colors <- pal[unique(gene_node_types)]
  if (!is.null(add_front)) {
    colors <- c(add_front, colors)
  }
  colors <- paste0('"', paste(colors, collapse = '", "'), '"')
  colors <- paste0("d3.scaleOrdinal().range([", colors, "]) ;")
  return(colors)
}

create_networkplot_data <- function(cp_df,
                                    geneList,
                                    pathway_categories = NULL,
                                    max_groups = 100) {
  if (is.null(pathway_categories)) {
    pathway_categories <- nrow(data.frame(cp_df))
  } else {
    pathway_categories <- min(c(nrow(data.frame(cp_df)), pathway_categories))
  }

  fortified <- cp_to_df(cp_df)[1:pathway_categories, ]

  if (class(cp_df)[1] == "gseaResult") {
    value <- "NES"
    fortified <- fortified %>%
      dplyr::rename(geneID = "core_enrichment")
  } else {
    value <- "GeneRatio"
  }

  process_nodes <- fortified %>%
    dplyr::select(ID, Description, contains(value)) %>%
    dplyr::rename(
      id = "ID",
      name = "Description",
      value = all_of(value)
    ) %>%
    dplyr::mutate(group = "process") %>%
    as_tibble()

  enriched_genes <- fortified %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    pull(geneID)

  gene_nodes <- get_fcs(cp_df, geneList) %>% setNames(c("ENTREZID", "id", "value"))
  gene_nodes <- gene_nodes %>%
    dplyr::select(-ENTREZID) %>%
    dplyr::mutate(name = id) %>%
    dplyr::select(id, name, value) %>%
    dplyr::filter(id %in% enriched_genes) %>%
    dplyr::arrange(value)

  # this makes the fc values symmetrical around 0, corresponding to the middle of interval
  # small epsilon needs to be added since the cut intervals are closed on one side
  max_val <- max(abs(min(gene_nodes$value)), abs(max(gene_nodes$value))) + 0.00001
  gene_nodes$group <- as.integer(cut(gene_nodes$value,
    breaks = seq(max_val, -max_val, length.out = max_groups)
  ))

  nodes <- rbind.data.frame(process_nodes, gene_nodes) %>%
    as.data.frame() %>%
    dplyr::mutate(id_num = row_number() - 1)

  edges <- fortified %>%
    dplyr::select(ID, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(from = "ID", to = "geneID") %>%
    dplyr::left_join(nodes %>% dplyr::select(id, id_num), by = c("from" = "id")) %>%
    dplyr::rename(from_num = "id_num") %>%
    dplyr::left_join(nodes %>% dplyr::select(id, id_num), by = c("to" = "id")) %>%
    dplyr::rename(to_num = "id_num")

  return(list(
    nodes = as.data.frame(nodes, stringsAsFactors = F),
    edges = as.data.frame(edges, stringsAsFactors = F)
  ))
}

create_emapplot_data <- function(cp_df, pathway_categories = NULL) {
  overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y)) / length(unique(c(x, y)))
  }

  fortified <- cp_to_df(cp_df) %>%
    dplyr::filter(!is.na(Description))

  if (is.null(pathway_categories)) {
    pathway_categories <- nrow(fortified)
  } else {
    pathway_categories <- min(nrow(fortified), pathway_categories)
  }

  geneSets <- clusterProfiler::geneInCategory(cp_df)[fortified$ID[1:pathway_categories]]

  nodes <- sapply(geneSets[names(geneSets)], length) %>%
    as_tibble(rownames = "id") %>%
    dplyr::left_join(fortified %>% dplyr::select(ID, Description), by = c("id" = "ID")) %>%
    dplyr::mutate(value = rescale(value, c(0, 50))) %>%
    dplyr::mutate(id_num = row_number() - 1) %>%
    dplyr::left_join(fortified %>% dplyr::select(ID, p.adjust), by = c("id" = "ID")) %>%
    dplyr::arrange(-p.adjust)
  nodes$group <- as.numeric(cut(nodes$p.adjust, 100))
  nodes <- nodes %>%
    dplyr::arrange(id_num)

  edges <- data.frame()

  if (length(fortified$ID) > 1) {
    geneSets_combinations <- combn(fortified$ID, 2, simplify = T) %>%
      t() %>%
      data.frame(stringsAsFactors = F) %>%
      dplyr::rename(from = "X1", to = "X2")
    geneSets_combinations$overlap <- pmap_dbl(
      geneSets_combinations,
      function(from, to, geneSets) {
        overlap_ratio(geneSets[[from]], geneSets[[to]])
      },
      geneSets = geneSets
    )

    edges <- geneSets_combinations %>%
      dplyr::left_join(nodes %>% dplyr::select(id, id_num), by = c("from" = "id")) %>%
      dplyr::rename(from_num = "id_num") %>%
      dplyr::left_join(nodes %>% dplyr::select(id, id_num), by = c("to" = "id")) %>%
      dplyr::rename(to_num = "id_num") %>%
      dplyr::mutate(width = rescale(overlap, c(0, 75), from = c(0.15, 1))) %>%
      dplyr::filter(overlap >= 0.2)
  }

  if (nrow(edges) == 0) {
    edges <- data.frame(
      from = character(1),
      to = character(1),
      overlap = numeric(1),
      from_num = numeric(1),
      to_num = numeric(1),
      width = numeric(1)
    )
  }

  return(list(
    nodes = as.data.frame(nodes, stringsAsFactors = F),
    edges = as.data.frame(edges, stringsAsFactors = F)
  ))
}

create_legend <- function(vals,
                          palette,
                          no_colors,
                          symmetric = TRUE,
                          rev = FALSE) {
  pal <- palette_to_hex(palette, no_colors, rev)

  if (symmetric) {
    max_val <- max(abs(min(vals)), abs(max(vals)))
    df <- data.frame(x = seq(max_val, -max_val, length.out = no_colors))
  } else {
    df <- data.frame(x = seq(min(vals), max(vals), length.out = no_colors))
  }

  tile_width <- max((max(df$x) - min(df$x)) / (no_colors - 1))

  ggplot(df, aes(fill = x)) +
    geom_tile(aes(
      x = x,
      y = 1,
      width = tile_width,
      height = 1
    )) +
    scale_fill_gradientn(colors = pal) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(min(df$x), max(df$x))) +
    theme_void() +
    theme(
      plot.margin = unit(c(20, 5, 5, 5), "point"),
      axis.text.x = element_text(),
      axis.ticks.length = unit(5, "points"),
      axis.ticks.x = element_line(),
      legend.position = "none"
    )
}

# visualization helpers --------------------------------------------------------
default_dt <- function(x, opts = NULL) {
  def_opts <- list(
    dom = "Blrtip",
    # specify content (search box, etc)
    deferRender = TRUE,
    scrollY = 600,
    scroller = TRUE,
    buttons = list(
      I("colvis"),
      "csv",
      "excel"
    )
  )

  DT::datatable(
    x,
    rownames = T,
    escape = F,
    filter = "top",
    extensions = c("Buttons", "Scroller"),
    style = "bootstrap",
    class = "compact stripe",
    width = "100%",
    options = c(def_opts, opts)
  )
}

# TODO verify if columns are present
format_gene_expr_dt <- function(data, sp_info) {
  data <- data %>%
    dplyr::mutate(
      ENTREZID = paste0(
        '<a href="https://www.ncbi.nlm.nih.gov/gene/',
        ENTREZID,
        '" target="_blank">',
        ENTREZID,
        "</a>"
      ),
      ENSEMBL = paste0(
        '<a href="https://www.ensembl.org/',
        gsub(" ", "_", sp_info["full_name"]),
        "/Gene/Summary?g=",
        ENSEMBL,
        '" target="_blank">',
        ENSEMBL,
        "</a>"
      ),
      SYMBOL = paste0(
        '<span title="',
        GENENAME, '">',
        SYMBOL,
        "</span>"
      )
    ) %>%
    dplyr::select(-GENENAME) %>%
    dplyr::select(
      ENTREZID,
      ENSEMBL,
      SYMBOL,
      !!fc_column_name,
      !!pval_column_name,
      !!padj_column_name,
      any_of(c("SYMBOL_HUMAN", "ENTREZID_HUMAN"))
    )
  return(data)
}

# upset related functions  -----------------------------------------------------
comb_string_to_lgl <- function(s) {
  l <- s %>%
    strsplit(split = "") %>%
    unlist() %>%
    as.numeric() %>%
    as.logical()
  return(l)
}

comb_string_to_degree <- function(s) {
  return(comb_string_to_lgl(s) %>% sum())
}

comb_string_to_set_name <- function(s, m) {
  l <- comb_string_to_lgl(s)
  return(set_name(m)[l])
}

set_contrast_to_comb_string <- function(sc, m, sep = "__int__") {
  tmp <- strsplit(sc, sep) %>% unlist()
  cs <- (set_name(m) %in% tmp) %>%
    as.numeric() %>%
    paste(collapse = "")
  return(cs)
}

filter_data <- function(data,
                        padj_column = "padj",
                        padj_threshold = 0.05,
                        fc_column = "log2FoldChange",
                        fc_threshold = 0.585) {
  data %>%
    dplyr::filter(!!as.symbol(padj_column) < padj_threshold) %>%
    dplyr::filter(abs(!!as.symbol(fc_column)) > fc_threshold)
}

create_set_intersection_name <- function(setid,
                                         m,
                                         remove_prefix = "SET.",
                                         sep = "__int__") {
  n <- paste(gsub(remove_prefix, "", setid) %>% comb_string_to_set_name(m),
    collapse = sep
  )
  return(n)
}

get_upset_list <- function(data_list, id_type = "ENTREZID", ...) {
  map(data_list, function(x) {
    x %>%
      filter_data(...) %>%
      pull(!!as.symbol(id_type))
  })
}

get_upset_sets <- function(data_list,
                           min_set_size,
                           multi = TRUE,
                           id_type = "ENTREZID",
                           set_sep = "__int__") {
  m <- get_upset_list(data_list, id_type = id_type) %>%
    make_comb_mat()

  if (multi == FALSE) {
    usets <- comb_name(m)[comb_degree(m) <= 1 &
      comb_size(m) >= min_set_size]
  } else {
    usets <- comb_name(m)[comb_degree(m) > 1 &
      comb_size(m) >= min_set_size]
  }

  usets_name <- map_chr(usets, ~ create_set_intersection_name(.x, m))

  return(
    data.frame(
      set_id = usets,
      set_name = usets_name,
      comb_size = comb_size(m)[usets] %>% unname(),
      comb_degree = comb_degree(m)[usets] %>% unname(),
      stringsAsFactors = F
    )
  )
}

get_unique_expr_data <- function(data_list,
                                 contrast,
                                 min_set_size = 50,
                                 id_type = "ENTREZID",
                                 set_sep = "__int__") {
  m <- get_upset_list(data_list, id_type = id_type) %>%
    make_comb_mat()

  sn <- set_contrast_to_comb_string(contrast, m, set_sep)
  gene_ids <- extract_comb(m, sn)

  if (length(gene_ids) >= min_set_size) {
    if (comb_string_to_degree(sn) == 1) {
      df <- data_list[[contrast]] %>%
        dplyr::filter(!!as.symbol(id_type) %in% gene_ids)
      return(df)
    } else {
      split_contrast <- strsplit(contrast, set_sep) %>% unlist()
      # TODO make nicer
      df <- dplyr::bind_rows(data_list[split_contrast]) %>%
        dplyr::filter(!!as.symbol(id_type) %in% gene_ids) %>%
        dplyr::group_by(ENTREZID, ENSEMBL, SYMBOL, GENENAME) %>%
        dplyr::summarise(
          across(everything(), ~ mean(.x, na.rm = TRUE)),
          .groups = "drop"
        ) %>%
        as.data.frame(stringsAsFactors = FALSE)
      return(df)
    }
  } else {
    return(data_list[[1]][FALSE, ])
  }
}

# main clusterProfiler functions  ----------------------------------------------
create_subpage_data <- function(expr_data, sp_info, params) {
  if (params$type %in% c("Disease", "ConsensusPathDB") &&
    sp_info["common"] != "human") {
    expr_data <- expr_data %>%
      add_human_ids(sp_info)
    select_entrez_id <- "ENTREZID_HUMAN"
  } else {
    select_entrez_id <- "ENTREZID"
  }

  if (params$analysis %in% c("ORA")) {
    expr_data <- filter_data(
      expr_data,
      padj_column = padj_column_name,
      padj_threshold = gene_padj_threshold,
      fc_column = fc_column_name,
      fc_threshold = gene_fc_threshold
    )
  } else if (params$analysis %in% c("GSEA")) {
    # TODO implement custom ordering
  }

  geneList <- structure(
    expr_data[, fc_column_name],
    names = as.character(expr_data[, select_entrez_id])
  ) %>%
    sort(decreasing = T)
  geneList <- filter_gene_list(geneList, params$gene_set)

  # if (length(filter_gene_list(geneList, params$gene_set)) < 1) {
  #   return(NULL)
  # }

  if (length(geneList) < 1) {
    return(NULL)
  }

  # Generate enrichResult or gseaResult dataframe to pass as data to subpage
  tryCatch(
    {
      df <- get_wrapper(params$type, params$analysis)(
        geneList,
        params$category,
        sp_info)
    },
    error = function(e) {
      df <- NULL
    }
  )

  if (nrow(data.frame(df)) == 0 | is.null(df)) {
    return(NULL)
  }
  if (df@readable == FALSE) {
    df <- setReadable(df, sp_info["orgdb"], keyType = "ENTREZID")
  }

  return(list(
    params = params,
    geneList = geneList,
    df = df,
    expr_data = expr_data
  ))
}

create_subpage_markdown <- function(subpage_data,
                                    subpage_markdown = "scripts/subpage.Rmd") {
  # Create temporary environment which we use for knitting subpages.RMD
  subpage_env <- new.env()

  assign(
    "subpage_item",
    paste0(
      subpage_data$params$analysis,
      " - ",
      subpage_data$params$gene_set
    ),
    subpage_env
  )
  assign(
    "subpage_category",
    subpage_data$params$heading,
    subpage_env
  )
  assign(
    "subpage_heading",
    paste0(
      subpage_data$params$heading,
      " - ",
      subpage_data$params$analysis,
      " - ",
      subpage_data$params$gene_set
    ),
    subpage_env
  )
  assign("params", subpage_data$params, subpage_env)
  assign("geneList", subpage_data$geneList, subpage_env)
  assign("subpage_df", subpage_data$df, subpage_env)

  return(knitr::knit_child(subpage_markdown, envir = subpage_env))
}

create_pathway_csv <- function(subpage_data, template, path) {
  df <- cp_to_df(subpage_data$df)
  params <- subpage_data$params

  if (class(subpage_data$df) == "gseaResult") {
    value <- "NES"
    df <- df %>% dplyr::rename(geneID = "core_enrichment")
  } else {
    value <- "GeneRatio"
  }

  fcs <- get_fcs(subpage_data$df, subpage_data$geneList) %>%
    setNames(c("ENTREZID", "geneID", "log2FoldChange"))

  df <- df %>%
    dplyr::select(ID, Description, contains(value), pvalue, p.adjust, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    as_tibble() %>%
    left_join(fcs, by = "geneID") %>%
    dplyr::rename(SYMBOL = geneID)

  filename <- paste0(
    gsub(
      " ",
      "_",
      paste(
        params$type,
        params$category,
        params$analysis,
        params$gene_set
      )
    ),
    "_", template, ".txt"
  )

  write.table(
    df,
    paste0(path, filename),
    quote = F,
    sep = "\t",
    row.names = F
  )
}

# TODO uses fc to order, might be useful to use custom ordering
# TODO add knitr options

run_cp <- function(data, sp_info, payloads, template, outdir, output_opts = list()) {
  map(payloads$type %>% unique(), function(x, data, sp_info) {
    payloads_sub <- payloads %>% dplyr::filter(type == x)
    render(cp_script_path,
      output_file = paste0(x, "_", template, ".html"),
      output_format = "all",
      output_dir = paste0(report_outdir, outdir),
      output_options = output_opts
    )
  }, data = data, sp_info = sp_info)
}
