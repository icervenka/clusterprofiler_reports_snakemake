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

group_order <- data.frame(
    group = str_split(snakemake@params[["group_order"]], ",")[[1]],
    stringsAsFactors = FALSE
)

sample_sheet <- read_delim(
    snakemake@input[["sample_sheet"]],
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
) %>%
    select(sample, group) %>%
    unique()

# ordering the metadata dataframe based on the grouporder
# specified on commandline
sample_sheet <- grouporder %>%
    left_join(sample_sheet, by = "group")

sample_expression <- read_delim(
    snakemake@input[["sample_expression"]],
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
) %>%
    mutate(SYMBOL = toupper(SYMBOL))

all_pathways <- map(snakemake@input[["cp_data"]], function(x) {
    read_delim(
        x,
        delim = "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        show_col_types = FALSE
    ) %>%
        # cast all the IDs to character, since some of the IDs are numeric
        dplyr::mutate(ID = as.character(ID))
})

filter_pathways <- function(pathways, heatmap_config) {
    # filter pathway ids based on specified files if there are any
    if (length(heatmap_config$pathway_id_files) > 0) {
        pathway_ids <- map_dfr(heatmap_config$pathway_id_files, function(x) {
            read_delim(
                x,
                delim = "\t",
                col_names = c("id"),
                escape_double = FALSE,
                trim_ws = TRUE,
                show_col_types = FALSE
            )
        }) %>%
            unique()
        pathways %>%
            filter(ID %in% pathway_ids$id)
    # if no files are specified filter according to config criteria
    } else {
        pathways %>%
            filter(type %in% heatmap_config$type$include) %>%
            filter(!(type %in% heatmap_config$type$exclude)) %>%
            filter(category %in% heatmap_config$category$include) %>%
            filter(!(category %in% heatmap_config$category$exclude)) %>%
            filter(analysis %in% heatmap_config$analysis)
    }
}

create_heatmap <- function(data, sample_sheet) {
    # generate heatmap, no column clustering
    p <- pheatmap::pheatmap(
        data,
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(
            n = 7, name = "RdBu"
        )))(100),
        cluster_cols = FALSE,
        scale = "row",
        cellwidth = 10,
        cellheight = 10,
        border_color = "white",
        treeheight_row = 10,
        annotation_col = sample_sheet %>%
            column_to_rownames("sample") %>%
            select(matches("group")),
        annotation_legend = FALSE,
        show_colnames = FALSE
    )
    return(p)
}

save_heatmap <- function(plot, outdir, filename, width, height) {
    # save as image, paste filename from pathway name
    # and original analysis name
    # dimensions are empirically determined
    ggplot2::ggsave(
        paste0(
            outdir,
            filename,
            ".png"
        ), # \\/:"*?<>|
        heatmap,
        units = "mm",
        width = width,
        height = height
    )
}

# walk through all the pathway IDs
cat("Exporting heatmap for:\n")
walk(all_pathways, function(x, samples, expression, heatmap_config) {
    y <- x %>%
        filter_pathways(heatmap_config)

    walk(y$ID %>% unique(), function(a) {
        z <- y %>%
            filter(ID == a)

        # code for merging with sample expression
        # TODO
        
        num_genes <- nrow(z)
        num_samples <- length(samples$sample)

        # code for creating file name
        # TODO
        filename <- "heatmap"

        p <- create_heatmap(z, samples)
        save_heatmap(
            plot = p,
            outdir = HEATMAP_OUTDIR,
            filename = filename,
            width = num_samples * 3.5 + 40,
            height = num_genes * 3.5 + 20
        )
    })
},
heatmap_config = heatmaps,
samples = sample_sheet,
expression = sample_expression
)
