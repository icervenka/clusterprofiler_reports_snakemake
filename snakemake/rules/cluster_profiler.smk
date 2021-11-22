def get_cp_output_files(wildcards):
    return expand(REPORT_OUTDIR + "{contrast}/{type}_{template}.html",
        contrast = Metadata.contrast, type = types, template = templates)

rule cp:
    input:
        file=get_file,
        type=TYPES_DIR + "{type}.csv"
    output:
        expand(REPORT_OUTDIR + "{{contrast}}/{{type}}_{template}.html",
            template = templates)
    params:
        contrast=get_contrast,
        input_dir=INPUT_DIR,
        report_outdir=REPORT_OUTDIR,
        csv_outdir=CSV_OUTDIR,
        metadata=config['metadata'],
        species=config["species"],
        self_contained=config["self_contained"],
        id_column_name=config["id_column_name"],
        id_type=config["id_type"],
        fc_column_name=config["fc_column_name"],
        padj_column_name=config["padj_column_name"],
        pval_column_name=config["pval_column_name"],
        gene_fc_threshold=config["gene_fc_threshold"],
        gene_padj_threshold=config["gene_padj_threshold"],
        enrich_minGS=config["enrich_minGS"],
        enrich_maxGS=config["enrich_maxGS"],
        enrich_pval_cutoff=config["enrich_pval_cutoff"],
        min_set_size=config["min_set_size"],
        #gse_nperm=config["gse_nperm"],
        gse_minGS=config["gse_minGS"],
        gse_maxGS=config["gse_maxGS"],
        gse_pval_cutoff=config["gse_pval_cutoff"],
        dotplot_categories=config["dotplot"]["categories"],
        networkplot_categories=config["networkplot"]["categories"],
        heatplot_categories=config["heatplot"]["categories"],
        emapplot_categories=config["enrichmap"]["categories"],
        hide_columns=config["table"]["hide_columns"]
    script:
        "../scripts/cp.R"

rule copy_config:
    input:
        config="config.yaml",
        metadata="metadata.tsv"
    output:
        config=OUTPUT_DIR + "analysis_params/config.yaml",
        metadata=OUTPUT_DIR + "analysis_params/metadata.tsv"
    shell:
        "cp {input.config} {output.config}; cp {input.metadata} {output.metadata}"
