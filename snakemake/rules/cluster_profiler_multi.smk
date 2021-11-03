def get_cp_multi_output_files(wildcards):
    return dynamic(expand(REPORT_OUTDIR + "{{contrast}}/{type}_multi.html", type = types))

rule cp_multi:
    input:
        type=TYPES_DIR + "{type}.csv"
    output:
        dynamic(REPORT_OUTDIR + "{contrast}/{type}" + "_multi.html"),
    params:
        input_dir=INPUT_DIR,
        report_outdir=REPORT_OUTDIR,
        csv_outdir=CSV_OUTDIR,
        metadata=config['metadata'],
        species=config["species"],
        self_contained=config["clusterProfiler"]["self_contained"],
        id_column_name=config["clusterProfiler"]["id_column_name"],
        id_type=config["clusterProfiler"]["id_type"],
        fc_column_name=config["clusterProfiler"]["fc_column_name"],
        padj_column_name=config["clusterProfiler"]["padj_column_name"],
        pval_column_name=config["clusterProfiler"]["pval_column_name"],
        gene_fc_threshold=config["clusterProfiler"]["gene_fc_threshold"],
        gene_padj_threshold=config["clusterProfiler"]["gene_padj_threshold"],
        enrich_minGS=config["clusterProfiler"]["enrich_minGS"],
        enrich_maxGS=config["clusterProfiler"]["enrich_maxGS"],
        enrich_pval_cutoff=config["clusterProfiler"]["enrich_pval_cutoff"],
        min_set_size=config["clusterProfiler"]["min_set_size"],
        gse_nperm=config["clusterProfiler"]["gse_nperm"],
        gse_minGS=config["clusterProfiler"]["gse_minGS"],
        gse_maxGS=config["clusterProfiler"]["gse_maxGS"],
        gse_pval_cutoff=config["clusterProfiler"]["gse_pval_cutoff"],
        dotplot_categories=config["clusterProfiler"]["dotplot"]["categories"],
        networkplot_categories=config["clusterProfiler"]["networkplot"]["categories"],
        heatplot_categories=config["clusterProfiler"]["heatplot"]["categories"],
        emapplot_categories=config["clusterProfiler"]["enrichmap"]["categories"],
        hide_columns=config["clusterProfiler"]["table"]["hide_columns"]
    script:
        "../scripts/cp_multi.R"
