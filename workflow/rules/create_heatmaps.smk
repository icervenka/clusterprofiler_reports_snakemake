def get_heatmap_output_files(wildcards):
    return expand(LOG_DIR + "{contrast}/{type}_{template}.done",
            contrast = Metadata.contrast, type = types, template = templates)

rule create_heatmaps:
    input:
        cp_data=rules.create_csvs.output.ora + rules.create_csvs.output.gsea,
        sample_sheet=get_sample_sheet,
        sample_expression=get_sample_expression
    output:
        expand(LOG_DIR + "{{contrast}}/{{type}}_{template}.done",
            template = templates)
    params:
        output_dir=HEATMAP_OUTDIR,
        group_order=get_group_order
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/pathways2heatmaps.R"