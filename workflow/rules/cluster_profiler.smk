def get_cp_output_files(wildcards):
    return expand(RDS_OUTDIR + "{contrast}_{template}/{type}.rds",
        contrast = Metadata.contrast, type = types, template = templates["contrast"])

rule cp:
    input:
        file=rules.preprocess_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        expand(RDS_OUTDIR +"{{contrast}}_{template}/{{type}}.rds",
            template = templates["contrast"])
    params:
        contrast=lambda wc: wc.get("contrast"),
        rds_outdir=RDS_OUTDIR,
        template=templates["contrast"]
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/cluster_profiler.R"
