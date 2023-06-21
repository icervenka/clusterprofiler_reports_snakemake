def get_cp_output_files(wildcards):
    return expand(RDS_OUTDIR + "{contrast}/{type}_{template}.rds",
        contrast = Metadata.contrast, type = types, template = templates)

rule cp:
    input:
        file=rules.preprocess_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        expand(RDS_OUTDIR + "{{contrast}}/{{type}}_{template}.rds",
            template = templates)
    params:
        contrast=lambda wc: wc.get("contrast"),
        rds_outdir=RDS_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yml"
    script:
        "../scripts/cp.R"
