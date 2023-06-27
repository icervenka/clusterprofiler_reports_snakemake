def get_cp_unique_output_files(wildcards):
    return expand(RDS_OUTDIR + "{contrast}_{template}/{type}.rds",
        contrast = Metadata.contrast, type = types, template = templates["unique"])

rule cp:
    input:
        file=rules.preprocess_data.output,
        all_data=rules.collate_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        RDS_OUTDIR + "{contrast}_{template}/{type}.rds"
    params:
        contrast=lambda wc: wc.get("contrast"),
        rds_outdir=RDS_OUTDIR,
        template=uniqe_template
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/cluster_profiler_unique.R"
