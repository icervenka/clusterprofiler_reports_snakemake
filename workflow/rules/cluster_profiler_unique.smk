def get_cp_unique_output(wildcards):
    return expand(RDS_OUTDIR + "{contrast}/{type}" + "_" + templates["unique"] + ".rds",
        contrast = Metadata.contrast, type = types)

rule cp_unique:
    input:
        file=rules.preprocess_data.output,
        all_data=rules.collate_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        RDS_OUTDIR + "{contrast}/{type}" + "_" + templates["unique"] + ".rds"
    params:
        contrast=lambda wc: wc.get("contrast"),
        template=templates["unique"],
        rds_outdir=RDS_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/cluster_profiler_unique.R"
