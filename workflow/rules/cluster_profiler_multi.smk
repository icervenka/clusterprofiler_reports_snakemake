def get_cp_multi_output_files(wildcards):
    return dynamic(expand(RDS_OUTDIR + "{{contrast}}/{type}_multi.rds", type = types))

rule cp_multi:
    input:
        type=TYPES_DIR + "{type}.csv",
        all_data=rules.collate_data.output,
    output:
        dynamic(RDS_OUTDIR + "{contrast}/{type}" + ".rds"),
    params:
        contrast=lambda wc: wc.get("contrast"),
        rds_outdir=RDS_OUTDIR,
        template=multi_template
        # input_dir=INPUT_DIR,
        # report_outdir=REPORT_OUTDIR,
        # csv_outdir=CSV_OUTDIR,
    conda:
        "../envs/clusterprofiler_reports.yml"
    script:
        "../scripts/cp_multi.R"
