

def get_report_output_files(wildcards):
    ANALYSES, = glob_wildcards(RDS_OUTDIR + "{analysis}/{type}.rds")
    return expand(REPORT_OUTDIR + "{analysis}/{type}.html",
        analysis =ANALYSES, type = types)

rule create_reports:
    input:
        cp_data=RDS_OUTDIR + "{analysis}/{type}.rds",
        diffexp_data=rules.preprocess_data.output,
        all_data=rules.collate_data.output,
    output:
        REPORT_OUTDIR + "{analysis}/{type}.html"
    params:
        contrast=lambda wc: wc.get("contrast"),
        report_outdir=REPORT_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_reports.R"