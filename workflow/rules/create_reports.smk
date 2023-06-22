def get_reports_output_files(wildcards):
    return expand(REPORT_OUTDIR + "{contrast}/{type}_{template}.html",
        contrast = Metadata.contrast, type = types, template = templates)

rule create_reports:
    input:
        cp_data=rules.cp.output,
        diffexp_data=rules.preprocess_data.output,
        all_data=rules.collate_data.output,
    output:
        expand(REPORT_OUTDIR + "{{contrast}}/{{type}}_{template}.html",
            template = templates)
    params:
        contrast=lambda wc: wc.get("contrast"),
        report_outdir=REPORT_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_reports.R"