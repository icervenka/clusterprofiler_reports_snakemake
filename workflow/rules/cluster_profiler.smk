def get_outfile_string(outdir, template, suffix):
    return outdir + "{contrast}/{type}" + "_" + template + suffix

# functions for generating lists of output files -------------------------------
def get_cp_output(wildcards):
    return expand(get_outfile_string(RDS_OUTDIR, templates["all"], ".rds"),
        contrast = Metadata.contrast, type = types)

def get_cp_report_output(wildcards):
    return expand(REPORT_OUTDIR + "{contrast}/{type}.html",
        contrast = Metadata.contrast, type = types)

def get_cp_csv_output(wildcards):
    return expand(CSV_OUTDIR + "{contrast}/{type}_ORA.txt",
        contrast = Metadata.contrast, type = types) + \
        expand(CSV_OUTDIR + "{contrast}/{type}_GSEA.txt",
        contrast = Metadata.contrast, type = types)

# standard cluster profiler workflow -------------------------------------------
rule cp:
    input:
        file=rules.preprocess_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        get_outfile_string(RDS_OUTDIR, templates["all"], ".rds")
    params:
        contrast=lambda wc: wc.get("contrast"),
        template=templates["all"],
        rds_outdir=RDS_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/cluster_profiler.R"

rule cp_report:
    input:
        cp_data=rules.cp.output,
        diffexp_data=rules.preprocess_data.output,
        all_data=rules.collate_data.output
    output:
        REPORT_OUTDIR + "{contrast}/{type}.html"
    params:
        contrast=lambda wc: wc.get("contrast"),
        template=templates["all"],
        report_outdir=REPORT_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_reports.R"

rule cp_csv:
    input:
        cp_data=rules.cp.output
    output:
        ora=CSV_OUTDIR + "{contrast}/{type}_ORA.txt",
        gsea=CSV_OUTDIR + "{contrast}/{type}_GSEA.txt"
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_csvs.R"

