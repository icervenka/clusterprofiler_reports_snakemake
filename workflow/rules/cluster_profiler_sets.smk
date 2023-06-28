# general use functions --------------------------------------------------------
def get_outfile_string(outdir, template, suffix):
    return outdir + "{contrast}/{type}" + "_" + template + suffix

# output functions -------------------------------------------------------------
def get_cp_unique_output(wildcards):
    return expand(get_outfile_string(RDS_OUTDIR, templates["unique"], ".rds"),
        contrast = Metadata.contrast, type = types)

def get_cp_unique_report_output(wildcards):
    return expand(get_outfile_string(REPORT_OUTDIR, templates["unique"], ".html"),
        contrast = Metadata.contrast, type = types)

def get_cp_unique_csv_output(wildcards):
    return expand(get_outfile_string(CSV_OUTDIR, templates["unique"], "_ORA.txt"),
        contrast = Metadata.contrast, type = types) + \
        expand(get_outfile_string(CSV_OUTDIR, templates["unique"], "_GSEA.txt"),
        contrast = Metadata.contrast, type = types)

def get_cp_set_output(wildcards):
    return expand(LOG_DIR + "{type}.done", type = types)

# rules for unique genes in contrast -------------------------------------------

rule cp_unique:
    input:
        file=rules.preprocess_data.output,
        all_data=rules.collate_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        get_outfile_string(RDS_OUTDIR, templates["unique"], ".rds")
    params:
        contrast=lambda wc: wc.get("contrast"),
        template=templates["unique"],
        rds_outdir=RDS_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/cluster_profiler.R"

rule cp_unique_report:
    input:
        cp_data=rules.cp_unique.output,
        diffexp_data=rules.preprocess_data.output,
        all_data=rules.collate_data.output
    output:
        get_outfile_string(REPORT_OUTDIR, templates["unique"], ".html")
    params:
        contrast=lambda wc: wc.get("contrast"),
        template=templates["unique"],
        report_outdir=REPORT_OUTDIR
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_reports.R"

rule cp_unique_csv:
    input:
        cp_data=rules.cp_unique.output,
    output:
        ora=get_outfile_string(CSV_OUTDIR, templates["unique"], "_ORA.txt"),
        gsea=get_outfile_string(CSV_OUTDIR, templates["unique"], "_GSEA.txt"),
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_csvs.R"

# rules for set combinations ---------------------------------------------------
rule cp_set:
    input:
        all_data=rules.collate_data.output,
        type=TYPES_DIR + "{type}.csv"
    output:
        LOG_DIR + "{type}.done"
    params:
        type=lambda wc: wc.get("type"),
        template=templates["set"],
        rds_outdir=RDS_OUTDIR,
        log_dir=LOG_DIR
    script:
        "../scripts/cluster_profiler_sets.R"

# rule: cp_set_report:
#     input:
#         aggregate_cp_set_input
#     output:
#          REPORT_OUTDIR + "{set_contrast}/{type}.html"
#     params:
#         contrast=lambda wc: wc.get("set_contrast"),
#         template=templates["all"],
#         report_outdir=REPORT_OUTDIR
#     conda:
#         "../envs/clusterprofiler_reports.yaml"
#     script:
#         "../scripts/create_reports.R"