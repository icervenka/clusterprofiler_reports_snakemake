def get_csv_output_files(wildcards):
    return expand(CSV_OUTDIR + "{contrast}/{type}_ORA.txt",
        contrast = Metadata.contrast, type = types) + \
        expand(CSV_OUTDIR + "{contrast}/{type}_GSEA.txt",
        contrast = Metadata.contrast, type = types)

rule create_csvs:
    input:
        cp_data=rules.cp.output
    output:
        ora=CSV_OUTDIR + "{contrast}/{type}_ORA.txt",
        gsea=CSV_OUTDIR + "{contrast}/{type}_GSEA.txt"
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_csvs.R"