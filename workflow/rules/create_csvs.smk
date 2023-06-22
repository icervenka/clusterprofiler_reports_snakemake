def get_csv_output_files(wildcards):
    return expand(CSV_OUTDIR + "{contrast}/{type}_{template}_ORA.txt",
        contrast = Metadata.contrast, type = types, template = templates) + \
        expand(CSV_OUTDIR + "{contrast}/{type}_{template}_GSEA.txt",
        contrast = Metadata.contrast, type = types, template = templates)

rule create_csvs:
    input:
        cp_data=rules.cp.output
    output:
        ora=expand(CSV_OUTDIR + "{{contrast}}/{{type}}_{template}_ORA.txt",
            template = templates),
        gsea=expand(CSV_OUTDIR + "{{contrast}}/{{type}}_{template}_GSEA.txt",
            template = templates)
    conda:
        "../envs/clusterprofiler_reports.yaml"
    script:
        "../scripts/create_csvs.R"