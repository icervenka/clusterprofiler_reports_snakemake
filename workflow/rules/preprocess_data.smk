def get_preprocess_data_files(wildcards):
    return [expand(RDS_OUTDIR + "{contrast}.rds", contrast=pd.unique(Metadata.contrast)) +
    RDS_OUTDIR + "all_data.rds"]

rule preprocess_data:
    input:
        get_file
    output:
        RDS_OUTDIR + "{contrast}.rds"
    params:
        input_contrast=lambda wc: wc.get("contrast")
    script:
        "../scripts/preprocess_data.R"

rule collate_data:
    input:
        files=expand(rules.preprocess_data.output, contrast=pd.unique(Metadata.contrast))
    output:
        RDS_OUTDIR + "all_data.rds"
    script:
        "../scripts/collate_data.R"