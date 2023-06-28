def get_preprocess_data_files(wildcards):
    return expand(RDS_OUTDIR + "input/{contrast}.rds", contrast=pd.unique(Metadata.contrast)) + [RDS_OUTDIR + "input/all_data.rds"]

rule preprocess_data:
    input:
        get_file
    output:
        RDS_OUTDIR + "input/{contrast}.rds"
    params:
        contrast=lambda wc: wc.get("contrast")
    script:
        "../scripts/preprocess_data.R"

rule collate_data:
    # I needed a named input to insert the contrast names as names of the 
    # resulting processed list
    input:
        **(dict(
            zip(
                pd.unique(Metadata.contrast), 
                expand(rules.preprocess_data.output, 
                    contrast=pd.unique(Metadata.contrast))
                )
            )
        )
    output:
        RDS_OUTDIR + "input/all_data.rds"
    script:
        "../scripts/collate_data.R"