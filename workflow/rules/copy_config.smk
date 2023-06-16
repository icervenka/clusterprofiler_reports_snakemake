def get_copy_config_files(wildcards):
    return [ANALYSIS_PARAMS_OUTDIR + "config.yaml", 
            ANALYSIS_PARAMS_OUTDIR + "metadata.tsv"]

rule copy_config:
    input:
        config=CONFIG_DIR+"config.yaml",
        metadata=CONFIG_DIR+"metadata.tsv"
    output:
        config=ANALYSIS_PARAMS_OUTDIR + "config.yaml",
        metadata=ANALYSIS_PARAMS_OUTDIR + "metadata.tsv"
    shell:
        "cp {input.config} {output.config}; cp {input.metadata} {output.metadata}"