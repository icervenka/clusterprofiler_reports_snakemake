rule create_analysis_types:
    input:
        "config/config.yaml"
    output:
        temp(expand(TYPES_DIR + "{type}.csv", type = config["types"]))
    params:
        outdir=TYPES_DIR
    shell:
        "workflow/scripts/export_analysis_types.R {params.outdir}"