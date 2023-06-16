def get_archive_output_files(wildcards):
    return ["pathway_archive.tar.gz"]

rule create_archive:
    input:
        get_cp_output_files,
        get_cp_multi_output_files,
    output:
        "pathway_archive.tar.gz"
    shell:
        "tar cvfz {output} {OUTPUT_DIR}"
