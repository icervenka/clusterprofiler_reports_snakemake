def get_result_archive_output_files(wildcards):
    return ["pathway_archive.tar.gz"]

rule create_archive:
    input:
        directory(OUTPUT_DIR)
    output:
        "pathway_archive.tar.gz"
    shell:
        "tar cvfz {input} {output}"
