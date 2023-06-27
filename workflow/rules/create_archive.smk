def get_archive_output_files(wildcards):
    return [ARCHIVE_OUTDIR + config["experiment_name"] + "_pathway_archive.tar.gz"]

rule create_archive:
    input:
        get_cp_output_files,
        get_cp_multi_output_files,
        get_report_output_files,
        get_csv_output_files,
        # get_heatmap_files,
        get_copy_config_files
    output:
        ARCHIVE_OUTDIR + config["experiment_name"] + "_pathway_archive.tar.gz"
    run:
        include_dirs = []
        for item in RESULT_ARCHIVE_DIRS:
            if(os.path.exists(item)):
                include_dirs = include_dirs + [item]
        shell(
            "tar "
            "-czf "    
            "{output} "
            "{include_dirs} "
        )
