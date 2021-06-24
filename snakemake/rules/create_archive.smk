rule create_archive:
    input:
        directory(REPORT_OUTDIR)
        directory(CSV_OUTDIR)
    output:
        "pathway_archive.tar.gz"
    shell:
        "tar cvf {input} {output}"
