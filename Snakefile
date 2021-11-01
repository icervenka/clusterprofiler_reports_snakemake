import pandas as pd
import glob
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.7.0")

##### load config and sample sheets #####
include: "snakemake/rules/common.smk"
include: "snakemake/rules/functions.smk"
configfile: "config.yaml"
# validate(config, schema="snakemake/schema/config.schema.yaml")

Metadata = pd.read_table(config["metadata"])
types = config["clusterProfiler"]["types"]
# validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")

rule all:
    input:
        expand(REPORT_OUTDIR +  "{contrast}/{type}_{template}.html",
            contrast = Metadata.contrast, type = types, template = ['contrast', 'unique']),
        # dynamic(expand(REPORT_OUTDIR + "{{contrast}}/{type}_multi.html", type = types)),
        REPORT_OUTDIR + "analysis_params/config.yaml",
        REPORT_OUTDIR + "analysis_params/metadata.tsv"

include: "snakemake/rules/cluster_profiler.smk"

if len(Metadata.contrast) > 1:
    include: "snakemake/rules/cluster_profiler_multi.smk"
