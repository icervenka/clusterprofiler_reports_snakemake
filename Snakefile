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
if len(pd.unique(Metadata.contrast))) == 1:
    templates = ['contrast']
else:
    templates = ['contrast', 'unique']
# validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")

rule all:
    input:
        get_cp_output_files,
        get_cp_multi_output_files,
        get_result_archive_output_files,
        OUTPUT_DIR + "analysis_params/config.yaml",
        OUTPUT_DIR + "analysis_params/metadata.tsv"

include: "snakemake/rules/cluster_profiler.smk"

if len(pd.unique(Metadata.contrast))) > 1:
    include: "snakemake/rules/cluster_profiler_multi.smk"
else:
    include: "snakemake/rules/skip_cluster_profiler_multi.smk"

include: "snakemake/rules/result_archive.smk"
