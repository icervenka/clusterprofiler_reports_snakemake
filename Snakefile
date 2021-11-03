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
if len(pd.unique(Metadata.contrast)) == 1:
    templates = ['contrast']
else:
    templates = ['contrast', 'unique']
# validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")

include: "snakemake/rules/cluster_profiler.smk"

if len(pd.unique(Metadata.contrast)) > 1:
    include: "snakemake/rules/cluster_profiler_multi.smk"
else:
    include: "snakemake/rules/skip_cluster_profiler_multi.smk"

include: "snakemake/rules/create_archive.smk"

# TODO implement go simplify cutoff
# TODO add ridgplots and upset plot
# TODO clustercompare
# TODO explore new ggplot graphs
# TODO enrichr libraries https://maayanlab.cloud/Enrichr/#libraries
rule all:
    input:
        directory(CSV_OUTDIR),
        get_cp_output_files,
        get_cp_multi_output_files,
        get_archive_output_files,
        OUTPUT_DIR + "analysis_params/config.yaml",
        OUTPUT_DIR + "analysis_params/metadata.tsv"
