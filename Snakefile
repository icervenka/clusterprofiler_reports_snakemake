import pandas as pd
import glob
import os
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.7.0")

##### load config and sample sheets #####
include: "snakemake/rules/common.smk"
include: "snakemake/rules/functions.smk"
configfile: "config.yaml"
validate(config, schema="snakemake/schema/config.schema.yaml")


try:
    os.mkdir(OUTPUT_DIR)
except OSError as e:
    print(e)

try:
    os.mkdir(CSV_OUTDIR)
except OSError as e:
    print(e)

Metadata = pd.read_table(config["metadata"])
validate(Metadata, schema="snakemake/schema/metadata.schema.yaml")

types = config["types"]
if len(pd.unique(Metadata.contrast)) == 1:
    templates = ['contrast']
else:
    templates = ['contrast', 'unique']


include: "snakemake/rules/cluster_profiler.smk"

if len(pd.unique(Metadata.contrast)) > 1:
    include: "snakemake/rules/cluster_profiler_multi.smk"
else:
    include: "snakemake/rules/skip_cluster_profiler_multi.smk"

include: "snakemake/rules/create_archive.smk"

# TODO add ridgeplots and upset plot
# TODO fix the overwriting color warning
# TODO fix the ggplotly warning
# TODO clustercompare
# TODO explore new ggplot graphs
# TODO enrichr libraries https://maayanlab.cloud/Enrichr/#libraries
rule all:
    input:
        get_cp_output_files,
        get_cp_multi_output_files,
        get_archive_output_files,
        OUTPUT_DIR + "analysis_params/config.yaml",
        OUTPUT_DIR + "analysis_params/metadata.tsv"
