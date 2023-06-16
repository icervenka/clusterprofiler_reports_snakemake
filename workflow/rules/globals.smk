# folder organization of the pipelline input and output files
RULES_DIR="workflow/rules/"
CONFIG_DIR="config/"
TYPES_DIR="resources/analysis_types/"
INPUT_DIR="resources/input_files/"
OUTPUT_DIR="results/"
REPORT_OUTDIR=OUTPUT_DIR+"reports/"
CSV_OUTDIR=OUTPUT_DIR+"csvs/"
HEATMAPS_OUTDIR=OUTPUT_DIR+"heatmaps/"
ANALYSIS_PARAMS_OUTDIR=OUTPUT_DIR+"analysis_params/"
LOG_DIR=OUTPUT_DIR+"logs/"
ARCHIVE_OUTDIR=OUTPUT_DIR+"archive/"

RESULT_ARCHIVE_DIRS = [REPORT_OUTDIR, CSV_OUTDIR, HEATMAPS_OUTDIR, LOG_DIR,
ANALYSIS_PARAMS_OUTDIR]