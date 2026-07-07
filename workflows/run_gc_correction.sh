#!/usr/bin/env bash
#
# run_gc_correction.sh
# --------------------
# Runs Griffin's GC-correction part (GC_counts -> GC_bias) for one or more
# samples straight from the terminal -- no Snakemake, no SLURM. Each sample
# uses CPU cores internally via the --CPU flag of griffin_GC_counts.py.
#
# This script is the GC-CORRECTION STAGE ONLY. It is deliberately self-contained
# so it can be run on its own (someone may only need GC correction). The
# nucleosome-profiling stage reads the samples.GC.yaml written at the end.
#
# Usage:
#   1. edit the CONFIG block below
#   2. create a tab-separated samples file:   <sample_name>\t<path/to/sample.bam>
#   3. ./run_gc_correction.sh samples.tsv
#
set -euo pipefail

# ============================ CONFIG (edit me) ============================
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "$0")/../gc_correction" 2>/dev/null && pwd)}"   # where the griffin_*.py live
OUT_DIR="${OUT_DIR:-./griffin_gc_output}"                     # results go here

REF_SEQ="${REF_SEQ:-/path/to/hg38.fa}"                        # reference fasta (indexed: .fai)
CHROM_SIZES="${CHROM_SIZES:-/path/to/hg38.standard.chrom.sizes}"
MAPPABLE_REGIONS="${MAPPABLE_REGIONS:-/path/to/k100_minus_exclusion_lists.mappable_regions.hg38.bed}"
GENOME_GC_FREQUENCY="${GENOME_GC_FREQUENCY:-/path/to/genome_GC_frequency}"  # folder from genome_GC_frequency step

MAP_Q="${MAP_Q:-20}"
SIZE_RANGE="${SIZE_RANGE:-15 500}"   # two ints: min max  (e.g. "15 500")
CPU="${CPU:-8}"
NO_PLOTS="${NO_PLOTS:-0}"            # set to 1 to skip the GC_bias QC PDFs (faster; .txt output unchanged)

PYTHON="${PYTHON:-python3}"
# =========================================================================

SAMPLES_FILE="${1:-samples.tsv}"

# mappable_name = basename with the final extension (.bed) removed,
# matching the original snakefile: rsplit('/',1)[1].rsplit('.',1)[0]
MAPPABLE_BASENAME="$(basename "$MAPPABLE_REGIONS")"
MAPPABLE_NAME="${MAPPABLE_BASENAME%.*}"

GC_COUNTS_SCRIPT="$SCRIPTS_DIR/griffin_GC_counts.py"
GC_BIAS_SCRIPT="$SCRIPTS_DIR/griffin_GC_bias.py"

echo "=== Griffin GC correction (terminal driver) ==="
echo "scripts_dir         = $SCRIPTS_DIR"
echo "out_dir             = $OUT_DIR"
echo "ref_seq             = $REF_SEQ"
echo "mappable_regions    = $MAPPABLE_REGIONS  (mappable_name=$MAPPABLE_NAME)"
echo "genome_GC_frequency = $GENOME_GC_FREQUENCY"
echo "map_q=$MAP_Q  size_range=$SIZE_RANGE  CPU=$CPU"
echo "samples_file        = $SAMPLES_FILE"
echo

# --- sanity checks ---
[[ -f "$SAMPLES_FILE" ]]      || { echo "ERROR: samples file not found: $SAMPLES_FILE" >&2; exit 1; }
[[ -f "$GC_COUNTS_SCRIPT" ]]  || { echo "ERROR: missing $GC_COUNTS_SCRIPT" >&2; exit 1; }
[[ -f "$GC_BIAS_SCRIPT" ]]    || { echo "ERROR: missing $GC_BIAS_SCRIPT" >&2; exit 1; }
[[ -f "$REF_SEQ" ]]          || { echo "ERROR: ref_seq not found: $REF_SEQ" >&2; exit 1; }

mkdir -p "$OUT_DIR"/{GC_counts,GC_bias,GC_plots}

# Collect (name, bam) pairs for the yaml written at the end.
declare -a NAMES=()
declare -a BAMS=()

# --- per-sample loop (sequential; each sample is multicore via --CPU) ---
while IFS=$'\t' read -r SAMPLE_NAME BAM_PATH || [[ -n "$SAMPLE_NAME" ]]; do
    # skip blank lines and comments
    [[ -z "${SAMPLE_NAME// }" ]] && continue
    [[ "$SAMPLE_NAME" == \#* ]] && continue

    echo ">>> [$SAMPLE_NAME] $BAM_PATH"
    [[ -f "$BAM_PATH" ]] || { echo "ERROR: bam not found: $BAM_PATH" >&2; exit 1; }

    echo "    step 1/2: GC_counts"
    "$PYTHON" "$GC_COUNTS_SCRIPT" \
        --bam_file "$BAM_PATH" \
        --bam_file_name "$SAMPLE_NAME" \
        --mappable_regions_path "$MAPPABLE_REGIONS" \
        --ref_seq "$REF_SEQ" \
        --chrom_sizes "$CHROM_SIZES" \
        --out_dir "$OUT_DIR" \
        --map_q "$MAP_Q" \
        --size_range $SIZE_RANGE \
        --CPU "$CPU"

    echo "    step 2/2: GC_bias"
    NO_PLOTS_FLAG=""
    [[ "$NO_PLOTS" == "1" ]] && NO_PLOTS_FLAG="--no_plots"
    "$PYTHON" "$GC_BIAS_SCRIPT" \
        --bam_file_name "$SAMPLE_NAME" \
        --mappable_name "$MAPPABLE_NAME" \
        --genome_GC_frequency "$GENOME_GC_FREQUENCY" \
        --out_dir "$OUT_DIR" \
        --size_range $SIZE_RANGE \
        $NO_PLOTS_FLAG

    NAMES+=("$SAMPLE_NAME")
    BAMS+=("$BAM_PATH")
    echo "    done [$SAMPLE_NAME]"
    echo
done < "$SAMPLES_FILE"

# --- write samples.GC.yaml (same format the snakefile produced) ---
YAML="$OUT_DIR/samples.GC.yaml"
{
    echo "samples:"
    for i in "${!NAMES[@]}"; do
        name="${NAMES[$i]}"
        bam="${BAMS[$i]}"
        gc_bias_path="$(cd "$OUT_DIR" && pwd)/GC_bias/${name}.GC_bias.txt"
        echo "  ${name}:"
        echo "    bam: ${bam}"
        echo "    GC_bias: ${gc_bias_path}"
    done
} > "$YAML"

echo "=== GC correction complete. Wrote $YAML ==="
