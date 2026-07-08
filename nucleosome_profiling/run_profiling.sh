#!/bin/bash
# =============================================================================
# run_profiling.sh  —  end-to-end nucleosome profiling for ONE sample.
#
#   02 (FCC count, if needed)  ->  per-chromosome profiling  ->  merge + features.
#   No Slurm, no Snakemake. Runs from a terminal or inside the container.
#
# USAGE:
#   bash run_profiling.sh --sample <SAMPLE> [--config config.yaml] [--tfs SEL]
#                         [--force-recount] [--threads N] [--results-root DIR]
#
#   --sample         (required) sample_id; must be a unique row in the sample sheet.
#   --config         path to config.yaml (default: ./config.yaml).
#   --tfs            ALL | FOXA1 | FOXA1,GRHL2,SPI1 | /path/to/tf_list.txt  (default ALL)
#   --force-recount  recount FCC even if it already exists.
#   --threads        override config threads for this run.
#
# All paths and parameters come from config.yaml + the sample sheet, read ONLY
# by config_helper.py (Python). This script parses no YAML itself.
# =============================================================================

set -euo pipefail

# ------------------------------- LOCATE SELF ---------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER="$SCRIPT_DIR/config_helper.py"
PROF_PY="$SCRIPT_DIR/04_profile_tf_sites.py"
MERGE_PY="$SCRIPT_DIR/04_merge_and_extract.py"
COUNT_SH="$SCRIPT_DIR/02_fcc_count.sh"
PY="${PY:-python}"

# ------------------------------- ARGS ----------------------------------------
CONFIG="./config.yaml"
SAMPLE=""
TF_SELECTION="ALL"
FORCE_RECOUNT=false
THREADS_OVERRIDE=""
RESULTS_ROOT_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)        CONFIG="$2";           shift 2 ;;
        --sample)        SAMPLE="$2";           shift 2 ;;
        --tfs)           TF_SELECTION="$2";     shift 2 ;;
        --force-recount) FORCE_RECOUNT=true;    shift 1 ;;
        --threads)       THREADS_OVERRIDE="$2"; shift 2 ;;
        --results-root)  RESULTS_ROOT_OVERRIDE="$2"; shift 2 ;;
        -h|--help)       sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; exit 1 ;;
    esac
done
[ -n "$SAMPLE" ] || { echo "ERROR: --sample is required" >&2; exit 1; }

# ------------------------------- PRE-FLIGHT ----------------------------------
# Fail fast on missing tools/files before doing any real work.
for tool in "$PY" gawk samtools parallel; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH" >&2; exit 1; }
done
for f in "$HELPER" "$PROF_PY" "$MERGE_PY" "$COUNT_SH" "$CONFIG"; do
    [ -e "$f" ] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

# ------------------------------- READ CONFIG ---------------------------------
# Driver scalars (deployment) + the method flags 02 needs. config_helper emits
# shell-safe UPPERCASE assignments; we eval them.
eval "$("$PY" "$HELPER" config "$CONFIG" \
        sample_sheet reference_fa mappable_bed chrom_sizes tf_dir results_root tmpdir threads \
        frag_len_min frag_len_max mapq gc_bias_min autosomes)"

# Resolve the sample row (errors nonzero if missing or duplicated).
eval "$("$PY" "$HELPER" sample "$SAMPLE_SHEET" "$SAMPLE")"
# -> BAM, GC_BIAS, FCC_DIR

# Optional threads override for this run.
[ -n "$THREADS_OVERRIDE" ] && THREADS="$THREADS_OVERRIDE"

# Optional results_root override (CLI beats config); redirects all Stage-2
# outputs (profiling/, logs/, features) under the given root. FCC location is
# unaffected here — it comes from the sample sheet's fcc_dir column.
[ -n "$RESULTS_ROOT_OVERRIDE" ] && RESULTS_ROOT="${RESULTS_ROOT_OVERRIDE%/}"

# Scratch for temp files (02 uses mktemp).
export TMPDIR="${TMPDIR:-/tmp}"

# Output layout.
OUTPUT_DIR="${RESULTS_ROOT%/}/${SAMPLE}/profiling"
RESULTS_DIR="${RESULTS_ROOT%/}/${SAMPLE}"
export TFP_LOG_DIR="${RESULTS_ROOT%/}/${SAMPLE}/logs"
mkdir -p "$OUTPUT_DIR" "$RESULTS_DIR" "$TFP_LOG_DIR"

# Autosomes as a bash array (safe: chromosome tokens have no spaces/specials).
read -r -a AUTOSOMES <<< "$AUTOSOMES"

echo "############################################################"
echo "# Nucleosome profiling — end to end"
echo "#   sample   : $SAMPLE"
echo "#   config   : $CONFIG"
echo "#   TF sel   : $TF_SELECTION"
echo "#   FCC dir  : $FCC_DIR"
echo "#   threads  : $THREADS"
echo "#   output   : $OUTPUT_DIR"
echo "#   features : $RESULTS_DIR/all_TF_features.tsv"
echo "############################################################"

# ------------------------------- FCC FILE HELPERS ----------------------------
fcc_path()  { echo "$FCC_DIR/${1}_counts.tsv.gz"; }
have_fcc () { ls "$FCC_DIR"/*_counts.tsv.gz >/dev/null 2>&1; }

# ------------------------------- STAGE 02: ENSURE FCC ------------------------
if have_fcc && [ "$FORCE_RECOUNT" != true ]; then
    echo "[FCC] found existing FCC in $FCC_DIR — reusing (use --force-recount to redo)."
else
    [ -f "$BAM" ]     || { echo "ERROR: BAM not found for counting: $BAM" >&2; exit 1; }
    [ -f "$GC_BIAS" ] || { echo "ERROR: GC-bias file not found: $GC_BIAS" >&2; exit 1; }
    echo "[FCC] counting from BAM -> $FCC_DIR"
    mkdir -p "$FCC_DIR"
    # 02 self-checks tools/paths and exits nonzero if any chromosome fails; set -e
    # aborts the driver on that. So a bad count can never fall through to profiling.
    bash "$COUNT_SH" \
        --bam "$BAM" --gc-bias "$GC_BIAS" --out "$FCC_DIR" \
        --ref "$REFERENCE_FA" --mappable "$MAPPABLE_BED" \
        --frag-min "$FRAG_LEN_MIN" --frag-max "$FRAG_LEN_MAX" \
        --mapq "$MAPQ" --gc-bias-min "$GC_BIAS_MIN" --threads "$THREADS"
fi

# --- completeness sweep: every autosome must have an FCC file ----------------
CHROMS_PRESENT=()
MISSING=()
for chrom in "${AUTOSOMES[@]}"; do
    if [ -f "$(fcc_path "$chrom")" ]; then
        CHROMS_PRESENT+=("$chrom")
    else
        MISSING+=("$chrom")
    fi
done
if [ "${#MISSING[@]}" -gt 0 ]; then
    echo "ERROR: FCC missing for autosome(s): ${MISSING[*]}" >&2
    echo "       (counting incomplete, or FCC_DIR points at a partial set)" >&2
    exit 1
fi
echo "[FCC] complete: ${#CHROMS_PRESENT[@]} autosomes present."

# ------------------------------- RESOLVE TF SELECTION ------------------------
# Warn on derived-name collisions (name before first '_'): two such beds MERGE
# under one TF. Checked on BOTH the ALL path and the subset path.
collision_check () {
    local dir="$1" dup
    dup=$(ls "$dir"/*.bed 2>/dev/null | xargs -n1 basename 2>/dev/null \
          | sed -E 's/_.*//; s/\.bed$//' | sort | uniq -d || true)
    [ -n "$dup" ] && {
        echo "  [WARN] TF-name collisions (will MERGE under one name): $(echo "$dup" | tr '\n' ' ')"
    }
    return 0
}

WORK="${OUTPUT_DIR}/_tf_selection"
rm -rf "$WORK"

link_bed () {  # $1 = TF name; link matching bed(s)
    local tf="$1" hits=() b
    for b in "$TF_DIR/$tf".bed "$TF_DIR/$tf"_*.bed; do
        [ -f "$b" ] && hits+=("$b")
    done
    if [ "${#hits[@]}" -eq 0 ]; then echo "  [WARN] no bed found for TF '$tf'" >&2; return 1; fi
    if [ "${#hits[@]}" -gt 1 ]; then
        echo "  [WARN] '$tf' matches ${#hits[@]} beds; keyed on name before first '_', so these MERGE: ${hits[*]##*/}" >&2
    fi
    for b in "${hits[@]}"; do ln -sf "$b" "$WORK/$(basename "$b")"; done
}

if [ "$TF_SELECTION" = "ALL" ]; then
    echo "[TF ] using ALL beds in $TF_DIR"
    USE_TF_DIR="$TF_DIR"
    collision_check "$TF_DIR"
else
    mkdir -p "$WORK"
    names=()
    if [ -f "$TF_SELECTION" ]; then
        mapfile -t names < <(grep -v '^[[:space:]]*$' "$TF_SELECTION")
        echo "[TF ] ${#names[@]} names from file $TF_SELECTION"
    else
        IFS=',' read -r -a names <<< "$TF_SELECTION"
        echo "[TF ] selecting: ${names[*]}"
    fi
    for tf in "${names[@]}"; do link_bed "$(echo "$tf" | xargs)" || true; done
    n_linked=$(ls "$WORK"/*.bed 2>/dev/null | wc -l | tr -d ' ')
    [ "$n_linked" -gt 0 ] || { echo "ERROR: no beds matched the TF selection" >&2; exit 1; }
    collision_check "$WORK"
    USE_TF_DIR="$WORK"
    echo "[TF ] $n_linked bed(s) selected"
fi

# ------------------------------- STAGE 04a: PROFILE --------------------------
# Wipe stale intermediates/composites so a re-run with a DIFFERENT TF selection
# can't leave old TFs in the merged output.
#rm -rf "$OUTPUT_DIR/intermediate" "$OUTPUT_DIR/composite"
rm -rf "$OUTPUT_DIR/intermediate"

PROF_JOBLOG="$OUTPUT_DIR/_prof_joblog.tsv"; rm -f "$PROF_JOBLOG"
echo "[PROF] profiling ${#CHROMS_PRESENT[@]} chromosomes (up to $THREADS in parallel)..."

# Proper pool: GNU parallel keeps all $THREADS busy regardless of chrom size,
# and --joblog records each job's exit status for the completeness gate below.
set +e
printf '%s\n' "${CHROMS_PRESENT[@]}" | \
    parallel --joblog "$PROF_JOBLOG" -j "$THREADS" \
        "$PY" "$PROF_PY" \
            --config "$CONFIG" \
            --fcc_file "$FCC_DIR/{}_counts.tsv.gz" \
            --tf_dir "$USE_TF_DIR" \
            --output_dir "$OUTPUT_DIR" \
            '>' "$TFP_LOG_DIR/prof_{}.out" '2>&1'
PROF_RC=$?
set -e

# Completeness gate: joblog present, no nonzero exits, one row per chromosome.
if [ ! -f "$PROF_JOBLOG" ]; then
    echo "ERROR: no profiling joblog produced" >&2; exit 1
fi
n_failed=$(awk 'NR>1 && $7 != 0' "$PROF_JOBLOG" | wc -l | tr -d ' ')
n_ran=$(awk 'NR>1' "$PROF_JOBLOG" | wc -l | tr -d ' ')
if [ "$n_failed" -gt 0 ] || [ "$PROF_RC" -ne 0 ] || [ "$n_ran" -ne "${#CHROMS_PRESENT[@]}" ]; then
    echo "ERROR: profiling incomplete ($n_failed failed, $n_ran/${#CHROMS_PRESENT[@]} ran). See:" >&2
    echo "       joblog: $PROF_JOBLOG   per-chrom logs: $TFP_LOG_DIR/prof_*.out" >&2
    awk 'NR==1 || $7 != 0' "$PROF_JOBLOG" >&2
    exit 1
fi
echo "[PROF] all ${#CHROMS_PRESENT[@]} chromosomes profiled OK."

# ------------------------------- STAGE 04b: MERGE + FEATURES -----------------
echo "[MERGE] compositing + extracting features (n_jobs=$THREADS)..."
"$PY" "$MERGE_PY" \
    --config "$CONFIG" \
    --output_dir "$OUTPUT_DIR" \
    --results_dir "$RESULTS_DIR" \
    --n_jobs "$THREADS"

# ------------------------------- DONE ----------------------------------------
N_PROF=$(find "$OUTPUT_DIR/composite" -name "*_profile.tsv.gz" 2>/dev/null | wc -l | tr -d ' ')
echo "############################################################"
echo "# DONE: $SAMPLE"
echo "#   composite profiles : $OUTPUT_DIR/composite  ($N_PROF TFs)"
echo "#   features table     : $RESULTS_DIR/all_TF_features.tsv"
echo "#   logs               : $TFP_LOG_DIR"
echo "############################################################"
