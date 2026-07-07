#!/usr/bin/env bash
# run_gc_correction.sh — Stage 1 driver (GC correction).
#
# For each sample in the sample sheet: run GC_counts (prefix-sum) then GC_bias,
# writing <results_root>/GC_counts/<name>.GC_counts.txt and
#         <results_root>/GC_bias/<name>.GC_bias.txt
# and recording the produced paths in <results_root>/gc_correction_manifest.tsv
# (the Stage 1 -> Stage 2 seam).
#
# NO Slurm, NO Snakemake: samples run sequentially in this process.
# The ENVIRONMENT IS THE CALLER'S JOB. This driver does not `conda activate`
# anything; it preflight-checks the env and fails loudly if it's wrong. Run it as:
#     conda activate griffin && bash run_gc_correction.sh          (interactive)
#     micromamba run -n griffin bash run_gc_correction.sh          (container)
#
# name convention (option a): bam_file_name = BAM basename minus '.bam'
#     e.g. SeCT-26_t2.sortByCoord.bam -> SeCT-26_t2.sortByCoord
set -euo pipefail

# ---------------- args ----------------
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CONFIG="$SCRIPT_DIR/config_gc.yaml"
FORCE_PLOTS=0
RESUME=1
RESULTS_ROOT_OVERRIDE=""
usage(){ echo "usage: $0 [-c config_gc.yaml] [--plots] [--no-resume] [--results-root DIR]"; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    -c|--config)      CONFIG="$2"; shift 2;;
    --plots)          FORCE_PLOTS=1; shift;;
    --no-resume)      RESUME=0; shift;;
    --results-root)   RESULTS_ROOT_OVERRIDE="$2"; shift 2;;
    -h|--help)        usage;;
    *) echo "unknown arg: $1"; usage;;
  esac
done

COUNTS="$SCRIPT_DIR/griffin_GC_counts_prefixsum.py"
BIAS="$SCRIPT_DIR/griffin_GC_bias.py"

# ---------------- read config (YAML via PyYAML) ----------------
[ -f "$CONFIG" ] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
CFG_VARS=$(python - "$CONFIG" <<'PY'
import sys, os, shlex, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
d, m = cfg['deployment'], cfg['method']
root = str(d['refs_root']).rstrip('/')
cfg_dir = os.path.dirname(os.path.abspath(sys.argv[1]))
def rel(p): return p if os.path.isabs(p) else os.path.join(root, p)
def sheet(p): return p if os.path.isabs(p) else os.path.join(cfg_dir, p)
def q(k, v): print(f'{k}={shlex.quote(str(v))}')
q('REF_FA',         rel(d['reference_fa']))
q('MAPPABLE_BED',   rel(d['mappable_bed']))
q('CHROM_SIZES',    rel(d['chrom_sizes']))
q('GENOME_GC_FREQ', rel(d['genome_gc_freq']))
q('SAMPLE_SHEET',   sheet(d['sample_sheet']))
q('RESULTS_ROOT',   str(d['results_root']).rstrip('/'))
q('THREADS',        d.get('threads', 4))
q('MAP_Q',          m['map_q'])
sr = m['size_range']; q('SIZE_LO', sr[0]); q('SIZE_HI', sr[1])
q('CFG_PLOTS',      1 if m.get('plots', False) else 0)
PY
) || { echo "ERROR: could not parse config (is PyYAML available in the env?)"; exit 1; }
eval "$CFG_VARS"

# CLI override beats config (keeps all Stage-1 outputs + manifest under one root)
[ -n "$RESULTS_ROOT_OVERRIDE" ] && RESULTS_ROOT="${RESULTS_ROOT_OVERRIDE%/}"

# mappable_name = mappable bed basename minus .bed  (derived, not a config key)
MAPPABLE_NAME=$(basename "$MAPPABLE_BED"); MAPPABLE_NAME=${MAPPABLE_NAME%.bed}
# plots: config default, overridable with --plots
PLOTS=$CFG_PLOTS; [ "$FORCE_PLOTS" = 1 ] && PLOTS=1

# ---------------- preflight (verify env; do NOT activate it) ----------------
echo "### preflight ###"
command -v python >/dev/null || { echo "ERROR: no python on PATH — activate your env first"; exit 1; }
python -c 'import numpy,pandas,pysam,yaml' 2>/dev/null \
  || { echo "ERROR: env missing numpy/pandas/pysam/yaml — activate the griffin env first"; exit 1; }
[ -f "$COUNTS" ] || { echo "ERROR: missing $COUNTS"; exit 1; }
[ -f "$BIAS" ]   || { echo "ERROR: missing $BIAS"; exit 1; }
for f in "$REF_FA" "$REF_FA.fai" "$MAPPABLE_BED" "$CHROM_SIZES"; do
  [ -f "$f" ] || { echo "ERROR: missing reference file: $f"; exit 1; }
done
[ -d "$GENOME_GC_FREQ" ] || { echo "ERROR: missing genome_GC_frequency dir: $GENOME_GC_FREQ"; exit 1; }
[ -f "$SAMPLE_SHEET" ]   || { echo "ERROR: missing sample sheet: $SAMPLE_SHEET"; exit 1; }
echo "  env + refs OK"
echo "  refs:        $REF_FA"
echo "  mappable:    $MAPPABLE_NAME"
echo "  gc_freq:     $GENOME_GC_FREQ"
echo "  results:     $RESULTS_ROOT"
echo "  params:      map_q=$MAP_Q  size_range=$SIZE_LO $SIZE_HI  threads=$THREADS  plots=$PLOTS  resume=$RESUME"

# ---------------- setup outputs + manifest ----------------
mkdir -p "$RESULTS_ROOT/GC_counts" "$RESULTS_ROOT/GC_bias" "$RESULTS_ROOT/logs"
MANIFEST="$RESULTS_ROOT/gc_correction_manifest.tsv"
printf 'sample_id\tbam\tgc_counts\tgc_bias\n' > "$MANIFEST"   # rewritten fresh each run

# ---------------- sample loop ----------------
n_ok=0; n_fail=0
# read sheet: expects a header line, columns: sample_id  bam  [extra cols ignored]
{
  read -r _header || true
  while IFS=$'\t' read -r sample_id bam _rest || [ -n "${sample_id:-}" ]; do
    sample_id=${sample_id%$'\r'}; bam=${bam%$'\r'}
    [ -z "${sample_id:-}" ] && continue
    case "$sample_id" in \#*) continue;; esac        # allow comment lines

    name=$(basename "$bam"); name=${name%.bam}         # option (a)
    counts_out="$RESULTS_ROOT/GC_counts/$name.GC_counts.txt"
    bias_out="$RESULTS_ROOT/GC_bias/$name.GC_bias.txt"
    log="$RESULTS_ROOT/logs/$name.log"

    echo; echo "==== $sample_id  ($name) ===="
    if [ ! -f "$bam" ]; then echo "  SKIP: bam not found: $bam"; n_fail=$((n_fail+1)); continue; fi

    # ---- GC_counts ----
    if [ "$RESUME" = 1 ] && [ -s "$counts_out" ]; then
      echo "  [counts] resume: exists, skipping"
    else
      echo "  [counts] running (prefix-sum, $THREADS threads) ..."
      if ! python "$COUNTS" \
            --bam_file "$bam" --bam_file_name "$name" \
            --mappable_regions_path "$MAPPABLE_BED" --ref_seq "$REF_FA" \
            --chrom_sizes "$CHROM_SIZES" --out_dir "$RESULTS_ROOT" \
            --map_q "$MAP_Q" --size_range "$SIZE_LO" "$SIZE_HI" --CPU "$THREADS" \
            > "$log" 2>&1; then
        echo "  [counts] FAILED — see $log"; n_fail=$((n_fail+1)); continue
      fi
    fi

    # ---- GC_bias ----
    if [ "$RESUME" = 1 ] && [ -s "$bias_out" ]; then
      echo "  [bias] resume: exists, skipping"
    else
      plots_flag=""; [ "$PLOTS" = 0 ] && plots_flag="--no_plots"
      echo "  [bias] running${plots_flag:+ (}${plots_flag}${plots_flag:+)} ..."
      if ! python "$BIAS" \
            --bam_file_name "$name" \
            --mappable_name "$MAPPABLE_NAME" \
            --genome_GC_frequency "$GENOME_GC_FREQ" \
            --out_dir "$RESULTS_ROOT" \
            --size_range "$SIZE_LO" "$SIZE_HI" \
            $plots_flag \
            >> "$log" 2>&1; then
        echo "  [bias] FAILED — see $log"; n_fail=$((n_fail+1)); continue
      fi
    fi

    printf '%s\t%s\t%s\t%s\n' "$sample_id" "$bam" "$counts_out" "$bias_out" >> "$MANIFEST"
    echo "  done -> $bias_out"
    n_ok=$((n_ok+1))
  done
} < "$SAMPLE_SHEET"

echo; echo "### summary ###"
echo "  ok:   $n_ok"
echo "  fail: $n_fail"
echo "  manifest: $MANIFEST"
[ "$n_fail" -eq 0 ] || exit 1
