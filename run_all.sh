#!/usr/bin/env bash
# run_all.sh — MASTER: Stage 1 (GC correction) -> Stage 2 (profiling), chained.
#
# The seam: Stage 1 writes gc_correction_manifest.tsv (sample -> real gc_bias
# path). This master turns that manifest into Stage 2's sample sheet (via
# build_stage2_sheet.py), so Stage 2 reads Stage 1's ACTUAL outputs instead of
# a hand-maintained sheet. Then it loops Stage 2 (one sample per invocation).
#
# The Stage-1 sample sheet is the single source of truth for which samples exist;
# everything downstream flows from it.
#
# ENV is the caller's job (conda activate / micromamba run). This script does
# not activate anything.
#
# fcc_root is read from the Stage-2 config (deployment.fcc_root), like every
# other path — no flag. fcc_dir = <fcc_root>/<sample_id>/fcc.
#
# USAGE:
#   bash run_all.sh [--gc-config F] [--prof-config F] [--tfs SEL] [--threads N] [--skip-stage1]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python}"

# Script/driver locations — override via env if your layout differs.
GC_DRIVER="${GC_DRIVER:-$SCRIPT_DIR/gc_correction/run_gc_correction.sh}"
PROF_DRIVER="${PROF_DRIVER:-$SCRIPT_DIR/nucleosome_profiling/run_profiling.sh}"
HELPER="${HELPER:-$SCRIPT_DIR/nucleosome_profiling/config_helper.py}"
GEN="${GEN:-$SCRIPT_DIR/build_stage2_sheet.py}"

# Defaults (override via flags).
GC_CONFIG="$SCRIPT_DIR/gc_correction/config_gc.yaml"
PROF_CONFIG="$SCRIPT_DIR/nucleosome_profiling/config.yaml"
TFS="ALL"
THREADS_OVERRIDE=""
SKIP_STAGE1=false
OUT_ROOT=""

usage(){ echo "usage: $0 [--out-root DIR] [--gc-config F] [--prof-config F] [--tfs SEL] [--threads N] [--skip-stage1]"; exit 1; }
while [ $# -gt 0 ]; do case "$1" in
  --out-root)    OUT_ROOT="$2"; shift 2;;
  --gc-config)   GC_CONFIG="$2"; shift 2;;
  --prof-config) PROF_CONFIG="$2"; shift 2;;
  --tfs)         TFS="$2"; shift 2;;
  --threads)     THREADS_OVERRIDE="$2"; shift 2;;
  --skip-stage1) SKIP_STAGE1=true; shift;;
  -h|--help)     usage;;
  *) echo "ERROR: unknown arg: $1" >&2; usage;;
esac; done

# --out-root redirects ALL outputs coherently, preserving the data/ vs results/
# split:  <OUT_ROOT>/data (Stage-1 + FCC)  and  <OUT_ROOT>/results/samples (Stage-2).
# When set, the three roots are derived here and threaded to the drivers + the
# sheet generator, so the manifest and generated sheet stay in lockstep.
GC_RESULTS_ROOT=""; FCC_ROOT_OVERRIDE=""; PROF_RESULTS_ROOT=""
gc_root_arg=(); prof_root_arg=()
if [ -n "$OUT_ROOT" ]; then
  OUT_ROOT="${OUT_ROOT%/}"
  GC_RESULTS_ROOT="$OUT_ROOT/data"
  FCC_ROOT_OVERRIDE="$OUT_ROOT/data/samples"
  PROF_RESULTS_ROOT="$OUT_ROOT/results/samples"
  gc_root_arg=(--results-root "$GC_RESULTS_ROOT")
  prof_root_arg=(--results-root "$PROF_RESULTS_ROOT")
fi

# preflight
for f in "$GC_DRIVER" "$PROF_DRIVER" "$HELPER" "$GEN" "$GC_CONFIG" "$PROF_CONFIG"; do
  [ -e "$f" ] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

# -------- Stage 1 --------
if [ "$SKIP_STAGE1" = true ]; then
  echo "### STAGE 1: skipped (--skip-stage1) ###"
else
  echo "### STAGE 1: GC correction ###"
  bash "$GC_DRIVER" -c "$GC_CONFIG" "${gc_root_arg[@]}"
fi

# -------- locate Stage-1 manifest --------
# If --out-root was given, the manifest is under the overridden Stage-1 root;
# otherwise read results_root from the gc config. A RELATIVE results_root is
# anchored to the GC config's own folder (same rule the drivers use), so the
# repo can be run from anywhere.
if [ -n "$OUT_ROOT" ]; then
  MANIFEST="$GC_RESULTS_ROOT/gc_correction_manifest.tsv"
else
  MANIFEST="$("$PY" - "$GC_CONFIG" <<'PY'
import sys, yaml, os
cfg_path = sys.argv[1]
cfg_dir = os.path.dirname(os.path.abspath(cfg_path))
c = yaml.safe_load(open(cfg_path))
rr = str(c['deployment']['results_root'])
rr = rr if os.path.isabs(rr) else os.path.normpath(os.path.join(cfg_dir, rr))
print(os.path.join(rr.rstrip('/'), 'gc_correction_manifest.tsv'))
PY
)"
fi
[ -f "$MANIFEST" ] || { echo "ERROR: Stage-1 manifest not found: $MANIFEST" >&2; exit 1; }
echo "### Stage-1 manifest: $MANIFEST ###"

# -------- build Stage-2 sheet from manifest, into the path Stage-2 config expects --------
# sample_sheet + fcc_root both come from the Stage-2 config (helper is strict:
# fcc_root must exist in deployment: or this errors out clearly). Both are returned
# already anchored/absolute by config_helper.py.
eval "$("$PY" "$HELPER" config "$PROF_CONFIG" sample_sheet fcc_root)"   # -> SAMPLE_SHEET, FCC_ROOT
# --out-root overrides the config fcc_root so FCC lands under the same base
[ -n "$OUT_ROOT" ] && FCC_ROOT="$FCC_ROOT_OVERRIDE"
echo "### building Stage-2 sheet -> $SAMPLE_SHEET  (fcc_root=$FCC_ROOT) ###"
"$PY" "$GEN" --gc_manifest "$MANIFEST" --fcc_root "$FCC_ROOT" --out "$SAMPLE_SHEET"

# -------- Stage 2: one invocation per sample --------
echo "### STAGE 2: profiling ###"
th_arg=(); [ -n "$THREADS_OVERRIDE" ] && th_arg=(--threads "$THREADS_OVERRIDE")
n_ok=0; n_fail=0
while IFS=$'\t' read -r sid _rest; do
  [ -z "${sid:-}" ] && continue
  echo; echo "==== STAGE 2: $sid ===="
  if bash "$PROF_DRIVER" --sample "$sid" --config "$PROF_CONFIG" --tfs "$TFS" "${th_arg[@]}" "${prof_root_arg[@]}"; then
    n_ok=$((n_ok+1))
  else
    echo "  STAGE 2 FAILED for $sid (continuing with remaining samples)" >&2
    n_fail=$((n_fail+1))
  fi
done < <(tail -n +2 "$MANIFEST" | cut -f1)

echo; echo "### ALL DONE — Stage 2 ok:$n_ok  fail:$n_fail ###"
[ "$n_fail" -eq 0 ] || exit 1
