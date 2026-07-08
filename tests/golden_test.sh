#!/bin/bash
# =============================================================================
# golden_test.sh  —  gate the stage-04 changes before swapping the edited scripts in.
#
# What it does (Part A):
#   Runs your CURRENT (original) 04_* scripts AND the EDITED 04_* scripts on the
#   SAME existing FCC input (a couple of small chromosomes, a TF subset), then
#   compares. The stage-04 changes only touch n_sites, so the expectation is:
#       * composite profiles      -> IDENTICAL
#       * every feature but n_sites-> IDENTICAL
#       * n_sites                  -> changed (the fix)
#   Any other movement => FAIL (we broke something).
#
# Optional (Part B, RUN_AWK_CHECK=true):
#   A synthetic micro-check of the 02_fcc_count.sh GC-counting change (S-fix +
#   seeding) — fast, no genome-wide re-run.
#
# This test needs NO Slurm — run it on an interactive node or login shell with
# tfp_env active.  It does NOT modify your originals.
# =============================================================================

set -euo pipefail

# ------------------------------- CONFIG --------------------------------------
# Edit these paths, then run:  bash golden_test.sh
ORIG_SCRIPTS="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/scripts"   # your CURRENT 04_*.py
NEW_SCRIPTS="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/scripts_edited"  # the EDITED 04_*.py (+ compare_outputs.py)

# An existing sample's FCC directory (per-chromosome *_counts.tsv.gz already computed)
FCC_DIR="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/data/samples/SeCT-26_t2/fcc"

TF_DIR="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/data/TFBS_10000ms"
CHROM_SIZES="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/refs/hg38.standard.chrom.sizes"

WORK="./golden_test_work"          # scratch dir for this test (safe to delete)

# --- Chromosomes -------------------------------------------------------------
# Leave TEST_CHROMS empty to auto-use ALL chromosomes present in FCC_DIR.
# Or list specific ones, e.g. TEST_CHROMS=("chr21" "chr22")
TEST_CHROMS=()

# --- TF selection ------------------------------------------------------------
# Option 1: name specific TFs (matched against the START of the bed filename,
#   i.e. the same "split before first _" rule the pipeline uses). Recommended.
#   e.g. TF_NAMES=("FOXA1" "GRHL2" "SPI1")
TF_NAMES=("FOXA1" "GRHL2" "SPI1" "CDX2")
# Option 2: if TF_NAMES is empty, fall back to the first N_TF beds (0 = all).
N_TF=20

PY="python"                        # interpreter for tfp_env (or full path to it)

RUN_AWK_CHECK=true                 # set false to skip the Part B synthetic GC check

# Feature/profile tolerances. Same FCC input => stage-04 changes are exact, so 0.
FEATURE_TOL=0
PROFILE_TOL=0
# -----------------------------------------------------------------------------

echo "############################################################"
echo "# Stage-04 golden test"
echo "#   original scripts : $ORIG_SCRIPTS"
echo "#   edited scripts   : $NEW_SCRIPTS"
echo "#   FCC input        : $FCC_DIR"
echo "#   chroms           : ${TEST_CHROMS[*]}   TF subset: $N_TF"
echo "############################################################"

# --- sanity checks -----------------------------------------------------------
for f in "$ORIG_SCRIPTS/04_profile_tf_sites.py" "$ORIG_SCRIPTS/04_merge_and_extract.py" \
         "$NEW_SCRIPTS/04_profile_tf_sites.py"  "$NEW_SCRIPTS/04_merge_and_extract.py" \
         "$NEW_SCRIPTS/compare_outputs.py" "$CHROM_SIZES"; do
    [ -e "$f" ] || { echo "ERROR: missing required file: $f"; exit 1; }
done
[ -d "$FCC_DIR" ] || { echo "ERROR: FCC_DIR not found: $FCC_DIR (set it to a real sample)"; exit 1; }
[ -d "$TF_DIR" ]  || { echo "ERROR: TF_DIR not found: $TF_DIR"; exit 1; }

# --- resolve chromosomes: if TEST_CHROMS is empty, use all present in FCC_DIR
if [ "${#TEST_CHROMS[@]}" -eq 0 ]; then
    mapfile -t TEST_CHROMS < <(ls "$FCC_DIR"/*_counts.tsv.gz "$FCC_DIR"/*.tsv.gz 2>/dev/null \
        | xargs -n1 basename 2>/dev/null \
        | sed -E 's/(_counts)?\.tsv\.gz$//' \
        | sort -u -V)
    [ "${#TEST_CHROMS[@]}" -gt 0 ] || { echo "ERROR: no FCC files found in $FCC_DIR"; exit 1; }
    echo "Auto-detected ${#TEST_CHROMS[@]} chromosomes in FCC_DIR: ${TEST_CHROMS[*]}"
fi

# Best-effort: the 04_*.py scripts open a hardcoded benchmark-log path at import.
# On the cluster this dir already exists; create it defensively so import can't fail.
mkdir -p /cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/scripts/final_griffin 2>/dev/null || true

rm -rf "$WORK"
mkdir -p "$WORK"

# --- build a small TF subset for speed ---------------------------------------
SUBSET_DIR="$WORK/tf_subset"
mkdir -p "$SUBSET_DIR"
if [ "${#TF_NAMES[@]}" -gt 0 ]; then
    # select by name: match beds whose filename starts with "<TF>" (before first _ or .)
    for tf in "${TF_NAMES[@]}"; do
        found=0
        for bed in "$TF_DIR/$tf".bed "$TF_DIR/$tf"_*.bed; do
            [ -f "$bed" ] || continue
            ln -sf "$bed" "$SUBSET_DIR/$(basename "$bed")"
            found=1
        done
        [ "$found" -eq 1 ] || echo "  [WARN] no bed found for TF '$tf' in $TF_DIR"
    done
    n_found=$(ls "$SUBSET_DIR"/*.bed 2>/dev/null | wc -l)
    [ "$n_found" -gt 0 ] || { echo "ERROR: none of the named TFs matched a bed file"; exit 1; }
    echo "Using named TF subset: $n_found bed(s) -> ${TF_NAMES[*]}"
    USE_TF_DIR="$SUBSET_DIR"
elif [ "$N_TF" -gt 0 ]; then
    ls "$TF_DIR"/*.bed | sort | head -n "$N_TF" | while read -r bed; do
        ln -sf "$bed" "$SUBSET_DIR/$(basename "$bed")"
    done
    echo "Using first-$N_TF TF subset: $(ls "$SUBSET_DIR"/*.bed | wc -l) beds"
    USE_TF_DIR="$SUBSET_DIR"
else
    echo "Using ALL TFs in $TF_DIR"
    USE_TF_DIR="$TF_DIR"
fi

# --- helper: resolve an FCC file for a chromosome (handles both naming schemes)
find_fcc () {
    local chrom="$1"
    for cand in "$FCC_DIR/${chrom}_counts.tsv.gz" "$FCC_DIR/${chrom}.tsv.gz"; do
        [ -f "$cand" ] && { echo "$cand"; return 0; }
    done
    return 1
}

# --- run one version (orig|new) end to end -----------------------------------
run_version () {
    local tag="$1" scripts="$2"
    local out="$WORK/$tag/prof" res="$WORK/$tag/res"
    mkdir -p "$out" "$res"
    echo ""
    echo ">>> [$tag] profiling ${TEST_CHROMS[*]}"
    for chrom in "${TEST_CHROMS[@]}"; do
        local fcc; fcc="$(find_fcc "$chrom")" || { echo "  skip $chrom (no FCC file)"; continue; }
        echo "    $chrom  <- $fcc"
        "$PY" "$scripts/04_profile_tf_sites.py" \
            --fcc_file "$fcc" \
            --tf_dir "$USE_TF_DIR" \
            --output_dir "$out" \
            --window_size 5000 \
            --chrom_sizes "$CHROM_SIZES" \
            --step 15 >/dev/null
    done
    echo ">>> [$tag] merge + features"
    "$PY" "$scripts/04_merge_and_extract.py" \
        --output_dir "$out" \
        --results_dir "$res" >/dev/null
    echo ">>> [$tag] done: profiles in $out/composite , features in $res/all_TF_features.tsv"
}

run_version orig "$ORIG_SCRIPTS"
run_version new  "$NEW_SCRIPTS"

# --- compare -----------------------------------------------------------------
echo ""
set +e
"$PY" "$NEW_SCRIPTS/compare_outputs.py" \
    --golden-composite "$WORK/orig/prof/composite" \
    --new-composite    "$WORK/new/prof/composite" \
    --golden-features  "$WORK/orig/res/all_TF_features.tsv" \
    --new-features     "$WORK/new/res/all_TF_features.tsv" \
    --feature-tol "$FEATURE_TOL" --profile-tol "$PROFILE_TOL"
STAGE04_RC=$?
set -e

# =============================================================================
# Part B (optional): synthetic check of the 02_fcc_count.sh GC change
# =============================================================================
if [ "$RUN_AWK_CHECK" = "true" ]; then
    echo ""
    echo "############################################################"
    echo "# Part B: GC-counting micro-check (S-fix + seeding, gawk)"
    echo "############################################################"

    if ! command -v gawk >/dev/null 2>&1; then
        echo "  gawk not found — skipping (the edited 02 requires gawk for reproducible seeding)."
    else
        # gc_old = original logic under whatever /usr/bin/awk is (often mawk).
        # gc_new = edited logic (S in GC class, seeded) under gawk.
        gc_old () { awk -v s="$1" 'BEGIN{seq=s; gc=gsub(/[GCgc]/,"",seq); n=gsub(/[NRYKMBHDVnrykmbhdv]/,"",seq); gc+=int(rand()*2*n); print gc}'; }
        gc_new () { gawk -v s="$1" -v seed="$2" 'BEGIN{seq=s; gc=gsub(/[GCScgs]/,"",seq); n=gsub(/[NRYKMBHDVnrykmbhdv]/,"",seq); srand(seed); gc+=int(rand()*2*n); print gc}'; }

        echo "  seq=ATGCATGC (no S/ambiguous): old=$(gc_old ATGCATGC)  new=$(gc_new ATGCATGC 100)  (expect 4 / 4)"
        echo "  seq=ATGCS    (one S)         : old=$(gc_old ATGCS)     new=$(gc_new ATGCS 100)     (expect 2 / 3  -> S now counts)"
        echo "  seq=ATGCW    (one W)         : old=$(gc_old ATGCW)     new=$(gc_new ATGCW 100)     (expect 2 / 2  -> W stays 0)"
        a=$(gc_new ATGCN 55); b=$(gc_new ATGCN 55)
        echo "  seq=ATGCN gawk seeded x2 (55): $a / $b  (must be EQUAL -> reproducible)"
        [ "$a" = "$b" ] && echo "    -> seeding reproducible under gawk: PASS" || echo "    -> NOT reproducible: FAIL (is this really gawk?)"
    fi
fi

# --- final verdict -----------------------------------------------------------
echo ""
echo "============================================================"
if [ "$STAGE04_RC" -eq 0 ]; then
    echo "STAGE-04 GOLDEN GATE: PASS  (only n_sites moved; profiles + other features identical)"
    echo "Safe to swap in the edited 04_* scripts."
else
    echo "STAGE-04 GOLDEN GATE: FAIL  (something other than n_sites changed — investigate above)"
fi
echo "Work dir kept at: $WORK  (delete when done)"
echo "============================================================"
exit "$STAGE04_RC"
