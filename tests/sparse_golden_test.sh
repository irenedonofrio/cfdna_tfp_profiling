#!/bin/bash
# =============================================================================
# sparse_golden_test.sh  —  prove that making the FCC files SPARSE does not change
# the profiling output, before you flip 02_fcc_count.sh to sparse in production.
#
# The trick: it needs NO BAM and NO re-count. A sparse FCC file is exactly the dense
# file with the "0 0" rows removed, so we build the sparse version by filtering your
# EXISTING dense FCC, then run 04 (profile + merge) on BOTH and diff the results.
#
#   dense FCC  --04-->  golden profiles + features
#   sparse FCC --04-->  new    profiles + features
#   expect: IDENTICAL (profiles, all features, and n_sites — same sites either way)
#
# No Slurm. Run on an interactive node with tfp_env active.
# =============================================================================

set -euo pipefail

# ------------------------------- CONFIG --------------------------------------
PY="python"
SCRIPTS="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/scripts_edited"   # edited 04_*.py + compare_outputs.py
DENSE_FCC_DIR="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/data/samples/SeCT-26_t2/fcc"
TF_DIR="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/data/TFBS_10000ms"
CHROM_SIZES="/cluster/work/medinfmk/cfDNA-SeCT/irene_tfp/refs/hg38.standard.chrom.sizes"

WORK="./sparse_test_work"
TF_NAMES=("FOXA1" "GRHL2" "SPI1" "CDX2")
WINDOW=5000
STEP=15

# OPTIONAL: if you have already run the SPARSE 02 on a sample, point this at its FCC
# dir to also verify the sparse-02 output equals the dense file's non-zero rows
# byte-for-byte. Leave empty to skip that check.
SPARSE_02_DIR=""
# -----------------------------------------------------------------------------

PROF_PY="$SCRIPTS/04_profile_tf_sites.py"
MERGE_PY="$SCRIPTS/04_merge_and_extract.py"
CMP_PY="$SCRIPTS/compare_outputs.py"
for f in "$PROF_PY" "$MERGE_PY" "$CMP_PY" "$CHROM_SIZES"; do
    [ -e "$f" ] || { echo "ERROR: missing: $f"; exit 1; }
done
[ -d "$DENSE_FCC_DIR" ] || { echo "ERROR: DENSE_FCC_DIR not found: $DENSE_FCC_DIR"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK"
export TFP_LOG_DIR="$WORK/logs"; mkdir -p "$TFP_LOG_DIR"

echo "############################################################"
echo "# Sparse-FCC golden test"
echo "#   dense FCC : $DENSE_FCC_DIR"
echo "#   TFs       : ${TF_NAMES[*]}"
echo "############################################################"

# --- build the sparse FCC by dropping zero rows ------------------------------
SPARSE_FCC_DIR="$WORK/fcc_sparse"; mkdir -p "$SPARSE_FCC_DIR"
echo ""
echo "[1] building sparse FCC (dropping rows where raw==0 and fcc==0)..."
total_dense=0; total_sparse=0
for f in "$DENSE_FCC_DIR"/*_counts.tsv.gz "$DENSE_FCC_DIR"/*.tsv.gz; do
    [ -f "$f" ] || continue
    base="$(basename "$f")"
    # keep only covered positions: raw != 0 OR fcc != 0
    nd=$(zcat "$f" | wc -l)
    zcat "$f" | awk -F'\t' '($3+0)!=0 || ($4+0)!=0' | gzip > "$SPARSE_FCC_DIR/$base"
    ns=$(zcat "$SPARSE_FCC_DIR/$base" | wc -l)
    total_dense=$((total_dense + nd)); total_sparse=$((total_sparse + ns))
    printf "    %-24s  %12d -> %10d rows\n" "$base" "$nd" "$ns"
done
echo "    TOTAL rows: $total_dense -> $total_sparse"
[ "$total_sparse" -gt 0 ] || { echo "ERROR: sparse FCC is empty — is the dense FCC tab-separated with raw in col 3?"; exit 1; }
if [ "$total_dense" -gt 0 ]; then
    pct=$(awk -v a="$total_sparse" -v b="$total_dense" 'BEGIN{printf "%.1f", 100.0*a/b}')
    echo "    (sparse kept ${pct}% of rows)"
fi

# --- OPTIONAL: verify a real sparse-02 output matches dense non-zero rows -----
if [ -n "$SPARSE_02_DIR" ]; then
    echo ""
    echo "[1b] checking sparse-02 output == dense non-zero rows, byte-for-byte..."
    mism=0
    for f in "$SPARSE_FCC_DIR"/*.gz; do
        base="$(basename "$f")"
        other="$SPARSE_02_DIR/$base"
        [ -f "$other" ] || { echo "    [WARN] $base not in SPARSE_02_DIR"; continue; }
        if diff <(zcat "$f") <(zcat "$other") >/dev/null; then
            echo "    $base: identical"
        else
            echo "    $base: DIFFERS"; mism=$((mism+1))
        fi
    done
    [ "$mism" -eq 0 ] && echo "    -> sparse-02 matches dense non-zero rows: PASS" \
                      || echo "    -> MISMATCH in $mism files: FAIL"
fi

# --- build TF subset ---------------------------------------------------------
SUBSET="$WORK/tf_subset"; mkdir -p "$SUBSET"
for tf in "${TF_NAMES[@]}"; do
    for bed in "$TF_DIR/$tf".bed "$TF_DIR/$tf"_*.bed; do
        [ -f "$bed" ] && ln -sf "$bed" "$SUBSET/$(basename "$bed")"
    done
done
[ "$(ls "$SUBSET"/*.bed 2>/dev/null | wc -l)" -gt 0 ] || { echo "ERROR: no TF beds matched"; exit 1; }

# --- run 04 on a given FCC dir -----------------------------------------------
find_fcc () { for c in "$1/${2}_counts.tsv.gz" "$1/${2}.tsv.gz"; do [ -f "$c" ] && { echo "$c"; return 0; }; done; return 1; }

run04 () {
    local tag="$1" fccdir="$2"
    local out="$WORK/$tag/prof" res="$WORK/$tag/res"
    mkdir -p "$out" "$res"
    echo ""
    echo "[2:$tag] profiling from $fccdir"
    # chromosomes present in this FCC dir
    mapfile -t chroms < <(ls "$fccdir"/*_counts.tsv.gz "$fccdir"/*.tsv.gz 2>/dev/null \
        | xargs -n1 basename | sed -E 's/(_counts)?\.tsv\.gz$//' | sort -u -V)
    for chrom in "${chroms[@]}"; do
        local fcc; fcc="$(find_fcc "$fccdir" "$chrom")" || continue
        "$PY" "$PROF_PY" --fcc_file "$fcc" --tf_dir "$SUBSET" --output_dir "$out" \
            --window_size "$WINDOW" --chrom_sizes "$CHROM_SIZES" --step "$STEP" >/dev/null
    done
    "$PY" "$MERGE_PY" --output_dir "$out" --results_dir "$res" --step "$STEP" >/dev/null
    echo "[2:$tag] done"
}

run04 dense  "$DENSE_FCC_DIR"
run04 sparse "$SPARSE_FCC_DIR"

# --- compare -----------------------------------------------------------------
echo ""
set +e
"$PY" "$CMP_PY" \
    --golden-composite "$WORK/dense/prof/composite" \
    --new-composite    "$WORK/sparse/prof/composite" \
    --golden-features  "$WORK/dense/res/all_TF_features.tsv" \
    --new-features     "$WORK/sparse/res/all_TF_features.tsv" \
    --feature-tol 0 --profile-tol 0
RC=$?
set -e

echo ""
echo "============================================================"
if [ "$RC" -eq 0 ]; then
    echo "SPARSE GATE: PASS — sparse FCC gives identical profiles & features."
    echo "Safe to flip 02_fcc_count.sh to sparse. (n_sites should show 0 changed here,"
    echo "since both runs use the same 04 scripts and the same sites.)"
else
    echo "SPARSE GATE: FAIL — something changed; do NOT flip to sparse yet. Investigate above."
fi
echo "Work dir: $WORK"
echo "============================================================"
exit "$RC"
