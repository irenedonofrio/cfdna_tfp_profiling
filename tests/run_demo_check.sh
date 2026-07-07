#!/usr/bin/env bash
#
# run_demo_check.sh
# -----------------
# Verifies the optimized griffin_GC_counts.py against the GROUND-TRUTH output
# that ships with the Griffin repo demo
# (demo/.../expected_results/Healthy_demo.GC_counts.txt).
#
# It uses the exact demo parameters (map_q=20, size_range=15 500, the shipped
# k100 mappable-regions file, sample name "Healthy_demo"), converts the demo
# CRAM to BAM, runs the optimized counter, and diffs the result against the
# expected file. PASS means byte-for-byte identical.
#
# Requirements: samtools, pysam, and a UCSC hg38.fa (the assembly the demo was
# aligned to). The repo does NOT ship hg38.fa -- point REF_SEQ at your copy.
#
# Usage:
#   REF_SEQ=/path/to/hg38.fa \
#   GRIFFIN_DIR=/path/to/Griffin \
#   OPT_SCRIPT=/path/to/griffin_GC_counts.py \
#   ./run_demo_check.sh
#
set -euo pipefail

_hash() {  # portable checksum (Linux: md5sum, macOS: shasum)
    if command -v md5sum >/dev/null 2>&1; then md5sum < "$1" | cut -d' ' -f1
    elif command -v shasum >/dev/null 2>&1; then shasum -a 256 < "$1" | cut -d' ' -f1
    else cksum < "$1" | cut -d' ' -f1; fi
}

GRIFFIN_DIR="${GRIFFIN_DIR:?set GRIFFIN_DIR to the cloned Griffin repo root}"
OPT_SCRIPT="${OPT_SCRIPT:?set OPT_SCRIPT to the optimized griffin_GC_counts.py}"
REF_SEQ="${REF_SEQ:?set REF_SEQ to your UCSC hg38.fa}"
PYTHON="${PYTHON:-python3}"
CPU="${CPU:-4}"

# --- fixed demo parameters (from the repo's GC-correction config) ---
NAME="Healthy_demo"
MAP_Q=20
SIZE_RANGE="15 500"
CRAM="$GRIFFIN_DIR/demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.cram"
MAPPABLE="$GRIFFIN_DIR/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed"
CHROM_SIZES="$GRIFFIN_DIR/Ref/hg38.standard.chrom.sizes"
EXPECTED="$GRIFFIN_DIR/demo/griffin_GC_correction_demo_files/expected_results/${NAME}.GC_counts.txt"

# --- checks ---
command -v samtools >/dev/null || { echo "ERROR: samtools not found" >&2; exit 1; }
for f in "$CRAM" "$MAPPABLE" "$CHROM_SIZES" "$EXPECTED" "$REF_SEQ" "$OPT_SCRIPT"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[[ -f "${REF_SEQ}.fai" ]] || { echo "indexing reference ..."; samtools faidx "$REF_SEQ"; }

TMP="$(mktemp -d)"
echo "tmp dir: $TMP"

# --- 1) CRAM -> BAM (needs the reference) ---
echo "=== converting demo CRAM -> BAM ==="
samtools view -b -T "$REF_SEQ" -o "$TMP/${NAME}.bam" "$CRAM"
samtools index "$TMP/${NAME}.bam"

# --- 2) run the optimized counter with the exact demo parameters ---
echo "=== running optimized griffin_GC_counts.py ==="
"$PYTHON" "$OPT_SCRIPT" \
    --bam_file "$TMP/${NAME}.bam" \
    --bam_file_name "$NAME" \
    --mappable_regions_path "$MAPPABLE" \
    --ref_seq "$REF_SEQ" \
    --chrom_sizes "$CHROM_SIZES" \
    --out_dir "$TMP" \
    --map_q "$MAP_Q" \
    --size_range $SIZE_RANGE \
    --CPU "$CPU"

GOT="$TMP/GC_counts/${NAME}.GC_counts.txt"

# --- 3) compare against the shipped ground truth ---
echo
echo "=== comparing to expected ==="
echo "expected : $EXPECTED   ($(_hash "$EXPECTED"))"
echo "optimized: $GOT   ($(_hash "$GOT"))"
echo
if diff -q "$EXPECTED" "$GOT" >/dev/null; then
    echo "PASS: optimized output matches the demo expected_results byte-for-byte."
    rm -rf "$TMP"
else
    echo "FAIL: outputs differ. Showing first differences (expected | optimized):"
    diff "$EXPECTED" "$GOT" | head -20
    echo
    echo "If many rows differ, the most likely cause is a different hg38 build"
    echo "than the one the demo was aligned to. Outputs kept in $TMP."
    exit 1
fi
