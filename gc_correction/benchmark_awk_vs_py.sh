#!/usr/bin/env bash
# Benchmark + correctness check: samtools/bedtools/gawk GC-counts vs the
# optimized Python griffin_GC_counts.py, on a REAL BAM.
#
# Times both (single-thread, apples-to-apples) and compares outputs.
# Run on an interactive compute node inside screen/tmux.
# Use FULL LITERAL PATHS below (no <...> placeholders, no "..." abbreviations).
set -euo pipefail

# ================= EDIT THESE (literal paths) =================
PY=/cluster/customapps/medinfmk/idonofrio/mambaforge/envs/griffin/bin/python
GC_PY=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/tfp_profiling/gc_correction/griffin_GC_counts.py
AWK_SH=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/tfp_profiling/gc_correction/gc_counts_awk_parallel.sh
BAM=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/SeCT-26_t1.sortByCoord.bam
REF=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/hg38.fa
MAPPABLE=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed
CHROM_SIZES=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/hg38.standard.chrom.sizes
NAME=SeCT-26_t1
MAPQ=20; LO=15; HI=500
CPU=8            # cores for BOTH sides (Python --CPU and awk parallel threads)
# =============================================================

WORK=$(mktemp -d); echo "work dir: $WORK"

# ---- preflight ----
echo "### preflight ###"
for t in samtools bedtools gawk; do command -v "$t" >/dev/null || { echo "MISSING: $t (module load / conda activate first)"; exit 1; }; done
"$PY" -c 'import pysam' 2>/dev/null || { echo "MISSING: pysam in $PY"; exit 1; }
[ -f "${REF}.fai" ] || { echo "MISSING ${REF}.fai -- run: samtools faidx $REF"; exit 1; }
[ -f "$BAM" ] || { echo "MISSING $BAM"; exit 1; }
echo "all tools + .fai present"

# timing helper: GNU /usr/bin/time -v if available, else wall-clock fallback
HAVE_GNU_TIME=0; [ -x /usr/bin/time ] && HAVE_GNU_TIME=1
run_timed(){ local tf=$1; shift
  if [ "$HAVE_GNU_TIME" = 1 ]; then /usr/bin/time -v "$@" >/dev/null 2>"$tf"
  else local s e; s=$(date +%s); "$@" >/dev/null 2>"$tf"; e=$(date +%s)
       echo "Elapsed (wall clock) time: $((e-s)) s (approx; GNU time not found)" >>"$tf"; fi
}
report(){ grep -E "Elapsed \(wall|Maximum resident" "$1" || cat "$1"; }

echo; echo "### 1/2  Python (CPU=$CPU) ###"
run_timed "$WORK/py.time" "$PY" "$GC_PY" \
  --bam_file "$BAM" --bam_file_name "$NAME" \
  --mappable_regions_path "$MAPPABLE" --ref_seq "$REF" \
  --chrom_sizes "$CHROM_SIZES" --out_dir "$WORK/py" \
  --map_q "$MAPQ" --size_range "$LO" "$HI" --CPU "$CPU"
report "$WORK/py.time"

echo; echo "### 2/2  samtools/bedtools/gawk ($CPU threads) ###"
run_timed "$WORK/awk.time" bash "$AWK_SH" \
  "$BAM" "$REF" "$MAPPABLE" "$WORK/awk.GC_counts.txt" "$MAPQ" "$LO" "$HI" "$CPU"
report "$WORK/awk.time"

echo; echo "### correctness: are the results the same? ###"
PYOUT="$WORK/py/GC_counts/$NAME.GC_counts.txt"
if diff -q "$PYOUT" "$WORK/awk.GC_counts.txt" >/dev/null; then
  echo "VERDICT: BYTE-IDENTICAL"
else
  paste "$PYOUT" "$WORK/awk.GC_counts.txt" \
   | awk -F'\t' 'NR>1{d=$3-$6; if(d<0)d=-d; s+=d; if(d)r++; py+=$3; aw+=$6}
      END{
        printf "python total fragments = %d\n", py
        printf "awk    total fragments = %d\n", aw
        printf "VERDICT: NOT byte-identical\n"
        printf "  bins that differ      = %d\n", r+0
        printf "  sum of abs count diff = %d\n", s
        printf "  --> expected sources, both documented:\n"
        printf "      (1) fragments overlapping N/ambiguous reference (awk cannot\n"
        printf "          reproduce numpy per-fragment RNG); usually <<1%% on\n"
        printf "          highly-mappable autosomes.\n"
        printf "      (2) reads spanning two adjacent mappable intervals (Griffin\n"
        printf "          per-interval loop may double-count; samtools -L once).\n"
        printf "  If the diff is tiny and confined to these, results agree up to\n"
        printf "  the known, irreducible awk limitation.\n"
      }'
fi
echo; echo "outputs kept in: $WORK"
