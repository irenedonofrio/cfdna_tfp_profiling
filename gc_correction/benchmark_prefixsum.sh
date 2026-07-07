#!/usr/bin/env bash
# Benchmark + byte-identity check: prefix-sum griffin_GC_counts vs your current
# optimized griffin_GC_counts, on a REAL BAM. Both are Python, same CLI, and
# MUST produce byte-identical output (prefix-sum preserves the seeded RNG).
# Reports /usr/bin/time -v for both -> Elapsed AND peak RSS (memory tradeoff).
set -euo pipefail

# ================= EDIT THESE (literal paths) =================
PY=/cluster/customapps/medinfmk/idonofrio/mambaforge/envs/griffin/bin/python
CUR=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/tfp_profiling/gc_correction/griffin_GC_counts.py
PFX=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/tfp_profiling/gc_correction/griffin_GC_counts_prefixsum.py
BAM=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/SeCT-26_t1.sortByCoord.bam
REF=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/hg38.fa
MAPPABLE=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed
CHROM_SIZES=/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/hg38.standard.chrom.sizes
NAME=SeCT-26_t1
MAPQ=20; LO=15; HI=500; CPU=8
# =============================================================

WORK=$(mktemp -d); echo "work dir: $WORK"
[ -f "${REF}.fai" ] || { echo "MISSING ${REF}.fai -- samtools faidx $REF"; exit 1; }
HAVE_GT=0; [ -x /usr/bin/time ] && HAVE_GT=1
run(){ local tf=$1; shift
  if [ "$HAVE_GT" = 1 ]; then /usr/bin/time -v "$@" >/dev/null 2>"$tf"
  else local s=$(date +%s); "$@" >/dev/null 2>"$tf"; echo "Elapsed (wall clock) time: $(( $(date +%s)-s )) s" >>"$tf"; fi; }
rep(){ grep -E "Elapsed \(wall|Maximum resident" "$1" || cat "$1"; }

echo; echo "### current optimized (str.count fetch), CPU=$CPU ###"
run "$WORK/cur.t" "$PY" "$CUR" --bam_file "$BAM" --bam_file_name "$NAME" \
  --mappable_regions_path "$MAPPABLE" --ref_seq "$REF" --chrom_sizes "$CHROM_SIZES" \
  --out_dir "$WORK/cur" --map_q "$MAPQ" --size_range "$LO" "$HI" --CPU "$CPU"
rep "$WORK/cur.t"

echo; echo "### prefix-sum, CPU=$CPU ###"
run "$WORK/pfx.t" "$PY" "$PFX" --bam_file "$BAM" --bam_file_name "$NAME" \
  --mappable_regions_path "$MAPPABLE" --ref_seq "$REF" --chrom_sizes "$CHROM_SIZES" \
  --out_dir "$WORK/pfx" --map_q "$MAPQ" --size_range "$LO" "$HI" --CPU "$CPU"
rep "$WORK/pfx.t"

echo; echo "### byte-identity (MUST be identical) ###"
if diff -q "$WORK/cur/GC_counts/$NAME.GC_counts.txt" "$WORK/pfx/GC_counts/$NAME.GC_counts.txt" >/dev/null; then
  echo "VERDICT: BYTE-IDENTICAL ✓  (prefix-sum preserves your guarantee)"
else
  echo "VERDICT: DIFFERS ✗  -- investigate (outputs kept in $WORK)"; exit 1
fi
echo; echo "outputs in: $WORK"
