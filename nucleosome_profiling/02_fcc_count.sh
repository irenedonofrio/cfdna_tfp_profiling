#!/bin/bash
# =============================================================================
# 02_fcc_count.sh  —  genome-wide GC-corrected midpoint (FCC) counting. SPARSE.
#
# INTERFACE (named flags; the driver builds this invocation from config.yaml
# + samples.tsv, so this script never reads YAML itself):
#
#   bash 02_fcc_count.sh \
#       --bam         <input.bam> \
#       --gc-bias     <GC_bias.tsv> \
#       --out         <fcc_output_dir> \
#       --ref         <reference.fa>          # .fai MUST sit beside it \
#       --mappable    <clean_regions.bed> \
#       --frag-min    120 \
#       --frag-max    200 \
#       --mapq        30 \
#       --gc-bias-min 0.05 \
#       --threads     4
#
# STANDALONE: runnable on its own, as long as tfp_env is ACTIVE (this script no
# longer activates conda) and gawk / samtools / parallel are on PATH.
#
# OUTPUT: SPARSE per-chromosome FCC — only covered positions are written.
# 04_profile_tf_sites.py zero-inits each window and overwrites only present
# positions, so "absent" == "present with 0". Validated bit-identical downstream.
#
# EXIT CODE: nonzero if ANY chromosome job fails (checked via parallel joblog).
# =============================================================================

set -euo pipefail

usage() {
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
}

# ------------------------------- DEFAULTS ------------------------------------
input_bam=""
gc_bias_file=""
output_base=""
reference_genome_file=""
clean_regions_bed=""
frag_min=120
frag_max=200
min_quality=30
gc_bias_min=0.05
threads=4

# ------------------------------- ARGS ----------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)         input_bam="$2";             shift 2 ;;
        --gc-bias)     gc_bias_file="$2";          shift 2 ;;
        --out)         output_base="$2";           shift 2 ;;
        --ref)         reference_genome_file="$2"; shift 2 ;;
        --mappable)    clean_regions_bed="$2";     shift 2 ;;
        --frag-min)    frag_min="$2";              shift 2 ;;
        --frag-max)    frag_max="$2";              shift 2 ;;
        --mapq)        min_quality="$2";           shift 2 ;;
        --gc-bias-min) gc_bias_min="$2";           shift 2 ;;
        --threads)     threads="$2";               shift 2 ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ------------------------------- VALIDATION ----------------------------------
for pair in "--bam:$input_bam" "--gc-bias:$gc_bias_file" "--out:$output_base" \
            "--ref:$reference_genome_file" "--mappable:$clean_regions_bed"; do
    flag="${pair%%:*}"; val="${pair#*:}"
    [ -n "$val" ] || { echo "ERROR: $flag is required" >&2; usage; exit 1; }
done
[ -f "$input_bam" ]             || { echo "ERROR: BAM not found: $input_bam" >&2; exit 1; }
[ -f "$gc_bias_file" ]          || { echo "ERROR: GC-bias file not found: $gc_bias_file" >&2; exit 1; }
[ -f "$reference_genome_file" ] || { echo "ERROR: reference not found: $reference_genome_file" >&2; exit 1; }
[ -f "${reference_genome_file}.fai" ] || { echo "ERROR: reference index missing: ${reference_genome_file}.fai (build it beside the .fa; a read-only mount cannot)" >&2; exit 1; }
[ -f "$clean_regions_bed" ]     || { echo "ERROR: mappable BED not found: $clean_regions_bed" >&2; exit 1; }

# gawk in particular: the ambiguous-base draw is seeded per fragment; mawk's
# srand(seed) is NOT deterministic across processes, so seeding silently fails
# under mawk. gawk's is deterministic.
for tool in gawk samtools parallel; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' is required but not on PATH" >&2; exit 1; }
done

output_base="${output_base%/}"
output_dir="${output_base}"
mkdir -p "$output_dir"
[ -d "$output_dir" ] || { echo "ERROR: failed to create output directory $output_dir" >&2; exit 1; }
echo "Output directory: $output_dir"

# GC bias lookup table (respects TMPDIR — set it in the container)
temp_gc_bias_file=$(mktemp)
awk '{print $1 "-" $2, $7}' "$gc_bias_file" > "$temp_gc_bias_file"

# ------------------------------- WORKER --------------------------------------
process_chromosome() {
    chrom=$1

    echo "Processing chromosome: $chrom (mapq>=$min_quality, frag ${frag_min}-${frag_max})"

    region_seq_file=$(mktemp)
    samtools faidx "$reference_genome_file" "$chrom" > "$region_seq_file"
    if [ ! -s "$region_seq_file" ]; then
        echo "Error: failed to fetch reference sequence for $chrom" >&2
        rm -f "$region_seq_file"
        exit 1
    fi

    intermediate_file=$(mktemp)

    # -L pre-filters to mappable regions (replaces Griffin's post-hoc exclusion).
    samtools view -q "$min_quality" -f 0x2 -F 0x410 -L "$clean_regions_bed" "$input_bam" "$chrom" | \
    gawk -v chrom="$chrom" \
         -v gc_bias_file="$temp_gc_bias_file" \
         -v ref_seq_file="$region_seq_file" \
         -v frag_min="$frag_min" \
         -v frag_max="$frag_max" \
         -v gc_bias_min="$gc_bias_min" '
    BEGIN {
        while ((getline line < gc_bias_file) > 0) {
            split(line, fields, " ");
            gc_bias[fields[1]] = fields[2];
        }
        close(gc_bias_file);

        while ((getline line < ref_seq_file) > 0) {
            if (line ~ /^>/) continue;
            ref_seq = ref_seq line;
        }
        close(ref_seq_file);
    }
    {
        fragment_length = $9;
        if (fragment_length >= frag_min && fragment_length <= frag_max) {
            fragment_start = $4 - 1;                       # 0-based
            fragment_end = fragment_start + fragment_length - 1;

            if (fragment_start < 0 || fragment_end >= length(ref_seq)) next;

            fragment_seq = substr(ref_seq, fragment_start + 1, fragment_length);

            # S (strong = G or C) counts as GC, matching Griffin. W (weak = A or T)
            # is deliberately NOT in this class -> counts as 0.
            gc_count = gsub(/[GCScgs]/, "", fragment_seq);

            # Seed per fragment, keyed to fragment_start (the SAME seed Griffin uses
            # via np.random.default_rng(reference_start)). REQUIRES gawk (guarded).
            # NOTE: draw is still uniform int(rand()*2*N), not a per-base coin-flip;
            # tiny effect (ambiguous bases rare in mappable regions). Left as-is.
            srand(fragment_start);
            ambiguous_count = gsub(/[NRYKMBHDVnrykmbhdv]/, "", fragment_seq);
            gc_count += int(rand() * 2 * ambiguous_count);

            key = int(fragment_length) "-" int(gc_count);

            if (key in gc_bias && gc_bias[key] >= gc_bias_min) {
                gc_correction = 1.0 / gc_bias[key];
                midpoint = int(($4 + $4 + fragment_length) / 2);
                print midpoint, 1, gc_correction;
            }
        }
    }' | \
    sort -k1,1n | \
    awk '{
        if ($1 == last_pos) {
            uncorrected_count += $2;
            corrected_count += $3;
        } else {
            if (NR > 1) {
                print last_pos, uncorrected_count, corrected_count;
            }
            last_pos = $1;
            uncorrected_count = $2;
            corrected_count = $3;
        }
    }
    END {
        if (NR > 0) {
            print last_pos, uncorrected_count, corrected_count;
        }
    }' > "$intermediate_file"

    # SPARSE output: write ONLY covered positions. 04_profile_tf_sites.py zero-inits
    # each window and overwrites only present positions (searchsorted), so "absent"
    # == "present with 0" to the reader. Non-zero rows are byte-identical to the old
    # dense output's non-zero rows, so composites/features are unchanged (validated).
    # This also removes the old dense zero-fill, whose END{} referenced an undefined
    # 'end' and never ran (dead code).
    awk -v chrom="$chrom" 'BEGIN {OFS="\t"} {
        print chrom, $1, $2, $3;
    }' "$intermediate_file" | gzip > "${output_dir}/${chrom}_counts.tsv.gz"

    rm -f "$intermediate_file" "$region_seq_file"
    echo "Finished processing chromosome: $chrom"
}

export -f process_chromosome
export input_bam temp_gc_bias_file reference_genome_file output_dir clean_regions_bed \
       min_quality frag_min frag_max gc_bias_min

# ------------------------------- RUN -----------------------------------------
# Chromosomes = those present in the mappable BED (unchanged from validated version).
joblog="${output_dir}/_joblog.tsv"
rm -f "$joblog"

set +e
cut -f1 "$clean_regions_bed" | sort -u | \
    parallel --joblog "$joblog" -j "$threads" process_chromosome {}
parallel_rc=$?
set -e

rm -f "$temp_gc_bias_file"

# ------------------------------- COMPLETENESS GATE ---------------------------
# Trust the joblog: fail loudly if any chromosome job returned nonzero.
if [ ! -f "$joblog" ]; then
    echo "ERROR: no joblog produced ($joblog) — cannot verify completeness" >&2
    exit 1
fi

n_failed=$(awk 'NR>1 && $7 != 0' "$joblog" | wc -l | tr -d ' ')
if [ "$n_failed" -gt 0 ] || [ "$parallel_rc" -ne 0 ]; then
    echo "ERROR: $n_failed chromosome job(s) failed. See joblog: $joblog" >&2
    awk 'NR==1 || $7 != 0' "$joblog" >&2
    exit 1
fi

echo "Processing completed. Output written to $output_dir"
echo "Joblog: $joblog"
