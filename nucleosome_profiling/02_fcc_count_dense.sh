#!/bin/bash

# Check if the input BAM file, GC bias file, and reference genome are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.bam> <GC_bias.tsv> <output_base> [min_quality]"
    exit 1
fi

# activate env
source ~/.bashrc
conda activate tfp_env

# CHANGED: require gawk. The ambiguous-base draw below is seeded per fragment for
# reproducibility, but mawk's srand(seed) is NOT deterministic across processes
# (it mixes in entropy), so under mawk the seeding silently fails. gawk's srand(seed)
# is deterministic. If your cluster/container 'awk' is mawk, this guard makes the
# requirement explicit instead of producing quietly non-reproducible output.
command -v gawk >/dev/null 2>&1 || { echo "ERROR: gawk is required (mawk's srand(seed) is non-deterministic)."; exit 1; }

# Input BAM file, GC bias file, and reference genome
input_bam="$1"
gc_bias_file="$2"
output_base="$3"
# Minimum mapping quality (default: 30)
min_quality=${4:-30}

reference_genome_file="/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/hg38.fa"
clean_regions_bed="/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed" #this is the k100_minus_exclusion_lists.mappable_regions.hg38.bed from Griffin repo


#get sample ID from bam file
base_name=$(basename "$input_bam")
base_name="${base_name%%.*}"   # keep only text before first dot

# Output directory
output_base="${output_base%/}"
output_dir="${output_base}"
mkdir -p "$output_dir"

# Debugging: Confirm output directory creation
if [ ! -d "$output_dir" ]; then
    echo "Error: Failed to create output directory $output_dir"
    exit 1
else
    echo "Output directory: $output_dir created successfully."
fi

# Create temporary file for GC bias lookup
temp_gc_bias_file=$(mktemp)
awk '{print $1 "-" $2, $7}' "$gc_bias_file" > "$temp_gc_bias_file"

# Function to process a single chromosome
process_chromosome() {
    chrom=$1
    min_quality=$2
    
    echo "Processing chromosome: $chrom with min_quality: $min_quality"

    # Fetch the reference sequence for the chromosome
    region_seq_file=$(mktemp)
    samtools faidx "$reference_genome_file" "$chrom" > "$region_seq_file"

    # Debugging: Confirm reference sequence fetch
    if [ ! -s "$region_seq_file" ]; then
        echo "Error: Failed to fetch reference sequence for $chrom"
        rm "$region_seq_file"
        exit 1
    fi

    intermediate_file=$(mktemp)

    # Process each fragment and calculate GC correction 
    # Note: here's you are also filtering out bad regions (encode blacklist, low mappability, etc.)
    samtools view -q "$min_quality" -f 0x2 -F 0x410  -L "$clean_regions_bed" "$input_bam" "$chrom" | \
    gawk -v chrom="$chrom" -v gc_bias_file="$temp_gc_bias_file" -v ref_seq_file="$region_seq_file" '
    BEGIN {
        # Load GC bias values into an array
        while ((getline line < gc_bias_file) > 0) {
            split(line, fields, " ");
            gc_bias[fields[1]] = fields[2];
        }
        close(gc_bias_file);

        # Load the reference sequence into a single string
        while ((getline line < ref_seq_file) > 0) {
            if (line ~ /^>/) continue; # Skip FASTA headers
            ref_seq = ref_seq line;
        }
        close(ref_seq_file);
    }
    {
        fragment_length = $9;
        if (fragment_length >= 120 && fragment_length <= 200) {
            # Calculate fragment positions
            fragment_start = $4 - 1; # 0-based index
            fragment_end = fragment_start + fragment_length - 1;

            # Skip invalid indices
            if (fragment_start < 0 || fragment_end >= length(ref_seq)) next;

            # Extract reference sequence
            fragment_seq = substr(ref_seq, fragment_start + 1, fragment_length);

            # Count GC bases in the reference sequence
            # CHANGED: added S/s so S (strong = G or C) is counted as GC, matching Griffin.
            # W (weak = A or T) is deliberately NOT in this class, so it counts as 0 (also Griffin).
            gc_count = gsub(/[GCScgs]/, "", fragment_seq);
            # CHANGED: seed the RNG per fragment, keyed to fragment_start (0-based, the SAME
            # seed value Griffin uses via np.random.default_rng(reference_start)). This makes the
            # ambiguous-base draws reproducible run-to-run instead of changing every run.
            # REQUIRES gawk (guarded at the top): under mawk, srand(seed) is non-deterministic
            # and this line would silently fail to reproduce.
            # NOTE (documented later cleanup, as agreed): the draw below is still uniform
            # int(rand()*2*N), NOT a per-base coin-flip like Griffin, so it can occasionally add
            # more GC than there are ambiguous bases. Effect is tiny (ambiguous bases are rare
            # inside the mappable regions), so left as-is for now; revisit if bit-matching Griffin.
            srand(fragment_start);
            # Count ambiguous bases and add random 0 or 1 for each
            ambiguous_count = gsub(/[NRYKMBHDVnrykmbhdv]/, "", fragment_seq);
            gc_count += int(rand() * 2 * ambiguous_count); # Randomly add 0 or 1 for each ambiguous base

            # Create the key with fragment length and GC count
            key = int(fragment_length) "-" int(gc_count);

            # Check GC bias and apply only if above threshold
            if (key in gc_bias && gc_bias[key] >= 0.05) {
                gc_correction = 1.0 / gc_bias[key];

                # Calculate midpoint
                midpoint = int(($4 + $4 + fragment_length) / 2);

                # Print to intermediate file (midpoint, uncorrected count, gc-corrected count)
                print midpoint, 1, gc_correction;
            }
        }
    }' | \
    sort -k1,1n | \
    awk '{
        # Accumulate counts for each position
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

    # Fill missing positions with zeros
    awk -v chrom="$chrom" 'BEGIN {OFS="\t"; pos = 1} {
        if (NR == 1) {
            while (pos < $1) {
                print chrom, pos, 0, 0;
                pos++;
            }
        }
        while (pos < $1) {
            print chrom, pos, 0, 0;
            pos++;
        }
        print chrom, $1, $2, $3;
        pos = $1 + 1;
    } END {
        while (pos <= end) {
            print chrom, pos, 0, 0;
            pos++;
        }
    }' "$intermediate_file" | gzip > "${output_dir}/${chrom}_counts.tsv.gz"

    # Clean up temporary files
    rm "$intermediate_file"
    rm "$region_seq_file"

    echo "Finished processing chromosome: $chrom"
}

export -f process_chromosome

# Export variables to make them available in the parallel environment
export input_bam temp_gc_bias_file reference_genome_file output_dir clean_regions_bed

# Run processing in parallel for each chromosome
#samtools idxstats "$input_bam" | awk '$3 > 0 {print $1}' | \
# process only canonical chromosomes (just autosomes)
cut -f1 "$clean_regions_bed" | sort -u | \
parallel -j 11 process_chromosome {} "$min_quality"

# Cleanup temporary GC bias file
rm "$temp_gc_bias_file"

echo "Processing completed. Output written to $output_dir"
