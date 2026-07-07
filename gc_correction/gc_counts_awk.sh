#!/usr/bin/env bash
# GC counts via samtools + bedtools + gawk (Griffin-faithful reconstruction).
#
# Reproduces the original griffin_GC_counts.py logic:
#   filters : paired, MAPQ>=map_q, not duplicate, not qcfail, size in [lo,hi]
#   fragment: forward -> [ref_start, ref_start+tlen)
#             reverse -> [ref_end+tlen, ref_end)  where ref_end = ref_start + ref_len(CIGAR)
#   GC      : counted from the REFERENCE over the whole fragment span
#
# DIVERGENCE from Griffin (documented, on purpose):
#   * ambiguous IUPAC bases (N,R,Y,K,M,B,D,H,V) are NOT randomly assigned here
#     (awk cannot reproduce numpy's per-fragment default_rng). They are counted
#     as non-GC. Affects only fragments that overlap N/ambiguous reference.
#   * samtools -L emits a boundary-spanning read once even if it overlaps two
#     intervals; Griffin's per-interval loop could count it twice. (Rare.)
set -euo pipefail

BAM=$1; REF=$2; MAPPABLE=$3; OUT=$4
MAPQ=${5:-20}; LO=${6:-15}; HI=${7:-500}; THREADS=${8:-1}

TMPD=$(mktemp -d); trap 'rm -rf "$TMPD"' EXIT

# Restrict to autosomes chr1..chr22, exactly as griffin_GC_counts.py does
# (mappable_intervals[...].isin(['chr1'..'chr22'])). The anchored regex keeps
# chr1..chr9, chr10..chr19, chr20..chr22 and nothing else (not chrX/Y/contigs).
gawk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$MAPPABLE" > "$TMPD/autosomes.bed"

# 0x1 paired ; exclude 0x400 dup + 0x200 qcfail = 1536
samtools view -@ "$THREADS" -q "$MAPQ" -f 1 -F 1536 -L "$TMPD/autosomes.bed" "$BAM" \
| gawk -v lo="$LO" -v hi="$HI" 'BEGIN{OFS="\t"}
    {
      flag=$2; chrom=$3; pos=$4+0; cigar=$6; tlen=$9+0
      rev = and(flag,16)
      # size filter (matches Griffin: forward uses +tlen, reverse uses -tlen)
      if (rev) { if (-tlen < lo || -tlen > hi) next }
      else     { if ( tlen < lo ||  tlen > hi) next }
      rstart = pos - 1                      # 0-based
      if (rev) {
        # reference length consumed by CIGAR ops M,D,N,=,X
        reflen=0; n=""
        for (i=1;i<=length(cigar);i++){
          ch=substr(cigar,i,1)
          if (ch ~ /[0-9]/){ n=n ch }
          else { if (ch ~ /[MDN=X]/) reflen+=n+0; n="" }
        }
        fe = rstart + reflen
        fs = fe + tlen                      # tlen negative
      } else {
        fs = rstart
        fe = rstart + tlen
      }
      if (fs < 0) next
      fraglen = (tlen<0)? -tlen : tlen
      print chrom, fs, fe, fraglen           # BED: 0-based, end-exclusive
    }' \
| bedtools getfasta -fi "$REF" -bed - -nameOnly -tab \
| gawk -v lo="$LO" -v hi="$HI" '
    {
      s=toupper($2)
      gc=gsub(/[CGS]/,"",s)                  # gsub returns #substitutions = GC count
      tally[$1 SUBSEP gc]++
    }
    END{
      print "length\tnum_GC\tnumber_of_fragments"
      for (L=lo; L<=hi; L++)
        for (g=0; g<=L; g++)
          print L"\t"g"\t"(tally[L SUBSEP g]+0)
    }' > "$OUT"
