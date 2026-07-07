#!/usr/bin/env bash
# Parallel samtools/bedtools/gawk GC-counts. Splits the (autosome-filtered)
# mappable regions into chunks, runs the full pipe per chunk across THREADS
# cores, then merges per-chunk (length,num_GC) tallies into Griffin's grid.
# Output is identical to the single-thread gc_counts_awk.sh.
set -euo pipefail

BAM=$1; REF=$2; MAPPABLE=$3; OUT=$4
MAPQ=${5:-20}; LO=${6:-15}; HI=${7:-500}; THREADS=${8:-8}

export BAM REF MAPQ LO HI
TMPD=$(mktemp -d); trap 'rm -rf "$TMPD"' EXIT

# autosomes only (chr1..chr22), matching griffin_GC_counts.py
gawk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' "$MAPPABLE" > "$TMPD/auto.bed"

# more chunks than cores => better load balance across uneven regions
NCHUNKS=$(( THREADS * 4 ))
split -n l/"$NCHUNKS" -d "$TMPD/auto.bed" "$TMPD/chunk."
ls "$TMPD"/chunk.* > "$TMPD/chunklist"

# per-chunk worker: emits sparse "length<TAB>num_GC<TAB>count" triples
worker(){
  local chunk=$1
  samtools view -q "$MAPQ" -f 1 -F 1536 -L "$chunk" "$BAM" \
  | gawk -v lo="$LO" -v hi="$HI" 'BEGIN{OFS="\t"}
      {
        flag=$2; chrom=$3; pos=$4+0; cigar=$6; tlen=$9+0
        rev=and(flag,16)
        if (rev){ if(-tlen<lo||-tlen>hi) next } else { if(tlen<lo||tlen>hi) next }
        rstart=pos-1
        if (rev){ reflen=0;n=""
          for(i=1;i<=length(cigar);i++){ch=substr(cigar,i,1)
            if(ch~/[0-9]/)n=n ch; else{if(ch~/[MDN=X]/)reflen+=n+0;n=""}}
          fe=rstart+reflen; fs=fe+tlen
        } else { fs=rstart; fe=rstart+tlen }
        if(fs<0) next
        fl=(tlen<0)?-tlen:tlen
        print chrom, fs, fe, fl
      }' \
  | bedtools getfasta -fi "$REF" -bed - -nameOnly -tab \
  | gawk '{gc=gsub(/[CGScgs]/,"",$2); t[$1 SUBSEP gc]++}
          END{for(k in t){split(k,a,SUBSEP); print a[1]"\t"a[2]"\t"t[k]}}'
}
export -f worker

# run chunks in parallel; concatenate all sparse triples
xargs -a "$TMPD/chunklist" -P "$THREADS" -I{} bash -c 'worker "$@"' _ {} \
  > "$TMPD/partials.tsv"

# merge: sum counts per (length,num_GC), emit full Griffin grid
gawk -v lo="$LO" -v hi="$HI" -F'\t' '
  {c[$1 SUBSEP $2]+=$3}
  END{
    print "length\tnum_GC\tnumber_of_fragments"
    for(L=lo;L<=hi;L++) for(g=0;g<=L;g++) print L"\t"g"\t"(c[L SUBSEP g]+0)
  }' "$TMPD/partials.tsv" > "$OUT"
