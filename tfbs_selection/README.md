# Transcription factor site selection

This folder builds two things the rest of the pipeline depends on:

1. the **list of transcription factors (TFs)** we profile, and
2. one **binding-site file per TF**, holding the genomic positions where that TF binds.

---

We start from a large public catalogue of where proteins bind DNA (GTRD). Not
everything in that catalogue is a real transcription factor, so we keep only the
entries that a curated reference (CIS-BP) recognises as genuine, sequence-specific
DNA binders (this logic follows the same logic applied by De Sarkar (2023). We then drop any TF that doesn't have enough binding sites (10000 as in Doebley (2022), N.B. in De Sarkar (2023) only 1000 were used instead)  to give a
stable signal, and for each TF that survives we keep only its strongest sites (the most supported by CHIP-seq experiment).
The output is a tidy set of per-TF site files, plus two bookkeeping files that
record exactly what was kept and what was thrown away.

---

## What goes in

**GTRD — meta-clusters (version 20.06).**
The file was downloaded from the GTRD site (http://gtrd.biouml.org:8888/downloads/20.06/intervals/chip-seq/Homo_sapiens_meta_clusters.zip).
GTRD collects thousands of ChIP-seq experiments (experiments that show where a
given protein sits on the genome) and merges them into "meta-clusters" — one
row per binding site, with a `peak.count` telling us how many experiments
supported that site. 

**CIS-BP (v2.00) — the "is this really a TF?" reference.**
GTRD will happily list any protein that was pulled down in a ChIP-seq
experiment, whether or not it truly binds a specific DNA sequence. CIS-BP is a
curated catalogue of TFs and their binding motifs, and we use it as the gate
that decides which GTRD entries count as real transcription factors.
This is where the data was dowloaded: https://cisbp.ccbr.utoronto.ca/bulk.php
---

## Workflow

**1. Read the sites.**
For each row we take the chromosome, start, end, the TF name, and the
`peak.count`. Columns are found *by their header name*, not by position — GTRD
reshuffled the column order between versions, so reading by name keeps this
robust across 19.10 / 20.06.

**2. Keep only the autosomes (chr1–22).**
We restrict to chr1–22, matching the reference methods
(Griffin, De Sarkar).

**3. The filter**
For every protein, we ask CIS-BP a simple question: *does this have a known
binding motif?* (In CIS-BP terms, `TF_Status` is `D` = motif directly measured,
`I` = motif inferred from a close relative, or `N` = no motif. We keep `D` and
`I`, and drop `N`.)

Why this matters: a lot of proteins land on chromatin without binding a specific
DNA sequence themselves — they're carried there by real TFs. Think of
co-factors and chromatin machinery like **HDAC6** or **CHD1**: they show up in
ChIP-seq, but they aren't transcription factors.

**4. Require at least 10,000 sites.**
The signal we ultimately measure is an *average* over many binding sites. We keep a TF
only if it has at least 10,000 binding sites on the autosomes — the same
threshold Doebley et al used (while De Sarkar used a lower 1,000-site threshold)

**5. Keep each TF's strongest sites.**
For every surviving TF we rank its sites by `peak.count` (how many experiments
backed each one) and keep the top 10,000. More experimental support means more
confidence the site is real, so we profile the best-evidenced ones.

**6. Pin down the exact position.**
Each site's anchor point is the midpoint of its start and end — the single
coordinate the downstream profiling centres its window on.

---

## What comes out

Inside the output folder:

- **`sites/<TF>.bed`** — one file per kept TF, listing its selected binding
  sites (chromosome, start, end, midpoint, `peak.count`). These feed straight
  into the GC-correction stage.
- **`tf_list.tsv`** — the roster. One row per TF, with how many autosomal sites
  it had, how many were written, and its status (`kept`, or the reason it was
  dropped). The TFs that actually made it are the ones marked `kept`.
- **`dropped_targets.tsv`** — the audit trail: every GTRD entry that failed the
  real-TF gate, with its site count and the reason. This is where the
  co-factors and non-TFs end up.

A quick way to read the final roster size:

```bash
awk -F'\t' '$4=="kept"' tf_list.tsv | wc -l
```

---

## Running it

```bash
python build_tf_sites.py \
    --metaclusters Homo_sapiens_meta_clusters.interval \
    --census-type cisbp --cisbp TF_Information.txt \
    --outdir ./tf_sites --min-sites 10000 --top-n 10000
```

When it starts, it prints the column mapping it inferred (e.g.
`tf='tfTitle' peakcount='peak.count'`) plus a few sample rows — worth a glance
to confirm it locked onto the right columns before trusting the output.

---

## Choices worth knowing about

**GTRD version (20.06).**
We use the current GTRD build. Griffin and De Sarkar used the older 19.10 version; the roster comes out a little larger
than their 377 / 338 — that's expected from the extra data, not a mistake.


**Which "real-TF" reference (CIS-BP vs Lambert).**
We gate on CIS-BP because that's what Griffin and De Sarkar used, which keeps us
comparable to them. An alternative is the Lambert 2018 human-TF census, which
judges TFs by whether they have a DNA-binding domain rather than by whether a
motif exists. The script supports both (`--census-type cisbp | lambert |
intersection`); the choice changes *which TFs* are on the list but never the
sites within a TF that passes.

**No mappability filter.**
De Sarkar additionally removed sites in poorly-mappable regions (Doebley et al did not). This step was skipped here.

**A note on gene names.**
The gate matches GTRD's TF names against CIS-BP's, so both sides need to agree
on the symbol. Most do. Occasionally a gene shows up under an old name — for
example brachyury appears as `T` (its legacy symbol) rather than `TBXT` — but
CIS-BP still lists it under `T` with a directly-measured motif, so it passes
correctly. If you ever see a surprising name in the roster, check it with an
exact-field match rather than a loose text search:

```bash
awk -F'\t' 'toupper($7)=="T"' TF_Information.txt
```


The BED files once created can then be moved into the `../refs/TFBS_10000ms/` folder for use in the profiling stage.

