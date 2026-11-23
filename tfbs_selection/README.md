# Transcription factor site selection

This folder builds two things the rest of the pipeline depends on:

1. the **list of transcription factors (TFs)** we profile, and
2. one **binding-site file per TF**, holding the genomic positions where that TF binds.

---

We start from a large public catalogue of where proteins bind DNA (GTRD). Not
everything in that catalogue is a real transcription factor, so we keep only the
entries that a curated reference (CIS-BP) lists with a known or inferred
DNA-binding motif. This is the same gate used by Griffin (Doebley et al., 2022)
and De Sarkar et al. (2023). We then drop any TF with fewer than 10,000 binding
sites on the autosomes (the threshold both papers used) and, for each TF that
survives, keep its 10,000 best-supported sites (the ones backed by the most
ChIP-seq evidence). The output is one site file per TF, plus bookkeeping files
that record what was kept, what was dropped, and what the input looked like.

---

## What goes in

**GTRD meta-clusters (version 20.06).**
Downloaded from
http://gtrd.biouml.org:8888/downloads/20.06/intervals/chip-seq/Homo_sapiens_meta_clusters.zip
(GTRD 20.06 is also served from gtrd20-06.biouml.org). Keep the zip exactly as
downloaded and pass it to the script: **do not extract it**. An earlier build
used a partial extraction that stopped in the middle of chr22, so every TF
silently lost all of its sites on chr3–chr9. The script now refuses incomplete
input (step 4 below).

GTRD collects thousands of ChIP-seq experiments (experiments that show where a
given protein sits on the genome) and merges their peaks into "meta-clusters":
one row per binding site. Two columns matter here:

- `peak.count`: how many individual peaks were merged into the site, counted
  over all experiments **and** all peak callers (GTRD runs several callers on
  each experiment). It measures support but is not the number of experiments.
  For example, the top CTCF site in GTRD 19.10 (as shipped with Griffin's demo
  files) has peaks from 541 experiments but a `peak.count` of 3,138.
- `exp.set`: the IDs of the experiments that support the site. Their number is
  used as the second ranking key. Tables with an `exp.count` column instead
  (as in the 19.10-derived site lists shipped with Griffin) also work.

**CIS-BP: the "is this really a TF?" reference.**
GTRD will happily list any protein that was pulled down in a ChIP-seq
experiment, whether or not it truly binds a specific DNA sequence. CIS-BP is a
curated catalogue of TFs and their binding motifs, and we use it as the gate
that decides which GTRD entries count as real transcription factors.
The human `TF_Information.txt` was downloaded from
https://cisbp.ccbr.utoronto.ca/bulk.php in July 2026
(`Homo_sapiens_2026_07_09_8_33_am/`). This is CIS-BP Database Build 3.10
(April 2026, as stated on the CIS-BP home page; every ID in the file carries
the `_3.10` suffix), a newer build than the v2.00 used by Griffin and
De Sarkar.

---

## Workflow

**1. Read the sites.**
The script streams the table straight out of the zip. For each row it takes the
chromosome, start, end, TF name, `peak.count` and the experiment set. Columns
are found *by their header name*, not by position: GTRD reshuffled the column
order between versions, so reading by name keeps this robust across 19.10 and
20.06.

**2. Keep only the autosomes (chr1–22).**
We restrict to chr1–22, matching the reference methods (Griffin, De Sarkar).

**3. The filter.**
For every protein, we ask CIS-BP a simple question: *does this have a known
binding motif?* In CIS-BP terms, `TF_Status` is `D` (motif directly measured),
`I` (motif inferred from a close relative) or `N` (no motif). We keep `D` and
`I`, and drop `N`.

Why this matters: a lot of proteins land on chromatin without binding a specific
DNA sequence themselves; they're carried there by real TFs. Think of
co-factors and chromatin machinery like **HDAC6** or **CHD1**: they show up in
ChIP-seq, but they aren't transcription factors. Every GTRD target removed here
is listed in `dropped_targets.tsv`.

**4. Check that the input is complete.**
Before anything is written, the script checks that:
- every autosome has sites;
- the furthest site on each autosome reaches at least 90% of its hg38 length;
- the file ends on a complete line.

If any check fails it stops with an error that names the problem. You can
override it with `--allow-incomplete`, but don't. The per-chromosome numbers are
written to `chrom_summary.tsv`. A corrupt or partially downloaded zip also stops
the run, because the zip's checksum is verified as it is read.

**5. Require at least 10,000 sites.**
The signal we ultimately measure is an *average* over many binding sites. We
keep a TF only if it has at least 10,000 binding sites on the autosomes, the
same threshold as Doebley et al. (2022) and De Sarkar et al. (2023).

**6. Keep each TF's best-supported sites.**
For every surviving TF we rank its sites and keep the top 10,000. The ranking
uses three keys in order:
1. `peak.count`, highest first;
2. the number of supporting experiments, highest first;
3. a hash of the TF name and site coordinates (seed set with `--seed`, default 0).

Many sites share the same `peak.count` and number of experiments, so the cutoff
often falls inside a large tie. The hash makes the choice among tied sites
reproducible and independent of the order of the input file. An earlier version
kept tied sites in file order instead. GTRD is sorted by chromosome, so that
favoured chr1, chr10, chr11 and so on; in that build, IRF9 had 9,213 of its
10,000 sites on chr1. `tf_list.tsv` reports, for each TF, how many kept sites
were chosen by the hash (`n_cutoff_kept`) and how many candidates they were
chosen from (`n_cutoff_total`).

**7. Pin down the exact position.**
Each site's anchor point is the midpoint of its start and end. This is the
single coordinate the profiling stage centres its window on.

---

## What comes out

Inside the output folder:

- **`sites/<TF>.bed`**: one file per kept TF, best site first. BED6 with no
  header: chromosome, midpoint, midpoint + 1, TF, `peak.count` (in the score
  column), `.` (strand). Keep it at exactly six columns: the profiling stage
  reads the file with six column names and fails on any extra column.
- **`tf_list.tsv`**: the roster. One row per TF that passed the CIS-BP gate:
  - `n_autosomal_sites`, `n_written` and `status` (`kept`, or the reason it was
    dropped);
  - for kept TFs, the number of chromosomes among the written sites
    (`n_chroms_written`);
  - for kept TFs, the cutoff diagnostics described in step 6
    (`cutoff_peak_count`, `cutoff_n_exp`, `n_cutoff_kept`, `n_cutoff_total`).
- **`dropped_targets.tsv`**: the audit trail. Every GTRD target that failed the
  real-TF gate, with its site count. This is where the co-factors and non-TFs
  end up.
- **`chrom_summary.tsv`**: per chromosome, the number of input rows, the
  furthest site and how far along the chromosome it lies, and the number of
  kept sites.

A quick way to read the final roster size:

```bash
awk -F'\t' '$4=="kept"' tf_list.tsv | wc -l
```

The end of the run also prints a one-line summary: GTRD targets → failed the
TF gate → passed → fewer than 10,000 sites → kept.

**Current build** (GTRD 20.06, CIS-BP 3.10): 110.5 M rows read; 1,358 GTRD
targets on the autosomes → 559 fail the TF gate → 799 pass → 369 have fewer
than 10,000 sites → **430 TFs kept** (4.3 M sites). All 430 have sites on all
22 autosomes. For every TF the cutoff falls inside a tie; 22.5% of the written
sites were chosen by the hash among equally supported candidates.
The 2,887 rows reported as skipped have no TF name: they are IgG control
experiments.

---

## Running it

```bash
python build_tf_sites.py \
    --metaclusters Homo_sapiens_meta_clusters.zip \
    --census-type cisbp --cisbp Homo_sapiens_2026_07_09_8_33_am/TF_Information.txt \
    --outdir ./tf_sites --min-sites 10000 --top-n 10000
```

When it starts, it prints the column mapping it inferred (e.g.
`tf='tfTitle' peakcount='peak.count' n_exp='exp.set'`) plus a few sample rows.
Glance at these to confirm it locked onto the right columns before trusting
the output. It then prints the result of the completeness check.

Memory stays modest: for each TF the script holds only its current top 10,000
candidates, not every site.

The BED files can then be copied into `../refs/TFBS_10000ms/` for use in the
profiling stage.

---

## Choices worth knowing about

**GTRD version (20.06).**
We use the most recent GTRD build. Griffin and De Sarkar used the older 19.10,
so our roster comes out larger: 430 TFs, versus their 377 and 338. That's expected from
the extra data, not a mistake.

**How this compares with Griffin and De Sarkar.**
- *Griffin (Doebley et al., 2022)*: GTRD 19.10, CIS-BP gate, at least 10,000
  autosomal sites (377 TFs), top 10,000 sites by `peak.count`, midpoint
  position. Our method is the same apart from the GTRD version and the
  tie-breaking.
- *De Sarkar et al. (2023)*: GTRD 19.10 and the same CIS-BP gate. Before
  counting sites, it first removes poorly mappable sites (mean hg38 Umap
  mappability over ±5 kb must be above 0.95). It then requires at least 10,000
  of the remaining autosomal sites (338 TFs) and profiles the top **1,000**
  sites per TF. The threshold is 10,000 in both papers; 1,000 is only how many
  sites De Sarkar profiled. The mappability filter applied before counting is
  the main reason their roster is shorter.

**No site-level mappability filter.**
As in Griffin, sites are not removed for low mappability. Instead, the profiling
stage only counts reads inside Griffin's mappable regions
(`k100_minus_exclusion_lists.mappable_regions.hg38.bed`).

**Which "real-TF" reference (CIS-BP vs Lambert).**
We gate on CIS-BP because that's what Griffin and De Sarkar used, which keeps us
comparable to them. An alternative is the Lambert et al. (2018) human-TF census,
which judges TFs by whether they have a DNA-binding domain rather than by
whether a motif exists. The script supports both
(`--census-type cisbp | lambert | intersection`). The choice changes *which TFs*
are on the list but never the sites within a TF that passes.

**Some chromatin proteins pass the gate.**
A protein passes CIS-BP if it has a DNA-binding domain with a known or inferred
motif. That lets through some chromatin enzymes, complex subunits and general
factors whose ChIP-seq peaks mostly reflect recruitment rather than their own
sequence-specific binding. Examples are KMT2A, KMT2B, KDM2B, CXXC1, KDM5B,
SETDB1, KAT7, SRCAP, TET2, TET3, BPTF, MBD2, MBD3, MECP2, GATAD2A, GATAD2B, HMG20A, HMGA1,
TBP, DRAP1, SP140 and NFYC. They are kept so the list stays comparable with
Griffin and De Sarkar. Expect their profiles to resemble a generic
CpG-island/promoter signal, and consider excluding them in sensitivity analyses.

**A note on gene names.**
The gate matches GTRD's TF names against CIS-BP's, so both sides need to agree
on the symbol. Most do. Occasionally a gene shows up under an old name: for
example, brachyury appears as `T` (its legacy symbol) rather than `TBXT`. CIS-BP
still lists it under `T` with a directly measured motif, so it passes
correctly. If you ever see a surprising name in the roster, check it with an
exact-field match rather than a loose text search:

```bash
awk -F'\t' 'toupper($7)=="T"' TF_Information.txt
```

---

## References

- Doebley A-L, Ko M, Liao H, et al. A framework for clinical cancer subtyping
  from nucleosome profiling of cell-free DNA. *Nat Commun* 13, 7475 (2022).
  doi:10.1038/s41467-022-35076-w
- De Sarkar N, Patton RD, Doebley A-L, et al. Nucleosome patterns in circulating
  tumor DNA reveal transcriptional regulation of advanced prostate cancer
  phenotypes. *Cancer Discov* 13(3):632–653 (2023). doi:10.1158/2159-8290.CD-22-0692
- Kolmykov S, Yevshin I, Kulyashov M, et al. GTRD: an integrated view of
  transcription regulation. *Nucleic Acids Res* 49(D1):D104–D111 (2021).
- Weirauch MT, Yang A, Albu M, et al. Determination and inference of eukaryotic
  transcription factor sequence specificity. *Cell* 158(6):1431–1443 (2014).
- Lambert SA, Jolma A, Campitelli LF, et al. The Human Transcription Factors.
  *Cell* 172(4):650–665 (2018).
