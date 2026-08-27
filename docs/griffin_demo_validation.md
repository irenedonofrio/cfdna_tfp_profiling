# Griffin demo validation — CTCF

**Date:** 2026-08-27
**Validates:** commit `7cc4fa7` — "Fix feature-window fencepost: match Griffin's half-open arange"
**Run on:** LeoMed, conda env `ae_env`

Provenance of every number below: **real data**, unless marked otherwise.

---

## Purpose

Confirm that `cfdna_tfp_profiling` reproduces Griffin on data where Griffin's
own expected output is published, after the feature-window fencepost fix.

This is an implementation check. It does **not** restate `COMPARISON.md` §11,
which is CDX2 on BH01 at a 120–200 fragment window and still requires
re-derivation.

---

## Setup

Griffin's published demo: `Healthy_GSM1833219_downsampled.sorted.mini.cram`,
0.26x, hg38, converted to BAM against `Griffin/Ref/hg38.fa`.

Sites: Griffin's own `CTCF.hg38.1000.txt` (1000 sites), converted to BED6 by
taking the `position` column as the summit and writing `start = position - 1`
(1-based → 0-based).

Parameters matched to Griffin for this run:

| | value |
|---|---|
| mapping quality | 20 (not the usual 30) |
| fragment length window | 100–200 (not the usual 120–200) |
| chromosomes | autosomes chr1–chr22 |

Stage-1 GC correction left at its validated defaults (`map_q: 20`,
`size_range: [15, 500]`) — these are byte-identity validated against Griffin
and were not touched.

Reference: `griffin_nucleosome_profiling_demo_files/expected_results/Healthy_demo.GC_corrected.coverage.tsv`,
which ships both the composite profile (132 position columns, −990…975) and
the feature values.

---

## Result — profile

| | |
|---|---|
| Pearson r | 0.990419 |
| max abs difference | 0.094506 |
| mean difference | +0.003235 |
| n_sites | 1000 both sides |

The residual is oscillatory and largest where the profile is steepest,
consistent with a small difference in the fragment population rather than a
systematic offset.

---

## Result — features

| feature | Griffin | this pipeline | difference |
|---|---|---|---|
| `mean_coverage` | 0.850336 | 0.853571 | +0.38% |
| `central_coverage` | 0.495717 | 0.493036 | −0.54% |
| `amplitude` | 13.958096 | 14.485981 | +3.78% |

---

## Result — feature window bin counts

| window | Griffin | this pipeline |
|---|---|---|
| `save_window` ±1000 | 132 | 132 |
| `center_window` ±30 | 4 | 4 |
| `fft_window` ±960 | 128 | 128 |

**This is the point of the exercise.** Before commit `7cc4fa7` these were
133 / 5 / 129 on this side, because the feature windows used closed intervals
(`>= -w & <= w`) where Griffin snaps each bound onto the step grid and then
builds a half-open `np.arange` (`griffin_merge_sites.py:226-228, 232-234`).
The bin counts now agree against Griffin's actual published output, not
against a reconstruction of it.

Features computed from the stored `all_TF_features.tsv` and independently
recomputed from the composite agreed exactly, so the comparison is not relying
on a reimplementation of `extract_features`.

---

## Residual attribution

Not measured here — **source code**:

- Mappability handling differs: this pipeline pre-filters via the mappable BED
  at counting; Griffin excludes per-bin at merge.
- Autosomes only.
- Three fixed bugs in Griffin's logic, which make this pipeline's output
  intentionally *more* correct than Griffin's rather than identical to it:
  `S` counted as GC, the ambiguous-base draw seeded per fragment, and
  `n_sites` read from `len(current_sites)`.

`amplitude` at +3.78% is the largest residual and the least explained.
Griffin's `outlier_cutoff` metadata column has not been checked against this
run's fallback (the outlier mask reported "Low coverage, resetting cutoff to 2",
expected at 0.26x). Worth checking before attributing the gap.

---

## Reproducing

```bash
conda activate ae_env
cd <repo>
bash run_all.sh \
  --gc-config   $PWD/gc_correction/config_gc_demo.yaml \
  --prof-config $PWD/nucleosome_profiling/config_demo.yaml \
  --out-root    <demo dir>/run \
  --tfs CTCF --threads 8
```

Comparison notebook: `compare_griffin_demo.ipynb`.

Note: `02_fcc_count.sh` requires `gawk`, and Stage 2 requires GNU `parallel` on
PATH. Both are present in `ae_env` on LeoMed.
