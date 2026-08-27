# TF Profiling — Comparison to upstream Griffin & change record

This file is to  record the divergences and exactly what was changed. It is the Stage-2 analogue of the Stage-1 `NOTICE`. It also records pipeline-wide infrastructure changes (environment consolidation, path handling) that affect both stages — see §10.

**Upstream reference.** Griffin (Doebley et al. 2022), commit
`b624c7a2ea7b8758b6eebd424066c815bfeece61`, Clear BSD License, Fred Hutchinson Cancer
Research Center. 

---

## 1. Architecture — a deliberate redesign, not a port

Griffin re-opens the BAM and fetches reads around **every site for every TF**
(`griffin_coverage.py` → per-site bigWig → `griffin_merge_sites.py`).

This pipeline instead computes a **genome-wide GC-corrected midpoint track once per sample**
(`02_fcc_count.sh`) and then **slices** windows out of that track
(`04_profile_tf_sites.py`), compositing and extracting features in
`04_merge_and_extract.py`. With hundreds of TF bed files this reads the BAM once instead of
hundreds of times.

**Consequence:** "bit-identical to Griffin" does **not** apply to this stage: Because this stage works differently from Griffin, its output is close to Griffin's but not an exact
(On CDX2 with the same sites, the two profiles match to r = 0.997 and the features to within a few percent; see §11.)

## 2. Script mapping

| This pipeline | Griffin equivalent | Role |
|---|---|---|
| `02_fcc_count.sh` | `griffin_coverage.py` | GC-corrected midpoint counting (genome-wide track vs per-site bigWig) |
| `04_profile_tf_sites.py` | `griffin_merge_sites.py` (`fetch_bw_values`, `sum_bins`) | Slice ±5 kb windows around TFBS, bin to 15 bp |
| `04_merge_and_extract.py` | `griffin_merge_sites.py` (`make_outlier_mask`, `normalize_and_smooth`, `calculate_features`) | Outlier mask → composite → smooth → normalize → features |

## 3. Stage-by-stage verdict

| Step | Griffin | This pipeline | Verdict                                                  |
|---|---|---|----------------------------------------------------------|
| Read filter | paired, mapq≥q, not dup, not qcfail, fwd + TLEN>0 | proper-pair, not dup, fwd (`-f 0x2 -F 0x410`), mapq≥30 | Kept (does not drop `qcfail 0x200` — negligible)         |
| Fragment length | 100–200 | 120–200 | **Intended** (aligns with Maria Lambrinos's filter)      |
| GC from reference over full fragment span | yes | yes (span matches exactly) | Kept                                                     |
| GC base classes | A/T/W→0, C/G/S→1, ambiguous→seeded per-base RNG | see §5 | **Was drift → fixed**                                    |
| GC-bias < 0.05 excluded | yes | yes (`>= 0.05`) | Kept                                                     |
| Midpoint | `floor((start+end)/2)`, 0-based | `int((2·POS+len)/2)`, 1-based | Formula kept; 1 bp coord offset in stage 04 (documented) |
| Windowing / bin sum | fetch ±5000, sum per 15 bp | slice ±5000, sum per 15 bp | Kept                                                     |
| Strand orientation | reverses `-` sites | strand ignored | **Intended** (GTRD sites unstranded, `.`)                |
| Region/mappability exclusion | per-bin at merge | pre-filtered via mappable BED at counting | **Intended**                                             |
| Outlier mask | z-score axis=None on *uncorrected*, `<10`, min_cutoff=2 fallback | identical | Kept (line for line)                                     |
| Order: average → smooth → normalize | yes | yes | Kept (aggregate-then-normalize, same as Griffin)         |
| Savgol window | `floor(165/15)=11`, odd, order 3 | same | Kept                                                     |
| Features (mean / central / FFT-index-10) | window = `np.arange(ceil(lo/step)*step, floor(hi/step)*step, step)` — **half-open**, right endpoint excluded (`griffin_merge_sites.py:226-228, 232-234`) | same windows & index **only after the fencepost fix**; before it the masks were closed (`>= -w & <= w`) | **Was drift → fixed**                                    |
| `n_sites` | `len(current_sites)` | (was left-edge bin) | **Was bug → fixed**                                      |

**Correction (feature windows).** The "same windows & index" verdict previously recorded on the
features row was **false**. Griffin snaps each window bound onto the step grid and then builds the
columns with a half-open `np.arange`, so the right endpoint is not a bin. This pipeline used closed
intervals (`>= -w & <= w`), which admitted one extra bin in every feature window at 15 bp bins:

| Window | This pipeline (pre-fix) | Griffin | Consequence |
|---|---|---|---|
| `save_window` ±1000 | 133 bins (−990…990) | 132 bins (−990…975) | `mean_coverage` averaged over one extra bin |
| `center_window` ±30 | 5 bins (−30…30) | 4 bins (−30…15) | `central_coverage` averaged over one extra bin |
| `fft_window` ±960 | 129 samples (−960…960) | 128 samples (−960…945) | `fft_index=10` read a **193.5 bp** period, not Griffin's **192.0 bp** |

Note that ±1000 is not a multiple of the 15 bp step, so a bare `<` does **not** reproduce Griffin's
save window — the bound must be floored onto the grid first. Fixed by matching Griffin's
construction exactly in `04_merge_and_extract.py:extract_features`, which now also takes `step` as a
required argument and rejects a `step` that disagrees with the composite's bin spacing.

## 4. Intended divergences (keep — these define the method)

- Genome-wide precompute + slice architecture.
- Fragment length 120–200.
- Mapping quality 30.
- Mappable-BED pre-filter replacing Griffin's post-hoc exclusion steps.
- Explicit autosome list (chr1–chr22).
- Corrected-track features only (Griffin also emits uncorrected/map-corrected).
- Strand ignored (GTRD meta-cluster beds are unstranded, so a strand flip would be a no-op).

## 5. Drift / bugs found

| # | Item | Effect | Status |
|---|---|---|---|
| D1 | `S` not counted as GC (Griffin counts S=1) | GC undercount on rare S bases → wrong bias key | **Fixed** |
| D2 | Ambiguous-base draw `int(rand()*2*N)` unseeded | non-reproducible run-to-run; can overcount | **Seeded (Fixed, requires gawk); distribution deferred** |
| D3 | `n_sites = count.iloc[0]` (left-edge bin) | undercounts sites → distorts Ulz noise correction | **Fixed** |
| D4 | 1 bp coord offset (1-based FCC vs 0-based BED summit) | sub-bin, no feature effect | Documented only |
| D5 | dead `end` var + redundant genome-wide zero-fill in `02` | inflated I/O & memory; trailing fill never runs | **Fixed (sparse output) & validated — see §9** |
| D6 | hardcoded benchmark-log paths in both `04_*.py` | crashes off-cluster / in container | **Fixed** (now `TFP_LOG_DIR`, default `./logs`) |
| D7 | `GC_corrected` not rounded to Griffin's 5 dp | only matters for byte-comparison to Griffin | Deferred |

## 6. Changes applied in this pass

Every edit is tagged `# CHANGED:` in the source; comment-only notes are tagged `# NOTE:`.

| Change | File | Track | Expected output effect |
|---|---|---|---|
| **S-fix**: GC class `[GCgc]` → `[GCScgs]` | `02_fcc_count.sh` | B | Tiny: only fragments overlapping `S` bases; shifts their bias key by one GC. |
| **Seed RNG**: `srand(fragment_start)` before ambiguous draw; **force gawk + guard** | `02_fcc_count.sh` | B→A | Makes ambiguous-base draws reproducible run-to-run. **Requires gawk** — mawk's `srand(seed)` is non-deterministic, so the script now calls `gawk` explicitly and errors if it is absent. |
| **`n_sites` = true count**: capture `len(tf_sites)` per chrom, store in npz, sum across chroms, use in features | `04_profile_tf_sites.py`, `04_merge_and_extract.py` | B | `n_sites` column changes to the true per-TF site total. No effect on profiles or coverage/amplitude features. |
| **Sparse FCC**: drop the genome-wide zero-fill; write only covered positions (+ remove dead trailing fill) | `02_fcc_count.sh` | A | **None** — non-zero rows are byte-identical; the reader treats absent as 0. Validated bit-identical downstream (§9). Big I/O/RAM reduction. |
| **Log path configurable**: `TFP_LOG_DIR` (default `./logs`), created safely | `04_profile_tf_sites.py`, `04_merge_and_extract.py` | A | None — only moves the benchmark log; outputs unchanged. Removes the off-cluster/container import crash. |
| Strand note | `04_profile_tf_sites.py` | — | none (documentation) |
| Coord-offset note | `04_profile_tf_sites.py` | — | none (documentation) |



## 7. Validation

Run `golden_test.sh`. It:
1. runs the **current** `04_*` scripts and the **edited** `04_*` scripts on the *same* existing
   FCC input (two small chromosomes, a TF subset), and
2. calls `compare_outputs.py`, which asserts the composite **profiles are identical** and every
   feature **except `n_sites` is identical**, and reports the `n_sites` old→new change.


`sparse_golden_test.sh` separately gates the sparse-FCC change (§9): it builds a sparse copy
from an existing dense FCC (no BAM/re-count needed), runs `04` on both, and asserts the
composite profiles and features are identical.

**Results (SeCT-26_t2, TFs FOXA1/GRHL2/SPI1).**

*Sparse gate:* PASS. Profiles and all features (incl. `n_sites`) bit-identical (max|Δ| = 0).
Re-verified against the stored `tests/sparse_test_work/` artifacts.

*Stage-04 gate:* **the result recorded here — "PASS, only `n_sites` moved" — was wrong.** The
stored artifacts of that run (`tests/golden_test_work/{orig,new}/`, both trees written 2026-07-01)
show that `central_coverage` **and** `amplitude` moved for all four TFs, with the composite
profiles bit-identical. CDX2: amplitude 0.126326 → 0.218102 (+73%), central_coverage 0.798939 →
0.799087. `compare_outputs.py` gates both of those columns at `--feature-tol 0`, so on those
artifacts it prints **FAIL (something else moved)** — not PASS. Stating it plainly: this section
recorded a pass for a gate that fails.

The cause is the feature-window fencepost corrected in §3. The stage-04 edit changed `center_mask`
and `fft_mask` from Griffin's half-open convention to closed intervals. `mean_coverage` did **not**
move in that run because its window was **already** closed before the edit — that fencepost
predates the stage-04 change and is present in **both** stored runs, which is why it left no trace
in this gate.

## 9. Dense vs sparse FCC

`02` previously wrote one row per base per chromosome (~248M rows for chr1), the vast majority
`chrom  pos  0  0`, since cfDNA midpoint coverage is sparse. `04_profile_tf_sites.py`
zero-initialises each window and overwrites only positions present in the file, so **absent ==
present-with-0** to the reader. The sparse `02` therefore emits only covered positions; the
non-zero rows are byte-identical to the old output's non-zero rows (same tab layout, same
`%.6g` values), so downstream results are unchanged — validated in §8.

Measured on `SeCT-26_t2` (whole sample):

| | Dense | Sparse | Ratio |
|---|---|---|---|
| On-disk (gzipped) | 7.2 GB | 412 MB | ~18× smaller (sparse ≈ 5.6%) |
| chr1 | 630 MB | 33 MB | |
| chr2 | 621 MB | 40 MB | |
| chr22 | 127 MB | 5.2 MB | |


**Caveats.** Do not feed sparse FCC into legacy dense-assuming scripts (Maria's
`smooth_all_files.sh`, `filter_blacklist_bedtools.sh`)

## 10. Pipeline infrastructure changes (both stages)


### 10.1 Environment consolidation — two envs → one `tfp`

**What.** The two runtime environments (`griffin` for Stage 1, `tfp_env` for Stage 2) were
replaced by a single minimal `environment.yml`, scoped to only what the scripts actually import
or call: `numpy`, `pandas`, `scipy`, `pysam`, `pyyaml`, `matplotlib-base`; `gawk`, `samtools`,
`parallel`. Everything else from the old exports (R/tidyverse, jupyter, snakemake, picard,
seaborn, scikit-learn, bedtools, whittaker-eilers, …) was dropped as unused. `gawk` was **added
explicitly** — it was absent from both prior exports and `02_fcc_count.sh` requires it.

**Version anchor.** The single env is pinned on the validated **Stage-2** stack
(pandas 2.3.3, numpy 2.3.4, scipy 1.16.2, pysam 0.23.3, samtools 1.22.1). Rationale:
1. Stage 2 is the more validation-sensitive half (float composite/FFT outputs), so it is left
   untouched → **no Stage-2 re-validation required**.
2. Stage 1's pandas surface is trivial (integer-tally I/O); its numerically sensitive parts
   (ambiguous-base RNG, cumulative-GC arrays) are numpy-driven and numpy stays on 2.x, so the
   RNG stream is unchanged.
3. `pysam 0.23.3` + `samtools 1.22.1` share one htslib in this stack; anchoring on Stage-1's
   `pysam 0.24.0` instead would want htslib 1.23 and conflict with `samtools 1.22.1`'s 1.22.

**Track.** A (output-preserving), pending the Stage-1 re-validation below.

**Validation.** Stage 1 was re-run on the new `tfp` env and diffed against its frozen golden
(`SeCT-26_t1`): `GC_counts` vs `original_SeCT-26_t1` (unoptimized Griffin) and `GC_bias` vs
`bench_bias_new`, both `diff -q` **byte-identical**. This validates the env move
(pandas 3.0→2.3, pysam 0.24→0.23, numpy 2.4→2.3) in one run. Stage 2's stack is unchanged, so no
Stage-2 numeric re-validation was needed; the end-to-end `run_all.sh` run
(`SeCT-26_t1`, FOXA1/GRHL2/SPI1) completed with all 22 chromosomes profiled and features written.


### 10.2 Relocatable paths — bundled refs + repo-relative outputs

**What.** Both configs now use a `refs_root` + relative ref names, and outputs are repo-relative,
so the whole project folder can be moved or mounted anywhere without editing configs. Relative
paths resolve against the **config file's own directory** (not the shell's CWD); an optional
`TFP_REFS` environment variable overrides `refs_root` for refs mounted outside the repo.

**Where the resolution lives.**
- Stage 1 driver `run_gc_correction.sh` resolves its own config (inline).
- Stage 2 `config_helper.py` resolves ref-typed keys onto `refs_root` and dir-typed keys
  (`sample_sheet`, `results_root`, `fcc_root`) onto the config dir.
- `run_all.sh` anchors its manifest-path lookup to the Stage-1 config dir.

**Fix B (`04_profile_tf_sites.py`).** This script reads `chrom_sizes` **directly** from the YAML
(bypassing `config_helper.py`), so a relative `chrom_sizes` failed with `FileNotFoundError` when
the script ran under a relative config. Fixed by adding `resolve_ref_path()`, which replicates
`config_helper.py`'s rule exactly (relative ref value joined onto `refs_root`; `refs_root` from
`TFP_REFS` or the config, anchored to the config dir; absolute values pass through). The
script's own `chrom_sizes` is resolved through it before the CLI > config > fallback precedence.

**Note — duplicated rule.** The anchoring logic now exists in two places (`config_helper.py` and
`04_profile_tf_sites.py::resolve_ref_path`). They are behaviourally identical by construction; if
the path-resolution rule ever changes, change both. This duplication was accepted deliberately so
that `04_profile_tf_sites.py` is self-sufficient with a relative config whether launched by the
driver or run standalone.

**Track.** A (output-preserving) — path handling only; no change to any computed value. Absolute
configs behave exactly as before (every resolver returns absolute inputs unchanged), so the
change is a strict superset. Confirmed by the end-to-end `run_all.sh` run in §10.1.


## 11. External validation vs upstream Griffin — CDX2, identical sites

Direct comparison against **upstream Griffin's actual output** on the **same
BAM and the same site list**, so any difference is attributable to the method, not the inputs.

**Setup.** Sample `SeCT-26_t1`, TF `CDX2`. To remove the confound of different bed files, CDX2 was
re-profiled on **Griffin's own site list**
(`/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/CDX2_binding_sites_without_blacklisted.tsv`, 286558 sites, `position` column used as the summit
verbatim), not the GTRD `TFBS_10000ms` bed. Griffin reference:
`/cluster/work/medinfmk/cfDNA-SeCT/results/BAM/gc_uncorrected/griffin/snakemakes/griffin_nucleosome_profiling/results/sample_name_12.GC_corrected.coverage.tsv`, `GC_corrected` +
smoothed row (the analogue of this pipeline's `normalized_FCC`). Compared on the overlapping
±990 bp window (132 bins, 15 bp grid).

**Result.**

| Quantity | Griffin | This pipeline | Δ |
|---|---|---|---|
| `n_sites` / `number_of_sites` | 286558 | 286558 | **0 (exact)** |
| Profile Pearson r | — | — | **0.9970** |
| mean ratio (yours/griffin) | — | — | 1.0031 |
| residual max\|Δ\| / rms | — | — | 0.0070 / 0.0035 |
| mean_coverage | 0.98296 | 0.98605 | +0.31% |
| central_coverage | 0.96807 | 0.97414 | +0.63% |
| amplitude | 0.05774 | 0.05407 | −6.3% |

> **All three scalar-feature rows above were computed with the pre-fix (closed-interval) feature
> windows and are superseded.** They are retained only as a record of what was measured. See the
> amplitude bullet below.

**Reading.**
- **`n_sites` exact.** On the identical bed, this pipeline's site loader reproduces Griffin's
  `len(current_sites)` to the site
- **r = 0.9970.** The composites are the same curve: both resolve the CDX2 central nucleosome
  depletion and flanking structure identically in shape.
- **The residual is smooth, sign-consistent, and centre-weighted.** Candidate contributing terms:
  the tighter fragment window (120–200 vs Griffin's 100–200), mapq 30, and the different
  normalization window (~0.3% global level, `mean ratio` 1.0031).
- **amplitude −6.3% — attribution NOT established. Requires re-derivation.** This row was computed
  with the closed-interval feature windows (§3): a 129-sample FFT window where Griffin uses 128, so
  `fft_index=10` read a 193.5 bp period against Griffin's 192.0 bp. That fencepost is not a small
  term. Measured on the stage-04 golden TFs, switching that one window between the two conventions
  moved amplitude by −0.5% (SPI1), +6.7% (GRHL2), +13.9% (FOXA1) and +73% (CDX2) — a different
  sample and site set, so not a transfer of magnitude to the CDX2/Griffin comparison, but enough to
  show the term can exceed 6.3% and can carry either sign. The −6.3% therefore cannot be read as
  the consequence of the fragment window and mapq. **No replacement number is offered here: this
  comparison must be re-run against Griffin's output after the fix.**

**Verdict.** On identical sites, the redesigned genome-wide-track-and-slice architecture
reproduces Griffin's CDX2 profile **shape** near-perfectly (r = 0.9970). That part stands: `r` is
computed on the overlapping ±990 bp window and does not depend on the feature-window fencepost.

The **scalar-feature** half of this verdict is withdrawn pending re-derivation. The claim that the
features agree "within a few percent, with every residual attributable to a documented intended
divergence (§4)" is not supported: the feature-window fencepost was neither documented nor
intended, and it is a larger term than the divergences it was attributed to.

**Caveat for the features table.** This CDX2 row was computed on Griffin's bed (286558 sites),
whereas the pipeline's default `TFBS_10000ms` CDX2 bed differs. If CDX2 is later re-profiled under
the default `tf_dir`, its row is replaced (§ upsert behaviour) and `n_sites` will no longer be
286558. Keep this Griffin-bed CDX2 result isolated (separate `--results-root`) if it must persist.
