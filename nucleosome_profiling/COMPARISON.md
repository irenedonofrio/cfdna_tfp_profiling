# TF Profiling — Comparison to upstream Griffin & change record

**Purpose.** This is the labelled record of how the nucleosome-profiling stage of this
pipeline relates to upstream Griffin, which divergences are *deliberate* (part of the
method) versus *accidental* (drift/bugs), and exactly what was changed in the
robustness/correctness pass. It is the Stage-2 analogue of the Stage-1 `NOTICE`.

**Upstream reference.** Griffin (Doebley et al. 2022), commit
`b624c7a2ea7b8758b6eebd424066c815bfeece61`, Clear BSD License, Fred Hutchinson Cancer
Research Center. Method influences: Ulz et al. 2019 (frequency decomposition / accessibility,
consumed downstream) and De Sarkar et al. 2023.

---

## 1. Architecture — a deliberate redesign, not a port

Griffin re-opens the BAM and fetches reads around **every site for every TF**
(`griffin_coverage.py` → per-site bigWig → `griffin_merge_sites.py`).

This pipeline instead computes a **genome-wide GC-corrected midpoint track once per sample**
(`02_fcc_count.sh`) and then **slices** windows out of that track
(`04_profile_tf_sites.py`), compositing and extracting features in
`04_merge_and_extract.py`. With hundreds of TF bed files this reads the BAM once instead of
hundreds of times.

**Consequence:** "bit-identical to Griffin" does **not** apply to this stage — it computes
different intermediate objects on purpose. The correctness baseline is this pipeline's *own
current output*, frozen as a golden file. Changes are therefore split into:
- **Track A — behaviour-preserving** (must reproduce the frozen output), and
- **Track B — behaviour-changing** (change output by design; justified individually).

## 2. Script mapping

| This pipeline | Griffin equivalent | Role |
|---|---|---|
| `02_fcc_count.sh` | `griffin_coverage.py` | GC-corrected midpoint counting (genome-wide track vs per-site bigWig) |
| `04_profile_tf_sites.py` | `griffin_merge_sites.py` (`fetch_bw_values`, `sum_bins`) | Slice ±5 kb windows around TFBS, bin to 15 bp |
| `04_merge_and_extract.py` | `griffin_merge_sites.py` (`make_outlier_mask`, `normalize_and_smooth`, `calculate_features`) | Outlier mask → composite → smooth → normalize → features |

## 3. Stage-by-stage verdict

| Step | Griffin | This pipeline | Verdict |
|---|---|---|---|
| Read filter | paired, mapq≥q, not dup, not qcfail, fwd + TLEN>0 | proper-pair, not dup, fwd (`-f 0x2 -F 0x410`), mapq≥30 | Kept (does not drop `qcfail 0x200` — negligible) |
| Fragment length | 100–200 | 120–200 | **Intended** (aligns with collaborator's filter) |
| GC from reference over full fragment span | yes | yes (span matches exactly) | Kept |
| GC base classes | A/T/W→0, C/G/S→1, ambiguous→seeded per-base RNG | see §5 | **Was drift → fixed** |
| GC-bias < 0.05 excluded | yes | yes (`>= 0.05`) | Kept |
| Midpoint | `floor((start+end)/2)`, 0-based | `int((2·POS+len)/2)`, 1-based | Formula kept; 1 bp coord offset in stage 04 (documented) |
| Windowing / bin sum | fetch ±5000, sum per 15 bp | slice ±5000, sum per 15 bp | Kept |
| Strand orientation | reverses `-` sites | strand ignored | **Intended** (GTRD sites unstranded, `.`) |
| Region/mappability exclusion | per-bin at merge | pre-filtered via mappable BED at counting | **Intended** |
| Outlier mask | z-score axis=None on *uncorrected*, `<10`, min_cutoff=2 fallback | identical | Kept (line for line) |
| Order: average → smooth → normalize | yes | yes | Kept (aggregate-then-normalize, same as Griffin) |
| Savgol window | `floor(165/15)=11`, odd, order 3 | same | Kept |
| Features (mean / central / FFT-index-10) | yes | yes (same windows & index) | Kept |
| `n_sites` | `len(current_sites)` | (was left-edge bin) | **Was bug → fixed** |

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

Because the S-fix and `n_sites` change output **by design**, freeze a golden output from the
*current* scripts first, then confirm only these move. See §8.

## 7. Deferred / documented-only

- **D2 distribution**: replace `int(rand()*2*N)` with a per-base coin-flip (sum of N draws) to
  match Griffin's binomial. Left for now; effect is tiny (ambiguous bases rare in mappable
  regions). Seeding already makes it deterministic **under gawk**.
- **GC counting in awk is inherently fragile** (S-class handling, RNG determinism, and
  awk-implementation dependence — mawk vs gawk). The robust long-term fix, consistent with
  Stage 1, is to compute fragment GC in the already-validated Python path and *share that one
  function* between the GC-bias builder and this corrector ("share the computation"), instead of
  maintaining a second awk implementation. Deferred, but this is the direction if the
  reproducibility guarantees need to tighten.
- **Remaining portability items** for the USZ container step: replace the Slurm array
  orchestrators with plain-bash drivers (done: `run_profiling.sh`); `02` still self-runs
  `conda activate tfp_env` internally (parameterise the interpreter); wire the drivers to
  `common/manifest.py`; ensure the container installs **gawk**.

## 8. Validation

Run `golden_test.sh` (see the script header for the config block). It:
1. runs the **current** `04_*` scripts and the **edited** `04_*` scripts on the *same* existing
   FCC input (two small chromosomes, a TF subset), and
2. calls `compare_outputs.py`, which asserts the composite **profiles are identical** and every
   feature **except `n_sites` is identical**, and reports the `n_sites` old→new change.

An optional awk micro-check demonstrates the S-fix and the seeding on synthetic sequences
without re-running the genome-wide counting.

`sparse_golden_test.sh` separately gates the sparse-FCC change (§9): it builds a sparse copy
from an existing dense FCC (no BAM/re-count needed), runs `04` on both, and asserts the
composite profiles and features are identical.

**Results (SeCT-26_t2, TFs FOXA1/GRHL2/SPI1):** both gates PASS. Stage-04 edits: profiles and
all coverage/amplitude features bit-identical, only `n_sites` moved. Sparse: profiles and all
features (incl. `n_sites`) bit-identical (max|Δ| = 0).

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

Note: gzip already compressed the long runs of zeros well, so the **disk** win understates the
real benefit. The zeros were expensive in **RAM**: `pd.read_csv` re-expands a dense file to its
full row count (~248M for chr1) regardless of gzip ratio, whereas the sparse file loads only the
covered rows. The memory/load-time reduction scales with the *row-count* drop (larger than 18×),
which is what removes the OOM risk when the USZ node runs many chromosomes in parallel. The
per-chromosome row counts printed by `sparse_golden_test.sh` are the figures to use for
container sizing.

**Caveats.** (1) Do not feed sparse FCC into legacy dense-assuming scripts (Maria's
`smooth_all_files.sh`, `filter_blacklist_bedtools.sh`). (2) The sparse gate proved the *format*
is safe by filtering existing dense FCC; when the sparse `02` is next used to recompute a sample,
run the optional `SPARSE_02_DIR` check in `sparse_golden_test.sh` to confirm the *script* emits
exactly those rows byte-for-byte.
