# Nucleosome profiling (Stage 2) — how to run

Infers transcription-factor accessibility from cfDNA fragment coverage (modified
Griffin). Three stages, chained by `run_profiling.sh` but each runnable on its own:

```
02_fcc_count.sh          genome-wide GC-corrected midpoint (FCC) track, per sample
      |   (FCC: one sparse *_counts.tsv.gz per chromosome)
      v
04_profile_tf_sites.py   slice ±5 kb windows around TFBS, bin to 15 bp, per chromosome
      |   (intermediate: <TF>_chr*.npz)
      v
04_merge_and_extract.py  outlier mask -> composite -> smooth -> normalize -> features
          (output: <TF>_profile.tsv.gz  +  all_TF_features.tsv)
```

You give it a **sample** and a **TF selection**; it produces composite profiles and a
features table.

---

## 1. Requirements

Activate the environment first — **the scripts do not activate conda themselves**
(this is deliberate, so they compose cleanly and work in the container):

```bash
conda activate tfp_env          # or: micromamba run -n tfp_env <cmd>
```

On PATH inside that env: `python` (with `pyyaml`), `gawk`, `samtools`, `parallel`.
The driver checks all four up front and fails in the first second if any is missing.

> **gawk, not mawk.** The ambiguous-base draw in `02` is seeded per fragment;
> mawk's `srand(seed)` is non-deterministic across processes, so under mawk the
> seeding silently fails. `02` errors if `gawk` is absent. Check with
> `readlink -f "$(command -v awk)"`.

---

## 2. One-time setup

Two files. The driver reads them **only via `config_helper.py`** (bash never parses
YAML itself).

### `config.yaml`

Two sections. Swap `deployment` when you move to the container; keep `method` frozen
so results stay comparable.

```yaml
deployment:
  sample_sheet:  /path/to/samples.tsv
  reference_fa:  /path/to/hg38.fa          # .fai MUST sit beside it
  mappable_bed:  /path/to/k100_minus_exclusion_lists.mappable_regions.hg38.bed
  chrom_sizes:   /path/to/hg38.standard.chrom.sizes
  tf_dir:        /path/to/TFBS_10000ms
  results_root:  /path/to/results
  tmpdir:        /tmp
  threads:       4

method:
  # Griffin-fixed — changing these breaks Griffin identity
  mapq:          30
  gc_bias_min:   0.05
  # deliberate divergence — must match Sara
  frag_len_min:  120
  frag_len_max:  200
  # profiling / feature params
  window_size:      5000
  step:             15
  smoothing_length: 165
  fft_window:       960
  fft_index:        10
  center_window:    30
  save_window:      1000
  autosomes: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12,
              chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22]
```

### `samples.tsv`

One row per sample. **`sample_id` must be unique** — a duplicate is a hard error.
Tab-separated, this exact header:

```
sample_id	bam	gc_bias	fcc_dir
SeCT-26_t2	/path/SeCT-26_t2.sortByCoord.bam	/path/SeCT-26_t2.GC_bias.txt	/path/SeCT-26_t2/fcc
```

`fcc_dir` is both where `02` writes FCC and where profiling reads it. If it already
holds `*_counts.tsv.gz`, counting is skipped (reused) unless you pass `--force-recount`.

---

## 3. Run the whole thing (usual case)

```bash
# all TFs
bash run_profiling.sh --sample SeCT-26_t2 --config config.yaml

# a few TFs
bash run_profiling.sh --sample SeCT-26_t2 --config config.yaml --tfs FOXA1,GRHL2,SPI1

# TFs from a file (one name per line)
bash run_profiling.sh --sample SeCT-26_t2 --config config.yaml --tfs /path/my_tfs.txt

# recount FCC even if it exists
bash run_profiling.sh --sample SeCT-26_t2 --config config.yaml --force-recount
```

Flags: `--sample` (required), `--config` (default `./config.yaml`), `--tfs`
(`ALL` | one name | comma-list | file; default `ALL`), `--force-recount`,
`--threads N` (override config for this run).

**TF names** are matched against the bed filename before the first `_`
(e.g. `FOXA1` matches `FOXA1.bed` and `FOXA1_meta.bed`). If two beds share that
prefix they **merge under one TF** — the driver warns when this happens, on both the
selective and the `ALL` path.

The driver refuses to proceed past any incomplete step: `02` must exit clean, every
autosome must have an FCC file, and every per-chromosome profiling job must succeed
before it merges. A partial run aborts loudly rather than reporting a false "DONE".

---

## 4. Run a stage standalone

Each stage stands alone — useful for debugging or partial reruns.

### `02` only (build the FCC track)

```bash
bash 02_fcc_count.sh \
  --bam      /path/SeCT-26_t2.sortByCoord.bam \
  --gc-bias  /path/SeCT-26_t2.GC_bias.txt \
  --out      /path/SeCT-26_t2/fcc \
  --ref      /path/hg38.fa \
  --mappable /path/mappable_regions.hg38.bed \
  --frag-min 120 --frag-max 200 --mapq 30 --gc-bias-min 0.05 --threads 4
```

Exits nonzero if any chromosome fails (checked via `parallel --joblog`, written to
`<out>/_joblog.tsv`).

### `04_profile` only (one chromosome)

```bash
python 04_profile_tf_sites.py \
  --config     config.yaml \
  --fcc_file   /path/SeCT-26_t2/fcc/chr1_counts.tsv.gz \
  --tf_dir     /path/TFBS_10000ms \
  --output_dir /path/results/SeCT-26_t2/profiling
```

`--config` supplies `window_size`, `step`, `chrom_sizes`. You can override any of
them with explicit flags (`--window_size`, `--step`, `--chrom_sizes`); explicit flags
win over config.

### `04_merge` only (composite + features)

```bash
python 04_merge_and_extract.py \
  --config      config.yaml \
  --output_dir  /path/results/SeCT-26_t2/profiling \
  --results_dir /path/results/SeCT-26_t2 \
  --n_jobs 4
```

Reads every `*_chr*.npz` under `<output_dir>/intermediate/`.

---

## 5. What you get

With `results_root: /path/results` and sample `SeCT-26_t2`:

```
/path/SeCT-26_t2/fcc/                         # FCC track (from the sample sheet)
    chr1_counts.tsv.gz ... chr22_counts.tsv.gz    cols: chrom  position  raw  fcc  (SPARSE)
    _joblog.tsv

/path/results/SeCT-26_t2/
    profiling/
        intermediate/   <TF>_chr*.npz          # per-chrom sliced windows (+ n_sites)
        composite/      <TF>_profile.tsv.gz    # cols: position mean_FCC smoothed_FCC normalized_FCC count TF
    all_TF_features.tsv                        # cols: TF mean_coverage central_coverage amplitude n_sites
    logs/                                      # benchmark logs + per-chrom stdout
```

- **Composite profiles** and the **features table** are the two deliverables.
- FCC is **sparse** (only covered positions written); this is bit-identical downstream
  to the old dense output — the reader zero-inits each window, so "absent" == "0".

---

## 6. Validate (when you change a script)

Stage 2 is a deliberate redesign of Griffin, not a port, so the correctness baseline
is this pipeline's **own** frozen output, not Griffin's. To confirm a change is
surgical:

```bash
bash golden_test.sh     # runs current vs edited 04_* on the same FCC; compares
```

Want: `STAGE-04 GOLDEN GATE: PASS` — profiles + all features identical, only the
intended column moved. See `COMPARISON.md` for the full Griffin-vs-this-pipeline
record and the list of deliberate divergences.

---

## 7. Troubleshooting

| Message | Cause / fix |
|---|---|
| `'gawk' is required` | mawk is your `awk`. Install gawk into the env / container. |
| `reference index missing: ...fai` | Build it beside the `.fa`: `samtools faidx hg38.fa`. A read-only mount can't build it at runtime — pre-build it. |
| `sample 'X' not found in ...` | `sample_id` not in `samples.tsv` (check spelling/tabs). |
| `sample 'X' appears N times` | Duplicate `sample_id` — must be unique. |
| `key 'Y' not found in config` | Missing key in `config.yaml`, or it's under the wrong section. |
| `FCC missing for autosome(s): ...` | Counting was incomplete, or `fcc_dir` points at a partial set. Rerun `02` / `--force-recount`. |
| `profiling incomplete (... failed ...)` | A chromosome job failed. See the joblog and `logs/prof_<chrom>.out` named in the error. |
| Log-file crash on import | Set `TFP_LOG_DIR` to a writable dir (the driver sets it automatically; only bites when running `04` by hand). |

---

## Notes / conventions

- FCC `position` is **1-based** (from BAM POS). BED `summit` is 0-based → a 1 bp
  offset that is sub-bin at 15 bp and averages out; documented, not corrected.
- Strand is **ignored** (GTRD meta-cluster beds are unstranded `.`).
- Autosomes are an **explicit chr1–chr22 allowlist**; `02` may still emit sex-chrom
  FCC if the mappable BED contains them, but the driver ignores anything outside the
  allowlist.
- Environment activation and log-dir belong to the **caller** (driver / container
  entrypoint), never the leaf scripts.
