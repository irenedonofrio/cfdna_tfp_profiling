# tfp_profiling
 
cfDNA transcription-factor profiling pipeline, built on a optimized Griffin.
 
It has two **independent** stages that hand off through one small file
(`samples.GC.yaml`):
 
1. **GC correction** (`gc_correction/`) — optimized + modernized Griffin GC
   correction. Can be run entirely on its own.
2. **Nucleosome profiling** (`nucleosome_profiling/`) — the modified profiling
   stage. Reads the manifest produced by stage 1.
 
## Layout
- `gc_correction/`        – `griffin_GC_counts.py` (optimized), `griffin_GC_bias.py` (modernized, `--no_plots`)
- `nucleosome_profiling/` – profiling scripts
- `common/`               – `manifest.py`: the single source of truth for the stage-1→stage-2 handoff
- `workflows/`            – `run_gc_correction.sh`: run stage 1 from the terminal (no scheduler)
- `tests/`                – correctness + regression + demo checks
- `env/`                  – conda environment (to be added)
- `container/`            – Apptainer/Singularity definition (to be added)
 
## Quick checks
    conda activate griffin
    python tests/test_gc_identity.py          # GC-count correctness, no genome needed
 
## Run GC correction (needs hg38.fa + the Griffin Ref/ files)
    # edit the CONFIG block in workflows/run_gc_correction.sh first
    ./workflows/run_gc_correction.sh samples.tsv
 
See NOTICE for upstream attribution and a summary of changes.
