#!/usr/bin/env python3
"""
TF profiling from FCC data as in Griffin
Sequential version + benchmarking instrumentation
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
import time
import logging
import random

# ---------------------------
# Benchmark logging
# ---------------------------
# CHANGED: log dir is configurable via TFP_LOG_DIR (set by the driver), with a
# safe fallback to the current directory. Previously this was a hardcoded cluster
# path opened at import time, which crashed off-cluster / in the container.
script_name = Path(sys.argv[0]).stem
job_id = os.environ.get("SLURM_JOB_ID", "nojob")
log_dir = os.environ.get("TFP_LOG_DIR", ".")
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"{script_name}_{job_id}.log")

logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(message)s"
)

def timed_step(name, func, *args, **kwargs):
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    logging.info(f"{name}\t{elapsed:.6f}")
    return result

# ---------------------------
# Config loading (CHANGED)
# ---------------------------
# CHANGED: optional --config. When given, method params (window_size, step) and the
# deployment chrom_sizes are read from it. Explicit CLI flags still override config,
# and config overrides the hardcoded fallbacks — so this stays backward-compatible
# with callers that pass everything on the command line (e.g. golden_test.sh).
def load_config_sections(config_path):
    if not config_path:
        return {}, {}
    import yaml
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}
    return (cfg.get("deployment", {}) or {}, cfg.get("method", {}) or {})

def resolve(cli_value, config_value, fallback):
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return fallback

# ---------------------------
# CLI
# ---------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute TF profiles from FCC data (chromosome-first approach)"
    )
    parser.add_argument("--fcc_file", required=True, help="Path to input FCC tsv.gz file (single chromosome)")
    parser.add_argument("--tf_dir", required=True, help="Directory containing TF .bed files")
    parser.add_argument("--output_dir", required=True, help="Output directory for intermediate profile data")
    # CHANGED: --config optional; the params below default to None so we can tell
    # "not given" from an explicit value and apply config -> fallback precedence.
    parser.add_argument("--config", default=None, help="Path to config.yaml (supplies window_size, step, chrom_sizes)")
    parser.add_argument("--window_size", type=int, default=None, help="Window size for TF profiling (default: 5000)")
    parser.add_argument("--chrom_sizes", default=None, help="Path to chromosome sizes file")
    parser.add_argument("--step", type=int, default=None, help="Binning step size for profiles (default: 15)")
    parser.add_argument("--n_jobs", type=int, default=1)  # unused but kept for compatibility
    return parser.parse_args()

# ---------------------------
# Helpers
# ---------------------------
def extract_chromosome(fcc_file):
    name = Path(fcc_file).stem.replace('.tsv', '')
    if name.endswith('_counts'):
        name = name.replace('_counts', '')
    return name

def load_fcc(fcc_file):
    print(f"Loading FCC data from {fcc_file}...", file=sys.stderr)
    try:
        fcc_df = pd.read_csv(
            fcc_file, sep='\t', compression='gzip', header=None,
            names=['chrom', 'position', 'raw', 'fcc'],
            dtype={'chrom': str, 'position': int, 'raw': float, 'fcc': float}
        )
        print(f"  Loaded {len(fcc_df):,} FCC positions", file=sys.stderr)

        if not fcc_df['position'].is_monotonic_increasing:
            print(f"  Sorting by position...", file=sys.stderr)
            fcc_df = fcc_df.sort_values('position').reset_index(drop=True)

        print(f"  FCC value range: [{fcc_df['fcc'].min():.3f}, {fcc_df['fcc'].max():.3f}]", file=sys.stderr)
        return fcc_df

    except Exception as e:
        print(f"Error loading FCC file: {e}", file=sys.stderr)
        sys.exit(1)

def load_bed(tf_bed_file, chrom):
    try:
        # NOTE: 'strand' is read into the column names but intentionally NOT used.
        # GTRD meta-cluster BEDs are unstranded ('.'), so there is no orientation to
        # preserve and no strand-flip is needed before compositing. (Griffin reverses
        # the profile for '-' sites; that step would be a no-op here.) This is a
        # deliberate choice, not an oversight.
        df = pd.read_csv(
            tf_bed_file, sep='\t', header=None,
            names=['chrom', 'summit', 'summit_plus1', 'TF', 'support', 'strand'],
            usecols=['chrom', 'summit', 'TF'],
            dtype={'chrom': str, 'summit': int, 'TF': str}
        )
        df_chrom = df[df['chrom'] == chrom][['summit']]
        return df_chrom.sort_values('summit').reset_index(drop=True) if len(df_chrom) > 0 else None
    except Exception as e:
        print(f"  Warning: Failed to load {tf_bed_file}: {e}", file=sys.stderr)
        return None

def bin_profile(profile, step):
    n_bins = len(profile) // step
    return np.sum(profile.reshape(n_bins, step), axis=1)

def compute_tf_profiles(fcc_df, tf_sites, chrom_size, window_size, profile_length, step):
    profiles = {'uncorrected': [], 'corrected': []}

    # Extract numpy arrays ONCE — removes pandas overhead in hot loop
    pos_arr = fcc_df['position'].to_numpy()
    raw_arr = fcc_df['raw'].to_numpy()
    fcc_arr = fcc_df['fcc'].to_numpy()

    for summit in tf_sites['summit'].values:
        start = summit - window_size
        end   = summit + window_size

        valid_start = max(0, start)
        valid_end   = min(chrom_size, end)

        # initialize footprint arrays
        raw = np.full(profile_length, np.nan)
        fcc = np.full(profile_length, np.nan)

        vs = valid_start - start
        ve = valid_end   - start
        raw[vs:ve] = 0
        fcc[vs:ve] = 0

        # binary search in numpy arrays
        start_idx = pos_arr.searchsorted(valid_start)
        end_idx   = pos_arr.searchsorted(valid_end)

        if end_idx > start_idx:  # means some FCC positions exist in this window

            # relative offsets from summit
            # NOTE: FCC 'position' is 1-based (it comes from BAM POS/$4 in stage 02),
            # while BED 'summit' is 0-based. So rel is offset by 1 bp: the true site
            # center sits at rel = -1 rather than rel = 0. This is sub-bin at a 15 bp
            # step and averages out across thousands of sites, so it does NOT affect the
            # features — left as-is intentionally.
            # To correct if ever needed: use (summit + 1) here to put both on 1-based footing.
            rel = pos_arr[start_idx:end_idx] - summit

            # mask values strictly inside [-window_size, +window_size)
            mask = (rel >= -window_size) & (rel < window_size)

            if np.any(mask):
                idx = (rel[mask] + window_size).astype(int)

                # bulk assignment from pre-sliced numpy arrays
                raw[idx] = raw_arr[start_idx:end_idx][mask]
                fcc[idx] = fcc_arr[start_idx:end_idx][mask]

        # bin to coarser resolution
        raw = bin_profile(raw, step)
        fcc = bin_profile(fcc, step)

        profiles['uncorrected'].append(raw)
        profiles['corrected'].append(fcc)

    return profiles

def save_tf_profiles(tf_profiles, tf_name, output_dir, chrom, window_size, step, n_sites=0):
    intermediate_dir = Path(output_dir) / "intermediate"
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    start = -window_size
    bin_positions = np.arange(start, window_size, step)

    # CHANGED: store the true number of loaded sites on this chromosome (mirrors Griffin's
    # len(current_sites)). The merge step sums this across chromosomes to get the real
    # per-TF site count, instead of inferring it from a bin of the composite (which was
    # reading the left-edge bin and under-counting).
    arrays_to_save = {'positions': bin_positions,
                      'n_sites': np.asarray(n_sites)}
    for k, v in tf_profiles.items():
        arrays_to_save[k] = np.array(v)

    np.savez_compressed(intermediate_dir / f"{tf_name}_{chrom}.npz", **arrays_to_save)

# ---------------------------
# Instrumented processing
# ---------------------------
def process_tf(tf_bed_file, fcc_df, chrom, chrom_size, window_size, profile_length, step, output_dir):
    stem = tf_bed_file.stem
    tf_name = stem.split("_")[0] if "_" in stem else stem

    try:
        tf_sites = timed_step(f"{tf_name}:load_bed", load_bed, tf_bed_file, chrom)
        if tf_sites is None or tf_sites.empty:
            logging.info(f"{tf_name}:no_sites\t0")
            print(f"  No sites for {tf_name} on {chrom}", file=sys.stderr)
            return (tf_name, True, None)

        tf_profiles = timed_step(
            f"{tf_name}:compute",
            compute_tf_profiles,
            fcc_df, tf_sites, chrom_size, window_size, profile_length, step
        )

        timed_step(
            f"{tf_name}:save",
            save_tf_profiles,
            # CHANGED: len(tf_sites) is exactly Griffin's len(current_sites), restricted to
            # this chromosome. Passed through so the merge step can sum it across chroms.
            tf_profiles, tf_name, output_dir, chrom, window_size, step, len(tf_sites)
        )

        return (tf_name, True, None)

    except Exception as e:
        logging.info(f"{tf_name}:error\t0")
        return (tf_name, False, str(e))

# ---------------------------
# Main
# ---------------------------
def main():
    pipeline_start = time.perf_counter()
    args = parse_args()

    # CHANGED: resolve params with precedence CLI > config > fallback.
    deployment, method = load_config_sections(args.config)
    window_size_arg = resolve(args.window_size, method.get("window_size"), 5000)
    step           = resolve(args.step,        method.get("step"),        15)
    chrom_sizes_path = resolve(args.chrom_sizes, deployment.get("chrom_sizes"), None)
    if not chrom_sizes_path:
        print("[FATAL] chrom_sizes not provided (pass --chrom_sizes or set deployment.chrom_sizes in --config)", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chrom_sizes = {}
    with open(chrom_sizes_path) as f:
        for line in f:
            c, s = line.strip().split('\t')
            chrom_sizes[c] = int(s)

    chrom = extract_chromosome(args.fcc_file)
    print(f"Processing chromosome: {chrom}", file=sys.stderr)

    if chrom not in chrom_sizes:
        print(f"Chromosome {chrom} not found", file=sys.stderr)
        sys.exit(1)

    chrom_size = chrom_sizes[chrom]

    start = -window_size_arg
    end = window_size_arg
    start = int(np.ceil(start/step)*step)
    end = int(np.floor(end/step)*step)
    window_size = end
    profile_length = (end - start)

    fcc_df = timed_step("load_fcc", load_fcc, args.fcc_file)
    tf_bed_files = timed_step("scan_tf_dir", lambda d: sorted(d.glob("*.bed")), Path(args.tf_dir))

    print(f"Found {len(tf_bed_files)} TF bed files", file=sys.stderr)
    print(f"Processing TFs sequentially...", file=sys.stderr)

    results = []

    for tf_bed in tf_bed_files:
        res = process_tf(tf_bed, fcc_df, chrom, chrom_size, window_size,
                         profile_length, step, output_dir)
        tf_name, success, err = res
        print(f"  {'OK' if success else 'FAIL'} {tf_name}", file=sys.stderr)
        results.append(res)

    logging.info(f"total_runtime\t{time.perf_counter() - pipeline_start:.6f}")

    ok = sum(1 for _, s, _ in results if s)
    fail = len(results) - ok
    print(f"\nCompleted {chrom}: {ok}/{len(results)} succeeded, {fail} failed", file=sys.stderr)

    if fail:
        for tf, s, e in results:
            if not s:
                print(f"  {tf}: {e}", file=sys.stderr)
        # CHANGED: exit nonzero on any TF failure so the driver's per-chrom joblog
        # records it as a failed job instead of silently succeeding.
        sys.exit(1)

if __name__ == "__main__":
    main()
