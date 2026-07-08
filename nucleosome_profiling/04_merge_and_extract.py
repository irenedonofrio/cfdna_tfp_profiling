#!/usr/bin/env python3
"""
Fast TF composite profile computation with parallel processing and NumPy vectorization.
"""

import argparse
import pandas as pd
import numpy as np
import scipy.stats
from concurrent.futures import ProcessPoolExecutor
from scipy.signal import savgol_filter
from pathlib import Path
import functools
import sys
import time
import logging
import os

# CHANGED: log dir configurable via TFP_LOG_DIR with a safe fallback (was a hardcoded
# cluster path opened at import time, which crashed off-cluster / in the container).
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

# -----------------------------------------------------------
# Config loading (CHANGED)
# -----------------------------------------------------------
# CHANGED: optional --config supplies method params. CLI > config > fallback, so this
# stays backward-compatible with callers that pass params explicitly (golden_test.sh).
def load_config_method(config_path):
    if not config_path:
        return {}
    import yaml
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg.get("method", {}) or {}

def resolve(cli_value, config_value, fallback):
    if cli_value is not None:
        return cli_value
    if config_value is not None:
        return config_value
    return fallback

# -----------------------------------------------------------
# Provenance (for the upsert into all_TF_features.tsv)
# -----------------------------------------------------------
# CHANGED: since the features table now accumulates rows across runs (upsert),
# guard against silently mixing rows computed under DIFFERENT method params.
# We stamp a small sidecar with a signature of the config's method: section and
# WARN (non-blocking) if a later run's signature differs. All of this is a no-op
# when --config is not given (e.g. golden_test.sh passes params on the CLI), so
# the fresh-write path stays byte-identical and no sidecar is created.
_PROVENANCE_SIDECAR = ".all_TF_features.method"

def _method_signature(config_path):
    if not config_path:
        return None
    import yaml, json, hashlib
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}
    method = cfg.get("method", {}) or {}
    blob = json.dumps(method, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:12]

def _stamp_provenance(results_dir, config_path):
    sig = _method_signature(config_path)
    if sig is None:
        return
    (Path(results_dir) / _PROVENANCE_SIDECAR).write_text(sig + "\n")

def _warn_provenance_mismatch(results_dir, config_path):
    sig = _method_signature(config_path)
    if sig is None:
        return
    sidecar = Path(results_dir) / _PROVENANCE_SIDECAR
    if sidecar.exists():
        prev = sidecar.read_text().strip()
        if prev != sig:
            print(f"[WARN] method params differ from the run that produced the existing "
                  f"all_TF_features.tsv (sig {prev} -> {sig}). The upserted table will MIX "
                  f"rows computed under different parameters. Use a fresh --results-root if "
                  f"this is a different analysis.", file=sys.stderr)
    _stamp_provenance(results_dir, config_path)  # update to current run's signature

# -----------------------------------------------------------
# Argument parser
# -----------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge TF profiles across chromosomes and compute composite profiles"
    )
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory containing intermediate profile data")
    parser.add_argument("--results_dir", type=str, required=True,
                        help="Directory where extracted features are saved")
    parser.add_argument("--config", default=None,
                        help="Path to config.yaml (supplies step/smoothing/window params)")
    parser.add_argument("--n_jobs", type=int, default=1,
                        help="Number of parallel processes to use (default: 1)")
    # CHANGED: method params default to None so config -> fallback precedence applies.
    parser.add_argument("--step", type=int, default=None, help="Bin size used in profiling (default: 15)")
    parser.add_argument("--smoothing_length", type=int, default=None, help="Fragment length for smoothing (default: 165)")
    parser.add_argument("--save_window", type=int, default=None, help="Half-width for mean coverage (default: 1000)")
    parser.add_argument("--center_window", type=int, default=None, help="Half-width for central coverage (default: 30)")
    parser.add_argument("--fft_window", type=int, default=None, help="Half-width for FFT window (default: 960)")
    parser.add_argument("--fft_index", type=int, default=None, help="FFT frequency index (default: 10)")
    return parser.parse_args()

# -----------------------------------------------------------
# Vectorized profile reconstruction
# -----------------------------------------------------------
def load_and_merge_profiles(intermediate_dir, tf_name):
    """Load and merge per-chromosome TF profiles into site-level profiles using NumPy vectorization.
    The intermediate files are npz fles, with 2d arrays of shape (n_tfbs_per_chrom,n_bins)
    """
    all_uncorrected = []
    all_corrected = []
    total_n_sites = 0   # CHANGED: accumulate the true site count across chromosomes

    tf_files = sorted(intermediate_dir.glob(f"{tf_name}_chr*.npz"))
    if not tf_files:
        raise FileNotFoundError(f"No files found for {tf_name}")

    first_data = np.load(tf_files[0])
    positions = first_data['positions']

    for tf_file in tf_files:
        data = np.load(tf_file)
        all_uncorrected.append(data['uncorrected'])
        all_corrected.append(data['corrected'])
        # CHANGED: sum the per-chromosome true site counts saved by stage 04 (mirrors
        # Griffin's len(current_sites)). Fall back to the row count for older npz files
        # generated before 'n_sites' was added, so we don't crash on stale intermediates.
        total_n_sites += int(data['n_sites']) if 'n_sites' in data.files else data['uncorrected'].shape[0]

    # vertical stack => (total_tfbs, n_bins)
    uncorrected_combined = np.vstack(all_uncorrected)
    corrected_combined = np.vstack(all_corrected)

    return uncorrected_combined, corrected_combined, positions, total_n_sites

# -----------------------------------------------------------
# Outlier masking
# -----------------------------------------------------------
def outlier_mask(uncorrected_combined, corrected_combined):
    """Mask outliers in the profiles using z-score thresholding.
    It is computed on the uncorrected profiles, but applied on all profiles."""
    scores = scipy.stats.zscore(uncorrected_combined, axis=None, nan_policy="omit")
    mask = np.where(scores < 10, 1, np.nan)

    # Apply mask to both
    uncorrected_masked = uncorrected_combined * mask
    corrected_masked = corrected_combined * mask

    # Handle low coverage - check UNCORRECTED
    min_cutoff = 2
    if np.nanmax(uncorrected_masked) < min_cutoff:
        print("Low coverage, resetting cutoff to", min_cutoff, file=sys.stderr)
        mask = np.where(uncorrected_combined <= min_cutoff, 1, np.nan)
        corrected_masked = corrected_combined * mask

    # Only return corrected (as you wanted)
    return corrected_masked


# -----------------------------------------------------------
# Compute composite profile
# -----------------------------------------------------------
def compute_composite_profile(masked_profile, positions,step=15,smoothing_length=165):
    """Compute composite FCC profile - smoothed and normalized"""

    # 1. Average across all TFBS => composite profile
    mean_profile = np.nanmean(masked_profile, axis=0) # this is the composite
    count_profile = np.sum(~np.isnan(masked_profile), axis=0)
    #sem_profile = scipy.stats.sem(masked_profile, axis=0, nan_policy="omit")

    # 2. Savitzky-Golay smoothing
    # savgol window should be approx one fragment length but it must be odd
    savgol_window = np.floor(smoothing_length / step)
    if savgol_window % 2 == 0:
        savgol_window = savgol_window + 1
    savgol_window = int(savgol_window)
    smoothed_profile = savgol_filter(mean_profile, savgol_window, polyorder=3)

    # 3. Scale profile to mean of 1
    profile_mean = np.nanmean(smoothed_profile)
    normalized_profile = smoothed_profile / profile_mean


    result_df = pd.DataFrame({
        "position": positions,
        "mean_FCC": mean_profile,
        "smoothed_FCC": smoothed_profile,
        "normalized_FCC": normalized_profile,
        "count": count_profile
    })
    return result_df

# -----------------------------------------------------------
# Extract features
# -----------------------------------------------------------

def extract_features(composite_df, save_window=1000, center_window=30, fft_window=960, fft_index=10, n_sites=None):
    """
    Extract coverage features from the composite TF profile.
    """
    # Sort by position
    composite_df = composite_df.sort_values('position').reset_index(drop=True)

    # 1. Mean Coverage: over ±1000bp window (or specified save_window)
    save_mask = (composite_df['position'] >= -save_window) & \
                (composite_df['position'] <= save_window)
    mean_coverage = composite_df.loc[save_mask, 'normalized_FCC'].mean()

    # 2. Central Coverage: over ±30bp window
    center_mask = (composite_df['position'] >= -center_window) & \
                  (composite_df['position'] <= center_window)
    central_coverage = composite_df.loc[center_mask, 'normalized_FCC'].mean()

    # 3. Amplitude: FFT over ±960bp window
    fft_mask = (composite_df['position'] >= -fft_window) & \
               (composite_df['position'] <= fft_window)
    fft_values = composite_df.loc[fft_mask, 'normalized_FCC'].values

    # Check we have enough data for FFT
    if len(fft_values) < 2:
        amplitude = np.nan
    else:
        fft_result = np.fft.fft(fft_values)
        amplitude = np.abs(fft_result[fft_index]) if fft_index < len(fft_result) else np.nan

    features = {
        'mean_coverage': mean_coverage,
        'central_coverage': central_coverage,
        'amplitude': amplitude,
        # CHANGED: use the true per-TF site count carried through from stage 04 (mirrors
        # Griffin's len(current_sites)). Previously this read composite_df['count'].iloc[0],
        # i.e. the LEFT-EDGE bin (-window) where the FEWEST sites contribute — a systematic
        # under-count of the site total.
        'n_sites': n_sites
    }

    return features

# -----------------------------------------------------------
# Per-TF processing (runs in parallel)
# -----------------------------------------------------------
def process_single_tf(tf_name, intermediate_dir, final_dir, step=15, smoothing_length=165,
                      save_window=1000, center_window=30, fft_window=960, fft_index=10):
    """Process one TF: load data, compute profile, extract features, save result."""
    try:
        print(f"[INFO] Processing {tf_name}...", file=sys.stderr)

        # CHANGED: also receive the true per-TF site count (summed across chroms).
        uncorrected, corrected, positions, n_sites = timed_step(
            f"{tf_name}:load_merge",
            load_and_merge_profiles,
            intermediate_dir, tf_name
        )

        if uncorrected.size == 0 or corrected.size == 0:
            print(f"  [WARN] No data found for {tf_name}", file=sys.stderr)
            return None

        # Apply outlier masking
        masked_profile = timed_step(
            f"{tf_name}:mask",
            outlier_mask,
            uncorrected, corrected
        )

        # Compute composite profile
        result_df = timed_step(
            f"{tf_name}:composite",
            compute_composite_profile,
            masked_profile,
            positions,
            step=step,
            smoothing_length=smoothing_length
        )

        if result_df is None:
            print(f"  [WARN] Could not compute profile for {tf_name}", file=sys.stderr)
            return None

        # Extract features
        features = timed_step(
            f"{tf_name}:features",
            extract_features,
            result_df,
            save_window=save_window,
            center_window=center_window,
            fft_window=fft_window,
            fft_index=fft_index,
            n_sites=n_sites   # CHANGED: pass the true site count through to the features
        )

        # Add TF name to profile
        result_df["TF"] = tf_name

        # Save profile
        output_file = final_dir / f"{tf_name}_profile.tsv.gz"
        result_df.to_csv(output_file, sep="\t", index=False, compression="gzip")

        print(f"  [OK] Saved {tf_name}", file=sys.stderr)

        features['TF'] = tf_name
        return features

    except Exception as e:
        print(f"[ERROR] {tf_name}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return None

# -----------------------------------------------------------
# Main orchestrator
# -----------------------------------------------------------
def main():
    pipeline_start = time.perf_counter()
    args = parse_args()

    # CHANGED: resolve method params with precedence CLI > config > fallback.
    method = load_config_method(args.config)
    step             = resolve(args.step,             method.get("step"),             15)
    smoothing_length = resolve(args.smoothing_length, method.get("smoothing_length"), 165)
    save_window      = resolve(args.save_window,      method.get("save_window"),      1000)
    center_window    = resolve(args.center_window,    method.get("center_window"),    30)
    fft_window       = resolve(args.fft_window,       method.get("fft_window"),       960)
    fft_index        = resolve(args.fft_index,        method.get("fft_index"),        10)

    output_dir = Path(args.output_dir)
    intermediate_dir = output_dir / "intermediate"
    final_dir = output_dir / "composite"
    final_dir.mkdir(parents=True, exist_ok=True)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    if not intermediate_dir.exists():
        print(f"[FATAL] Intermediate directory not found: {intermediate_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect TF names
    tf_files = list(intermediate_dir.glob("*_chr*.npz"))

    if not tf_files:
        print(f"[FATAL] No NPZ files found in {intermediate_dir}", file=sys.stderr)
        sys.exit(1)

    tf_names = sorted(set(f.stem.rsplit("_chr", 1)[0] for f in tf_files))
    print(f"[INFO] Found {len(tf_names)} TFs to process", file=sys.stderr)

    # Prepare function for executor
    process_func = functools.partial(
        process_single_tf,
        intermediate_dir=intermediate_dir,
        final_dir=final_dir,
        step=step,
        smoothing_length=smoothing_length,
        save_window=save_window,
        center_window=center_window,
        fft_window=fft_window,
        fft_index=fft_index
    )

    # Process TFs
    if args.n_jobs > 1:
        print(f"[INFO] Using {args.n_jobs} parallel jobs", file=sys.stderr)
        with ProcessPoolExecutor(max_workers=args.n_jobs) as executor:
            results = list(executor.map(process_func, tf_names))
    else:
        results = [process_func(tf) for tf in tf_names]

    # Collect features from successful runs
    features_list = [r for r in results if r is not None]

    # Create features summary DataFrame
    if features_list:
        run_df = pd.DataFrame(features_list)
        cols = ['TF', 'mean_coverage', 'central_coverage', 'amplitude', 'n_sites']
        run_df = run_df[cols]
        features_file = results_dir / "all_TF_features.tsv"

        # CHANGED: UPSERT this run's TF rows into an existing table instead of
        # overwriting the whole file. Rows for TFs profiled in a previous run are
        # preserved; a TF profiled again replaces its own row. If no table exists
        # yet, the file is written exactly as before (byte-identical fresh output).
        if features_file.exists():
            _warn_provenance_mismatch(results_dir, args.config)  # non-blocking
            prev = pd.read_csv(features_file, sep="\t")
            keep = prev[~prev['TF'].isin(run_df['TF'])] if 'TF' in prev.columns else prev.iloc[0:0]
            combined = pd.concat([keep, run_df], ignore_index=True)
            combined = combined[cols].sort_values('TF').reset_index(drop=True)
            combined.to_csv(features_file, sep="\t", index=False)
            print(f"[INFO] Upserted {len(run_df)} TF(s) into {features_file} "
                  f"({len(combined)} total after merge)", file=sys.stderr)
        else:
            run_df.to_csv(features_file, sep="\t", index=False)
            _stamp_provenance(results_dir, args.config)
            print(f"[INFO] Saved features summary: {features_file}", file=sys.stderr)

    successful = len(features_list)
    print(f"[INFO] Completed {successful}/{len(tf_names)} TFs successfully", file=sys.stderr)
    print(f"[DONE] Profiles saved in: {final_dir}", file=sys.stderr)
    print(f"[DONE] Features summary: {results_dir}/all_TF_features.tsv", file=sys.stderr)

    logging.info(f"total_runtime\t{time.perf_counter() - pipeline_start:.6f}")

    # CHANGED: exit nonzero if any TF failed to produce features, so the driver can
    # trust the exit code instead of parsing stderr.
    if successful < len(tf_names):
        print(f"[FATAL] {len(tf_names) - successful} TF(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()