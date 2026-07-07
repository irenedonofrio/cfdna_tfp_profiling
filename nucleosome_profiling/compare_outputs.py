#!/usr/bin/env python3
"""
Categorized comparison of two TF-profiling output trees (golden = current scripts,
new = edited scripts), run on the SAME FCC input.

Because only stage-04 changed here (n_sites plumbing), the expectation is:
  * composite profiles           -> IDENTICAL (n_sites is not stored in the profile files)
  * feature 'n_sites'            -> EXPECTED to change (that was the fix)
  * every other feature column   -> IDENTICAL within --feature-tol

Exit code 0 = all gates pass, 1 = a gate failed (something moved that should not have).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _fmt(x):
    return "nan" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3e}"


def compare_features(golden_f, new_f, tol):
    g = pd.read_csv(golden_f, sep="\t").set_index("TF").sort_index()
    n = pd.read_csv(new_f, sep="\t").set_index("TF").sort_index()

    only_g = sorted(set(g.index) - set(n.index))
    only_n = sorted(set(n.index) - set(g.index))
    if only_g or only_n:
        print("  [WARN] TF set differs between runs")
        if only_g:
            print("         only in golden:", ", ".join(only_g[:10]), "..." if len(only_g) > 10 else "")
        if only_n:
            print("         only in new   :", ", ".join(only_n[:10]), "..." if len(only_n) > 10 else "")

    common = g.index.intersection(n.index)
    g, n = g.loc[common], n.loc[common]

    passed = True

    # --- columns that must NOT change ---
    stable_cols = [c for c in ["mean_coverage", "central_coverage", "amplitude"] if c in g.columns]
    print("\n  Features that must be UNCHANGED (tol = {:g}):".format(tol))
    for c in stable_cols:
        a, b = g[c].to_numpy(float), n[c].to_numpy(float)
        nan_match = np.array_equal(np.isnan(a), np.isnan(b))
        diff = np.nanmax(np.abs(a - b)) if np.isfinite(a).any() else 0.0
        ok = nan_match and (diff <= tol)
        passed &= ok
        print(f"    {c:18s} max|Δ| = {_fmt(diff):>10s}   nan-mask {'ok' if nan_match else 'MISMATCH'}   -> {'PASS' if ok else 'FAIL'}")

    # --- the column we EXPECT to change ---
    if "n_sites" in g.columns:
        changed = (g["n_sites"].to_numpy() != n["n_sites"].to_numpy())
        n_changed = int(changed.sum())
        print(f"\n  n_sites (EXPECTED to change): {n_changed}/{len(common)} TFs differ")
        show = pd.DataFrame({"golden": g["n_sites"], "new": n["n_sites"]})
        show = show[show["golden"] != show["new"]]
        if len(show):
            show["Δ"] = show["new"] - show["golden"]
            with pd.option_context("display.max_rows", 25):
                print(show.head(25).to_string())
            if len(show) > 25:
                print(f"    ... ({len(show) - 25} more)")
        else:
            print("    (no n_sites changed — check that the edited scripts were actually used)")

    return passed


def compare_profiles(golden_dir, new_dir, tol):
    golden_dir, new_dir = Path(golden_dir), Path(new_dir)
    g_files = {p.name: p for p in golden_dir.glob("*_profile.tsv.gz")}
    n_files = {p.name: p for p in new_dir.glob("*_profile.tsv.gz")}
    common = sorted(set(g_files) & set(n_files))

    missing = sorted(set(g_files) ^ set(n_files))
    if missing:
        print("  [WARN] profile files present in only one run:", ", ".join(missing[:10]))

    cols = ["mean_FCC", "smoothed_FCC", "normalized_FCC", "count"]
    worst = 0.0
    worst_where = ""
    nan_ok = True
    n_compared = 0

    for name in common:
        gg = pd.read_csv(g_files[name], sep="\t")
        nn = pd.read_csv(n_files[name], sep="\t")
        if len(gg) != len(nn):
            nan_ok = False
            worst_where = f"{name}: row count {len(gg)} vs {len(nn)}"
            worst = np.inf
            break
        for c in cols:
            if c not in gg.columns or c not in nn.columns:
                continue
            a, b = gg[c].to_numpy(float), nn[c].to_numpy(float)
            if not np.array_equal(np.isnan(a), np.isnan(b)):
                nan_ok = False
                worst_where = f"{name}:{c} (nan-mask mismatch)"
                worst = np.inf
                break
            d = np.nanmax(np.abs(a - b)) if np.isfinite(a).any() else 0.0
            if d > worst:
                worst, worst_where = d, f"{name}:{c}"
        n_compared += 1

    print(f"\n  Composite profiles compared: {n_compared} TFs")
    print(f"    max|Δ| across all profile columns = {_fmt(worst)}   ({worst_where or 'n/a'})")
    ok = nan_ok and (worst <= tol)
    print(f"    nan-masks {'ok' if nan_ok else 'MISMATCH'}   -> {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--golden-composite", required=True)
    ap.add_argument("--new-composite", required=True)
    ap.add_argument("--golden-features", required=True)
    ap.add_argument("--new-features", required=True)
    ap.add_argument("--feature-tol", type=float, default=0.0,
                    help="max allowed |Δ| for unchanged feature columns (default 0 = exact)")
    ap.add_argument("--profile-tol", type=float, default=0.0,
                    help="max allowed |Δ| for composite profile columns (default 0 = exact)")
    args = ap.parse_args()

    print("=" * 70)
    print("GOLDEN COMPARISON  (golden = current scripts, new = edited scripts)")
    print("=" * 70)

    prof_ok = compare_profiles(args.golden_composite, args.new_composite, args.profile_tol)
    feat_ok = compare_features(args.golden_features, args.new_features, args.feature_tol)

    print("\n" + "-" * 70)
    verdict = prof_ok and feat_ok
    print(f"PROFILES : {'PASS' if prof_ok else 'FAIL'}")
    print(f"FEATURES : {'PASS (only n_sites moved)' if feat_ok else 'FAIL (something else moved)'}")
    print(f"OVERALL  : {'PASS' if verdict else 'FAIL'}")
    print("-" * 70)
    sys.exit(0 if verdict else 1)


if __name__ == "__main__":
    main()
