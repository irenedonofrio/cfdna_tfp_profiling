#!/usr/bin/env python3
"""
test_feature_windows.py — pin the three feature windows in 04_merge_and_extract.py
against Griffin's own window construction.

No data needed: runs on a synthetic composite with a known bin grid.

Griffin builds each feature window in two steps (griffin_merge_sites.py):
    :226-228   bounds snapped onto the step grid
               window = [ceil(lo/step)*step, floor(hi/step)*step]
    :232-234   columns = np.arange(window[0], window[1], step)     <-- HALF-OPEN
`arange` excludes the right endpoint, so a closed interval (>= -w & <= w) admits one
extra bin per window. At w=960 / step=15 that is 129 samples instead of 128, and
fft_index=10 then reads a 193.5 bp period instead of Griffin's 192.0 bp.

This test asserts, for save / center / fft:
  A) the number of bins the shipped code actually selects == len(Griffin's arange).
     Read out exactly: on an all-ones composite, np.fft.fft(x)[0] == len(x), so
     `amplitude` with fft_index=0 IS the window length. Each width is fed through
     fft_window in turn, which exercises the snapping rule at all three widths --
     including save_window=1000, which is NOT a multiple of the 15 bp step.
  B) mean_coverage and central_coverage equal a reference computed directly on
     Griffin's arange columns. (A) alone would not catch save_mask or center_mask
     drifting away from fft_mask, since they are three separate expressions.

Exit 0 = all gates pass, 1 = a gate failed.

The module under test defaults to ../nucleosome_profiling/04_merge_and_extract.py;
set TFP_MERGE_MODULE to point elsewhere (used to A/B a patched copy).
"""

import importlib.util
import inspect
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

STEP = 15
WINDOW_HALF = 4995          # composite half-width: floor(5000/15)*15, per 04_profile_tf_sites.py:304-308
SAVE_W, CENTER_W, FFT_W = 1000, 30, 960


def griffin_columns(w, step):
    """Verbatim Griffin construction for a +/-w window: merge_sites.py:226-228 then :232-234."""
    lo = int(np.ceil(-w / step) * step)
    hi = int(np.floor(w / step) * step)
    return np.arange(lo, hi, step)


def load_module():
    default = Path(__file__).resolve().parent.parent / "nucleosome_profiling" / "04_merge_and_extract.py"
    path = Path(os.environ.get("TFP_MERGE_MODULE", default))
    if not path.is_file():
        print(f"[FATAL] module under test not found: {path}", file=sys.stderr)
        sys.exit(1)
    # the module opens a log file at import time; keep it out of the tree
    os.environ.setdefault("TFP_LOG_DIR", tempfile.mkdtemp(prefix="tfp_test_log_"))
    spec = importlib.util.spec_from_file_location("merge_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print(f"module under test: {path}")
    return mod


def call(mod, **kw):
    """Call extract_features, dropping `step` if this version does not accept it.

    Kept deliberately: a revert of the window fix must surface here as a WINDOW LENGTH
    failure, not as an import-time TypeError that hides which gate broke.
    """
    if "step" not in inspect.signature(mod.extract_features).parameters:
        kw.pop("step", None)
    return mod.extract_features(**kw)


def make_composite(values):
    """A composite on the real bin grid: np.arange(-4995, 4995, 15)."""
    pos = np.arange(-WINDOW_HALF, WINDOW_HALF, STEP)
    return pd.DataFrame({
        "position": pos,
        "mean_FCC": values(pos),
        "smoothed_FCC": values(pos),
        "normalized_FCC": values(pos),
        "count": np.ones(len(pos), dtype=int),
    })


def main():
    mod = load_module()
    ones = make_composite(lambda p: np.ones(len(p), dtype=float))
    ramp = make_composite(lambda p: p.astype(float))

    passed = True

    # ---- A) window length, read out through fft_index=0 on an all-ones composite ----
    print("\n  A) window length vs Griffin's np.arange (fft_index=0 on ones -> |FFT[0]| == N)")
    print(f"     {'window':8s} {'w':>5s} {'observed N':>11s} {'Griffin N':>10s}   result")
    for name, w in [("save", SAVE_W), ("center", CENTER_W), ("fft", FFT_W)]:
        feats = call(mod, composite_df=ones, save_window=SAVE_W, center_window=CENTER_W,
                     fft_window=w, fft_index=0, n_sites=None, step=STEP)
        observed = int(round(feats["amplitude"]))
        expected = len(griffin_columns(w, STEP))
        ok = observed == expected
        passed &= ok
        print(f"     {name:8s} {w:5d} {observed:11d} {expected:10d}   {'PASS' if ok else 'FAIL'}")

    # ---- B) the save_mask / center_mask expressions themselves ----
    print("\n  B) mean_coverage / central_coverage vs a reference on Griffin's columns")
    feats = call(mod, composite_df=ramp, save_window=SAVE_W, center_window=CENTER_W,
                 fft_window=FFT_W, fft_index=10, n_sites=None, step=STEP)
    pos = ramp["position"].to_numpy()
    val = ramp["normalized_FCC"].to_numpy()
    for label, key, w in [("mean_coverage", "mean_coverage", SAVE_W),
                          ("central_coverage", "central_coverage", CENTER_W)]:
        ref = val[np.isin(pos, griffin_columns(w, STEP))].mean()
        got = feats[key]
        ok = bool(np.isclose(got, ref, rtol=1e-12, atol=0))
        passed &= ok
        print(f"     {label:18s} got {got:12.6f}   Griffin ref {ref:12.6f}   {'PASS' if ok else 'FAIL'}")

    # ---- C) step is required and is checked against the data ----
    print("\n  C) step guard")
    try:
        mod.extract_features(ones, n_sites=None)
    except TypeError:
        print("     omitting step            -> TypeError                        PASS")
    else:
        passed = False
        print("     omitting step            -> accepted (silent default)        FAIL")
    try:
        mod.extract_features(ones, n_sites=None, step=10)
    except ValueError:
        print("     wrong step (10 vs 15)    -> ValueError                       PASS")
    except TypeError:
        passed = False
        print("     wrong step (10 vs 15)    -> TypeError (step not accepted)    FAIL")
    else:
        passed = False
        print("     wrong step (10 vs 15)    -> accepted silently                FAIL")

    print("\n" + "-" * 70)
    print(f"  OVERALL: {'PASS' if passed else 'FAIL'}")
    print("-" * 70)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
