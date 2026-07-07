#!/usr/bin/env python3
"""Prove the optimized GC counting is bit-identical to Griffin's original.

This needs NO bam / reference / data. It loads count_fragment_gc() from the
optimized griffin_GC_counts.py and compares it against a verbatim copy of the
ORIGINAL per-fragment logic over a large set of sequences, deliberately
including every ambiguous IUPAC code and edge cases.

Run:  python3 test_gc_identity.py
"""
import os
import random
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Find the optimized griffin_GC_counts.py. Checked in order:
#   1. the GRIFFIN_GC_COUNTS environment variable, if set
#   2. the same folder as this test
#   3. ../gc_correction/   (the tfp_profiling layout)
_candidates = [
    os.environ.get("GRIFFIN_GC_COUNTS"),
    os.path.join(HERE, "griffin_GC_counts.py"),
    os.path.join(HERE, "..", "gc_correction", "griffin_GC_counts.py"),
]
script_path = next((p for p in _candidates if p and os.path.exists(p)), None)
if script_path is None:
    raise SystemExit(
        "Could not find griffin_GC_counts.py.\n"
        "Put it next to this test or in ../gc_correction/, "
        "or set GRIFFIN_GC_COUNTS=/full/path/to/griffin_GC_counts.py")

SPEC = importlib.util.spec_from_file_location("griffin_gc", script_path)
gc_mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gc_mod)
gc_fast = gc_mod.count_fragment_gc
print("testing script:", os.path.abspath(script_path))

AMBIG = ['N', 'R', 'Y', 'K', 'M', 'B', 'D', 'H', 'V']


def gc_slow(fragment_seq, fragment_start):
    """Verbatim copy of the original griffin_GC_counts.py inner computation."""
    fragment_seq = np.array(list(fragment_seq.upper()))
    fragment_seq[np.isin(fragment_seq, ['A', 'T', 'W'])] = 0
    fragment_seq[np.isin(fragment_seq, ['C', 'G', 'S'])] = 1
    rng = np.random.default_rng(fragment_start)
    fragment_seq[np.isin(fragment_seq, AMBIG)] = rng.integers(
        2, size=len(fragment_seq[np.isin(fragment_seq, AMBIG)]))
    fragment_seq = fragment_seq.astype(float)
    return int(fragment_seq.sum())


def check(seq, fs):
    a, b = gc_slow(seq, fs), gc_fast(seq, fs)
    assert a == b, f"MISMATCH seq={seq!r} start={fs}: original={a} optimized={b}"


def main():
    random.seed(1234)
    n = 0

    # 1) lots of pure-ACGT fragments of varying length and start
    for _ in range(100_000):
        L = random.randint(1, 250)
        seq = ''.join(random.choice('ACGT') for _ in range(L))
        check(seq, random.randint(0, 10**9)); n += 1

    # 2) fragments with soft-masked (lowercase) bases
    for _ in range(20_000):
        L = random.randint(1, 250)
        seq = ''.join(random.choice('acgtACGT') for _ in range(L))
        check(seq, random.randint(0, 10**9)); n += 1

    # 3) fragments containing 1+ ambiguous bases (exercise the RNG fallback)
    for _ in range(50_000):
        L = random.randint(1, 250)
        seq = list(random.choice('ACGT') for _ in range(L))
        for _ in range(random.randint(1, 5)):
            seq[random.randrange(L)] = random.choice(AMBIG + ['S', 'W'])
        check(''.join(seq), random.randint(0, 10**9)); n += 1

    # 4) every ambiguous code on its own, across many seeds
    for code in AMBIG + ['S', 'W']:
        for fs in range(0, 500):
            check(code * random.randint(1, 10), fs); n += 1

    # 5) edge cases
    check('', 0); n += 1
    check('N', 0); n += 1
    check('S', 999); n += 1
    check('ACGTSWNRYKMBDHV', 42); n += 1

    print(f"PASS \u2014 {n:,} fragments, optimized == original on every one")


if __name__ == "__main__":
    main()
