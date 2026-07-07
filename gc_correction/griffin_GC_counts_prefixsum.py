#!/usr/bin/env python
# coding: utf-8
#
# Prefix-sum drop-in replacement for Griffin's griffin_GC_counts.py
# ----------------------------------------------------------------
# Same CLI and same output file (<out_dir>/GC_counts/<name>.GC_counts.txt),
# same filters / fragment-coordinate math / output ordering as the original.
#
# What changed (and why results stay BIT-IDENTICAL):
#   Instead of ref_seq.fetch() + counting the sequence for EVERY fragment, we
#   build, ONCE per chromosome, a cumulative-GC array:
#         cumGC[i] = number of GC bases (C/G/S) in reference[0:i]
#   Then a fragment's solid GC count is a single subtraction:
#         gc = cumGC[end] - cumGC[start]              (no per-fragment fetch)
#   Ambiguous IUPAC bases (anything not A/C/G/T/S/W) still need Griffin's
#   per-fragment seeded coin flip, which a fixed array cannot encode. So we
#   also keep the sorted positions of ambiguous bases and, for the fragment,
#   count how many fall in its span:
#         n_amb = (#ambiguous positions in [start,end))
#   That count + the fragment_start seed are ALL the original RNG depends on,
#   so we reproduce it exactly:  gc += default_rng(start).integers(2, n_amb).sum()
#   For the overwhelming majority of fragments n_amb==0 and nothing extra runs.
#
#   The per-fragment lookup is vectorized over batches, so once the array is
#   built the counting is essentially array arithmetic.
#
# Parallelism is per-chromosome (each worker builds its own chromosome's array
# and processes that chromosome's intervals).

import os
import sys
import time
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from multiprocessing import Pool

# ASCII codes (we upper-case the reference, matching the original's .upper())
_A, _C, _G, _T, _S, _W = 65, 67, 71, 84, 83, 87

_BATCH = 2_000_000   # fragments per vectorized flush (caps memory)


def _init_worker(bam_path, ref_path, mq, srange, intervals_by_chrom):
    global _BAM, _REF, _MAPQ, _SR, _IVLS
    _BAM, _REF, _MAPQ, _SR, _IVLS = bam_path, ref_path, mq, srange, intervals_by_chrom


def _build_prefix(ref_path, chrom):
    """cumGC (uint32, len N+1) and sorted ambiguous positions for one chromosome."""
    import pysam
    ref = pysam.FastaFile(ref_path)
    seq = ref.fetch(chrom).upper()
    ref.close()
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    N = arr.size

    is_gc = (arr == _C) | (arr == _G) | (arr == _S)
    cumGC = np.zeros(N + 1, dtype=np.uint32)
    np.cumsum(is_gc, dtype=np.uint32, out=cumGC[1:])
    del is_gc

    known = (arr == _A) | (arr == _C) | (arr == _G) | (arr == _T) | (arr == _S) | (arr == _W)
    amb_pos = np.flatnonzero(~known).astype(np.int64)
    del known, arr
    return cumGC, amb_pos, N


def process_chrom(chrom):
    import pysam
    lo, hi = _SR[0], _SR[1]
    map_q = _MAPQ

    t0 = time.time()
    cumGC, amb_pos, N = _build_prefix(_REF, chrom)
    build_t = time.time() - t0

    # length x num_GC tally (dense; tiny: ~ (hi-lo+1) x (hi+1))
    n_len = hi - lo + 1
    M = np.zeros((n_len, hi + 1), dtype=np.int64)
    has_amb = amb_pos.size > 0

    # batch buffers
    b_fs = np.empty(_BATCH, dtype=np.int64)
    b_fe = np.empty(_BATCH, dtype=np.int64)
    b_ln = np.empty(_BATCH, dtype=np.int64)
    n = 0
    count_t = 0.0

    def flush(k):
        nonlocal count_t
        if k == 0:
            return
        tc = time.time()
        fs = b_fs[:k]; fe = b_fe[:k]; ln = b_ln[:k]
        fs_c = np.clip(fs, 0, N)          # match pysam fetch clamping at chrom ends
        fe_c = np.clip(fe, 0, N)
        gc = cumGC[fe_c].astype(np.int64) - cumGC[fs_c].astype(np.int64)
        if has_amb:
            n_amb = (np.searchsorted(amb_pos, fe_c) - np.searchsorted(amb_pos, fs_c))
            for j in np.flatnonzero(n_amb > 0):
                seed = int(fs[j])          # raw fragment_start, exactly as original
                gc[j] += int(np.random.default_rng(seed).integers(2, size=int(n_amb[j])).sum())
        np.add.at(M, (ln - lo, gc), 1)
        count_t += time.time() - tc

    bam = pysam.AlignmentFile(_BAM, "rb")
    for (c, start, end) in _IVLS[chrom]:
        for read in bam.fetch(c, start, end):
            tlen = read.template_length
            if read.is_reverse:
                if not (lo <= -tlen <= hi):
                    continue
            else:
                if not (lo <= tlen <= hi):
                    continue
            if not (read.is_paired
                    and read.mapping_quality >= map_q
                    and not read.is_duplicate
                    and not read.is_qcfail):
                continue
            if read.is_reverse:
                fragment_end = read.reference_start + read.reference_length
                fragment_start = fragment_end + tlen
            else:
                fragment_start = read.reference_start
                fragment_end = read.reference_start + tlen

            b_fs[n] = fragment_start
            b_fe[n] = fragment_end
            b_ln[n] = abs(tlen)
            n += 1
            if n == _BATCH:
                flush(n); n = 0
    bam.close()
    flush(n)
    del cumGC, amb_pos
    return M, build_t, count_t


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_file', required=True)
    parser.add_argument('--bam_file_name', required=True)
    parser.add_argument('--mappable_regions_path', required=True)
    parser.add_argument('--ref_seq', required=True)
    parser.add_argument('--chrom_sizes', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--map_q', type=int, required=True)
    parser.add_argument('--size_range', nargs=2, type=int, required=True)
    parser.add_argument('--CPU', type=int, required=True)
    args = parser.parse_args()

    bam_file_name = args.bam_file_name
    ref_seq_path = args.ref_seq
    out_dir = args.out_dir
    map_q = args.map_q
    size_range = args.size_range
    CPU = args.CPU

    print('arguments provided:')
    print('\tbam_file_path = "' + args.bam_file + '"')
    print('\tbam_file_name = "' + bam_file_name + '"')
    print('\tmappable_regions_path = "' + args.mappable_regions_path + '"')
    print('\tref_seq_path = "' + ref_seq_path + '"')
    print('\tout_dir = "' + out_dir + '"')
    print('\tmap_q = ' + str(map_q))
    print('\tsize_range = ' + str(size_range))
    print('\tCPU = ' + str(CPU))

    out_file = out_dir + '/GC_counts/' + bam_file_name + '.GC_counts.txt'
    print('out_file', out_file)
    os.makedirs(out_dir + '/GC_counts/', exist_ok=True)

    # mappable intervals, autosomes only (chr1..chr22), exactly as the original
    mappable_intervals = pd.read_csv(args.mappable_regions_path, sep='\t', header=None)
    chroms = ['chr' + str(m) for m in range(1, 23)]
    mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]
    print('number_of_intervals:', len(mappable_intervals))
    sys.stdout.flush()

    intervals_by_chrom = defaultdict(list)
    for (c, s, e) in mappable_intervals.iloc[:, 0:3].itertuples(index=False, name=None):
        intervals_by_chrom[c].append((c, int(s), int(e)))
    work_chroms = [c for c in chroms if c in intervals_by_chrom]

    start_time = time.time()
    lo, hi = size_range[0], size_range[1]
    n_len = hi - lo + 1
    master = np.zeros((n_len, hi + 1), dtype=np.int64)
    tot_build = tot_count = 0.0

    with Pool(processes=CPU,
              initializer=_init_worker,
              initargs=(args.bam_file, ref_seq_path, map_q, size_range,
                        dict(intervals_by_chrom))) as p:
        for M, bt, ct in p.map(process_chrom, work_chroms):
            master += M
            tot_build = max(tot_build, bt)   # build runs in parallel; report the longest
            tot_count += ct

    print('counting done, aggregating ...', np.round(time.time() - start_time), 's')
    print('  (slowest per-chromosome prefix build: %.1fs; summed vectorized-lookup time: %.1fs)'
          % (tot_build, tot_count))
    sys.stdout.flush()

    # ---- output: identical format/ordering to the original ----
    lengths, gcs, counts = [], [], []
    for length in range(lo, hi + 1):
        row = master[length - lo]
        for g in range(0, length + 1):
            lengths.append(length)
            gcs.append(g)
            counts.append(int(row[g]))

    all_GC_df = pd.DataFrame({'length': lengths,
                              'num_GC': gcs,
                              'number_of_fragments': counts})
    all_GC_df.to_csv(out_file, sep='\t', index=False)
    print('done')
