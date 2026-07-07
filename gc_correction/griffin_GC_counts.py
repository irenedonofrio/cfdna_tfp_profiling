#!/usr/bin/env python
# coding: utf-8
#
# Optimized drop-in replacement for Griffin's griffin_GC_counts.py
# ----------------------------------------------------------------
# Same command-line interface and same output file
# (<out_dir>/GC_counts/<name>.GC_counts.txt) as the original.
#
# What changed (and why results stay BIT-IDENTICAL):
#   1. Per-fragment GC counting now uses C-level str.count for the common
#      case (fragments with no ambiguous bases) instead of building a NumPy
#      string array, running np.isin three times, and constructing a fresh
#      np.random.default_rng PER FRAGMENT. The original random assignment of
#      ambiguous bases (N/R/Y/K/M/B/D/H/V) is reproduced exactly, and only on
#      the rare fragments that actually contain one. See count_fragment_gc().
#   2. Interval iteration uses itertuples instead of per-row .iloc.
#   3. Work is split into more chunks than CPUs for better load balance
#      (counts are summed back together, so the partition cannot change totals).
#   4. Aggregation uses dict-merge + a single DataFrame build instead of
#      DataFrame.append (which is O(n^2) and was REMOVED in pandas >= 2.0).
#
# Filtering, fragment-coordinate computation, the BAM fetch, and the output
# format/ordering are unchanged from the original.

import os
import sys
import time
import argparse

import numpy as np
import pandas as pd

from multiprocessing import Pool


# Ambiguous IUPAC codes that the original assigns a random 0/1 to.
# (A/T/W -> 0  i.e. "not GC";  C/G/S -> 1  i.e. "GC")
_AMBIGUOUS = ("N", "R", "Y", "K", "M", "B", "D", "H", "V")


def count_fragment_gc(fragment_seq, fragment_start):
    """Number of GC bases in a fragment, bit-identical to Griffin's original.

    Original logic, per base:
        A, T, W            -> 0
        C, G, S            -> 1
        N,R,Y,K,M,B,D,H,V  -> random 0/1, drawn from
                              np.random.default_rng(fragment_start)

    The total GC count therefore equals
        (#C + #G + #S)  +  sum of `n_ambiguous` random draws (each 0 or 1)
    where the draws come from a generator seeded with `fragment_start`.
    The sum of the draws does not depend on which position gets which value,
    so we can compute it directly. For the overwhelmingly common case of a
    fragment with no ambiguous bases, no RNG is constructed at all.
    """
    s = fragment_seq.upper()
    gc = s.count("C") + s.count("G") + s.count("S")
    at = s.count("A") + s.count("T") + s.count("W")
    n_ambiguous = len(s) - gc - at
    if n_ambiguous == 0:
        return gc
    # Rare path: reproduce the original random assignment exactly.
    rng = np.random.default_rng(fragment_start)
    return gc + int(rng.integers(2, size=n_ambiguous).sum())


def _init_worker(bam_path, ref_path, mq, srange):
    """Pool initializer: set per-worker config explicitly.

    This makes the workers correct under BOTH multiprocessing start methods:
    'fork' (Linux / inside the container) and 'spawn' (macOS default), where
    workers re-import the module and would not otherwise inherit these values.
    """
    global _BAM_PATH, _REF_PATH, _MAP_Q, _SIZE_RANGE
    _BAM_PATH, _REF_PATH, _MAP_Q, _SIZE_RANGE = bam_path, ref_path, mq, srange


def collect_reads(sublist):
    """Count fragments by (length, GC content) over a chunk of intervals.

    Config comes from the worker globals set by _init_worker(), so this is safe
    regardless of the multiprocessing start method.
    """
    map_q = _MAP_Q
    size_range = _SIZE_RANGE
    lo, hi = size_range[0], size_range[1]

    # dict[length][num_GC] = count
    GC_dict = {length: {g: 0 for g in range(0, length + 1)}
               for length in range(lo, hi + 1)}

    import pysam  # imported here so the GC helper can be used without pysam installed
    # Opened inside the worker (avoids the truncated-file warning across forks).
    bam_file = pysam.AlignmentFile(_BAM_PATH, "rb")
    ref_seq = pysam.FastaFile(_REF_PATH)

    for chrom, start, end in sublist:
        for read in bam_file.fetch(chrom, start, end):
            tlen = read.template_length
            # size filter: forward reads use +tlen, reverse reads use -tlen
            if read.is_reverse:
                if not (lo <= -tlen <= hi):
                    continue
            else:
                if not (lo <= tlen <= hi):
                    continue
            # qc filters (kept identical to the original)
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

            fragment_seq = ref_seq.fetch(read.reference_name,
                                         fragment_start, fragment_end)
            num_GC = count_fragment_gc(fragment_seq, fragment_start)
            GC_dict[abs(tlen)][num_GC] += 1

    bam_file.close()
    ref_seq.close()
    return GC_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_file', help='sample_bam_file', required=True)
    parser.add_argument('--bam_file_name', help='sample name (does not need to match actual file name)', required=True)
    parser.add_argument('--mappable_regions_path', help='highly mappable regions to be used in GC correction, bedGraph or bed foramt', required=True)
    parser.add_argument('--ref_seq', help='reference sequence (fasta format)', required=True)
    parser.add_argument('--chrom_sizes', help='path to chromosome sizes for the reference seq', required=True)
    parser.add_argument('--out_dir', help='folder for GC bias results', required=True)
    parser.add_argument('--map_q', help='minimum mapping quality for reads to be considered', type=int, required=True)
    parser.add_argument('--size_range', help='range of read sizes to be included', nargs=2, type=int, required=True)
    parser.add_argument('--CPU', help='number of CPU for parallelizing', type=int, required=True)
    args = parser.parse_args()

    bam_file_path = args.bam_file
    bam_file_name = args.bam_file_name
    mappable_regions_path = args.mappable_regions_path
    ref_seq_path = args.ref_seq
    chrom_sizes_path = args.chrom_sizes
    out_dir = args.out_dir
    map_q = args.map_q
    size_range = args.size_range
    CPU = args.CPU

    print('arguments provided:')
    print('\tbam_file_path = "' + bam_file_path + '"')
    print('\tbam_file_name = "' + bam_file_name + '"')
    print('\tmappable_regions_path = "' + mappable_regions_path + '"')
    print('\tref_seq_path = "' + ref_seq_path + '"')
    print('\tchrom_sizes_path = "' + chrom_sizes_path + '"')
    print('\tout_dir = "' + out_dir + '"')
    print('\tmap_q = ' + str(map_q))
    print('\tsize_range = ' + str(size_range))
    print('\tCPU = ' + str(CPU))

    out_file = out_dir + '/GC_counts/' + bam_file_name + '.GC_counts.txt'
    print('out_file', out_file)
    os.makedirs(out_dir + '/GC_counts/', exist_ok=True)

    # import filter
    mappable_intervals = pd.read_csv(mappable_regions_path, sep='\t', header=None)
    # keep autosomes only (chr1..chr22), exactly as the original
    chroms = ['chr' + str(m) for m in range(1, 23)]
    mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]
    print('chroms:', chroms)
    print('number_of_intervals:', len(mappable_intervals))
    sys.stdout.flush()

    # ---- parallel counting ----
    start_time = time.time()

    # Pull the three needed columns into a plain list of (chrom, start, end)
    # tuples once. (Splitting the DataFrame directly is not safe: newer numpy
    # turns np.array_split(df, ...) into a bare ndarray with no .itertuples.)
    intervals = list(mappable_intervals.iloc[:, 0:3].itertuples(index=False, name=None))
    n_intervals = len(intervals)

    # more chunks than CPUs => better load balance across uneven regions.
    # this is just a finer partition, so summed counts are unchanged.
    n_chunks = max(1, min(n_intervals, CPU * 8))
    k, m = divmod(n_intervals, n_chunks)
    sublists = [intervals[i * k + min(i, m):(i + 1) * k + min(i + 1, m)]
                for i in range(n_chunks)]

    with Pool(processes=CPU,
              initializer=_init_worker,
              initargs=(bam_file_path, ref_seq_path, map_q, size_range)) as p:
        GC_dict_list = p.map(collect_reads, sublists)

    print('counting done, aggregating ...', np.round(time.time() - start_time), 's')
    sys.stdout.flush()

    # ---- aggregate (dict merge, identical ordering/format to the original) ----
    lo, hi = size_range[0], size_range[1]
    master = {length: {g: 0 for g in range(0, length + 1)}
              for length in range(lo, hi + 1)}
    for d in GC_dict_list:
        for length, sub in d.items():
            m = master[length]
            for g, c in sub.items():
                m[g] += c

    lengths, gcs, counts = [], [], []
    for length in range(lo, hi + 1):
        for g in range(0, length + 1):
            lengths.append(length)
            gcs.append(g)
            counts.append(master[length][g])

    all_GC_df = pd.DataFrame({'length': lengths,
                              'num_GC': gcs,
                              'number_of_fragments': counts})
    all_GC_df.to_csv(out_file, sep='\t', index=False)

    print('done')
