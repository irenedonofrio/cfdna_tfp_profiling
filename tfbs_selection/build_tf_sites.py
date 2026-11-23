#!/usr/bin/env python3
"""
build_tf_sites.py  -  Definitive TF list + per-TF site files from GTRD.

Upstream site-generation stage for cfdna_tfp_profiling. Produces the TF
universe and the per-TF site files that the profiling stage reads from
refs/TFBS_10000ms/.

Pipeline (single streaming pass over the meta-clusters file):
  1. Open GTRD meta-clusters: the .zip exactly as downloaded (preferred), or an
     extracted .interval / .interval.gz. Reading the .zip directly avoids
     silently using a partial extraction.
  2. Locate columns BY HEADER NAME (GTRD reorders columns between releases).
  3. Drop non-TF targets using a "real TF" gate (--census-type):
       cisbp        keep CIS-BP TF_Status != N   (default; matches Griffin/De Sarkar)
       lambert      keep Lambert 'Is TF?' == Yes (DBD-based census)
       intersection keep only symbols in BOTH, and report how much they disagree
     This is the ONLY "real TF" gate -- no site-level mappability filtering
     (mappability is handled in the profiling stage, as in Griffin).
  4. Keep autosomes (chr1-22); site position = (Start + End) // 2 (Griffin rule).
  5. Integrity check: every autosome must be present, reach (close to) its
     hg38 end, and the file must end on a complete line. Otherwise stop
     (override with --allow-incomplete).
  6. Threshold: keep TFs with >= --min-sites autosomal sites.
  7. Per surviving TF keep the --top-n sites, ranked by
         peak.count (desc)  ->  number of experiments (desc)  ->  seeded hash
     peak.count is the number of individual peaks (experiments x peak callers)
     merged into the meta-cluster. The hash breaks the remaining ties in a way
     that is reproducible and independent of input order, so tied sites are
     not taken in file order (which is sorted by chromosome and would bias
     the selection towards chr1, chr10, chr11, ...).

Outputs (under --outdir):
  sites/<TF>.bed       chrom  pos  pos+1  TF  peak.count  .   (BED6, no header;
                       the profiling stage reads exactly six columns)
  tf_list.tsv          one row per gated TF: site counts, status, cutoff diagnostics
  dropped_targets.tsv  GTRD targets that failed the TF gate (audit trail)
  chrom_summary.tsv    per chromosome: input rows, furthest site, kept sites

tf_list.tsv columns:
  tf, n_autosomal_sites, n_written, status,
  n_chroms_written     autosomes represented among the written sites
  cutoff_peak_count    peak.count of the lowest-ranked written site
  cutoff_n_exp         its number of experiments
  n_cutoff_kept        written sites with exactly that (peak.count, n_exp)
  n_cutoff_total       all sites of the TF with that (peak.count, n_exp)
  When n_cutoff_total > n_cutoff_kept, the hash chose which of those tied
  sites were written; n_cutoff_kept is how many sites that choice affects.

--------------------------------------------------------------------------
CHECK before trusting the output: the script prints the header->column
mapping it inferred and a few sample rows. If auto-detection guesses wrong,
override with --tf-col / --chrom-col / --start-col / --end-col /
--peakcount-col / --nexp-col.
--------------------------------------------------------------------------

Example (default CIS-BP gate):
  python build_tf_sites.py \\
      --metaclusters Homo_sapiens_meta_clusters.zip \\
      --census-type cisbp --cisbp TF_Information.txt \\
      --outdir ./tf_sites --min-sites 10000 --top-n 10000

  # Lambert gate:                    --census-type lambert --lambert DatabaseExtract_v_1.01.csv
  # both, with a divergence report:  --census-type intersection --cisbp ... --lambert ...
"""

import argparse
import csv
import gzip
import heapq
import io
import os
import sys
import zipfile
import zlib
from collections import Counter, defaultdict


# ---- header/column detection -------------------------------------------------

# Candidate names per logical column, in preference order. Matched
# case-insensitively against the header after stripping a leading '#'.
META_ALIASES = {
    "chrom":     ["chrom", "chr", "chromosome", "seqnames"],
    "start":     ["start", "chromstart", "begin"],
    "end":       ["end", "chromend", "stop", "finish"],
    "tf":        ["tftitle", "tf", "target", "tfname", "gene", "genesymbol", "uniprotid"],
    "peakcount": ["peak.count", "peakcount", "peak_count", "peaks"],
}
# Optional: number of experiments supporting a site. GTRD 20.06 has 'exp.set'
# (';'-separated experiment IDs); tables with an 'exp.count' column also work.
META_OPTIONAL = {
    "nexp": ["exp.count", "exp_count", "exp.set", "exp_set"],
}
LAMBERT_ALIASES = {
    "symbol": ["hgnc symbol", "hgncsymbol", "symbol", "gene symbol", "genesymbol", "gene"],
    "is_tf":  ["is tf?", "is tf", "istf", "is_tf"],
}
CISBP_ALIASES = {   # species handled separately (optional column)
    "symbol": ["tf_name", "gene_name", "tf name", "gene name"],
    "status": ["tf_status", "tf status", "status"],
}

# hg38 autosome lengths, used only for the completeness check.
HG38_AUTOSOMES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}
MIN_END_FRACTION = 0.90   # furthest site on each autosome must reach 90% of its length


def _norm(s):
    return s.strip().lstrip("#").strip().lower()


def resolve_columns(header_fields, aliases, overrides):
    """Map each logical column to an index in header_fields."""
    lut = {_norm(name): i for i, name in enumerate(header_fields)}
    resolved = {}
    for key, candidates in aliases.items():
        if overrides.get(key):  # explicit --*-col wins
            want = _norm(overrides[key])
            if want not in lut:
                sys.exit(f"[error] requested column '{overrides[key]}' for '{key}' "
                         f"not in header: {header_fields}")
            resolved[key] = lut[want]
            continue
        hit = next((lut[c] for c in candidates if c in lut), None)
        if hit is None:
            sys.exit(f"[error] could not auto-detect column '{key}'. "
                     f"Tried {candidates}. Header was: {header_fields}\n"
                     f"        Override with --{key}-col.")
        resolved[key] = hit
    return resolved


def resolve_optional(header_fields, candidates, override):
    """Index of an optional column, or None if absent."""
    lut = {_norm(name): i for i, name in enumerate(header_fields)}
    if override:
        want = _norm(override)
        if want not in lut:
            sys.exit(f"[error] requested column '{override}' not in header: {header_fields}")
        return lut[want]
    return next((lut[c] for c in candidates if c in lut), None)


def open_meta(path, member=None):
    """Text stream over the meta-clusters table: .zip, .gz or plain text."""
    if path.endswith(".zip"):
        try:
            zf = zipfile.ZipFile(path)
        except zipfile.BadZipFile as e:
            sys.exit(f"[error] {path} is not a readable zip ({e}). "
                     f"The download is probably incomplete -- download it again.")
        files = [i for i in zf.infolist() if not i.is_dir()]
        if member:
            info = zf.getinfo(member)
        else:
            cands = [i for i in files if i.filename.endswith(".interval")] or files
            if len(cands) != 1:
                sys.exit(f"[error] {path} has {len(cands)} candidate members; "
                         f"pick one with --zip-member: {[i.filename for i in cands]}")
            info = cands[0]
        print(f"[meta] reading '{info.filename}' from zip "
              f"({info.file_size / 1e9:.1f} GB uncompressed)", file=sys.stderr)
        # zipfile verifies the member's CRC when the end is reached, so a
        # corrupt archive raises instead of yielding a short table.
        return io.TextIOWrapper(zf.open(info), encoding="utf-8")
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


# ---- census ------------------------------------------------------------------

def load_lambert(path):
    """UPPERCASED symbols with 'Is TF?' == Yes (Lambert 2018 census CSV)."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        cols = resolve_columns(header, LAMBERT_ALIASES, {})
        si, ti = cols["symbol"], cols["is_tf"]
        tfs = {row[si].strip().upper() for row in reader
               if len(row) > max(si, ti)
               and row[ti].strip().lower() == "yes" and row[si].strip()}
    if not tfs:
        sys.exit("[error] Lambert census produced 0 TFs -- check the file/columns.")
    print(f"[lambert] {len(tfs)} symbols with 'Is TF?' == Yes", file=sys.stderr)
    return tfs


def load_cisbp(path):
    """UPPERCASED symbols with TF_Status != N (CIS-BP TF_Information.txt, TSV).

    If a species column is present, rows are restricted to Homo_sapiens.
    """
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cols = resolve_columns(header, CISBP_ALIASES, {})
        si, ki = cols["symbol"], cols["status"]
        lut = {_norm(h): i for i, h in enumerate(header)}
        pi = next((lut[c] for c in ("tf_species", "species") if c in lut), None)
        tfs = set()
        for row in reader:
            if len(row) <= max(si, ki):
                continue
            if pi is not None and len(row) > pi \
               and row[pi].strip().lower().replace(" ", "_") != "homo_sapiens":
                continue
            if row[ki].strip().upper() != "N":
                sym = row[si].strip().upper()
                if sym:
                    tfs.add(sym)
    if not tfs:
        sys.exit("[error] CIS-BP produced 0 TFs -- check the file/columns.")
    print(f"[cisbp] {len(tfs)} symbols with TF_Status != N", file=sys.stderr)
    return tfs


def resolve_tf_set(args):
    """Pick the 'real TF' gate. Default cisbp matches Griffin/De Sarkar."""
    if args.census_type == "lambert":
        return load_lambert(args.lambert)
    if args.census_type == "cisbp":
        return load_cisbp(args.cisbp)
    # intersection: gate on the overlap, and report how much the two disagree
    L, C = load_lambert(args.lambert), load_cisbp(args.cisbp)
    inter = L & C
    only_l, only_c = sorted(L - C), sorted(C - L)
    print(f"[compare] lambert={len(L)} cisbp={len(C)} intersection={len(inter)}",
          file=sys.stderr)
    print(f"[compare] lambert-only n={len(only_l)} e.g. {only_l[:8]}", file=sys.stderr)
    print(f"[compare] cisbp-only   n={len(only_c)} e.g. {only_c[:8]}", file=sys.stderr)
    return inter


# ---- meta_clusters -----------------------------------------------------------

AUTOSOMES = {str(i) for i in range(1, 23)}


def canon_chrom(chrom):
    """'chr7' / 'Chr7' / '7' -> 'chr7'."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "chr" + c


def is_autosome(chrom):
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return c in AUTOSOMES


def stream_sites(meta_path, tf_symbols, overrides, nexp_override, top_n, seed,
                 member=None):
    """
    One pass over the meta-clusters table. Non-TF targets and non-autosomes are
    dropped immediately. For each TF only the current top-N candidates are held
    (min-heap), so memory stays at ~ n_TFs x top_n instead of every site.

    Rank key (higher is better), a total order:
        (peak.count, n_experiments, crc32(seed|TF|chrom|start|end), chrom, start, end)
    """
    heaps = defaultdict(list)            # TF -> min-heap of rank keys (size <= top_n)
    n_sites = defaultdict(int)           # TF -> autosomal sites
    pair_hist = defaultdict(Counter)     # TF -> Counter((peak, n_exp)) for rows that
                                         #       could still be at the final cutoff
    dropped_counts = defaultdict(int)    # target -> autosomal sites, for non-TF targets
    chrom_rows = Counter()               # all rows (any target) per chromosome
    chrom_max_end = {}                   # furthest END seen per chromosome
    chrom_canon = {}                     # share one str object per chromosome name
    skipped_bad = 0
    last_line_complete = True
    n_lines = 0

    with open_meta(meta_path, member) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = resolve_columns(header, META_ALIASES, overrides)
        ci, si, ei, ti, pi = (cols["chrom"], cols["start"], cols["end"],
                              cols["tf"], cols["peakcount"])
        xi = resolve_optional(header, META_OPTIONAL["nexp"], nexp_override)
        x_is_set = xi is not None and _norm(header[xi]).endswith("set")
        maxcol = max(ci, si, ei, ti, pi, xi if xi is not None else 0)
        print(f"[meta] column map: chrom={header[ci]!r} start={header[si]!r} "
              f"end={header[ei]!r} tf={header[ti]!r} peakcount={header[pi]!r} "
              f"n_exp={(header[xi] if xi is not None else None)!r}"
              f"{' (counted from the set)' if x_is_set else ''}",
              file=sys.stderr)
        if xi is None:
            print("[meta] no experiment column found: ties on peak.count go "
                  "straight to the seeded hash", file=sys.stderr)

        shown = 0
        for line in fh:
            n_lines += 1
            last_line_complete = line.endswith("\n")
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= maxcol:
                skipped_bad += 1
                continue
            chrom = parts[ci]
            chrom = chrom_canon.setdefault(chrom, chrom)
            chrom_rows[chrom] += 1
            try:
                start = int(parts[si]); end = int(parts[ei])
                peak = int(float(parts[pi]))
            except ValueError:
                skipped_bad += 1
                continue
            if end > chrom_max_end.get(chrom, -1):
                chrom_max_end[chrom] = end
            if not is_autosome(chrom):
                continue
            tf = parts[ti].strip().upper()
            if not tf:
                skipped_bad += 1
                continue

            if tf not in tf_symbols:          # not a real TF -> drop now
                dropped_counts[tf] += 1
                continue

            n_sites[tf] += 1
            heap = heaps[tf]
            # A site below the current N-th best peak.count can never be
            # written nor sit at the final cutoff (the cutoff only rises).
            if len(heap) >= top_n and peak < heap[0][0]:
                continue

            if xi is None:
                n_exp = 0
            elif x_is_set:
                s = parts[xi].strip()
                n_exp = s.count(";") + 1 if s else 0
            else:
                try:
                    n_exp = int(float(parts[xi]))
                except ValueError:
                    n_exp = 0
            h = zlib.crc32(f"{seed}|{tf}|{chrom}|{start}|{end}".encode())
            key = (peak, n_exp, h, chrom, start, end)
            pair_hist[tf][(peak, n_exp)] += 1
            if len(heap) < top_n:
                heapq.heappush(heap, key)
            elif key > heap[0]:
                heapq.heapreplace(heap, key)

            if shown < 3:                     # sample rows for the user to eyeball
                print(f"[meta] sample: {chrom} {start}-{end} pos={(start + end) // 2} "
                      f"tf={tf} peak.count={peak} n_exp={n_exp}", file=sys.stderr)
                shown += 1

    print(f"[meta] {n_lines:,} data rows read", file=sys.stderr)
    if skipped_bad:
        print(f"[meta] skipped {skipped_bad} malformed/blank rows", file=sys.stderr)
    return {
        "heaps": heaps, "n_sites": n_sites, "pair_hist": pair_hist,
        "dropped_counts": dropped_counts, "chrom_rows": chrom_rows,
        "chrom_max_end": chrom_max_end, "last_line_complete": last_line_complete,
    }


def check_completeness(res, allow_incomplete):
    """Stop if the input is visibly truncated (missing/short autosomes, cut last line)."""
    rows = {canon_chrom(c): n for c, n in res["chrom_rows"].items()}
    max_end = {canon_chrom(c): e for c, e in res["chrom_max_end"].items()}
    problems = []
    missing = [c for c in HG38_AUTOSOMES if rows.get(c, 0) == 0]
    if missing:
        problems.append(f"no rows at all on {', '.join(missing)}")
    short = [f"{c} (furthest site {max_end[c] / 1e6:.1f} of {L / 1e6:.1f} Mb)"
             for c, L in HG38_AUTOSOMES.items()
             if c in max_end and max_end[c] < MIN_END_FRACTION * L]
    if short:
        problems.append("sites stop well before the chromosome end on " + ", ".join(short))
    if not res["last_line_complete"]:
        problems.append("the last line has no newline (file cut off mid-line)")

    if not problems:
        print(f"[check] input complete: all 22 autosomes present and reach "
              f">= {MIN_END_FRACTION:.0%} of their length; last line intact",
              file=sys.stderr)
        return
    msg = ("meta-clusters input looks INCOMPLETE:\n  - " + "\n  - ".join(problems) +
           "\n  Download the zip again and pass it directly (do not extract it).")
    if allow_incomplete:
        print("[warning] " + msg + "\n  Continuing because --allow-incomplete was given.",
              file=sys.stderr)
    else:
        sys.exit("[error] " + msg + "\n  (Override with --allow-incomplete.)")


# ---- write -------------------------------------------------------------------

def write_outputs(res, outdir, min_sites):
    heaps, n_sites, pair_hist = res["heaps"], res["n_sites"], res["pair_hist"]
    dropped_counts = res["dropped_counts"]
    sites_dir = os.path.join(outdir, "sites")
    os.makedirs(sites_dir, exist_ok=True)

    tf_rows, kept = [], 0
    kept_per_chrom = Counter()
    for tf in sorted(n_sites):
        n_total = n_sites[tf]
        if n_total < min_sites:
            tf_rows.append((tf, n_total, 0, f"dropped: <{min_sites} autosomal sites",
                            "", "", "", "", ""))
            continue
        top = sorted(heaps[tf], reverse=True)      # best first
        with open(os.path.join(sites_dir, f"{tf}.bed"), "w") as out:
            for peak, _n_exp, _h, chrom, start, end in top:
                pos = (start + end) // 2           # Griffin: mean of Start/End
                # BED6: chrom, start, end, name, score=peak.count, strand(.)
                # Exactly six columns: 04_profile_tf_sites.py reads the file with six
                # column names and fails on any extra column.
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{tf}\t{peak}\t.\n")
                kept_per_chrom[canon_chrom(chrom)] += 1
        cut_peak, cut_exp = top[-1][0], top[-1][1]
        n_cut_kept = sum(1 for k in top if k[0] == cut_peak and k[1] == cut_exp)
        n_cut_total = pair_hist[tf][(cut_peak, cut_exp)]
        n_chroms = len({canon_chrom(k[3]) for k in top})
        tf_rows.append((tf, n_total, len(top), "kept", n_chroms,
                        cut_peak, cut_exp, n_cut_kept, n_cut_total))
        kept += 1

    with open(os.path.join(outdir, "tf_list.tsv"), "w") as out:
        out.write("tf\tn_autosomal_sites\tn_written\tstatus\tn_chroms_written\t"
                  "cutoff_peak_count\tcutoff_n_exp\tn_cutoff_kept\tn_cutoff_total\n")
        for row in tf_rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    # audit: targets present in GTRD but absent from the census TF list
    with open(os.path.join(outdir, "dropped_targets.tsv"), "w") as out:
        out.write("target\tn_autosomal_sites\treason\n")
        for tgt in sorted(dropped_counts, key=lambda t: (-dropped_counts[t], t)):
            out.write(f"{tgt}\t{dropped_counts[tgt]}\tnot in TF gate set\n")

    # per-chromosome QC
    rows = Counter()
    max_end = {}
    for c, n in res["chrom_rows"].items():
        rows[canon_chrom(c)] += n
    for c, e in res["chrom_max_end"].items():
        max_end[canon_chrom(c)] = max(e, max_end.get(canon_chrom(c), -1))

    def chrom_order(c):
        x = c[3:]
        return (0, int(x), "") if x.isdigit() else (1, 0, x)

    with open(os.path.join(outdir, "chrom_summary.tsv"), "w") as out:
        out.write("chrom\tlength_hg38\tinput_rows\tfurthest_site_end\t"
                  "frac_of_length\tkept_sites\n")
        for c in sorted(set(rows) | set(HG38_AUTOSOMES), key=chrom_order):
            L = HG38_AUTOSOMES.get(c)
            e = max_end.get(c)
            frac = f"{e / L:.3f}" if (L and e is not None) else ""
            out.write(f"{c}\t{L or ''}\t{rows.get(c, 0)}\t{'' if e is None else e}\t"
                      f"{frac}\t{kept_per_chrom.get(c, 0)}\n")

    n_gated = len(n_sites)
    print(f"[done] GTRD targets on autosomes: {n_gated + len(dropped_counts)}  |  "
          f"failed TF gate: {len(dropped_counts)}  |  passed: {n_gated}  |  "
          f"<{min_sites} sites: {n_gated - kept}  |  kept: {kept}", file=sys.stderr)
    print(f"[done] site files in {sites_dir}", file=sys.stderr)


# ---- main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Build definitive TF list + site files from GTRD.")
    p.add_argument("--metaclusters", required=True,
                   help="GTRD meta-clusters: the downloaded .zip (preferred) or .interval(.gz)")
    p.add_argument("--zip-member", help="member to read if the zip holds several files")
    p.add_argument("--outdir", required=True)
    p.add_argument("--census-type", choices=["cisbp", "lambert", "intersection"],
                   default="cisbp",
                   help="'real TF' gate. cisbp (default) matches Griffin/De Sarkar; "
                        "lambert is the DBD-based census; intersection uses both.")
    p.add_argument("--cisbp", help="CIS-BP TF_Information.txt (TSV)")
    p.add_argument("--lambert", help="Lambert 2018 census CSV")
    p.add_argument("--min-sites", type=int, default=10000)
    p.add_argument("--top-n", type=int, default=10000)
    p.add_argument("--seed", type=int, default=0,
                   help="seed for the hash that breaks remaining ties (default 0)")
    p.add_argument("--allow-incomplete", action="store_true",
                   help="continue even if the input looks truncated (not recommended)")
    # meta_clusters column overrides (optional; auto-detected by default)
    p.add_argument("--tf-col"); p.add_argument("--chrom-col")
    p.add_argument("--start-col"); p.add_argument("--end-col")
    p.add_argument("--peakcount-col"); p.add_argument("--nexp-col")
    args = p.parse_args()

    need = {"cisbp": ["cisbp"], "lambert": ["lambert"],
            "intersection": ["cisbp", "lambert"]}[args.census_type]
    missing = [f"--{x}" for x in need if getattr(args, x) is None]
    if missing:
        p.error(f"--census-type {args.census_type} requires: {', '.join(missing)}")
    if args.top_n < 1 or args.min_sites < 0:
        p.error("--top-n must be >= 1 and --min-sites >= 0")

    meta_ovr = {"chrom": args.chrom_col, "start": args.start_col, "end": args.end_col,
                "tf": args.tf_col, "peakcount": args.peakcount_col}

    tf_symbols = resolve_tf_set(args)
    try:
        res = stream_sites(args.metaclusters, tf_symbols, meta_ovr, args.nexp_col,
                           args.top_n, args.seed, args.zip_member)
    except (zipfile.BadZipFile, zlib.error, EOFError, gzip.BadGzipFile) as e:
        sys.exit(f"[error] reading {args.metaclusters} failed: {e}\n"
                 f"        The file is corrupt or incomplete -- download it again.")
    check_completeness(res, args.allow_incomplete)
    os.makedirs(args.outdir, exist_ok=True)
    write_outputs(res, args.outdir, args.min_sites)


if __name__ == "__main__":
    main()
