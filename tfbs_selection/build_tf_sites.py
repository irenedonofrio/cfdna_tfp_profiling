#!/usr/bin/env python3
"""
build_tf_sites.py  -  Definitive TF list + per-TF site files from GTRD.

Upstream site-generation stage for cfdna_tfp_profiling. Produces the TF
universe and per-TF site files that feed Stage 1 (GC correction).

Pipeline (single streaming pass over the meta_clusters file):
  2. Parse GTRD meta_clusters.interval(.gz), locating columns BY HEADER NAME.
  3. Drop non-TF targets on the fly using a "real TF" gate (--census-type):
       cisbp        keep CIS-BP TF_Status != N   (default; matches Griffin/De Sarkar)
       lambert      keep Lambert 'Is TF?' == Yes (stricter, DBD-based)
       intersection keep only symbols in BOTH, and report how much they disagree
     This is the ONLY "real TF" gate -- no mappability filtering.
  4. Keep autosomes (chr1-22); site position = mean(Start, End) (Griffin rule).
  5. Threshold: keep TFs with >= --min-sites autosomal sites.
  6. Per surviving TF: rank sites by 'peak.count' (desc), take --top-n, write.

Outputs (under --outdir):
  sites/<TF>.txt      Chrom  Start  End  position  peak.count   (tab, w/ header)
  tf_list.tsv         tf  n_autosomal_sites  n_written  status
  dropped_targets.tsv target  n_autosomal_sites  reason   (audit trail)

--------------------------------------------------------------------------
TWO THINGS TO CONFIRM before trusting the output:
  1. COLUMN DETECTION. The script prints the header->column mapping it inferred
     and a couple of sample rows. Eyeball them. GTRD renamed/repackaged files
     across 18.01 / 19.10 / 20.06; if auto-detect guesses wrong, override with
     --tf-col / --chrom-col / --start-col / --end-col / --peakcount-col.
  2. GRIFFIN BYTE-IDENTITY. Two choices here mirror Griffin but are worth
     checking against Griffin's published beds during validation:
       - position = (Start + End) // 2  (integer floor of the mean)
       - top-N tie-breaking is STABLE (input order kept on equal peak.count)
     If your CDX2 comparison diverges on a handful of sites, look here first.
--------------------------------------------------------------------------

Example (default CIS-BP gate):
  python build_tf_sites.py \
      --metaclusters "Homo sapiens_meta_clusters.interval.gz" \
      --census-type cisbp --cisbp TF_Information.txt \
      --outdir ./tf_sites --min-sites 10000 --top-n 10000

  # stricter Lambert gate:      --census-type lambert  --lambert DatabaseExtract_v_1.01.csv
  # both, with a divergence report: --census-type intersection --cisbp ... --lambert ...
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict


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
LAMBERT_ALIASES = {
    "symbol": ["hgnc symbol", "hgncsymbol", "symbol", "gene symbol", "genesymbol", "gene"],
    "is_tf":  ["is tf?", "is tf", "istf", "is_tf"],
}
CISBP_ALIASES = {   # species handled separately (optional column)
    "symbol": ["tf_name", "gene_name", "tf name", "gene name"],
    "status": ["tf_status", "tf status", "status"],
}


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
                     f"        Override with --{key.replace('peakcount','peakcount')}-col.")
        resolved[key] = hit
    return resolved


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


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


def is_autosome(chrom):
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return c in AUTOSOMES


def stream_sites(meta_path, tf_symbols, overrides):
    """
    One pass. Drop non-TF targets and non-autosomes immediately.
    Returns dict: TF -> list of (chrom, start, end, position, peak_count),
    in file order (so a later stable sort keeps input order on ties).
    """
    per_tf = defaultdict(list)
    dropped_counts = defaultdict(int)   # target -> autosomal-site count, for targets NOT in census
    skipped_bad = 0

    with open_maybe_gz(meta_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = resolve_columns(header, META_ALIASES, overrides)
        ci, si, ei, ti, pi = (cols["chrom"], cols["start"], cols["end"],
                              cols["tf"], cols["peakcount"])
        print(f"[meta] column map: chrom={header[ci]!r} start={header[si]!r} "
              f"end={header[ei]!r} tf={header[ti]!r} peakcount={header[pi]!r}",
              file=sys.stderr)

        shown = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ci, si, ei, ti, pi):
                skipped_bad += 1
                continue
            chrom = parts[ci]
            if not is_autosome(chrom):
                continue
            tf = parts[ti].strip().upper()
            if not tf:
                skipped_bad += 1
                continue
            try:
                start = int(parts[si]); end = int(parts[ei])
                peak = int(float(parts[pi]))
            except ValueError:
                skipped_bad += 1
                continue

            if tf not in tf_symbols:          # not a real TF -> drop now
                dropped_counts[tf] += 1
                continue

            pos = (start + end) // 2          # Griffin: mean of Start/End
            per_tf[tf].append((chrom, start, end, pos, peak))

            if shown < 3:                     # sample rows for the user to eyeball
                print(f"[meta] sample: {chrom} {start}-{end} pos={pos} "
                      f"tf={tf} peak.count={peak}", file=sys.stderr)
                shown += 1

    if skipped_bad:
        print(f"[meta] skipped {skipped_bad} malformed/blank rows", file=sys.stderr)
    return per_tf, dropped_counts


# ---- write -------------------------------------------------------------------

def write_outputs(per_tf, dropped_counts, outdir, min_sites, top_n):
    sites_dir = os.path.join(outdir, "sites")
    os.makedirs(sites_dir, exist_ok=True)

    tf_rows, kept = [], 0
    for tf in sorted(per_tf):
        sites = per_tf[tf]
        n_total = len(sites)
        if n_total < min_sites:
            tf_rows.append((tf, n_total, 0, f"dropped: <{min_sites} autosomal sites"))
            continue
        # stable sort by peak.count desc -> ties keep file order
        top = sorted(sites, key=lambda r: r[4], reverse=True)[:top_n]
        path = os.path.join(sites_dir, f"{tf}.txt")
        path = os.path.join(sites_dir, f"{tf}.bed")
        with open(path, "w") as out:
            for chrom, start, end, pos, peak in top:
                # BED6+1: chrom, start, end, name, score(.), strand(.), peak.count
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{tf}\t.\t.\t{peak}\n")
        tf_rows.append((tf, n_total, len(top), "kept"))
        kept += 1

    with open(os.path.join(outdir, "tf_list.tsv"), "w") as out:
        out.write("tf\tn_autosomal_sites\tn_written\tstatus\n")
        for tf, n_total, n_written, status in tf_rows:
            out.write(f"{tf}\t{n_total}\t{n_written}\t{status}\n")

    # audit: targets present in GTRD but absent from the census TF list
    with open(os.path.join(outdir, "dropped_targets.tsv"), "w") as out:
        out.write("target\tn_autosomal_sites\treason\n")
        for tgt in sorted(dropped_counts, key=lambda t: -dropped_counts[t]):
            out.write(f"{tgt}\t{dropped_counts[tgt]}\tnot in TF gate set\n")

    print(f"[done] {kept} TFs kept, site files in {sites_dir}", file=sys.stderr)
    print(f"[done] {len(dropped_counts)} non-TF targets logged to dropped_targets.tsv",
          file=sys.stderr)


# ---- main --------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Build definitive TF list + site files from GTRD.")
    p.add_argument("--metaclusters", required=True, help="GTRD meta_clusters.interval(.gz)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--census-type", choices=["cisbp", "lambert", "intersection"],
                   default="cisbp",
                   help="'real TF' gate. cisbp (default) matches Griffin/De Sarkar; "
                        "lambert is the stricter DBD-based test; intersection uses both.")
    p.add_argument("--cisbp", help="CIS-BP TF_Information.txt (TSV)")
    p.add_argument("--lambert", help="Lambert 2018 census CSV")
    p.add_argument("--min-sites", type=int, default=10000)
    p.add_argument("--top-n", type=int, default=10000)
    # meta_clusters column overrides (optional; auto-detected by default)
    p.add_argument("--tf-col"); p.add_argument("--chrom-col")
    p.add_argument("--start-col"); p.add_argument("--end-col")
    p.add_argument("--peakcount-col")
    args = p.parse_args()

    need = {"cisbp": ["cisbp"], "lambert": ["lambert"],
            "intersection": ["cisbp", "lambert"]}[args.census_type]
    missing = [f"--{x}" for x in need if getattr(args, x) is None]
    if missing:
        p.error(f"--census-type {args.census_type} requires: {', '.join(missing)}")

    meta_ovr = {"chrom": args.chrom_col, "start": args.start_col, "end": args.end_col,
                "tf": args.tf_col, "peakcount": args.peakcount_col}

    os.makedirs(args.outdir, exist_ok=True)
    tf_symbols = resolve_tf_set(args)
    per_tf, dropped = stream_sites(args.metaclusters, tf_symbols, meta_ovr)
    write_outputs(per_tf, dropped, args.outdir, args.min_sites, args.top_n)


if __name__ == "__main__":
    main()