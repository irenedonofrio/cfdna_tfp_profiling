#!/usr/bin/env python3
"""
config_helper.py — the ONLY YAML/sheet reader on the bash side of the pipeline.

The bash driver never parses YAML itself; it calls this helper and `eval`s the
shell-safe `KEY=value` assignments it prints. Python does all YAML/TSV reading
(here and in the 04_*.py scripts); bash reads neither.

Two modes:

  config <config.yaml> KEY [KEY ...]
      Look up scalar (or list-of-scalars) keys across the merged
      deployment + method sections and print `<UPPER_KEY>='value'`.

      PATH RESOLUTION (so the whole project folder is relocatable):
        * refs_root       : anchored to THIS config file's folder; the TFP_REFS
                            env var overrides it (for refs mounted outside the repo).
        * ref-typed keys  : reference_fa, mappable_bed, chrom_sizes, tf_dir,
                            genome_gc_freq -> joined onto refs_root when relative.
        * dir-typed keys  : sample_sheet, results_root, fcc_root -> anchored to
                            THIS config file's folder when relative.
        * everything else : emitted verbatim (threads, mapq, autosomes, tmpdir, ...).
      Absolute values ALWAYS pass through unchanged (back-compat with old configs).

      Lists are space-joined. Errors on unknown keys, cross-section name
      collisions, or nested structures.

  sample <samples.tsv> <sample_id>
      Print BAM=, GC_BIAS=, FCC_DIR= for the single matching row.
      Exits nonzero if the sample_id is missing OR appears more than once
      (sample_id must be unique — duplicates are a hard error).

Usage from the driver:
    eval "$(python config_helper.py config config.yaml threads results_root tf_dir chrom_sizes)"
    eval "$(python config_helper.py sample samples.tsv "$SAMPLE")"

Dependencies: pyyaml (already in the env). No other dependencies.
"""

import os
import sys
import csv
import shlex


def die(msg):
    print(f"config_helper: {msg}", file=sys.stderr)
    sys.exit(1)


def emit(key, value):
    """Print a shell-safe assignment `KEY='value'` (KEY uppercased)."""
    if isinstance(value, dict):
        die(f"key '{key}' is a mapping; not shell-emittable")
    if isinstance(value, (list, tuple)):
        if any(isinstance(x, (list, dict)) for x in value):
            die(f"key '{key}' is a nested structure; not shell-emittable")
        value = " ".join(str(x) for x in value)
    print(f"{key.upper()}={shlex.quote(str(value))}")


# --------------------------------------------------------------------------
# config mode
# --------------------------------------------------------------------------
# Keys whose relative values are joined onto refs_root:
_REF_RELATIVE = {"reference_fa", "mappable_bed", "chrom_sizes", "tf_dir", "genome_gc_freq"}
# Keys whose relative values are anchored to the config file's own directory:
_DIR_RELATIVE = {"sample_sheet", "results_root", "fcc_root", "refs_root"}


def cmd_config(args):
    if len(args) < 2:
        die("usage: config <config.yaml> KEY [KEY ...]")
    cfg_path, keys = args[0], args[1:]

    try:
        import yaml
    except ImportError:
        die("pyyaml not available in this environment")

    try:
        with open(cfg_path) as fh:
            cfg = yaml.safe_load(fh) or {}
    except FileNotFoundError:
        die(f"config not found: {cfg_path}")

    deployment = cfg.get("deployment", {}) or {}
    method = cfg.get("method", {}) or {}

    # Merge the two sections into one lookup namespace; a key present in both
    # is a config error (the two sections must not overlap).
    merged = {}
    for section_name, section in (("deployment", deployment), ("method", method)):
        for k, v in section.items():
            if k in merged:
                die(f"key '{k}' present in both deployment and method sections")
            merged[k] = v

    # ---- path anchoring ------------------------------------------------------
    # Relative paths are resolved so the project folder can be moved/mounted
    # anywhere. Absolute paths are returned unchanged.
    cfg_dir = os.path.dirname(os.path.abspath(cfg_path))

    def anchor(p):
        p = str(p)
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(cfg_dir, p))

    # refs_root: TFP_REFS env var wins; else the config value; anchored to cfg_dir.
    refs_root_raw = os.environ.get("TFP_REFS") or merged.get("refs_root")
    refs_root = anchor(refs_root_raw) if refs_root_raw is not None else None

    def resolve(k, v):
        if k in _REF_RELATIVE:
            if os.path.isabs(str(v)):
                return v
            if refs_root is None:
                die(f"key '{k}' is relative but no refs_root (and no TFP_REFS) is set")
            return os.path.join(refs_root, str(v))
        if k == "refs_root":
            if refs_root is None:
                die("key 'refs_root' requested but neither refs_root nor TFP_REFS is set")
            return refs_root
        if k in _DIR_RELATIVE:
            return anchor(v)
        return v  # non-path keys: verbatim (threads, mapq, autosomes, tmpdir, ...)

    for k in keys:
        # refs_root may come solely from TFP_REFS, so allow it even if absent in the file.
        if k not in merged and k != "refs_root":
            die(f"key '{k}' not found in config (looked in deployment + method)")
        emit(k, resolve(k, merged.get(k)))


# --------------------------------------------------------------------------
# sample mode
# --------------------------------------------------------------------------
def cmd_sample(args):
    if len(args) != 2:
        die("usage: sample <samples.tsv> <sample_id>")
    sheet_path, sample_id = args

    try:
        with open(sheet_path, newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except FileNotFoundError:
        die(f"sample sheet not found: {sheet_path}")

    # Drop blank lines and comment lines
    rows = [r for r in rows if r and not r[0].lstrip().startswith("#")]
    if not rows:
        die(f"sample sheet is empty: {sheet_path}")

    header = [c.strip() for c in rows[0]]
    idx = {name: i for i, name in enumerate(header)}

    required = ["sample_id", "bam", "gc_bias", "fcc_dir"]
    missing = [c for c in required if c not in idx]
    if missing:
        die(f"sample sheet missing column(s): {', '.join(missing)} "
            f"(need: {', '.join(required)})")

    si = idx["sample_id"]
    matches = [r for r in rows[1:]
               if len(r) > si and r[si].strip() == sample_id]

    if len(matches) == 0:
        die(f"sample '{sample_id}' not found in {sheet_path}")
    if len(matches) > 1:
        die(f"sample '{sample_id}' appears {len(matches)} times in {sheet_path} "
            f"(sample_id must be unique)")

    row = matches[0]

    def cell(col):
        i = idx[col]
        if i >= len(row) or not row[i].strip():
            die(f"sample '{sample_id}': empty value for column '{col}'")
        return row[i].strip()

    emit("bam", cell("bam"))
    emit("gc_bias", cell("gc_bias"))
    emit("fcc_dir", cell("fcc_dir"))


# --------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        die("usage: config_helper.py <config|sample> ...")
    mode, rest = sys.argv[1], sys.argv[2:]
    if mode == "config":
        cmd_config(rest)
    elif mode == "sample":
        cmd_sample(rest)
    else:
        die(f"unknown mode '{mode}' (use 'config' or 'sample')")


if __name__ == "__main__":
    main()
