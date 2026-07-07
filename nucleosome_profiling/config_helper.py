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
      Lists are space-joined (only used for simple token lists such as
      chromosome names). Errors on unknown keys, cross-section name
      collisions, or nested structures.

  sample <samples.tsv> <sample_id>
      Print BAM=, GC_BIAS=, FCC_DIR= for the single matching row.
      Exits nonzero if the sample_id is missing OR appears more than once
      (sample_id must be unique — duplicates are a hard error).

Usage from the driver:
    eval "$(python config_helper.py config config.yaml threads results_root tf_dir chrom_sizes)"
    eval "$(python config_helper.py sample samples.tsv "$SAMPLE")"

Dependencies: pyyaml (already in tfp_env). No other dependencies.
"""

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

    for k in keys:
        if k not in merged:
            die(f"key '{k}' not found in config (looked in deployment + method)")
        emit(k, merged[k])


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
