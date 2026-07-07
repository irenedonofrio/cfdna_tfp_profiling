#!/usr/bin/env python3
"""samples.GC.yaml -- the handoff between the two stages.

This is the ONE piece of shared logic worth centralizing: GC correction
(stage 1) WRITES this file, nucleosome profiling (stage 2) READS it. Keeping a
single source of truth for its format means the two separable halves can never
silently disagree about it.

The write format is byte-compatible with the original Griffin snakefile output:

    samples:
      <name>:
        bam: <bam path>
        GC_bias: <absolute path to out_dir/GC_bias/<name>.GC_bias.txt>

Use as a library (write_manifest / read_manifest) or from the command line:

    python manifest.py write --samples samples.tsv --out_dir results
    python manifest.py read  --manifest results/samples.GC.yaml
"""
import os


def gc_bias_path(out_dir, sample_name):
    """Absolute path to a sample's GC_bias.txt, as recorded in the manifest."""
    return os.path.abspath(os.path.join(out_dir, 'GC_bias', sample_name + '.GC_bias.txt'))


def mappability_bias_path(out_dir, sample_name):
    """Absolute path to a sample's mappability_bias.txt (only if mappability is on)."""
    return os.path.abspath(os.path.join(out_dir, 'mappability_bias',
                                        sample_name + '.mappability_bias.txt'))


def write_manifest(path, samples, out_dir, mappability=False):
    """Write samples.GC.yaml.

    path       : output file path
    samples    : iterable of (sample_name, bam_path) pairs (order is preserved)
    out_dir    : the GC-correction output directory (where GC_bias/ lives)
    mappability: if True, also record the mappability_bias path (off by default,
                 matching the standard Griffin GC-only setup)

    Output is identical to what the original Griffin snakefile produced.
    """
    with open(path, 'w') as f:
        f.write('samples:\n')
        for name, bam in samples:
            f.write('  ' + name + ':\n')
            f.write('    bam: ' + bam + '\n')
            f.write('    GC_bias: ' + gc_bias_path(out_dir, name) + '\n')
            if mappability:
                f.write('    mappability_bias: ' + mappability_bias_path(out_dir, name) + '\n')
    return path


def read_manifest(path):
    """Parse samples.GC.yaml -> {name: {'bam': ..., 'GC_bias': ..., ...}}."""
    try:
        import yaml
    except ImportError as e:  # PyYAML ships with essentially every bioinfo env
        raise ImportError("read_manifest needs PyYAML (conda/pip install pyyaml)") from e
    with open(path) as f:
        data = yaml.safe_load(f)
    if not data or 'samples' not in data:
        raise ValueError(f"{path}: not a valid samples manifest (missing 'samples:')")
    return data['samples']


def _read_samples_tsv(tsv_path):
    """Read a <name>\\t<bam> TSV (blank lines and #comments ignored)."""
    samples = []
    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2 or not parts[0].strip() or not parts[1].strip():
                raise ValueError(f"bad samples line (need name<TAB>bam): {line!r}")
            samples.append((parts[0].strip(), parts[1].strip()))
    return samples


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description="read/write the samples.GC.yaml manifest")
    sub = ap.add_subparsers(dest='cmd', required=True)

    w = sub.add_parser('write', help='write a manifest from a samples.tsv')
    w.add_argument('--samples', required=True, help='TSV: <name>\\t<bam> per line')
    w.add_argument('--out_dir', required=True, help='GC-correction output dir')
    w.add_argument('--manifest', default=None,
                   help='output path (default: <out_dir>/samples.GC.yaml)')
    w.add_argument('--mappability', action='store_true',
                   help='also record mappability_bias paths')

    r = sub.add_parser('read', help='parse and print a manifest')
    r.add_argument('--manifest', required=True)

    args = ap.parse_args()

    if args.cmd == 'write':
        out = args.manifest or os.path.join(args.out_dir, 'samples.GC.yaml')
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        write_manifest(out, _read_samples_tsv(args.samples), args.out_dir,
                       mappability=args.mappability)
        print('wrote', out)
    elif args.cmd == 'read':
        for name, info in read_manifest(args.manifest).items():
            print(name, info)
