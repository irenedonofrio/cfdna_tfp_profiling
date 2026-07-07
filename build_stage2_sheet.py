#!/usr/bin/env python
# build_stage2_sheet.py — turn the Stage-1 GC-correction manifest into the
# Stage-2 profiling sample sheet, so Stage 2 reads Stage 1's ACTUAL outputs
# instead of re-deriving filenames.
#
# Stage-1 manifest columns : sample_id  bam  gc_counts  gc_bias
# Stage-2 sheet columns     : sample_id  bam  gc_bias    fcc_dir
#   - gc_bias : taken from the manifest (the file Stage 1 actually wrote)
#   - fcc_dir : derived as <fcc_root>/<sample_id>/fcc  (Stage-2 output location,
#               keyed on the BARE sample_id, matching run_profiling.sh)
import argparse, csv, os, sys

ap = argparse.ArgumentParser()
ap.add_argument('--gc_manifest', required=True, help='Stage-1 gc_correction_manifest.tsv')
ap.add_argument('--fcc_root', required=True, help='fcc_dir = <fcc_root>/<sample_id>/fcc')
ap.add_argument('--out', required=True, help='Stage-2 sample sheet to write')
a = ap.parse_args()

rows, missing = [], []
with open(a.gc_manifest) as f:
    rd = csv.DictReader(f, delimiter='\t')
    need = {'sample_id', 'bam', 'gc_bias'}
    if not need.issubset(rd.fieldnames or []):
        sys.exit(f"ERROR: manifest missing columns {need - set(rd.fieldnames or [])}: {a.gc_manifest}")
    for row in rd:
        sid, bam, gc_bias = row['sample_id'], row['bam'], row['gc_bias']
        if not sid:
            continue
        if not os.path.isfile(gc_bias):        # seam integrity: fail loud if Stage 1 didn't produce it
            missing.append((sid, gc_bias)); continue
        fcc_dir = os.path.join(a.fcc_root.rstrip('/'), sid, 'fcc')
        rows.append((sid, bam, gc_bias, fcc_dir))

if missing:
    sys.stderr.write("ERROR: gc_bias file(s) missing — Stage 1 did not complete for:\n")
    for sid, p in missing:
        sys.stderr.write(f"  {sid}: {p}\n")
    sys.exit(1)

out_dir = os.path.dirname(os.path.abspath(a.out))
os.makedirs(out_dir, exist_ok=True)
if os.path.exists(a.out):                       # keep a backup of any existing sheet
    os.replace(a.out, a.out + '.bak')
with open(a.out, 'w') as f:
    f.write('sample_id\tbam\tgc_bias\tfcc_dir\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')

print(f"wrote {len(rows)} sample(s) -> {a.out}")
if os.path.exists(a.out + '.bak'):
    print(f"(previous sheet backed up -> {a.out}.bak)")
