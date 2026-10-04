#!/usr/bin/env python3
"""Turn a PacBio per-barcode FL count matrix into a membership table for the pipeline.

Iso-Seq runs are usually barcoded and pooled before polishing, so the HQ transcript FASTA has no
per-sample identity — but the run reports a matrix of full-length read counts per transcript per
biosample. This converts that matrix into the `sequence_id<TAB>sample` table that
`--longread_membership` expects, so a pooled long-read set still yields per-sample statistics.

A transcript supported by reads from several samples appears once per sample: it was genuinely
detected in each of them. This lives outside the pipeline because the matrix is a platform
artifact, not something the workflow produces.

Usage:
  pacbio_fl_counts_to_membership.py --counts <*.fl_counts.csv> --out membership.tsv \\
      [--names Control,UVB,UVR] [--min-reads 1]
"""
import argparse
import csv

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True, help="CSV: name,<sample columns...>")
ap.add_argument("--out", required=True)
ap.add_argument("--names", default=None,
                help="comma-separated sample names replacing the matrix's column headers, in order")
ap.add_argument("--min-reads", type=int, default=1, help="FL reads needed to call a transcript present")
ap.add_argument("--prefix", default="", help="prepended to each sequence id (e.g. 'transcript/')")
a = ap.parse_args()

with open(a.counts, newline="") as fh:
    rd = csv.reader(fh)
    header = next(rd)
    cols = header[1:]
    if a.names:
        names = [n.strip() for n in a.names.split(",")]
        if len(names) != len(cols):
            raise SystemExit(f"--names has {len(names)} entries but the matrix has {len(cols)} sample columns: {cols}")
        cols = names
    n_rows = n_pairs = 0
    per = {c: 0 for c in cols}
    with open(a.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        for row in rd:
            if not row:
                continue
            n_rows += 1
            seq = a.prefix + row[0]
            for name, v in zip(cols, row[1:]):
                try:
                    if int(v) >= a.min_reads:
                        w.writerow([seq, name])
                        per[name] += 1
                        n_pairs += 1
                except ValueError:
                    pass

print(f"[membership] {n_rows} sequences -> {n_pairs} sequence-sample pairs")
for c in cols:
    print(f"    {c:<20} {per[c]:>8,} sequences detected")
print(f"    (pairs beyond {n_rows} are sequences detected in more than one sample)")
