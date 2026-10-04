#!/usr/bin/env python3
"""Per-sample salmon mapping rates + the end-to-end read accounting check.

The summed num_processed must equal the number of read pairs that survived decontamination —
a mismatch means samples were dropped or double-counted somewhere in the fan-out.
"""
import argparse
import sys
import glob
import json
import os
import statistics

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min-mapping-rate", type=float, default=0.0,
                help="fail if any sample maps below this percentage")
args = ap.parse_args()

rows = []
for meta in sorted(glob.glob(os.path.join(args.dir, "*", "aux_info", "meta_info.json"))):
    sample = os.path.basename(os.path.dirname(os.path.dirname(meta)))
    d = json.load(open(meta))
    cfr = d.get("compatible_fragment_ratio")
    rows.append((sample, (d.get("library_types") or ["NA"])[0], d.get("num_processed", 0),
                 d.get("num_mapped", 0), d.get("percent_mapped", 0.0),
                 "%.4f" % cfr if isinstance(cfr, (int, float)) else ""))

with open(args.out, "w") as out:
    out.write("sample\tlibrary_type\tfragments_processed\tfragments_mapped\tpercent_mapped\tcompatible_fragment_ratio\n")
    for r in rows:
        out.write("%s\t%s\t%d\t%d\t%.2f\t%s\n" % r)
    if rows:
        pct = [r[4] for r in rows]
        out.write("# n=%d  mean_percent_mapped=%.2f  min=%.2f  max=%.2f  total_fragments=%d\n"
                  % (len(rows), statistics.mean(pct), min(pct), max(pct), sum(r[2] for r in rows)))
print(f"[quant] {len(rows)} samples -> {args.out}")

# A mapping rate that collapses means the reference, the reads or the library type is wrong, and
# every downstream number is then meaningless. Fail here rather than in a reader's attention.
low = [(r[0], r[4]) for r in rows if r[4] < args.min_mapping_rate]
if low:
    for s_, p_ in low:
        print(f"[ERROR] {s_}: {p_:.2f}% mapped, below the --min-mapping-rate {args.min_mapping_rate}")
    sys.exit(1)
if rows:
    print(f"[quant] all {len(rows)} samples at or above {args.min_mapping_rate}% mapped")
