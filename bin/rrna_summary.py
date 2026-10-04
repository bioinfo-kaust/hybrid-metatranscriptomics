#!/usr/bin/env python3
"""One-line per-sample rRNA depletion summary parsed from a bbduk log.

bbduk counts READS, so on paired input every pair is two of them. --paired halves the counts and
reports read PAIRS, which is the unit every other per-sample table in this pipeline uses — mixing
the two would make a stacked composition chart silently wrong.

Usage: rrna_summary.py SAMPLE bbduk.log CONDITION [--paired]
"""
import re
import sys

args = [a for a in sys.argv[1:] if not a.startswith("--")]
paired = "--paired" in sys.argv
sample, log, cond = args[0], args[1], args[2]
txt = open(log, errors="replace").read()


def grab(tag):
    m = re.search(tag + r":\s+([\d,]+) reads", txt)
    n = int(m.group(1).replace(",", "")) if m else 0
    return n // 2 if paired else n


total = grab("Input")
removed = grab("Total Removed")
pct = 100.0 * removed / total if total else 0.0
# the unit travels with the numbers so nothing downstream has to assume it
unit = "pairs" if paired else "sequences"
print("sample\tcondition\treads_in\trrna_removed\treads_kept\trrna_pct\tunit")
print(f"{sample}\t{cond}\t{total}\t{removed}\t{total - removed}\t{pct:.2f}\t{unit}")
