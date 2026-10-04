#!/usr/bin/env python3
"""One human-readable summary of the run, assembled from the tables the pipeline produced.

Each section prints the numbers and the file they came from, so every figure quoted downstream
can be re-derived from a single place.
"""
import argparse
import csv
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--primary-map", default=None,
                help="which gene map the headline numbers come from: string | locus | corset")
args = ap.parse_args()

L = []
ACCOUNTING_FAILED = []


def section(title):
    L.append("")
    L.append(title)
    L.append("-" * len(title))


def read_tsv(name):
    path = os.path.join(args.dir, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]


def read_csv(name):
    path = os.path.join(args.dir, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return list(csv.reader(fh))


L.append("Run summary")
L.append("=" * 64)

t = read_tsv("rrna_per_sample.tsv")
if t and len(t) > 1:
    section("rRNA depletion")
    pct = [float(r[5]) for r in t[1:] if len(r) > 5]
    kept = sum(int(r[4]) for r in t[1:] if len(r) > 4)
    L.append(f"samples {len(t) - 1}  mean rRNA {sum(pct) / len(pct):.2f}%  "
             f"(min {min(pct):.2f}, max {max(pct):.2f})  read pairs kept {kept:,}")
    L.append("  source: rrna_per_sample.tsv")

t = read_tsv("decontam_per_sample.tsv")
if t and len(t) > 1:
    section("Contaminant removal (read level, before assembly)")
    rm = sum(int(r[3]) for r in t[1:])
    kept = sum(int(r[4]) for r in t[1:])
    pct = [float(r[5]) for r in t[1:]]
    L.append(f"samples {len(t) - 1}  removed {rm:,} pairs  kept {kept:,} pairs")
    L.append(f"  mean {sum(pct) / len(pct):.2f}% per sample, pooled {100 * rm / (rm + kept):.2f}%")
    L.append("  source: decontam_per_sample.tsv")

t = read_tsv("longread_decontam.tsv")
if t:
    section("Contaminant removal (long reads)")
    for r in t:
        L.append("  " + "\t".join(r))
    L.append("  source: longread_decontam.tsv")

t = read_tsv("route_contribution.tsv")
if t and len(t) > 1:
    section("Reference ladder")
    widths = [max(len(r[i]) for r in t) for i in range(len(t[0]))]
    for r in t:
        L.append("  " + "  ".join(c.ljust(widths[i]) for i, c in enumerate(r)))
    L.append("  source: route_contribution.tsv")
    L.append("  NOTE: gene counts from string-stripped transcript ids are an UPPER BOUND on")
    L.append("        distinct loci — see locus_rebuild_report.txt for the coordinate-based map.")

import glob as _glob
hits = sorted(_glob.glob(os.path.join(args.dir, "*summary*.txt")))
hits = [h for h in hits if "annotation" in os.path.basename(h)]
if hits:
    section("Annotation")
    L += ["  " + l for l in open(hits[0]).read().strip().split("\n")]
    L.append(f"  source: {os.path.basename(hits[0])}")

t = read_tsv("mapping_rates.tsv")
if t and len(t) > 1:
    section("Quantification")
    pct = [float(r[4]) for r in t[1:]]
    frags = sum(int(r[2]) for r in t[1:])
    libs = {r[1] for r in t[1:]}
    L.append(f"samples {len(t) - 1}  mean mapping {sum(pct) / len(pct):.2f}% "
             f"(min {min(pct):.2f}, max {max(pct):.2f})")
    L.append(f"library types detected: {', '.join(sorted(libs))}")
    L.append(f"fragments processed: {frags:,}  <- must equal the pairs kept after decontamination")
    L.append("  source: mapping_rates.tsv")
    # assert it rather than inviting the reader to check: a mismatch means reads were lost or
    # double-counted between decontamination and quantification, which no plot would reveal
    d = read_tsv("decontam_per_sample.tsv")
    if d and len(d) > 1:
        kept = sum(int(r[4]) for r in d[1:])
        if kept == frags:
            L.append(f"  CHECK ok: fragments processed == pairs kept after decontamination ({kept:,})")
        else:
            L.append(f"  CHECK FAILED: pairs kept {kept:,} != fragments processed {frags:,} "
                     f"(difference {abs(kept - frags):,})")
            ACCOUNTING_FAILED.append((kept, frags))

import glob as _dg


def map_kind(dirname):
    return "locus" if dirname.endswith("_locus") else "corset" if dirname.endswith("_corset") else "string"


def map_label(dirname):
    """Which gene map a dge_* directory holds, and whether it is the headline one."""
    label = {"locus": "locus rebuild (genomic coordinates)",
             "corset": "Corset (shared reads + expression)",
             "string": "string-stripped transcript ids"}[map_kind(dirname)]
    if args.primary_map and map_kind(dirname) == args.primary_map:
        label += "  [PRIMARY]"
    return label


# One section per gene map. They are separate answers to "what is a gene", not one result, so
# they are printed side by side and never merged.
de_dirs = sorted(d for d in _dg.glob(os.path.join(args.dir, "dge_*")) if os.path.isdir(d))
for d in de_dirs:
    rows = read_csv(os.path.join(os.path.basename(d), "DE_summary.csv"))
    if not rows or len(rows) < 2:
        continue
    section(f"Differential expression — {map_label(os.path.basename(d))}")
    for r in rows:
        L.append("  " + "\t".join(r))
    L.append(f"  source: {os.path.basename(d)}/DE_summary.csv")
    L.append("  NOTE: sig = padj cut AND fold-change cut (= up + down); sig_padj drops the")
    L.append("        fold-change cut, so it is the larger number.")
    L.append("        log2FoldChange is shrunken; the raw estimate is log2FoldChange_unshrunken.")

    frows = read_csv(os.path.join(os.path.basename(d), "filtering_summary.csv"))
    if frows and len(frows) > 1:
        for r in frows:
            L.append("  " + "\t".join(r))
        L.append(f"  source: {os.path.basename(d)}/filtering_summary.csv")
        L.append("  NOTE: independent filtering trims low-count genes at a per-contrast baseMean")
        L.append("        threshold; Cook's outliers leave the test entirely (n=5 is below the")
        L.append("        replacement threshold of 7).")

if len(de_dirs) > 1:
    section("Gene maps — read this before quoting a DE number")
    L.append("  The counts above come from DIFFERENT definitions of a gene over the SAME")
    L.append("  quantifications. They are a sensitivity analysis, not three measurements to")
    if args.primary_map:
        L.append(f"  The headline map for this run is: {args.primary_map} (--dge_primary_map).")
    else:
        L.append("  No primary map was named, so none of them is the headline result.")

import glob as _g

def dump(title, pattern, note=None):
    """Print the first file matching pattern, if the stage that writes it ran."""
    hits = sorted(_g.glob(os.path.join(args.dir, pattern)))
    if not hits:
        return
    section(title)
    L.extend("  " + l for l in open(hits[0]).read().strip().split("\n"))
    L.append(f"  source: {os.path.basename(hits[0])}")
    if note:
        L.append(f"  NOTE: {note}")


# filter_by_expression.py prints the "rescued by ORF+hom" line only when --rescue-annotated is on,
# so the report itself says which note applies.
_filt = sorted(_g.glob(os.path.join(args.dir, "*_expressed.filter.txt")))
_rescue_on = bool(_filt) and "rescued by ORF+hom" in open(_filt[0]).read()
dump("Reference: expression filter", "*_expressed.filter.txt",
     "this cuts the reference for QUANTIFICATION AND TESTING, not the assembly. The ladder above\n"
     "        describes the set it came from, and that set is what the annotation and the reported\n"
     "        BUSCO completeness are measured on -- a dropped transcript keeps its annotation row.\n"
     + ("        The ORF+homology rescue is ON: rescued transcripts are quantified but not\n"
        "        testable. DE and ORA run on genes with >=1 expressed transcript, which keep the\n"
        "        counts of their rescued isoforms."
        if _rescue_on else
        "        The ORF+homology rescue is off, so the quantified set and the DE-eligible set are\n"
        "        the same transcripts."))
dump("Long-read models: SQANTI3 structural categories", "sqanti3_*.summary.txt",
     "FSM/ISM are supported by the short-read annotation; NNC and antisense are the ones to check.")
dump("Gene identity — Corset (read sharing)", "*_corset.summary.txt",
     "a third map beside the string-stripped ids and the genomic loci; compare, do not average.")

path = os.path.join(args.dir, "locus_rebuild_report.txt")
if os.path.exists(path):
    section("Gene identity — locus rebuild")
    L += ["  " + l for l in open(path).read().split("\n")[:20]]
    L.append("  source: locus_rebuild_report.txt (full report)")

L.append("")
open(args.out, "w").write("\n".join(L) + "\n")

if ACCOUNTING_FAILED:
    kept, frags = ACCOUNTING_FAILED[0]
    sys.exit(f"[ERROR] read accounting: {kept:,} pairs kept after decontamination but "
             f"{frags:,} fragments quantified — reads were lost or double-counted between "
             f"the two stages.")
