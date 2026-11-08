#!/usr/bin/env python3
"""One human-readable summary of the run, assembled from the tables the pipeline produced.

Each section prints the numbers and the file they came from, so every figure quoted downstream
can be re-derived from a single place.
"""
import argparse
import csv
import glob
import os
import re
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

# The expression-filtered set is a ladder rung too, and the one every reported number below
# describes. route_contribution.tsv is written before the filter runs, so it cannot carry it --
# SEQ_STATS_EXPR and BUSCO_EXPR publish for it instead. Appending it here puts the analysed set and
# the assembly it was cut from on one line each, with completeness for both, so the plain-text
# summary says the same thing the dashboard does.
if t and len(t) > 1:
    for _st in sorted(glob.glob(os.path.join(args.dir, "*_expressed.stats.tsv"))):
        _name = os.path.basename(_st)[:-len(".stats.tsv")]
        if any(r and r[0] == _name for r in t[1:]):
            continue
        _rows = read_tsv(os.path.basename(_st))
        if not _rows or len(_rows) < 2:
            continue
        _c = dict(zip(_rows[0], _rows[1]))
        _busco = "-"
        _bp = os.path.join(args.dir, f"{_name}.busco.txt")
        if os.path.exists(_bp):
            _m = re.search(r"C:[\d.]+%\[S:[\d.]+%,D:[\d.]+%\],F:[\d.]+%,M:[\d.]+%",
                           open(_bp).read())
            if _m:
                _busco = _m.group(0)
        _genes = "-"
        _fp = os.path.join(args.dir, f"{_name}.filter.txt")
        if os.path.exists(_fp):
            _g = re.search(r"genes\s+[\d,]+\s*->\s*([\d,]+)", open(_fp).read())
            if _g:
                _genes = _g.group(1).replace(",", "")
        t.append([_name, "expression-filtered (THE ANALYSED SET: annotated, quantified, tested)",
                  _c.get("num_seqs", "-"), _c.get("sum_len", "-"), "-", _genes, _busco])

if t and len(t) > 1:
    section("Reference ladder")
    widths = [max(len(r[i]) for r in t) for i in range(len(t[0]))]
    for r in t:
        L.append("  " + "  ".join(c.ljust(widths[i]) for i, c in enumerate(r)))
    L.append("  source: route_contribution.tsv")
    L.append("  NOTE: gene counts from string-stripped transcript ids are an UPPER BOUND on")
    L.append("        distinct loci — see locus_rebuild_report.txt for the coordinate-based map.")
    L.append("  NOTE: the last row is the analysed set — everything reported below it (annotation,")
    L.append("        quantification, DE, enrichment) describes that set, not the full assembly.")
    L.append("        Both carry their own BUSCO: the assembly's says how complete the")
    L.append("        reconstruction is, the analysed set's how much of the conserved gene space")
    L.append("        this experiment actually expressed.")

import glob as _glob
hits = sorted(_glob.glob(os.path.join(args.dir, "*summary*.txt")))
hits = [h for h in hits if "annotation" in os.path.basename(h)]
if hits:
    section("Annotation")
    # The annotation module is handed the string-stripped tx2gene, because the locus rebuild needs
    # the genome placement and runs after it. Every other line of its summary counts TRANSCRIPTS
    # and so is the same under any map; this one line is not, and unlabelled it reads as the
    # headline gene count. Label it and put the locus figure beside it.
    _locus_after = None
    _lr = os.path.join(args.dir, "locus_rebuild_report.txt")
    if os.path.exists(_lr):
        _m = re.search(r"gene count AFTER\s*\(genomic loci\)\s*:\s*([\d,]+)", open(_lr).read())
        if _m:
            # the module's own rows are unformatted, so match them rather than the report's commas
            _locus_after = _m.group(1).replace(",", "")
    for line in open(hits[0]).read().strip().split("\n"):
        _g = re.match(r"^(genes)(\s+)([\d,]+)\s*$", line)
        if _g and _locus_after:
            width = len(_g.group(1)) + len(_g.group(2))
            L.append("  " + "genes (string-stripped)".ljust(width) + _g.group(3))
            L.append("  " + "genes (locus rebuild)".ljust(width) + _locus_after)
        else:
            L.append("  " + line)
    if _locus_after:
        L.append("  NOTE: every other row counts TRANSCRIPTS and is therefore the same under any")
        L.append("        gene map. Only the gene row moves, and the headline map for this run")
        L.append(f"        ({args.primary_map}) is the second of the two.")
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


# The filter decides the analysed set, so this note has to be read before any count below it.
dump("Reference: the expression filter", "*_expressed.filter.txt",
     "this is the set that was ANNOTATED, quantified and tested. Everything after this section\n"
     "        describes it and nothing else. What the filter does not touch is the assembly: the\n"
     "        ladder above keeps its own transcript counts and BUSCO completeness, measured on the\n"
     "        full deduplicated set, so those two numbers describe more transcripts than any other\n"
     "        number in this report.")
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
