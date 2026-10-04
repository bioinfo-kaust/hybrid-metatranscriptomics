#!/usr/bin/env python3
"""Build the reference ladder table: one row per reference variant.

Reads three directories of per-variant files keyed on the variant label:
    <stats>/<label>.stats.tsv     seqkit stats -T output
    <tx2gene>/tx2gene.<label>.tsv transcript -> gene map
    <busco>/<label>.busco.txt     the BUSCO completeness string (or NA)

Any variant missing a file still gets a row, with '-' in that column: the ladder is meant to be
readable even when an optional rung was skipped.
"""
import argparse
import glob
import os

ap = argparse.ArgumentParser()
ap.add_argument("--stats", required=True)
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--busco", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

ORDER = [
    "isoseq",
    "isoseq_stringtie",
    "isoseq_stringtie_trinity",
    "isoseq_stringtie_lrfree",
    "isoseq_stringtie_trinity_lrfree",
    "all_sources_dedup",
]
DESCRIPTION = {
    "isoseq":                          "Iso-Seq models on the genome (Route A)",
    "isoseq_stringtie":                "+ genome-guided short-read loci (StringTie)",
    "isoseq_stringtie_trinity":        "+ de novo short-read loci (Trinity)",
    "isoseq_stringtie_lrfree":         "+ reference-free long-read loci only",
    "isoseq_stringtie_trinity_lrfree": "+ de novo short-read + reference-free long-read loci",
    "all_sources_dedup":               "cross-source deduplicated (the assembly: annotated and scored here)",
}


def read_stats(path):
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh]
    if len(rows) < 2:
        return None
    return dict(zip(rows[0], rows[1]))


labels = []
for f in sorted(glob.glob(os.path.join(args.stats, "*.stats.tsv"))):
    labels.append(os.path.basename(f)[: -len(".stats.tsv")])
labels.sort(key=lambda l: ORDER.index(l) if l in ORDER else len(ORDER))

with open(args.out, "w") as out:
    out.write("reference\tbasis\tn_transcripts\ttotal_bp\tn50\tn_genes\tBUSCO\n")
    for label in labels:
        st = read_stats(os.path.join(args.stats, f"{label}.stats.tsv")) or {}
        t2g = os.path.join(args.tx2gene, f"tx2gene.{label}.tsv")
        n_genes = "-"
        if os.path.exists(t2g):
            n_genes = str(len({l.split("\t")[1].strip() for l in open(t2g) if "\t" in l}))
        bus = os.path.join(args.busco, f"{label}.busco.txt")
        busco = open(bus).read().strip() if os.path.exists(bus) else "-"
        out.write("\t".join([
            label,
            DESCRIPTION.get(label, "-"),
            st.get("num_seqs", "-"),
            st.get("sum_len", "-"),
            st.get("N50", "-"),
            n_genes,
            busco or "-",
        ]) + "\n")
print(f"[ladder] {len(labels)} reference variants -> {args.out}")
