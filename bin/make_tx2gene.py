#!/usr/bin/env python3
"""Derive a transcript -> gene map by stripping the isoform suffix from transcript ids.

Rules are given as 'prefix:regex' pairs separated by ';', applied in order; the first rule whose
prefix matches the transcript id wins, and a rule with an empty regex means "no grouping"
(gene == transcript). A rule with an empty prefix is the fallback for everything else.

    --rules 'TR_:_i[0-9]+$;PBR_:;:\\.[0-9]+$'

This is deliberately weak: it groups isoforms WITHIN one assembly source and cannot reconcile
gene identity ACROSS sources, so the gene count it produces is an upper bound on distinct loci.
See locus_tx2gene.py for the coordinate-based alternative.
"""
import argparse
import re

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--rules", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

rules = []
for chunk in args.rules.split(";"):
    if not chunk:
        continue
    prefix, _, pattern = chunk.partition(":")
    rules.append((prefix, re.compile(pattern) if pattern else None))

n_tx = 0
genes = set()
with open(args.out, "w") as out:
    for line in open(args.fasta):
        if not line.startswith(">"):
            continue
        tx = line[1:].split()[0].split("|")[0]
        gene = tx
        for prefix, pattern in rules:
            if prefix and not tx.startswith(prefix):
                continue
            gene = pattern.sub("", tx) if pattern else tx
            break
        out.write(f"{tx}\t{gene}\n")
        n_tx += 1
        genes.add(gene)
print(f"[tx2gene] {n_tx} transcripts -> {len(genes)} genes ({n_tx / max(1, len(genes)):.2f} tx/gene)")
