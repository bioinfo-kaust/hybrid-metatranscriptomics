#!/usr/bin/env python3
"""Best spliced alignment per transcript -> exon-block BED (stdout) + splice junctions (file).

The transcript's genomic SPAN includes its introns, so merging spans fuses any gene nested in
another gene's intron. Split the CIGAR on N (intron) operations instead and emit one BED
interval per exon, plus one junction record per intron.

Usage: paf_to_exons.py placed.paf MINCOV junctions.tsv > best.bed
"""
import re
import sys

paf, mincov, jout = sys.argv[1], float(sys.argv[2]), sys.argv[3]

best = {}
for line in open(paf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
    strand, target, tstart, nmatch = f[4], f[5], int(f[7]), int(f[9])
    if (qe - qs) / qlen < mincov:
        continue
    if q in best and nmatch <= best[q][0]:
        continue
    cigar = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
    best[q] = (nmatch, target, tstart, strand, cigar)

with open(jout, "w") as jf:
    for q, (_, target, tstart, strand, cigar) in best.items():
        if not cigar:                       # no CIGAR: fall back to a point interval
            print(f"{target}\t{tstart}\t{tstart + 1}\t{q}\t0\t{strand}")
            continue
        pos = blk_start = tstart
        blk_len = 0
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            n = int(n)
            if op in "M=XD":
                blk_len += n
                pos += n
            elif op == "N":                 # intron: close the exon, record the junction
                if blk_len:
                    print(f"{target}\t{blk_start}\t{blk_start + blk_len}\t{q}\t0\t{strand}")
                jf.write(f"{q}\t{target}\t{strand}\t{pos}\t{pos + n}\n")
                pos += n
                blk_start = pos
                blk_len = 0
        if blk_len:
            print(f"{target}\t{blk_start}\t{blk_start + blk_len}\t{q}\t0\t{strand}")
