#!/usr/bin/env python3
"""Classify sequences as host vs contaminant from two minimap2 PAF files.

Per query, score = max residue matches (PAF column 10) over all alignments to that target set.
Decision (default margin 1.1):
    host == contam == 0    -> unmapped
    host  >= contam*margin -> host
    contam >= host*margin  -> contaminant
    otherwise              -> ambiguous

Only queries that appear in at least one PAF are enumerated: sequences aligning to neither
reference are absent from every id list, and the caller must take the complement to keep them.

Usage: classify_host_contaminant.py host.paf contam.paf out_prefix [margin]
Writes: <prefix>.classification.tsv and <prefix>.{host,contaminant,ambiguous,unmapped}.ids
"""
import sys
from collections import defaultdict


def best_matches(paf):
    best, qlen = defaultdict(int), {}
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, ql, nmatch = f[0], int(f[1]), int(f[9])
            qlen[q] = ql
            if nmatch > best[q]:
                best[q] = nmatch
    return best, qlen


def main():
    host_paf, contam_paf, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
    margin = float(sys.argv[4]) if len(sys.argv) > 4 else 1.1

    hb, hql = best_matches(host_paf)
    cb, cql = best_matches(contam_paf)
    queries = set(hb) | set(cb) | set(hql) | set(cql)

    classes = ("host", "contaminant", "ambiguous", "unmapped")
    counts = defaultdict(int)
    fout = {k: open(f"{prefix}.{k}.ids", "w") for k in classes}
    with open(f"{prefix}.classification.tsv", "w") as tab:
        tab.write("query\tqlen\thost_matches\tcontaminant_matches\tclass\n")
        for q in sorted(queries):
            h, c = hb.get(q, 0), cb.get(q, 0)
            ql = hql.get(q, cql.get(q, 0))
            if h == 0 and c == 0:
                cls = "unmapped"
            elif h >= c * margin:
                cls = "host"
            elif c >= h * margin:
                cls = "contaminant"
            else:
                cls = "ambiguous"
            counts[cls] += 1
            tab.write(f"{q}\t{ql}\t{h}\t{c}\t{cls}\n")
            fout[cls].write(q + "\n")
    for f in fout.values():
        f.close()
    sys.stderr.write("Classification counts: "
                     + ", ".join(f"{k}={counts[k]}" for k in classes) + "\n")


if __name__ == "__main__":
    main()
