#!/usr/bin/env python3
"""Bucket this run's genes against a published reference transcriptome.

Deliberately NOT part of the pipeline: comparing to a specific paper is an interpretation step
that depends on which reference you consider authoritative, and it should not silently change
when the pipeline is re-run.

Buckets, per gene:
  matched                  at least one isoform hits the published set convincingly
  annotation_gap           hits it, but the published entry is uncharacterised while we annotate it
  novel_locus_annotated    no hit, and we have functional evidence
  novel_locus_unannotated  no hit, and no functional evidence

Three caveats that materially change the numbers — read them before quoting a novelty rate:
 1. The reference is usually a TRANSCRIPTOME, so genes the authors annotated but did not express
    are invisible and ours are over-called as novel.
 2. Coverage is measured on OUR transcript, so short de novo fragments fail the coverage test
    even when they are real. Check --report for how many rejections that accounts for.
 3. `blastn -max_target_seqs 1` does NOT return the best hit; use 25 and let bitscore decide
    (compare_to_published.sh does this).
"""
import argparse
import csv
import glob
import os
import re
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--blast", required=True, help="outfmt 6: qseqid sseqid pident length qlen evalue bitscore")
ap.add_argument("--annotation", required=True)
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--dge-dir", default=None, help="flag genes significant in any contrast")
ap.add_argument("--reference-annotation", default=None,
                help="optional TSV: <reference sequence id>\\t<description>, for the annotation-gap call")
ap.add_argument("--min-identity", type=float, default=80.0)
ap.add_argument("--min-coverage", type=float, default=0.4)
ap.add_argument("--max-evalue", type=float, default=1e-10)
ap.add_argument("--padj", type=float, default=0.05)
ap.add_argument("--lfc", type=float, default=1.0)
ap.add_argument("--out", required=True)
ap.add_argument("--report", default=None)
a = ap.parse_args()

ref_desc = {}
if a.reference_annotation:
    for line in open(a.reference_annotation):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            ref_desc[f[0]] = f[1]

ann = {}
with open(a.annotation) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ann[r["transcript"]] = r

tx2gene = {}
for line in open(a.tx2gene):
    f = line.split()
    if len(f) >= 2:
        tx2gene[f[0]] = f[1]

de_genes = set()
if a.dge_dir:
    for f in glob.glob(os.path.join(a.dge_dir, "DE_*.csv")):
        if "DE_summary" in f:
            continue
        with open(f) as fh:
            rd = csv.DictReader(fh)
            for r in rd:
                gene = r[rd.fieldnames[0]]
                try:
                    if r["padj"] not in ("NA", "") and float(r["padj"]) < a.padj \
                       and abs(float(r["log2FoldChange"])) >= a.lfc:
                        de_genes.add(gene)
                except (KeyError, ValueError):
                    pass

best = {}
n_rejected_cov = 0
for line in open(a.blast):
    f = line.rstrip("\n").split("\t")
    if len(f) < 7:
        continue
    q, s = f[0], f[1]
    pid, ln, qlen, ev, bit = float(f[2]), int(f[3]), int(f[4]), float(f[5]), float(f[6])
    cov = ln / qlen if qlen else 0.0
    if q not in best or bit > best[q][-1]:
        best[q] = (s, pid, cov, ev, bit)


def is_match(h):
    if not h:
        return False
    _, pid, cov, ev, _ = h
    return ev <= a.max_evalue and pid >= a.min_identity and cov >= a.min_coverage


for q, h in best.items():
    if h[1] >= a.min_identity and h[3] <= a.max_evalue and h[2] < a.min_coverage:
        n_rejected_cov += 1

UNCHAR = re.compile(r"uncharacter|hypothetical|unknown|predicted protein", re.I)

genes = defaultdict(lambda: dict(iso=0, matched=0, ref_ann=set(), ours=set(), has_go=False))
for t, g in tx2gene.items():
    G = genes[g]
    G["iso"] += 1
    h = best.get(t)
    if is_match(h):
        G["matched"] += 1
        if h[0] in ref_desc:
            G["ref_ann"].add(ref_desc[h[0]])
    r = ann.get(t)
    if r:
        if r.get("protein_desc"):
            G["ours"].add(r["protein_desc"])
        if r.get("GOs"):
            G["has_go"] = True

counts = defaultdict(int)
with open(a.out, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["gene", "n_isoforms", "n_matched", "bucket", "differentially_expressed",
                "our_protein_hit", "our_has_GO", "reference_annotation"])
    for g, G in genes.items():
        our_annot = bool(G["ours"]) or G["has_go"]
        if G["matched"] == 0:
            bucket = "novel_locus_annotated" if our_annot else "novel_locus_unannotated"
        elif our_annot and ref_desc and (not G["ref_ann"] or all(UNCHAR.search(x) for x in G["ref_ann"])):
            bucket = "annotation_gap"
        else:
            bucket = "matched"
        counts[bucket] += 1
        w.writerow([g, G["iso"], G["matched"], bucket,
                    "yes" if g in de_genes else "",
                    " | ".join(sorted(G["ours"]))[:120],
                    "yes" if G["has_go"] else "",
                    " | ".join(sorted(G["ref_ann"]))[:120]])

tot = sum(counts.values()) or 1
lines = [f"=== {sum(counts.values())} genes compared against the published reference ==="]
for k in ("matched", "annotation_gap", "novel_locus_annotated", "novel_locus_unannotated"):
    lines.append(f"  {k}: {counts[k]} ({100 * counts[k] / tot:.1f}%)")
nov_de = sum(1 for g, G in genes.items() if G["matched"] == 0 and g in de_genes)
lines.append(f"  novel loci that are differentially expressed: {nov_de}")
lines.append(f"  hits rejected on coverage < {a.min_coverage} alone: {n_rejected_cov}")
if not ref_desc:
    lines.append("  (no --reference-annotation given: the annotation_gap bucket is not evaluated)")
text = "\n".join(lines)
print(text)
if a.report:
    open(a.report, "w").write(text + "\n")
