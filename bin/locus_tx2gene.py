#!/usr/bin/env python3
"""Group transcripts into genes by genomic locus, and validate the result honestly.

A gene is a CONNECTED COMPONENT of transcripts linked by shared splice junctions (or shared
exons). Component grouping is required rather than "assign each transcript to one cluster":
a multi-exon gene spans several exon clusters, so picking one would be arbitrary.

Locus merging trades under-merging for over-merging. Two metrics pull in opposite directions and
both are reported:
  * genes per protein accession > 1.00  -> one protein split across several genes (under-merged)
  * loci carrying >1 distinct accession -> distinct genes fused into one locus (over-merged)
Transcripts with no confident placement keep their own gene id, flagged with a LOCUSU_ prefix.
"""
import argparse
import collections
import csv
import re

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True)
ap.add_argument("--pairs", required=True)
ap.add_argument("--bed", required=True)
ap.add_argument("--old-tx2gene", required=True)
ap.add_argument("--annotation", required=True)
ap.add_argument("--mincov", required=True)
ap.add_argument("--groupby", required=True)
ap.add_argument("--out-tx2gene", required=True)
ap.add_argument("--out-annotation", required=True)
ap.add_argument("--report", required=True)
ap.add_argument("--hit-column", default="protein_hit",
                help="annotation column holding the protein accession used for validation")
ap.add_argument("--desc-column", default="protein_desc")
args = ap.parse_args()

all_tx = [l[1:].split()[0] for l in open(args.fasta) if l.startswith(">")]
placed = {l.split("\t")[3] for l in open(args.bed)}

parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


n_pairs = 0
for line in open(args.pairs):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 2:
        continue
    union(parts[0], parts[1])
    n_pairs += 1
for t in all_tx:                      # every placed transcript is at least its own component
    if t in placed:
        find(t)

comp = collections.defaultdict(list)
for t in list(parent):
    comp[find(t)].append(t)

tx2locus, members = {}, {}
for i, (_, mem) in enumerate(sorted(comp.items(), key=lambda kv: -len(kv[1])), 1):
    lid = f"LOCUS_{i:06d}"
    members[lid] = mem
    for t in mem:
        tx2locus[t] = lid
unplaced = [t for t in all_tx if t not in tx2locus]
for t in unplaced:
    tx2locus[t] = f"LOCUSU_{t}"

with open(args.out_tx2gene, "w") as fo:
    for t in all_tx:
        fo.write(f"{t}\t{tx2locus[t]}\n")

# ---- remap the annotation's gene column (rows and evidence unchanged) --------------
with open(args.annotation) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    cols = rd.fieldnames
    rows = list(rd)
with open(args.out_annotation, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        r["gene"] = tx2locus.get(r["transcript"], r["gene"])
        w.writerow(r)

# ---- validation -------------------------------------------------------------------
old = {}
for line in open(args.old_tx2gene):
    a = line.split()
    if len(a) >= 2:
        old[a[0]] = a[1]

hit_col = args.hit_column if args.hit_column in (cols or []) else None
desc_col = args.desc_column if args.desc_column in (cols or []) else None


def genes_per_accession(gene_of):
    acc = collections.defaultdict(set)
    if hit_col:
        for r in rows:
            h = r.get(hit_col, "")
            if h in ("", "NA", "-"):
                continue
            acc[h].add(gene_of(r["transcript"]))
    if not acc:
        return 0, 0, 0.0
    # mean genes per accession. Summing per accession is the point: a gene shared by several
    # accessions must be counted once for each, or over-merging would lower this number too and
    # the two checks would stop trading off against each other.
    n_g = sum(len(v) for v in acc.values())
    multi = sum(1 for v in acc.values() if len(v) > 1)
    return len(acc), multi, n_g / len(acc)


def family_of_accession():
    """UniProt gene symbol per accession, from the description's GN= field.

    Two accessions with the same symbol are the same gene hit in different species, so a locus
    carrying both is an orthologue group, not a fusion. Only accessions with DIFFERENT symbols
    make a locus unambiguously wrong.
    """
    fam = {}
    if not (hit_col and desc_col):
        return fam
    for r in rows:
        h = r.get(hit_col, "")
        if h in ("", "NA", "-"):
            continue
        d = r.get(desc_col, "") or ""
        m = re.search(r"\bGN=(\S+)", d)
        if m:
            fam[h] = m.group(1).lower()
        else:                                   # no symbol: fall back to the protein name
            name = re.split(r"\s+(?:OS|OX|GN|PE|SV)=", re.sub(r"^\S+\s+", "", d))[0]
            fam[h] = re.sub(r"\W+", " ", name.lower()).strip() or h
    return fam


def acc_per_gene(gene_of):
    d = collections.defaultdict(set)
    if hit_col:
        for r in rows:
            h = r.get(hit_col, "")
            if h not in ("", "NA", "-"):
                d[gene_of(r["transcript"])].add(h)
    return d


o_acc, o_multi, o_ratio = genes_per_accession(lambda t: old.get(t, t))
n_acc, n_multi, n_ratio = genes_per_accession(lambda t: tx2locus[t])
n_old, n_new = len(set(old.values())), len(set(tx2locus.values()))

L = []
L.append("Locus-based gene rebuild — validation report")
L.append("=" * 64)
L.append(f"transcripts                        : {len(all_tx):,}")
L.append(f"placed on the genome (cov >= {args.mincov}) : {len(all_tx) - len(unplaced):,}")
L.append(f"unplaced (kept as own gene)        : {len(unplaced):,}")
L.append("")
L.append(f"gene count BEFORE (string-stripped): {n_old:,}")
L.append(f"gene count AFTER  (genomic loci)   : {n_new:,}")
if n_old:
    L.append(f"  reduction                        : {n_old - n_new:,}  ({100 * (n_old - n_new) / n_old:.1f}%)")
L.append("")
if hit_col:
    L.append(f"Redundancy check — genes sharing one {hit_col} accession")
    L.append(f"  BEFORE: {o_acc:,} accessions, {o_multi:,} hit by >1 gene "
             f"({100 * o_multi / max(1, o_acc):.1f}%), {o_ratio:.2f} genes/accession")
    L.append(f"  AFTER : {n_acc:,} accessions, {n_multi:,} hit by >1 gene "
             f"({100 * n_multi / max(1, n_acc):.1f}%), {n_ratio:.2f} genes/accession")
    L.append("  (under-merging pushes this ABOVE 1.00; BELOW 1.00 means loci carry several")
    L.append("   different proteins, i.e. distinct genes have been fused — see the next check)")
    L.append("")
    o_apg, n_apg = acc_per_gene(lambda t: old.get(t, t)), acc_per_gene(lambda t: tx2locus[t])
    o_bad = sum(1 for v in o_apg.values() if len(v) > 1)
    n_bad = sum(1 for v in n_apg.values() if len(v) > 1)
    L.append("OVER-MERGE check — genes carrying >1 DISTINCT accession")
    L.append(f"  BEFORE: {o_bad:,} / {len(o_apg):,} annotated genes ({100 * o_bad / max(1, len(o_apg)):.1f}%)")
    L.append(f"  AFTER : {n_bad:,} / {len(n_apg):,} annotated loci  ({100 * n_bad / max(1, len(n_apg)):.1f}%)")
    L.append(f"  worst locus carries {max((len(v) for v in n_apg.values()), default=0)} distinct accessions")
    fam = family_of_accession()
    if fam:
        def true_over(apg):
            bad = [v for v in apg.values() if len(v) > 1]
            return sum(1 for v in bad if len({fam.get(a, a) for a in v}) > 1)
        o_true, n_true = true_over(o_apg), true_over(n_apg)
        L.append("  of those, carrying accessions with DIFFERENT gene symbols (a real fusion;")
        L.append("  the rest are one gene matched in several species, which is not a wrong merge)")
        L.append(f"  BEFORE: {o_true:,} ({100 * o_true / max(1, len(o_apg)):.1f}% of annotated genes)")
        L.append(f"  AFTER : {n_true:,} ({100 * n_true / max(1, len(n_apg)):.1f}% of annotated loci)")
    if n_bad / max(1, len(n_apg)) > 1.5 * (o_bad / max(1, len(o_apg))):
        L.append("  >> WARNING: over-merging rose sharply. Raise --locus_mincov, or keep the")
        L.append("     string-stripped map for anything that depends on gene identity.")
    L.append("")

L.append(f"grouping rule                      : {args.groupby}")
L.append(f"linking pairs                      : {n_pairs:,}")
L.append(f"loci (connected components)        : {len(members):,}")
L.append("")
L.append("LARGEST LOCI — inspect these for over-merging before adopting the map")
L.append("-" * 64)
desc = {}
if desc_col:
    for r in rows:
        d = r.get(desc_col, "")
        if d not in ("", "NA", "-"):
            desc.setdefault(tx2locus[r["transcript"]], collections.Counter())[d.split(" OS=")[0]] += 1
for lid, mem in sorted(members.items(), key=lambda kv: -len(kv[1]))[:20]:
    top = desc.get(lid, collections.Counter()).most_common(2)
    tops = "; ".join(f"{d[:52]} x{c}" for d, c in top) or "(no protein hit)"
    L.append(f"{lid}  {len(mem):4d} tx")
    L.append(f"    {tops}")
L.append("")
L.append("A locus whose top descriptions disagree is a likely over-merge.")
L.append("If many do, re-run with a higher --locus_mincov, or keep the string-stripped map.")

open(args.report, "w").write("\n".join(L) + "\n")
print("\n".join(L[:24]))
print(f"\n[written] {args.out_tx2gene}\n[written] {args.out_annotation}\n[written] {args.report}")
