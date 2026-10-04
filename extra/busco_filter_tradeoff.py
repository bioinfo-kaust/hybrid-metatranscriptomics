#!/usr/bin/env python3
"""What the expression filter costs in BUSCO terms, at several thresholds.

Run from the project root (the directory holding results_nf/). Re-running real BUSCO per candidate
threshold would cost hours; this reads the BUSCO run that already exists for the unfiltered
reference and asks which BUSCOs survive each filter. Validated against the real filtered run:
predicted C:70.7% S:47.9% D:22.8% against measured C:70.7% [S:48.0% D:22.8%].

BUSCO renames input headers to t1,t2,... in FASTA order, so its full_table can be mapped back to
transcript ids by position. A BUSCO survives a filter when at least one of the transcripts it
matched survives; single vs duplicated is then recounted from how many survive.
"""
import csv, glob, os, sys
from collections import defaultdict

R = "results_nf"
fasta = f"{R}/transcriptome/all_sources_dedup.fasta"
table = glob.glob(f"{R}/transcriptome/qc/busco_all_sources_dedup/run_*/full_table.tsv")[0]
ann   = f"{R}/annotation/annotation.all_sources_dedup.tsv"

# t<N> -> transcript id, by position in the FASTA
order = []
with open(fasta) as fh:
    for line in fh:
        if line.startswith(">"):
            order.append(line[1:].split()[0])
tname = {f"t{i+1}": t for i, t in enumerate(order)}
print(f"reference: {len(order):,} transcripts")

# BUSCO id -> status, matched transcripts
status, hits = {}, defaultdict(set)
with open(table) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        bid, st = f[0], f[1]
        status[bid] = st
        if len(f) > 2 and f[2]:
            seq = f[2].split(":")[0]
            if seq in tname:
                hits[bid].add(tname[seq])
total = len(status)
print(f"BUSCOs in lineage: {total:,}")

# per-transcript: in how many libraries does it reach a TPM, and does it have ORF+homology
quants = sorted(glob.glob("nextflow_work/*/*/*/quant.sf"))
quants = [q for q in quants if sum(1 for _ in open(q)) - 1 == len(order)]
byname = {}
for q in quants:
    byname.setdefault(os.path.basename(os.path.dirname(q)), q)
print(f"first-pass quantifications: {len(byname)}")

nsamp = defaultdict(int); nsamp05 = defaultdict(int)
for s, q in byname.items():
    with open(q) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            tpm = float(r["TPM"])
            if tpm >= 1.0: nsamp[r["Name"]] += 1
            if tpm >= 0.5: nsamp05[r["Name"]] += 1

coding_hom = set()
with open(ann) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        if r.get("coding_status") == "coding" and (r.get("protein_hit") or r.get("uniprot_sprot_hit") or r.get("InterPro")):
            coding_hom.add(r["transcript"])

def busco_of(keep):
    s = d = f = m = 0
    for bid, st in status.items():
        surv = hits[bid] & keep
        if not surv:
            m += 1
        elif st == "Fragmented":
            f += 1
        elif len(surv) > 1:
            d += 1
        else:
            s += 1
    n = total
    return (100*(s+d)/n, 100*s/n, 100*d/n, 100*f/n, 100*m/n, len(keep))

scenarios = [
    ("unfiltered (all)",                 set(order)),
    ("TPM>=1 in >=3  (current)",          {t for t, n in nsamp.items() if n >= 3}),
    ("TPM>=1 in >=2",                     {t for t, n in nsamp.items() if n >= 2}),
    ("TPM>=1 in >=1",                     {t for t, n in nsamp.items() if n >= 1}),
    ("TPM>=0.5 in >=3",                   {t for t, n in nsamp05.items() if n >= 3}),
    ("TPM>=1 in >=3  OR ORF+homology",    {t for t, n in nsamp.items() if n >= 3} | coding_hom),
    ("TPM>=1 in >=2  OR ORF+homology",    {t for t, n in nsamp.items() if n >= 2} | coding_hom),
]
print(f"\n{'scenario':<34}{'transcripts':>12}{'C':>8}{'S':>7}{'D':>7}{'F':>7}{'M':>7}")
for name, keep in scenarios:
    c, s_, d_, f_, m_, n = busco_of(keep)
    print(f"{name:<34}{n:>12,}{c:>7.1f}%{s_:>6.1f}%{d_:>6.1f}%{f_:>6.1f}%{m_:>6.1f}%")
