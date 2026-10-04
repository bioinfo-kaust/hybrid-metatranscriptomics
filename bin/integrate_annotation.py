#!/usr/bin/env python3
"""Integrate every line of annotation evidence into one per-transcript table.

Inputs are whatever is available: one or more DIAMOND result files (diamond/<db>.tsv), an
eggNOG-mapper annotations file, and an InterProScan TSV. The first DIAMOND database in
alphabetical order — or the one named by --primary-db — supplies the generic protein_* columns
that downstream steps key on; every other database contributes <name>_hit / <name>_desc.

Column 3 (has_ORF) is yes/no for EVERY row, so counting non-empty values there returns the whole
table; the ORF count is coding_status == "coding".
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--pep", required=True)
ap.add_argument("--diamond-dir", required=True)
ap.add_argument("--primary-db", default=None)
ap.add_argument("--emapper", default=None)
ap.add_argument("--interpro", default=None)
ap.add_argument("--interest-regex", default=None)
ap.add_argument("--interest-label", default="gene_family_of_interest")
ap.add_argument("--out", required=True)
ap.add_argument("--summary", default=None)
args = ap.parse_args()


def strip_orf(pid):          # PB.1.1.p1 -> PB.1.1
    return re.sub(r"\.p\d+$", "", pid)


tx2gene = {}
for line in open(args.tx2gene):
    parts = line.rstrip("\n").split("\t")
    if len(parts) >= 2:
        tx2gene[parts[0]] = parts[1]
transcripts = list(tx2gene)

has_orf = set()
for line in open(args.pep):
    if line.startswith(">"):
        has_orf.add(strip_orf(line[1:].split()[0]))


def top_hits(path):
    """qseqid -> (sseqid, pident, evalue, stitle), keeping the best bitscore.

    outfmt: qseqid sseqid pident length evalue bitscore stitle
    """
    best = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        q = strip_orf(f[0])
        bit = float(f[5])
        if q not in best or bit > best[q][-1]:
            title = f[6] if len(f) > 6 else f[1]
            best[q] = (f[1], f[2], f[4], title, bit)
    return best


dbs = {}
for path in sorted(glob.glob(os.path.join(args.diamond_dir, "*.tsv"))):
    name = os.path.basename(path)
    name = re.sub(r"^diamond_", "", name)
    name = re.sub(r"\.tsv$", "", name)
    dbs[name] = top_hits(path)

if args.primary_db and args.primary_db not in dbs:
    print(f"[WARN] --primary-db '{args.primary_db}' is not among the DIAMOND results "
          f"({', '.join(sorted(dbs)) or 'none'}); falling back to alphabetical order",
          file=sys.stderr)
primary = args.primary_db if args.primary_db in dbs else (sorted(dbs)[0] if dbs else None)
secondary = [d for d in sorted(dbs) if d != primary]

egg = {}
if args.emapper:
    hdr = None
    for line in open(args.emapper):
        if line.startswith("##"):
            continue
        if line.startswith("#query"):
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if not f or not f[0] or hdr is None:
            continue
        egg[strip_orf(f[0])] = dict(zip(hdr, f))

interpro = defaultdict(set)
interpro_go = defaultdict(set)
if args.interpro:
    for line in open(args.interpro):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 12 and f[11].startswith("IPR"):
            interpro[strip_orf(f[0])].add(f[11] + ":" + (f[12] if len(f) > 12 else ""))
        # InterProScan runs with -goterms, so column 14 carries the GO annotations it inferred
        # from the matched signatures. Those are a second, independent GO source next to
        # eggNOG's — the enrichment uses the union, so keep them in their own column.
        if len(f) >= 14:
            for g in f[13].split("|"):
                g = re.sub(r"\(.*\)$", "", g.strip())   # newer releases append "(InterPro)"
                if g.startswith("GO:"):
                    interpro_go[strip_orf(f[0])].add(g)

INTEREST = re.compile(args.interest_regex, re.I) if args.interest_regex else None


def egg_field(d, *keys):
    for k in keys:
        if d and d.get(k) and d[k] != "-":
            return d[k]
    return ""


header = ["transcript", "gene", "has_ORF", "coding_status",
          "protein_hit", "protein_pident", "protein_evalue", "protein_desc"]
for name in secondary:
    header += [f"{name}_hit", f"{name}_desc"]
header += ["eggnog_desc", "COG_cat", "GOs", "InterPro_GOs", "KEGG_ko", "KEGG_Pathway",
           "InterPro", args.interest_label]

n = orf_n = prot_n = go_n = kegg_n = ipr_n = interest_n = nc_n = 0
ipr_go_n = go_any_n = path_n = 0
sec_n = defaultdict(int)

with open(args.out, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    for t in transcripts:
        orf = t in has_orf
        pr = dbs.get(primary, {}).get(t) if primary else None
        ed = egg.get(t)
        pr_desc = pr[3] if pr else ""
        egg_desc = egg_field(ed, "Description", "Preferred_name")
        gos = egg_field(ed, "GOs")
        kegg = egg_field(ed, "KEGG_ko")
        path = egg_field(ed, "KEGG_Pathway")
        cog = egg_field(ed, "COG_category")
        ipr = ";".join(sorted(interpro.get(t, []))) if args.interpro else ""
        ipr_go = ",".join(sorted(interpro_go.get(t, []))) if args.interpro else ""
        has_func = bool(pr or ed or ipr)
        coding = "coding" if orf else ("noncoding_candidate" if not has_func else "noncoding_with_homology")
        flag = "yes" if (INTEREST and INTEREST.search(" ".join([pr_desc, egg_desc, ipr]))) else ""

        row = [t, tx2gene[t], "yes" if orf else "no", coding,
               pr[0] if pr else "", pr[1] if pr else "", pr[2] if pr else "", pr_desc]
        for name in secondary:
            h = dbs[name].get(t)
            row += [h[0] if h else "", h[3] if h else ""]
            if h:
                sec_n[name] += 1
        row += [egg_desc, cog, gos, ipr_go, kegg, path, ipr, flag]
        w.writerow(row)

        n += 1
        orf_n += bool(orf)
        prot_n += bool(pr)
        go_n += bool(gos)
        ipr_go_n += bool(ipr_go)
        go_any_n += bool(gos or ipr_go)
        kegg_n += bool(kegg)
        path_n += bool(path)
        ipr_n += bool(ipr)
        interest_n += bool(flag)
        nc_n += coding == "noncoding_candidate"

genes = len(set(tx2gene.values()))
def line(label, value):
    return f"{label:<32}{value:>10}"


lines = [
    line("transcripts", n),
    line("genes", genes),
    line("with ORF (coding)", orf_n),
    line(f"protein hit ({primary})", prot_n),
]
for name in secondary:
    lines.append(line(f"hit in {name}", sec_n[name]))
lines += [
    line("InterPro domains", ipr_n),
    line("GO terms (eggNOG)", go_n),
    line("GO terms (InterPro)", ipr_go_n),
    line("GO terms (either source)", go_any_n),
    line("KEGG KO", kegg_n),
    line("KEGG pathways", path_n),
    line(args.interest_label, interest_n),
    line("noncoding candidates", nc_n),
]
text = "\n".join(lines) + "\n"
if args.summary:
    open(args.summary, "w").write(text)
print(text, end="")
print(f"[annotation] wrote {args.out}")
