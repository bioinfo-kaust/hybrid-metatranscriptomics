#!/usr/bin/env python3
"""Per-sample accounting for the long-read side, from the pooled classification.

Long reads are pooled before they are mapped — one pass over a multi-gigabase contaminant
reference instead of one per sample — so the per-sample numbers are recovered afterwards from a
membership map (sequence id -> sample) rather than from separate runs. Two cases, one code path:

  * several long-read samples in the sample sheet -> one membership file per sample, written by
    the pipeline from each sample's own reads;
  * one already-pooled file (a polished Iso-Seq set, say) plus an external membership table
    (--longread_membership) -> the same accounting, with sequences attributed to the samples that
    contributed reads to them.

In the second case a sequence can belong to more than one sample. That is real, not a bug: the
per-sample columns then say "detected in", they sum to more than the pooled total, and the note
at the bottom of the table says so.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--membership", required=True,
                help="directory of, or a single, 'sequence_id<TAB>sample' TSV")
ap.add_argument("--classification", required=True,
                help="classifier output: query, qlen, host, contaminant, class")
ap.add_argument("--input-ids", default=None, help="FASTA or id list that entered decontamination")
ap.add_argument("--rrna-dir", default=None, help="directory of *.rrna_hits.fasta (removed by the rRNA screen)")
ap.add_argument("--conditions", default=None, help="'sample<TAB>condition' TSV")
ap.add_argument("--keep-ambiguous", action="store_true",
                help="count ambiguous as kept (matches --keep_ambiguous_long)")
ap.add_argument("--out", required=True)
args = ap.parse_args()


def read_ids(path):
    """Ids from a FASTA (its headers) or from a plain one-per-line list.

    The format is decided by the first character of the file rather than per line, so a sequence
    line that happens to look like an identifier can never be mistaken for one.
    """
    out = set()
    with open(path, errors="replace") as fh:
        first = fh.read(1)
        if not first:
            return out
        fh.seek(0)
        is_fasta = first == ">"
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if is_fasta:
                if line.startswith(">"):
                    out.add(line[1:].split()[0])
            else:
                out.add(line.split()[0])
    return out


# ---- membership: sequence -> the samples it belongs to -----------------------------
paths = sorted(glob.glob(os.path.join(args.membership, "*"))) if os.path.isdir(args.membership) \
    else [args.membership]
members = defaultdict(set)          # sample -> {seq ids}
samples = []
for p in paths:
    if os.path.isdir(p):
        continue
    with open(p, errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or f[0] in ("sequence_id", "id", "name"):
                continue
            members[f[1]].add(f[0])
for s in members:
    samples.append(s)
samples.sort()
if not samples:
    # nothing to attribute: leave the output absent rather than writing an empty table that
    # would show up downstream as a per-sample section with no samples in it
    sys.stderr.write("[longread] no membership rows — no per-sample table written\n")
    raise SystemExit(0)

# ---- what the rRNA screen removed --------------------------------------------------
rrna_removed = set()
if args.rrna_dir and os.path.isdir(args.rrna_dir):
    for p in sorted(glob.glob(os.path.join(args.rrna_dir, "*"))):
        rrna_removed |= read_ids(p)

# Restrict membership to sequences this run actually saw. --input-ids is the set that entered
# DECONTAMINATION, i.e. after the rRNA screen, so the rRNA-removed ids must be added back in or
# the per-sample rRNA column would read zero and the input counts would be understated.
if args.input_ids and os.path.exists(args.input_ids):
    universe = read_ids(args.input_ids) | rrna_removed
    if universe:
        for s in samples:
            members[s] &= universe

# ---- the classifier's verdict ------------------------------------------------------
verdict = {}
with open(args.classification, errors="replace") as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd, None)
    for row in rd:
        if len(row) >= 5:
            verdict[row[0]] = row[4]

conditions = {}
if args.conditions and os.path.exists(args.conditions):
    for line in open(args.conditions):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            conditions[f[0]] = f[1]

DROPPED = {"contaminant"} if args.keep_ambiguous else {"contaminant", "ambiguous"}
CLASSES = ("host", "no_alignment", "contaminant", "ambiguous")
# older classifier outputs name the removal class after the organism rather than the role
ALIAS = {"symbiont": "contaminant", "sym": "contaminant"}
unknown = set()


def klass(seq_id):
    """The classifier only enumerates sequences seen in at least one alignment; the rest aligned
    to neither reference and are kept as host candidates."""
    c = verdict.get(seq_id)
    if c is None:
        return "no_alignment"
    c = ALIAS.get(c, c)
    if c not in CLASSES:
        unknown.add(c)
        return "contaminant"          # an unrecognised class is never silently kept
    return c

rows = []
for s in samples:
    ids = members[s]
    n_rrna = len(ids & rrna_removed)
    entering = ids - rrna_removed
    counts = dict.fromkeys(CLASSES, 0)
    for i in entering:
        counts[klass(i)] += 1
    if sum(counts.values()) != len(entering):
        sys.stderr.write(f"[longread] {s}: class counts do not sum to the sequences entering "
                         f"decontamination — the accounting is wrong\n")
    removed = sum(counts[c] for c in DROPPED)
    kept = len(entering) - removed
    rows.append({
        "sample": s, "condition": conditions.get(s, ""),
        "isoforms_in": len(ids), "rrna_removed": n_rrna,
        "host": counts["host"], "no_alignment": counts["no_alignment"],
        "contaminant": counts["contaminant"], "ambiguous": counts["ambiguous"],
        "contaminant_removed": removed, "kept": kept,
        "pct_removed": round(100.0 * removed / len(entering), 2) if entering else 0.0,
    })

cols = ["sample", "condition", "isoforms_in", "rrna_removed", "host", "no_alignment",
        "contaminant", "ambiguous", "contaminant_removed", "kept", "pct_removed"]
with open(args.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    pooled = len(set().union(*(members[s] for s in samples))) if samples else 0
    summed = sum(r["isoforms_in"] for r in rows)
    if summed > pooled:
        fh.write(f"# {summed - pooled:,} sequence-sample pairs are shared: a sequence detected in "
                 f"several samples is counted once per sample, so these columns sum to more than "
                 f"the {pooled:,} pooled sequences\n")

if unknown:
    sys.stderr.write(f"[longread] classes treated as contaminant: {', '.join(sorted(unknown))}\n")
print(f"[longread] {len(rows)} samples -> {args.out}")
for r in rows:
    print(f"    {r['sample']:<22} in {r['isoforms_in']:>8,}  rRNA {r['rrna_removed']:>6,}  "
          f"contaminant {r['contaminant_removed']:>7,} ({r['pct_removed']:.2f}%)  kept {r['kept']:>8,}")
