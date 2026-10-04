#!/usr/bin/env python3
"""Drop the transcripts that are neither expressed nor protein-coding.

Nothing else in the ladder filters on expression, so the quantification reference carries every
near-duplicate and every unexpressed fragment the four assembly sources produced. A transcript is
kept when it reaches --min-tpm in at least --min-samples libraries, measured on the first
quantification pass over the unfiltered set.

The threshold is a PRESENCE rule, not a mean: a transcript expressed in one treatment group and
absent elsewhere is exactly what a UV experiment is looking for, and averaging would bury it.

--rescue-annotated adds a second way to survive: an ORF *and* homology support (a protein hit or
an InterPro domain). It is OFF by default. It was measured on this dataset and does not earn its
place in the DE reference:

  DE genes gained      0 -- no rescued-only gene is significant in any gene map or contrast
  mapping rate gained  +0.127 pts (92.14% -> 92.27%), identical to within 0.07 pts across all 15
                       libraries and all three groups, so it cannot move a contrast
  BUSCO                it makes the DE reference WORSE by the measure that matters here:
                         expression only  48,224 tx  C:70.6% S:48.2% D:22.4%
                         + rescue         56,599 tx  C:72.1% S:45.5% D:26.6%
                       48 more complete markers bought with 144 more duplicated ones. Expression
                       only gives the highest single-copy fraction of any reference in the ladder.

The rescue used to be justified by that 1.5-point completeness gap, but that was the wrong number
to defend: completeness is a property of the ASSEMBLY, and the assembly is the unfiltered ladder
tip (all_sources_dedup, 79,835 tx, C:72.1%), which is what BUSCO reports on and what the paper
quotes. What this script produces is the reference for quantification and testing, and it has no
obligation to reproduce the assembly's completeness.

Nor are the rescued transcripts thrown away. ANNOTATE runs upstream of this filter, so all 79,835
keep their annotation in annotation.<reference>.tsv -- including the 8,375 the rescue would have
added back, which carry 24% of the reference's protein hits, 25% of its InterPro domains and 34%
of its UV_stress_family calls. Roughly half of them sit in loci with no expressed isoform at all:
real gene models, dormant in this experiment. They stay in the annotation as "this species encodes
X"; they are simply not quantified and not tested.

Two sets leave this script. With the rescue off they are identical, and the split is a safety
interlock rather than a live distinction:

  kept = expressed | rescued   -> --out-fasta / --out-tx2gene / --out-annotation
      The quantification reference. Salmon indexes it and every transcript-level number for the
      DE analysis comes from it.

  expressed                    -> --out-expressed-tx
      The DE-eligible set. deseq2.R keeps only genes with at least one transcript in this list,
      which makes the tested universe "genes seen above threshold in the experiment" while the
      counts for those genes still include reads landing on any rescued isoforms. GO/KEGG ORA
      then inherits that universe automatically, because enrichment.R takes its background from
      the genes present in the DESeq2 results tables. Re-enabling --rescue-annotated therefore
      widens the index without silently widening what gets tested.
"""
import argparse
import csv
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--quant-dir", required=True, help="directory of <sample>/quant.sf")
ap.add_argument("--fasta", required=True)
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--min-tpm", type=float, default=1.0)
ap.add_argument("--min-samples", type=int, default=3)
ap.add_argument("--out-fasta", required=True)
ap.add_argument("--out-tx2gene", required=True)
ap.add_argument("--annotation", default=None, help="per-transcript annotation to restrict too")
ap.add_argument("--out-annotation", default=None)
ap.add_argument("--out-expressed-tx", default=None,
                help="write the transcript ids that passed the TPM rule (rescued ones "
                     "excluded). This is the DE-eligible set; see the module docstring.")
ap.add_argument("--rescue-annotated", action="store_true",
                help="also keep a transcript carrying an ORF AND homology support "
                     "(needs --annotation). Off by default; see the module docstring.")
ap.add_argument("--log", default=None)
args = ap.parse_args()

quants = sorted(
    os.path.join(args.quant_dir, d, "quant.sf")
    for d in os.listdir(args.quant_dir)
    if os.path.exists(os.path.join(args.quant_dir, d, "quant.sf"))
)
if not quants:
    sys.exit(f"[ERROR] no <sample>/quant.sf under {args.quant_dir}")

n_pass = {}
for path in quants:
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if float(row["TPM"]) >= args.min_tpm:
                n_pass[row["Name"]] = n_pass.get(row["Name"], 0) + 1

expressed = {t for t, n in n_pass.items() if n >= args.min_samples}

# Second route in: coding evidence. "Homology" is a protein hit OR an InterPro domain — either is
# independent evidence the ORF is a real gene rather than a chance open frame.
rescued = set()
if args.rescue_annotated:
    if not args.annotation:
        sys.exit("[ERROR] --rescue-annotated needs --annotation")
    with open(args.annotation) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = {"transcript", "has_ORF", "protein_hit", "InterPro"} - set(rd.fieldnames or [])
        if missing:
            sys.exit(f"[ERROR] --rescue-annotated needs columns {sorted(missing)} in {args.annotation}")
        for row in rd:
            orf = (row.get("has_ORF") or "").strip().lower() in ("yes", "true", "1")
            homology = bool((row.get("protein_hit") or "").strip()) or \
                       bool((row.get("InterPro") or "").strip())
            if orf and homology:
                rescued.add(row["transcript"])

keep = expressed | rescued

# The two sets have different jobs downstream and must not be conflated:
#   keep      -> the quantification reference. A rescued transcript is a real gene model with
#                an ORF and homology support, so leaving it out of the index pushes its reads onto
#                whatever near-duplicate remained. Measured here that costs 0.127 pts of mapping
#                rate, evenly across every library -- which is why the rescue is off by default.
#   expressed -> the DE-eligible set. A transcript the experiment never saw above threshold
#                carries no evidence about treatment, and testing it only costs multiple-testing
#                burden. Genes are restricted to this set in deseq2.R, NOT here, so that reads
#                landing on a rescued transcript still count toward its gene.
if args.out_expressed_tx:
    with open(args.out_expressed_tx, "w") as out:
        for t in sorted(expressed):
            out.write(t + "\n")

kept_seqs = dropped_seqs = 0
with open(args.fasta) as fh, open(args.out_fasta, "w") as out:
    writing = False
    for line in fh:
        if line.startswith(">"):
            tid = line[1:].split()[0]
            writing = tid in keep
            kept_seqs += writing
            dropped_seqs += not writing
        if writing:
            out.write(line)

kept_tx = kept_genes = 0
genes_before = set()
genes_after = set()
with open(args.tx2gene) as fh, open(args.out_tx2gene, "w") as out:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        genes_before.add(parts[1])
        if parts[0] in keep:
            out.write(line)
            kept_tx += 1
            genes_after.add(parts[1])
kept_genes = len(genes_after)

ann_rows = 0
if args.annotation and args.out_annotation:
    # keep the annotation in step with the reference: leaving rows for dropped transcripts would
    # let them keep their old gene ids and inflate every gene-level count computed from this table
    with open(args.annotation) as fh, open(args.out_annotation, "w", newline="") as out:
        rd = csv.reader(fh, delimiter="\t")
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        header = next(rd, None)
        if header:
            w.writerow(header)
            for row in rd:
                if row and row[0] in keep:
                    w.writerow(row)
                    ann_rows += 1

rule = f"TPM >= {args.min_tpm} in >= {args.min_samples} of {len(quants)} libraries"
if args.rescue_annotated:
    rule += " OR (ORF AND homology)"
text = "\n".join([
    f"expression filter: {rule}",
    f"  transcripts kept    : {kept_seqs:,}",]
+ ([f"    expressed         : {len(expressed):,}",
    f"    rescued by ORF+hom: {len(rescued - expressed):,}"] if args.rescue_annotated else []) +
   [
    f"  transcripts dropped : {dropped_seqs:,}",
    f"  genes  {len(genes_before):,} -> {kept_genes:,}",
    f"  tx2gene rows written: {kept_tx:,}",
] + ([f"  annotation rows kept: {ann_rows:,}"] if args.annotation else [])
  + ([f"  DE-eligible transcripts (expressed only): {len(expressed):,}"]
     if args.out_expressed_tx else [])) + "\n"
if args.log:
    open(args.log, "w").write(text)
print(text, end="")

if kept_seqs == 0:
    sys.exit("[ERROR] expression filter kept nothing — check --min-tpm / --min-samples")
