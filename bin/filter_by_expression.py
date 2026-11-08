#!/usr/bin/env python3
"""Cut the reference down to the transcripts the experiment actually saw.

Nothing else in the ladder filters on expression, so the ladder's cross-source union carries every
near-duplicate and every unexpressed fragment the four assembly sources produced. A transcript is
kept when it reaches --min-tpm (and/or --min-reads) in at least --min-samples libraries, measured
on a first quantification pass over the unfiltered set.

Two presence tests are available and either or both can be on. --min-tpm is length-normalised, so
it treats a short transcript and a long one with the same molar abundance alike; --min-reads is
raw evidence, so it will not admit a transcript that only looks abundant because it is short.
They disagree most on the extremes of the length distribution, which is exactly where spurious
assembly fragments live, so which one is used changes what the analysed set contains.

The threshold is a PRESENCE rule, not a mean: a transcript expressed in one treatment group and
absent elsewhere is exactly what a UV experiment is looking for, and averaging would bury it.

What leaves this script is the ANALYSED set. Everything after it — annotation, the reported salmon
quantification, DESeq2, enrichment, the figures and both reports — runs on this FASTA and this
tx2gene map, so every published number describes one set of transcripts and none has to be
qualified with which reference it came from.

What this does NOT cut is the assembly. Every ladder rung keeps its own FASTA, seqkit stats and
BUSCO summary, so the transcript count and the completeness quoted for the transcriptome remain
properties of the full deduplicated set (which is what a paper quotes); this script's output gets
the same two measurements, so the assembly and the analysed set can be read side by side.

Two files are written for the reference, and with no second route into it they list the same
transcripts:

  --out-fasta / --out-tx2gene   the quantification reference. Salmon indexes it and every
                                transcript-level number comes from it.
  --out-expressed-tx            the DE-eligible set. deseq2.R keeps only genes with at least one
                                transcript in this list, so the tested universe is "genes seen
                                above threshold in the experiment"; GO/KEGG ORA inherits it,
                                because enrichment.R takes its background from the genes present
                                in the DESeq2 results tables. The split is an interlock: a future
                                second way into the index (homology, say) would widen what is
                                quantified without silently widening what is tested.

An earlier version of this pipeline annotated the unfiltered set and offered --rescue-annotated,
a second way in for a transcript with an ORF and homology support. Both are gone. The rescue was
measured on this dataset and earned nothing — 0 DE genes gained, +0.127 pts of mapping rate
spread evenly over all 15 libraries (so it cannot move a contrast), and a DE reference that was
worse where it counts: 56,599 tx at C:72.1% S:45.5% D:26.6% with it, against 48,224 tx at
C:70.6% S:48.2% D:22.4% without, i.e. 48 more complete markers bought with 144 more duplicated
ones. Expression only gives the highest single-copy fraction of any reference in the ladder.
"""
import argparse
import csv
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("--quant-dir", required=True, help="directory of <sample>/quant.sf")
ap.add_argument("--fasta", required=True)
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--min-tpm", type=float, default=1.0,
                help="TPM a transcript must reach in a library to count as present there. "
                     "0 disables the TPM test (then --min-reads must be set).")
ap.add_argument("--min-reads", type=float, default=0.0,
                help="Estimated reads (salmon NumReads) a transcript must reach in a library to "
                     "count as present there. 0 disables the read test. When both tests are on, a "
                     "library counts only if the transcript passes BOTH in that same library.")
ap.add_argument("--min-samples", type=int, default=3)
ap.add_argument("--out-fasta", required=True)
ap.add_argument("--out-tx2gene", required=True)
ap.add_argument("--out-expressed-tx", default=None,
                help="write the transcript ids that passed the TPM rule. This is the "
                     "DE-eligible set; see the module docstring.")
ap.add_argument("--log", default=None)
args = ap.parse_args()

quants = sorted(
    os.path.join(args.quant_dir, d, "quant.sf")
    for d in os.listdir(args.quant_dir)
    if os.path.exists(os.path.join(args.quant_dir, d, "quant.sf"))
)
if not quants:
    sys.exit(f"[ERROR] no <sample>/quant.sf under {args.quant_dir}")

if args.min_tpm <= 0 and args.min_reads <= 0:
    sys.exit("[ERROR] both --min-tpm and --min-reads are 0, so nothing would be filtered; "
             "set at least one, or switch the filter off upstream")

# A library counts as "saw this transcript" only if it passes every ENABLED test in that same
# library. Evaluating the tests per library rather than on a per-transcript summary is what keeps
# this a presence rule: a transcript on in one treatment group and off elsewhere still survives.
n_pass = {}
for path in quants:
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if args.min_tpm > 0 and float(row["TPM"]) < args.min_tpm:
                continue
            if args.min_reads > 0 and float(row["NumReads"]) < args.min_reads:
                continue
            n_pass[row["Name"]] = n_pass.get(row["Name"], 0) + 1

expressed = {t for t, n in n_pass.items() if n >= args.min_samples}
keep = expressed

# The reference and the DE-eligible set are written separately although they are currently the
# same transcripts: the reference decides what salmon may assign a read to, the DE-eligible list
# decides what deseq2.R is allowed to test. Keeping them distinct is what makes a future second
# route into the index (homology, say) widen the first without silently widening the second.
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

_tests = []
if args.min_tpm > 0:
    _tests.append(f"TPM >= {args.min_tpm}")
if args.min_reads > 0:
    _tests.append(f"reads >= {args.min_reads:g}")
rule = f"{' AND '.join(_tests)} in >= {args.min_samples} of {len(quants)} libraries"
text = "\n".join([
    f"expression filter: {rule}",
    f"  transcripts kept    : {kept_seqs:,}",
    f"  transcripts dropped : {dropped_seqs:,}",
    f"  genes  {len(genes_before):,} -> {kept_genes:,}",
    f"  tx2gene rows written: {kept_tx:,}",
    "  this set is what gets annotated, quantified and tested",
] + ([f"  DE-eligible transcripts: {len(expressed):,}"]
     if args.out_expressed_tx else [])) + "\n"
if args.log:
    open(args.log, "w").write(text)
print(text, end="")

if kept_seqs == 0:
    sys.exit("[ERROR] expression filter kept nothing — check --min-tpm / --min-samples")
