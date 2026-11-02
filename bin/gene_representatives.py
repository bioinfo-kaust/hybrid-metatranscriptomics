#!/usr/bin/env python3
"""One sequence per GENE, for a reference that only has one sequence per transcript.

Everything downstream of a gene map works on counts: the tx2gene file says which transcripts
belong to a gene, DESeq2 sums them, and no sequence for the gene itself ever exists. Anything that
needs the gene as a sequence — BLAST against another catalogue, orthology, a protein-mode BUSCO,
primer design for a DE locus, a supplementary FASTA for the paper — then has to be rebuilt by hand
from the transcript FASTA and the map. This writes it once, next to the map it came from.

A locus-defined gene has no single true sequence, so the choice is stated rather than hidden:

  nucleotide   the LONGEST transcript of the gene. Longest carries the most sequence into a
               homology search, which is what a gene-level FASTA is mostly used for. It is a
               representative, not a consensus: an alternative isoform's unique exons are not in
               it, so transcript-level questions stay with the transcript FASTA.
  peptide      the longest TransDecoder ORF found on ANY transcript of the gene (--pep). Taking
               the ORF of the representative transcript instead would drop a gene whose longest
               transcript happens to be the one TransDecoder found no ORF on. A gene with no ORF
               anywhere is absent from the peptide file and carries rep_orf=NA in the table.

The table lists every gene either way, so the peptide file's shortfall is readable as a number
rather than inferred from two line counts.
"""
import argparse
import collections


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def wrap(out, seq, width=60):
    for i in range(0, len(seq), width):
        out.write(seq[i:i + width] + "\n")


ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True, help="transcript FASTA of the analysed reference")
ap.add_argument("--tx2gene", required=True)
ap.add_argument("--pep", default=None,
                help="TransDecoder peptides (<transcript>.p<N> ids). Omit or pass the NO_FILE "
                     "sentinel to write nucleotide representatives only.")
ap.add_argument("--out-fasta", required=True)
ap.add_argument("--out-pep", required=True)
ap.add_argument("--out-table", required=True)
args = ap.parse_args()

t2g = {}
with open(args.tx2gene) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            t2g[parts[0]] = parts[1]

# longest transcript per gene
best = {}
n_tx = collections.Counter()
missing_map = 0
for header, seq in read_fasta(args.fasta):
    tid = header.split()[0]
    gene = t2g.get(tid)
    if gene is None:
        missing_map += 1
        continue
    n_tx[gene] += 1
    if gene not in best or len(seq) > best[gene][0]:
        best[gene] = (len(seq), tid, seq)

if not best:
    raise SystemExit("[ERROR] no transcript in --fasta appears in --tx2gene")

# longest ORF per gene, over all of its transcripts
best_pep = {}
if args.pep and not args.pep.startswith("NO_FILE"):
    for header, seq in read_fasta(args.pep):
        orf_id = header.split()[0]
        gene = t2g.get(orf_id.rsplit(".p", 1)[0])
        if gene is None:
            continue
        pep = seq.rstrip("*")
        if gene not in best_pep or len(pep) > len(best_pep[gene][1]):
            best_pep[gene] = (orf_id, pep)

genes = sorted(best)
with open(args.out_fasta, "w") as out:
    for gene in genes:
        length, tid, seq = best[gene]
        out.write(f">{gene} rep_transcript={tid} len={length} n_tx={n_tx[gene]}\n")
        wrap(out, seq)

with open(args.out_pep, "w") as out:
    for gene in genes:
        if gene in best_pep:
            orf_id, pep = best_pep[gene]
            out.write(f">{gene} rep_orf={orf_id} len={len(pep)}\n")
            wrap(out, pep)

with open(args.out_table, "w") as out:
    out.write("gene\tn_transcripts\trep_transcript\trep_tx_len\trep_orf\trep_orf_aa\n")
    for gene in genes:
        length, tid, _seq = best[gene]
        orf_id, pep = best_pep.get(gene, ("NA", ""))
        out.write(f"{gene}\t{n_tx[gene]}\t{tid}\t{length}\t{orf_id}\t{len(pep)}\n")

lengths = sorted((v[0] for v in best.values()), reverse=True)
total = sum(lengths)
running = 0
n50 = next(length for length in lengths if (running := running + length) >= total / 2)
print("\n".join([
    f"gene representatives: {len(genes):,} genes from {sum(n_tx.values()):,} transcripts",
    f"  with a peptide    : {len(best_pep):,}",
    f"  total nt          : {total:,}",
    f"  N50 / median      : {n50:,} / {lengths[len(lengths) // 2]:,}",
] + ([f"  [WARN] not in the map: {missing_map:,} transcripts"] if missing_map else [])))
