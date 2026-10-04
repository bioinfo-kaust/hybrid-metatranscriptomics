#!/usr/bin/env bash
# Compare a finished pipeline run against a published reference transcriptome.
#
# This lives OUTSIDE the workflow on purpose. Which published reference counts as authoritative
# is an interpretation choice, it changes as new references appear, and it must not silently
# alter pipeline outputs on a re-run. Run it whenever you like, against as many references as
# you like, without touching the pipeline.
#
# Usage:
#   extra/compare_to_published.sh <pipeline_outdir> <published.fasta> [out_dir] [reference_annotation.tsv]
#
# Example (the Cassiopea run vs the Gacesa et al. TSA assembly):
#   extra/compare_to_published.sh results_nf \
#       refs/gacesa/Cassiopea_andromeda_TSA_GJJJ.fasta \
#       results_nf/comparison_gacesa
set -euo pipefail

OUTDIR=${1:?usage: compare_to_published.sh <pipeline_outdir> <published.fasta> [out_dir] [ref_annotation.tsv]}
PUBLISHED=${2:?published reference FASTA}
DEST=${3:-$OUTDIR/comparison}
REF_ANN=${4:-}
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
SET=${SET:-final_nr}

QUERY=$OUTDIR/transcriptome/${SET}.fasta
TX2GENE=$OUTDIR/transcriptome/tx2gene.${SET}.tsv
ANNOT=$(ls "$OUTDIR"/annotation/annotation.*.tsv 2>/dev/null | head -1)
DGE=$(ls -d "$OUTDIR"/dge/dge_* 2>/dev/null | head -1)

for f in "$QUERY" "$TX2GENE" "$ANNOT"; do
  [ -s "$f" ] || { echo "[ERROR] missing pipeline output: $f" >&2; exit 1; }
done
mkdir -p "$DEST"

# `-max_target_seqs 1` does NOT return the best hit — it returns the first one found above the
# threshold. Ask for 25 and let the bitscore selection in the classifier decide.
if [ ! -s "$DEST/hits.tsv" ]; then
  makeblastdb -in "$PUBLISHED" -dbtype nucl -out "$DEST/published_db" > /dev/null
  blastn -query "$QUERY" -db "$DEST/published_db" -evalue 1e-5 -max_target_seqs 25 \
         -num_threads "$THREADS" \
         -outfmt "6 qseqid sseqid pident length qlen evalue bitscore" > "$DEST/hits.tsv"
fi
echo "blast hits: $(wc -l < "$DEST/hits.tsv")"

python3 "$(dirname "$0")/novelty_vs_reference.py" \
    --blast "$DEST/hits.tsv" \
    --annotation "$ANNOT" \
    --tx2gene "$TX2GENE" \
    ${DGE:+--dge-dir "$DGE"} \
    ${REF_ANN:+--reference-annotation "$REF_ANN"} \
    --out "$DEST/novelty_per_gene.tsv" \
    --report "$DEST/novelty_summary.txt"

# side-by-side metrics
{
  printf 'metric\tthis run\tpublished reference\n'
  printf 'transcripts\t%s\t%s\n' "$(grep -c '^>' "$QUERY")" "$(grep -c '^>' "$PUBLISHED")"
  printf 'genes\t%s\t%s\n' "$(cut -f2 "$TX2GENE" | sort -u | wc -l)" "-"
  printf 'transcripts with a protein hit\t%s\t%s\n' \
    "$(awk -F'\t' 'NR>1 && $5!=""' "$ANNOT" | wc -l)" "-"
} > "$DEST/metrics.tsv"
column -t -s$'\t' "$DEST/metrics.tsv"
echo "[done] -> $DEST/"
