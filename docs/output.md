# bioinfo-kaust/hybrid-metatranscriptomics: Output

Everything is written under `--outdir`. Start with `dashboard/index.html`, or `run_summary.txt` for
the plain-text version; both name the files their numbers come from. `<set>` below is a rung of
the reference ladder (`isoseq` … `all_sources_dedup`, and `<dge_reference>_expressed` after the
expression filter); `<sample>` is a library from `--input` or `--assembly_input`.

The expression filter runs **before** the annotation, so `<dge_reference>_expressed` is the set the
annotation, the published quantification, DESeq2, enrichment and the figures all describe. The
ladder's own FASTAs, seqkit stats and BUSCO summaries still cover every rung, so the transcript
count and the completeness quoted for the assembly describe the larger, unfiltered set.

Directories for steps that did not run (a missing genome, no contaminant genomes, `--skip_*`) are
simply absent.

## Reports

| Path | Content |
| --- | --- |
| `dashboard/index.html` | Self-contained interactive report of this run (README §3.1) |
| `run_summary.txt` | Headline numbers, the read-accounting check, and the source file of each number |
| `qc/multiqc/multiqc_report.html` | MultiQC over fastp, FastQC, STAR, BUSCO, salmon and the pipeline's own tables |
| `figures_<set>[_locus\|_corset]/` | Publication figures (PDF + 300 dpi PNG), one directory per gene map |
| `pipeline_info/` | Nextflow execution report, timeline, trace and DAG (timestamped), `params_*.json`, `hybridmetatranscriptomics_software_mqc_versions.yml`, `run_manifest.json` (every parameter plus the effective step switches; the dashboard reads it) |

## References — `references/`

Only what the pipeline downloaded or built: the genome and contaminant FASTAs (when given by
accession), `rrna_db.fasta`, the DIAMOND databases (`<name>.fasta.gz`, `<name>.dmnd`) and
`star_index/`. Staged inputs are not copied.

## Read QC — `qc/`

| Path | Content |
| --- | --- |
| `fastp/<sample>.fastp.{html,json,log}` | Adapter and quality trimming |
| `rrna/<sample>.bbduk.log`, `<sample>.rrna.tsv` | BBDuk rRNA screen, and its per-sample fraction |
| `fastqc/<sample>_{1,2}_fastqc.{html,zip}` | FastQC on the reads that go into assembly and quantification |
| `decontam/<sample>.{tsv,log}` | Per-sample decontamination accounting (pairs kept / removed) |
| `star/<sample>.Log.final.out` | STAR mapping summary (genome-guided route) |
| `rrna_per_sample.tsv`, `decontam_per_sample.tsv`, `*_assembly_only.tsv`, `longread_*.tsv` | The per-sample tables above, aggregated |

`--save_decontaminated_reads` additionally writes the host-only reads to `reads/host_only/`.

## Decontamination — `decontam/`

`decon_index.log` (the combined host + contaminant index), the long-read classification and
host-only long reads, and `audit/` — the residual-contamination audit of each assembly source
(`routeA`, `lrfree`, `denovo`): per-sequence classification, the kept ids and FASTA, and counts.

## Assembly — `assembly/`

| Path | Content |
| --- | --- |
| `isoseq/` | Route A: `collapsed.{gff,gtf,fasta}`, the genome-sequenced models (`collapsed.genome.fasta`), stats, `longreads_aligned.flagstat.txt` |
| `stringtie/` | Route C: per-sample `<sample>.gtf`, `stringtie_merged.gtf`, `gffcmp.stats`, the novel loci (`stringtie_novel.fasta`, `stringtie_novel_ids.txt`) |
| `denovo/` | Route D: `trinity.fa.gz`, `trinity.log` |
| `clustering/` | Novelty logs of the de novo and reference-free sets against the ladder |

## Reference transcriptome — `transcriptome/`

| Path | Content |
| --- | --- |
| `<set>.fasta` | One FASTA per ladder rung; `<dge_reference>_expressed.fasta` is what is annotated, quantified and tested |
| `tx2gene.<set>.tsv` | String-stripped gene map (an upper bound on distinct loci) |
| `tx2gene.<set>_locus.tsv`, `locus_rebuild_report.txt` | Genomic-locus gene map and its validation |
| `tx2gene.<set>_corset.tsv` | Corset read-sharing gene map |
| `annotation.<set>[_locus\|_corset].tsv` | The annotation re-keyed to each gene map |
| `<set>[_locus].genes.fasta`, `.genes.pep` | One representative sequence per gene of that map: its longest transcript, and the longest ORF on any of its transcripts |
| `<set>[_locus].gene_rep.tsv` | Which transcript and ORF represents each gene, with its transcript count (`rep_orf = NA` when the gene has no ORF) |
| `<set>.filter.txt`, `<set>.de_eligible_tx.txt` | Expression-filter log and the transcripts it kept (the analysed set) |
| `qc/<set>.stats.tsv` | seqkit statistics per rung |
| `qc/<set>.busco.txt`, `qc/short_summary.*`, `qc/*-busco.batch_summary.txt` | BUSCO per rung (`NA` when BUSCO failed) |
| `qc/route_contribution.tsv` | What each assembly source added to the ladder |
| `qc/sqanti3_isoseq*` | SQANTI3 structural classification of the long-read models |

## Annotation — `annotation/`

`annotation.<set>.tsv` (one row per transcript: ORF, protein hits, eggNOG, InterPro, GO, the
`--interest_label` flag) and its `.summary.txt`; TransDecoder's `<set>.pep`, `<set>.cds`,
`<set>.transdecoder.{gff3,bed}`; the raw searches in `diamond/diamond_<db>.tsv`,
`eggnog/eggnog.emapper.annotations` and `interproscan/interpro.tsv`.

`<set>` here is `<dge_reference>_expressed`: the annotation covers the transcripts the experiment
expressed, the same ones the DE tables are keyed on. Unexpressed gene models from the assembly are
not annotated, so this table is not a catalogue of everything the species encodes.

## Quantification — `quantification/`

`<sample>/` is salmon's output (`quant.sf`, `aux_info/` with the equivalence classes, bootstraps)
for the analysed set; the pre-filter pass that only produced TPMs for the filter is not published.
`mapping_rates.tsv`; Corset's clusters and counts (`corset_<set>/`,
`<set>_corset.summary.txt`).

## Differential expression

| Path | Content |
| --- | --- |
| `dge/dge_<set>/` | DESeq2 on the string-stripped map: one table per contrast, `DE_summary.csv`, `filtering_summary.csv`, normalised counts, PCA |
| `dge_locus/dge_<set>_locus/` | The same on the locus map (the default primary map) |
| `dge_corset/dge_<set>_corset/` | The same on the Corset map |
| `enrichment[_locus\|_corset]/enrichment_<set>*/` | GO / KEGG over-representation and GSEA per contrast, `enrichment_summary.csv` |
| `isoforms/` | Differential isoform usage (only without `--skip_dtu`) |
