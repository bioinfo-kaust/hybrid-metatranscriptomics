# hybrid-metatranscriptomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-11-08

The pipeline now follows the nf-core template: most of what changes is how it is laid out,
configured and provisioned, with every step, argument and threshold carried over. The one
exception is the order of the expression filter and the annotation, which does change what is
analysed and is listed under **Changed** below.

### Added

- **`COUNT_MATRIX`: the count matrices themselves, one set per gene map, in
  `quantification/gene_level/`.** The run published per-sample transcript `quant.sf` and, from
  DESeq2, normalised and VST counts — but those are transformed, and they cover only the genes
  that passed the presence pre-filter, so the raw gene-level matrix existed nowhere outside the
  DESeq2 task and the gene FASTA shipped without the counts belonging to it. The pipeline does
  not quantify genes directly: it quantifies transcripts and sums them with tximport, and this
  writes that aggregation out, unfiltered — one row per gene in the map, including the ones DE
  later drops. Published per map: `gene_counts` (what DESeq2 modelled), `gene_length` (the
  `avgTxLength` offset it is paired with, without which the counts are a different analysis),
  `gene_counts_lengthScaledTPM` (that correction folded in, for tools that cannot take an
  offset) and `gene_tpm`, plus the three transcript-level matrices once.
  The process is handed each map's `dds.rds` and **fails rather than publishing** unless
  `round(counts)` is cell-for-cell identical to the count assay DESeq2 modelled and the length
  matrix matches its offset exactly, so the published matrix is the DE input rather than a
  plausible reconstruction of it; `gene_level/checks.txt` records the comparison. `README.txt`
  also warns off `<sample>/quant.genes.sf`, which salmon writes from the string-stripped map
  because the locus rebuild does not exist yet at quantification time, and which no reported
  result uses.
- **Every gene count on the dashboard and in `run_summary.txt` now names the map it comes from.**
  Three stages report a gene count from the string-stripped tx2gene because that is the map they
  were handed — the expression filter, the annotation module, and the dedup step feeding the
  methods diagram — and the headline quotes the primary (locus) map, so an unlabelled figure read
  as the headline and overstated the gene count by every split the locus rebuild reconciles
  (34,853 against 29,709 on this dataset). The expression-filter and annotation cards now carry
  both counts, labelled, because the analysed set is exactly the set the locus rebuild is defined
  on; the assembly-wide count stays string-stripped and says why (the rebuild needs genome
  placement). Every other row of the annotation summary counts transcripts and is unaffected by
  the choice of map, which the tab now states rather than leaving a reader to work out.
- **The DE filtering card now accounts for the whole gap** between the headline gene count and
  `tested`, not just DESeq2's share of it: the presence pre-filter applied before `DESeq2()` is
  called, then Cook's-distance outliers (a property of the design, so the same genes leave every
  contrast), then per-contrast independent filtering. Only the last is contrast-specific, and it
  is the only one the card used to mention.
- **An Explorer tab on the dashboard, carrying the expression matrix itself** rather than a
  summary of it. The page ships gene-level normalised counts for every gene map, transcript-level
  salmon TPM, and per-contrast fold changes both as reported (apeglm-shrunken) and as DESeq2
  estimated them before shrinkage. From it: search genes or transcripts by id or annotation,
  filter by contrast / significance / direction, select rows and plot exactly that selection
  (condition means with SD, per-library dots, per-library bars, profiles, z-scored heatmap),
  export any chart as PNG or SVG and the underlying values as TSV, break a gene down into its
  transcripts to see which isoform carries the signal, and compare the two fold-change estimates
  gene by gene with live counts of what each one calls significant.
  The matrices are quantised (uint16 `log2(x+1) * 1000` for expression, int32 at 1e-4 for fold
  changes and 1e-5 for `-log10(padj)`), gzipped, base64'd and inflated in the browser on first
  use, so the page stays one self-contained file that works from `file://` and a reader who never
  opens the tab pays nothing for it. The expression quantisation is ~0.04% of value, which the
  tab states; the DE quantisation is fine enough that the tab's live significance counts reproduce
  `DE_summary.csv` exactly — verified for all three gene maps × three contrasts on `tested`,
  `sig_padj`, `sig`, `up` and `down`. This cost ~7 MB on a 12 MB page.
  The fold-change cut is applied strictly (`|log2FC| > cut`) to match `deseq2.R`; at 1/100 of a
  log2 unit the counts had disagreed with `DE_summary.csv` by a gene on two contrasts, which is
  why the DE arrays are int32 rather than int16.
- The dashboard is now staged the per-sample salmon directories, which is where the
  transcript-level TPM matrix comes from.
- `transcriptome/qc/<set>.full_table.tsv` is published for every BUSCO run: which ortholog each
  sequence carries. This is the label-free truth the gene-map identity evaluation scores the three
  maps against, and until it was published that evaluation could only be run out of a BUSCO task's
  work directory, which tied it to one particular run of the pipeline. The ids in it are the
  serial ids `BUSCO_RENAME` wrote (`t1..tN` in file order), resolved back through the rung's own
  FASTA order.
  `BUSCO_SUMMARY` copies it out rather than `BUSCO_BUSCO` publishing it directly, and that is not
  a stylistic choice: `BUSCO_BUSCO` declares the whole `*-busco` directory as an output as well as
  the table nested inside it, and Nextflow then treats the directory as the unit of publication,
  so a `publishDir` pattern matching the nested file publishes nothing at all. Verified against
  that exact output shape — an explicit glob, a `**` glob and a `saveAs`-only filter all published
  nothing, while publishing the directory did carry the table but would cost 123 MB per rung for a
  512 KB file.
- `--expr_filter_reads`: a read-count presence test for the expression filter, alongside
  `--expr_filter_tpm`. TPM is length-normalised, so a long transcript at low molar abundance can
  fall below TPM 1 while carrying plenty of reads, and a short one can clear TPM 1 on very few;
  the two therefore disagree most at the extremes of the length distribution, which is where
  spurious assembly fragments live. Either test may be used alone, or both together (a library
  must then pass both for that library to count). Both at 0 switches the filter off. Default 0,
  so existing behaviour is unchanged. Measured on the Cassiopea data, swapping TPM >= 1 for
  reads >= 5 (same ">= 3 libraries" rule) keeps 49,768 transcripts instead of 44,461: it admits
  5,564 with median length 1,712 and drops 257 with median length 338.
- `--genemap_eval`: a directory of externally-computed gene-map evaluation tables
  (`genemap_eval.tsv`, `genemap_null.tsv`, `genemap_verdict.txt`). Which transcript-to-gene map is
  primary cannot be decided from a run's own outputs, because a map that calls more DE genes is
  not thereby the better map; deciding needs independent truth and a label-permutation null. When
  this is set the dashboard renders a "Why <map> is the primary map" card comparing the maps on
  F1 against each truth, under- and over-merge rates and false DE calls per label shuffle. When it
  is unset the card is simply absent, so the pipeline stays generic.

- nf-core template layout: `main.nf` → `workflows/hybridmetatranscriptomics.nf` →
  `subworkflows/local/*/main.nf`, one directory per module, `nextflow_schema.json`, nf-schema
  validation and `--help`, completion email/summary, institutional configs (`-profile <site>`)
- nf-core modules for fastp, BBDuk, FastQC, cat/fastq, STAR (index + align), StringTie (+ merge),
  GffCompare, pbmm2, samtools flagstat, Trinity, vsearch, BUSCO, seqkit (stats, replace),
  TransDecoder, DIAMOND (makedb + blastp), eggNOG-mapper, InterProScan, salmon (index + quant),
  gunzip and MultiQC; versions pinned in `modules.json`
- Container and conda environment on every process: `-profile docker`, `singularity`,
  `apptainer`, `podman`, `shifter`, `charliecloud`, `wave`, `conda`, `mamba`. The pipeline's own
  modules use Seqera Containers built from their `environment.yml`
- Samplesheet schemas (`assets/schema_input.json`, `assets/schema_longread.json`): columns and
  file existence are checked at start-up
- `GENE_FASTA`: a gene-level FASTA (`<set>.genes.fasta`) and peptide file (`.genes.pep`) for the
  string-stripped and the genomic-locus gene map, plus `<set>.gene_rep.tsv` naming each gene's
  representative. A gene map alone leaves a gene as a sum of counts with no sequence to BLAST,
  score in protein mode or attach to a paper; the representative is the longest transcript of the
  gene, and for the peptide the longest ORF on any of its transcripts
- `--interproscan_db`: an InterProScan `data/` directory (default: the installation's own)
- `--slurm_account`, `--slurm_partition`, `--slurm_extra` declared in the schema (used by `-profile hpc`)
- `assets/test/contaminant.fa`, so `-profile test` can exercise decontamination

### Fixed

- **ORA was gated on the number of significant genes, not on the number that carry a term**, so a
  near-null contrast could return a page of "enriched" GO terms driven by one gene. Measured on
  the Cassiopea read-filter run: `UVAUVB_vs_UVB` had 5 significant genes of which 1 was annotated,
  and `enrichment.R` emitted 38 GO terms and 1 KEGG pathway as significant — every GO row at
  `GeneRatio 1/1`, all 38 sharing a single `p.adjust`, because a one-gene query makes every term
  that gene carries look enriched against a 13,272-gene background and BH over 38 tests adjusts
  almost nothing away. The gate now counts the significant genes present in the term map (GO and
  KEGG separately) and logs the reason when it declines. Re-running the affected contrast drops
  the two spurious files and leaves every other contrast byte-identical.

### Changed

- **The expression filter now runs before the annotation, and defines the analysed set.** The
  order was ladder → annotate the full deduplicated set → salmon pass 1 → filter → salmon pass 2
  → DE; it is now ladder → salmon pass 1 → filter → annotate the filtered set → salmon pass 2 →
  DE. The annotation, the published quantification, DESeq2, enrichment and the figures therefore
  all describe one set of transcripts (on the Cassiopea run: 44,454 expressed, from 149,385 in
  `all_sources_dedup`), so an annotation count and a DE count are comparable without a footnote.
  Consequences: `annotation.<dge_reference>.tsv` becomes
  `annotation.<dge_reference>_expressed.tsv` and no longer carries rows for transcripts the
  experiment did not express, so it is not a catalogue of everything the species encodes; the
  annotation searches (TransDecoder, DIAMOND, eggNOG-mapper, InterProScan) run on roughly a third
  of the transcripts and start later, after the filter. The ASSEMBLY is untouched: every ladder
  rung still publishes its own FASTA, seqkit stats and BUSCO summary, and the full deduplicated
  set keeps its own count and completeness score on the dashboard's Reference tab, beside the
  analysed set rather than in place of it.
- **The dashboard now leads with the analysed set, and quotes BUSCO for it.** The headline figure
  was the assembly's transcript count (149,385 on the Cassiopea run) with the assembly's BUSCO;
  it is now the expression-filtered count (44,454) with that set's gene count and its own BUSCO
  (78.0% complete vs 81.0% on the assembly), and names the assembly it was filtered out of. The
  `Reference transcriptome` card leads with `transcripts analysed` and its completeness, and
  carries the assembly as a `the assembly, before the filter` row underneath. The reasoning: the
  analysed set is the denominator every other number on the page shares, so heading the page with
  the larger count invited readers to scale the DE and annotation counts against transcripts that
  were never tested. BUSCO on the two sets answers different questions and both are still shown —
  the assembly's is how complete the reconstruction is, the analysed set's is how much of the
  conserved gene space the experiment actually expressed.
  What this does NOT change: the DE reference, the quantification and the DESeq2 results, because
  the old filter already cut the annotation table to the kept transcripts (44,454 rows) before
  enrichment saw it. The homology searches are per query, so their hits are unaffected; the one
  genuine difference is that TransDecoder trains its coding-likelihood model on the set it is
  given, so a small number of ORF calls — and the GO terms that follow from them — can move
- New subworkflow `subworkflows/local/expression_filter` holds salmon pass 1 and the filter; the
  pre-filter pass is no longer published (it only feeds the filter), so `quantification/` now
  holds one quant.sf per sample instead of the second pass overwriting the first
- **Parameters renamed:** `--samplesheet` → `--input`, `--longread_samplesheet` →
  `--longread_input`, `--assembly_samplesheet` → `--assembly_input`
- `--outdir` is required (no default)
- Per-task caps are `process.resourceLimits` (conf/base.config, conf/hpc.config) instead of
  `--max_cpus` / `--max_memory` / `--max_time`
- Every task runs under `bash -e -u -o pipefail`: a failing command inside a pipe now fails the task
- Downloads use `curl -sfL`: an HTTP error fails the task instead of saving the error page
- `--r_libs`, `--corset_bin` and `--extra_path` default to null; they are only for `-profile
  envmodules`, and `--r_libs` must not be set with containers
- MultiQC: default config is `assets/multiqc_config.yml`; set the report title with `--multiqc_title`
- Trinity's `--max_memory` is 80% of the task memory (200G at the default 250 GB) instead of a parameter
- Published outputs that changed name or form: `assembly/denovo/trinity.fa.gz` + `trinity.log`
  (was `trinity.Trinity.fasta` + `trinity.stats.tsv`); BUSCO publishes `short_summary.*` and
  `*-busco.batch_summary.txt` in `transcriptome/qc/` (was one `busco_<set>/` directory per set);
  `pipeline_info/` file names carry a timestamp, and the software versions file is
  `hybridmetatranscriptomics_software_mqc_versions.yml`; FastQC reports are `<sample>_1_fastqc.*` /
  `<sample>_2_fastqc.*` (was `<sample>_R1.rrna_clean_fastqc.*`)
- STAR alignments carry an `@RG` line, StringTie writes gene abundances (`-A`), and salmon is given
  the tx2gene map (`--geneMap`): none of these changes a count

### Removed

- `--expr_filter_rescue` (the ORF+homology second way into the DE reference). It needed the
  annotation to decide what to keep, which is no longer available when the filter runs, and it had
  already been measured on this dataset as buying nothing: 0 DE genes gained, +0.127 pts of
  mapping rate spread evenly over all 15 libraries, and a worse reference where it counts (56,599
  tx at C:72.1% S:45.5% D:26.6% with it, against 48,224 tx at C:70.6% S:48.2% D:22.4% without).
  Expression alone gives the highest single-copy BUSCO fraction of any reference in the ladder
- `--trinity_max_memory` (see above), `--max_cpus`, `--max_memory`, `--max_time`, `--tracedir`
- `params.conda_pkg`, `--use_conda_pkgs` and `conf/conda.config`: `-profile conda` uses each
  module's `environment.yml`
- `lib/Inputs.groovy`: the input helpers live in `subworkflows/local/utils_*`
- No longer published: the vsearch cluster logs (`assembly/clustering/*.cluster.log`) and BBDuk's
  per-database `*.bbduk_stats.txt` / `*.bbduk_refstats.txt` (the rRNA fraction comes from the BBDuk log)

### Known issues

- Not yet validated on real data: stubs, nf-core lint and `nextflow lint` pass; a subsampled
  Cassiopea run (`-profile hpc,envmodules`) reached annotation and stopped at
  `TRANSDECODER_PREDICT` (its symlink staging of the LONGORF directory fails on the compute
  node). A full run must then be compared with `results_nf_29Sep26/run_summary.txt`.
- `LOCUS_TX2GENE` has no container yet (placeholder tag `LOCUS_*_PENDING`; Wave was rate-limited):
  it runs under `-profile envmodules` only until the image is built.

## 1.1.0dev — optional inputs (merged into 2.0.0)

### Changed

- `--genome` / `--genome_accession` are optional. Without a genome the genome-guided long-read
  route, StringTie, SQANTI3 and the locus gene map are skipped (`dge_primary_map locus` falls back
  to `string`); Trinity and the reference-free long reads build the reference on their own, salmon
  indexes without decoys, and decontamination removes every read that places on a contaminant
- The run manifest records the effective step switches, so the dashboard reports steps that a
  missing input turned off as skipped

### Fixed

- `TRINITY` memory: the nf-core module passes 80% of `task.memory` as `--max_memory`, so the
  ported 250 GB request gave Trinity 200G where the local module had hard-coded 240G. Measured
  cost of the 40 GB: 15 h 23 m against 9 h 18 m for the same assembly. Now requested at 300 GB
- No contaminant genomes with `skip_decontam` left at `false` stalled the run silently after
  FastQC (no assembly, no DGE, still reported success); decontamination is now skipped instead
- `--longread_samplesheet null` was opened as a file called `null`
- A run left with no assembly route now fails at start-up instead of producing nothing
- A missing `corset` executable (no environment module provides it) now fails at start-up, naming
  `--corset_bin`, `-profile conda` and `--skip_corset`, instead of failing at the CORSET step

## 1.0.0 - [2026-10-04]

First public release. Builds a host-only reference transcriptome from long + short reads,
removing contaminant (e.g. symbiont) reads _before_ assembly, then annotates, quantifies and
tests it for differential expression. Developed on a _Cassiopea andromeda_ UV-radiation
experiment.

### Added

- Competitive read-level decontamination: host and contaminant genomes in one minimap2 index
- Four assembly routes (Iso-Seq on genome, reference-free long reads, StringTie, Trinity),
  layered into a reference ladder with BUSCO on every rung
- Annotation: TransDecoder, DIAMOND, eggNOG-mapper, InterProScan
- salmon → DESeq2 → GO/KEGG enrichment, on three gene maps (genomic locus by default)
- Interactive HTML dashboard, MultiQC and a plain-text run summary
- Profiles: `hpc` / `slurm` (executor), `envmodules` / `conda` (software), `test` (stub run)
- `assets/params_cassiopea_uv.yml` reproduces the Cassiopea run
