# hybrid-metatranscriptomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2.0.0dev

The pipeline now follows the nf-core template: most of what changes is how it is laid out,
configured and provisioned, with every step, argument and threshold carried over. The one
exception is the order of the expression filter and the annotation, which does change what is
analysed and is listed under **Changed** below.

### Added

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
  rung still publishes its own FASTA, seqkit stats and BUSCO summary, so the transcript count and
  the completeness quoted for the transcriptome remain properties of the full deduplicated set.
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

## 1.1.0dev — optional inputs (merged into 2.0.0dev)

### Changed

- `--genome` / `--genome_accession` are optional. Without a genome the genome-guided long-read
  route, StringTie, SQANTI3 and the locus gene map are skipped (`dge_primary_map locus` falls back
  to `string`); Trinity and the reference-free long reads build the reference on their own, salmon
  indexes without decoys, and decontamination removes every read that places on a contaminant
- The run manifest records the effective step switches, so the dashboard reports steps that a
  missing input turned off as skipped

### Fixed

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
