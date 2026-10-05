# hybrid-metatranscriptomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
