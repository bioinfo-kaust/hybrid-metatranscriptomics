# hybrid-metatranscriptomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
