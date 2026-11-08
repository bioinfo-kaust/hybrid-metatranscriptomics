#!/usr/bin/env Rscript
# Gene-level and transcript-level count matrices, one set per gene map.
#
# WHY THIS EXISTS: the pipeline publishes per-sample salmon quant.sf (transcript level) and, from
# DESeq2, normalized_counts.csv / vst_counts.csv -- both transformed, and both covering only the
# genes that passed the presence pre-filter. Nothing published a raw gene-level count matrix, so
# the gene FASTA shipped without the counts that belong to it and anyone wanting them had to
# re-derive the aggregation by hand. The pipeline does not quantify genes directly: it quantifies
# transcripts and sums them with tximport, which is exactly the step this writes out.
#
# IT IS NOT A RE-QUANTIFICATION. It makes the same tximport call deseq2.R makes, on the same
# quantifications, and -- when the matching dds.rds is available -- refuses to write unless
# round(counts) is cell-for-cell identical to the count assay DESeq2 actually modelled. The
# published matrix is therefore the DE input itself rather than a reconstruction of it.
#
# Difference from the DE input: no filter. Every gene in the map gets a row, including the ones the
# presence pre-filter later drops, because a count matrix that silently omits rows is a trap.
#
# Two ways to call it:
#   staged (the pipeline)  --quant-dir quants --tx2gene a.tsv,b.tsv --gene-fasta ... --dge-dir ...
#   published results dir  --results <outdir>     paths derived from the published layout
suppressMessages(library(tximport))

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
splitopt <- function(flag) {
  v <- getopt(flag)
  if (is.null(v) || !nzchar(v)) return(character(0))
  v <- unlist(strsplit(v, "[,[:space:]]+"))
  v[nzchar(v)]
}

res         <- getopt("--results")
samplesheet <- getopt("--samplesheet")
primary     <- getopt("--primary-map", "locus")
outdir      <- getopt("--outdir")
if (is.null(samplesheet)) stop("need --samplesheet")

# the gene maps, keyed by the suffix their tx2gene file carries
LABEL <- c(locus  = "locus rebuild (genomic coordinates)",
           corset = "Corset (shared reads + expression)",
           string = "string-stripped transcript ids")
map_key <- function(set) {
  ifelse(grepl("_locus$", set), "locus", ifelse(grepl("_corset$", set), "corset", "string"))
}

if (!is.null(res)) {
  # ---- published layout -------------------------------------------------------------------------
  TX    <- file.path(res, "transcriptome")
  QUANT <- file.path(res, "quantification")
  if (is.null(outdir)) outdir <- file.path(QUANT, "gene_level")
  t2g_files <- Sys.glob(file.path(TX, "tx2gene.*.tsv"))
  # only the maps of the analysed set: the ladder's earlier rungs have tx2gene files too, and they
  # describe references these quantifications were never made against
  t2g_files <- t2g_files[grepl("_expressed(_locus|_corset)?\\.tsv$", basename(t2g_files))]
  fa_files  <- Sys.glob(file.path(TX, "*.genes.fasta"))
  dge_files <- Sys.glob(file.path(res, "dge*", "dge_*", "dds.rds"))
  tx_fasta  <- Sys.glob(file.path(TX, "*_expressed.fasta"))
} else {
  # ---- staged layout (nextflow work dir, everything flat) ---------------------------------------
  QUANT <- getopt("--quant-dir", "quants")
  if (is.null(outdir)) outdir <- "gene_level"
  t2g_files <- splitopt("--tx2gene")
  fa_files  <- splitopt("--gene-fasta")
  dge_files <- file.path(splitopt("--dge-dir"), "dds.rds")
  dge_files <- dge_files[file.exists(dge_files)]
  tx_fasta  <- splitopt("--transcript-fasta")
}
if (!length(t2g_files)) stop("no tx2gene files given")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sets <- sub("\\.tsv$", "", sub("^tx2gene\\.", "", basename(t2g_files)))
# primary first, so the headline matrix is the first thing written and the first thing logged
ord  <- order(map_key(sets) != primary, sets)
t2g_files <- t2g_files[ord]
sets      <- sets[ord]

# ---- sample order: built the same way deseq2.R builds it, so the columns line up ---------------
sheet <- read.csv(samplesheet, stringsAsFactors = FALSE)
meta  <- unique(sheet[, setdiff(colnames(sheet), c("fastq_1", "fastq_2")), drop = FALSE])
meta  <- meta[!duplicated(meta$sample), , drop = FALSE]
meta  <- meta[order(meta$sample), , drop = FALSE]
files <- file.path(QUANT, meta$sample, "quant.sf")
names(files) <- meta$sample
if (any(!file.exists(files)))
  stop("missing quantifications: ", paste(names(files)[!file.exists(files)], collapse = ", "))
cat(sprintf("[counts] %d samples: %s\n", length(files), paste(names(files), collapse = ", ")))

fasta_ids <- function(path) {
  if (is.null(path) || !length(path) || !file.exists(path[1])) return(NULL)
  h <- grep("^>", readLines(path[1], warn = FALSE), value = TRUE)
  sub("^>([^ \t]+).*$", "\\1", h)
}

write_mat <- function(m, path, id_col, digits = 4) {
  df <- data.frame(round(m, digits), check.names = FALSE)
  df <- cbind(setNames(data.frame(rownames(m), stringsAsFactors = FALSE), id_col), df)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[counts]   %-60s %d x %d\n", basename(path), nrow(m), ncol(m)))
}

REPORT <- c("Count matrices: provenance and checks", strrep("=", 64), "",
            sprintf("generated      : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            sprintf("R              : %s", R.version.string),
            sprintf("tximport       : %s", as.character(packageVersion("tximport"))),
            sprintf("samples        : %d (%s)", length(files),
                    paste(names(files), collapse = ", ")),
            sprintf("primary map    : %s", primary), "")

# dropInfReps: the bootstrap replicates are a separate element of the tximport object and do not
# touch counts/abundance/length, so skipping them changes nothing here and saves reading the
# bootstrap matrices once per map.
TXI <- function(...) tximport(files, type = "salmon", dropInfReps = TRUE, ...)

# ---- transcript level, once ---------------------------------------------------------------------
# Published for the same reason as the gene matrices: quant.sf exists per sample, but no matrix.
# This is the input side of the aggregation, so a reader can check the sum themselves.
txi_tx <- TXI(txOut = TRUE)
write_mat(txi_tx$counts,    file.path(outdir, "transcript_counts.tsv"),           "transcript")
write_mat(txi_tx$abundance, file.path(outdir, "transcript_tpm.tsv"),              "transcript")
write_mat(txi_tx$length,    file.path(outdir, "transcript_effective_length.tsv"), "transcript")
tx_fa <- fasta_ids(tx_fasta)
REPORT <- c(REPORT, "transcript level", strrep("-", 16),
            sprintf("  rows %d", nrow(txi_tx$counts)),
            if (!is.null(tx_fa))
              sprintf("  CHECK %s rows == %s (%d records)",
                      if (setequal(rownames(txi_tx$counts), tx_fa)) "ok:" else "FAILED:",
                      basename(tx_fasta[1]), length(tx_fa)),
            "")
if (!is.null(tx_fa) && !setequal(rownames(txi_tx$counts), tx_fa))
  stop("transcript matrix does not match the reference FASTA")

# ---- gene level, per map ------------------------------------------------------------------------
for (i in seq_along(sets)) {
  set <- sets[i]
  key <- map_key(set)
  lab <- paste0(LABEL[[key]], if (key == primary) "  -- PRIMARY" else "")
  cat(sprintf("\n[counts] %s -- %s\n", key, lab))

  t2g <- read.table(t2g_files[i], sep = "\t", col.names = c("tx", "gene"),
                    stringsAsFactors = FALSE)
  txi <- TXI(tx2gene = t2g, ignoreTxVersion = FALSE)

  lines <- c(sprintf("gene level -- %s  [%s]", key, lab), strrep("-", 16),
             sprintf("  tx2gene   %s (%d transcripts -> %d genes)", basename(t2g_files[i]),
                     nrow(t2g), length(unique(t2g$gene))),
             sprintf("  rows      %d", nrow(txi$counts)))

  # the FASTA is the deliverable these counts belong to: one row per record, no more, no less
  gfa_f <- fa_files[basename(fa_files) == sprintf("%s.genes.fasta", set)]
  gfa   <- fasta_ids(gfa_f)
  if (!is.null(gfa)) {
    ok <- setequal(rownames(txi$counts), gfa)
    lines <- c(lines, sprintf("  CHECK %s rows == %s.genes.fasta (%d records)%s",
                              if (ok) "ok:" else "FAILED:", set, length(gfa),
                              if (ok) "" else
                                sprintf("  [only-in-matrix %d, only-in-fasta %d]",
                                        length(setdiff(rownames(txi$counts), gfa)),
                                        length(setdiff(gfa, rownames(txi$counts))))))
    if (!ok) stop("gene matrix does not match ", set, ".genes.fasta")
  }

  # The load-bearing check: is this the matrix DESeq2 got? DESeqDataSetFromTximport stores
  # round(txi$counts) for the genes that passed the pre-filter, and txi$length as the avgTxLength
  # offset, so every cell of both must agree. A mismatch means the published matrix describes a
  # different quantification from the published DE tables, which is the one failure here that
  # would be invisible downstream -- so it stops the run rather than warning.
  dds_f <- dge_files[basename(dirname(dge_files)) == sprintf("dge_%s", set)]
  if (length(dds_f)) {
    suppressMessages(library(SummarizedExperiment))
    dds <- readRDS(dds_f[1])
    dc  <- SummarizedExperiment::assay(dds, "counts")
    nbad <- sum(round(txi$counts[rownames(dc), colnames(dc), drop = FALSE]) != dc)
    lines <- c(lines, sprintf(
      "  CHECK %s round(counts) identical to dge_%s/dds.rds over its %d genes x %d samples%s",
      if (nbad == 0) "ok:" else "FAILED:", set, nrow(dc), ncol(dc),
      if (nbad == 0) "" else sprintf("  [%d cells differ]", nbad)))
    if (nbad != 0) stop("counts do not match the DE input for map ", key)
    al <- tryCatch(SummarizedExperiment::assay(dds, "avgTxLength"), error = function(e) NULL)
    if (!is.null(al)) {
      d <- max(abs(txi$length[rownames(al), colnames(al), drop = FALSE] - al))
      lines <- c(lines, sprintf(
        "  CHECK %s length matrix == the avgTxLength offset DESeq2 used (max abs diff %.3g)",
        if (d < 1e-6) "ok:" else "FAILED:", d))
      if (d >= 1e-6) stop("length matrix does not match the DESeq2 offset for map ", key)
    }
  } else {
    lines <- c(lines, sprintf("  note: no dge_%s/dds.rds given, matrix unverified against a DE run",
                              set))
  }

  sfx <- sprintf(".%s.tsv", set)
  write_mat(txi$counts,    file.path(outdir, paste0("gene_counts", sfx)), "gene")
  write_mat(txi$abundance, file.path(outdir, paste0("gene_tpm", sfx)),    "gene")
  write_mat(txi$length,    file.path(outdir, paste0("gene_length", sfx)), "gene")

  # A second counts matrix that carries the length correction INSIDE the numbers. DESeq2 got the
  # raw counts plus the length matrix as an offset; a tool that cannot take an offset (edgeR as
  # usually driven, or anyone loading the TSV by hand) would silently ignore gene_length and
  # compare genes whose effective length moved between samples. lengthScaledTPM is for that reader.
  txi_ls <- TXI(tx2gene = t2g, ignoreTxVersion = FALSE, countsFromAbundance = "lengthScaledTPM")
  write_mat(txi_ls$counts, file.path(outdir, paste0("gene_counts_lengthScaledTPM", sfx)), "gene")

  raw_tot <- colSums(txi$counts)
  ls_tot  <- colSums(txi_ls$counts)
  lines <- c(lines,
             sprintf("  library sizes (raw counts)  min %s  max %s",
                     format(round(min(raw_tot)), big.mark = ","),
                     format(round(max(raw_tot)), big.mark = ",")),
             sprintf("  genes with zero counts in every sample: %d", sum(rowSums(txi$counts) == 0)),
             sprintf("  CHECK %s lengthScaledTPM preserves every library size",
                     if (max(abs(ls_tot - raw_tot) / raw_tot) < 1e-6) "ok:" else "FAILED:"),
             "")
  REPORT <- c(REPORT, lines)
}

README <- c(
"Gene-level and transcript-level count matrices",
strrep("=", 64),
"",
"WHAT THESE ARE",
"  gene_counts.<set>.tsv",
"      Estimated reads per gene per sample: salmon's transcript NumReads summed over the",
"      transcripts of each gene by tximport. Fractional, because salmon's estimates are.",
"      THIS IS THE MATRIX DESeq2 WAS GIVEN -- the pipeline does not quantify genes directly, it",
"      quantifies transcripts and aggregates. Rounded to integers it is identical to the count",
"      assay of the matching dge_*/dds.rds; checks.txt records that comparison, and the step",
"      fails rather than publishing if it does not hold.",
"",
"  gene_length.<set>.tsv",
"      Average transcript length per gene per sample, abundance-weighted. DESeq2 takes this as a",
"      per-cell offset (assay 'avgTxLength'), which is how a gene whose isoform mix shifts",
"      between conditions stays comparable. It is not decoration: using gene_counts without it",
"      is a different analysis from the published one.",
"",
"  gene_counts_lengthScaledTPM.<set>.tsv",
"      The same counts with that length correction folded into the numbers, and library sizes",
"      preserved. Use this one in any tool that cannot take a length offset, rather than",
"      ignoring gene_length.",
"",
"  gene_tpm.<set>.tsv",
"      Transcript TPM summed per gene. Within-sample composition; not for differential testing.",
"",
"  transcript_{counts,tpm,effective_length}.tsv",
"      The same three quantities at transcript level -- the input side of the aggregation, i.e.",
"      the per-sample quant.sf files assembled into matrices.",
"",
"WHICH MAP",
"  The gene maps are different answers to 'what is a gene' over the SAME quantifications, so they",
"  are published side by side and must not be averaged or mixed. The <set> in each filename names",
"  the map; the one marked PRIMARY in checks.txt is the headline map for this run (and the map the",
"  published gene FASTA and PEP use).",
"",
"  Do not mistake quantification/<sample>/quant.genes.sf for these. Salmon writes those with the",
"  --geneMap it was handed at quantification time, which is the string-stripped map -- the locus",
"  rebuild does not exist yet at that point in the pipeline. They are a different gene set and no",
"  reported result uses them.",
"",
"UNFILTERED",
"  Every gene in the map has a row, including genes the DE presence pre-filter later drops",
"  (de_min_count in at least de_min_samples libraries). The DE tables therefore cover fewer genes",
"  than these matrices; that is the filter, not a discrepancy.",
"",
"PROVENANCE",
"  bin/count_matrix.R, from the COUNT_MATRIX process, on the same quantifications and in the same",
"  R environment as the DESEQ2 process. See checks.txt.")

writeLines(README, file.path(outdir, "README.txt"))
writeLines(REPORT, file.path(outdir, "checks.txt"))
cat(sprintf("\n[counts] wrote %s\n", outdir))
cat(paste(grep("CHECK", REPORT, value = TRUE), collapse = "\n"), "\n")
