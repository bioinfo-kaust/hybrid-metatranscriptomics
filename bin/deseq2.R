#!/usr/bin/env Rscript
# Gene-level differential expression from salmon quantifications.
#
# The design formula, the reference level and the contrast list are all arguments, so this
# serves any factorial design — not just the three-treatment UV experiment it was written for.
#
# log2FoldChange is SHRUNKEN (apeglm where the contrast is a model coefficient, ashr otherwise);
# the raw MLE estimate is kept beside it as log2FoldChange_unshrunken, and `stat` is the Wald
# statistic from the unshrunken fit — that is what the enrichment ranks on.
suppressMessages({library(tximport); library(DESeq2)})

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
quant_dir   <- getopt("--quant-dir")
tx2gene_f   <- getopt("--tx2gene")
samplesheet <- getopt("--samplesheet")
cond_col    <- getopt("--condition-column", "condition")
ref_level   <- getopt("--reference-level", NULL)
design_str  <- getopt("--design", paste("~", cond_col))
contrast_f  <- getopt("--contrasts", NULL)
outdir      <- getopt("--outdir", "dge")
min_count   <- as.numeric(getopt("--min-count", "5"))
min_samples <- as.numeric(getopt("--min-samples", "3"))
alpha       <- as.numeric(getopt("--alpha", "0.05"))
lfc_cut     <- as.numeric(getopt("--lfc", "1"))
shrink_type <- getopt("--shrink", "apeglm")
expr_tx_f   <- getopt("--expressed-tx", NULL)
if (!is.null(expr_tx_f) && (expr_tx_f == "NO_FILE" || !file.exists(expr_tx_f))) expr_tx_f <- NULL
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[deseq2] quants=%s\n[deseq2] tx2gene=%s\n[deseq2] design=%s\n[deseq2] out=%s\n",
            quant_dir, tx2gene_f, design_str, outdir))

tx2gene <- read.table(tx2gene_f, sep = "\t", col.names = c("tx", "gene"), stringsAsFactors = FALSE)

# sample metadata: one row per sample (the sheet may carry one row per unit/lane)
sheet <- read.csv(samplesheet, stringsAsFactors = FALSE)
stopifnot(all(c("sample", cond_col) %in% colnames(sheet)))
meta <- unique(sheet[, setdiff(colnames(sheet), c("fastq_1", "fastq_2")), drop = FALSE])
meta <- meta[!duplicated(meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample
meta <- meta[order(meta$sample), , drop = FALSE]

files <- file.path(quant_dir, meta$sample, "quant.sf")
names(files) <- meta$sample
missing <- files[!file.exists(files)]
if (length(missing)) stop("missing quantifications: ", paste(names(missing), collapse = ", "))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)

meta[[cond_col]] <- factor(meta[[cond_col]])
if (!is.null(ref_level)) {
  if (!ref_level %in% levels(meta[[cond_col]]))
    stop("reference level '", ref_level, "' not among: ", paste(levels(meta[[cond_col]]), collapse = ", "))
  meta[[cond_col]] <- relevel(meta[[cond_col]], ref = ref_level)
}
for (cl in setdiff(colnames(meta), c("sample", cond_col))) {
  if (is.character(meta[[cl]])) meta[[cl]] <- factor(meta[[cl]])
}

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = as.formula(design_str))

# ---- DE-eligible genes ---------------------------------------------------------------
# The quantification reference is whatever the expression filter kept. With expr_filter_rescue
# off (the default) that is exactly the expressed set and this restriction is a no-op. With it on,
# the kept set is a UNION — transcripts the experiment saw above threshold, plus unexpressed ones
# rescued on ORF + homology — and the two parts are not equivalent: a gene made only of transcripts
# the experiment never saw above threshold carries no evidence about treatment and costs
# multiple-testing burden for nothing, even though its model belongs in the index.
#
# So the restriction is applied HERE, to genes, and AFTER tximport rather than by cutting the
# tx2gene map: a gene that qualifies keeps the counts of all its transcripts, rescued isoforms
# included, so no reads are discarded. A gene qualifies if at least one of its transcripts is in
# the expressed list.
#
# This is also what sets the ORA background: enrichment.R takes its universe from the genes in
# these results tables, so restricting the model restricts the background in the same step.
if (!is.null(expr_tx_f)) {
  expr_tx <- readLines(expr_tx_f)
  expr_tx <- expr_tx[nzchar(expr_tx)]
  eligible <- unique(tx2gene$gene[tx2gene$tx %in% expr_tx])
  n_all <- nrow(dds)
  hit <- rownames(dds) %in% eligible
  if (!any(hit)) stop("--expressed-tx matched no gene in the count matrix: ",
                      "wrong transcript ids for this gene map?")
  dds <- dds[hit, ]
  cat(sprintf(paste0("[deseq2] DE-eligible: %d expressed transcripts -> %d of %d genes kept ",
                     "(%d dropped as unexpressed-but-annotated)\n"),
              length(expr_tx), nrow(dds), n_all, n_all - nrow(dds)))
} else {
  cat("[deseq2] no --expressed-tx: every gene in the reference is DE-eligible\n")
}

# Presence filter, not a total: a sum across all libraries passes a gene carried by one sample,
# which on a redundant reference is usually the EM splitting reads between near-duplicates.
n_before <- nrow(dds)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
cat(sprintf("[deseq2] filter: >= %g counts in >= %g samples -> %d of %d genes kept\n",
            min_count, min_samples, nrow(dds), n_before))

dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "dds.rds"))

write.csv(counts(dds, normalized = TRUE), file.path(outdir, "normalized_counts.csv"))
vsd <- vst(dds, blind = TRUE)
write.csv(assay(vsd), file.path(outdir, "vst_counts.csv"))
pc <- plotPCA(vsd, intgroup = cond_col, returnData = TRUE)
write.csv(pc, file.path(outdir, "pca_data.csv"), row.names = FALSE)
# the percentages live in an attribute, so they never reached the figure without this
write.csv(data.frame(PC = c("PC1", "PC2"), percentVar = round(100 * attr(pc, "percentVar"), 1)),
          file.path(outdir, "pca_variance.csv"), row.names = FALSE)

# sample-to-sample distances — the honest check on whether the factor is the dominant axis
sd_mat <- as.matrix(dist(t(assay(vsd))))
write.csv(sd_mat, file.path(outdir, "sample_distances.csv"))

# dispersion fit, for the diagnostic panel: a fit that misses the cloud invalidates every p-value
disp <- data.frame(baseMean = mcols(dds)$baseMean, dispGeneEst = mcols(dds)$dispGeneEst,
                   dispFit = mcols(dds)$dispFit, dispersion = dispersions(dds))
write.csv(disp[!is.na(disp$baseMean) & disp$baseMean > 0, ],
          file.path(outdir, "dispersion_data.csv"), row.names = FALSE)

# ---- contrasts ---------------------------------------------------------------------
if (!is.null(contrast_f) && file.exists(contrast_f)) {
  ct <- read.csv(contrast_f, stringsAsFactors = FALSE)
  stopifnot(all(c("id", "factor", "level", "reference") %in% colnames(ct)))
} else {
  # default: every level of the condition column against the reference level
  lv <- levels(meta[[cond_col]])
  ct <- data.frame(id = paste0(lv[-1], "_vs_", lv[1]), factor = cond_col,
                   level = lv[-1], reference = lv[1], stringsAsFactors = FALSE)
  cat("[deseq2] no contrast file: testing every level against", lv[1], "\n")
}

# apeglm shrinks a model COEFFICIENT, so out of the box it only serves contrasts against the
# reference level. Mixing estimators across contrasts would make the |log2FC| cut mean different
# things in different columns of the same summary table, so for a contrast that is not a
# coefficient the factor is releveled and only the Wald test is refitted — dispersions do not
# depend on the reference level, so this costs seconds, not a second DESeq() run. ashr is the
# fallback if that fails. The MLE estimate is preserved in its own column either way.
shrink <- function(res, f, lvl, ref) {
  if (identical(shrink_type, "none")) return(list(res = res, how = "none"))
  coefname <- paste0(f, "_", lvl, "_vs_", ref)
  if (shrink_type == "apeglm") {
    d <- dds
    if (!coefname %in% resultsNames(d)) {
      d <- tryCatch({
        colData(d)[[f]] <- relevel(factor(colData(d)[[f]]), ref = ref)
        nbinomWaldTest(d, quiet = TRUE)
      }, error = function(e) NULL)
    }
    if (!is.null(d) && coefname %in% resultsNames(d)) {
      out <- tryCatch(lfcShrink(d, coef = coefname, type = "apeglm", res = res, quiet = TRUE),
                      error = function(e) NULL)
      if (!is.null(out)) return(list(res = out, how = "apeglm"))
    }
  }
  fallback <- if (shrink_type %in% c("ashr", "normal")) shrink_type
              else if (requireNamespace("ashr", quietly = TRUE)) "ashr" else "normal"
  out <- tryCatch(lfcShrink(dds, contrast = c(f, lvl, ref), type = fallback, res = res, quiet = TRUE),
                  error = function(e) NULL)
  if (is.null(out)) return(list(res = res, how = "none (shrinkage failed)"))
  list(res = out, how = fallback)
}

summary_tab <- data.frame()
filter_tab  <- data.frame()
for (i in seq_len(nrow(ct))) {
  nm <- ct$id[i]
  raw <- tryCatch(results(dds, contrast = c(ct$factor[i], ct$level[i], ct$reference[i]), alpha = alpha),
                  error = function(e) { cat("[SKIP]", nm, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(raw)) next
  sh  <- shrink(raw, ct$factor[i], ct$level[i], ct$reference[i])
  res <- sh$res
  cat("[deseq2]", nm, "- shrinkage:", sh$how, "\n")

  out <- as.data.frame(res)
  out$log2FoldChange_unshrunken <- raw$log2FoldChange[match(rownames(out), rownames(raw))]
  if (!"stat" %in% colnames(out)) out$stat <- raw$stat[match(rownames(out), rownames(raw))]
  out <- out[order(out$padj), ]
  write.csv(out, file.path(outdir, paste0("DE_", nm, ".csv")))

  up <- sum(out$padj < alpha & out$log2FoldChange >  lfc_cut, na.rm = TRUE)
  dn <- sum(out$padj < alpha & out$log2FoldChange < -lfc_cut, na.rm = TRUE)
  summary_tab <- rbind(summary_tab, data.frame(
    contrast = nm, tested = sum(!is.na(out$padj)),
    sig_padj = sum(out$padj < alpha, na.rm = TRUE),   # no fold-change filter (kept for continuity)
    sig      = up + dn,                                # padj AND |log2FC| — the headline number
    up = up, down = dn))

  # what the defaults threw away, per contrast: independent filtering trims low-count genes at a
  # threshold chosen per contrast, and Cook's distance blanks outlier genes outright (with n=5 it
  # cannot replace them, so they simply leave the test).
  thr <- tryCatch(as.numeric(metadata(raw)$filterThreshold), error = function(e) NA_real_)
  filter_tab <- rbind(filter_tab, data.frame(
    contrast = nm, genes_in_model = nrow(raw), tested = sum(!is.na(raw$padj)),
    filter_threshold_basemean = if (length(thr)) round(thr[1], 3) else NA_real_,
    dropped_independent_filter = sum(!is.na(raw$pvalue) & is.na(raw$padj)),
    dropped_cooks_outlier = sum(is.na(raw$pvalue) & raw$baseMean > 0),
    shrinkage = sh$how))
}
write.csv(summary_tab, file.path(outdir, "DE_summary.csv"), row.names = FALSE)
# NOT named DE_*: every consumer globs DE_*.csv for contrast tables
write.csv(filter_tab, file.path(outdir, "filtering_summary.csv"), row.names = FALSE)
cat(sprintf("=== DE summary (sig = padj < %g AND |log2FC| > %g; sig_padj drops the fold-change cut) ===\n",
            alpha, lfc_cut))
print(summary_tab)
cat("=== what the filters removed ===\n"); print(filter_tab)
cat("Samples per level:\n"); print(table(meta[[cond_col]]))
