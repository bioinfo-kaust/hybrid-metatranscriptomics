#!/usr/bin/env Rscript
# Publication figure set (PDF + 300 dpi PNG). Every panel is computed from this run's own
# outputs; a panel whose input is missing is skipped with a message rather than failing the run.
suppressMessages({library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
dge_dir  <- getopt("--dge-dir")
enr_dir  <- getopt("--enrichment-dir")
ann_file <- getopt("--annotation")
qc_dir   <- getopt("--qc-dir")
busco_d  <- getopt("--busco-dir")
ladder_f <- getopt("--ladder")
cond_col <- getopt("--condition-column", "condition")
fig      <- getopt("--outdir", "figures")
# same thresholds the tables were made with — passed in, never re-declared here
padj_cut <- as.numeric(getopt("--padj", "0.05"))
lfc_cut  <- as.numeric(getopt("--lfc", "1"))
sig_lab  <- sprintf("padj < %g, |log2FC| > %g", padj_cut, lfc_cut)
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
save2 <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
ok <- function(expr, nm) tryCatch({expr; cat("  [ok]", nm, "\n")},
                                  error = function(e) cat("  [SKIP]", nm, ":", conditionMessage(e), "\n"))
qc <- function(f) { p <- file.path(qc_dir, f); if (!file.exists(p)) stop("missing ", f); read.delim(p) }

# F1 — rRNA depletion per sample
ok({
  rr <- qc("rrna_per_sample.tsv"); rr$sample <- factor(rr$sample, levels = rr$sample)
  p <- ggplot(rr, aes(sample, rrna_pct, fill = .data[[cond_col]])) + geom_col() + th +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       labs(title = "rRNA content per sample", y = "% reads rRNA", x = NULL, fill = cond_col)
  save2(p, "F1_rRNA_per_sample", 8, 5)
}, "F1 rRNA")

# F2 — contaminant read removal per sample
ok({
  sr <- qc("decontam_per_sample.tsv")
  d <- data.frame(sample = rep(sr$sample, 2),
                  fraction = rep(c("Host (kept)", "Contaminant (removed)"), each = nrow(sr)),
                  pairs = c(sr$pairs_kept, sr$contaminant_removed))
  d$sample <- factor(d$sample, levels = sr$sample)
  d$fraction <- factor(d$fraction, levels = c("Host (kept)", "Contaminant (removed)"))
  p <- ggplot(d, aes(sample, pairs / 1e6, fill = fraction)) + geom_col() +
       scale_fill_manual(values = c(`Host (kept)` = "#127A8A", `Contaminant (removed)` = "#CC6611")) +
       th + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       labs(title = "Contaminant read removal (competitive placement)",
            y = "read pairs (millions)", x = NULL, fill = NULL)
  save2(p, "F2_contaminant_removal", 8.5, 5)
}, "F2 decontamination")

# F3 — BUSCO completeness across the reference ladder
ok({
  files <- list.files(busco_d, pattern = "\\.busco\\.txt$", full.names = TRUE)
  if (!length(files)) stop("no BUSCO summaries")
  rows <- do.call(rbind, lapply(files, function(f) {
    l <- readLines(f, warn = FALSE)[1]
    g <- regmatches(l, regexec("C:([0-9.]+)%\\[S:([0-9.]+)%,D:([0-9.]+)%\\],F:([0-9.]+)%,M:([0-9.]+)%", l))[[1]]
    if (length(g) < 6) return(NULL)
    data.frame(set = sub("\\.busco\\.txt$", "", basename(f)),
               Single = as.numeric(g[3]), Duplicated = as.numeric(g[4]),
               Fragmented = as.numeric(g[5]), Missing = as.numeric(g[6]))
  }))
  d <- reshape(rows, direction = "long", varying = list(2:5), v.names = "pct",
               timevar = "category", times = c("Single", "Duplicated", "Fragmented", "Missing"),
               idvar = "set")
  d$category <- factor(d$category, levels = c("Single", "Duplicated", "Fragmented", "Missing"))
  p <- ggplot(d, aes(set, pct, fill = category)) + geom_col() + coord_flip() + th +
       scale_fill_manual(values = c(Single = "#3377BB", Duplicated = "#88BBDD",
                                    Fragmented = "#EEBB55", Missing = "#CC5544")) +
       labs(title = "BUSCO completeness by reference variant", x = NULL, y = "% of lineage set", fill = NULL)
  save2(p, "F3_busco_ladder", 8, 5)
}, "F3 BUSCO")

# F4 — PCA of the variance-stabilised counts, with the variance each axis actually carries
ok({
  pc <- read.csv(file.path(dge_dir, "pca_data.csv"))
  grp <- if (cond_col %in% colnames(pc)) cond_col else "group"
  vfile <- file.path(dge_dir, "pca_variance.csv")
  ax <- if (file.exists(vfile)) {
    v <- read.csv(vfile); c(sprintf("PC1 (%.1f%%)", v$percentVar[1]), sprintf("PC2 (%.1f%%)", v$percentVar[2]))
  } else c("PC1", "PC2")
  p <- ggplot(pc, aes(PC1, PC2, colour = .data[[grp]])) + geom_point(size = 3) + th +
       labs(title = "PCA of variance-stabilised counts", x = ax[1], y = ax[2], colour = cond_col)
  save2(p, "F4_PCA", 6.5, 5)
}, "F4 PCA")

# F5 — DE counts per contrast
ok({
  s <- read.csv(file.path(dge_dir, "DE_summary.csv"))
  d <- data.frame(contrast = rep(s$contrast, 2),
                  direction = rep(c("up", "down"), each = nrow(s)),
                  n = c(s$up, s$down))
  p <- ggplot(d, aes(contrast, n, fill = direction)) + geom_col(position = "dodge") + th +
       scale_fill_manual(values = c(up = "#CC5544", down = "#3377BB")) +
       theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
       labs(title = "Differentially expressed genes", y = paste0("genes (", sig_lab, ")"), x = NULL)
  save2(p, "F5_DE_counts", 7, 5)
}, "F5 DE counts")

# F6 — volcano per contrast
ok({
  for (f in list.files(dge_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)) {
    if (grepl("DE_summary", f)) next
    nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
    de <- read.csv(f, row.names = 1)
    de <- de[!is.na(de$padj) & !is.na(de$log2FoldChange), ]
    de$sig <- de$padj < padj_cut & abs(de$log2FoldChange) > lfc_cut
    xlab <- if ("log2FoldChange_unshrunken" %in% colnames(de)) "log2 fold change (shrunken)"
            else "log2 fold change (unshrunken)"
    p <- ggplot(de, aes(log2FoldChange, -log10(padj), colour = sig)) +
         geom_point(size = 0.6, alpha = 0.5) + th +
         scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#CC5544"), guide = "none") +
         labs(title = paste("Volcano:", nm), subtitle = sig_lab, x = xlab, y = "-log10 padj")
    save2(p, paste0("F6_volcano_", nm), 6, 5)
  }
}, "F6 volcanoes")

# F7 — annotation evidence
ok({
  a <- read.delim(ann_file, stringsAsFactors = FALSE)
  # count a column only if it exists: a missing column would otherwise plot a silent zero
  n_of <- function(col, pred = function(x) x != "") {
    if (!col %in% colnames(a)) { cat("    [note] no column", col, "\n"); return(NA_integer_) }
    sum(pred(a[[col]]), na.rm = TRUE)
  }
  # GO comes from two sources now; count a transcript once if either carries a term
  go_any <- if (all(c("GOs", "InterPro_GOs") %in% colnames(a))) {
    sum(a$GOs != "" | a$InterPro_GOs != "", na.rm = TRUE)
  } else n_of("GOs")
  ev <- data.frame(
    evidence = c("ORF (coding)", "protein hit", "InterPro", "GO (eggNOG or InterPro)", "KEGG pathway"),
    n = c(n_of("coding_status", function(x) x == "coding"), n_of("protein_hit"),
          n_of("InterPro"), go_any, n_of("KEGG_Pathway")))
  ev <- ev[!is.na(ev$n), ]
  p <- ggplot(ev, aes(reorder(evidence, n), n)) + geom_col(fill = "#127A8A") + coord_flip() + th +
       labs(title = "Annotation evidence per transcript", x = NULL, y = "transcripts")
  save2(p, "F7_annotation_evidence", 6.5, 4)
}, "F7 annotation")

# F8 — the reference ladder as a table figure
ok({
  lad <- read.delim(ladder_f, stringsAsFactors = FALSE)
  lad$n_transcripts <- suppressWarnings(as.numeric(gsub(",", "", lad$n_transcripts)))
  p <- ggplot(lad, aes(reorder(reference, n_transcripts), n_transcripts)) +
       geom_col(fill = "#3377BB") + coord_flip() + th +
       labs(title = "Reference ladder: transcripts per variant", x = NULL, y = "transcripts")
  save2(p, "F8_reference_ladder", 7.5, 4)
}, "F8 ladder")

# F9 - strongest enriched terms per contrast (whatever enrichment produced)
ok({
  files <- list.files(enr_dir, pattern = "_(GO|KEGG)_(ORA|GSEA)\\.csv$", full.names = TRUE)
  if (!length(files)) stop("no enrichment tables")
  d <- do.call(rbind, lapply(files, function(f) {
    x <- read.csv(f, stringsAsFactors = FALSE)
    if (!nrow(x) || !"p.adjust" %in% colnames(x)) return(NULL)   # read p.adjust BY NAME
    x <- x[order(x$p.adjust), ]
    x <- head(x[x$p.adjust < 0.05, ], 10)
    if (!nrow(x)) return(NULL)
    lab <- if ("Description" %in% colnames(x) && any(x$Description != "")) x$Description else x$ID
    data.frame(analysis = sub("\\.csv$", "", basename(f)), term = lab, padj = x$p.adjust)
  }))
  if (is.null(d) || !nrow(d)) stop("no terms with padj < 0.05")
  d$term <- substr(d$term, 1, 55)
  p <- ggplot(d, aes(reorder(term, -padj), -log10(padj))) +
       geom_col(fill = "#127A8A") + coord_flip() + facet_wrap(~analysis, scales = "free_y") + th +
       labs(title = sprintf("Top enriched terms (padj < %g)", padj_cut), x = NULL, y = "-log10 padj")
  save2(p, "F9_enrichment", 11, 2 + 0.35 * nrow(d))
}, "F9 enrichment")

# ---------------------------------------------------------------------------------------------
# Diagnostics. These decide whether the tables above can be believed, so they are figures in the
# set rather than numbers in a log nobody opens.

# F10 — sample-to-sample distance heatmap: does the treatment group samples at all?
ok({
  m <- as.matrix(read.csv(file.path(dge_dir, "sample_distances.csv"), row.names = 1))
  d <- data.frame(a = rep(rownames(m), ncol(m)),
                  b = rep(colnames(m), each = nrow(m)),
                  dist = as.vector(m))
  d$a <- factor(d$a, levels = rownames(m)); d$b <- factor(d$b, levels = rownames(m))
  p <- ggplot(d, aes(a, b, fill = dist)) + geom_tile() + th +
       scale_fill_gradient(low = "#127A8A", high = "grey92") +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       labs(title = "Sample-to-sample distance (VST, Euclidean)", x = NULL, y = NULL, fill = "distance")
  save2(p, "F10_sample_distances", 7.5, 6)
}, "F10 sample distances")

# F11 — dispersion fit: a fit that misses the cloud invalidates every p-value downstream
ok({
  d <- read.csv(file.path(dge_dir, "dispersion_data.csv"))
  d <- d[d$baseMean > 0, ]
  p <- ggplot(d, aes(baseMean)) +
       geom_point(aes(y = dispGeneEst), size = 0.3, alpha = 0.25, colour = "grey55") +
       geom_point(aes(y = dispersion), size = 0.3, alpha = 0.35, colour = "#127A8A") +
       geom_line(aes(y = dispFit), colour = "#CC5544", linewidth = 0.6) +
       scale_x_log10() + scale_y_log10() + th +
       labs(title = "Dispersion estimates (grey: per-gene, teal: final, red: fit)",
            x = "mean of normalised counts", y = "dispersion")
  save2(p, "F11_dispersion", 6.5, 5)
}, "F11 dispersion")

# F12 — p-value histograms: a flat bulk with a spike near zero is what a healthy test looks like
ok({
  d <- do.call(rbind, lapply(list.files(dge_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE), function(f) {
    if (grepl("DE_summary", f)) return(NULL)
    x <- read.csv(f, row.names = 1)
    x <- x[!is.na(x$pvalue), ]
    if (!nrow(x)) return(NULL)
    data.frame(contrast = sub("^DE_", "", sub("\\.csv$", "", basename(f))), pvalue = x$pvalue)
  }))
  if (is.null(d) || !nrow(d)) stop("no p-values")
  p <- ggplot(d, aes(pvalue)) + geom_histogram(binwidth = 0.025, fill = "#3377BB", colour = "white") +
       facet_wrap(~contrast, scales = "free_y") + th +
       labs(title = "p-value distribution per contrast", x = "raw p-value", y = "genes")
  save2(p, "F12_pvalue_histograms", 9, 3.5)
}, "F12 p-value histograms")

# F13 — MA per contrast: where the significant genes sit in expression
ok({
  for (f in list.files(dge_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)) {
    if (grepl("DE_summary", f)) next
    nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
    de <- read.csv(f, row.names = 1)
    de <- de[!is.na(de$log2FoldChange) & de$baseMean > 0, ]
    de$sig <- !is.na(de$padj) & de$padj < padj_cut & abs(de$log2FoldChange) > lfc_cut
    p <- ggplot(de, aes(baseMean, log2FoldChange, colour = sig)) +
         geom_point(size = 0.4, alpha = 0.4) + scale_x_log10() + th +
         geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
         scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#CC5544"), guide = "none") +
         labs(title = paste("MA:", nm), subtitle = sig_lab,
              x = "mean of normalised counts", y = "log2 fold change")
    save2(p, paste0("F13_MA_", nm), 6, 5)
  }
}, "F13 MA plots")

cat("[figures] written to", fig, "\n")
