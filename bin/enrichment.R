#!/usr/bin/env Rscript
# GO + KEGG over-representation (ORA) and GO GSEA per contrast, using the pipeline's own
# gene -> GO / KO annotation as TERM2GENE. No model-organism OrgDb is required.
#
# Read p.adjust BY NAME in any downstream check: core_enrichment contains quoted commas, so an
# awk test on the last field silently reports every row as significant.
suppressMessages({library(clusterProfiler)})

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
ann_file <- getopt("--annotation")
dge_dir  <- getopt("--dge-dir")
outdir   <- getopt("--outdir", "enrichment")
padj_cut <- as.numeric(getopt("--padj", "0.05"))
lfc_cut  <- as.numeric(getopt("--lfc", "1"))
seed     <- as.integer(getopt("--seed", "42"))
# fgsea's multilevel algorithm is randomised: without a fixed seed the same inputs give a
# different set of significant terms on every run (measured on this project: 7 vs 13 vs 19
# significant GO terms for one contrast across three runs). Seed it, and pass seed = TRUE so
# clusterProfiler seeds fgsea's own RNG as well.
set.seed(seed)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[enrich] ann=%s\n[enrich] dge=%s\n[enrich] out=%s\n", ann_file, dge_dir, outdir))

ann <- read.delim(ann_file, stringsAsFactors = FALSE)

expand_map <- function(genecol, termcol, sep = ",") {
  if (!termcol %in% colnames(ann)) return(NULL)
  d <- ann[ann[[termcol]] != "" & !is.na(ann[[termcol]]), c(genecol, termcol)]
  if (!nrow(d)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    ts <- trimws(strsplit(d[[termcol]][i], sep)[[1]])
    ts <- ts[ts != "" & ts != "-"]
    if (length(ts)) data.frame(term = ts, gene = d[[genecol]][i], stringsAsFactors = FALSE) else NULL
  }))
}

# GO from BOTH sources. eggNOG transfers GO by orthology, InterProScan infers it from matched
# signatures (it already runs with -goterms); on a non-model cnidarian each finds terms the other
# misses, and using only eggNOG discarded the second set entirely.
go_egg <- expand_map("gene", "GOs")
go_ipr <- expand_map("gene", "InterPro_GOs")
go2gene <- unique(rbind(go_egg, go_ipr))
cat("[enrich] GO mappings: eggNOG", if (is.null(go_egg)) 0 else nrow(go_egg),
    "+ InterPro", if (is.null(go_ipr)) 0 else nrow(go_ipr),
    "-> union", if (is.null(go2gene)) 0 else nrow(go2gene), "\n")

# KEGG at PATHWAY level. The KEGG_ko column holds orthologue ids, which are genes and not gene
# sets: testing them asks whether a KO is over-represented among genes carrying that KO, which is
# a tautology. eggNOG's KEGG_Pathway column is the actual pathway membership.
ko2gene <- unique(expand_map("gene", "KEGG_Pathway"))
if (!is.null(ko2gene)) {
  ko2gene$term <- sub("^ko", "map", ko2gene$term)     # ko00010 and map00010 are the same pathway
  ko2gene <- unique(ko2gene[grepl("^map[0-9]+$", ko2gene$term), ])
}
cat("[enrich] KEGG pathway mappings:", if (is.null(ko2gene)) 0 else nrow(ko2gene), "\n")

t2name <- NULL
try({
  library(GO.db)
  ids <- unique(go2gene$term)
  nm <- AnnotationDbi::select(GO.db, keys = ids, columns = "TERM", keytype = "GOID")
  t2name <- data.frame(term = nm$GOID, name = nm$TERM, stringsAsFactors = FALSE)
}, silent = TRUE)
if (is.null(t2name)) cat("[enrich] NOTE: GO.db unavailable, terms will carry bare GO ids\n")

# pathway names, so the KEGG table is readable without a browser open. Needs network; a failure
# costs the names, not the analysis.
k2name <- NULL
try({
  kk <- clusterProfiler::download_KEGG("ko")$KEGGPATHID2NAME
  k2name <- data.frame(term = sub("^ko", "map", kk[[1]]), name = kk[[2]], stringsAsFactors = FALSE)
}, silent = TRUE)
if (is.null(k2name)) cat("[enrich] NOTE: KEGG pathway names unavailable, ids only\n")

de_files <- list.files(dge_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("DE_summary", de_files)]
if (!length(de_files)) stop("no DE_*.csv files in ", dge_dir)

# ONE universe for every contrast: every gene that entered the model. Taking it per contrast from
# the genes with a non-NA padj made the background differ between contrasts and excluded exactly
# the low-expression genes that independent filtering removes, which biases ORA in a direction
# nobody had measured.
is_results <- function(f) all(c("padj", "log2FoldChange") %in% colnames(read.csv(f, nrows = 1)))
de_files <- Filter(is_results, de_files)
if (!length(de_files)) stop("no DESeq2 results tables in ", dge_dir)
universe <- unique(unlist(lapply(de_files, function(f) rownames(read.csv(f, row.names = 1)))))
cat("[enrich] universe (genes in the model, shared by all contrasts):", length(universe), "\n")
cat("[enrich] of those, annotated with >=1 GO term:",
    length(intersect(universe, unique(go2gene$gene))), "\n")

# Tables are written in full (pvalueCutoff = 1) and significance is applied once, at --padj, by the
# summary below and by the figures. Previously ORA wrote at p<0.1/q<0.2, GSEA at p<0.25 and the
# summary counted p.adjust<0.05, so three different thresholds were in play at once.
write_res <- function(obj, nm, kind) {
  if (is.null(obj)) return(invisible(NULL))
  d <- as.data.frame(obj)
  write.csv(d, file.path(outdir, paste0(nm, "_", kind, ".csv")), row.names = FALSE)
  cat("   ", kind, ":", nrow(d), "terms returned,",
      if ("p.adjust" %in% colnames(d)) sum(d$p.adjust < padj_cut, na.rm = TRUE) else 0, "significant\n")
}

for (f in de_files) {
  nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de <- read.csv(f, row.names = 1)
  det <- de[!is.na(de$padj), ]
  sig <- rownames(det)[det$padj < padj_cut & abs(det$log2FoldChange) >= lfc_cut]
  cat(nm, "- sig genes:", length(sig), "of", nrow(det), "tested\n")

  if (length(sig) >= 5 && !is.null(go2gene)) {
    write_res(tryCatch(enricher(sig, TERM2GENE = go2gene[, c("term", "gene")], TERM2NAME = t2name,
                                universe = universe, pvalueCutoff = 1, qvalueCutoff = 1),
                       error = function(e) NULL), nm, "GO_ORA")
  }
  if (length(sig) >= 5 && !is.null(ko2gene)) {
    write_res(tryCatch(enricher(sig, TERM2GENE = ko2gene[, c("term", "gene")], TERM2NAME = k2name,
                                universe = universe, pvalueCutoff = 1, qvalueCutoff = 1),
                       error = function(e) NULL), nm, "KEGG_ORA")
  }
  if (!is.null(go2gene)) {
    # rank on the Wald statistic: it carries the evidence, where a fold change alone is dominated
    # by low-count genes with huge ratios. Falls back to the fold change if `stat` is absent.
    rankcol <- if ("stat" %in% colnames(de)) "stat" else "log2FoldChange"
    gl <- de[[rankcol]]; names(gl) <- rownames(de)
    gl <- sort(gl[!is.na(gl)], decreasing = TRUE)
    cat("    GSEA ranked on:", rankcol, "\n")
    set.seed(seed)
    write_res(tryCatch(GSEA(gl, TERM2GENE = go2gene[, c("term", "gene")], TERM2NAME = t2name,
                            pvalueCutoff = 1, minGSSize = 10, seed = TRUE),
                       error = function(e) NULL), nm, "GO_GSEA")
  }
}

# significant-term counts, read by name
out <- do.call(rbind, lapply(list.files(outdir, pattern = "csv$", full.names = TRUE), function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  n <- if ("p.adjust" %in% colnames(d)) sum(d$p.adjust < padj_cut, na.rm = TRUE) else NA
  data.frame(file = basename(f), rows = nrow(d), sig_padj = n)
}))
if (!is.null(out)) {
  write.csv(out, file.path(outdir, "enrichment_summary.csv"), row.names = FALSE)
  print(out)
}
