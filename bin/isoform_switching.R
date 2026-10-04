#!/usr/bin/env Rscript
# Differential isoform usage: structure from the long-read collapsed models, abundances from the
# short-read replicate quantification (IsoformSwitchAnalyzeR + DEXSeq).
#
# CAVEAT: this tests only the isoforms present in the supplied GTF. On a cross-source deduplicated
# reference roughly half the long-read models can be dropped, so a gene's isoform set may be
# incomplete and results at different reference levels are NOT interchangeable.
suppressMessages({library(IsoformSwitchAnalyzeR)})

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
quant_dir   <- getopt("--quant-dir")
iso_fasta   <- getopt("--fasta")
gtf         <- getopt("--gtf")
samplesheet <- getopt("--samplesheet")
cond_col    <- getopt("--condition-column", "condition")
outdir      <- getopt("--outdir", "isoforms")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

salmon <- importIsoformExpression(parentDir = quant_dir, addIsofomIdAsColumn = TRUE)

# keep only the isoforms the GTF describes — the rest have no structure to test
gtf_ids <- unique(sub('.*transcript_id "([^"]+)".*', "\\1",
                      grep("transcript_id", readLines(gtf), value = TRUE)))
keep <- salmon$counts$isoform_id %in% gtf_ids
salmon$counts    <- salmon$counts[keep, ]
salmon$abundance <- salmon$abundance[keep, ]
cat("isoforms with a model in the GTF:", sum(keep), "\n")

sheet <- read.csv(samplesheet, stringsAsFactors = FALSE)
sheet <- sheet[!duplicated(sheet$sample), ]
samples <- colnames(salmon$abundance)[-1]
design <- data.frame(sampleID = samples,
                     condition = sheet[[cond_col]][match(samples, sheet$sample)],
                     stringsAsFactors = FALSE)
stopifnot(!any(is.na(design$condition)))

lv <- unique(design$condition)
cmp <- as.data.frame(t(combn(lv, 2)), stringsAsFactors = FALSE)
colnames(cmp) <- c("condition_1", "condition_2")

sw <- importRdata(
  isoformCountMatrix   = salmon$counts,
  isoformRepExpression = salmon$abundance,
  designMatrix         = design,
  comparisonsToMake    = cmp,
  isoformExonAnnoation = gtf,
  isoformNtFasta       = iso_fasta,
  showProgress = FALSE, ignoreAfterPeriod = FALSE, fixStringTieAnnotationProblem = FALSE)

sw <- preFilter(sw, geneExpressionCutoff = 1, isoformExpressionCutoff = 0, removeSingleIsoformGenes = TRUE)
sw <- isoformSwitchTestDEXSeq(sw, reduceToSwitchingGenes = TRUE, alpha = 0.05, dIFcutoff = 0.1)
sw <- analyzeORF(sw, genomeObject = NULL, orfMethod = "longest", showProgress = FALSE)
sw <- analyzeIntronRetention(sw)

saveRDS(sw, file.path(outdir, "switchlist.rds"))
summ <- extractSwitchSummary(sw)
write.csv(summ, file.path(outdir, "switch_summary.csv"), row.names = FALSE)
cat("=== isoform-switch summary (DTU genes per contrast) ===\n"); print(summ)

cons <- extractTopSwitches(sw, n = Inf, sortByQvals = TRUE, inEachComparison = TRUE)
write.csv(cons, file.path(outdir, "condition_specific_isoforms.csv"), row.names = FALSE)
cat("Top switching genes written:", nrow(cons), "rows\n")
