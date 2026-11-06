# dmrseqx() on the whole genome against the checkpoints of a run of the R implementation, with the
# labels of that run.
# Usage: Rscript check_genome.R <bsseq.rds> <checkpoint dir> <threads> <out dir> [lib] [fit]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 5) .libPaths(c(args[5], .libPaths()))
fit <- if (length(args) >= 6) args[6] else "dmrseq"
suppressPackageStartupMessages({ library(dmrseqx); library(bsseq); library(data.table); library(GenomicRanges) })
BS <- args[1]; CK <- args[2]; threads <- as.integer(args[3]); OUT <- args[4]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
options(width = 170, digits = 4, dmrseqx.keepNull = TRUE)

t0 <- Sys.time()
bs <- readRDS(BS)
cov <- as.matrix(getCoverage(bs, type = "Cov"))
bs <- bs[rowSums(cov == 0) == 0]
rm(cov); invisible(gc())
log_msg("loaded ", nrow(bs), " CpGs, ", ncol(bs), " samples, ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

fs <- list.files(CK, pattern = "\\.rds$")
lab <- sub("_chr.*", "", fs)
tab <- table(lab)
obs <- paste(as.integer(pData(bs)$group == "case"), collapse = "")
done <- setdiff(names(tab)[tab == length(unique(seqnames(bs)))], obs)
perms <- sapply(done, function(s) as.integer(strsplit(s, "")[[1]]))
log_msg("permutations finished in the reference: ", length(done))

set.seed(1)
t0 <- Sys.time()
res <- dmrseqx(bs, testCovariate = "group", cutoff = 0.1, threads = threads, perms = perms, verbose = TRUE, fit = fit)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
log_msg("dmrseqx: ", round(el, 1), " s with ", threads, " thread(s), ", length(res), " regions, null ", metadata(res)$dmrseqx$nullStatistics)
A <- as.data.table(as.data.frame(res))
null <- as.data.table(metadata(res)$null)
fwrite(A[, .(chr = seqnames, start, end, L, area, beta, stat, pval, qval)], file.path(OUT, "regions.tsv.gz"), sep = "\t")
fwrite(null[, .(pass, chr, start, end, L, area, beta, stat)], file.path(OUT, "null.tsv.gz"), sep = "\t")

labels <- metadata(res)$dmrseqx$labels
rows <- list()
for (k in seq_along(labels)) {
  ref <- rbindlist(lapply(fs[lab == labels[k]], function(f) as.data.table(readRDS(file.path(CK, f)))))
  me <- if (k == 1) A[, .(chr = as.character(seqnames), start, end, L, area, beta, stat)] else null[pass == k - 1L, .(chr, start, end, L, area, beta, stat)]
  m <- merge(ref[, .(chr, start, end, L, area, beta, stat)], me[, .(chr, start, end, area_x = area, beta_x = beta, stat_x = stat)], by = c("chr", "start", "end"))
  m[, dt := abs(stat_x - stat)]
  rows[[k]] <- data.table(pass = k - 1L, ref = nrow(ref), mine = nrow(me), same = nrow(m), na_ref = sum(is.na(ref$stat)), na_mine = sum(is.na(me$stat)),
                          long_n = sum(m$L >= 40), long_max = max(c(0, m[L >= 40]$dt), na.rm = TRUE), short_n = sum(m$L < 40),
                          short_1e3 = mean(m[L < 40]$dt < 1e-3, na.rm = TRUE), short_q99 = quantile(m[L < 40]$dt, 0.99, na.rm = TRUE), short_max = max(c(0, m[L < 40]$dt), na.rm = TRUE))
  if (k == 1) {
    r0 <- ref
    m0 <- m
  }
}
print(rbindlist(rows))

# calls with the null of the reference for the reference, and our own null for ours
refNull <- rbindlist(lapply(fs[lab %in% done], function(f) as.data.table(readRDS(file.path(CK, f)))[, .(stat)]))
ns <- sort(abs(refNull$stat[!is.na(refNull$stat)]))
p <- (1 + length(ns) - findInterval(abs(r0$stat), ns, left.open = TRUE)) / (1 + length(ns))
r0[, q := p.adjust(p, "BH")]
cmp <- merge(r0[, .(chr, start, end, L, stat, q)], A[, .(chr = as.character(seqnames), start, end, stat_x = stat, q_x = qval)], by = c("chr", "start", "end"))
a <- cmp$q <= 0.05; b <- cmp$q_x <= 0.05
cat(sprintf("calls q <= 0.05: reference %d, dmrseqx %d, both %d, only reference %d, only dmrseqx %d, of them L < 40: %d, %d\n", sum(a), sum(b), sum(a & b), sum(a & !b), sum(b & !a), sum(a & !b & cmp$L < 40), sum(b & !a & cmp$L < 40)))
cat(sprintf("share of calls that differ: %.3f%%\n", 100 * sum(a != b) / sum(a | b)))
log_msg("GENOME_DONE")
