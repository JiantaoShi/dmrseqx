# dmrseqx() against dmrseq() on small data, with one thread each. Covers what the ESCC comparison
# doesn't: missing coverage, adjustment and matching covariates, unbalanced groups, small groups.
# Usage: Rscript check_dmrseq.R [lib]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({
  library(dmrseqx)
  library(dmrseq); library(bsseq); library(GenomicRanges); library(BiocParallel); library(data.table)
})
options(width = 170, digits = 4)
data("BS.chr21", package = "dmrseq")
keep_groups <- function(bs, group) {
  cv <- as.matrix(getCoverage(bs, type = "Cov"))
  ok <- Reduce(`&`, lapply(unique(group), function(l) rowSums(cv[, group == l, drop = FALSE]) > 0))
  ok <- ok & rowSums(cv == 0) < ncol(bs) - 1
  bs[ok]
}
compare <- function(name, bs, ...) {
  set.seed(11)
  t0 <- Sys.time()
  a <- suppressWarnings(suppressMessages(dmrseq(bs, BPPARAM = SerialParam(), verbose = FALSE, ...)))
  ta <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  set.seed(11)
  t0 <- Sys.time()
  b <- suppressMessages(dmrseqx(bs, verbose = FALSE, threads = 1L, ...))
  tb <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(a) || is.null(b)) { cat(name, ": a", length(a), "b", length(b), "\n"); return(invisible(NULL)) }
  A <- as.data.table(as.data.frame(a)); B <- as.data.table(as.data.frame(b))
  m <- merge(A[, .(seqnames, start, end, L, area, beta, stat, pval, qval, index.start)], B[, .(seqnames, start, end, area_x = area, beta_x = beta, stat_x = stat, pval_x = pval, qval_x = qval, index_x = index.start)], by = c("seqnames", "start", "end"))
  m[, dt := abs(stat_x - stat)]
  called <- function(q) !is.na(q) & q <= 0.05
  cat(sprintf("%-22s regions %5d %5d same %5d | area %.1e index %d | NA %d %d | long max dt %.2e | short within 1e-3 %.4f max dt %.2f | |dp| max %.2e | calls %d %d both %d | time %.1f s %.2f s\n",
              name, nrow(A), nrow(B), nrow(m), max(abs(m$area - m$area_x), na.rm = TRUE), sum(m$index.start != m$index_x), sum(is.na(A$stat)), sum(is.na(B$stat)),
              max(c(0, m[L >= 40]$dt), na.rm = TRUE), mean(m[L < 40]$dt < 1e-3, na.rm = TRUE), max(c(0, m[L < 40]$dt), na.rm = TRUE),
              max(abs(m$pval - m$pval_x), na.rm = TRUE), sum(called(A$qval)), sum(called(B$qval)), sum(called(m$qval) & called(m$qval_x)), ta, tb))
  invisible(m)
}

# A, B: real data of two cell types with two replicates each, parts of chr21
slice <- BS.chr21[240001:250000, ]
ct <- pData(slice)$CellType
cv <- as.matrix(getCoverage(slice, type = "Cov"))
compare("BS.chr21 complete", slice[rowSums(cv == 0) == 0], testCovariate = "CellType", cutoff = 0.05)
m_b <- compare("BS.chr21 missing", keep_groups(slice, ct), testCovariate = "CellType", cutoff = 0.05)

# Simulated counts on the positions of the slice
simulate <- function(n1, n2, missing = 0, batch_effect = 0, seed = 3) {
  set.seed(seed)
  pos <- start(slice)[seq_len(6000)]; L <- length(pos); N <- n1 + n2
  base <- pmin(pmax(as.numeric(stats::filter(rbinom(L, 1, 0.7), rep(1 / 25, 25), sides = 2)), 0.05), 0.95)
  base[is.na(base)] <- 0.7
  dmr <- rep(0, L)
  for (s in sample(L - 200, 20)) dmr[s:(s + sample(8:150, 1))] <- sample(c(-0.3, 0.3), 1)
  group <- rep(c("A", "B"), c(n1, n2))
  batch <- rep(c("x", "y"), length.out = N)
  M <- Cov <- matrix(0L, L, N)
  for (j in seq_len(N)) {
    cvj <- rpois(L, 15) + 1L
    if (missing > 0) cvj[runif(L) < missing] <- 0L
    p <- pmin(pmax(base + (group[j] == "B") * dmr + (batch[j] == "y") * batch_effect, 0.01), 0.99)
    Cov[, j] <- cvj
    M[, j] <- rbinom(L, cvj, p)
  }
  ids <- paste0("S", seq_len(N))
  bs <- BSseq(M = M, Cov = Cov, pos = pos, chr = rep("chr21", L), sampleNames = ids,
              pData = data.frame(group = group, batch = batch, id = paste0("s", c(seq_len(n1), seq_len(n2))), row.names = ids))
  bs
}
compare("sim 3v3", simulate(3, 3), testCovariate = "group")
compare("sim 3v4 unbalanced", simulate(3, 4), testCovariate = "group")
s <- simulate(3, 3, missing = 0.1); compare("sim 3v3 missing 10%", keep_groups(s, pData(s)$group), testCovariate = "group")
s <- simulate(4, 3, missing = 0.3, seed = 4); compare("sim 4v3 missing 30%", keep_groups(s, pData(s)$group), testCovariate = "group")
compare("sim 3v3 adjust", simulate(3, 3, batch_effect = 0.1), testCovariate = "group", adjustCovariate = "batch")
s <- simulate(3, 3, missing = 0.1, batch_effect = 0.1, seed = 5); compare("sim adjust missing", keep_groups(s, pData(s)$group), testCovariate = "group", adjustCovariate = "batch")
compare("sim 3v3 match", simulate(3, 3), testCovariate = "group", matchCovariate = "id")
compare("sim minNumRegion 3", simulate(3, 3), testCovariate = "group", minNumRegion = 3, cutoff = 0.2)
compare("sim no smoothing", simulate(3, 3), testCovariate = "group", smooth = FALSE, cutoff = 0.2)
cat("CHECK_DMRSEQ_DONE\n")
