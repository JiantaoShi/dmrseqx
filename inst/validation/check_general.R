# Several groups and a continuous covariate, pass by pass against dmrseq:::bumphunt with the same
# permutations, one thread each. Usage: Rscript check_general.R sim [lib] [regex of the cases]
# or Rscript check_general.R chr21 [lib] [bsseq.rds]
args <- commandArgs(trailingOnly = TRUE)
what <- args[1]
if (length(args) >= 2) .libPaths(c(args[2], .libPaths()))
suppressPackageStartupMessages({
  library(dmrseqx); library(dmrseq); library(bsseq); library(GenomicRanges); library(BiocParallel); library(data.table)
})
options(width = 170, digits = 4, dmrseqx.keepNull = TRUE)
BiocParallel::register(SerialParam())

# dmrseq stops when a permutation leaves a group without coverage at a CpG: the design of the samples
# with coverage isn't of full rank there and backsolve() fails. As a reference for these passes, this
# gives NA at such CpGs and leaves everything else to dmrseq, one call for each pattern of missing
# samples, since a single pattern on a chromosome stops dmrseq too. dmrseqx gives no effect either.
getEstimate0 <- dmrseq:::getEstimate
singularCpGs <- 0L
getEstimateNA <- function(mat, design, coeff, coeff.adj) {
  m <- as.matrix(mat)
  na <- is.na(m)
  pattern <- drop(na %*% 2^(seq_len(ncol(m)) - 1))
  diff <- sd <- rep(NA_real_, nrow(m))
  for (p in unique(pattern)) {
    rows <- which(pattern == p)
    if (qr(design[!na[rows[1], ], , drop = FALSE])$rank < ncol(design)) {
      singularCpGs <<- singularCpGs + length(rows)
      next
    }
    out <- getEstimate0(m[rows, , drop = FALSE], design, coeff, coeff.adj)
    diff[rows] <- out$meth.diff
    sd[rows] <- out$sd.meth.diff
  }
  list(meth.diff = diff, sd.meth.diff = sd, df.residual = ncol(m) - ncol(design))
}

compare <- function(name, bs, test, adjust = NULL, nPerm = 2, cutoff = 0.1, seed = 1) {
  pD <- as.data.frame(pData(bs))
  testIdx <- match(test, colnames(pD))
  adjIdx <- if (is.null(adjust)) NULL else match(adjust, colnames(pD))
  N <- ncol(bs)
  testCov <- pD[[testIdx]]
  if (is.factor(testCov)) testCov <- droplevels(testCov)
  k <- length(unique(testCov))
  continuous <- is.numeric(testCov) && k > 2
  coeff <- if (continuous) 2L else seq(2L, 2L + k - 2L)
  if (is.null(adjust)) {
    design <- model.matrix(~testCov)
    coeff.adj <- NULL
  } else {
    mmdat <- cbind(data.frame(testCov = testCov), pD[, adjIdx, drop = FALSE])
    design <- model.matrix(as.formula(paste0("~", paste0(colnames(mmdat), collapse = " + "))), data = mmdat)
    coeff.adj <- (max(coeff) + 1):ncol(design)
  }
  colnames(design)[coeff] <- test
  set.seed(seed)
  perms <- replicate(nPerm, sample(N))
  t0 <- Sys.time()
  x <- suppressMessages(dmrseqx(bs, testCovariate = test, adjustCovariate = adjust, cutoff = cutoff, perms = perms, threads = 1L, verbose = FALSE))
  tx <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  A <- as.data.table(as.data.frame(x))
  A <- A[, !c("seqnames", "width", "strand", "pval", "qval", "index.start", "index.end", "index.width"), with = FALSE]
  A[, pass := 0L]
  null <- as.data.table(metadata(x)$null)
  all <- rbind(A, null[, colnames(A), with = FALSE])
  tr <- 0
  for (j in 0:nPerm) {
    designr <- design
    bsr <- bs
    if (j > 0) {
      designr[, coeff] <- designr[perms[, j], coeff]
      pData(bsr) <- pData(bs)[perms[, j], , drop = FALSE]
    }
    reference <- function() suppressWarnings(suppressMessages(dmrseq:::bumphunt(bs = bsr, design = designr, coeff = coeff, coeff.adj = coeff.adj,
      minInSpan = 30, minNumRegion = 5, cutoff = cutoff, maxGap = 1000, maxGapSmooth = 2500, smooth = TRUE, bpSpan = 1000,
      verbose = FALSE, parallel = FALSE, block = FALSE, blockSize = 5000, chrsPerChunk = 1, fact = !continuous, adjustCovariate = adjIdx)))
    t0 <- Sys.time()
    ref <- tryCatch(reference(), error = identity)
    note <- ""
    if (inherits(ref, "error")) {
      message(name, " pass ", j, ": dmrseq stops: ", conditionMessage(ref))
      singularCpGs <<- 0L
      assignInNamespace("getEstimate", getEstimateNA, ns = "dmrseq")
      ref <- reference()
      assignInNamespace("getEstimate", getEstimate0, ns = "dmrseq")
      note <- sprintf(" | dmrseq stops, reference has NA at %d CpGs", singularCpGs)
    }
    tr <- tr + as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ref <- as.data.table(ref)
    me <- all[pass == j]
    if (nrow(ref) == 0) { cat(sprintf("%-31s pass %d: no regions in dmrseq, %d in dmrseqx%s\n", name, j, nrow(me), note)); next }
    m <- merge(ref, me, by = c("start", "end"), suffixes = c("", "_x"))
    m[, dt := abs(stat_x - stat)]
    bcols <- grep("^beta", colnames(ref), value = TRUE)
    db <- max(sapply(bcols, function(b) max(abs(m[L >= 40][[b]] - m[L >= 40][[paste0(b, "_x")]]), 0, na.rm = TRUE)))
    cat(sprintf("%-31s pass %d: regions %5d %5d same %5d | area %.1e | NA %d %d | long n %4d max dt %.2e max dbeta %.2e | short n %4d within 1e-3 %.3f max dt %.2f%s\n",
                name, j, nrow(ref), nrow(me), nrow(m), max(abs(m$area - m$area_x), na.rm = TRUE), sum(is.na(ref$stat)), sum(is.na(me$stat)),
                sum(m$L >= 40), max(c(0, m[L >= 40]$dt), na.rm = TRUE), db, sum(m$L < 40), mean(m[L < 40]$dt < 1e-3, na.rm = TRUE), max(c(0, m[L < 40]$dt), na.rm = TRUE), note))
  }
  cat(sprintf("%-31s time dmrseq %.1f s, dmrseqx %.2f s, %d passes; columns: %s\n", name, tr, tx, nPerm + 1, paste(colnames(A), collapse = " ")))
}

if (what == "sim") {
  data("BS.chr21", package = "dmrseq")
  positions <- start(BS.chr21)[240001:246000]
  simulate <- function(cov_of_sample, missing = 0, seed = 3, batch = NULL) {
    # cov_of_sample: a list with the effect multiplier of each sample
    set.seed(seed)
    L <- length(positions); N <- length(cov_of_sample$effect)
    base <- pmin(pmax(as.numeric(stats::filter(rbinom(L, 1, 0.7), rep(1 / 25, 25), sides = 2)), 0.05), 0.95)
    base[is.na(base)] <- 0.7
    dmr <- rep(0, L)
    for (s in sample(L - 200, 20)) dmr[s:(s + sample(8:150, 1))] <- sample(c(-0.3, 0.3), 1)
    M <- Cov <- matrix(0L, L, N)
    for (j in seq_len(N)) {
      cvj <- rpois(L, 15) + 1L
      if (missing > 0) cvj[runif(L) < missing] <- 0L
      p <- pmin(pmax(base + cov_of_sample$effect[j] * dmr, 0.01), 0.99)
      Cov[, j] <- cvj
      M[, j] <- rbinom(L, cvj, p)
    }
    ids <- paste0("S", seq_len(N))
    pd <- data.frame(x = cov_of_sample$x, row.names = ids)
    if (!is.null(batch)) pd$batch <- batch
    BSseq(M = M, Cov = Cov, pos = positions, chr = rep("chr21", L), sampleNames = ids, pData = pd)
  }
  # no_empty_group: also without CpGs where a permutation could leave a group without coverage
  keep_rows <- function(bs, no_empty_group = FALSE) {
    cv <- as.matrix(getCoverage(bs, type = "Cov"))
    x <- pData(bs)$x
    ok <- rowSums(cv == 0) < ncol(bs) - 1
    if (!is.numeric(x)) {
      ok <- ok & Reduce(`&`, lapply(unique(x), function(l) rowSums(cv[, x == l, drop = FALSE]) > 0))
      if (no_empty_group) ok <- ok & rowSums(cv == 0) < min(table(x))
    }
    bs[ok]
  }
  three <- list(x = rep(c("A", "B", "C"), each = 3), effect = rep(c(0, 1, -1), each = 3))
  four <- list(x = factor(rep(c("d", "a", "c", "b"), each = 2), levels = c("d", "a", "c", "b")), effect = rep(c(0, 1, -0.5, 0.5), each = 2))
  set.seed(9)
  cont <- list(x = c(-1, -0.6, -0.3, 0, 0.2, 0.5, 0.8, 1.1), effect = c(-1, -0.6, -0.3, 0, 0.2, 0.5, 0.8, 1.1))
  scaled <- list(x = 10 * c(1, 1.4, 1.7, 2, 2.2, 2.5, 2.8, 3.1), effect = c(-1, -0.6, -0.3, 0, 0.2, 0.5, 0.8, 1.1))
  only <- if (length(args) >= 3) args[3] else "."
  run <- function(name, ...) if (grepl(only, name)) compare(name, ...)
  run("3 groups", simulate(three), "x")
  # about 5 h of dmrseq: small groups, so that noise reaches the cutoff and regions have 2,000 CpGs
  run("4 groups, factor levels", simulate(four), "x")
  run("3 groups, missing 10%", keep_rows(simulate(three, missing = 0.1, seed = 4)), "x")
  run("3 groups, missing 30%", keep_rows(simulate(three, missing = 0.3, seed = 5)), "x")
  run("3 groups, missing 10% filtered", keep_rows(simulate(three, missing = 0.1, seed = 4), TRUE), "x")
  run("3 groups, missing 30% filtered", keep_rows(simulate(three, missing = 0.3, seed = 5), TRUE), "x")
  run("3 groups, adjust", simulate(three, batch = rep(c("u", "v", "w"), 3)), "x", adjust = "batch")
  run("continuous", simulate(cont), "x")
  # the effect between the values 0 and 1 is an extrapolation here, and small
  run("continuous, large values", simulate(scaled), "x", cutoff = 0.005)
  run("continuous, missing 10%", keep_rows(simulate(cont, missing = 0.1, seed = 6)), "x")
  run("continuous, missing 30%", keep_rows(simulate(cont, missing = 0.3, seed = 7)), "x")
  run("continuous, adjust", simulate(cont, batch = rep(c("u", "v"), 4)), "x", adjust = "batch")
} else {
  bs <- readRDS(args[3])
  bs <- bs[as.character(seqnames(bs)) == "chr21"]
  samples <- c("ctrl_N1", "ctrl_N2", "ctrl_N6", "case_T1", "case_T2", "case_T6", "case_T7", "case_T8", "case_T9")
  bs <- bs[, samples]
  cv <- as.matrix(getCoverage(bs, type = "Cov"))
  bs <- bs[rowSums(cv == 0) == 0]
  cat("ESCC chr21:", nrow(bs), "CpGs,", ncol(bs), "samples\n")
  # cutoff 0.2: at 0.1 the longest regions have more than 2,000 CpGs, many hours each for the dense fit of dmrseq
  b3 <- bs
  pData(b3)$x <- rep(c("normal", "tumorA", "tumorB"), each = 3)
  compare("ESCC chr21, 3 groups", b3, "x", nPerm = 1, cutoff = 0.2)
  b4 <- bs
  # a continuous covariate along the tumor status
  pData(b4)$x <- c(0, 0.1, 0.2, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
  compare("ESCC chr21, continuous", b4, "x", nPerm = 1, cutoff = 0.2)
}
cat("GENERAL_DONE\n")
