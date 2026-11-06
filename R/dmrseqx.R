

dmrseqx <- function(bs, testCovariate, adjustCovariate = NULL, cutoff = 0.1,
                    minNumRegion = 5, smooth = TRUE, bpSpan = 1000,
                    minInSpan = 30, maxGapSmooth = 2500, maxGap = 1000,
                    verbose = TRUE, maxPerms = 10, matchCovariate = NULL,
                    stat = "stat", threads = 1L, fit = c("dmrseq", "global"),
                    perms = NULL) {
  stopifnot(methods::is(bs, "BSseq"))
  fit <- match.arg(fit)
  if (!is.numeric(cutoff) || length(cutoff) != 1L || abs(cutoff) > 1 || cutoff <= 0) {
    stop("'cutoff' must be one number above 0 and at most 1.")
  }
  if (minNumRegion < 3) stop("minNumRegion must be at least 3")
  if (!(stat %in% c("L", "area", "beta", "stat", "avg"))) {
    stop("Specified '", stat, "' as the test statistic which is not in the results. ",
         "Please specify a valid name from one of L, area, beta, stat, or avg")
  }
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("'threads' must be one positive integer.")
  }

  pD <- as.data.frame(bsseq::pData(bs))
  cols <- .resolve_covariates(pD, testCovariate, adjustCovariate, matchCovariate)
  testCovariate <- cols$test
  adjustCovariate <- cols$adjust
  mC <- cols$match

  if (is.unsorted(bs)) stop("'bs' must be sorted. Use 'sort(bs)'.")

  setup <- .make_design(pD, testCovariate, adjustCovariate, verbose)
  design <- setup$design
  coeff <- setup$coeff
  sampleSize <- setup$sampleSize
  twoGroup <- setup$kind == "two groups"
  if (!is.null(matchCovariate) && !twoGroup) {
    stop("Covariate matching can only be carried out for 2-group comparisons")
  }
  if (!is.null(matchCovariate) && choose(nrow(design), min(sampleSize)) >= 5e5) {
    stop("matchCovariate can't be used when the sample size is large enough to yield more ",
         "than 500000 possible permutations. Perhaps you'd like to add an adjustCovariate instead?")
  }
  if (stat == "beta" && length(coeff) > 1) {
    stop("Specified 'beta' as the test statistic, but several groups have several coefficients.")
  }
  if (verbose) {
    if (twoGroup) {
      message("Condition: ", unique(pD[, testCovariate][which(design[, coeff] == 1)]), " vs ",
              unique(pD[, testCovariate][which(design[, coeff] == 0)]))
    }
    if (!is.null(adjustCovariate)) {
      message("Adjusting for covariate (s): ", paste(colnames(pD)[adjustCovariate], collapse = ", "))
    }
    if (!is.null(matchCovariate)) {
      message("Matching permutations on covariate: ", colnames(pD)[mC])
    }
  }

  # counts of all samples in memory, as integers
  M <- as.matrix(bsseq::getCoverage(bs, type = "M"))
  Cov <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
  if (any(M > .Machine$integer.max) || any(Cov > .Machine$integer.max)) {
    stop("Counts above the largest integer are not supported.")
  }
  storage.mode(M) <- "integer"
  storage.mode(Cov) <- "integer"
  dimnames(M) <- dimnames(Cov) <- NULL

  if (setup$kind == "continuous") {
    uncovered <- matrixStats::rowSums2(Cov == 0L) >= ncol(Cov) - 1
    if (any(uncovered)) {
      stop(sum(uncovered), " loci have zero coverage in at least ", ncol(Cov) - 1,
           " samples. Please remove these loci before running dmrseqx")
    }
  } else {
    zeroInGroup <- Reduce(`|`, lapply(unique(pD[[testCovariate]]), function(level) {
      matrixStats::rowSums2(Cov, cols = which(pD[[testCovariate]] == level)) == 0
    }))
    if (any(zeroInGroup)) {
      stop(sum(zeroInGroup), " loci have zero coverage in all samples of at least one condition. ",
           "Please remove these loci before running dmrseqx")
    }
  }

  rl <- GenomicRanges::seqnames(bs)
  chrs <- as.character(S4Vectors::runValue(rl))
  if (anyDuplicated(chrs)) stop("'bs' must be sorted. Use 'sort(bs)'.")
  chrStart <- c(0, cumsum(as.numeric(S4Vectors::runLength(rl))))
  pos <- as.integer(GenomicRanges::start(bs))

  bpSpan2 <- mean(vapply(seq_along(chrs), function(i) {
    first <- chrStart[i] + 1
    last <- chrStart[i + 1]
    minInSpan * (pos[last] - pos[first] + 1) / (last - first + 1)
  }, numeric(1)), na.rm = TRUE)
  covQ75 <- stats::quantile(matrixStats::rowMeans2(Cov), 0.75, names = FALSE)

  passes <- .labellings(pD, design, coeff, testCovariate, adjustCovariate, mC, sampleSize,
                        maxPerms, perms, verbose, twoGroup)
  nPass <- length(passes$design)
  designs <- do.call(cbind, passes$design)
  X <- do.call(cbind, passes$X)
  storage.mode(designs) <- "double"
  storage.mode(X) <- "double"

  params <- c(cutoff, bpSpan, minInSpan, bpSpan2, maxGapSmooth, maxGap, covQ75,
              minNumRegion, as.numeric(isTRUE(smooth)), as.numeric(fit == "global"),
              as.numeric(twoGroup), as.numeric(setup$kind != "continuous"), length(coeff), ncol(design))
  seed <- sample.int(.Machine$integer.max, 1L)
  if (verbose) {
    message("Detecting candidate regions with coefficient larger than ", cutoff,
            " in magnitude, using ", threads, " thread(s).")
  }
  res <- .Call(C_dmrseqx_run, M, Cov, pos, chrStart, designs, X, as.double(params),
               .t_critical(), .grubbs_critical(), threads, as.double(seed), isTRUE(verbose))
  betas <- matrix(res$beta, ncol = length(coeff))
  colnames(betas) <- passes$betaNames
  res$beta <- NULL
  res <- cbind(as.data.frame(res), as.data.frame(betas))
  res$chr <- chrs[res$chrom]

  OBS <- res[res$pass == 0L, , drop = FALSE]
  FLIP <- res[res$pass > 0L, , drop = FALSE]
  if (nrow(OBS) == 0L) {
    message("No candidate regions pass the cutoff of ", cutoff)
    return(NULL)
  }
  .finish(OBS, FLIP, stat, passes$betaNames, fit, seed, passes$labels)
}

.resolve_covariates <- function(pD, testCovariate, adjustCovariate, matchCovariate) {
  if (is.character(testCovariate)) {
    if (length(testCovariate) > 1) stop("Only one testCovariate can be specified")
    if (is.character(adjustCovariate) && sum(testCovariate %in% adjustCovariate) > 0) {
      stop("adjustCovariate can't contain testCovariate")
    }
    if (is.character(matchCovariate) && sum(testCovariate %in% matchCovariate) > 0) {
      stop("matchCovariate can't contain testCovariate")
    }
    testCovariate <- which(colnames(pD) == testCovariate)
    if (length(testCovariate) == 0) {
      stop("testCovariate not found in pData(). Please specify a valid testCovariate")
    }
  }
  if (is.character(adjustCovariate)) {
    if (is.character(matchCovariate) && any(matchCovariate == adjustCovariate)) {
      stop("matchCovariate can't be identical to adjustCovariate")
    }
    adjustCovariate <- which(colnames(pD) %in% adjustCovariate)
    if (length(adjustCovariate) == 0) {
      stop("adjustCovariate not found in pData(). Please specify a valid adjustCovariate")
    }
  }
  if (ncol(pD) < max(testCovariate, adjustCovariate)) {
    stop("pData(bs) has too few columns. Please specify valid covariates to use in the analysis")
  }
  mC <- NULL
  if (!is.null(matchCovariate)) {
    if (length(matchCovariate) > 1) stop("Covariate matching can only be carried out for one covariate")
    if (is.character(matchCovariate)) {
      mC <- grep(matchCovariate, colnames(pD))
      if (length(mC) == 0) stop("no column in pData() found that matches the matchCovariate")
      if (length(mC) > 1) stop("matchCovariate matches more than one column in pData()")
    } else {
      stopifnot(matchCovariate <= ncol(pD))
      mC <- matchCovariate
    }
  }
  list(test = testCovariate, adjust = adjustCovariate, match = mC)
}

.make_design <- function(pD, testCovariate, adjustCovariate, verbose) {
  coeff <- 2L
  testCov <- pD[, testCovariate]
  if (is.factor(testCov)) testCov <- droplevels(testCov)
  if (any(is.na(testCov))) {
    stop("Missing values found in the testCovariate. Please remove these samples from the analysis.")
  }
  nValues <- length(unique(testCov))
  sampleSize <- table(testCov)
  if (nValues == 1L) {
    stop("testCovariate has only one value.")
  } else if (nValues > 2L && !is.numeric(testCov)) {
    if (nValues > 60L) stop("testCovariate has too many groups.")
    if (verbose) {
      message("Performing a global test of H0: no difference among ", nValues, " groups (assuming the test ",
              "covariate ", colnames(pD)[testCovariate], " is a factor).")
    }
    coeff <- seq(coeff, coeff + nValues - 2L)
    kind <- "several groups"
  } else if (nValues > 2L) {
    if (verbose) message("Assuming the test covariate ", colnames(pD)[testCovariate], " is continuous.")
    kind <- "continuous"
  } else {
    if (verbose) message("Assuming the test covariate ", colnames(pD)[testCovariate], " is a factor.")
    if (min(sampleSize) < 2) {
      stop("At least one group has only one sample! Replicates are required to run dmrseqx.")
    }
    testCov <- as.factor(testCov)
    kind <- "two groups"
  }
  if (!is.null(adjustCovariate)) {
    mmdat <- data.frame(testCov = testCov)
    adjustCov <- pD[, adjustCovariate, drop = FALSE]
    if (any(is.na(adjustCov))) {
      stop("Missing values found in the adjustCovariate. Please remove these samples from the analysis.")
    }
    if (any(vapply(adjustCov, function(x) length(unique(x)), 1L) < 2)) {
      stop("At least one adjust covariate is constant across samples. ",
           "Please remove this covariate from the model and try again.")
    }
    for (f in seq_len(ncol(adjustCov))) {
      if (is.factor(adjustCov[, f])) adjustCov[, f] <- droplevels(adjustCov[, f])
    }
    mmdat <- cbind(mmdat, adjustCov)
    design <- stats::model.matrix(stats::as.formula(paste0("~", paste0(colnames(mmdat), collapse = " + "))),
                                  data = mmdat)
  } else {
    design <- stats::model.matrix(~testCov)
  }
  e <- eigen(crossprod(design), symmetric = TRUE, only.values = TRUE)$values
  if (!(e[1] > 0 && abs(e[length(e)] / e[1]) > 1e-13)) stop("Design matrix is not full rank")
  if (nrow(pD) < ncol(design) + 1) {
    stop("Not enough degrees of freedom to estimate ", ncol(design) - 1, " covariates using ",
         nrow(pD), " samples. Please use a larger number of samples, or specify fewer adjust covariates.")
  }
  list(design = design, coeff = coeff, sampleSize = sampleSize, kind = kind)
}

# The permutation scheme of dmrseq. For two groups and not too many samples the columns of the
# result hold the samples of one group, otherwise they are permutations of the samples.
.draw_permutations <- function(design, coeff, sampleSize, maxPerms, twoGroup) {
  n <- nrow(design)
  if (twoGroup && choose(n, min(sampleSize)) < 5e5) {
    perms <- utils::combn(seq_len(n), min(sampleSize))
    if (length(unique(table(design[, coeff]))) == 1) {
      perms <- perms[, seq_len(ncol(perms) / 2), drop = FALSE]
    }
    rmv <- which(apply(perms, 2, function(p) length(unique(design[p, coeff])) == 1))
    if (length(rmv) > 0) perms <- perms[, -rmv, drop = FALSE]
    if (maxPerms < ncol(perms)) {
      # preference to those with the least similarity to the original partition
      similarity <- apply(perms, 2, function(x) max(table(design[x, coeff])))
      all_perms <- perms
      perms <- NULL
      levs <- sort(unique(similarity))
      l <- 1
      num <- 0
      while (!(num == maxPerms) && l <= length(levs)) {
        candidates <- which(similarity == levs[l])
        keep <- candidates[sample.int(length(candidates), min(maxPerms - num, length(candidates)))]
        perms <- cbind(perms, all_perms[, keep, drop = FALSE])
        l <- l + 1
        num <- ncol(perms)
      }
    }
    return(perms)
  }
  # unrestricted permutations of the samples
  perms <- as.matrix(seq_len(n))
  for (p in seq_len(maxPerms)) {
    tries <- 0
    candidate <- sample(seq_len(n), n)
    duplicate <- function(x) {
      sum(apply(perms, 2, function(z) isTRUE(all.equal(z, x)))) > 0 ||
        sum(apply(perms, 2, function(z) isTRUE(all.equal(z, rev(x))))) > 0
    }
    while (duplicate(candidate) && tries <= 20) {
      candidate <- sample(seq_len(n), n)
      tries <- tries + 1
    }
    if (tries <= 20) perms <- cbind(perms, candidate)
  }
  perms[, -1, drop = FALSE]
}

# The covariates of one labelling without intercept, the covariate of interest first
.fit_covariates <- function(pDat, testCovariate, adjustCovariate) {
  g.fac <- pDat[, testCovariate]
  if (is.factor(g.fac)) g.fac <- droplevels(g.fac)
  if (length(unique(g.fac)) == 2L) g.fac <- as.factor(g.fac)
  mm <- data.frame(g.fac = g.fac)
  if (length(adjustCovariate) > 0) {
    adjust <- pDat[, adjustCovariate, drop = FALSE]
    for (f in seq_len(ncol(adjust))) {
      if (is.factor(adjust[, f])) adjust[, f] <- droplevels(adjust[, f])
    }
    mm <- cbind(mm, adjust)
  }
  stats::model.matrix(stats::as.formula(paste0("~", paste0(colnames(mm), collapse = " + "))), data = mm)[, -1L, drop = FALSE]
}

# One entry for each labelling: the observed one, then the permutations
.labellings <- function(pD, design, coeff, testCovariate, adjustCovariate, mC, sampleSize,
                        maxPerms, perms, verbose, twoGroup) {
  firstX <- .fit_covariates(pD, testCovariate, adjustCovariate)
  betaNames <- if (length(coeff) > 1) paste0("beta_", sub("^g\\.fac", "", colnames(firstX)[seq_along(coeff)])) else "beta"
  out_design <- list(unname(design))
  out_X <- list(unname(firstX))
  given <- !is.null(perms)
  if (given) {
    perms <- as.matrix(perms)
    if (twoGroup) {
      # one column per permutation, with the 0 or 1 of each sample
      if (nrow(perms) != nrow(design) || !all(perms %in% c(0, 1))) {
        stop("'perms' must have one row per sample holding 0 or 1.")
      }
    } else if (nrow(perms) != nrow(design) || !all(apply(perms, 2, function(p) identical(sort(as.integer(p)), seq_len(nrow(design)))))) {
      # one column per permutation of the samples
      stop("'perms' must have one row per sample, and each column must be a permutation of the samples.")
    }
  } else {
    if (verbose) {
      message("Performing permutations of the condition across samples ",
              "to generate a null distribution of region test statistics")
    }
    perms <- .draw_permutations(design, coeff, sampleSize, maxPerms, twoGroup)
  }
  levs <- unique(pD[[testCovariate]])
  groupOf <- as.factor(if (is.factor(pD[[testCovariate]])) droplevels(pD[[testCovariate]]) else pD[[testCovariate]])
  adjustFollows <- FALSE
  for (j in seq_len(ncol(perms))) {
    designr <- design
    pDr <- pD
    if (given && twoGroup) {
      designr[, coeff] <- perms[, j]
      pDr[[testCovariate]][perms[, j] == 1] <- levels(groupOf)[2]
      pDr[[testCovariate]][perms[, j] == 0] <- levels(groupOf)[1]
      adjustFollows <- TRUE
    } else if (nrow(perms) != nrow(designr)) {
      reorder <- perms[, j]
      designr[, coeff] <- 0
      designr[reorder, coeff] <- 1
      pDr[[testCovariate]] <- rep(levs[1], nrow(pDr))
      pDr[[testCovariate]][reorder] <- levs[2]
      if (!all(sort(pD[[testCovariate]]) == sort(pDr[[testCovariate]]))) {
        designr[, coeff] <- 1
        designr[reorder, coeff] <- 0
        pDr[[testCovariate]] <- rep(levs[2], nrow(pDr))
        pDr[[testCovariate]][reorder] <- levs[1]
      }
      adjustFollows <- TRUE
    } else {
      reorder <- perms[, j]
      designr[, coeff] <- designr[reorder, coeff]
      pDr <- pD[reorder, , drop = FALSE]
      adjustFollows <- FALSE
    }
    if (adjustFollows && length(adjustCovariate) > 0) {
      # the adjustment covariates follow the samples of the same group
      xr <- NULL
      for (rd in seq_len(nrow(pD))) {
        match <- which(pD[[testCovariate]] %in% pDr[rd, testCovariate])
        taken <- which(match %in% xr)
        if (length(taken) > 0) match <- match[-taken]
        if (length(match) > 0) xr <- c(xr, match[1])
      }
      pDr[, adjustCovariate] <- pD[xr, adjustCovariate]
    }
    if (!is.null(mC)) {
      c1 <- pDr[designr[, coeff] == 1, mC]
      c2 <- pDr[designr[, coeff] == 0, mC]
      if (!(sum(c1 %in% c2) > 0 && sum(c2 %in% c1) > 0)) {
        if (verbose) message("Skipping permutation ", j)
        next
      }
    }
    out_design[[length(out_design) + 1L]] <- unname(designr)
    out_X[[length(out_X) + 1L]] <- unname(.fit_covariates(pDr, testCovariate, adjustCovariate))
  }
  labels <- vapply(out_design, function(d) paste0(d[, coeff[1]], collapse = if (twoGroup) "" else ","), "")
  list(design = out_design, X = out_X, labels = labels, betaNames = betaNames)
}

# Pooled null, empirical p-values, Benjamini-Hochberg q-values, GRanges
.finish <- function(OBS, FLIP, stat, betaNames, fit, seed, labels) {
  columns <- c("chr", "start", "end", "indexStart", "indexEnd", "L", "area", betaNames, "stat")
  OBS <- OBS[, columns, drop = FALSE]
  if (nrow(FLIP) == 0L) {
    warning("No candidate regions found in permutation, so inference can't be carried out. ",
            "Try decreasing the cutoff, or running on a larger dataset if you are currently using a subset.")
    OBS$pval <- NA_real_
    OBS$qval <- NA_real_
  } else {
    if (nrow(FLIP) < 0.05 * nrow(OBS)) {
      message("Note: Very few null candidate regions were found. For more accurate and sensitive inference, ",
              "try decreasing the cutoff, or running on a larger dataset if you are currently using a subset.")
    }
    if (nrow(FLIP) > 1e+06) FLIP <- FLIP[sample.int(nrow(FLIP), 1e+06), , drop = FALSE]
    if (stat == "avg") {
      OBS$avg <- OBS$area / OBS$L
      FLIP$avg <- FLIP$area / FLIP$L
    }
    nullStats <- sort(abs(FLIP[[stat]]))
    observed <- abs(OBS[[stat]])
    # (1 + number of null statistics at least as large) / (1 + number of null statistics)
    pval <- (1 + length(nullStats) - findInterval(observed, nullStats, left.open = TRUE)) / (1 + length(nullStats))
    OBS$pval <- pval
    OBS$qval <- stats::p.adjust(pval, method = "BH")
    OBS$avg <- NULL
  }
  index <- IRanges::IRanges(OBS$indexStart, OBS$indexEnd)
  gr <- GenomicRanges::makeGRangesFromDataFrame(OBS[, setdiff(colnames(OBS), c("indexStart", "indexEnd")), drop = FALSE],
                                                keep.extra.columns = TRUE)
  gr$index <- index
  names(gr) <- NULL
  S4Vectors::metadata(gr)$dmrseqx <- list(version = as.character(utils::packageVersion("dmrseqx")),
                                          fit = fit, labels = labels, seed = seed,
                                          nullStatistics = nrow(FLIP))
  # the regions of the permutations, for checks
  if (isTRUE(getOption("dmrseqx.keepNull"))) S4Vectors::metadata(gr)$null <- FLIP
  gr[order(gr$pval, -abs(gr$stat))]
}

# qt(0.995, df) for the slope tests of the trimming of regions
.t_critical <- function() stats::qt(0.995, seq_len(20000L))

# qt(1 - 0.01 / n, n - 2) for the test for one outlier, n = 3, ..., 10
.grubbs_critical <- function() vapply(3:10, function(n) stats::qt(1 - 0.01 / n, n - 2), numeric(1))

# The fit of one region: for tests and checks
.fit_region <- function(M, Cov, pos, X, weighted = TRUE, fit = c("dmrseq", "global"), seed = 1, nTest = 1L) {
  fit <- match.arg(fit)
  storage.mode(M) <- "integer"
  storage.mode(Cov) <- "integer"
  params <- as.double(c(0.1, 1000, 30, 3000, 2500, 1000, 20, 5, 1, as.numeric(fit == "global"), 0, 1, nTest, nTest + 1))
  out <- .Call(C_fit_one, M, Cov, as.integer(pos), matrix(as.double(X), nrow = ncol(M)), isTRUE(weighted),
               params, .t_critical(), .grubbs_critical(), as.double(seed))
  stats::setNames(out, c("stat", "phi", "constant", if (nTest > 1) paste0("beta", seq_len(nTest)) else "beta"))
}

# The smoother on one block: for tests and checks
.smooth_block <- function(x, y, w, nn, h = 1000, cut = 0.8, at = x) {
  .Call(C_smooth_one, as.double(x), as.double(y), as.double(w), as.double(nn), as.double(h),
        as.double(cut), as.double(at))
}
