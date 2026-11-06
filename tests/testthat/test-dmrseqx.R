simulate_bs <- function(n1 = 3, n2 = 3, L = 3000, missing = 0, seed = 1) {
  set.seed(seed)
  pos <- cumsum(c(1e5, rpois(L - 1, 60) + 2))
  chr <- rep(c("chr1", "chr2"), each = L / 2)
  pos[chr == "chr2"] <- pos[chr == "chr2"] - pos[L / 2]
  N <- n1 + n2
  group <- rep(c("A", "B"), c(n1, n2))
  base <- rep(rbeta(L / 25 + 1, 2, 2), each = 25, length.out = L)
  effect <- rep(0, L)
  for (s in seq(100, L - 200, by = 280)) effect[s:(s + 40)] <- 0.35 * ifelse(base[s] > 0.5, -1, 1)
  M <- Cov <- matrix(0L, L, N)
  for (j in seq_len(N)) {
    Cov[, j] <- rpois(L, 15) + 1L
    if (missing > 0) Cov[runif(L) < missing, j] <- 0L
    M[, j] <- rbinom(L, Cov[, j], pmin(pmax(base + effect * (group[j] == "B"), 0.02), 0.98))
  }
  ids <- paste0("s", seq_len(N))
  bsseq::BSseq(M = M, Cov = Cov, pos = as.integer(pos), chr = chr, sampleNames = ids,
               pData = data.frame(group = factor(group, levels = c("A", "B")),
                                  batch = rep(c("x", "y"), length.out = N), row.names = ids))
}

test_that("the regions have the columns of dmrseq and find the simulated differences", {
  bs <- simulate_bs()
  r <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE))
  expect_s4_class(r, "GRanges")
  expect_true(all(c("L", "area", "beta", "stat", "pval", "qval", "index") %in% colnames(S4Vectors::mcols(r))))
  expect_gt(length(r), 5)
  expect_false(is.unsorted(r$pval, na.rm = TRUE))
  expect_equal(r$L, BiocGenerics::width(r$index))
  expect_equal(GenomicRanges::start(bs)[BiocGenerics::start(r$index)], GenomicRanges::start(r))
  expect_true(all(r$pval > 0 & r$pval <= 1))
})

# without loci that have no coverage in a group
covered <- function(bs) {
  cv <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
  g <- bsseq::pData(bs)$group
  bs[rowSums(cv[, g == "A"]) > 0 & rowSums(cv[, g == "B"]) > 0 & rowSums(cv == 0) < ncol(bs) - 1]
}

test_that("results don't depend on the number of threads and follow the seed", {
  bs <- covered(simulate_bs(missing = 0.1, seed = 2))
  set.seed(5)
  a <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE, threads = 1L))
  set.seed(5)
  b <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE, threads = 3L))
  expect_identical(as.data.frame(a), as.data.frame(b))
  set.seed(5)
  d <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE, threads = 2L))
  expect_identical(as.data.frame(a), as.data.frame(d))
})

test_that("given permutations are used", {
  bs <- simulate_bs()
  labels <- cbind(c(1, 0, 1, 0, 1, 0), c(0, 1, 1, 0, 0, 1))
  r <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE, perms = labels))
  expect_equal(S4Vectors::metadata(r)$dmrseqx$labels, c("000111", "101010", "011001"))
  expect_error(dmrseqx(bs, testCovariate = "group", verbose = FALSE, perms = labels[1:3, ]), "one row per sample")
})

test_that("agrees with dmrseq", {
  skip_if_not_installed("dmrseq")
  skip_on_cran()
  for (cfg in list(list(n1 = 3, n2 = 3, missing = 0, seed = 3),
                   list(n1 = 2, n2 = 3, missing = 0.1, seed = 4))) {
    bs <- covered(do.call(simulate_bs, cfg))
    a <- suppressWarnings(suppressMessages(dmrseq::dmrseq(bs, testCovariate = "group", verbose = FALSE,
                                                          BPPARAM = BiocParallel::SerialParam())))
    b <- suppressMessages(dmrseqx(bs, testCovariate = "group", verbose = FALSE))
    A <- as.data.frame(a)
    B <- as.data.frame(b)
    key <- function(x) paste(x$seqnames, x$start, x$end)
    expect_setequal(key(A), key(B))
    o <- match(key(A), key(B))
    expect_equal(A$area, B$area[o], tolerance = 1e-6)
    long <- A$L >= 40
    if (any(long)) expect_equal(A$stat[long], B$stat[o][long], tolerance = 1e-3)
    expect_gt(mean(abs(A$stat - B$stat[o]) < 1e-3, na.rm = TRUE), 0.9)
  }
})

# Three groups of three samples, with an order of the levels that isn't alphabetical
three_groups <- function(seed = 6, L = 3000) {
  bs <- simulate_bs(n1 = 3, n2 = 6, L = L, seed = seed)
  bsseq::pData(bs)$group <- factor(rep(c("c", "a", "b"), each = 3), levels = c("c", "a", "b"))
  bs
}

# Eight samples along a continuous covariate, in decades around 50 years: the effect of a CpG is
# the change of methylation between the values 0 and 1 of the covariate
along_age <- function(seed = 7, L = 3000) {
  bs <- simulate_bs(n1 = 4, n2 = 4, L = L, seed = seed)
  bsseq::pData(bs)$age <- (c(31, 35, 42, 46, 52, 60, 66, 71) - 50) / 10
  bs
}

test_that("several groups: one coefficient per level, statistics without direction", {
  bs <- three_groups()
  perms <- cbind(c(2, 1, 5, 4, 3, 6, 9, 8, 7), 9:1)
  r <- suppressMessages(dmrseqx(bs, testCovariate = "group", cutoff = 0.3, perms = perms, verbose = FALSE))
  expect_s4_class(r, "GRanges")
  expect_true(all(c("L", "area", "beta_a", "beta_b", "stat", "pval", "qval", "index") %in% colnames(S4Vectors::mcols(r))))
  expect_false("beta" %in% colnames(S4Vectors::mcols(r)))
  expect_true(all(r$stat >= 0, na.rm = TRUE))
  expect_gt(length(r), 5)
  expect_length(S4Vectors::metadata(r)$dmrseqx$labels, 3)
  expect_true(all(r$pval > 0 & r$pval <= 1))
})

test_that("several groups: refused input, and permutations follow the seed and not the threads", {
  bs <- three_groups()
  expect_error(dmrseqx(bs, testCovariate = "group", stat = "beta", verbose = FALSE), "several coefficients")
  expect_error(dmrseqx(bs, testCovariate = "group", matchCovariate = "batch", verbose = FALSE), "2-group")
  expect_error(dmrseqx(bs, testCovariate = "group", perms = cbind(rep(1, 9)), verbose = FALSE), "permutation of the samples")
  set.seed(5)
  a <- suppressMessages(dmrseqx(bs, testCovariate = "group", cutoff = 0.3, maxPerms = 2, threads = 1L, verbose = FALSE))
  set.seed(5)
  b <- suppressMessages(dmrseqx(bs, testCovariate = "group", cutoff = 0.3, maxPerms = 2, threads = 3L, verbose = FALSE))
  expect_identical(as.data.frame(a), as.data.frame(b))
  expect_length(S4Vectors::metadata(a)$dmrseqx$labels, 3)
})

test_that("several groups: a permutation may leave a group without coverage at a CpG", {
  bs <- three_groups(seed = 10, L = 1500)
  cov0 <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
  m0 <- as.matrix(bsseq::getCoverage(bs, type = "M"))
  # One sample of each group lacks these CpGs, inside a simulated difference, which the observed
  # labels allow; this permutation gives the three samples the same label
  cov0[675:680, c(1, 4, 7)] <- 0L
  m0[675:680, c(1, 4, 7)] <- 0L
  thin <- bsseq::BSseq(M = m0, Cov = cov0, gr = GenomicRanges::granges(bs), pData = bsseq::pData(bs), sampleNames = bsseq::sampleNames(bs))
  old <- options(dmrseqx.keepNull = TRUE)
  on.exit(options(old))
  r <- suppressMessages(dmrseqx(thin, testCovariate = "group", cutoff = 0.15, perms = cbind(c(4, 1, 2, 5, 3, 7, 6, 8, 9)), verbose = FALSE))
  null <- as.data.frame(S4Vectors::metadata(r)$null)
  expect_gt(length(r), 5)
  expect_false(anyNA(r$stat))
  expect_false(anyNA(null$stat))
  # the region of the permutation passes over these CpGs
  expect_equal(sum(null$indexStart <= 675 & null$indexEnd >= 680), 1L)
})

test_that("a continuous covariate: the columns of dmrseq, and a coefficient with a sign", {
  bs <- along_age()
  perms <- cbind(c(2, 1, 4, 3, 6, 5, 8, 7), 8:1)
  r <- suppressMessages(dmrseqx(bs, testCovariate = "age", cutoff = 0.03, perms = perms, verbose = FALSE))
  expect_true(all(c("L", "area", "beta", "stat", "pval", "qval", "index") %in% colnames(S4Vectors::mcols(r))))
  expect_gt(length(r), 5)
  expect_true(all(sign(r$stat) == sign(r$beta), na.rm = TRUE))
  expect_true(any(r$stat < 0, na.rm = TRUE) && any(r$stat > 0, na.rm = TRUE))
})

test_that("a continuous covariate: refused input", {
  bs <- along_age(L = 1000)
  cov0 <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
  m0 <- as.matrix(bsseq::getCoverage(bs, type = "M"))
  cov0[5, 1:7] <- 0L
  m0[5, 1:7] <- 0L
  thin <- bsseq::BSseq(M = m0, Cov = cov0, gr = GenomicRanges::granges(bs), pData = bsseq::pData(bs), sampleNames = bsseq::sampleNames(bs))
  expect_error(dmrseqx(thin, testCovariate = "age", verbose = FALSE), "zero coverage in at least")
  expect_error(dmrseqx(bs, testCovariate = "age", matchCovariate = "batch", verbose = FALSE), "2-group")
  expect_error(dmrseqx(bs, testCovariate = "age", perms = cbind(c(1, 1, 2, 3, 4, 5, 6, 7)), verbose = FALSE), "permutation of the samples")
})

test_that("several groups and a continuous covariate agree with dmrseq", {
  skip_if_not_installed("dmrseq")
  skip_on_cran()
  cases <- list(list(bs = three_groups(seed = 8, L = 2000), test = "group", cutoff = 0.3),
                list(bs = along_age(seed = 9, L = 2000), test = "age", cutoff = 0.03))
  for (cs in cases) {
    a <- suppressWarnings(suppressMessages(dmrseq::dmrseq(cs$bs, testCovariate = cs$test, cutoff = cs$cutoff, maxPerms = 2,
                                                          verbose = FALSE, BPPARAM = BiocParallel::SerialParam())))
    b <- suppressMessages(dmrseqx(cs$bs, testCovariate = cs$test, cutoff = cs$cutoff, maxPerms = 2, verbose = FALSE))
    A <- as.data.frame(a)
    B <- as.data.frame(b)
    key <- function(x) paste(x$seqnames, x$start, x$end)
    expect_setequal(key(A), key(B))
    o <- match(key(A), key(B))
    expect_equal(A$area, B$area[o], tolerance = 1e-6)
    long <- A$L >= 40
    if (any(long)) expect_equal(A$stat[long], B$stat[o][long], tolerance = 1e-3)
    expect_gt(mean(abs(A$stat - B$stat[o]) < 1e-3, na.rm = TRUE), 0.9)
  }
})

test_that("unsupported input is refused", {
  bs <- simulate_bs()
  one <- bs
  bsseq::pData(one)$group <- "A"
  expect_error(dmrseqx(one, testCovariate = "group", verbose = FALSE), "only one value")
  many <- bs
  bsseq::pData(many)$group <- c("A", "A", "B", "B", "C", "D")
  expect_error(dmrseqx(many, testCovariate = "group", verbose = FALSE, perms = cbind(1:5)), "one row per sample")
  expect_error(dmrseqx(bs, testCovariate = "nothing", verbose = FALSE), "not found")
  expect_error(dmrseqx(bs, testCovariate = "group", cutoff = 2, verbose = FALSE), "cutoff")
  expect_error(dmrseqx(bs, testCovariate = "group", minNumRegion = 2, verbose = FALSE), "minNumRegion")
  empty <- bs
  cov0 <- bsseq::getCoverage(empty, type = "Cov")
  m0 <- bsseq::getCoverage(empty, type = "M")
  cov0[5, 1:3] <- 0L
  m0[5, 1:3] <- 0L
  empty <- bsseq::BSseq(M = as.matrix(m0), Cov = as.matrix(cov0), gr = GenomicRanges::granges(bs),
                        pData = bsseq::pData(bs), sampleNames = bsseq::sampleNames(bs))
  expect_error(dmrseqx(empty, testCovariate = "group", verbose = FALSE), "zero coverage")
})
