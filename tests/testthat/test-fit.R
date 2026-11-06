simulate_region <- function(L, N = 8, seed = 1, gap = NULL, zero = 0, effect = 0.3) {
  set.seed(seed)
  pos <- cumsum(c(1000, if (is.null(gap)) rpois(L - 1, 40) + 2 else rep(gap, L - 1)))
  group <- rep(0:1, length.out = N)
  # correlated between neighbouring CpGs within samples
  Cov <- matrix(rpois(L * N, 18) + 1L, L, N)
  if (zero > 0) Cov[runif(L * N) < zero] <- 0L
  M <- Cov
  for (s in seq_len(N)) {
    e <- as.numeric(stats::arima.sim(list(ar = 0.6), L)) * 0.08
    p <- pmin(pmax(0.5 + effect * group[s] + e, 0.02), 0.98)
    M[, s] <- rbinom(L, Cov[, s], p)
  }
  list(M = M, Cov = Cov, pos = as.integer(pos), group = group)
}

gls_stat <- function(r, adjust = NULL, weighted = TRUE) {
  L <- nrow(r$M)
  N <- ncol(r$M)
  dat <- data.frame(g.fac = factor(rep(r$group, each = L)), sample = factor(rep(seq_len(N), each = L)),
                    meth = as.vector(r$M), cov = as.vector(r$Cov), L = rep(as.numeric(r$pos), N))
  if (!is.null(adjust)) dat$adj <- factor(rep(adjust, each = L))
  dat$pos <- as.numeric(factor(dat$L))
  dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, stats::median)), N)
  dat <- dat[dat$cov > 0, ]
  dat$Z <- asin(2 * (dat$meth + 0.05) / (dat$cov + 0.1) - 1)
  corr <- if (L < 40) nlme::corCAR1(form = ~L | sample) else nlme::corAR1(form = ~1 | sample)
  w <- if (weighted && all(dat$MedCov > 0)) nlme::varPower(form = ~1 / MedCov, fixed = 0.5) else NULL
  frm <- if (is.null(adjust)) Z ~ g.fac + factor(L) else Z ~ g.fac + factor(L) + adj
  fit <- summary(nlme::gls(frm, weights = w, data = dat, correlation = corr))
  c(beta = unname(fit$tTable[2, 1]), stat = unname(fit$tTable[2, 3]),
    phi = unname(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE)))
}

test_that("long regions agree with gls", {
  skip_if_not_installed("nlme")
  for (cfg in list(list(L = 40, seed = 1), list(L = 90, seed = 2), list(L = 60, N = 11, seed = 3))) {
    r <- do.call(simulate_region, cfg)
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group)
    b <- gls_stat(r)
    expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
    expect_equal(a[["beta"]], b[["beta"]], tolerance = 1e-4)
  }
})

test_that("long regions with missing coverage agree with gls", {
  skip_if_not_installed("nlme")
  for (seed in 1:3) {
    r <- simulate_region(70, seed = seed, zero = 0.2)
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group)
    b <- gls_stat(r)
    expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
    expect_equal(a[["beta"]], b[["beta"]], tolerance = 1e-4)
  }
})

test_that("short regions with close CpGs reach the optimum of gls", {
  skip_if_not_installed("nlme")
  for (seed in 1:4) {
    r <- simulate_region(20, seed = seed, gap = 2)
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group)
    b <- gls_stat(r)
    expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
  }
})

test_that("short regions with distant CpGs stay at the start value, as gls does", {
  skip_if_not_installed("nlme")
  for (seed in 1:4) {
    r <- simulate_region(20, seed = seed, gap = 30)
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group)
    b <- gls_stat(r)
    expect_equal(a[["phi"]], 0.2)
    expect_equal(b[["phi"]], 0.2, tolerance = 1e-2)
    expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
  }
})

test_that("fit = global finds the true optimum in the same regions", {
  r <- simulate_region(20, seed = 1, gap = 30)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, fit = "global")
  expect_gt(a[["phi"]], 0.2)
})

test_that("the sparse path gives the optimum of gls with one row missing", {
  skip_if_not_installed("nlme")
  r <- simulate_region(30, seed = 8)
  r$Cov[5, 2] <- 0L
  r$M[5, 2] <- 0L
  gap <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, fit = "global")
  expect_equal(gap[["stat"]], {
    L <- nrow(r$M); N <- ncol(r$M)
    dat <- data.frame(g.fac = factor(rep(r$group, each = L)), sample = factor(rep(seq_len(N), each = L)),
                      meth = as.vector(r$M), cov = as.vector(r$Cov), L = rep(as.numeric(r$pos), N))
    dat$pos <- as.numeric(factor(dat$L))
    dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, stats::median)), N)
    dat <- dat[dat$cov > 0, ]
    dat$Z <- asin(2 * (dat$meth + 0.05) / (dat$cov + 0.1) - 1)
    # the optimum, from a start close to it
    ctl <- nlme::corCAR1(value = gap[["phi"]], form = ~L | sample)
    fit <- summary(nlme::gls(Z ~ g.fac + factor(L), weights = nlme::varPower(form = ~1 / MedCov, fixed = 0.5),
                             data = dat, correlation = ctl))
    unname(fit$tTable[2, 3])
  }, tolerance = 1e-3)
})

test_that("adjustment covariates agree with gls", {
  skip_if_not_installed("nlme")
  r <- simulate_region(50, N = 10, seed = 5)
  batch <- rep(c(0, 1), each = 5)
  X <- cbind(r$group, batch)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, X)
  b <- gls_stat(r, adjust = batch)
  expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
  expect_equal(a[["beta"]], b[["beta"]], tolerance = 1e-4)
})

# The gls of dmrseq for a covariate of interest x of the samples, a factor or numeric.
gls_covariate <- function(r, x) {
  L <- nrow(r$M)
  N <- ncol(r$M)
  dat <- data.frame(x = rep(x, each = L), sample = factor(rep(seq_len(N), each = L)),
                    meth = as.vector(r$M), cov = as.vector(r$Cov), L = rep(as.numeric(r$pos), N))
  dat$pos <- as.numeric(factor(dat$L))
  dat$MedCov <- rep(as.numeric(by(dat$cov, dat$pos, stats::median)), N)
  dat <- dat[dat$cov > 0, ]
  dat$Z <- asin(2 * (dat$meth + 0.05) / (dat$cov + 0.1) - 1)
  corr <- if (L < 40) nlme::corCAR1(form = ~L | sample) else nlme::corAR1(form = ~1 | sample)
  nlme::gls(Z ~ x + factor(L), weights = nlme::varPower(form = ~1 / MedCov, fixed = 0.5), data = dat, correlation = corr)
}

test_that("several groups: the root of the sequential F test of gls and one coefficient per level", {
  skip_if_not_installed("nlme")
  for (cfg in list(list(L = 60, N = 9, seed = 11), list(L = 45, N = 12, seed = 12, zero = 0.15), list(L = 20, N = 9, seed = 13, gap = 2))) {
    r <- do.call(simulate_region, cfg)
    g <- factor(rep(c("a", "b", "c"), length.out = ncol(r$M)))
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, stats::model.matrix(~g)[, -1], nTest = 2L)
    fit <- gls_covariate(r, g)
    expect_equal(a[["stat"]], sqrt(stats::anova(fit)[["F-value"]][2]), tolerance = 1e-3)
    expect_equal(unname(a[c("beta1", "beta2")]), unname(stats::coef(fit)[2:3]), tolerance = 1e-4)
  }
})

test_that("a continuous covariate agrees with gls", {
  skip_if_not_installed("nlme")
  for (cfg in list(list(L = 70, N = 8, seed = 21), list(L = 50, N = 10, seed = 22, zero = 0.15), list(L = 20, N = 8, seed = 23, gap = 2))) {
    r <- do.call(simulate_region, cfg)
    set.seed(cfg$seed)
    x <- round(stats::runif(ncol(r$M), 20, 80), 1)
    a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, x)
    fit <- summary(gls_covariate(r, x))
    expect_equal(a[["stat"]], unname(fit$tTable[2, 3]), tolerance = 1e-3)
    expect_equal(a[["beta"]], unname(fit$tTable[2, 1]), tolerance = 1e-4)
  }
})

test_that("the statistics don't depend on the units of a continuous covariate", {
  r <- simulate_region(50, N = 8, seed = 24)
  x <- c(21, 25, 33, 38, 44, 52, 61, 70)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, x)
  b <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, x / 100)
  expect_equal(a[["stat"]], b[["stat"]])
  expect_equal(a[["beta"]], b[["beta"]] / 100)
})

test_that("the refit without variance weights agrees with gls", {
  skip_if_not_installed("nlme")
  r <- simulate_region(60, seed = 7)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, weighted = FALSE)
  b <- gls_stat(r, weighted = FALSE)
  expect_equal(a[["stat"]], b[["stat"]], tolerance = 1e-3)
})

test_that("a coverage outlier in a short region removes the variance weights", {
  skip_if_not_installed("outliers")
  r <- simulate_region(8, seed = 3, gap = 2)
  r$Cov[4, ] <- 500L
  r$M[4, ] <- rbinom(ncol(r$M), 500L, 0.5)
  expect_lt(outliers::grubbs.test(apply(r$Cov, 1, stats::median))$p.value, 0.01)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, weighted = TRUE)
  b <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, weighted = FALSE)
  expect_equal(a[["stat"]], b[["stat"]])
})

test_that("constant regions and models that can't be fit are recognized", {
  r <- simulate_region(10, seed = 1)
  zero <- r$M * 0L
  expect_equal(dmrseqx:::.fit_region(zero, r$Cov, r$pos, r$group)[["constant"]], 1)
  expect_equal(dmrseqx:::.fit_region(r$Cov, r$Cov, r$pos, r$group)[["constant"]], 1)
  # the covariate of interest is confounded with the loci
  r$Cov[1:5, r$group == 1] <- 0L
  r$Cov[6:10, r$group == 0] <- 0L
  r$M <- pmin(r$M, r$Cov)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group)
  expect_true(is.na(a[["stat"]]))
  expect_equal(a[["constant"]], 0)
})

test_that("zero median coverage leads to the fit without weights", {
  r <- simulate_region(50, seed = 2)
  a <- dmrseqx:::.fit_region(r$M, r$Cov, r$pos, r$group, weighted = FALSE)
  r0 <- r
  # A row where most samples have no coverage: no weights, and a different fit
  r0$Cov[10, 1:6] <- 0L
  r0$M[10, 1:6] <- 0L
  b <- dmrseqx:::.fit_region(r0$M, r0$Cov, r0$pos, r0$group, weighted = TRUE)
  d <- dmrseqx:::.fit_region(r0$M, r0$Cov, r0$pos, r0$group, weighted = FALSE)
  expect_equal(b[["stat"]], d[["stat"]])
  expect_false(is.na(a[["stat"]]))
})
