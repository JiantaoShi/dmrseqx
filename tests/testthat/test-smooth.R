locfit_fit <- function(x, y, w, nn, h = 1000, at = x) {
  sdata <- data.frame(posi = x, yi = y, weightsi = w)
  fit <- suppressWarnings(locfit::locfit(yi ~ locfit::lp(posi, nn = nn, h = h), data = sdata,
                                         weights = weightsi, family = "gaussian", maxk = 10000))
  pp <- suppressWarnings(stats::preplot(fit, where = "data", band = "local",
                                        newdata = data.frame(posi = at)))
  pp$trans(pp$fit)
}

simulate_block <- function(n, seed, spacing = 60) {
  set.seed(seed)
  x <- cumsum(c(1e6, rpois(n - 1, spacing) + 1))
  list(x = as.numeric(x), y = sin(x / 2000) / 4 + rnorm(n, 0, 0.1), w = runif(n, 0.5, 20))
}

test_that("the smoother reproduces locfit on the tree of vertices", {
  skip_if_not_installed("locfit")
  for (seed in 1:3) {
    b <- simulate_block(800, seed)
    nn <- min(3000 / (max(b$x) - min(b$x) + 1), 30 / length(b$x), 0.75)
    expect_equal(dmrseqx:::.smooth_block(b$x, b$y, b$w, nn),
                 locfit_fit(b$x, b$y, b$w, nn), tolerance = 1e-6)
  }
})

test_that("the smoother reproduces locfit with a wide bandwidth and a sparse block", {
  skip_if_not_installed("locfit")
  b <- simulate_block(60, 5, spacing = 700)
  nn <- min(3000 / (max(b$x) - min(b$x) + 1), 30 / length(b$x), 0.75)
  expect_equal(dmrseqx:::.smooth_block(b$x, b$y, b$w, nn),
               locfit_fit(b$x, b$y, b$w, nn), tolerance = 1e-6)
})

test_that("positions outside the fit points are close to locfit, and inside equal", {
  skip_if_not_installed("locfit")
  b <- simulate_block(200, 6)
  inner <- 4:195
  nn <- 0.2
  a <- dmrseqx:::.smooth_block(b$x[inner], b$y[inner], b$w[inner], nn, at = b$x)
  e <- locfit_fit(b$x[inner], b$y[inner], b$w[inner], nn, at = b$x)
  expect_equal(a[inner], e[inner], tolerance = 1e-6)
  expect_lt(max(abs(a - e)), 0.01)
})

test_that("nearly singular local fits follow locfit", {
  skip_if_not_installed("locfit")
  y <- c(0.1, 0.3, -0.2, 0.4, 0.2)
  w <- c(5, 4, 3, 6, 2)
  for (x in list(c(10000, 12000, 14000, 16000, 18000),
                 c(10000, 10010, 12000, 12010, 14000),
                 c(10000, 10002, 10004, 12400, 12402))) {
    nn <- min(3221.396 / (max(x) - min(x) + 1), 30 / 5, 0.75)
    expect_equal(dmrseqx:::.smooth_block(x, y, w, nn), locfit_fit(x, y, w, nn), tolerance = 1e-6)
  }
})
