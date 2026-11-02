sim_car <- function(n, seed) {
  set.seed(seed)
  S <- matrix(c(1, 0.3, 0.3, 0.6), 2)
  W <- plogis(MASS::mvrnorm(n, c(0.6, -0.8), S))
  X <- runif(n, 0.05, 0.95)
  list(X = X, Y = X * W[, 1] + (1 - X) * W[, 2], W = W)
}

test_that("Gibbs parametric CAR recovers in-sample means under a true CAR model", {
  s <- sim_car(400, 11)
  f <- eiDPM(s$X, s$Y, K = 1, n.draws = 1500, burnin = 500, seed = 2)
  a <- dpmAggregates(f)
  expect_lt(abs(mean(a[, "mean_W1"]) - mean(s$W[, 1])), 0.03)
  expect_lt(abs(mean(a[, "mean_W2"]) - mean(s$W[, 2])), 0.03)
})

test_that("Gibbs DP CAR recovers in-sample means and keeps the truncation mass small", {
  s <- sim_car(400, 12)
  f <- eiDPM(s$X, s$Y, K = 20, n.draws = 1500, burnin = 500, seed = 3)
  a <- dpmAggregates(f)
  expect_lt(abs(mean(a[, "mean_W1"]) - mean(s$W[, 1])), 0.04)
  expect_lt(abs(mean(a[, "mean_W2"]) - mean(s$W[, 2])), 0.04)
  expect_lt(stats::quantile(f$diagnostics$pi_last, 0.99), 0.01)
})

test_that("Gibbs parametric NCAR runs and gives rates on the lines", {
  set.seed(4)
  n <- 300
  Z <- MASS::mvrnorm(n, c(0.5, -0.5, 0), matrix(c(1, .2, .4, .2, 1, -.3, .4, -.3, 1), 3))
  W <- plogis(Z[, 1:2]); X <- plogis(Z[, 3]); Y <- X * W[, 1] + (1 - X) * W[, 2]
  f <- eiDPM(X, Y, context = TRUE, K = 1, n.draws = 1500, burnin = 500, seed = 5)
  a <- dpmAggregates(f)
  expect_lt(abs(mean(a[, "mean_W1"]) - mean(W[, 1])), 0.04)
  expect_lt(abs(mean(a[, "mean_W2"]) - mean(W[, 2])), 0.04)
})

test_that("Gibbs agrees with eco on the reg data (parametric CAR)", {
  skip_on_cran()
  skip_if_not_installed("eco")
  data(reg, package = "eco")
  f <- eiDPM(reg$X, reg$Y, reg$N, K = 1, n.draws = 3000, burnin = 1000, seed = 6)
  set.seed(6)
  e <- eco::eco(Y ~ X, N = N, data = reg, n.draws = 4000, burnin = 1000)
  expect_lt(abs(mean(rowMeans(f$W1)) - mean(rowMeans(e$W[, 1, ]))), 0.015)
  expect_gt(cor(colMeans(f$W1), colMeans(e$W[, 1, ])), 0.98)
})

test_that("Stan parametric CAR agrees with the Gibbs sampler", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if(Sys.getenv("EIDPM_TEST_STAN") != "true", "set EIDPM_TEST_STAN=true to compile and test the Stan model")
  s <- sim_car(150, 13)
  g <- eiDPM(s$X, s$Y, K = 1, n.draws = 2000, burnin = 500, seed = 7)
  h <- eiDPM(s$X, s$Y, K = 1, method = "stan", n.draws = 500, burnin = 500, chains = 2, seed = 7)
  expect_lt(abs(mean(rowMeans(g$W1)) - mean(rowMeans(h$W1))), 0.02)
  expect_equal(h$diagnostics$divergences, 0)
})

test_that("exact DP (collapsed) sampler recovers in-sample means and agrees with the blocked sampler", {
  s <- sim_car(300, 21)
  f <- eiDPM(s$X, s$Y, n.draws = 3000, burnin = 500, seed = 8)          # K = Inf
  expect_true(is.infinite(f$K))
  a <- dpmAggregates(f)
  expect_lt(abs(mean(a[, "mean_W1"]) - mean(s$W[, 1])), 0.04)
  expect_lt(abs(mean(a[, "mean_W2"]) - mean(s$W[, 2])), 0.04)
  expect_true(all(f$diagnostics$n_clusters >= 1))
  b <- eiDPM(s$X, s$Y, K = 30, n.draws = 3000, burnin = 500, seed = 9)
  expect_gt(cor(colMeans(f$W1), colMeans(b$W1)), 0.98)
})

test_that("exact DP sampler handles fully and partially observed units (eco's rule)", {
  set.seed(5)
  X <- c(runif(40, 0.1, 0.9), 1, 0, 0.5); Y <- c(runif(40, 0.2, 0.8), 0.7, 0.4, 1)
  f <- eiDPM(X, Y, n.draws = 300, burnin = 100, seed = 1)
  expect_equal(unname(f$W1[1, 41]), 0.7); expect_equal(unname(f$W2[1, 42]), 0.4)
  expect_equal(unname(f$W1[1, 43]), 1 - 1e-4)
  dev <- abs(sweep(f$W1[, 1:40], 2, X[1:40], `*`) + sweep(f$W2[, 1:40], 2, 1 - X[1:40], `*`) - rep(Y[1:40], each = 300))
  expect_lt(max(dev), 1e-10)
  expect_error(eiDPM(X, Y, method = "stan"), "finite K")
})
