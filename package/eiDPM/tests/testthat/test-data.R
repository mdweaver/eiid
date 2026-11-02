test_that("dpmData applies eco's boundary rule", {
  X <- c(0.3, 1, 0, 0.5, 0.5, 0.4)
  Y <- c(0.6, 0.7, 0.2, 1, 0, 0.5)
  d <- dpmData(X, Y)
  expect_equal(d$lat, c(1, 6))
  expect_equal(d$x1, 2); expect_equal(d$x0, 3); expect_equal(d$yb, c(4, 5))
  expect_equal(d$Wyb, c(1 - 1e-4, 1e-4))
  expect_equal(unname(d$Wobs[1, ]), qlogis(c(1 - 1e-4, 1e-4)))
  # bounds: Duncan-Davis, shrunk inward
  L <- pmax(0, (X[d$lat] + Y[d$lat] - 1) / X[d$lat]); U <- pmin(1, Y[d$lat] / X[d$lat])
  expect_true(all(d$L > L & d$U < U & d$L < d$U))
  expect_error(dpmData(X, Y, context = TRUE), "0 < X < 1")
  expect_error(dpmData(c(0.2, 1.2), c(0.5, 0.5)), "\\[0, 1\\]")
})

test_that("rates stay on the tomography line and inside the bounds", {
  set.seed(3)
  n <- 60; X <- runif(n, 0.1, 0.9)
  W <- plogis(MASS::mvrnorm(n, c(0.5, -0.5), diag(0.5, 2)))
  Y <- X * W[, 1] + (1 - X) * W[, 2]
  f <- eiDPM(X, Y, K = 5, n.draws = 200, burnin = 100, seed = 1)
  g <- eiDPM(X, Y, n.draws = 200, burnin = 100, seed = 1)
  expect_equal(dim(g$W1), c(200, n))
  expect_equal(dim(f$W1), c(200, n))
  dev <- abs(sweep(f$W1, 2, X, `*`) + sweep(f$W2, 2, 1 - X, `*`) - rep(Y, each = 200))
  expect_lt(max(dev), 1e-10)
  L <- pmax(0, (X + Y - 1) / X); U <- pmin(1, Y / X)
  expect_true(all(sweep(f$W1, 2, L, `>=`) & sweep(f$W1, 2, U, `<=`)))
})
