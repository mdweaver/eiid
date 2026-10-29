sim2x3 <- function(n = 400, seed = 1, conc = 50) {
  set.seed(seed)
  x <- stats::runif(n, 0.05, 0.95)
  X <- cbind(x, 1 - x)
  m1 <- cbind(0.70 + 0.10 * x, 0.10 - 0.05 * x, 0.20 - 0.05 * x)      # row 1 means (sum 1)
  m2 <- cbind(0.20 + 0.05 * x, 0.65 - 0.10 * x, 0.15 + 0.05 * x)      # row 2 means
  draw <- function(M) t(apply(M, 1, function(m) { g <- stats::rgamma(3, conc * m); g / sum(g) }))
  b1 <- draw(m1); b2 <- draw(m2)
  T <- X[, 1] * b1 + X[, 2] * b2
  N <- rep(500, n)
  list(X = X, T = T, N = N, b1 = b1, m1 = m1)
}

test_that("the three input interfaces agree", {
  s <- sim2x3()
  df <- data.frame(r1 = s$X[, 1] * s$N, r2 = s$X[, 2] * s$N,
                   c1 = s$T[, 1] * s$N, c2 = s$T[, 2] * s$N, c3 = s$T[, 3] * s$N)
  a <- eiBounds(cbind(c1, c2, c3) ~ cbind(r1, r2), data = df, ci = "none", tests = FALSE)
  b <- eiBounds(data = df, rows = c("r1", "r2"), cols = c("c1", "c2", "c3"), ci = "none", tests = FALSE)
  cc <- eiBounds(x = s$X, t = s$T, n = s$N, ci = "none", tests = FALSE)
  expect_equal(a$table$estimate_lower, b$table$estimate_lower)
  expect_equal(a$table$estimate_lower, cc$table$estimate_lower)
  expect_equal(a$table$estimate_upper, cc$table$estimate_upper)
})

test_that("mismatched totals are refused unless rescaling is requested", {
  df <- data.frame(r1 = c(10, 20, 30, 5), r2 = c(10, 20, 30, 15), c1 = c(5, 10, 20, 8), c2 = c(14, 30, 40, 12))
  expect_error(eiBounds(data = df, rows = c("r1", "r2"), cols = c("c1", "c2"), ci = "none", tests = FALSE),
               "totals differ")
})

test_that("extrapolating to pure units reproduces the closed form (Phase 2a Prop. 5)", {
  s <- sim2x3(seed = 2)
  fit <- suppressWarnings(eiBounds(x = s$X, t = s$T, n = s$N, range = "extrapolate", ci = "none", tests = FALSE,
                                   targets = list(a = c(1, 1), b = c(1, 2))))
  Z <- cbind(s$X[, 1]^2, s$X[, 1] * s$X[, 2], s$X[, 2]^2)
  w <- s$N * s$X[, 1] / sum(s$N * s$X[, 1])
  mu11 <- sum(w * s$X[, 1]); mu12 <- sum(w * s$X[, 2])
  for (k in 1:2) {
    b <- stats::lm.fit(Z, s$T[, k])$coefficients
    Q11 <- b[1]; Q12 <- b[2] / 2
    mm <- min(Q12, 1 - Q12)
    expect_equal(fit$table$estimate_lower[k], unname(Q11 * mu11 + (Q12 - mm) * mu12), tolerance = 1e-9)
    expect_equal(fit$table$estimate_upper[k], unname(Q11 * mu11 + (Q12 + mm) * mu12), tolerance = 1e-9)
  }
})

test_that("Goodman (L6) equals the no-intercept regression", {
  s <- sim2x3(seed = 3)
  fit <- eiBounds(x = s$X, t = s$T, n = s$N, model = "L6", targets = list(a = c(1, 1)), tests = FALSE, seed = 1)
  g <- stats::lm.fit(s$X, s$T[, 1])$coefficients
  expect_equal(fit$table$estimate_lower, unname(g[1]), tolerance = 1e-10)
  expect_equal(fit$table$estimate_upper, unname(g[1]), tolerance = 1e-10)
})

test_that("restrictions can only narrow the bound", {
  s <- sim2x3(seed = 4)
  base <- eiBounds(x = s$X, t = s$T, n = s$N, targets = list(a = c(1, 1)), ci = "none", tests = FALSE)
  o <- eiBounds(x = s$X, t = s$T, n = s$N, targets = list(a = c(1, 1)), ci = "none", tests = FALSE,
                restrictions = list(order = list(c(high = 1, low = 2, col = 1)),
                                    slope = 0.5,
                                    external = list(list(target = c(1, 1), lower = 0.70, upper = 0.80))))
  expect_true(o$table$estimate_lower >= base$table$estimate_lower - 1e-12)
  expect_true(o$table$estimate_upper <= base$table$estimate_upper + 1e-12)
  expect_true(o$table$estimate_upper <= 0.80 + 1e-12)
})

test_that("a structural zero that holds is compatible with the data", {
  set.seed(5)
  x <- stats::runif(300, 0.05, 0.95); X <- cbind(x, 1 - x)
  m1 <- cbind(0.7 + 0.1 * x, 0.3 - 0.1 * x, 0)          # row 1 never goes to column 3
  m2 <- cbind(0.2, 0.5, 0.3)[rep(1, 300), ]
  T <- X[, 1] * m1 + X[, 2] * m2
  fit <- eiBounds(x = X, t = T, n = rep(100, 300), ci = "none", tests = FALSE,
                  targets = list(z = c(1, 3)), restrictions = list(zeros = list(c(row = 1, col = 3))))
  expect_false(fit$table$empty)
  expect_equal(fit$table$estimate_upper, 0, tolerance = 1e-8)
})

test_that("the variance falls back from HC3 when leverage is one", {
  set.seed(6)
  g <- rep(1:3, c(40, 40, 3))
  x <- stats::runif(length(g), 0.1, 0.9); X <- cbind(x, 1 - x)
  t1 <- 0.5 + 0.1 * x + stats::rnorm(length(g), 0, 0.02)
  fit <- eiBounds(x = X, t = cbind(t1, 1 - t1), n = rep(100, length(g)), group = g, model = "M2",
                  targets = list(a = c(1, 1)), ci = "none", tests = FALSE)
  expect_equal(fit$vcov_type, "HC1")
  expect_true(any(grepl("fell back", fit$notes)))
})

test_that("confidence intervals cover the truth in a small simulation", {
  cover_model <- cover_actual <- logical(0)
  for (k in 1:8) {
    s <- sim2x3(n = 300, seed = 100 + k)
    fit <- eiBounds(x = s$X, t = s$T, n = s$N, targets = list(a = c(1, 1)), tests = FALSE, seed = k,
                    ndraw = 40, nz_draw = 1500)
    w <- s$N * s$X[, 1] / sum(s$N * s$X[, 1])
    truth_model <- sum(w * s$m1[, 1]); truth_actual <- sum(w * s$b1[, 1])
    cover_model <- c(cover_model, fit$table$ci_model_lower <= truth_model && truth_model <= fit$table$ci_model_upper)
    cover_actual <- c(cover_actual, fit$table$ci_actual_lower <= truth_actual && truth_actual <= fit$table$ci_actual_upper)
  }
  expect_true(mean(cover_model) >= 0.75)
  expect_true(all(cover_actual))
})

test_that("sensitivity bounds widen with eta and breakdown is monotone", {
  s <- sim2x3(seed = 7)
  fit <- eiBounds(x = s$X, t = s$T, n = s$N, targets = list(a = c(1, 1)), tests = FALSE, seed = 1,
                  ndraw = 30, nz_draw = 1000)
  sv <- sensitivity(fit, eta = c(0, 0.02, 0.05))
  expect_true(all(diff(sv$estimate_upper) >= -1e-12))
  expect_true(all(diff(sv$estimate_lower) <= 1e-12))
  e1 <- breakdown(fit, "a", threshold = 0.5, which = "estimate")
  e2 <- breakdown(fit, "a", threshold = 0.6, which = "estimate")
  expect_true(e1 >= e2)
})
