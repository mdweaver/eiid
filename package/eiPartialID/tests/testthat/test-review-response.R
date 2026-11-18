# Package changes made in response to the external review (docs/extension/REVIEW_RESPONSE_PLAN.md, WS5).

sim_groups <- function(G = 6, ng = 60, zeta_sd = 0, seed = 1, conc = 20) {
  set.seed(seed)
  H <- rbind(c(0.70, 0.45), c(0.30, 0.20))           # column-1 rates (H x)_r
  z <- if (zeta_sd > 0) matrix(stats::runif(2 * G, -zeta_sd, zeta_sd), G) else matrix(0, G, 2)
  z[1, ] <- 0
  g <- rep(seq_len(G), each = ng)
  a <- stats::runif(G, 0.1, 0.5)
  x1 <- unlist(lapply(seq_len(G), function(k) stats::runif(ng, a[k], a[k] + 0.4)))
  X <- cbind(x1, 1 - x1)
  m1 <- H[1, 1] * x1 + H[1, 2] * (1 - x1) + z[g, 1]
  m2 <- H[2, 1] * x1 + H[2, 2] * (1 - x1) + z[g, 2]
  b1 <- stats::rbeta(length(x1), m1 * conc, (1 - m1) * conc)
  b2 <- stats::rbeta(length(x1), m2 * conc, (1 - m2) * conc)
  t1 <- x1 * b1 + (1 - x1) * b2
  N <- sample(200:1000, length(x1), replace = TRUE)
  list(X = X, T = cbind(t1, 1 - t1), N = N, g = g, b1 = b1, m1 = m1)
}

test_that("the cluster argument is gone (no clustered variance; P3 section 0.4)", {
  d <- sim_groups()
  expect_error(eiBounds(x = d$X, t = d$T, n = d$N, model = "M0", cluster = d$g), "unused argument")
})

test_that("CLR fits also report projection intervals, which contain the estimate", {
  d <- sim_groups(seed = 3)
  fit <- eiBounds(x = d$X, t = d$T, n = d$N, model = "M0", seed = 1, targets = list(a = c(1, 1)))
  tb <- fit$table
  expect_true(all(c("proj_model_lower", "proj_model_upper", "proj_actual_lower", "proj_actual_upper") %in% names(tb)))
  expect_true(tb$proj_model_lower <= tb$estimate_lower + 1e-9 && tb$proj_model_upper >= tb$estimate_upper - 1e-9)
  expect_true(tb$proj_actual_lower >= tb$dd_lower - 1e-9 && tb$proj_actual_upper <= tb$dd_upper + 1e-9)
  expect_true(any(grepl("algorithmic", fit$notes)))
  fc <- eiBounds(x = d$X, t = d$T, n = d$N, model = "M0", seed = 1, ci = "conservative", targets = list(a = c(1, 1)))
  expect_true(is.na(fc$table$proj_model_lower))
})

test_that("projection interval for the model-implied rate covers in a small simulation", {
  cover <- logical(0)
  for (s in 1:30) {
    d <- sim_groups(G = 1, ng = 300, seed = 100 + s)
    fit <- eiBounds(x = d$X, t = d$T, n = d$N, model = "M0", seed = s, targets = list(a = c(1, 1)), tests = FALSE)
    w <- d$N * d$X[, 1] / sum(d$N * d$X[, 1])
    truth <- sum(w * d$m1)
    cover <- c(cover, fit$table$proj_model_lower <= truth && truth <= fit$table$proj_model_upper)
  }
  expect_gte(mean(cover), 0.9)
})

test_that("pooling test: P3.6(c), df = (G-1)(C-1)R; holds size under zeta = 0 and rejects otherwise", {
  d <- sim_groups(zeta_sd = 0.08, seed = 7)
  fit <- eiBounds(x = d$X, t = d$T, n = d$N, group = d$g, model = "M1", seed = 1, targets = list(a = c(1, 1)))
  expect_equal(fit$tests$pool$df, (6 - 1) * (2 - 1) * 2)
  expect_equal(fit$tests$m1$df, (6 - 1) * (2 - 1) * 2 * 1 / 2)
  expect_lt(fit$tests$pool$p.value, 0.01)
  expect_true(fit$tests$pool$spanning)
  rej <- vapply(1:40, function(s) {
    d0 <- sim_groups(zeta_sd = 0, seed = 200 + s)
    eiPartialID:::.rxc_pool_test(d0$X, d0$T, d0$N, d0$g)$p.value < 0.05
  }, logical(1))
  expect_lte(mean(rej), 0.15)
})

test_that("rank failures name the offending groups", {
  set.seed(5)
  X <- rbind(t(replicate(40, { v <- stats::rexp(3); v / sum(v) })),
             t(replicate(40, { u <- stats::runif(1); c(u, 1 - u, 0) * 0.9 + c(0, 0, 0.1) })))
  T <- X %*% matrix(c(0.6, 0.4, 0.3, 0.7, 0.5, 0.5), 3, byrow = TRUE)
  g <- rep(c("A", "B"), each = 40)
  expect_error(eiBounds(x = X, t = T, n = rep(500, 80), group = g, model = "M2", ci = "none", tests = FALSE),
               "groups without .*B")
})

test_that("Wald df of M1's testable part is (G-1)(C-1)R(R-1)/2 for a 3x3 table", {
  d <- read.csv(test_path("rxc_example_3x3.csv"))
  X <- as.matrix(d[, c("X1", "X2", "X3")]); T <- as.matrix(d[, c("T1", "T2", "T3")])
  G <- length(unique(d$group))
  tt <- eiPartialID:::.rxc_m1_test(X, T, d$N, d$group)
  expect_equal(tt$df, (G - 1) * 2 * 3)
})
