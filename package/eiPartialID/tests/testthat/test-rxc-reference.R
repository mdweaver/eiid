# Reference values computed independently in Python (docs/extension/verify_phase2b.py machinery)
# for a noise-free 3x3 example with 6 groups (tests/testthat/rxc_example_3x3.csv).
ref <- list(
  loy = list(rc = c(1, 2), M0 = c(0.405474371136, 0.930814566222), M2 = c(0.357607056303, 0.935056634124),
             M1 = c(0.419172246347, 0.894376343064), M1delta = c(0.382661680118, 0.933029434672)),
  def = list(rc = c(1, 3), M0 = c(0.031646907934, 0.355394880682), M2 = c(0.022321546465, 0.42021825439),
             M1 = c(0.046447984358, 0.352629674547), M1delta = c(0.022321546465, 0.397786442847)),
  mob = list(rc = c(3, 2), M0 = c(0.061561650806, 0.455298411589), M2 = c(0.055724651297, 0.523017645223),
             M1 = c(0.081640431546, 0.453034510662), M1delta = c(0.055938558834, 0.48251977772)))

d <- utils::read.csv(test_path("rxc_example_3x3.csv"))
X <- as.matrix(d[, c("X1", "X2", "X3")]); T <- as.matrix(d[, c("T1", "T2", "T3")])

test_that("plug-in bounds match the independent Python implementation (M0, M2, M1, M1delta)", {
  for (m in c("M0", "M2", "M1", "M1delta")) {
    tg <- lapply(ref, function(r) r$rc)
    fit <- eiBounds(x = X, t = T, n = d$N, group = d$group, model = m, delta = 0.05,
                    targets = tg, ci = "none", tests = FALSE)
    for (k in names(ref)) {
      row <- fit$table[fit$table$target == k, ]
      expect_equal(c(row$estimate_lower, row$estimate_upper), ref[[k]][[m]], tolerance = 1e-9,
                   label = paste(m, k))
    }
  }
})

test_that("the ladder is nested: M2 contains M1(delta) contains M1", {
  tg <- list(loy = c(1, 2), mob = c(3, 2))
  fits <- list(
    M2 = eiBounds(x = X, t = T, n = d$N, group = d$group, model = "M2", targets = tg, ci = "none", tests = FALSE),
    d10 = eiBounds(x = X, t = T, n = d$N, group = d$group, model = "M1delta", delta = 0.10, targets = tg, ci = "none", tests = FALSE),
    d02 = eiBounds(x = X, t = T, n = d$N, group = d$group, model = "M1delta", delta = 0.02, targets = tg, ci = "none", tests = FALSE),
    M1 = eiBounds(x = X, t = T, n = d$N, group = d$group, model = "M1", targets = tg, ci = "none", tests = FALSE))
  for (i in 2:length(fits)) {
    a <- fits[[i - 1]]$table; b <- fits[[i]]$table
    expect_true(all(b$estimate_lower >= a$estimate_lower - 1e-9))
    expect_true(all(b$estimate_upper <= a$estimate_upper + 1e-9))
  }
})

test_that("the test of M1's testable part does not reject when M1 holds and rejects when it fails", {
  set.seed(8)
  G <- 8; ng <- 60
  g <- rep(seq_len(G), each = ng)
  x <- stats::runif(G * ng, 0.05, 0.95)
  shift <- stats::rnorm(G, 0, 0.05)[g]
  make <- function(slope_w) {
    bw <- 0.55 + shift + slope_w * x; bb <- 0.35 + shift + 0.20 * x
    t1 <- x * bb + (1 - x) * bw + stats::rnorm(length(x), 0, 0.01)
    cbind(t1, 1 - t1)
  }
  ok <- eiBounds(x = cbind(x, 1 - x), t = make(-0.10), n = rep(100, length(x)), group = g, model = "M1",
                 targets = list(a = c(1, 1)), ci = "none", tests = TRUE)
  expect_gt(ok$tests$m1$p.value, 0.01)
  slopes <- seq(-0.6, 0.6, length.out = G)[g]                  # contextual effect differs by group
  bad <- eiBounds(x = cbind(x, 1 - x), t = make(slopes), n = rep(100, length(x)), group = g, model = "M1",
                  targets = list(a = c(1, 1)), ci = "none", tests = TRUE)
  expect_lt(bad$tests$m1$p.value, 0.01)
})
