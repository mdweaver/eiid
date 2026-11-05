# Step 1 close-out: module registry (meta) and the slack dial (S9) for sign / CR modules.

sim <- function(n = 200, seed = 7, ab = c(0.40, 0.20, 0.55, -0.05)) {
  set.seed(seed)
  X <- pmin(pmax(stats::rbeta(n, 2, 2) * 0.8 + 0.1, 0.02), 0.98); N <- sample(300:3000, n, replace = TRUE)
  m0 <- ab[1] + ab[2] * X; m1 <- ab[3] + ab[4] * X
  t <- (1 - X) * stats::rbeta(n, m0 * 15, (1 - m0) * 15) + X * stats::rbeta(n, m1 * 15, (1 - m1) * 15)
  list(x = X, t = t, n = N)
}
Dg <- list(D = matrix(c(1, -1, 0, 0), 2, 2))

test_that("registry: every restriction type has meta; fits report their active modules", {
  m <- eiModules()
  expect_setequal(m$type, c("zeros", "order", "sign", "slope", "linear", "external", "assoc", "cr"))
  expect_true(all(c("source", "object", "target", "refutable", "form", "sharp", "dial", "proof_status") %in% names(m)))
  expect_true(all(grepl("^proved", m$proof_status)))
  s <- sim()
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE, diagnostics = FALSE,
                restrictions = list(slope = 2, cr = list(list(row = 1, slack = 0.01)),
                                    assoc = list(list(type = "dW", row = 1, sign = 1))))
  expect_equal(f$modules$module, c("slope", "assoc[1]", "cr[1]"))
  expect_equal(f$modules$slack, c(NA, 0, 0.01))
  expect_match(f$modules$restriction[3], "slack 0.01")
  expect_error(eiBounds(x = s$x, t = s$t, n = s$n, ci = "none", tests = FALSE,
                        restrictions = list(cr = list(list(row = 1, slack = -1)))), "nonnegative")
})

test_that("slack: 2x2 closed forms shift ER / NM by slack * scale / variance", {
  s <- sim()
  fb <- function(res = NULL) eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE,
                                      diagnostics = FALSE, restrictions = res)
  f0 <- fb(); st <- f0$setup
  J <- unlist(f0$table[1, c("estimate_lower", "estimate_upper")])
  ds <- eiPartialID:::.rxc_descriptive(st, 1, f0$rows, f0$cols)
  dz <- eiPartialID:::.rxc_design_r(st, 1)
  sB <- eiPartialID:::.rxc_assoc_scale(st, "dB", "W1", 1); sW <- eiPartialID:::.rxc_assoc_scale(st, "dW", "W1", 1)
  # scales are the natural weights: E[X(1-X)] and mu Var_1(X) + (1-mu) Var_0(X)
  p <- s$n / sum(s$n); mu <- sum(p * s$x)
  wv <- function(w) { w <- w / sum(w); sum(w * s$x^2) - sum(w * s$x)^2 }
  expect_equal(sB, sum(p * s$x * (1 - s$x))); expect_equal(sW, mu * wv(p * s$x) + (1 - mu) * wv(p * (1 - s$x)))
  for (eta in c(0, 0.02, 0.1)) {
    a <- eta * sB / (dz$mu * (1 - dz$mu)); b <- eta * sW / dz$var
    cases <- list(list("dW", 1, c(J[1], min(J[2], ds$ER + b))), list("dW", -1, c(max(J[1], ds$ER - b), J[2])),
                  list("dB", 1, c(max(J[1], ds$NM - a), J[2])), list("dB", -1, c(J[1], min(J[2], ds$NM + a))))
    for (cs in cases) {
      f <- fb(list(assoc = list(list(type = cs[[1]], row = 1, sign = cs[[2]], slack = eta))))
      if (cs[[3]][1] > cs[[3]][2]) expect_true(f$table$empty) else
        expect_equal(unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])), unname(cs[[3]]), tolerance = 1e-9)
    }
    # CR with slack: one interval [min(NM - a, ER - b), max(ER + b, NM + a)] within J ...
    f <- fb(list(cr = list(list(row = 1, slack = eta))))
    ex <- c(max(J[1], min(ds$NM - a, ds$ER - b)), min(J[2], max(ds$ER + b, ds$NM + a)))
    expect_equal(unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])), unname(ex), tolerance = 1e-9)
    # ... equal to the union of the two slack branches
    cons <- eiPartialID:::.rxc_constraints(f$setup)
    br <- eiPartialID:::.rxc_systems(cons, "branch")
    u <- eiPartialID:::.rxc_sys_interval(br, eiPartialID:::.rxc_target(f$setup, Dg$D), f$setup$beta)$union
    expect_equal(nrow(u), 1L)
    expect_equal(unname(u[1, ]), unname(ex), tolerance = 1e-9)
  }
})

test_that("slack dial widens monotonically (estimate and CIs); breakdown on the slack dial", {
  s <- sim(seed = 3)
  for (ci in c("conservative", "clr")) {
    f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = ci, seed = 1, diagnostics = FALSE,
                  restrictions = list(cr = list(list(row = 1)),
                                      assoc = list(list(type = "dW", row = 1, sign = 1))))
    sv <- sensitivity(f, eta = c(0, 0.01, 0.05, 0.2), dial = "slack")
    expect_equal(sv$slack, c(0, 0.01, 0.05, 0.2))
    expect_true(all(diff(sv$estimate_lower) <= 1e-12) && all(diff(sv$estimate_upper) >= -1e-12))
    expect_equal(unname(unlist(sv[1, c("estimate_lower", "estimate_upper")])),
                 unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])))
    if (ci == "conservative")    # the same draws at each value: exact nesting
      expect_true(all(diff(sv$ci_model_lower) <= 1e-12) && all(diff(sv$ci_model_upper) >= -1e-12))
    # only one module varied
    sv1 <- sensitivity(f, eta = c(0, 0.2), dial = "slack", modules = "assoc[1]")
    expect_true(sv1$estimate_upper[2] >= sv1$estimate_upper[1])
  }
  # breakdown: "D > 0" under CR; how much wrong-sign slack before it fails
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE, diagnostics = FALSE,
                restrictions = list(cr = list(list(row = 1))))
  if (f$table$estimate_lower > 0) {
    b <- breakdown(f, "D", threshold = 0, which = "estimate", dial = "slack", max_eta = 0.5, tol = 1e-4)
    expect_gt(b, 0)
    lo <- function(e) sensitivity(f, eta = e, dial = "slack")$estimate_lower
    expect_gt(lo(b * 0.98), 0); expect_lte(lo(min(0.5, b + 2e-4)), 1e-9)
  }
  expect_error(sensitivity(eiBounds(x = s$x, t = s$t, n = s$n, ci = "none", tests = FALSE), eta = 0.1, dial = "slack"),
               "needs assoc or cr")
})

test_that("3x2 gap case: slack keeps the union structure and only widens", {
  ref <- utils::read.csv(test_path("msign_reference_3x2.csv"))
  gp <- utils::read.csv(test_path("msign_gap_3x2.csv"))
  od <- ref[ref$case == "gap" & ref$name == "order", ]
  run <- function(sl) eiBounds(x = as.matrix(gp[, 1:3]), t = as.matrix(gp[, 5:6]), n = gp$N, targets = list(B2 = c(2, 1)),
                               ci = "none", tests = FALSE, diagnostics = FALSE,
                               restrictions = list(order = list(list(high = od$r, low = od$s, col = 1)),
                                                   cr = list(list(row = 1, version = "W2", slack = sl))))
  u0 <- run(0)$pieces$B2$estimate; u1 <- run(0.005)$pieces$B2$estimate
  inside <- function(a, b) all(apply(a, 1, function(r) any(b[, 1] <= r[1] + 1e-9 & r[2] <= b[, 2] + 1e-9)))
  expect_true(inside(u0, u1))
})
