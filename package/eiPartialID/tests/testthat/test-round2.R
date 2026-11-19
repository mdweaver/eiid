# Round-2 referee reports (Nov 2026): EI report items 1, 4 and "also"; Phases 1-3 report item 3.
# Plan: docs/extension/reviews/ROUND2_PACKAGE_PLAN.md (decisions D1 = A, D2 = projection column only).

sim2r <- function(n = 200, seed = 7, ab = c(0.40, 0.20, 0.55, -0.05), conc = 15) {
  set.seed(seed)
  X <- pmin(pmax(stats::rbeta(n, 2, 2) * 0.8 + 0.1, 0.02), 0.98); N <- sample(300:3000, n, replace = TRUE)
  m0 <- ab[1] + ab[2] * X; m1 <- ab[3] + ab[4] * X
  t <- (1 - X) * stats::rbeta(n, m0 * conc, (1 - m0) * conc) + X * stats::rbeta(n, m1 * conc, (1 - m1) * conc)
  list(x = X, t = t, n = N)
}
simg2 <- function(G = 6, ng = 60, seed = 4, conc = 20) {
  set.seed(seed)
  g <- rep(seq_len(G), each = ng)
  a <- stats::runif(G, 0.1, 0.5)
  x1 <- unlist(lapply(seq_len(G), function(k) stats::runif(ng, a[k], a[k] + 0.4)))
  m1 <- 0.70 * x1 + 0.45 * (1 - x1); m2 <- 0.30 * x1 + 0.20 * (1 - x1)
  t1 <- x1 * stats::rbeta(length(x1), m1 * conc, (1 - m1) * conc) + (1 - x1) * stats::rbeta(length(x1), m2 * conc, (1 - m2) * conc)
  list(X = cbind(x1, 1 - x1), T = cbind(t1, 1 - t1), N = sample(200:1000, length(x1), replace = TRUE), g = g)
}
Dg2 <- list(D = matrix(c(1, -1, 0, 0), 2, 2))
cr1 <- list(cr = list(list(row = 1, col = 1)))

test_that("EI-1: a positive W1 sign does not mean both groups' rates rise (two-precinct example)", {
  # A's share .25 / .75; A's rate .40 -> .60, the rest's .55 -> .45 (equal precincts)
  st <- list(X = cbind(c(.25, .75), c(.75, .25)), N = c(1, 1), R = 2)
  W <- eiPartialID:::.rxc_assoc_W(st, "dW", "W1", 1)
  partA <- sum(W[, 1] * c(.40, .60)); partR <- sum(W[, 2] * c(.55, .45))
  expect_equal(partA + partR, 0.0046875)
  expect_equal(partA, 0.009375)
  expect_equal(partR, -0.0046875)
})

test_that("EI-1: the pairwise CR label names both sides of the pairing", {
  ex <- utils::read.csv(test_path("msign_gap_3x2.csv"))
  f <- eiBounds(x = as.matrix(ex[, c("X1", "X2", "X3")]), t = as.matrix(ex[, c("T1", "T2")]), n = ex$N,
                targets = list(B1 = c(1, 1)), ci = "none", tests = FALSE, diagnostics = FALSE, base = FALSE,
                restrictions = list(cr = list(list(row = 1, col = 1, version = "pair", other = 2))))
  expect_match(f$modules$restriction, "\\[pair\\] and Delta_W\\(X1\\) \\[W2, pooled over all groups\\]")
})

test_that("EI-4: print states which result covers each fit's intervals", {
  s <- sim2r(seed = 3)
  f0 <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg2, seed = 1, tests = FALSE)
  o0 <- capture.output(print(f0))
  expect_true(any(grepl("asymptotic guarantee \\(P3.5\\(b2\\)\\(ii\\)\\), pooled model; V-hat rank 3 of 3 \\(diagnostic\\)", o0)))
  fc <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg2, seed = 1, tests = FALSE, restrictions = cr1)
  expect_true(any(grepl("with sign/CR rows PHASE5a R10, pooled M0", capture.output(print(fc)))))
  fk <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg2, seed = 1, tests = FALSE, ci = "conservative")
  expect_true(any(grepl("ci_\\* \\(conservative\\): asymptotic guarantee \\(P3.5\\(a\\)\\)", capture.output(print(fk)))))
  d <- simg2()
  g1 <- eiBounds(x = d$X, t = d$T, n = d$N, group = d$g, model = "M1", targets = Dg2, seed = 1, tests = FALSE)
  expect_true(any(grepl("only with few, large groups \\(P3 regime R1\\)", capture.output(print(g1)))))
  gc <- eiBounds(x = d$X, t = d$T, n = d$N, group = d$g, model = "M1", targets = Dg2, seed = 1, tests = FALSE,
                 restrictions = cr1)
  expect_true(any(grepl("no coverage theorem for sign/CR restrictions under grouped rungs \\(R16 open\\)",
                        capture.output(print(gc)))))
  expect_true(any(grepl("nor are its confidence intervals", gc$notes)))
  expect_equal(gc$inference$p, ncol(gc$setup$V))
})

test_that("projection critical value uses q = dim(beta); the rank is a diagnostic only", {
  st <- list(beta = c(0, 0, 0))
  cu <- list(G = rbind(c(1, 0, 0)), ic = 1); cl <- list(G = rbind(c(1, 0, 0)), ic = -1)
  pr <- eiPartialID:::.rxc_projection(st, cu, cl, list(V = diag(3)), 0.05)
  expect_equal(as.numeric(pr), c(-1, 1) + c(-1, 1) * sqrt(stats::qchisq(0.95, 3)))
  expect_equal(eiPartialID:::.rxc_vrank(diag(c(1, 1, 0))), 2L)
})

test_that("EI-also: with several restrictions the base block is labelled as their joint effect", {
  s <- sim2r(seed = 11)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg2, ci = "none", tests = FALSE, diagnostics = FALSE,
                restrictions = c(cr1, list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1)))))
  out <- capture.output(print(f))
  expect_true(any(grepl("all 2 optional restrictions removed: their joint effect", out)))
  expect_true(any(grepl("added-last gains in eiRedundancy", out)))
  f1 <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg2, ci = "none", tests = FALSE, diagnostics = FALSE,
                 restrictions = cr1)
  expect_true(any(grepl("without the optional restriction\\)", capture.output(print(f1)))))
})

test_that("II: proj_actual takes its Bhatia-Davis bound from a rung-only projection interval; ci_actual unchanged", {
  s <- sim2r(seed = 3)
  res <- list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1)))
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE,
                restrictions = res, diagnostics = FALSE, base = FALSE)
  st <- f$setup
  cons <- eiPartialID:::.rxc_constraints(st)
  st0 <- st; st0$extra <- NULL; cons0 <- eiPartialID:::.rxc_constraints(st0)
  set.seed(5); Zs <- eiPartialID:::.rxc_zdraws(st$V, 2000)
  set.seed(6); new <- eiPartialID:::.rxc_fpc_sd(st, cons, 1, 1, "bhatia-davis", 0.005, Zs, 50, cons_rung = cons0)
  r1 <- stats::runif(1)
  set.seed(6); old <- eiPartialID:::.rxc_fpc_sd(st, cons, 1, 1, "bhatia-davis", 0.005, Zs, 50)
  r2 <- stats::runif(1)
  expect_equal(r1, r2)                                   # the extra rung-only call leaves the stream intact
  expect_equal(new[1], old[1])                           # CLR bound for ci_actual unchanged (D2)
  set.seed(6); rung <- eiPartialID:::.rxc_fpc_sd(st0, cons0, 1, 1, "bhatia-davis", 0.005, Zs, 50)
  expect_equal(new[2], rung[2])                          # projection bound = the rung alone's
  out <- capture.output(print(f))
  expect_true(any(grepl("projection interval for the rung alone", out)))
  expect_true(any(grepl("CLR \\(algorithmic\\); its Bhatia-Davis bound also uses a CLR interval", out)))
  # the Bhatia-Davis sd bound never exceeds Popoviciu's (the columns themselves need not nest:
  # with Popoviciu the variance-bound share of alpha goes to estimation)
  expect_lte(new[2], 0.5 * sqrt(sum(st$w[, 1]^2)) + 1e-12)
  fp <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE,
                 restrictions = res, diagnostics = FALSE, base = FALSE, fpc = "popoviciu")
  expect_true(any(grepl("Popoviciu bound 1/4", capture.output(print(fp)))))
})

test_that("II: near-boundary note when the weighted mean rate's interval lies within [0, .05] (P3.5(c1))", {
  s <- sim2r(seed = 2, ab = c(0.01, 0, 0.01, 0), conc = 300)    # every rate about .01
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE)
  expect_true(f$inference$bd_boundary)
  expect_true(any(grepl("bounded away from 0 and 1 \\(P3.5\\(c1\\)\\)", f$notes)))
  fp <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE, fpc = "popoviciu")
  expect_false(any(grepl("P3.5\\(c1\\)", fp$notes)))
  s2 <- sim2r(seed = 3)
  g <- eiBounds(x = s2$x, t = s2$t, n = s2$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE)
  expect_false(g$inference$bd_boundary)
})
