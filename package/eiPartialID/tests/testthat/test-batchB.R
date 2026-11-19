# Package changes made in response to the two batch-B referee reports (M-SIGN; release 0.3.0).
# Plan: docs/extension/reviews/B_PI_RESPONSE_PLAN.md (WS-B0, WS-B3 code part, WS-B9 subset, WS6 items).

sim2 <- function(n = 200, seed = 7, ab = c(0.40, 0.20, 0.55, -0.05), conc = 15) {
  set.seed(seed)
  X <- pmin(pmax(stats::rbeta(n, 2, 2) * 0.8 + 0.1, 0.02), 0.98); N <- sample(300:3000, n, replace = TRUE)
  m0 <- ab[1] + ab[2] * X; m1 <- ab[3] + ab[4] * X
  t <- (1 - X) * stats::rbeta(n, m0 * conc, (1 - m0) * conc) + X * stats::rbeta(n, m1 * conc, (1 - m1) * conc)
  list(x = X, t = t, n = N)
}
simg <- function(G = 6, ng = 60, zeta_sd = 0, seed = 1, conc = 20) {
  set.seed(seed)
  z <- if (zeta_sd > 0) matrix(stats::runif(2 * G, -zeta_sd, zeta_sd), G) else matrix(0, G, 2)
  z[1, ] <- 0
  g <- rep(seq_len(G), each = ng)
  a <- stats::runif(G, 0.1, 0.5)
  x1 <- unlist(lapply(seq_len(G), function(k) stats::runif(ng, a[k], a[k] + 0.4)))
  m1 <- 0.70 * x1 + 0.45 * (1 - x1) + z[g, 1]
  m2 <- 0.30 * x1 + 0.20 * (1 - x1) + z[g, 2]
  b1 <- stats::rbeta(length(x1), m1 * conc, (1 - m1) * conc)
  b2 <- stats::rbeta(length(x1), m2 * conc, (1 - m2) * conc)
  t1 <- x1 * b1 + (1 - x1) * b2
  list(X = cbind(x1, 1 - x1), T = cbind(t1, 1 - t1), N = sample(200:1000, length(x1), replace = TRUE), g = g)
}
Dg <- list(D = matrix(c(1, -1, 0, 0), 2, 2))
cr1 <- list(cr = list(list(row = 1, col = 1)))

# ---------------------------------------------------------------------------------------------
# WS-B0: labels and wording (no numeric change)
# ---------------------------------------------------------------------------------------------
test_that("M-SIGN is labelled experimental and print() shows the EXPERIMENTAL flag", {
  m <- eiModules()
  expect_match(m$proof_status[m$type == "assoc"], "^experimental")
  expect_match(m$proof_status[m$type == "cr"], "^experimental")
  s <- sim2()
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE, restrictions = cr1)
  out <- capture.output(print(f))
  expect_true(any(grepl("EXPERIMENTAL module CR", out)))
  expect_true(any(grepl("sharp for the regression \\(PHASE5a R3\\(b\\), pooled model\\)", out)))
})

test_that("no unqualified 'never refuted' claim in notes, module meta or help pages", {
  bad <- "never be refuted|cannot be refuted by aggregates on their own|never alone with 2 rows"
  expect_false(any(grepl(bad, unlist(eiModules()))))
  rd <- tools::Rd_db("eiPartialID")
  txt <- paste(unlist(lapply(rd, function(x) paste(as.character(x), collapse = ""))), collapse = " ")
  expect_false(grepl(bad, txt))
  expect_true(grepl("continuous unit tables", txt))
  v <- utils::read.csv(test_path("voter_vax_MonotoneEI.csv"))
  f <- eiBounds(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes, weights = "n", targets = Dg,
                ci = "conservative", seed = 1,
                restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1),
                                                 list(type = "dB", row = 1, col = 1, sign = 1))))
  expect_false(any(grepl(bad, f$notes)))
  expect_true(any(grepl("sample diagnostic, not a test", f$notes)))
  expect_true(any(grepl("rung alone is not rejected", f$notes)))       # beside the joint rejection
  expect_true(any(grepl("rung alone \\(M0\\) has a non-empty estimate", f$notes)))
  expect_false(any(grepl("refuted", f$sign$verdicts$verdict_first)))
  expect_false(any(grepl("refuted", f$sign$verdicts$verdict_last)))
  expect_true(all(f$sign$verdicts$verdict_last == "empty at the point estimate (see the emptiness test)"))
})

test_that("face test includes the target intercept (P-B18)", {
  sq <- list(s = c(0, 1, 0, 1), A = rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1)))
  tg <- list(cvec = numeric(0), c0 = 2, gamma = c(1, 0))          # phi = 2 + z1, U = 3
  keep <- list(list(s = 1.2, A = rbind(c(-1, -1))))                  # meets the face z1 = 1
  cut <- list(list(s = 0.8, A = rbind(c(-1, -1))))                   # removes it
  expect_false(eiPartialID:::.rxc_face_moves(sq, tg, numeric(0), keep, "upper", 3, 1e-7))
  expect_true(eiPartialID:::.rxc_face_moves(sq, tg, numeric(0), cut, "upper", 3, 1e-7))
  # target-specific zero gain (R2(e), P-B21): the row is not implied, yet does not move z1's end
  Fo <- list(s = c(sq$s, 0.5), A = rbind(sq$A, c(-1, 0)))
  N <- list(list(s = 0.8, A = rbind(c(-1, -1))))
  expect_false(eiPartialID:::.rxc_implied(Fo, N, 1e-7))
  expect_false(eiPartialID:::.rxc_face_moves(Fo, list(cvec = numeric(0), c0 = 0, gamma = c(1, 0)),
                                             numeric(0), N, "upper", 0.5, 1e-7))
  expect_true(eiPartialID:::.rxc_face_moves(Fo, list(cvec = numeric(0), c0 = 0, gamma = c(1, 1)),
                                            numeric(0), N, "upper", 1.5, 1e-7))
})

test_that("grouped rungs: fit-time note (one pooled restriction) and no R3(b) sharpness claim (U2(a))", {
  d <- simg(seed = 4)
  f <- eiBounds(x = d$X, t = d$T, n = d$N, group = d$g, model = "M1", targets = Dg, ci = "none",
                tests = FALSE, restrictions = cr1)
  expect_true(any(grepl("one pooled restriction across groups", f$notes)))
  expect_true(any(grepl("experimental", f$notes)))
  out <- capture.output(print(f))
  expect_false(any(grepl("R3\\(b\\), pooled model", out)))
  expect_true(any(grepl("sharpness under grouped rungs is not yet proved", out)))
  expect_true(any(grepl("EXPERIMENTAL", out)))
  # the M1 rows are the district functionals: each unit uses its own group's H (.rxc_mlin), so the
  # W1 identity delta_W = Var(X)(ER - D) holds with the district-fitted tau at every k
  st <- f$setup
  dz <- eiPartialID:::.rxc_design_r(st, 1)
  fw <- eiPartialID:::.rxc_assoc_form(st, "dW", "W1", 1, 1)
  fD <- eiPartialID:::.rxc_assoc_form(st, "D", "W1", 1, 1)
  fE <- eiPartialID:::.rxc_ER_form(st, 1, 1)
  set.seed(1)
  for (k in list(rep(0, st$nz), stats::rnorm(st$nz))) {
    lhs <- eiPartialID:::.rxc_val(fw, st$beta) + sum(fw$a * k)
    rhs <- dz$var * (eiPartialID:::.rxc_val(fE, st$beta) - eiPartialID:::.rxc_val(fD, st$beta) - sum(fD$a * k))
    expect_equal(lhs, rhs, tolerance = 1e-9)
  }
})

test_that("W2 is documented and reported as one pooled sum (P-E1, U8(a))", {
  ex <- utils::read.csv(test_path("msign_gap_3x2.csv"))
  f <- eiBounds(x = as.matrix(ex[, c("X1", "X2", "X3")]), t = as.matrix(ex[, c("T1", "T2")]), n = ex$N,
                targets = list(B1 = c(1, 1)), ci = "none", tests = FALSE,
                restrictions = list(cr = list(list(row = 1, col = 1, version = "W2"))))
  expect_true(any(grepl("W2 signs one pooled sum", f$notes)))
  expect_true(any(grepl("W2 is one pooled sum", capture.output(print(f)))))
  expect_match(eiModules()$object[eiModules()$type == "cr"], "W2 pooled over groups")
})

# ---------------------------------------------------------------------------------------------
# WS-B3 (code): projection rank guard and zero-variance candidates
# ---------------------------------------------------------------------------------------------
test_that("projection rank guard: candidates using a null direction of V-hat are dropped (P-B40)", {
  V <- diag(c(1, 0))
  st <- list(beta = c(0, 0))
  Zs <- list(V = V)
  cu <- list(G = rbind(c(1, 0), c(0, 1)), ic = c(1, 0.5))   # second candidate: sigma = 0, nonzero gradient
  cl <- list(G = rbind(c(1, 0)), ic = -1)
  k <- sqrt(stats::qchisq(0.95, 2))                          # q = dim(beta) = 2 (0.3.1), not rank 1
  pr <- eiPartialID:::.rxc_projection(st, cu, cl, Zs, 0.05)
  expect_equal(as.numeric(pr), c(-1 - k, 1 + k))             # old code: upper end 0.5 (sigma 0 taken as exact)
  expect_equal(attr(pr, "dropped"), 1)
  pr2 <- eiPartialID:::.rxc_projection(st, list(G = rbind(c(0, 1)), ic = 0.5), cl, Zs, 0.05)
  expect_equal(pr2[2], Inf)                                    # nothing left: no information on that side
  # nonsingular V: unchanged; V identically zero (exact fit): every candidate is exact
  expect_equal(attr(eiPartialID:::.rxc_projection(st, cu, cl, list(V = diag(2)), 0.05), "dropped"), 0)
  expect_false(any(eiPartialID:::.rxc_null_use(cu$G, matrix(0, 2, 2))))
  # a nonzero gradient with a small component in the null space is caught; a zero gradient is not
  expect_equal(eiPartialID:::.rxc_null_use(rbind(c(1, 1e-3), c(0, 0), c(1, 0)), V), c(TRUE, FALSE, FALSE))
})

test_that("CLR: a zero-variance candidate is exact only if its gradient is zero (P-B23)", {
  set.seed(1)
  V <- diag(c(1, 0))
  Zs <- list(V = V, draws = matrix(stats::rnorm(4000), 2000) %*% sqrt(V))
  # candidate 2 has sigma = 0 but a nonzero gradient: dropped, so the bound is not 0.5
  u <- eiPartialID:::.rxc_clr_min(c(1, 0.5), rbind(c(1, 0), c(0, 1)), Zs, 0.05, 200)
  expect_equal(u, 1 + stats::qnorm(0.95))
  # a candidate with a zero gradient is a known constant: kept as exact
  u0 <- eiPartialID:::.rxc_clr_min(c(1, 0.5), rbind(c(1, 0), c(0, 0)), Zs, 0.05, 200)
  expect_lte(u0, 0.5 + 1e-12)
  expect_equal(eiPartialID:::.rxc_clr_min(0.5, rbind(c(0, 1)), Zs, 0.05, 200), Inf)
})

test_that("a fit on a full-rank design never fires the rank guard", {
  s <- sim2(seed = 3)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, seed = 1, tests = FALSE)
  expect_false(any(grepl("numerically singular", f$notes)))
})

# ---------------------------------------------------------------------------------------------
# WS6 items: pooling test calibration; proj_actual covered by P3.5(b2)
# ---------------------------------------------------------------------------------------------
test_that("pooling test: wild bootstrap p-value by default, chi^2 and F reported, RNG stream untouched", {
  d <- simg(zeta_sd = 0.08, seed = 7)
  set.seed(11); r1 <- stats::runif(1)
  set.seed(11); pt <- eiPartialID:::.rxc_pool_test(d$X, d$T, d$N, d$g, B = 99); r2 <- stats::runif(1)
  expect_equal(r1, r2)
  expect_equal(pt$method, "wild")
  expect_equal(pt$df, (6 - 1) * (2 - 1) * 2)
  expect_true(all(c("p.value_chisq", "p.value_F") %in% names(pt)))
  expect_equal(pt$p.value, 1 / 100)                            # shifts present: never exceeded
  pc <- eiPartialID:::.rxc_pool_test(d$X, d$T, d$N, d$g, method = "chisq")
  expect_equal(pc$p.value, pc$p.value_chisq)
  expect_equal(pc$statistic, pt$statistic)
  # the bootstrap statistic at the observed data equals the reported Wald statistic
  f <- eiBounds(x = d$X, t = d$T, n = d$N, group = d$g, model = "M1", seed = 1, ci = "none", targets = Dg)
  expect_output(print(f), "wild bootstrap, B = 399")
})

test_that("pooling test: size of the wild bootstrap in a small null design", {
  rej <- vapply(1:30, function(s) {
    d0 <- simg(zeta_sd = 0, seed = 300 + s)
    eiPartialID:::.rxc_pool_test(d0$X, d0$T, d0$N, d0$g, B = 99)$p.value <= 0.05
  }, logical(1))
  expect_lte(mean(rej), 0.2)
})

test_that("proj_actual takes its Bhatia-Davis bound from a projection interval (WS6(b))", {
  s <- sim2(seed = 3)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE)
  st <- f$setup
  cons <- eiPartialID:::.rxc_constraints(st)
  set.seed(5); Zs <- eiPartialID:::.rxc_zdraws(st$V, 2000)
  set.seed(6); sd2 <- eiPartialID:::.rxc_fpc_sd(st, cons, 1, 1, "bhatia-davis", 0.005, Zs, 50)
  expect_length(sd2, 2)
  # recompute both from the same draws
  wr <- st$w[, 1]; wpi <- st$w; wpi[, 1] <- wr^2 / sum(wr^2)
  tg <- eiPartialID:::.rxc_target(st, matrix(c(1, 0, 0, 0), 2), w = wpi)
  set.seed(6); ci <- eiPartialID:::.rxc_ci_clr(st, cons, tg, 0.0025, 0.0025, Zs, 50, proj_alpha = 0.005)
  bd <- function(iv) { lo <- max(0, min(1, iv[1])); hi <- max(0, min(1, iv[2]))
    sqrt(sum(wr^2) * (if (lo <= 0.5 && hi >= 0.5) 0.25 else max(lo * (1 - lo), hi * (1 - hi)))) }
  expect_equal(as.numeric(sd2), c(bd(ci[1:2]), bd(attr(ci, "projection"))))
  expect_equal(eiPartialID:::.rxc_fpc_sd(st, cons, 1, 1, "popoviciu", 0.005, Zs, 50), rep(0.5 * sqrt(sum(wr^2)), 2))
  tb <- f$table
  expect_true(tb$proj_actual_lower <= tb$proj_model_lower && tb$proj_actual_upper >= tb$proj_model_upper)
})

# ---------------------------------------------------------------------------------------------
# U6(a), U7(a): CR interval labels; every union piece; base column; weight concentration
# ---------------------------------------------------------------------------------------------
test_that("CR fits print 'algorithmic' on the CLR lines and the projection beside them (U6(a))", {
  s <- sim2(seed = 11)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, seed = 2, tests = FALSE, restrictions = cr1)
  out <- capture.output(print(f))
  expect_true(any(grepl("proj_model", out)))
  expect_true(any(grepl("ci_model: .*", out)) && any(grepl("CLR: algorithmic", out)))
  expect_true(any(grepl("CLR \\(algorithmic\\); its Bhatia-Davis bound also uses a CLR interval", out)))
  expect_true(any(grepl("union of the CLR branch intervals \\(algorithmic\\)", out)))
  expect_false(is.na(f$table$proj_model_lower))
})

test_that("print shows every piece of a union with three or more pieces (P-E15)", {
  s <- sim2(seed = 11)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE, restrictions = cr1,
                diagnostics = FALSE)
  # mock a three-piece estimate (as several W2 / pair CRs can produce)
  u <- rbind(c(-0.5, -0.4), c(-0.3, -0.2), c(-0.1, 0))
  f$pieces$D$estimate <- u
  f$table$estimate_lower <- -0.5; f$table$estimate_upper <- 0
  f$table$estimate_gap_lower <- -0.4; f$table$estimate_gap_upper <- -0.3; f$table$n_pieces <- 3L
  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "[-0.500, -0.400] U [-0.300, -0.200] U [-0.100, 0.000]", fixed = TRUE)
  f$settings$union <- "hull"
  f$table$estimate_gap_lower <- NA; f$table$estimate_gap_upper <- NA
  expect_match(paste(capture.output(print(f)), collapse = "\n"), "[-0.500, 0.000]", fixed = TRUE)
})

test_that("base column: the rung-only result is reported beside the restricted one (U7(a))", {
  s <- sim2(seed = 11)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, seed = 2, restrictions = cr1)
  g <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, seed = 2)
  tb <- f$table
  expect_true(all(c("base_estimate_lower", "base_ci_model_lower", "base_ci_actual_upper", "base_empty") %in% names(tb)))
  expect_equal(tb$base_estimate_lower, g$table$estimate_lower)
  expect_equal(tb$base_ci_model_lower, g$table$ci_model_lower)        # same draws as a separate fit
  expect_equal(tb$base_ci_actual_upper, g$table$ci_actual_upper)
  expect_output(print(f), "The rung alone \\(M0, without the optional restriction\\)")
  expect_false("base_estimate_lower" %in% names(g$table))
  expect_false("base_estimate_lower" %in% names(eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none",
                                                          tests = FALSE, restrictions = cr1, base = FALSE)$table))
})

test_that("weight concentration per signed association (P-E25, U7(a))", {
  s <- sim2(seed = 11)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dg, ci = "none", tests = FALSE, restrictions = cr1)
  wc <- f$sign$concentration
  expect_equal(wc$association, c("delta_B(X vs rest)", "delta_W(X vs rest)"))
  # 2 x 2 delta_B: unit weights proportional to N_i x_i (1 - x_i)
  u <- s$n * s$x * (1 - s$x); u <- u / sum(u)
  expect_equal(wc$max_weight[1], max(u)); expect_equal(wc$unit_max[1], which.max(u))
  expect_equal(wc$n_eff[1], 1 / sum(u^2))
  expect_output(print(f$sign), "Weight concentration")
})

# ---------------------------------------------------------------------------------------------
# Decisions 3(a) and 3(b) (Nov 2026): conservative ci_actual and the M1-vs-M2 test
# ---------------------------------------------------------------------------------------------
test_that("conservative ci_actual takes its Bhatia-Davis bound from the conservative interval (3(a))", {
  s <- sim2(seed = 3)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE)
  st <- f$setup
  cons <- eiPartialID:::.rxc_constraints(st)
  set.seed(5); Zs <- eiPartialID:::.rxc_zdraws(st$V, 2000)
  sdc <- eiPartialID:::.rxc_fpc_sd(st, cons, 1, 1, "bhatia-davis", 0.005, Zs, 50, ci = "conservative")
  wr <- st$w[, 1]; wpi <- st$w; wpi[, 1] <- wr^2 / sum(wr^2)
  tg <- eiPartialID:::.rxc_target(st, matrix(c(1, 0, 0, 0), 2), w = wpi)
  cc <- eiPartialID:::.rxc_ci_conservative(st, cons, tg, 0.005, Zs)$ci
  lo <- max(0, min(1, cc[1])); hi <- max(0, min(1, cc[2]))
  expect_equal(as.numeric(sdc), rep(sqrt(sum(wr^2) * (if (lo <= 0.5 && hi >= 0.5) 0.25 else max(lo * (1 - lo), hi * (1 - hi)))), 2))
  fc <- eiBounds(x = s$x, t = s$t, n = s$n, targets = list(a = c(1, 1)), seed = 1, tests = FALSE, ci = "conservative")
  expect_true(fc$table$ci_actual_lower <= fc$table$ci_model_lower + 1e-9)
  expect_true(fc$table$ci_actual_lower >= fc$table$dd_lower - 1e-9 && fc$table$ci_actual_upper <= fc$table$dd_upper + 1e-9)
})

test_that("M1-vs-M2 test: wild bootstrap by default, chi^2 and F reported, RNG untouched (3(b))", {
  d <- simg(zeta_sd = 0.08, seed = 7)
  set.seed(11); r1 <- stats::runif(1)
  set.seed(11); mt <- eiPartialID:::.rxc_m1_test(d$X, d$T, d$N, d$g, B = 99); r2 <- stats::runif(1)
  expect_equal(r1, r2)
  expect_equal(mt$method, "wild")
  expect_equal(mt$df, (6 - 1) * (2 - 1) * 2 * 1 / 2)
  expect_true(all(c("p.value_chisq", "p.value_F") %in% names(mt)))
  mc <- eiPartialID:::.rxc_m1_test(d$X, d$T, d$N, d$g, method = "chisq")
  expect_equal(mc$statistic, mt$statistic)
  expect_equal(mc$p.value, mc$p.value_chisq)
  rej <- vapply(1:30, function(s) {
    d0 <- simg(zeta_sd = 0.08, seed = 400 + s)                  # M1 holds (shifts, common slopes)
    eiPartialID:::.rxc_m1_test(d0$X, d0$T, d0$N, d0$g, B = 99)$p.value <= 0.05
  }, logical(1))
  expect_lte(mean(rej), 0.2)
})
