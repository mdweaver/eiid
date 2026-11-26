# eiChange / eiATT: brute force over integer two-election tables for every descriptive rung (2e
# Props E3, E5, E6, E26, Thm E27), closed forms against the LPs (E5, E15, E17, E20, E21), the
# labelling rule (Def. 5.1), unions (E5(f), E18, E24), Prop. E25, and the gap register (Sec. 8).

two_el <- function(a, b0, b1, units = paste0("U", seq_len(nrow(a)))) {
  mk <- function(b) data.frame(u = units, g1 = a[, 1], g2 = a[, 2], yes = b, no = rowSums(a) - b)
  list(e0 = mk(b0), e1 = mk(b1))
}
# every integer two-election 2 x 2 table: y[t, i] = n_{i, g1, yes, t}
brute <- function(a, b) {
  rg <- list()
  for (t in 1:2) for (i in seq_len(nrow(a))) {
    N <- sum(a[i, ]); rg[[length(rg) + 1]] <- seq(max(0, a[i, 1] + b[[t]][i] - N), min(a[i, 1], b[[t]][i]))
  }
  g <- as.matrix(expand.grid(rg)); p <- nrow(a)
  list(y0 = g[, 1:p, drop = FALSE], y1 = g[, p + 1:p, drop = FALSE])
}
pieces <- function(a, b, bt) {
  N <- rowSums(a); pw <- N / sum(N); x <- a[, 1] / N; A1 <- sum(a[, 1]); A2 <- sum(a[, 2])
  dYt <- sum(b[[2]] - b[[1]])
  dY1 <- rowSums(bt$y1) - rowSums(bt$y0)
  D <- dY1 / A1 - (dYt - dY1) / A2
  dT <- (b[[2]] - b[[1]]) / N
  mu <- sum(pw * x); vx <- sum(pw * x^2) - mu^2; cv <- sum(pw * x * dT) - mu * sum(pw * dT)
  ER <- cv / vx; NM <- cv / (mu * (1 - mu))
  list(D = D, dW = vx * (ER - D), dB = mu * (1 - mu) * (D - NM))
}
designs <- function(n = 12, seed = 7) {
  set.seed(seed)
  lapply(seq_len(n), function(k) {
    p <- 3; N <- sample(4:7, p, TRUE); a1 <- vapply(N, function(m) sample(1:(m - 1), 1), 0)
    a <- cbind(a1, N - a1); b0 <- vapply(N, function(m) sample(0:m, 1), 0); b1 <- vapply(N, function(m) sample(0:m, 1), 0)
    list(a = a, b = list(b0, b1))
  })
}
fitd <- function(d, rows = NULL, integer = "auto", targets = NULL)
  eiChange(two_el(d$a, d$b[[1]], d$b[[2]]), cbind(yes, no) ~ cbind(g1, g2), unit = "u", row = "g1", col = "yes",
           integer = integer, restrictions = rows, targets = targets)
rngs <- function(v) if (!length(v)) c(NA, NA) else range(v)

test_that("brute force: margins, W1 signs, CR on changes (integer sets of Props E3, E5, E26)", {
  n_cr_empty <- 0
  for (d in designs()) {
    bt <- brute(d$a, d$b); pc <- pieces(d$a, d$b, bt)
    tol <- 1e-9
    f <- fitd(d)
    expect_equal(c(f$table$lower, f$table$upper), range(pc$D), tolerance = 1e-9)
    for (rw in list(list(type = "dW", sign = 1, keep = pc$dW >= -tol), list(type = "dW", sign = -1, keep = pc$dW <= tol),
                    list(type = "dB", sign = 1, keep = pc$dB >= -tol), list(type = "dB", sign = -1, keep = pc$dB <= tol))) {
      f <- fitd(d, list(sign_change = list(type = rw$type, sign = rw$sign, row = "g1", col = "yes", from = "e0", to = "e1")))
      v <- pc$D[rw$keep]
      if (!length(v)) expect_true(is.na(f$table$lower)) else expect_equal(c(f$table$lower, f$table$upper), range(v), tolerance = 1e-9)
      expect_match(f$table$integer_status[1], "E26|infeasible")
      # the general TU path (Thm E27, R = 2) agrees with the E26 lattice
      g <- eiPartialID:::.ec_range(f$P, f$built, f$targets[[1]], "auto")
      if (length(v)) expect_equal(c(g$lower, g$upper), range(v), tolerance = 1e-9)
    }
    keep <- (pc$dW >= -tol & pc$dB >= -tol) | (pc$dW <= tol & pc$dB <= tol)
    f <- fitd(d, list(cr_change = list(row = "g1", col = "yes", from = "e0", to = "e1")))
    v <- pc$D[keep]
    if (!length(v)) { expect_true(is.na(f$table$lower)); n_cr_empty <- n_cr_empty + 1 } else
      expect_equal(c(f$table$lower, f$table$upper), range(v), tolerance = 1e-9)
    # continuous CR is never empty (E5(d))
    expect_true(f$table$continuous_feasible)
  }
})

test_that("brute force: unit orderings of changes and the count contrast on changes (MILP)", {
  for (d in designs(6, seed = 11)) {
    bt <- brute(d$a, d$b); pc <- pieces(d$a, d$b, bt); N <- rowSums(d$a)
    # order_change vs rest in every unit: n_t - n_u >= a Db / N
    ok <- apply(bt$y1 - bt$y0, 1, function(z) all(z >= d$a[, 1] * (d$b[[2]] - d$b[[1]]) / N - 1e-9))
    f <- fitd(d, list(order_change = list(row = "g1", col = "yes", from = "e0", to = "e1")))
    v <- pc$D[ok]
    if (!length(v)) expect_false(isTRUE(f$table$integer_feasible)) else expect_equal(c(f$table$lower, f$table$upper), range(v), tolerance = 1e-9)
    # contrast on changes: bins {U1} and {U2, U3} ordered by the share of g1
    x <- d$a[, 1] / N
    sb <- function(S) sum(d$a[S, 1]) / sum(N[S])
    L <- if (sb(1) < sb(2:3)) 1 else 2:3; H <- setdiff(1:3, L)
    if (abs(sb(L) - sb(H)) < 1e-9) next
    bdid <- function(b) {
      d1 <- rowSums(bt$y1[, b, drop = FALSE]) - rowSums(bt$y0[, b, drop = FALSE]); dt <- sum(d$b[[2]][b] - d$b[[1]][b])
      d1 / sum(d$a[b, 1]) - (dt - d1) / sum(d$a[b, 2])
    }
    if (any(d$a[L, 1] == 0) && sum(d$a[L, 1]) == 0) next
    Dbar <- 0.5 * (bdid(L) + bdid(H))
    lam <- sb(H) - sb(L)
    dtau <- sum(d$b[[2]][H] - d$b[[1]][H]) / sum(N[H]) - sum(d$b[[2]][L] - d$b[[1]][L]) / sum(N[L])
    keep <- Dbar >= -1e-12 & lam * Dbar <= dtau + 1e-12
    uL <- paste0("U", L); uH <- paste0("U", H)
    f <- fitd(d, list(contrast_change = list(col = "yes", from = "e0", to = "e1", low = uL, high = uH)),
              targets = list(bin_did("yes", "e0", "e1", low = uL, high = uH)))
    expect_equal(c(f$table$dd_lower, f$table$dd_upper), range(Dbar), tolerance = 1e-9)
    v <- Dbar[keep]
    if (!length(v)) expect_false(isTRUE(f$table$integer_feasible)) else expect_equal(c(f$table$lower, f$table$upper), range(v), tolerance = 1e-9)
  }
})

test_that("Theorem E27: change rows on one column, two consecutive pairs: TU class, LP = MILP", {
  set.seed(5)
  p <- 4; N <- sample(20:40, p, TRUE); a1 <- round(N * runif(p, .2, .8)); a <- cbind(a1, N - a1)
  b <- lapply(1:3, function(t) vapply(N, function(m) sample(0:m, 1), 0))
  mk <- function(bb) data.frame(u = paste0("U", 1:p), g1 = a[, 1], g2 = a[, 2], yes = bb, no = N - bb)
  dl <- list(e0 = mk(b[[1]]), e1 = mk(b[[2]]), e2 = mk(b[[3]]))
  rows <- list(order_change = list(list(row = "g1", col = "yes", from = "e0", to = "e1", units = c("U1", "U2"), sign = 1)),
               sign_change = list(type = "dW", sign = 1, row = "g1", col = "yes", from = "e1", to = "e2"))
  tg <- list(did("g1", "g2", "yes", "e0", "e2"), rate_at("g1", "yes", "e1"))
  f <- eiChange(dl, cbind(yes, no) ~ cbind(g1, g2), unit = "u", base = "e0", targets = tg, restrictions = rows)
  if (all(f$table$integer_feasible %in% TRUE)) {
    expect_true(all(f$table$tu_class == "E27"))
    fx <- eiChange(dl, cbind(yes, no) ~ cbind(g1, g2), unit = "u", base = "e0", targets = tg, restrictions = rows, integer = "exact")
    expect_true(all(grepl("MILP", fx$table$integer_status)))
    expect_equal(fx$table$lower, f$table$lower, tolerance = 1e-12); expect_equal(fx$table$upper, f$table$upper, tolerance = 1e-12)
    w <- witness(f)
    expect_true(all(vapply(w, function(z) z$margins_ok && z$rows_ok && z$integer, TRUE)))
  }
  # three pairs sharing an election (degree 3) are outside E27: MILP
  cls <- eiPartialID:::.ed_classify
  P <- f$P; P$T <- 4
  mkrow <- function(pr) list(kind = "change", cell = c(1, 1), pair = pr, units = 1:2)
  expect_true(is.na(cls(P, list(mkrow(c(1, 2)), mkrow(c(2, 3)), mkrow(c(2, 4))))$class))
  expect_equal(cls(P, list(mkrow(c(1, 2)), mkrow(c(2, 3)), mkrow(c(3, 4))))$class, "E27")
})

test_that("Prop. E26 for R = 3 (one-vs-rest): lattice = MILP; E5 closed form = LP", {
  set.seed(2)
  p <- 3; A <- matrix(sample(1:3, 9, TRUE), p); N <- rowSums(A)
  b0 <- vapply(N, function(m) sample(0:m, 1), 0); b1 <- vapply(N, function(m) sample(0:m, 1), 0)
  mk <- function(bb) data.frame(u = 1:p, r1 = A[, 1], r2 = A[, 2], r3 = A[, 3], yes = bb, no = N - bb)
  dl <- list(e0 = mk(b0), e1 = mk(b1))
  for (rows in list(NULL, list(sign_change = list(type = "dW", sign = 1, row = "r1", col = "yes", from = "e0", to = "e1")),
                    list(cr_change = list(row = "r1", col = "yes", from = "e0", to = "e1")))) {
    f <- eiChange(dl, cbind(yes, no) ~ cbind(r1, r2, r3), unit = "u", row = "r1", col = "yes", restrictions = rows)
    if (is.na(f$table$lower)) next
    g <- eiPartialID:::.ec_range(f$P, f$built, f$targets[[1]], "exact")
    expect_equal(c(g$lower, g$upper), c(f$table$lower, f$table$upper), tolerance = 1e-9)
    fl <- eiChange(dl, cbind(yes, no) ~ cbind(r1, r2, r3), unit = "u", row = "r1", col = "yes", restrictions = rows, integer = "lp")
    expect_equal(c(fl$table$closed_lower, fl$table$closed_upper), c(fl$table$lp_lower, fl$table$lp_upper), tolerance = 1e-8)
  }
})

test_that("Prop. E3(c): the joint event-study set is not the box; eiJoint", {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10)); yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  f <- eiChange(list(early = mk(1), base = mk(2), post = mk(3)), cbind(yes, no) ~ cbind(vet, oth), unit = "p",
                base = "base", row = "vet", col = "yes", integer = "lp")
  expect_false(f$event_study$corner$margins_only)
  tb <- f$table
  expect_false(eiJoint(f, c(tb$upper[1], tb$lower[2])))
  expect_true(eiJoint(f, c(tb$upper[1], tb$upper[2])))       # both at their max: the base gap at its minimum serves both? check E3
  # each DiD's range is the Minkowski sum of the election gap ranges (E3(b))
  G <- f$event_study$G
  expect_equal(tb$dd_upper[2], G$G_upper[3] - G$G_lower[2], tolerance = 1e-9)
  expect_equal(tb$dd_lower[1], G$G_lower[1] - G$G_upper[2], tolerance = 1e-9)
})

test_that("unions: CR in two pairs gives 2^2 branches; RM is a union of 2L LPs; the union is reported", {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10)); yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  cr <- function(a, b) list(row = "vet", col = "yes", from = a, to = b)
  f <- eiChange(list(early = mk(1), base = mk(2), post = mk(3)), cbind(yes, no) ~ cbind(vet, oth), unit = "p",
                base = "base", row = "vet", col = "yes", integer = "lp", restrictions = list(cr_change = list(cr("early", "base"), cr("base", "post"))))
  br <- f$solutions[[2]]$branches
  expect_equal(nrow(br), 4)
  ok <- !is.na(br$lower)
  expect_equal(f$table$lower[2], min(br$lower[ok])); expect_equal(f$table$upper[2], max(br$upper[ok]))
  a <- eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = rel_mag(1), spillover = 0, no_anticipation = TRUE)
  b <- a$results[[1]]$branches
  expect_equal(nrow(b), 2 * 4)                                  # 2L RM branches x 4 count branches
  ok <- !is.na(b$lower)
  expect_equal(a$results[[1]]$lower, min(b$lower[ok]), tolerance = 1e-12)
  expect_equal(a$results[[1]]$upper, max(b$upper[ok]), tolerance = 1e-12)
})

# random margins for causal closed-form checks
rnd <- function(seed, p = 4, R = 2, T = 2, C = 2) {
  set.seed(seed)
  N <- sample(30:80, p, TRUE)
  A <- t(vapply(N, function(m) as.vector(stats::rmultinom(1, m, rep(1, R))), numeric(R)))
  dl <- lapply(seq_len(T), function(t) {
    B <- t(vapply(N, function(m) as.vector(stats::rmultinom(1, m, rep(1, C))), numeric(C)))
    d <- data.frame(u = seq_len(p)); for (r in seq_len(R)) d[[paste0("g", r)]] <- A[, r]
    for (c in seq_len(C)) d[[paste0("c", c)]] <- B[, c]; d })
  names(dl) <- paste0("e", seq_len(T))
  fo <- stats::as.formula(sprintf("cbind(%s) ~ cbind(%s)", paste0("c", seq_len(C), collapse = ","), paste0("g", seq_len(R), collapse = ",")))
  eiChange(dl, fo, unit = "u", base = if (T > 2) "e2" else NULL, row = "g1", col = "c1", integer = "lp")
}
A1 <- function(f, tr, comp = "g2", L = 0, ...) eiATT(f, treated = "g1", comparison = comp, population = "fixed", trend = tr,
                                                   spillover = L, no_anticipation = TRUE, ...)$results[[1]]

test_that("causal closed forms equal the LPs: E15, E17 (margins only), E20, E21", {
  for (sd in 1:8) {
    f <- rnd(sd)
    e <- A1(f, no_trend())
    expect_equal(c(e$lower, e$upper), c(e$B_r[1] - 1, e$B_r[2]), tolerance = 1e-9)   # E15
    expect_true(e$lower <= 0 && e$upper >= 0)
    for (ML in list(c(0.05, 0), c(0.05, 0.03), c(0.6, 0))) {
      e <- A1(f, trend_dial(ML[1]), L = ML[2])
      d <- e$dial
      expect_true(e$lower >= d$outer[1] - 1e-9 && e$upper <= d$outer[2] + 1e-9)    # E17(a)
      expect_equal(d$criterion_upper, abs(e$upper - d$outer[2]) < 1e-7)            # E17(b): criterion <=> equality
      expect_equal(d$criterion_lower, abs(e$lower - d$outer[1]) < 1e-7)
    }
    e <- A1(f, mtr_nd_ps(), L = 0.03)
    if (!is.na(e$lower)) expect_equal(c(e$lower, e$upper), e$e21, tolerance = 1e-9)                  # E21(c), named row
    f3 <- rnd(100 + sd, R = 3)
    e3 <- A1(f3, mtr_nd_ps(), comp = "rest", L = 0.03)
    if (!is.na(e3$lower)) expect_equal(c(e3$lower, e3$upper), e3$e21, tolerance = 1e-9)              # E21(c), rest
    e3 <- A1(f3, mtr_nd_ps(), comp = "g2")
    if (!is.na(e3$lower)) expect_equal(c(e3$lower, e3$upper), e3$e21, tolerance = 1e-9)
    fT <- rnd(200 + sd, T = 3)
    eS <- A1(fT, smoothness(0.05), L = 0.02, rate_limits = FALSE)
    expect_equal(c(eS$lower, eS$upper), eS$e20, tolerance = 1e-9)                                    # E20
    eR <- A1(fT, rel_mag(1))
    eD <- A1(fT, trend_dial(eR$rm$M_implied))
    expect_true(eR$lower >= eD$lower - 1e-9 && eR$upper <= eD$upper + 1e-9)                          # E19(a)
    expect_true(eR$rm$artefact <= eR$rm$M_implied + 1e-9)                                            # E19(c2)
  }
})

test_that("Def. 5.1: eiATT refuses a causal label without the required statements", {
  f <- rnd(1)
  expect_error(eiATT(f, treated = "g1", comparison = "g2", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE), "population")
  expect_error(eiATT(f, treated = "g1", comparison = "g2", population = "fixed", spillover = 0, no_anticipation = TRUE), "trend assumption")
  expect_error(eiATT(f, treated = "g1", population = "fixed", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE), "comparison group")
  expect_error(eiATT(f, treated = "g1", population = "fixed", trend = mtr_nd_ps(), spillover = 0, no_anticipation = TRUE), "comparison group")
  expect_error(eiATT(f, treated = "g1", comparison = "g2", population = "fixed", trend = parallel_trends(), no_anticipation = TRUE), "spillover")
  expect_error(eiATT(f, treated = "g1", comparison = "g2", population = "fixed", trend = parallel_trends(), spillover = 0), "no anticipation")
  expect_error(eiATT(f, treated = "g1", comparison = "g2", population = "fixed", trend = parallel_trends(), spillover = 0, no_anticipation = FALSE), "no anticipation")
  expect_error(eiATT(f, comparison = "g2", population = "fixed", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE), "treated group")
  expect_error(eiATT(f, treated = "g1", comparison = "g2", population = "stayers", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE), "population")
  # MTR + ND alone need no comparison group (E15 remark)
  e <- eiATT(f, treated = "g1", population = "fixed", trend = mtr_nd_ps(ps = FALSE), spillover = 0, no_anticipation = TRUE)
  expect_true(is.na(e$results[[1]]$lower) || abs(e$results[[1]]$lower) < 1e-12)
  expect_output(print(e), "CAUSAL")
  expect_output(print(f), "DESCRIPTIVE")
})

test_that("Prop. E22: the spillover allowance widens the set; it is always printed", {
  f <- rnd(3)
  e0 <- A1(f, parallel_trends()); e1 <- A1(f, parallel_trends(), L = 0.05)
  expect_true(e1$lower <= e0$lower + 1e-12 && e1$upper >= e0$upper - 1e-12)
  a <- eiATT(f, treated = "g1", comparison = "g2", population = "fixed", trend = parallel_trends(), spillover = 0.05, no_anticipation = TRUE)
  expect_output(print(a), "no information on spillover")
})

sim3 <- function(p = 200, seed = 3) {
  set.seed(seed)
  x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
  t0 <- .3 + .2 * x + rnorm(p, 0, .03); t1 <- t0 + .02 + .1 * x + rnorm(p, 0, .02); t2 <- t1 + .05 + .15 * x + rnorm(p, 0, .02)
  mkT <- function(l) round(N * pmin(pmax(l, 0), 1))
  mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = mkT(t), no = N - mkT(t))
  list(m1 = mk(t0), b = mk(t1), p = mk(t2))
}

test_that("model layer: FD (E7, E9, E23), per-period LCE (E11), cross-layer (E13)", {
  dl <- sim3()
  f <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                model = fd(rows = list(list(type = "didS"))))
  m <- f$model$table
  e9 <- f$model$extra[[2]]$e9
  expect_equal(e9$interval, c(m$plugin_lower[2], m$plugin_upper[2]), tolerance = 1e-8)   # E9(a) closed form = LP
  expect_true(all(m$model_lower <= m$plugin_lower + 1e-9 & m$model_upper >= m$plugin_upper - 1e-9))
  expect_true(all(m$lower >= m$count_lower - 1e-12 & m$upper <= m$count_upper + 1e-12))  # intersection
  expect_true(all(grepl("pooled FD", m$status)))
  # R1 point-identifies the model DiD (2d FD2)
  f1 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                 model = fd(rows = list(list(type = "R1"))))
  expect_equal(f1$model$table$plugin_lower, f1$model$table$plugin_upper, tolerance = 1e-9)
  # the allowance for a DiD is rho_i = 2 (w_ir + w_is) (Lemma E12)
  A <- dl$b$v; O <- dl$b$o
  S <- 0.5 * sqrt(sum((2 * (A / sum(A) + O / sum(O)))^2))
  expect_equal(f1$model$table$S_phi[2], S, tolerance = 1e-10)
  # non-consecutive pairs are refused (E25 (H1))
  expect_error(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                        model = fd(pairs = list(c("m1", "p")))), "consecutive")
  # per-period LCE with model CR on changes: union of two branches; pooled status
  f3 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                 model = lce_periods(rows = list(list(type = "cr", from = "b", to = "p"))))
  expect_equal(nrow(f3$model$extra[[2]]$branches), 2)
  expect_true(all(grepl("Prop. E11", f3$model$table$status)))
  # level targets have no model value under FD
  f4 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", integer = "lp",
                 targets = list(rate_at("v", "yes", "p")), model = fd())
  expect_match(f4$model$table$note, "levels free")
})

test_that("Prop. E25: model-based pre-trends; FD alone is OUTER (no rate limits), LCE imposes them", {
  dl <- sim3()
  fF <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                 model = fd(rows = list(list(type = "R1"))))
  aF <- eiATT(fF, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(1), spillover = 0,
              no_anticipation = TRUE, pretrends = "model")
  expect_true(aF$model$outer)
  expect_output(print(aF), "OUTER")
  expect_equal(nrow(aF$results[[1]]$branches), 2)
  # the box radii are Lemma E12 radii; the projection share and the box shares add to alpha
  expect_equal(aF$model$alpha_proj + length(aF$model$radii) * aF$model$alpha_q, 0.05, tolerance = 1e-12)
  fL <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                 model = lce_periods())
  aL <- eiATT(fL, treated = "v", comparison = "o", population = "fixed", trend = trend_dial(0.02), spillover = 0,
              no_anticipation = TRUE, pretrends = "model")
  expect_false(aL$model$outer)
  e <- aL$results[[1]]
  expect_true(e$lower >= -1 && e$upper <= 1)                    # clipped to the ATT's logical range
  # a larger dial can only widen the set
  aL2 <- eiATT(fL, treated = "v", comparison = "o", population = "fixed", trend = trend_dial(0.1), spillover = 0,
               no_anticipation = TRUE, pretrends = "model")
  expect_true(aL2$results[[1]]$lower <= e$lower + 1e-9 && aL2$results[[1]]$upper >= e$upper - 1e-9)
  expect_error(eiATT(rnd(1), treated = "g1", comparison = "g2", population = "fixed", trend = parallel_trends(), spillover = 0,
                     no_anticipation = TRUE, pretrends = "model"), "model layer")
})

test_that("gap register: changing populations, units and bases are refused; estimated margins warn", {
  d <- two_el(cbind(c(3, 4), c(2, 3)), c(2, 3), c(3, 4))
  d2 <- d; d2$e1$g1[1] <- 4; d2$e1$no[1] <- d2$e1$no[1] + 1
  expect_error(eiChange(d2, cbind(yes, no) ~ cbind(g1, g2), unit = "u"), "G1")
  d3 <- d; d3$e1$u[2] <- "Z"
  expect_error(eiChange(d3, cbind(yes, no) ~ cbind(g1, g2), unit = "u"), "G6")
  d4 <- d; d4$e1$no[1] <- d4$e1$no[1] - 1
  expect_error(eiChange(d4, cbind(yes, no) ~ cbind(g1, g2), unit = "u"), "common base")
  expect_warning(f <- eiChange(d, cbind(yes, no) ~ cbind(g1, g2), unit = "u", margins = "estimated"), "G7")
  expect_false(f$P$integer)
  expect_match(f$table$integer_status, "not applicable")
  d5 <- d; d5$e0$yes <- d5$e0$yes + 0.5; d5$e0$no <- d5$e0$no - 0.5
  f5 <- eiChange(d5, cbind(yes, no) ~ cbind(g1, g2), unit = "u")
  expect_false(f5$P$integer)
  expect_output(print(f5), "G7")
  expect_error(eiChange(d, cbind(yes, no) ~ cbind(g1, g2), unit = "u",
                        restrictions = list(sign_change = list(row = "g1", col = "yes", from = "e0", to = "e1", object = "model"))), "model layer")
})
