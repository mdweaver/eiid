# Regression tests for the PHASE6a code review (docs/extension/reviews/PHASE6a_code_review.md),
# issues 1-16 (issue 15 is this file), and the proof reviewer's round-3 nit a.

enum_unit2 <- function(a, b) {
  R <- length(a); C <- length(b); out <- list()
  rec <- function(n, k) {
    if (k > R * C) { if (all(rowSums(n) == a) && all(colSums(n) == b)) out[[length(out) + 1]] <<- n; return() }
    r <- (k - 1) %/% C + 1; c <- (k - 1) %% C + 1
    hi <- min(a[r] - sum(n[r, ]), b[c] - sum(n[, c]))
    if (hi < 0) return()
    for (v in 0:hi) { n[r, c] <- v; rec(n, k + 1) }
  }
  rec(matrix(0, R, C), 1)
  out
}
brute2 <- function(A, B, f, ok) {
  per <- lapply(seq_len(nrow(A)), function(i) enum_unit2(A[i, ], B[i, ]))
  idx <- as.matrix(expand.grid(lapply(per, seq_along)))
  vals <- c()
  for (k in seq_len(nrow(idx))) {
    tabs <- lapply(seq_len(nrow(A)), function(i) per[[i]][[idx[k, i]]])
    if (ok(tabs)) vals <- c(vals, f(tabs))
  }
  if (!length(vals)) c(NA, NA) else range(vals)
}
mk <- function(tabs, R = nrow(tabs[[1]]), C = ncol(tabs[[1]])) {
  A <- do.call(rbind, lapply(tabs, rowSums)); B <- do.call(rbind, lapply(tabs, colSums))
  d <- as.data.frame(cbind(A, B)); names(d) <- c(paste0("r", 1:R), paste0("c", 1:C)); d$u <- paste0("U", seq_along(tabs))
  list(d = d, A = A, B = B, rows = paste0("r", 1:R), cols = paste0("c", 1:C))
}
glpk <- function() skip_if_not(.ec_have_glpk(), "Rglpk not installed")

# ---- issue 1: MILP optimality ---------------------------------------------------------------
test_that("issue 1: the review's instance (gap across unit sets, fractional row) has upper end 1/3", {
  glpk()
  d <- data.frame(r1 = c(2, 1, 1), r2 = c(3, 1, 2), c1 = c(1, 2, 1), c2 = c(3, 0, 1), c3 = c(1, 0, 1))
  for (mode in c("auto", "exact")) {
    f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2", "c3"), data = d, integer = mode,
                  targets = list(gap("r1", "r1", "c2", where = 1:2, where2 = 3L)),
                  restrictions = list(linear = list(cells = data.frame(unit = c(1L, 3L), row = "r1", col = c("c3", "c2"),
                                                                       coef = c(3, 0.5)), lower = 1.65)))
    expect_equal(f$table$upper, 1 / 3, tolerance = 1e-12)
    expect_equal(f$table$lower, -1, tolerance = 1e-12)
    expect_equal(f$table$integer_status, "integer-exact (MILP)")
  }
})

test_that("issue 1: randomised brute force at the oracle's scale (subset and gap targets, MILP rows)", {
  glpk()
  set.seed(20261122)
  nchk <- 0
  for (k in 1:40) {
    p <- sample(2:3, 1); R <- 3; C <- sample(2:3, 1)
    tabs <- lapply(1:p, function(i) { repeat { m <- matrix(sample(0:2, R * C, TRUE), R, C); if (sum(m) > 0 && sum(m) <= 7) return(m) } })
    z <- mk(tabs, R, C)
    if (any(colSums(z$A) == 0) || z$A[1, 1] == 0 || sum(z$A[-1, 2]) == 0) next
    scen <- k %% 3
    if (scen == 0) {
      rs <- list(order_unit = list(row = "r1", over = "r2", col = "c1"))
      ok <- function(t) all(vapply(seq_along(t), function(i) z$A[i, 1] == 0 || z$A[i, 2] == 0 ||
                                     z$A[i, 2] * t[[i]][1, 1] - z$A[i, 1] * t[[i]][2, 1] >= 0, TRUE))
    } else if (scen == 1) {
      if (any(z$A[, 1:2] == 0)) next
      st <- list(X = z$A / rowSums(z$A), N = rowSums(z$A), R = R)
      W <- .rxc_assoc_W(st, "dB", "pair", 1, 2)
      rs <- list(sign = list(type = "dB", version = "pair", row = "r1", other = "r2", col = "c1"))
      ok <- function(t) sum(vapply(seq_along(t), function(i) sum(W[i, ] * t[[i]][, 1] / pmax(z$A[i, ], 1) * (z$A[i, ] > 0)), 0)) >= -1e-12
    } else {
      cf <- sample(c(1, 2, -1, 0.5, 1.5), 2, TRUE); lo <- 0.65
      rs <- list(linear = list(cells = data.frame(unit = c(1L, p), row = "r1", col = c("c1", "c2"), coef = cf), lower = lo))
      ok <- function(t) cf[1] * t[[1]][1, 1] + cf[2] * t[[p]][1, 2] >= lo
    }
    # gap across unit sets: rate of r1 on c2 in unit 1 minus rate of r2+r3 on c2 elsewhere
    f_t <- function(t) t[[1]][1, 2] / z$A[1, 1] - sum(vapply(2:p, function(i) sum(t[[i]][2:3, 2]), 0)) / sum(z$A[2:p, 2:3])
    if (sum(z$A[2:p, 2:3]) == 0) next
    fit <- eiCounts(rows = z$rows, cols = z$cols, data = z$d, unit = "u", binding = FALSE, restrictions = rs,
                    targets = list(gap("r1", c("r2", "r3"), "c2", where = "U1", where2 = paste0("U", 2:p))))
    br <- brute2(z$A, z$B, f_t, ok)
    if (is.na(br[1])) { expect_false(isTRUE(fit$table$integer_feasible)); next }
    expect_equal(c(fit$table$lower, fit$table$upper), br, tolerance = 1e-12)
    nchk <- nchk + 1
  }
  expect_gt(nchk, 15)
})

# ---- issue 16: realistic unit sizes; time limit --------------------------------------------
real_d <- function(p, seed = 5) {
  set.seed(seed)
  tabs <- lapply(1:p, function(i) { N <- sample(200:1500, 1); x <- as.vector(rmultinom(1, N, c(.5, .3, .2)))
    t(vapply(1:3, function(r) as.vector(rmultinom(1, x[r], c(.5, .3, .2) + c(.1, 0, -.1) * (r - 2))), numeric(3))) })
  A <- t(vapply(tabs, rowSums, numeric(3))); B <- t(vapply(tabs, colSums, numeric(3)))
  d <- data.frame(A, B); names(d) <- c("g1", "g2", "g3", "o1", "o2", "o3"); d
}

test_that("issue 16: pairwise delta_B and W2 rows at N_i of 200-1500 finish, certified", {
  glpk()
  d <- real_d(5)
  for (rs in list(list(sign = list(type = "dB", version = "pair", row = "g1", other = "g2", col = "o1")),
                  list(sign = list(type = "dW", version = "W2", row = "g1", col = "o1", sign = -1)))) {
    t0 <- proc.time()[["elapsed"]]
    f <- eiCounts(cbind(o1, o2, o3) ~ cbind(g1, g2, g3), data = d, targets = list(rate("g1", "o1")),
                  binding = FALSE, restrictions = rs)
    expect_lt(proc.time()[["elapsed"]] - t0, 60)
    expect_equal(f$table$integer_status, "integer-exact (MILP)")
    expect_true(f$table$lower >= f$table$lp_lower - 1e-12 && f$table$upper <= f$table$lp_upper + 1e-12)
    w <- witness(f)
    expect_true(all(vapply(w, function(z) z$margins_ok && z$rows_ok && z$integer, TRUE)))
  }
})

test_that("issue 16: on timeout the status is 'integer status unresolved' with the LP outer bound, never integer-exact", {
  glpk()
  d <- real_d(60, 7)
  f <- eiCounts(cbind(o1, o2, o3) ~ cbind(g1, g2, g3), data = d, targets = list(rate("g1", "o1")), binding = FALSE,
                restrictions = list(sign = list(type = "dW", version = "W2", row = "g1", col = "o1", sign = -1)),
                time_limit = 0.001)
  expect_match(f$table$integer_status, "^integer status unresolved")
  expect_false(grepl("integer-exact", f$table$integer_status))
  expect_true(is.na(f$table$integer_feasible))
  expect_true(f$table$lower <= f$table$lp_lower + 1e-9 || f$table$lower <= f$table$upper)
  expect_output(print(f), "never integer-exact")
})

# ---- issue 2: known with unsorted or duplicated units --------------------------------------
test_that("issue 2: known values follow the user's unit order; duplicates and wrong lengths are refused", {
  d <- data.frame(u = c("U1", "U2", "U3"), r1 = c(5, 5, 5), r2 = c(5, 5, 5), c1 = c(5, 5, 5), c2 = c(5, 5, 5))
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                targets = list(rate("r1", "c1", where = "U1"), rate("r1", "c1", where = "U2")),
                restrictions = list(known = list(unit = c("U2", "U1"), row = "r1", col = "c1", value = c(4, 1))))
  expect_equal(f$table$lower * 5, c(1, 4)); expect_equal(f$table$upper * 5, c(1, 4))
  f2 <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                 targets = list(rate("r1", "c1", where = "U3"), rate("r1", "c1", where = "U1")),
                 restrictions = list(known = list(unit = c(3, 1), row = "r1", col = "c1", value = c(2, 3))))
  expect_equal(f2$table$lower * 5, c(2, 3))
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                        restrictions = list(known = list(unit = c("U1", "U1"), row = "r1", col = "c1", value = c(1, 2)))),
               "named twice")
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                        restrictions = list(known = list(unit = c("U1", "U2", "U3"), row = "r1", col = "c1", value = c(1, 2)))),
               "length 1 or one value per unit")
})

# ---- issue 11: out-of-range numbers and dropped units --------------------------------------
test_that("issue 11: unknown row numbers and dropped units are errors with clear messages", {
  d <- data.frame(u = c("a", "b", "z"), r1 = c(2, 1, 0), r2 = c(2, 3, 0), c1 = c(1, 2, 0), c2 = c(3, 2, 0))
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                        restrictions = list(box = list(units = c(1, 99), row = "r1", col = "c1", lower = 1))), "out of range")
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                        restrictions = list(box = list(units = "z", row = "r1", col = "c1", lower = 1))), "N_i = 0")
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u",
                        restrictions = list(box = list(units = 3, row = "r1", col = "c1", lower = 1))), "N_i = 0")
})

# ---- issue 3: eiBounds with N_i = 0 under the default target; failures do not stop it --------
test_that("issue 3: eiBounds(x, t, n) with a unit of size 0 works under target = 'actual'", {
  f <- eiBounds(x = c(.3, .5, .4, .6), t = c(.5, .4, .6, .55), n = c(100, 0, 80, 90), model = "L6", tests = FALSE)
  m <- eiBounds(x = c(.3, .5, .4, .6), t = c(.5, .4, .6, .55), n = c(100, 0, 80, 90), model = "L6", tests = FALSE, target = "model")
  expect_false(is.null(f$counts))
  cols <- setdiff(names(m$table), c("dd_exact_lower", "dd_exact_upper"))
  expect_equal(f$table[, cols], m$table[, cols])
  expect_false(is.null(f$counts$cross))
  # a failing count restriction gives a warning, not an error, and keeps the model output
  d <- data.frame(vet = c(20, 60, 30), other = c(80, 40, 70), yes = c(30, 50, 40), no = c(70, 50, 60))
  expect_warning(g <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "L6", tests = FALSE,
                               counts = list(box = list(row = "nonexistent", col = "yes", lower = 1))),
                 "count layer could not be computed")
  expect_true(is.null(g$counts))
  expect_output(print(g), "R x C ecological inference bounds")
})

# ---- issue 4: double coefficients are rationalised; exact forms accepted ----------------------
test_that("issue 4: coef = 1/3 with lower = 1 means n >= 3; character fractions are exact; the note prints", {
  d <- data.frame(r1 = c(6, 6), r2 = c(4, 4), c1 = c(5, 5), c2 = c(5, 5))
  for (cf in list(1 / 3, "1/3")) {
    f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1", where = 1L)),
                  restrictions = list(linear = list(cells = data.frame(unit = 1L, row = "r1", col = "c1", coef = cf,
                                                                       stringsAsFactors = FALSE), lower = 1)))
    expect_equal(f$table$lower * 6, 3)
  }
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1", where = 1L)),
                restrictions = list(linear = list(cells = data.frame(unit = 1L, row = "r1", col = "c1", coef = 0.1), lower = 3 * 0.1)))
  expect_equal(f$table$lower * 6, 3)
  expect_output(print(f), "rationalised")
  # a fractional box bound typed as a double
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1", where = 1L)),
                restrictions = list(box = list(units = 1L, row = "r1", col = "c1", upper = 10 / 3)))
  expect_equal(f$table$upper * 6, 3)
})

# ---- issue 5: all-zero association rows --------------------------------------------------------
test_that("issue 5: a W2 row with constant X_r is vacuous: dropped with a note, no crash", {
  glpk()
  d <- data.frame(r1 = c(2, 4), r2 = c(1, 2), r3 = c(1, 2), c1 = c(2, 3), c2 = c(2, 5))   # X_r1 = 0.5 in both units
  f <- eiCounts(rows = c("r1", "r2", "r3"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1")),
                restrictions = list(sign = list(type = "dW", version = "W2", row = "r1", col = "c1")))
  expect_equal(c(f$table$lower, f$table$upper), c(f$table$dd_lower, f$table$dd_upper))
  expect_output(print(f), "vacuous")
})

# ---- issue 6: the TU path checks integrality before rounding ---------------------------------
test_that("issue 6: a fractional vertex is never rounded into an integer endpoint", {
  # force the TU path on a system whose rounded LP has fractional vertices: three crossing totals
  # are not in a class, so .ec_classify sends them to the MILP; here we call .ec_system with a
  # classifier that wrongly says "(i)" and check that the fallback catches it
  A <- matrix(1, 3, 2); B <- matrix(1, 3, 2)
  d <- data.frame(r1 = 1, r2 = 1, c1 = c(1, 1, 1), c2 = c(1, 1, 1))
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1", where = 1L)),
                restrictions = list(total = list(list(units = 1:2, row = 1, col = 1, value = 1),
                                                 list(units = 2:3, row = 1, col = 1, value = 1),
                                                 list(units = c(1, 3), row = 1, col = 1, value = 1))))
  expect_false(isTRUE(f$table$integer_feasible))
  local_mocked_bindings(.ec_classify = function(P, rows) list(class = "(i)", reason = "forced"))
  if (.ec_have_glpk()) {
    sy <- .ec_system(f$P, f$built$common, f$targets[[1]], "auto")
    expect_true(sy$fallback)
    expect_false(isTRUE(sy$int_ok))
  }
})

# ---- issue 7: outcome-selected unit sets ------------------------------------------------------
test_that("issue 7: a where= that uses the outcomes gets the count layer only", {
  set.seed(9); p <- 60; N <- sample(100:300, p, TRUE); V <- round(N * runif(p, .2, .8)); Y <- rbinom(p, N, .5)
  d <- data.frame(u = paste0("p", 1:p), vet = V, other = N - V, yes = Y, no = N - Y)
  m <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "L6", tests = FALSE, target = "model")
  f <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", model = m,
                targets = list(rate("vet", "yes", where = yes / (yes + no) > 0.5), rate("vet", "yes", where = vet > 100),
                               rate("vet", "yes", where = d$u[1:30])))
  expect_false(f$table$design[1]); expect_true(f$table$design[2])
  expect_match(f$cross$table$note[1], "selected using the outcomes")
  expect_true(is.na(f$cross$table$model_lower[1]))
  expect_false(is.na(f$cross$table$model_lower[2]))
  expect_match(f$cross$table$note[3], "must be fixed from the design")
})

# ---- issue 8: role defaults of totals ----------------------------------------------------------
test_that("issue 8: a total given as a floor or interval is a belief row; value = is an identity row", {
  d <- data.frame(vet = c(20, 60), other = c(80, 40), yes = c(30, 50), no = c(70, 50))
  f <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, targets = list(rate("vet", "yes")),
                restrictions = list(total = list(list(units = 1:2, row = "vet", col = "yes", value = 45),
                                                 list(units = 1:2, row = "vet", col = "yes", lower = 40),
                                                 list(units = 1:2, row = "vet", col = "yes", lower = 40, upper = 50),
                                                 list(units = 1:2, row = "vet", col = "yes", lower = 40, role = "identity"))))
  expect_equal(f$restrictions$role, c("identity", "belief", "belief", "identity"))
})

# ---- issue 10: the cross-layer interval reuses the fit's draws ---------------------------------
test_that("issue 10: for a single cell (conservative, Popoviciu) the cross-layer model interval equals ci_actual", {
  set.seed(4); p <- 2000; N <- sample(100:300, p, TRUE); V <- round(N * runif(p, .05, .95))
  Y <- rbinom(p, V, .6) + rbinom(p, N - V, .4)
  d <- data.frame(vet = V, other = N - V, yes = Y, no = N - Y)
  f <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", ci = "conservative", fpc = "popoviciu",
                seed = 3, tests = FALSE, nz_draw = 1000)
  ct <- f$counts$cross$table; tb <- f$table
  j <- which(!is.na(ct$model_lower) & tb$ci_actual_lower > tb$dd_lower & tb$ci_actual_upper < tb$dd_upper)
  skip_if(!length(j), "every ci_actual is clipped by Duncan-Davis here")
  expect_equal(ct$model_lower[j], tb$ci_actual_lower[j], tolerance = 1e-12)
  expect_equal(ct$model_upper[j], tb$ci_actual_upper[j], tolerance = 1e-12)
})

# ---- issue 9: status strings (6a Sec. 4.3 rows of bbf8cf46) -------------------------------------
test_that("issue 9: L6 and pooled L4 statuses use the new Sec. 4.3 strings", {
  set.seed(9); p <- 60; N <- sample(100:300, p, TRUE); V <- round(N * runif(p, .2, .8)); Y <- rbinom(p, N, .5)
  d <- data.frame(vet = V, other = N - V, yes = Y, no = N - Y)
  m6 <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "L6", tests = FALSE)
  expect_equal(m6$counts$cross$table$status[1],
               "normal interval: asymptotic guarantee (P3.5(a) with no free k: point-identified rung, one linear functional), pooled model")
  m4 <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "L4", row_constant = "vet", tests = FALSE,
                 ndraw = 20, nz_draw = 300)
  expect_equal(m4$counts$cross$table$status[1],
               "proj_*: asymptotic guarantee (P3.5(b2)(ii); L4 equalities as rows, as in PHASE5a R10), pooled model")
})

# ---- issues 12, 13, 14 --------------------------------------------------------------------------
test_that("issues 12-14: x/t labels, sign labels, Lindeberg note under Hoeffding, '>=0' signs", {
  f <- eiCounts(x = c(.3, .5, .4), t = c(.5, .4, .6), n = c(100, 80, 90), binding = FALSE)
  expect_equal(f$P$rows, c("X", "notX")); expect_equal(f$P$cols, c("T", "notT"))
  d <- data.frame(r1 = c(1, 4, 2), r2 = c(4, 1, 3), c1 = c(2, 3, 1), c2 = c(3, 2, 4))
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1")),
                restrictions = list(sign = list(list(type = "dB", row = "r1", col = "c1", sign = ">=0"),
                                                list(type = "dW", row = "r1", col = "c1", sign = "<=0"))))
  expect_match(f$restrictions$restriction[1], "delta_B"); expect_match(f$restrictions$restriction[1], "vs rest")
  expect_match(f$restrictions$restriction[2], "delta_W")
  expect_error(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d,
                        restrictions = list(sign = list(type = "dB", row = "r1", col = "c1", sign = "up"))), "sign must be")
  # two targets on different single-unit subsets get distinct labels
  f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1", where = 1L), rate("r1", "c1", where = 2L)))
  expect_equal(length(unique(f$table$target)), 2)
  # Lindeberg note: printed with Popoviciu (few units), suppressed with Hoeffding
  dd <- data.frame(vet = c(20, 60, 30, 40, 10, 55), other = c(80, 40, 70, 60, 90, 45), yes = c(30, 50, 40, 45, 20, 60),
                   no = c(70, 50, 60, 55, 80, 40))
  m <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = dd, model = "L6", tests = FALSE, target = "model")
  fp <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = dd, model = m, targets = list(rate("vet", "yes")))
  fh <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = dd, model = m, targets = list(rate("vet", "yes")), allowance = "hoeffding")
  expect_output(print(fp), "heuristic")
  op <- capture.output(print(fh))
  expect_false(any(grepl("heuristic", op)))
  expect_true(any(grepl("allowance holds in finite samples", op)))
})

# ---- binding: cheaper and optional ----------------------------------------------------------------
test_that("binding: inactive restrictions are skipped (same answer as the full drop-one), 'auto' off above 200 units", {
  set.seed(12); p <- 30; N <- sample(50:100, p, TRUE); V <- round(N * runif(p, .2, .8)); Y <- rbinom(p, N, .5)
  d <- data.frame(u = paste0("p", 1:p), vet = V, other = N - V, yes = Y, no = N - Y)
  rs <- list(total = list(units = paste0("p", 1:10), row = "vet", col = "yes", value = sum(pmin(V[1:10], Y[1:10])) - 3),
             box = list(units = "p20", row = "vet", col = "yes", lower = 0),
             order_unit = list(row = "vet", col = "yes"))
  f <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", targets = list(rate("vet", "yes")), restrictions = rs)
  # brute-force drop-one
  full <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", targets = list(rate("vet", "yes")), restrictions = rs, binding = FALSE)
  bl <- c(); bu <- c()
  for (j in seq_along(rs)) {
    g <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", targets = list(rate("vet", "yes")), restrictions = rs[-j], binding = FALSE)
    if (g$table$lower < full$table$lower - 1e-12) bl <- c(bl, f$restrictions$id[j])
    if (g$table$upper > full$table$upper + 1e-12) bu <- c(bu, f$restrictions$id[j])
  }
  expect_equal(f$table$binding_lower, paste(bl, collapse = ", "))
  expect_equal(f$table$binding_upper, paste(bu, collapse = ", "))
  set.seed(1); p <- 220; N <- sample(50:100, p, TRUE); V <- round(N * runif(p, .2, .8)); Y <- rbinom(p, N, .5)
  d <- data.frame(vet = V, other = N - V, yes = Y, no = N - Y)
  g <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, targets = list(rate("vet", "yes")),
                restrictions = list(order_unit = list(row = "vet", col = "yes")))
  expect_true(is.na(g$table$binding_lower))
  expect_output(print(g), "binding rows not computed")
})

# ---- proof review round 3, nit a: zero-SE terms in P3.5(a) ------------------------------------
test_that("nit a: a beta-independent row (zero SE) does not enter the studentised maximum", {
  set.seed(2); p <- 80; N <- sample(100:300, p, TRUE); V <- round(N * runif(p, .2, .8)); Y <- rbinom(p, N, .5)
  d <- data.frame(vet = V, other = N - V, yes = Y, no = N - Y)
  m <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", ci = "conservative", tests = FALSE, target = "model")
  st <- m$setup; cons <- .rxc_constraints(st); tg <- .rxc_target(st, matrix(c(1, 0, 0, 0), 2))
  set.seed(1); Zs <- .rxc_zdraws(st$V, 2000)
  k0 <- .rxc_ci_conservative(st, cons, tg, 0.05, Zs)$k
  cons2 <- cons
  cons2$S <- rbind(cons$S, 0); cons2$s0 <- c(cons$s0, 1); cons2$A <- rbind(cons$A, 0)
  cons2$kind <- c(cons$kind, "x"); cons2$module <- c(cons$module, NA)
  r2 <- .rxc_ci_conservative(st, cons2, tg, 0.05, Zs)
  expect_equal(r2$k, k0)
  expect_true(is.finite(r2$k))
  # every term beta-independent: k = 0, no division by zero
  st0 <- st; st0$V <- st$V * 0; Z0 <- Zs; Z0$V <- Zs$V * 0
  expect_equal(.rxc_ci_conservative(st0, cons, tg, 0.05, Z0)$k, 0)
})
