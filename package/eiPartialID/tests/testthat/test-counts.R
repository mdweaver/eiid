# The count layer (M-COUNT): eiCounts(), witness(), the cross-layer block and the eiBounds()
# embedding. Theory and check names: docs/extension/PHASE6a_count_layer.md (contract C10).

# ---- helpers ---------------------------------------------------------------------------
# all nonnegative integer R x C tables with row sums a and column sums b
enum_unit <- function(a, b) {
  R <- length(a); C <- length(b); out <- list()
  rec <- function(n, k) {
    if (k > R * C) { if (all(rowSums(n) == a) && all(colSums(n) == b)) out[[length(out) + 1]] <<- n; return() }
    r <- (k - 1) %/% C + 1; c <- (k - 1) %% C + 1
    hi <- min(a[r] - sum(n[r, ]), b[c] - sum(n[, c]))
    if (c == C) { v <- a[r] - sum(n[r, ]); if (v >= 0 && v <= b[c] - sum(n[, c])) { n[r, c] <- v; rec(n, k + 1) }; return() }
    if (r == R) { v <- b[c] - sum(n[, c]); if (v >= 0 && v <= hi) { n[r, c] <- v; rec(n, k + 1) }; return() }
    for (v in 0:max(hi, -1)) { n[r, c] <- v; rec(n, k + 1) }
  }
  rec(matrix(0, R, C), 1)
  out
}
# brute force: range of target f(tabs) over all integer tables with ok(tabs) TRUE
brute <- function(A, B, f, ok = function(tabs) TRUE) {
  per <- lapply(seq_len(nrow(A)), function(i) enum_unit(A[i, ], B[i, ]))
  idx <- as.matrix(expand.grid(lapply(per, seq_along)))
  vals <- c()
  for (k in seq_len(nrow(idx))) {
    tabs <- lapply(seq_len(nrow(A)), function(i) per[[i]][[idx[k, i]]])
    if (ok(tabs)) vals <- c(vals, f(tabs))
  }
  if (!length(vals)) return(c(NA, NA))
  range(vals)
}
mkdata <- function(A, B) {
  d <- as.data.frame(cbind(A, B)); R <- ncol(A); C <- ncol(B)
  names(d) <- c(paste0("r", seq_len(R)), paste0("c", seq_len(C))); d$u <- paste0("U", seq_len(nrow(A))); d
}
fitc <- function(A, B, ...) {
  d <- mkdata(A, B)
  eiCounts(rows = paste0("r", seq_len(ncol(A))), cols = paste0("c", seq_len(ncol(B))), data = d, unit = "u", ...)
}
rand_tabs <- function(I, R, C, maxv = 3) lapply(seq_len(I), function(i) {
  repeat { m <- matrix(sample(0:maxv, R * C, TRUE), R, C); if (sum(m) > 0) return(m) } })
AB <- function(tabs) list(A = do.call(rbind, lapply(tabs, rowSums)), B = do.call(rbind, lapply(tabs, colSums)))
wit_ok <- function(fit) {
  w <- witness(fit)
  all(vapply(w, function(z) z$margins_ok && z$rows_ok && abs(z$value - z$reported) < 1e-9 &&
               (!grepl("integer-exact", fit$table$integer_status[match(z$target, fit$table$target)]) || z$integer), TRUE))
}

# ---- D1 / D1b: margins only = Duncan-Davis ---------------------------------------------
test_that("D1b: margins-only range reproduces the package's Duncan-Davis to 1e-12 (cells, one-column contrasts)", {
  for (f in c("rxc_example_3x3.csv", "msign_example_3x2.csv")) {
    d <- read.csv(test_path(f))
    X <- as.matrix(d[, grep("^X", names(d))]); T <- as.matrix(d[, grep("^T", names(d))])
    R <- ncol(X); C <- ncol(T)
    dd <- .rxc_dd_cells(X, T, d$N)
    tg <- list(); ref <- list()
    for (r in 1:R) for (c in 1:C) { a <- matrix(0, R, C); a[r, c] <- 1; tg[[paste(r, c)]] <- a
      ref[[paste(r, c)]] <- c(dd$lower[r, c], dd$upper[r, c]) }
    for (c in 1:C) { a <- matrix(0, R, C); a[1, c] <- 1; a[2, c] <- -1; tg[[paste("g", c)]] <- a
      ref[[paste("g", c)]] <- .rxc_interval_lin(a, dd$lower, dd$upper) }
    fit <- eiCounts(x = X, t = T, n = d$N, targets = tg, binding = FALSE)
    expect_equal(unname(fit$table$dd_lower), unname(vapply(ref, `[`, 0, 1)), tolerance = 1e-12)
    expect_equal(unname(fit$table$dd_upper), unname(vapply(ref, `[`, 0, 2)), tolerance = 1e-12)
    expect_true(all(fit$table$integer_status == "not applicable: non-integer margins"))
    # the general LP path gives the same (two-row contrast over two columns forces the LP)
    a <- matrix(0, R, C); a[1, 1] <- 1; a[2, 1] <- -1
    expect_true(wit_ok(eiCounts(x = X, t = T, n = d$N, targets = list(g = a), binding = FALSE)))
  }
})

test_that("D1: random integer data, single cells = Frechet sums; LP path agrees", {
  set.seed(11)
  for (k in 1:20) {
    tabs <- rand_tabs(sample(2:4, 1), sample(2:3, 1), sample(2:3, 1), 6); ab <- AB(tabs)
    R <- ncol(ab$A); C <- ncol(ab$B); N <- rowSums(ab$A)
    r <- sample(R, 1); c <- sample(C, 1)
    if (sum(ab$A[, r]) == 0) next
    fit <- fitc(ab$A, ab$B, targets = list(rate(r, c)), binding = FALSE)
    expect_equal(fit$table$dd_lower, sum(pmax(0, ab$A[, r] + ab$B[, c] - N)) / sum(ab$A[, r]), tolerance = 1e-14)
    expect_equal(fit$table$dd_upper, sum(pmin(ab$A[, r], ab$B[, c])) / sum(ab$A[, r]), tolerance = 1e-14)
    # through the LP: a nonbinding box forces the solver
    f2 <- fitc(ab$A, ab$B, targets = list(rate(r, c)), binding = FALSE,
               restrictions = list(box = list(row = 1, col = 1, lower = 0)))
    expect_equal(f2$table$lower, fit$table$dd_lower, tolerance = 1e-12)
    expect_equal(f2$table$upper, fit$table$dd_upper, tolerance = 1e-12)
  }
})

test_that("Lemma 2.1(c) / RC6: the one-vs-rest gap is exact = single-cell DD mapped, narrower than cellwise", {
  set.seed(3)
  nwide <- 0
  for (k in 1:15) {
    tabs <- rand_tabs(sample(2:4, 1), 3, 2, 8); ab <- AB(tabs)
    if (any(colSums(ab$A) == 0)) next
    fit <- fitc(ab$A, ab$B, targets = list(ovr_gap(1, 1), rate(1, 1)), binding = FALSE)
    N <- sum(ab$A); Ar <- sum(ab$A[, 1]); Yc <- sum(ab$B[, 1]); kap <- N / (Ar * (N - Ar))
    Ylo <- fit$table$dd_lower[2] * Ar; Yhi <- fit$table$dd_upper[2] * Ar
    expect_equal(fit$table$dd_lower[1], kap * Ylo - Yc / (N - Ar), tolerance = 1e-12)
    expect_equal(fit$table$dd_upper[1], kap * Yhi - Yc / (N - Ar), tolerance = 1e-12)
    X <- ab$A / rowSums(ab$A); T <- ab$B / rowSums(ab$B); dd <- .rxc_dd_cells(X, T, rowSums(ab$A))
    a <- matrix(0, 3, 2); a[1, 1] <- 1; a[2:3, 1] <- -colSums(ab$A)[2:3] / (N - Ar)
    cw <- .rxc_interval_lin(a, dd$lower, dd$upper)
    expect_true(cw[1] <= fit$table$dd_lower[1] + 1e-12 && cw[2] >= fit$table$dd_upper[1] - 1e-12)
    nwide <- nwide + (cw[2] - cw[1] > fit$table$dd_upper[1] - fit$table$dd_lower[1] + 1e-9)
  }
  expect_gt(nwide, 0)
})

# ---- worked examples --------------------------------------------------------------------
d2 <- data.frame(precinct = c("P1", "P2"), vet = c(20, 60), other = c(80, 40), yes = c(30, 50), no = c(70, 50))

test_that("E1: realised ordering raises the floors to 6 and 30; witnesses are whole tables", {
  fit <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d2, unit = "precinct",
                  targets = list(rate("vet", "yes")),
                  restrictions = list(order_unit = list(row = "vet", col = "yes")))
  expect_equal(c(fit$table$dd_lower, fit$table$dd_upper), c(0.125, 0.875))
  expect_equal(c(fit$table$lower, fit$table$upper), c(0.45, 0.875))
  expect_equal(fit$table$integer_status, "integer-exact (TU class)")
  w <- witness(fit)
  expect_equal(c(w[[1]]$table["P1", "vet", "yes"], w[[1]]$table["P2", "vet", "yes"]), c(6, 30))
  expect_equal(c(w[[2]]$table["P1", "vet", "yes"], w[[2]]$table["P2", "vet", "yes"]), c(20, 50))
  expect_true(wit_ok(fit))
  expect_equal(fit$table$binding_lower, "order_unit[1]")
})

test_that("E2: a published total and where= targets", {
  fit <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d2, unit = "precinct",
                  targets = list(rate("vet", "yes", where = precinct == "P2"), rate("vet", "yes")),
                  restrictions = list(total = list(units = c("P1", "P2"), row = "vet", col = "yes", value = 45)))
  expect_equal(c(fit$table$lower[1], fit$table$upper[1]) * 60, c(25, 45))
  expect_equal(c(fit$table$lower[2], fit$table$upper[2]), c(0.5625, 0.5625))
  expect_equal(fit$restrictions$role, "identity")
})

test_that("N17: odds-ratio dial [1, 4] in P2 gives the veteran rate [0.500, 0.633]", {
  fit <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d2, unit = "precinct",
                  targets = list(rate("vet", "yes", where = precinct == "P2")),
                  restrictions = list(odds = list(psi = c(1, 4), units = "P2")))
  # NS N17's [0.500, 0.633] is the continuous range; with whole people the upper end is 37/60:
  # the odds ratio at 38 veteran yes votes is 38 * 28 / (22 * 12) = 4.03 > 4
  expect_equal(fit$table$lower, 0.5)
  expect_equal(round(fit$table$lp_upper, 3), 0.633)
  expect_equal(fit$table$upper, 37 / 60)
  expect_true(wit_ok(fit))
  expect_lte(37 * (100 - 60 - 50 + 37) / ((60 - 37) * (50 - 37)), 4)
  expect_gt(38 * (100 - 60 - 50 + 38) / ((60 - 38) * (50 - 38)), 4)
  expect_output(print(fit), "untestable")
})

# ---- T5, S5, RC11: round the rows, not the endpoints ------------------------------------
A5 <- rbind(c(1, 4), c(1, 4)); B5 <- rbind(c(1, 4), c(1, 4))

test_that("T5: rounding the LP endpoint is not enough; the rows are rounded", {
  fit <- fitc(A5, B5, targets = list(rate(1, 1)), restrictions = list(order_unit = list(row = 1, col = 1)))
  expect_equal(fit$table$lp_lower * 2, 0.4)          # continuous lower total 0.4; its ceiling 1 is wrong
  expect_equal(fit$table$lower * 2, 2)               # integer answer 2
  expect_equal(fit$table$integer_status, "integer-exact (TU class)")
  f2 <- fitc(A5, B5, targets = list(rate(1, 1)), restrictions = list(order_unit = list(row = 1, col = 1)), integer = "lp")
  expect_equal(f2$table$lower * 2, 0.4)
  expect_equal(f2$table$integer_status, "continuous only")
  f3 <- fitc(A5, B5, targets = list(rate(1, 1)), restrictions = list(order_unit = list(row = 1, col = 1)), integer = "exact")
  expect_equal(f3$table$lower * 2, 2)
  expect_equal(f3$table$integer_status, "integer-exact (MILP)")
})

test_that("S5: with fractional boxes the integer values are not the lattice in the continuous range", {
  fit <- fitc(A5, B5, targets = list(rate(1, 1)),
              restrictions = list(order_unit = list(row = 1, col = 1), sign = list(type = "dB", row = 1, col = 1)))
  expect_equal(c(fit$table$lp_lower, fit$table$lp_upper) * 2, c(0.4, 2))
  expect_equal(c(fit$table$lower, fit$table$upper) * 2, c(2, 2))      # integer values {2}; 1 is not attained
  expect_true(wit_ok(fit))
})

test_that("RC11: fractional totals also shrink the integer set; every integer of [2, 4] is attained", {
  A <- do.call(rbind, rep(list(c(1, 4)), 4)); B <- A
  rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, lower = 0.2),
                          list(units = c("U3", "U4"), row = 1, col = 1, lower = 0.2)),
             sign = list(type = "dB", row = 1, col = 1))
  fit <- fitc(A, B, targets = list(rate(1, 1)), restrictions = rs)
  expect_equal(c(fit$table$lp_lower, fit$table$lp_upper) * 4, c(0.8, 4))
  expect_equal(c(fit$table$lower, fit$table$upper) * 4, c(2, 4))
  expect_equal(fit$table$tu_class, "(i)")
  att <- vapply(0:4, function(k) {
    f <- fitc(A, B, targets = list(rate(1, 1)), binding = FALSE,
              restrictions = c(rs, list(linear = list(cells = data.frame(unit = paste0("U", 1:4), row = 1, col = 1, coef = 1),
                                                  lower = k, upper = k))))
    isTRUE(f$table$integer_feasible) }, TRUE)
  expect_equal(att, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  # brute force
  ok <- function(t) { y <- vapply(t, function(m) m[1, 1], 0); y[1] + y[2] >= 0.2 && y[3] + y[4] >= 0.2 && sum(y) >= 0.8 }
  expect_equal(brute(A, B, function(t) sum(vapply(t, function(m) m[1, 1], 0)), ok), c(2, 4))
})

# ---- brute force per TU class, mixtures and MILP rows -----------------------------------
check_brute <- function(tabs, rs, ok, f_tg, target, expect_class = NULL) {
  ab <- AB(tabs)
  if (any(colSums(ab$A) == 0)) return(invisible(NULL))
  fit <- fitc(ab$A, ab$B, targets = list(target), restrictions = rs, binding = FALSE)
  br <- brute(ab$A, ab$B, f_tg, ok)
  if (is.na(br[1])) expect_false(isTRUE(fit$table$integer_feasible)) else {
    expect_equal(c(fit$table$lower, fit$table$upper), br, tolerance = 1e-12)
    expect_true(wit_ok(fit))
  }
  if (!is.null(expect_class) && isTRUE(fit$table$continuous_feasible)) {
    if (is.na(expect_class)) expect_true(is.na(fit$table$tu_class) || fit$table$integer_status != "integer-exact (TU class)")
    else expect_equal(fit$table$tu_class, expect_class)
  }
  invisible(fit)
}
cnt <- function(t, i, r, c) t[[i]][r, c]

test_that("brute force, class (i): one cell, two laminar (crossing) families, fractional box, known cell", {
  set.seed(21)
  for (k in 1:12) {
    tabs <- rand_tabs(3, 2, 2, 2); y <- vapply(tabs, function(m) m[1, 1], 0)
    t12 <- y[1] + y[2] + sample(c(-0.5, 0.5), 1); t23 <- y[2] + y[3] - 0.3
    rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, lower = t12),
                            list(units = c("U2", "U3"), row = 1, col = 1, upper = t23 + 1)),
               box = list(row = 1, col = 1, units = "U3", lower = 0.4),
               known = list(unit = "U1", row = 2, col = 2, value = tabs[[1]][2, 2]))
    ok <- function(t) { yy <- vapply(t, function(m) m[1, 1], 0)
      yy[1] + yy[2] >= t12 && yy[2] + yy[3] <= t23 + 1 && yy[3] >= 0.4 && t[[1]][2, 2] == tabs[[1]][2, 2] }
    ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    check_brute(tabs, rs, ok, function(t) sum(vapply(t, function(m) m[1, 1], 0)) / sum(ab$A[, 1]) -
                  sum(vapply(t, function(m) m[2, 1], 0)) / sum(ab$A[, 2]), gap(1, 2, 1), "(i)")
  }
})

test_that("brute force, class (ii): totals on cells of one column, one laminar family each", {
  set.seed(22)
  for (k in 1:10) {
    tabs <- rand_tabs(3, 3, 2, 2)
    v1 <- tabs[[1]][1, 1] + tabs[[2]][1, 1] + 0.5; v2 <- tabs[[2]][2, 1] + tabs[[3]][2, 1] - 0.5
    rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, upper = v1),
                            list(units = c("U2", "U3"), row = 2, col = 1, lower = v2)))
    ok <- function(t) t[[1]][1, 1] + t[[2]][1, 1] <= v1 && t[[2]][2, 1] + t[[3]][2, 1] >= v2
    ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    check_brute(tabs, rs, ok, function(t) sum(vapply(t, function(m) m[1, 1] / sum(ab$A[, 1]) + m[3, 1] / sum(ab$A[, 3]), 0)),
                combo(rate(1, 1), rate(3, 1), weights = c(1, 1)), "(ii)")
  }
})

test_that("brute force, class (iii): R = 2, totals in several columns (RC1 design)", {
  set.seed(23)
  for (k in 1:10) {
    tabs <- rand_tabs(3, 2, 2, 2)
    s1 <- sum(vapply(tabs, function(m) m[1, 1], 0)) - 0.4; s2 <- tabs[[1]][2, 1] + tabs[[2]][2, 1] + 0.3
    s3 <- tabs[[2]][1, 2] + tabs[[3]][1, 2] - 0.2
    rs <- list(total = list(list(units = c("U1", "U2", "U3"), row = 1, col = 1, lower = s1),
                            list(units = c("U1", "U2"), row = 2, col = 1, upper = s2),
                            list(units = c("U2", "U3"), row = 1, col = 2, lower = s3)))
    ok <- function(t) sum(vapply(t, function(m) m[1, 1], 0)) >= s1 && t[[1]][2, 1] + t[[2]][2, 1] <= s2 &&
      t[[2]][1, 2] + t[[3]][1, 2] >= s3
    ab <- AB(tabs); if (any(colSums(ab$A) == 0) || ab$A[1, 1] == 0 || ab$A[3, 2] == 0) next
    check_brute(tabs, rs, ok, function(t) t[[1]][1, 2] - t[[3]][2, 1], combo(rate(1, 2, where = "U1"), rate(2, 1, where = "U3"),
                weights = c(sum(ab$A[1, 1]), -sum(ab$A[3, 2]))), "(iii)")
  }
})

test_that("C1: the four row types outside the classes are routed to the MILP and are exact", {
  set.seed(24)
  for (k in 1:6) {    # pairwise unit ordering, R = 3
    tabs <- rand_tabs(2, 3, 2, 3); ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    rs <- list(order_unit = list(row = 1, over = 2, col = 1))
    ok <- function(t) all(vapply(seq_along(t), function(i) ab$A[i, 2] * t[[i]][1, 1] - ab$A[i, 1] * t[[i]][2, 1] >= 0, TRUE))
    f <- check_brute(tabs, rs, ok, function(t) sum(vapply(t, function(m) m[1, 1], 0)) / sum(ab$A[, 1]), rate(1, 1), NA)
    expect_equal(f$table$integer_status, "integer-exact (MILP)")
  }
  for (k in 1:6) {    # pairwise realised delta_B(r, s), R = 3 (checked against 2c Definition 1 with realised rates)
    tabs <- rand_tabs(3, 3, 2, 3); ab <- AB(tabs); if (any(ab$A == 0)) next
    X <- ab$A / rowSums(ab$A); st <- list(X = X, N = rowSums(ab$A), R = 3)
    W <- .rxc_assoc_W(st, "dB", "pair", 1, 2)
    ok <- function(t) sum(vapply(seq_along(t), function(i) sum(W[i, ] * t[[i]][, 1] / ab$A[i, ]), 0)) >= -1e-12
    f <- check_brute(tabs, list(sign = list(type = "dB", version = "pair", row = 1, other = 2, col = 1)), ok,
                     function(t) sum(vapply(t, function(m) m[2, 1], 0)) / sum(ab$A[, 2]), rate(2, 1), NA)
    expect_equal(f$table$integer_status, "integer-exact (MILP)")
  }
  for (k in 1:5) {    # totals on cells in different rows and columns, R = C = 3
    tabs <- rand_tabs(2, 3, 3, 2); ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    v1 <- tabs[[1]][1, 1] + tabs[[2]][1, 1]; v2 <- tabs[[1]][2, 2] + tabs[[2]][2, 2]
    rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, value = v1),
                            list(units = c("U1", "U2"), row = 2, col = 2, value = v2)))
    ok <- function(t) t[[1]][1, 1] + t[[2]][1, 1] == v1 && t[[1]][2, 2] + t[[2]][2, 2] == v2
    f <- check_brute(tabs, rs, ok, function(t) t[[1]][3, 3] / ab$A[1, 3], rate(3, 3, where = "U1"))
    if (ab$A[1, 3] > 0) expect_true(is.na(f$table$tu_class) && f$table$integer_status == "integer-exact (MILP)")
  }
  # three pairwise-crossing totals: LP feasible, no integer table (T4 / C1)
  A <- matrix(1, 3, 2); B <- matrix(1, 3, 2)
  rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, value = 1),
                          list(units = c("U2", "U3"), row = 1, col = 1, value = 1),
                          list(units = c("U1", "U3"), row = 1, col = 1, value = 1)))
  fit <- fitc(A, B, targets = list(rate(1, 1)), restrictions = rs)
  expect_true(fit$table$continuous_feasible)
  expect_false(fit$table$integer_feasible)
  expect_match(fit$table$integer_status, "no whole-number table satisfies these rows")
  expect_output(print(fit), "integer-infeasible")
})

test_that("RC4: a mixture of classes is sent to the MILP (whole row set, not row by row)", {
  set.seed(25)
  for (k in 1:8) {
    tabs <- rand_tabs(4, 3, 2, 2); ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    y <- vapply(tabs, function(m) m[1, 1], 0); z <- vapply(tabs, function(m) m[2, 1], 0)
    rs <- list(total = list(list(units = c("U1", "U2"), row = 1, col = 1, value = y[1] + y[2]),
                            list(units = c("U2", "U3"), row = 1, col = 1, value = y[2] + y[3]),
                            list(units = c("U3", "U4"), row = 2, col = 1, value = z[3] + z[4])))
    ok <- function(t) { yy <- vapply(t, function(m) m[1, 1], 0); zz <- vapply(t, function(m) m[2, 1], 0)
      yy[1] + yy[2] == y[1] + y[2] && yy[2] + yy[3] == y[2] + y[3] && zz[3] + zz[4] == z[3] + z[4] }
    f <- check_brute(tabs, rs, ok, function(t) sum(vapply(t, function(m) m[1, 2] - m[3, 1], 0)),
                     combo(rate(1, 2), rate(3, 1), weights = c(sum(ab$A[, 1]), -sum(ab$A[, 3]))), NA)
    expect_equal(f$table$integer_status, "integer-exact (MILP)")
  }
  P <- .ec_problem(diag(2) + 1, diag(2) + 1, c("a", "b"), c("x", "y"), 1:2)
  expect_true(.ec_two_laminar(list(1:2, 2:3, 3:4)))          # a path: bipartite
  expect_false(.ec_two_laminar(list(1:2, 2:3, c(1, 3))))     # a triangle of crossings
})

test_that("N10c: LP feasible, rows rounded, still no integer table: reported as integer-infeasible", {
  Ns <- c(13, 10, 14); Vs <- c(2, 8, 9); Ys <- c(13, 8, 4)
  A <- cbind(Vs, Ns - Vs); B <- cbind(Ys, Ns - Ys)
  rows <- lapply(list(c(1, 1, 1), Ns - Vs), function(w) {
    cf <- w * (1 / Vs + 1 / (Ns - Vs)); rhs <- sum(w * Ys / (Ns - Vs))
    # scaled to integers (any rational row can be; exact rounding needs exact coefficients)
    s <- 2 * 8 * 9 * 11 * 2 * 5
    list(cells = data.frame(unit = paste0("U", 1:3), row = 1, col = 1, coef = round(cf * s)), lower = rhs * s)
  })
  rs <- list(total = list(units = paste0("U", 1:3), row = 1, col = 1, value = 11), linear = rows)
  fit <- fitc(A, B, targets = list(rate(1, 1)), restrictions = rs)
  expect_true(fit$table$continuous_feasible)
  expect_false(fit$table$integer_feasible)
  expect_match(fit$table$integer_status, "integer-infeasible")
})

test_that("T1-T3: random TU-class instances, rounded LP = MILP", {
  set.seed(26)
  for (k in 1:10) {
    tabs <- rand_tabs(4, 2, 3, 4); ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    y <- vapply(tabs, function(m) m[1, 1], 0)
    rs <- list(total = list(list(units = paste0("U", 1:4), row = 1, col = 1, lower = sum(y) - 1.5),
                            list(units = c("U1", "U2"), row = 1, col = 1, upper = y[1] + y[2] + 0.7)),
               box = list(row = 1, col = 2, units = "U3", lower = 0.3), order_unit = list(row = 2, col = 3))
    tg <- combo(rate(1, 1), rate(2, 3), weights = c(2, -1))
    fa <- fitc(ab$A, ab$B, targets = list(tg), restrictions = rs, binding = FALSE)
    fe <- fitc(ab$A, ab$B, targets = list(tg), restrictions = rs, binding = FALSE, integer = "exact")
    if (!isTRUE(fa$table$integer_feasible)) { expect_false(isTRUE(fe$table$integer_feasible)); next }
    expect_equal(fa$table$integer_status, "integer-exact (TU class)")
    expect_equal(c(fa$table$lower, fa$table$upper), c(fe$table$lower, fe$table$upper), tolerance = 1e-12)
  }
})

# ---- 2c M2 / M3 closed forms in 2 x 2 (count coordinates) -------------------------------
test_that("2c M3 in counts: between sign lower end max(Y-, ceil Y^NM); within sign emptiness iff Y^ER < Y-", {
  set.seed(27)
  for (k in 1:25) {
    I <- sample(2:4, 1); Ns <- sample(4:12, I, TRUE); Vs <- vapply(Ns, function(n) sample(1:(n - 1), 1), 0)
    Ys <- vapply(Ns, function(n) sample(0:n, 1), 0)
    A <- cbind(Vs, Ns - Vs); B <- cbind(Ys, Ns - Ys)
    X <- Vs / Ns; if (var(X) == 0) next
    N <- sum(Ns); Ar <- sum(Vs); Yc <- sum(Ys); p <- Ns / N
    YNM <- sum(Vs * Ys / Ns); Ylo <- sum(pmax(0, Vs + Ys - Ns)); Yhi <- sum(pmin(Vs, Ys))
    mu <- sum(p * X); vx <- sum(p * X^2) - mu^2; ER <- (sum(p * X * Ys / Ns) - mu * sum(p * Ys / Ns)) / vx
    kap <- N / (Ar * (N - Ar)); YER <- (ER + Yc / (N - Ar)) / kap
    fB <- fitc(A, B, targets = list(rate(1, 1)), restrictions = list(sign = list(type = "dB", row = 1, col = 1)))
    expect_equal(fB$table$lp_lower * Ar, max(Ylo, YNM), tolerance = 1e-9)
    expect_equal(fB$table$lower * Ar, max(Ylo, ceiling(YNM - 1e-9)), tolerance = 1e-9)
    expect_equal(fB$table$upper * Ar, Yhi, tolerance = 1e-9)
    fW <- fitc(A, B, targets = list(ovr_gap(1, 1)), restrictions = list(sign = list(type = "dW", row = 1, col = 1)))
    if (YER < Ylo - 1e-9) expect_false(fW$table$continuous_feasible) else if (YER > Ylo + 1e-9) {
      expect_equal(fW$table$lp_upper, kap * min(Yhi, YER) - Yc / (N - Ar), tolerance = 1e-9)
      expect_equal(fW$table$upper, kap * min(Yhi, floor(YER + 1e-9)) - Yc / (N - Ar), tolerance = 1e-9)
    }
  }
})

test_that("2c Remark M2.1 examples: continuous CR image {0} but no integer table; between-sign lower gap 2/15", {
  A <- cbind(c(1, 2, 3), c(3, 2, 1)); B <- cbind(c(0, 1, 0), c(4, 3, 4))
  f <- fitc(A, B, targets = list(ovr_gap(1, 1)), restrictions = list(cr = list(row = 1, col = 1)))
  expect_true(f$table$continuous_feasible)
  expect_equal(c(f$table$lp_lower, f$table$lp_upper), c(0, 0), tolerance = 1e-12)
  expect_false(f$table$integer_feasible)
  A <- cbind(c(1, 2), c(3, 2)); B <- cbind(c(1, 1), c(3, 3))
  f <- fitc(A, B, targets = list(ovr_gap(1, 1)), restrictions = list(sign = list(type = "dB", row = 1, col = 1)))
  expect_equal(f$table$lower, 2 / 15, tolerance = 1e-12)
  expect_equal(f$table$lp_lower, 0, tolerance = 1e-12)
  expect_true(wit_ok(f))
})

test_that("S2: W1 sign rows are totals (class (i)); the generic W1 form gives the same range", {
  set.seed(28)
  for (k in 1:8) {
    tabs <- rand_tabs(3, 3, 2, 4); ab <- AB(tabs); if (any(ab$A == 0)) next
    f <- fitc(ab$A, ab$B, targets = list(rate(1, 1)), binding = FALSE,
              restrictions = list(sign = list(type = "dB", row = 1, col = 1), box = list(row = 2, col = 1, units = "U1", lower = 0.5)))
    if (!isTRUE(f$table$integer_feasible)) next
    expect_equal(f$table$tu_class, "(i)")
    P <- f$P
    g <- .ec_assoc_gen(P, 1:3, "dB", "W1", 1, 1)
    th <- .ec_ovr_thresholds(P, 1:3, 1, 1)
    # g'n = (Y - Y^NM) / N for every table: compare on the two witnesses
    for (w in witness(f)) {
      x <- as.vector(aperm(w$table, c(3, 2, 1)))
      expect_equal(as.double(sum(g$coef * gmp::as.bigq(x[g$idx]))), (sum(x[.ec_vix(P, 1:3, 1, 1)]) - as.double(th$YNM)) / as.double(th$N), tolerance = 1e-12)
    }
  }
})

# ---- exact arithmetic (RC9) ---------------------------------------------------------------
test_that("RC9: thresholds that are integers in rational arithmetic are rounded to themselves", {
  set.seed(29); n_int <- 0
  for (k in 1:3000) {
    I <- sample(2:4, 1); Ns <- sample(3:9, I, TRUE); Vs <- vapply(Ns, function(n) sample(1:(n - 1), 1), 0)
    Ys <- vapply(Ns, function(n) sample(0:n, 1), 0)
    if (var(Vs / Ns) == 0) next
    P <- .ec_problem(cbind(Vs, Ns - Vs), cbind(Ys, Ns - Ys), c("a", "b"), c("y", "n"), seq_len(I))
    th <- .ec_ovr_thresholds(P, seq_len(I), 1, 1)
    if (is.null(th$ER) || !gmp::is.whole(th$YER_hi)) next
    n_int <- n_int + 1
    rw <- .ec_round_row(.ec_row(P, "s", "total", .ec_vix(P, seq_len(I), 1, 1), rep(1, I), NULL, th$YER_hi))
    expect_identical(rw$hi_rd, as.double(th$YER_hi))
  }
  expect_gt(n_int, 5)
  # margins count as integers only if exactly integer in the input
  P <- .ec_problem(cbind(c(1, 2) + 1e-13, c(3, 2)), cbind(c(2, 2) + 1e-13, c(2, 2)), c("a", "b"), c("y", "n"), 1:2)
  expect_false(P$integer)
})

# ---- restrictions only narrow; witnesses -------------------------------------------------
test_that("adding a row never widens a range, and every witness satisfies every row", {
  set.seed(30)
  for (k in 1:8) {
    tabs <- rand_tabs(3, 2, 2, 5); ab <- AB(tabs); if (any(colSums(ab$A) == 0)) next
    rs1 <- list(total = list(units = c("U1", "U2"), row = 1, col = 1, value = tabs[[1]][1, 1] + tabs[[2]][1, 1]))
    rs2 <- c(rs1, list(order_unit = list(row = 1, col = 1, sign = sample(c(-1, 1), 1))))
    rs3 <- c(rs2, list(cr = list(row = 1, col = 1)))
    ff <- lapply(list(NULL, rs1, rs2, rs3), function(r) fitc(ab$A, ab$B, targets = list(rate(1, 1), ovr_gap(1, 2)), restrictions = r))
    for (j in 2:4) for (t in 1:2) {
      a <- ff[[j - 1]]$table[t, ]; b <- ff[[j]]$table[t, ]
      lo_a <- if (j == 2) a$dd_lower else a$lower; hi_a <- if (j == 2) a$dd_upper else a$upper
      if (!is.na(b$lower)) expect_true(b$lower >= lo_a - 1e-12 && b$upper <= hi_a + 1e-12)
    }
    for (f in ff) expect_true(wit_ok(f))
  }
})

# ---- non-integer margins (R19(e)) ----------------------------------------------------------
test_that("non-integer margins: integer check off, the note printed, integer = 'exact' refused", {
  d <- data.frame(vet = c(59.4, 20.2), other = c(40.6, 79.8), yes = c(50, 30), no = c(50, 70))
  fit <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, targets = list(rate("vet", "yes")),
                  restrictions = list(order_unit = list(row = "vet", col = "yes")))
  expect_equal(fit$table$integer_status, "not applicable: non-integer margins")
  expect_output(print(fit), "exact only for the stated margins", fixed = TRUE)
  expect_true(wit_ok(fit))
  expect_error(eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, integer = "exact"), "integer margins")
  # the lower end is attained by a real table with exactly those margins (E5)
  w <- witness(fit)[[1]]
  expect_equal(sum(w$table[, "vet", ]), sum(d$vet))
})

test_that("units with N_i = 0 are dropped; degenerate units in orderings and dials are skipped", {
  d <- data.frame(u = c("a", "b", "z"), v = c(2, 0, 0), o = c(3, 4, 0), y = c(1, 2, 0), n = c(4, 2, 0))
  fit <- eiCounts(cbind(y, n) ~ cbind(v, o), data = d, unit = "u", targets = list(rate("v", "y")),
                  restrictions = list(order_unit = list(row = "v", col = "y"), odds = list(psi = c(1, 2))))
  expect_equal(fit$P$p, 2)
  expect_equal(fit$P$n_dropped, 1)
  expect_true(wit_ok(fit))
  expect_output(print(fit), "1 with N_i = 0 dropped")
})

test_that("CR: branches are reported and the union may be two pieces", {
  A <- cbind(c(1, 2, 3), c(3, 2, 1)); B <- cbind(c(1, 2, 2), c(3, 2, 2))
  f <- fitc(A, B, targets = list(ovr_gap(1, 1)), restrictions = list(cr = list(row = 1, col = 1, slack = 0.3)))
  expect_output(print(f), "CR branches")
  expect_true(wit_ok(f))
  expect_true(f$table$n_pieces %in% 1:2)
})

# ---- object tag (C9) ---------------------------------------------------------------------
test_that("C9: object = 'model' and object = c('actual', 'model') are refused in the count layer", {
  expect_error(eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d2,
                        restrictions = list(sign = list(type = "dB", row = "vet", col = "yes", object = "model"))),
               "eiBounds")
  expect_error(eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d2,
                        restrictions = list(sign = list(type = "dB", row = "vet", col = "yes", object = c("actual", "model")))),
               "one label in both layers")
})

# ---- cross-layer block (R27, C8) ---------------------------------------------------------
sim_counts <- function(n = 80, seed = 5) {
  set.seed(seed)
  N <- sample(50:150, n, TRUE); x <- runif(n, 0.1, 0.9); V <- round(N * x)
  b <- 0.3 + 0.3 * V / N; yv <- rbinom(n, V, b); yo <- rbinom(n, N - V, 0.4)
  data.frame(u = paste0("p", 1:n), vet = V, other = N - V, yes = yv + yo, no = N - yv - yo)
}

test_that("C8: identity rows only -> coverage wording; belief rows -> the belief wording, never 'given the count rows'", {
  d <- sim_counts(800)
  m <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", seed = 1, ndraw = 30, nz_draw = 500,
                tests = FALSE, target = "model")
  f1 <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", model = m,
                 targets = list(rate("vet", "yes"), gap("vet", "other", "yes"), ovr_gap("vet", "yes"), share("vet", "yes")),
                 restrictions = NULL)
  o1 <- paste(capture.output(print(f1)), collapse = "\n")
  expect_true(grepl(.ec_txt("cov_id"), o1, fixed = TRUE))
  expect_false(grepl("given the count rows", o1))
  expect_match(f1$cross$table$note[4], "count layer only")
  f2 <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", model = m,
                 targets = list(rate("vet", "yes")),
                 restrictions = list(order_unit = list(row = "vet", col = "yes")))
  o2 <- paste(capture.output(print(f2)), collapse = "\n")
  expect_true(grepl(.ec_txt("cov_bel"), o2, fixed = TRUE))
  expect_false(grepl("given the count rows", o2))
  # a false belief row makes the sets disjoint: the belief wording of the rejection
  f3 <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", model = m,
                 targets = list(rate("vet", "yes")),
                 restrictions = list(total = list(row = "vet", col = "yes", role = "belief",
                                                  upper = sum(pmax(0, d$yes - d$other)) + 2)))
  o3 <- paste(capture.output(print(f3)), collapse = "\n")
  expect_true(f3$cross$table$disjoint)
  expect_true(grepl(.ec_txt("rej_bel"), o3, fixed = TRUE))
  # identity rows only and disjoint: a published total far from the truth
  tot <- sum(d$vet)
  f4 <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "u", model = m,
                 targets = list(rate("vet", "yes")),
                 restrictions = list(total = list(row = "vet", col = "yes",
                                                  value = sum(pmax(0, d$yes - d$other)) + 2)))
  expect_true(f4$table$integer_feasible)
  expect_true(f4$cross$table$disjoint)
  o4 <- paste(capture.output(print(f4)), collapse = "\n")
  expect_true(grepl(.ec_txt("rej_id"), o4, fixed = TRUE))
  expect_false(grepl(.ec_txt("cov_bel"), o4, fixed = TRUE))
  # the intersection never lies outside the count set
  for (f in list(f1, f2, f3, f4)) {
    tb <- f$cross$table; j <- !is.na(tb$lower)
    expect_true(all(tb$lower[j] >= tb$count_lower[j] - 1e-12 & tb$upper[j] <= tb$count_upper[j] + 1e-12))
  }
  # Lemma 4.1: S_phi for a single cell is P3.5(a)'s S1; for a one-column gap rho = w_r + w_s
  w <- d$vet / sum(d$vet); w2 <- d$other / sum(d$other)
  expect_equal(f1$cross$table$S_phi[1], 0.5 * sqrt(sum(w^2)), tolerance = 1e-12)
  expect_equal(f1$cross$table$S_phi[2], 0.5 * sqrt(sum((w + w2)^2)), tolerance = 1e-12)
  # refusals: a model estimate (ci = "none") and trimmed fits
  m0 <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", ci = "none", tests = FALSE, target = "model")
  expect_error(eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, model = m0), "not an actual-target")
  expect_error(eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, model = list(a = 1)), "eiBounds")
})

# ---- eiBounds embedding (D1/D2, C1) ------------------------------------------------------
test_that("eiBounds prints the actual-rate block first (DD first); target = 'model' restores the old order", {
  d <- sim_counts(60, 7)
  fa <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", seed = 1, ndraw = 30, nz_draw = 500, tests = FALSE)
  fm <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", seed = 1, ndraw = 30, nz_draw = 500, tests = FALSE,
                 target = "model")
  oa <- capture.output(print(fa)); om <- capture.output(print(fm))
  expect_match(oa[1], "Actual \\(realised\\) rates")
  expect_true(grep("margins only = Duncan-Davis", oa)[1] < grep("R x C ecological inference bounds", oa)[1])
  expect_equal(om[1], "R x C ecological inference bounds")
  expect_false(any(grepl("count layer", om)))
  # the model results do not depend on the target argument
  cols <- setdiff(names(fm$table), c("dd_exact_lower", "dd_exact_upper"))
  expect_equal(fa$table[, cols], fm$table[, cols])
  # exact DD = cellwise DD for single cells
  expect_equal(fa$table$dd_exact_lower, fa$table$dd_lower, tolerance = 1e-12)
  expect_equal(fa$table$dd_exact_upper, fa$table$dd_upper, tolerance = 1e-12)
  expect_true(wit_ok(fa$counts))
  expect_s3_class(witness(fa), "ecWitness")
})

test_that("eiBounds: assoc with object = 'actual' goes to the count layer; both objects under one label are refused", {
  d <- sim_counts(60, 8)
  f <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", seed = 1, ndraw = 30, nz_draw = 500, tests = FALSE,
                restrictions = list(assoc = list(list(type = "dB", row = "vet", col = "yes", object = "actual"))))
  expect_null(f$setup$extra$assoc)
  expect_equal(f$counts$restrictions$id, "sign[1]")
  o <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(o, "EXPERIMENTAL")
  expect_true(grepl(.ec_txt("cov_bel"), o, fixed = TRUE))
  expect_error(eiBounds(cbind(yes, no) ~ cbind(vet, other), data = d, model = "M0", tests = FALSE,
                        restrictions = list(assoc = list(list(type = "dB", row = "vet", col = "yes", object = c("actual", "model"))))),
               "one label in both layers")
})
