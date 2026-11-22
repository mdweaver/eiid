# Code review of eiCounts() / the eiBounds() count block (PHASE6a, M-COUNT, step 1).
# Independent checks; see docs/extension/reviews/PHASE6a_code_review.md.
#
#   Rscript phase6a_code_review_checks.R oracle INSTANCES.json RESULTS.json
#       run eiCounts() on the oracle instances written by phase6a_code_review_oracle.py
#   Rscript phase6a_code_review_checks.R all
#       the named checks K1-K.. below (prints PASS / FAIL / FINDING per check)
suppressPackageStartupMessages({ library(eiPartialID); library(jsonlite) })
ns <- asNamespace("eiPartialID")
args <- commandArgs(TRUE)

# ---------------------------------------------------------------------------
# oracle driver
# ---------------------------------------------------------------------------
to_spec <- function(s) {
  t <- s$type; o <- s; o$type <- NULL
  if (t == "sign") { o$type <- if (is.null(s$stype)) "dB" else s$stype; o$stype <- NULL }
  if (t == "linear") {
    cl <- do.call(rbind, lapply(s$cells, function(z) as.numeric(unlist(z))))
    o$cells <- data.frame(unit = as.integer(cl[, 1]), row = as.integer(cl[, 2]), col = as.integer(cl[, 3]), coef = cl[, 4])
  }
  if (!is.null(o$units)) o$units <- as.integer(unlist(o$units))
  if (!is.null(o$psi)) o$psi <- as.numeric(unlist(o$psi))
  list(type = t, spec = o)
}
to_target <- function(t, rl, cl) {
  w <- if (is.null(t$units)) NULL else as.integer(unlist(t$units))
  switch(t$kind,
         rate = eiPartialID::rate(rl[unlist(t$row)], cl[unlist(t$col)], where = w),
         gap = { w2 <- if (is.null(t$units2)) w else as.integer(unlist(t$units2))
                 eiPartialID::gap(rl[t$row1], rl[t$row2], cl[t$col], where = w, where2 = w2) },
         ovr_gap = eiPartialID::ovr_gap(rl[t$row], cl[t$col], where = w),
         share = eiPartialID::share(rl[t$row], cl[t$col], where = w))
}
run_oracle <- function(fin, fout) {
  insts <- fromJSON(fin, simplifyVector = FALSE)
  res <- vector("list", length(insts))
  t0 <- proc.time()[3]
  for (k in seq_along(insts)) {
    I <- insts[[k]]
    A <- do.call(rbind, lapply(I$A, function(z) as.numeric(unlist(z))))
    B <- do.call(rbind, lapply(I$B, function(z) as.numeric(unlist(z))))
    rl <- paste0("r", seq_len(ncol(A))); cl <- paste0("c", seq_len(ncol(B)))
    d <- data.frame(A, B); names(d) <- c(rl, cl)
    restr <- list()
    for (s in I$specs) { z <- to_spec(s); restr[[z$type]] <- c(restr[[z$type]], list(z$spec)) }
    # labels rather than indices for rows / columns
    for (nm in names(restr)) restr[[nm]] <- lapply(restr[[nm]], function(s) {
      for (f in c("row", "over", "other")) if (!is.null(s[[f]])) s[[f]] <- rl[s[[f]]]
      if (!is.null(s$col)) s$col <- cl[s$col]
      if (!is.null(s$cells)) { s$cells$row <- rl[s$cells$row]; s$cells$col <- cl[s$cells$col] }
      s })
    tg <- lapply(I$targets, to_target, rl = rl, cl = cl)
    out <- list(id = I$id)
    fit <- tryCatch(eiCounts(rows = rl, cols = cl, data = d, targets = tg,
                             restrictions = if (length(restr)) restr else NULL, binding = FALSE),
                    error = function(e) e)
    if (inherits(fit, "error")) { out$error <- conditionMessage(fit); res[[k]] <- out; next }
    tb <- fit$table
    na2null <- function(v) lapply(v, function(z) if (is.na(z)) NULL else z)
    out$lower <- na2null(tb$lower); out$upper <- na2null(tb$upper)
    out$lp_lower <- na2null(tb$lp_lower); out$lp_upper <- na2null(tb$lp_upper)
    out$status <- tb$integer_status; out$class <- ifelse(is.na(tb$tu_class), "", tb$tu_class)
    out$witness <- lapply(seq_len(nrow(tb)), function(j) {
      ws <- witness(fit, target = j)
      lapply(c("lower", "upper"), function(e) {
        z <- Filter(function(q) q$end == e, ws)
        if (!length(z)) return(NULL)
        arr <- z[[1]]$table                    # units x rows x cols -> (i, r, c) row-major
        list(x = as.vector(aperm(arr, c(3, 2, 1))), margins_ok = z[[1]]$margins_ok, rows_ok = z[[1]]$rows_ok)
      })
    })
    res[[k]] <- out
  }
  cat(sprintf("eiCounts on %d instances: %.1f s\n", length(insts), proc.time()[3] - t0))
  writeLines(toJSON(res, auto_unbox = TRUE, null = "null", digits = NA), fout)
}

if (nzchar(Sys.getenv("PATCH"))) source(Sys.getenv("PATCH"))
if (length(args) && args[1] == "oracle") { run_oracle(args[2], args[3]); quit(save = "no") }

# ---------------------------------------------------------------------------
# named checks
# ---------------------------------------------------------------------------
RES <- list()
report <- function(id, ok, msg, finding = FALSE) {
  tag <- if (finding) "FINDING" else if (isTRUE(ok)) "PASS" else "FAIL"
  RES[[id]] <<- tag
  cat(sprintf("[%s] %s: %s\n", tag, id, msg))
}
mk <- function(A, B) {
  d <- as.data.frame(cbind(A, B)); names(d) <- c(paste0("r", seq_len(ncol(A))), paste0("c", seq_len(ncol(B))))
  d$u <- paste0("U", seq_len(nrow(A))); d
}
fc <- function(A, B, ...) eiCounts(rows = paste0("r", seq_len(ncol(A))), cols = paste0("c", seq_len(ncol(B))),
                                   data = mk(A, B), unit = "u", ...)
tpath <- function(f) file.path(Sys.getenv("PKG", "/home/user/eiid/package/eiPartialID"), "tests", "testthat", f)

if (length(args) && args[1] == "all") {

# K1 T5 ------------------------------------------------------------------------------------
A5 <- rbind(c(1, 4), c(1, 4))
f <- fc(A5, A5, targets = list(rate("r1", "c1")), restrictions = list(order_unit = list(row = "r1", col = "c1")))
report("K1 T5", abs(f$table$lp_lower * 2 - 0.4) < 1e-12 && f$table$lower * 2 == 2 && f$table$integer_status == "integer-exact (TU class)",
       sprintf("LP lower total %.1f, integer lower total %g (%s)", f$table$lp_lower * 2, f$table$lower * 2, f$table$integer_status))

# K2 S5 ------------------------------------------------------------------------------------
f <- fc(A5, A5, targets = list(rate("r1", "c1")),
        restrictions = list(order_unit = list(row = "r1", col = "c1"), sign = list(type = "dB", row = "r1", col = "c1")))
report("K2 S5", all(abs(c(f$table$lp_lower, f$table$lp_upper) * 2 - c(0.4, 2)) < 1e-12) && all(c(f$table$lower, f$table$upper) * 2 == 2),
       sprintf("continuous [%.1f, %.1f], integer [%g, %g]", f$table$lp_lower * 2, f$table$lp_upper * 2, f$table$lower * 2, f$table$upper * 2))

# K3 RC11 ----------------------------------------------------------------------------------
A4 <- do.call(rbind, rep(list(c(1, 4)), 4))
rs <- list(total = list(list(units = c("U1", "U2"), row = "r1", col = "c1", lower = 0.2),
                        list(units = c("U3", "U4"), row = "r1", col = "c1", lower = 0.2)),
           sign = list(type = "dB", row = "r1", col = "c1"))
f <- fc(A4, A4, targets = list(rate("r1", "c1")), restrictions = rs)
rs2 <- rs
att <- vapply(0:4, function(k) { rs2$total <- c(rs$total, list(list(units = paste0("U", 1:4), row = "r1", col = "c1", value = k)))
  isTRUE(fc(A4, A4, targets = list(rate("r1", "c1")), binding = FALSE, restrictions = rs2)$table$integer_feasible) }, TRUE)
report("K3 RC11", all(abs(c(f$table$lp_lower, f$table$lp_upper) * 4 - c(0.8, 4)) < 1e-12) && all(c(f$table$lower, f$table$upper) * 4 == c(2, 4)) &&
         identical(att, c(FALSE, FALSE, TRUE, TRUE, TRUE)),
       sprintf("continuous [%.1f, %g], integer [%g, %g], attained integers %s (class %s)", f$table$lp_lower * 4, f$table$lp_upper * 4,
               f$table$lower * 4, f$table$upper * 4, paste((0:4)[att], collapse = ","), f$table$tu_class))

# K4 N10c, with the rows given in their natural (non-integer) form ------------------------------
Ns <- c(13, 10, 14); Vs <- c(2, 8, 9); Ys <- c(13, 8, 4)
A <- cbind(Vs, Ns - Vs); B <- cbind(Ys, Ns - Ys)
lin <- lapply(list(c(1, 1, 1), Ns - Vs), function(w) list(cells = data.frame(unit = paste0("U", 1:3), row = "r1", col = "c1",
  coef = w * (1 / Vs + 1 / (Ns - Vs))), lower = sum(w * Ys / (Ns - Vs))))
f <- fc(A, B, targets = list(rate("r1", "c1")), restrictions = list(total = list(units = paste0("U", 1:3), row = "r1", col = "c1", value = 11), linear = lin))
report("K4 N10c", isTRUE(f$table$continuous_feasible) && identical(f$table$integer_feasible, FALSE) && grepl("integer-infeasible", f$table$integer_status),
       sprintf("continuous feasible %s, status '%s'", f$table$continuous_feasible, f$table$integer_status))

# K5 RC4 mixture -> MILP; crossing triple -> MILP; class (iii) -> TU ----------------------------
set.seed(4)
tabs <- lapply(1:3, function(i) matrix(sample(1:3, 9, TRUE), 3))
A <- do.call(rbind, lapply(tabs, rowSums)); B <- do.call(rbind, lapply(tabs, colSums))
f1 <- fc(A, B, targets = list(rate("r1", "c1")), binding = FALSE, restrictions = list(total = list(
  list(units = c("U1", "U2"), row = "r1", col = "c1", lower = 1.5), list(units = c("U2", "U3"), row = "r2", col = "c1", lower = 1.5),
  list(units = c("U1", "U2", "U3"), row = "r1", col = "c1", upper = 6.5), list(units = c("U2", "U3"), row = "r1", col = "c1", upper = 5.5))))
f2 <- fc(A, B, targets = list(rate("r1", "c1")), binding = FALSE, restrictions = list(total = list(
  list(units = c("U1", "U2"), row = "r1", col = "c1", lower = 1), list(units = c("U2", "U3"), row = "r1", col = "c1", lower = 1),
  list(units = c("U1", "U3"), row = "r1", col = "c1", lower = 1))))
report("K5 RC4 routing", f1$table$integer_status == "integer-exact (MILP)" && f2$table$integer_status == "integer-exact (MILP)",
       sprintf("two families on cell (1,1) + cell (2,1): %s; three crossing sets on one cell: %s", f1$table$integer_status, f2$table$integer_status))

# K6 exact rounding of thresholds (RC9) ---------------------------------------------------------
# random populations where Y^ER is an integer in exact arithmetic: the integer upper end of Y_rc
# under delta_W >= 0 must be min(Y^+, Y^ER), computed here independently with gmp
set.seed(9); ntry <- 0; nbad <- 0
while (ntry < 150) {
  I <- sample(2:4, 1); Nn <- sample(2:11, I, TRUE); V <- vapply(Nn, function(n) sample(1:(n - 1), 1), 0); Y <- vapply(Nn, function(n) sample(0:n, 1), 0)
  q <- gmp::as.bigq; N <- sum(Nn); p <- q(Nn) / N; X <- q(V) / q(Nn); Tt <- q(Y) / q(Nn)
  mu <- sum(p * X); vr <- sum(p * (X - mu)^2); if (vr == 0) next
  ER <- sum(p * (Tt - sum(p * Tt)) * (X - mu)) / vr
  Av <- sum(V); kap <- q(1, Av) + q(1, N - Av); YER <- (ER + q(sum(Y), N - Av)) / kap
  if (gmp::denominator(YER) != 1) next
  ntry <- ntry + 1
  A <- cbind(V, Nn - V); B <- cbind(Y, Nn - Y)
  Yp <- sum(pmin(V, Y)); Ym <- sum(pmax(0, V + Y - Nn))
  f <- fc(A, B, targets = list(rate("r1", "c1")), binding = FALSE, restrictions = list(sign = list(type = "dW", row = "r1", col = "c1")))
  exp_hi <- min(Yp, as.double(YER)); exp_ok <- as.double(YER) >= Ym
  got <- f$table$upper * Av
  if (exp_ok) { if (is.na(got) || abs(got - exp_hi) > 1e-9) nbad <- nbad + 1 } else if (isTRUE(f$table$integer_feasible)) nbad <- nbad + 1
}
report("K6 RC9 exact thresholds", nbad == 0, sprintf("%d populations with integer Y^ER: %d wrong integer upper ends", ntry, nbad))

# K7 linear rows: coefficients as doubles (deviation 2) ------------------------------------------
d <- data.frame(u = c("U1", "U2"), r1 = c(5, 5), r2 = c(5, 5), c1 = c(5, 5), c2 = c(5, 5))
f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u", targets = list(rate("r1", "c1", where = "U1")), binding = FALSE,
              restrictions = list(linear = list(cells = data.frame(unit = "U1", row = "r1", col = "c1", coef = 1/3), lower = 1)))
report("K7 linear row 1/3 n >= 1", f$table$lower * 5 == 3,
       sprintf("intended n >= 3; integer lower %g (LP %g), status %s: the double 1/3 is taken as an exact rational < 1/3",
               f$table$lower * 5, f$table$lp_lower * 5, f$table$integer_status), finding = f$table$lower * 5 != 3)

# K8 non-integer margins switch the integer check off ------------------------------------------
A <- rbind(c(59.4, 40.6), c(30, 70)); B <- rbind(c(50, 50), c(40, 60))
f <- fc(A, B, targets = list(rate("r1", "c1")), restrictions = list(order_unit = list(row = "r1", col = "c1")))
o <- capture.output(print(f))
e <- tryCatch(fc(A, B, targets = list(rate("r1", "c1")), integer = "exact"), error = function(e) "refused")
A2 <- rbind(c(3, 7), c(4, 6)) ; B2 <- rbind(c(5, 5), c(5, 5))
f2 <- fc(A2 + 1e-12, B2 + c(1e-12, 0), targets = list(rate("r1", "c1")))
f3 <- eiCounts(x = A2[, 1] / 10, t = B2[, 1] / 10, n = c(10, 10), targets = list(rate(1, 1)))
report("K8b x/t vectors: row and column labels", identical(f3$P$rows, c("X", "notX")),
       sprintf("rows %s, cols %s (documented X/notX, T/notT; cbind(x, 1 - x) names them 'x' and '')",
               paste0("'", f3$P$rows, "'", collapse = ","), paste0("'", f3$P$cols, "'", collapse = ",")), finding = !identical(f3$P$rows, c("X", "notX")))
report("K8 non-integer margins", f$table$integer_status == "not applicable: non-integer margins" && any(grepl("exact only for the stated margins", o)) &&
         identical(e, "refused") && f2$table$integer_status == "not applicable: non-integer margins" && f3$table$integer_status == "not applicable: non-integer margins",
       sprintf("59.4: '%s'; note printed %s; exact refused %s; 3+1e-12: '%s'; shares input: '%s'", f$table$integer_status,
               any(grepl("exact only for the stated margins", o)), identical(e, "refused"), f2$table$integer_status, f3$table$integer_status))

# K9 DD to 1e-12 on the reference data (single cells, one-column contrasts) --------------------
mx <- 0; mxc <- 0; nt <- 0; narrower <- 0
for (fn in c("rxc_example_3x3.csv", "msign_example_3x2.csv", "msign_gap_3x2.csv")) {
  dd <- read.csv(tpath(fn)); X <- as.matrix(dd[, grep("^X", names(dd))]); Tm <- as.matrix(dd[, grep("^T", names(dd))])
  R <- ncol(X); C <- ncol(Tm); tg <- list()
  for (r in 1:R) for (c in 1:C) { a <- matrix(0, R, C); a[r, c] <- 1; tg[[sprintf("cell%d%d", r, c)]] <- a }
  for (c in 1:C) for (r in 1:R) for (s in 1:R) if (r < s) { a <- matrix(0, R, C); a[r, c] <- 1; a[s, c] <- -1; tg[[sprintf("g%d%d_%d", r, s, c)]] <- a }
  a <- matrix(0, R, C); a[1, ] <- 1; tg[["rowsum"]] <- a                     # combination: exact LP narrower or equal
  fit <- eiBounds(x = X, t = Tm, n = dd$N, targets = tg, model = "M0", ci = "none")
  tb <- fit$table; single <- !(tb$target %in% "rowsum")
  mx <- max(mx, abs(tb$dd_lower - tb$dd_exact_lower)[single], abs(tb$dd_upper - tb$dd_exact_upper)[single]); nt <- nt + sum(single)
  narrower <- narrower + sum(tb$dd_exact_lower[!single] >= tb$dd_lower[!single] - 1e-12 & tb$dd_exact_upper[!single] <= tb$dd_upper[!single] + 1e-12)
}
report("K9 DD = exact LP for cells and one-column contrasts", mx <= 1e-12 && narrower == 3,
       sprintf("%d targets on 3 reference datasets: max |cellwise DD - exact LP| = %.1e; combination targets inside cellwise DD: %d/3", nt, mx, narrower))

} # end 'all' (part 1)

if (length(args) && args[1] %in% c("all", "r27")) {
# ---------------------------------------------------------------------------
# R27 code: alpha split, S_phi / t_E, representation, statuses, wording, refusals
# ---------------------------------------------------------------------------
# simulated 2 x 2 population (pooled LCE true: rates with a common mean); used for K11 and K16
sim_pop <- function(p, seed, m = c(0.7, 0.3), sd = 0.12) {
  set.seed(seed)
  Nn <- sample(150:600, p, TRUE); x <- stats::rbeta(p, 2, 2.5)
  a1 <- stats::rbinom(p, Nn, x); a2 <- Nn - a1
  b1 <- pmin(pmax(stats::rnorm(p, m[1], sd), 0), 1); b2 <- pmin(pmax(stats::rnorm(p, m[2], sd), 0), 1)
  y1 <- stats::rbinom(p, a1, b1); y2 <- stats::rbinom(p, a2, b2)
  data.frame(u = paste0("P", 1:p), g1 = a1, g2 = a2, yes = y1 + y2, no = Nn - y1 - y2, n11 = y1, n21 = y2)
}
dd <- sim_pop(1500, 1)
dd <- dd[dd$g1 > 0 & dd$g2 > 0, ]
fitP <- eiBounds(cbind(yes, no) ~ cbind(g1, g2), data = dd, model = "M0", fpc = "popoviciu", seed = 11,
                 targets = list(c11 = c("g1", "yes"), c21 = c("g2", "yes")))
cr <- fitP$counts$cross$table; tb <- fitP$table
cl_lo <- abs(tb$proj_actual_lower - tb$dd_lower) < 1e-12; cl_hi <- abs(tb$proj_actual_upper - tb$dd_upper) < 1e-12
dif <- max(abs(cr$model_lower - tb$proj_actual_lower)[!cl_lo], abs(cr$model_upper - tb$proj_actual_upper)[!cl_hi])
report("K11 single-cell model interval = eiBounds proj_actual (fpc = popoviciu)", sum(!cl_lo) + sum(!cl_hi) > 0 && dif < 1e-10,
       sprintf("%d endpoints not clipped by DD (1500 simulated units): max difference %.1e (alpha split 0.9 / 0.1 alpha as in eiBounds)", sum(!cl_lo) + sum(!cl_hi), dif))
fitC <- eiBounds(cbind(yes, no) ~ cbind(g1, g2), data = dd, model = "M0", fpc = "popoviciu", ci = "conservative", seed = 11,
                 targets = list(c11 = c("g1", "yes"), c21 = c("g2", "yes")))
crc <- fitC$counts$cross$table; tbc <- fitC$table
cl_lo <- abs(tbc$ci_actual_lower - tbc$dd_lower) < 1e-12; cl_hi <- abs(tbc$ci_actual_upper - tbc$dd_upper) < 1e-12
difc <- max(abs(crc$model_lower - tbc$ci_actual_lower)[!cl_lo], abs(crc$model_upper - tbc$ci_actual_upper)[!cl_hi])
report("K11b conservative interval = eiBounds ci_actual (ci = conservative, fpc = popoviciu)", difc < 1e-10,
       sprintf("max difference %.1e: the cross-layer block redraws the simulation draws (Zs) instead of reusing the fit's, so its model interval is not the one printed in the model table", difc),
       finding = difc >= 1e-10)

# the reference population for K12-K14: msign_example_3x2 rounded to integer counts
d2 <- read.csv(tpath("msign_example_3x2.csv"))
Nn <- d2$N; X <- as.matrix(d2[, c("X1", "X2", "X3")]); Tm <- as.matrix(d2[, c("T1", "T2")])
lr <- function(N, s) { v <- floor(N * s); k <- N - sum(v); if (k > 0) { o <- order(-(N * s - v)); v[o[1:k]] <- v[o[1:k]] + 1 }; v }
A <- t(vapply(seq_along(Nn), function(i) lr(Nn[i], X[i, ]), numeric(3))); B <- t(vapply(seq_along(Nn), function(i) lr(Nn[i], Tm[i, ]), numeric(2)))
dd <- data.frame(u = paste0("P", seq_along(Nn)), g1 = A[, 1], g2 = A[, 2], g3 = A[, 3], yes = B[, 1], no = B[, 2])
fitP <- eiBounds(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = "M0", fpc = "popoviciu", seed = 11,
                 targets = list(c11 = c("g1", "yes"), c21 = c("g2", "yes"), c31 = c("g3", "yes")))

# K12 S_phi and t_E recomputed independently ---------------------------------------------------
Aw <- A; p <- nrow(A)
rho_rate <- function(rows, cols, S = 1:p) { den <- sum(Aw[S, rows]); r <- numeric(p)
  if (length(cols) < 2) r[S] <- rowSums(Aw[S, rows, drop = FALSE]) / den else r[S] <- 0   # all columns -> max - min = 0
  r }
tg <- list(rate("g1", "yes"), gap("g1", "g2", "yes"), ovr_gap("g1", "yes"), rate(c("g1", "g2"), "yes"),
           rate("g1", c("yes", "no")), rate("g1", "yes", where = 1:5), combo(rate("g1", "yes"), rate("g3", "yes"), weights = c(1, -2)))
ec2 <- eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, unit = "u", model = fitP, targets = tg)
ec2h <- eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, unit = "u", model = fitP, targets = tg, allowance = "hoeffding")
A1 <- sum(A[, 1]); A2 <- sum(A[, 2]); A3 <- sum(A[, 3]); N <- sum(A)
rhos <- list(A[, 1] / A1, A[, 1] / A1 + A[, 2] / A2,
             NULL,                                                        # ovr_gap: min over two representations
             (A[, 1] + A[, 2]) / (A1 + A2), rep(0, p), c(A[1:5, 1] / sum(A[1:5, 1]), rep(0, p - 5)), A[, 1] / A1 + 2 * A[, 3] / A3)
r_cells <- A[, 1] / A1 + (A[, 2] + A[, 3]) / (N - A1)                     # D over cells: range per unit
r_kappa <- A[, 1] * N / (A1 * (N - A1))                                  # kappa Y_rc - Y_c / (N - A_r)
rhos[[3]] <- if (sum(r_cells^2) <= sum(r_kappa^2)) r_cells else r_kappa
S <- vapply(rhos, function(r) 0.5 * sqrt(sum(r^2)), 0)
aE <- 0.1 * 0.05
tE <- vapply(rhos, function(r) sqrt(0.5 * log(2 / aE) * sum(r^2)), 0)
c2 <- ec2$cross$table
wid_est <- function(tab) tab$model_upper - tab$model_lower
dS <- max(abs(c2$S_phi - S))
# Hoeffding vs Popoviciu: the width difference is 2 (t_E - z S)
dW <- max(abs((wid_est(ec2h$cross$table) - wid_est(c2)) - 2 * (tE - qnorm(1 - aE / 2) * S)))
report("K12 S_phi, t_E, representation", dS < 1e-12 && dW < 1e-9,
       sprintf("7 targets (cell, gap, ovr_gap, merged rows, all columns, subset, combo): max |S_phi - independent| = %.1e; Hoeffding-minus-Popoviciu width vs 2(t_E - z S): %.1e; ovr_gap representation chosen: %s (S %.4f vs %.4f)",
               dS, dW, if (sum(r_cells^2) <= sum(r_kappa^2)) "cells" else "kappa", 0.5 * sqrt(sum(r_cells^2)), 0.5 * sqrt(sum(r_kappa^2))))

# K13 refusals and count-only targets ------------------------------------------------------------
fitN <- eiBounds(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = "M0", ci = "none")
r1 <- tryCatch(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = fitN), error = function(e) conditionMessage(e))
r2 <- tryCatch(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = list(table = fitP$table)), error = function(e) conditionMessage(e))
ecs <- eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = fitP, targets = list(share("g1", "yes")))
r3 <- tryCatch(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd[-1, ], model = fitP), error = function(e) conditionMessage(e))
report("K13 refusals", grepl("not an actual-target", r1) && grepl("never with a model estimate", r2) && is.na(ecs$cross$table$model_lower) &&
         grepl("outcome margins", ecs$cross$table$note) && grepl("not on the same units", r3),
       "ci = 'none' refused; a non-eiBounds object refused; share() gets the count layer only; different units refused")

# K14 wording (C8) -------------------------------------------------------------------------------
src <- unlist(lapply(list.files(file.path(Sys.getenv("PKG", "/home/user/eiid/package/eiPartialID"), "R"), full.names = TRUE), readLines))
given <- grep("given the count rows", src, value = TRUE)
given <- given[!grepl("Never print|never print", given)]
o_id <- capture.output(print(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = fitP, targets = list(rate("g1", "yes")),
                                      restrictions = list(known = list(unit = 1, row = "g1", col = "yes", value = A[1, 1] - 0)))))
o_bel <- capture.output(print(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = fitP, targets = list(rate("g1", "yes")),
                                       restrictions = list(order_unit = list(row = "g1", col = "yes")))))
o_tot <- capture.output(print(eiCounts(cbind(yes, no) ~ cbind(g1, g2, g3), data = dd, model = fitP, targets = list(rate("g1", "yes")),
                                       restrictions = list(total = list(row = "g1", col = "yes", lower = 10)))))
cov_bel <- "coverage at least 1 − α − P(belief rows false), which is unknown; level 1 − α only if the belief rows are taken as certain"
cov_bel_a <- "coverage at least 1 - alpha - P(belief rows false), which is unknown; level 1 - alpha only if the belief rows are taken as certain"
hasb <- function(o) any(grepl(cov_bel, o, fixed = TRUE) | grepl(cov_bel_a, o, fixed = TRUE))
report("K14 C8 wording", !length(given) && !hasb(o_id) && hasb(o_bel),
       sprintf("'given the count rows' in printed strings: %d; identity-only print uses the coverage line: %s; belief print has the C8 sentence verbatim: %s",
               length(given), any(grepl("identity rows only", o_id)), hasb(o_bel)))
report("K14b a stated floor as a 'total' defaults to an identity row", hasb(o_tot),
       sprintf("total with lower = 10 only: belief wording printed: %s (C8 lists stated floors among belief rows)", hasb(o_tot)), finding = !hasb(o_tot))

# K15 inherited status strings -------------------------------------------------------------------
d3 <- read.csv(tpath("rxc_example_3x3.csv"))
lr3 <- function(i, s, N) lr(N[i], s[i, ])
X3 <- as.matrix(d3[, c("X1", "X2", "X3")]); T3 <- as.matrix(d3[, c("T1", "T2", "T3")])
A3m <- t(vapply(seq_len(nrow(d3)), function(i) lr3(i, X3, d3$N), numeric(3))); B3m <- t(vapply(seq_len(nrow(d3)), function(i) lr3(i, T3, d3$N), numeric(3)))
dd3 <- data.frame(g1 = A3m[, 1], g2 = A3m[, 2], g3 = A3m[, 3], o1 = B3m[, 1], o2 = B3m[, 2], o3 = B3m[, 3], grp = d3$group)
st <- list()
for (m in c("M0", "M1", "M2", "L4", "L5", "L6")) {
  f <- suppressWarnings(eiBounds(cbind(o1, o2, o3) ~ cbind(g1, g2, g3), data = dd3, group = if (m %in% c("M1", "M2", "L5")) dd3$grp else NULL,
                                 model = m, row_constant = if (m == "L4") "g3" else NULL, seed = 1, targets = list(c11 = c("g1", "o1"))))
  st[[m]] <- f$counts$cross$table$status
}
fs <- eiBounds(cbind(o1, o2, o3) ~ cbind(g1, g2, g3), data = dd3, group = dd3$grp, model = "M1", seed = 1, targets = list(c11 = c("g1", "o1")),
               restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1))))
st[["M1+sign"]] <- fs$counts$cross$table$status
for (m in names(st)) cat(sprintf("    %-8s %s\n", m, st[[m]]))
report("K15 status strings follow the 6a table", grepl("^P \\(asymptotic\\)$", st$M0) && grepl("R16", st[["M1+sign"]]) &&
         grepl("regime R1", st$M1) && grepl("point-identified", st$L6),
       "M0: P (asymptotic); M1/M2: P only in R1, simulation otherwise; grouped + sign rows: no theorem (R16); L5/L6: labelled 'P (asymptotic; point-identified rung ...)' though 6a §4.3 has no row for it (developer point b)",
       finding = TRUE)
}

if (length(args) && args[1] %in% c("sim")) {
# ---------------------------------------------------------------------------
# K16: coverage and size of the R27 intersection (identity rows; belief rows)
# Goodman (L6) model true: unit rates drawn around constant means; target rate(g1, yes).
# ---------------------------------------------------------------------------
sim_pop <- function(p, seed, m = c(0.7, 0.3), sd = 0.12) {
  set.seed(seed)
  Nn <- sample(150:600, p, TRUE); x <- stats::rbeta(p, 2, 2.5)
  a1 <- stats::rbinom(p, Nn, x); a2 <- Nn - a1
  b1 <- pmin(pmax(stats::rnorm(p, m[1], sd), 0), 1); b2 <- pmin(pmax(stats::rnorm(p, m[2], sd), 0), 1)
  y1 <- stats::rbinom(p, a1, b1); y2 <- stats::rbinom(p, a2, b2)
  data.frame(u = paste0("P", 1:p), g1 = a1, g2 = a2, yes = y1 + y2, no = Nn - y1 - y2, n11 = y1, n21 = y2)
}
nrep <- as.integer(Sys.getenv("NREP", "300")); p <- 300; alpha <- 0.05
set.seed(2026); design <- sim_pop(p, 1)            # fixed design (N, a) across replications
known_units <- sort(sample(p, 30))                  # chosen from the design only (sorted: see issue on known-cell order)
res <- matrix(NA, nrep, 9, dimnames = list(NULL, c("cov_model", "cov_I", "rej_I", "cov_B1", "rej_B1", "row_B1", "cov_B2", "rej_B2", "row_B2")))
t_floor <- NA
for (k in seq_len(nrep)) {
  d <- design
  set.seed(10000 + k)
  b1 <- pmin(pmax(stats::rnorm(p, 0.7, 0.12), 0), 1); b2 <- pmin(pmax(stats::rnorm(p, 0.3, 0.12), 0), 1)
  d$n11 <- stats::rbinom(p, d$g1, b1); d$n21 <- stats::rbinom(p, d$g2, b2)
  d$yes <- d$n11 + d$n21; d$no <- d$g1 + d$g2 - d$yes
  truth <- sum(d$n11) / sum(d$g1)
  if (is.na(t_floor)) t_floor <- round(0.7 * sum(d$g1) + 0.2 * sqrt(sum(d$g1)))   # a floor that is false in a sizeable share of replications
  fit <- eiBounds(cbind(yes, no) ~ cbind(g1, g2), data = d, model = "L6", target = "model")
  run <- function(rs) {
    e <- eiCounts(cbind(yes, no) ~ cbind(g1, g2), data = d, unit = "u", model = fit, targets = list(rate("g1", "yes")),
                  restrictions = rs, binding = FALSE)
    ct <- e$cross$table
    inI <- !is.na(ct$lower) && truth >= ct$lower - 1e-12 && truth <= ct$upper + 1e-12
    c(cov = inI, rej = isTRUE(ct$disjoint), mcov = truth >= ct$model_lower && truth <= ct$model_upper)
  }
  rI <- run(list(known = list(unit = known_units, row = "g1", col = "yes", value = d$n11[known_units])))
  rB1 <- run(list(known = list(unit = known_units, row = "g1", col = "yes", value = d$n11[known_units]),
                  total = list(row = "g1", col = "yes", lower = t_floor, role = "belief")))
  rB2 <- run(list(sign = list(type = "dB", row = "g1", col = "yes")))
  res[k, ] <- c(rI["mcov"], rI["cov"], rI["rej"], rB1["cov"], rB1["rej"], sum(d$n11) >= t_floor, rB2["cov"], rB2["rej"],
                sum(d$n11) / sum(d$g1) - sum(d$n21) / sum(d$g2) >= 0)   # B2 row: delta_B >= 0 holds iff Y >= Y^NM; recorded below
}
m <- colMeans(res)
pf1 <- 1 - m["row_B1"]
rej_true1 <- mean(res[, "rej_B1"] & res[, "row_B1"] == 1)
cat(sprintf("    %d replications, %d units, L6 (Goodman) true, alpha = %.2f; known cells in %d design-chosen units\n", nrep, p, alpha, length(known_units)))
cat(sprintf("    model interval coverage %.3f; identity rows: intersection coverage %.3f, rejection rate %.3f\n", m["cov_model"], m["cov_I"], m["rej_I"]))
cat(sprintf("    belief floor (Y_11 >= %d, false in %.3f): coverage %.3f (bound 1 - alpha - P(false) = %.3f); P(reject and row true) = %.3f\n",
            t_floor, pf1, m["cov_B1"], 1 - alpha - pf1, rej_true1))
cat(sprintf("    belief sign delta_B >= 0: coverage %.3f, rejection rate %.3f\n", m["cov_B2"], m["rej_B2"]))
se <- function(q) sqrt(q * (1 - q) / nrep)
report("K16 R27 simulation", m["cov_I"] == m["cov_model"] | m["cov_I"] >= m["cov_model"] &&
         m["cov_I"] >= 1 - alpha - 2 * se(1 - alpha) && m["rej_I"] <= alpha + 2 * se(alpha) &&
         m["cov_B1"] >= 1 - alpha - pf1 - 2 * se(1 - alpha) && rej_true1 <= alpha + 2 * se(alpha),
       "identity rows: coverage = the model interval's, size <= alpha; belief rows: (*) bounds hold")
}

if (length(args) && args[1] %in% c("all", "issues")) {
# ---------------------------------------------------------------------------
# reproductions of the numbered issues of PHASE6a_code_review.md
# ---------------------------------------------------------------------------
tr <- function(expr) tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err"))
# issue 1: lpSolve MILP returns a suboptimal integer endpoint (oracle instance; brute force max = 1/3)
d <- data.frame(r1 = c(2, 1, 1), r2 = c(3, 1, 2), c1 = c(1, 2, 1), c2 = c(3, 0, 1), c3 = c(1, 0, 1))
f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2", "c3"), data = d, binding = FALSE,
              targets = list(gap("r1", "r1", "c2", where = 1:2, where2 = 3L)),
              restrictions = list(linear = list(cells = data.frame(unit = c(1L, 3L), row = "r1", col = c("c3", "c2"), coef = c(3, 0.5)), lower = 1.65)))
x <- c(0, 1, 1, 1, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0)   # an integer table attaining 1/3 (unit-major, r, c)
ns <- asNamespace("eiPartialID")
ok_x <- ns$.ec_check(f$P, f$built$common, x, TRUE)
report("I1 MILP optimality (lpSolve)", abs(f$table$upper - 1/3) < 1e-12,
       sprintf("integer upper reported %.4f (%s); the table x with value 1/3 satisfies every row exactly: %s", f$table$upper, f$table$integer_status, ok_x),
       finding = abs(f$table$upper - 1/3) > 1e-12)
# issue 2: known cells with units not in sorted order
d <- data.frame(u = c("U1", "U2"), r1 = c(5, 5), r2 = c(5, 5), c1 = c(5, 5), c2 = c(5, 5))
f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u", binding = FALSE,
              targets = list(u1 = rate("r1", "c1", where = "U1"), u2 = rate("r1", "c1", where = "U2")),
              restrictions = list(known = list(unit = c("U2", "U1"), row = "r1", col = "c1", value = c(4, 1))))
report("I2 known cells keep the user's unit order", all(abs(f$table$lower * 5 - c(1, 4)) < 1e-12),
       sprintf("unit = c('U2','U1'), value = c(4, 1): fitted U1 = %g, U2 = %g", f$table$lower[1] * 5, f$table$lower[2] * 5),
       finding = !all(abs(f$table$lower * 5 - c(1, 4)) < 1e-12))
# issue 3: eiBounds(x, t, n) with an N_i = 0 unit fails by default (0.3.1 and target = "model" work)
e1 <- tr(eiBounds(x = c(.3, .5, .4), t = c(.5, .4, .6), n = c(100, 0, 80), model = "L6"))
e2 <- tr(eiBounds(x = c(.3, .5, .4), t = c(.5, .4, .6), n = c(100, 0, 80), model = "L6", target = "model"))
report("I3 eiBounds with an N = 0 unit (x, t, n input)", !inherits(e1, "err"),
       sprintf("default target: %s; target = 'model': %s", if (inherits(e1, "err")) paste("ERROR", e1) else "ok", if (inherits(e2, "err")) "ERROR" else "ok"),
       finding = inherits(e1, "err"))
# issue 4: K7 above (linear coefficients as doubles)
# issue 5: W2 row with all-zero coefficients
d7 <- data.frame(r1 = c(2, 4), r2 = c(2, 4), r3 = c(1, 2), c1 = c(2, 5), c2 = c(3, 5))
e5 <- tr(eiCounts(rows = c("r1", "r2", "r3"), cols = c("c1", "c2"), data = d7, targets = list(rate("r1", "c1")),
                  restrictions = list(sign = list(type = "dW", version = "W2", row = "r3", col = "c1"))))
report("I5 W2 sign with constant X_r", !inherits(e5, "err"), sprintf("%s", if (inherits(e5, "err")) paste("ERROR:", e5) else "ok"), finding = inherits(e5, "err"))
# issue 7: where= on outcomes is treated as design-measurable
set.seed(1); p <- 300; N <- sample(100:400, p, TRUE); v <- rbinom(p, N, runif(p, .1, .6)); y11 <- rbinom(p, v, .7); y21 <- rbinom(p, N - v, .35)
dw <- data.frame(u = paste0("P", 1:p), vet = v, other = N - v, yes = y11 + y21, no = N - y11 - y21)
m <- eiBounds(cbind(yes, no) ~ cbind(vet, other), data = dw, model = "L6", target = "model")
e <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = dw, unit = "u", model = m, targets = list(rate("vet", "yes", where = yes / (yes + no) > 0.5)))
report("I7 outcome-selected unit set gets no cross-layer interval", !isTRUE(e$table$design),
       sprintf("where = yes/(yes+no) > 0.5: design flag %s, model interval [%.3f, %.3f] printed as covering", e$table$design,
               e$cross$table$model_lower, e$cross$table$model_upper), finding = isTRUE(e$table$design))
# issue 11: numeric units out of range are dropped silently
f <- tr(eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, unit = "u", targets = list(rate("r1", "c1")),
                 restrictions = list(box = list(units = c(1, 99), row = "r1", col = "c1", lower = 1))))
report("I11 unknown numeric unit refused", inherits(f, "err"), if (inherits(f, "err")) f else paste("accepted:", f$restrictions$restriction),
       finding = !inherits(f, "err"))
# issue 13: dB sign rows are labelled delta_W
d <- data.frame(r1 = c(1, 4, 2), r2 = c(4, 1, 3), c1 = c(2, 3, 1), c2 = c(3, 2, 4))
f <- eiCounts(rows = c("r1", "r2"), cols = c("c1", "c2"), data = d, targets = list(rate("r1", "c1")),
              restrictions = list(sign = list(type = "dB", row = "r1", col = "c1")))
report("I13 sign row label", grepl("delta_B", f$restrictions$restriction), f$restrictions$restriction, finding = !grepl("delta_B", f$restrictions$restriction))
}

if (length(args) && args[1] == "hang") {
# issue 16: pairwise delta_B / W2 rows at realistic unit sizes. Run under `timeout 120`; the
# 3-unit case returns in ~1-2 s, the 5-unit case does not return (lpSolve MILP, no time limit).
# Usage: timeout 120 Rscript phase6a_code_review_checks.R hang [pair|W2]
kind <- if (length(args) > 1) args[2] else "pair"
for (p in c(3, 5)) {
  set.seed(5); R <- 3; C <- 3
  tabs <- lapply(1:p, function(i) { N <- sample(200:1500, 1); x <- as.vector(rmultinom(1, N, c(.5, .3, .2)))
    t(vapply(1:R, function(r) as.vector(rmultinom(1, x[r], c(.5, .3, .2) + c(.1, 0, -.1) * (r - 2))), numeric(C))) })
  A <- t(vapply(tabs, rowSums, numeric(R))); B <- t(vapply(tabs, colSums, numeric(C)))
  d <- data.frame(A, B); names(d) <- c("g1", "g2", "g3", "o1", "o2", "o3")
  rs <- if (kind == "pair") list(sign = list(type = "dB", version = "pair", row = "g1", other = "g2", col = "o1")) else
    list(sign = list(type = "dW", version = "W2", row = "g1", col = "o1", sign = -1))
  t0 <- proc.time()[3]
  f <- eiCounts(cbind(o1, o2, o3) ~ cbind(g1, g2, g3), data = d, targets = list(rate("g1", "o1")), binding = FALSE, restrictions = rs)
  cat(sprintf("p = %d (%s): %.1f s, %s\n", p, kind, proc.time()[3] - t0, f$table$integer_status))
}
}
