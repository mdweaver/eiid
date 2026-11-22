# The count layer (M-COUNT), internal core: margins, rows on the realised table, exact
# rounding, total-unimodularity classes, LP / MILP and witnesses.
# Theory: docs/extension/PHASE6a_count_layer.md (R19, Lemma 2.1, Lemma C1, Prop. C2).
#
# Unknowns: realised cell counts n[i, r, c] >= 0 with unit row sums A[i, r] and column sums
# B[i, c]. Variable index k(i, r, c) = (i - 1) R C + (r - 1) C + c.
#
# Every count restriction becomes one or more rows  lo <= g' n <= hi  (lo = -Inf or hi = Inf
# for one-sided rows). A row has a 'kind':
#   "bound"  one cell (known cell, structural zero, box, one-vs-rest ordering, odds-ratio dial);
#   "total"  coefficient 1 on cell (r, c) over a unit set (published totals, W1 sign/CR rows);
#   "gen"    anything else (pairwise orderings, pairwise and W2 signs).
# With exactly integer margins every threshold is computed in exact rational arithmetic
# (package gmp) and rounded inward (6a §1.4, contract C4, review check RC9); with non-integer
# margins thresholds are doubles and no rounding is done (R19(e)).

# ---------------------------------------------------------------------------
# numbers: exact (bigq) or double
# ---------------------------------------------------------------------------
.ec_ceil <- function(x) if (inherits(x, c("bigq", "bigz"))) -floor(-x) else ceiling(x)
.ec_dbl <- function(x) if (is.null(x)) NULL else as.double(x)
.ec_is_q <- function(x) inherits(x, c("bigq", "bigz"))

.ec_vix <- function(P, i, r, c) (i - 1L) * P$R * P$C + (r - 1L) * P$C + c
.ec_unit_of <- function(P, k) (k - 1L) %/% (P$R * P$C) + 1L

# ---------------------------------------------------------------------------
# the problem: margins (units with N_i = 0 dropped, 6a §1.1)
# ---------------------------------------------------------------------------
.ec_problem <- function(A, B, rl, cl, units, integer_input = TRUE, data = NULL, tol = 1e-9) {
  A <- as.matrix(A); B <- as.matrix(B)
  storage.mode(A) <- "double"; storage.mode(B) <- "double"
  if (any(!is.finite(A)) || any(!is.finite(B))) stop("missing or non-finite counts")
  if (any(A < 0) || any(B < 0)) stop("counts must be nonnegative")
  NA_ <- rowSums(A); NB <- rowSums(B)
  if (any(abs(NA_ - NB) > tol * pmax(NA_, NB, 1)))
    stop("row and column totals differ in some units: the count layer needs the same population in both margins")
  keep <- NA_ > 0
  orig <- which(keep)
  A <- A[keep, , drop = FALSE]; B <- B[keep, , drop = FALSE]
  N <- rowSums(A)
  # integers only if exactly integer in the input (not N * X recomputed in floating point)
  int <- isTRUE(integer_input) && all(A == round(A)) && all(B == round(B)) && sum(N) < 2^53 &&
    all(rowSums(A) == rowSums(B))
  p <- nrow(A); R <- ncol(A); C <- ncol(B)
  P <- list(A = A, B = B, N = N, p = p, R = R, C = C, nv = p * R * C, rows = rl, cols = cl,
            units = as.character(units)[keep], orig = orig, n_orig = length(keep),
            n_dropped = sum(!keep), integer = int, data = data)
  P$num <- if (int) gmp::as.bigq else identity
  P$margin_all <- .ec_margin_trip(P)
  P
}

# margins as triplets (row id, var, value) with right-hand sides; equality rows
.ec_margin_trip <- function(P, units = seq_len(P$p), local = FALSE) {
  R <- P$R; C <- P$C
  trip <- vector("list", length(units) * (R + C)); rhs <- numeric(length(units) * (R + C)); j <- 0
  for (u in seq_along(units)) {
    i <- units[u]
    base <- if (local) 0L else (i - 1L) * R * C
    for (r in seq_len(R)) {
      j <- j + 1; ks <- base + (r - 1L) * C + seq_len(C)
      trip[[j]] <- cbind(j, ks, 1); rhs[j] <- P$A[i, r]
    }
    for (c in seq_len(C)) {
      j <- j + 1; ks <- base + (seq_len(R) - 1L) * C + c
      trip[[j]] <- cbind(j, ks, 1); rhs[j] <- P$B[i, c]
    }
  }
  list(trip = do.call(rbind, trip), rhs = rhs, n = j)
}

# ---------------------------------------------------------------------------
# rows
# ---------------------------------------------------------------------------
# lo, hi may be exact (bigq) or double; NULL = no bound on that side
.ec_row <- function(P, id, kind, idx, coef, lo = NULL, hi = NULL, cell = NULL, units = NULL,
                    label = NULL) {
  o <- order(idx); idx <- as.integer(idx[o])
  coefq <- if (.ec_is_q(coef)) coef[o] else if (P$integer) gmp::as.bigq(coef[o]) else NULL
  coefd <- as.double(if (.ec_is_q(coef)) coef[o] else coef[o])
  q <- function(x) if (is.null(x)) NULL else if (P$integer) gmp::as.bigq(x) else as.double(x)
  list(id = id, kind = kind, idx = idx, coef = coefd, coefq = coefq,
       loq = q(lo), hiq = q(hi),
       lo = if (is.null(lo)) -Inf else as.double(lo), hi = if (is.null(hi)) Inf else as.double(hi),
       cell = cell, units = units, label = label)
}

# inward rounding of one row (6a §1.4; [Chv]). Integer margins only.
# bound / total rows: lo -> ceiling, hi -> floor. gen rows: scale to coprime integers first.
.ec_round_row <- function(row) {
  if (is.null(row$coefq)) return(row)
  cq <- row$coefq
  if (row$kind %in% c("bound", "total")) s <- gmp::as.bigq(1) else {
    den <- gmp::denominator(cq); num <- abs(gmp::numerator(cq))
    l <- den[1]; if (length(den) > 1) for (j in 2:length(den)) l <- gmp::lcm.bigz(l, den[j])
    g <- num[1]; if (length(num) > 1) for (j in 2:length(num)) g <- gmp::gcd(g, num[j])
    s <- gmp::as.bigq(l) / gmp::as.bigq(g)
  }
  ci <- cq * s
  if (max(abs(as.double(ci))) > 2^40) {           # too large to pass to the solver as exact integers
    row$rounded <- FALSE; row$coef_rd <- row$coef; row$lo_rd <- row$lo; row$hi_rd <- row$hi
    return(row)
  }
  row$rounded <- TRUE
  row$coef_rd <- as.double(ci)
  row$lo_rd <- if (is.null(row$loq)) -Inf else as.double(.ec_ceil(row$loq * s))
  row$hi_rd <- if (is.null(row$hiq)) Inf else as.double(floor(row$hiq * s))
  row
}

# does point x satisfy the row? exact for integer margins and integer x; tolerance otherwise
.ec_row_ok <- function(row, x, P, exact, tol = 1e-7) {
  xv <- x[row$idx]
  if (exact && !is.null(row$coefq)) {
    v <- sum(row$coefq * gmp::as.bigq(xv))
    ok <- TRUE
    if (!is.null(row$loq)) ok <- ok && (v >= row$loq)
    if (!is.null(row$hiq)) ok <- ok && (v <= row$hiq)
    return(isTRUE(ok))
  }
  v <- sum(row$coef * xv); sc <- max(1, sum(abs(row$coef * xv)), abs(row$lo[is.finite(row$lo)]),
                                     abs(row$hi[is.finite(row$hi)]))
  (v >= row$lo - tol * sc) && (v <= row$hi + tol * sc)
}

.ec_margins_ok <- function(P, x, exact, tol = 1e-7) {
  n <- array(x, c(P$C, P$R, P$p))          # [c, r, i]
  rs <- matrix(apply(n, c(2, 3), sum), P$R)   # R x p
  cs <- matrix(apply(n, c(1, 3), sum), P$C)   # C x p
  if (exact) return(all(x >= 0) && all(rs == t(P$A)) && all(cs == t(P$B)))
  sc <- max(1, P$N)
  all(x >= -tol * sc) && all(abs(rs - t(P$A)) <= tol * sc) && all(abs(cs - t(P$B)) <= tol * sc)
}

# ---------------------------------------------------------------------------
# thresholds (exact with integer margins): one-vs-rest associations over a unit set S
# (6a Lemma C1; 2c Definition 1). slack in the package's probability units (dial S9):
# sign * delta >= -slack * scale, scale as in .rxc_assoc_scale.
# ---------------------------------------------------------------------------
.ec_design_S <- function(P, S, r) {
  nm <- P$num
  Ni <- nm(P$N[S]); a <- nm(P$A[S, r]); N <- sum(Ni)
  p <- Ni / N; X <- a / Ni
  mu <- sum(p * X); varX <- sum(p * X * X) - mu * mu
  list(Ni = Ni, a = a, N = N, p = p, X = X, mu = mu, varX = varX, Ar = sum(a))
}

.ec_wvar <- function(w, x) { w <- w / sum(w); sum(w * x * x) - sum(w * x)^2 }

.ec_ovr_thresholds <- function(P, S, r, c, slack_B = 0, slack_W = 0) {
  nm <- P$num
  d <- .ec_design_S(P, S, r)
  b <- nm(P$B[S, c]); Yc <- sum(b); Tc <- b / d$Ni
  out <- list(N = d$N, Ar = d$Ar, Yc = Yc, varX = d$varX, mu = d$mu)
  out$YNM <- sum(d$a * b / d$Ni)
  out$scale_B <- sum(d$p * d$X * (1 - d$X))
  if (d$Ar > 0 && d$Ar < d$N) {
    out$kappa <- d$N / (d$Ar * (d$N - d$Ar))
    if (as.double(d$varX) > 0) {
      out$ER <- (sum(d$p * d$X * Tc) - d$mu * sum(d$p * Tc)) / d$varX
      out$scale_W <- d$mu * .ec_wvar(d$p * d$X, d$X) + (1 - d$mu) * .ec_wvar(d$p * (1 - d$X), d$X)
    }
  }
  epsB <- nm(slack_B) * out$scale_B
  out$YNM_lo <- out$YNM - d$N * epsB          # delta_B >= -eps  <=>  Y >= YNM - N eps
  out$YNM_hi <- out$YNM + d$N * epsB          # delta_B <=  eps  <=>  Y <= YNM + N eps
  if (!is.null(out$ER)) {
    epsW <- nm(slack_W) * out$scale_W
    base <- out$ER + out$Yc / (d$N - d$Ar)
    out$YER_hi <- (base + epsW / d$varX) / out$kappa   # delta_W >= -eps  <=>  Y <= YER_hi
    out$YER_lo <- (base - epsW / d$varX) / out$kappa   # delta_W <=  eps  <=>  Y >= YER_lo
  }
  out
}

# generic (pairwise / W2) realised association rows: coefficients on n, plus scale
.ec_assoc_gen <- function(P, S, type, version, r, c, other = NULL) {
  nm <- P$num
  d <- .ec_design_S(P, S, r)
  idx <- integer(0); coef <- NULL
  addc <- function(k, v) { idx <<- c(idx, k); coef <<- if (is.null(coef)) v else c(coef, v) }
  if (type == "dB" && version == "pair") {
    s <- other; as_ <- nm(P$A[S, s])
    for (j in seq_along(S)) {
      i <- S[j]
      addc(.ec_vix(P, i, r, c), as_[j] / (d$N * d$Ni[j]))
      addc(.ec_vix(P, i, s, c), -d$a[j] / (d$N * d$Ni[j]))
    }
    scale <- sum(d$p * d$X * (as_ / d$Ni))
  } else if (type == "dB" && version == "W1") {
    for (j in seq_along(S)) {
      i <- S[j]
      for (s in setdiff(seq_len(P$R), r)) {
        as_ <- nm(P$A[i, s])
        addc(.ec_vix(P, i, r, c), as_ / (d$N * d$Ni[j]))
        addc(.ec_vix(P, i, s, c), -d$a[j] / (d$N * d$Ni[j]))
      }
    }
    scale <- sum(d$p * d$X * (1 - d$X))
  } else if (type == "dW" && version == "W2") {
    scale <- NULL
    for (g in seq_len(P$R)) {
      ag <- nm(P$A[S, g]); Xg <- ag / d$Ni; mug <- sum(d$p * Xg)
      if (as.double(mug) <= 0) next
      xb <- sum(d$p * Xg * d$X) / mug
      for (j in seq_along(S)) if (P$A[S[j], g] > 0) addc(.ec_vix(P, S[j], g, c), (d$X[j] - xb) / d$N)
      sc_g <- mug * .ec_wvar(d$p * Xg, d$X)
      scale <- if (is.null(scale)) sc_g else scale + sc_g
    }
  } else stop("internal: unknown generic association ", type, "/", version)
  # merge duplicate indices
  u <- sort(unique(idx))
  cm <- NULL
  for (k in u) { v <- sum(coef[idx == k]); cm <- if (is.null(cm)) v else c(cm, v) }
  list(idx = u, coef = cm, scale = scale)
}

# odds-ratio dial on a 2x2 unit: integer (or continuous) range of n11 with
# psi_L <= n11 n22 / (n12 n21) <= psi_U, written k n22 >= psi_L n12 n21 (6a §1.3 convention)
.ec_odds_box <- function(P, i, psiL, psiU) {
  a <- P$A[i, 1]; b <- P$B[i, 1]; N <- P$N[i]
  if (a == 0 || a == N || b == 0 || b == N) return(NULL)       # cell fixed by the margins: ignored
  L <- max(0, a + b - N); U <- min(a, b)
  g <- function(k, psi) k * (N - a - b + k) - psi * (a - k) * (b - k)
  root <- function(psi) {            # g increasing on [L, U]; g(L) <= 0 <= g(U)
    if (psi == 1) return(a * b / N)
    qa <- 1 - psi; qb <- N - a - b + psi * (a + b); qc <- -psi * a * b
    disc <- sqrt(max(0, qb^2 - 4 * qa * qc))
    rts <- c((-qb + disc) / (2 * qa), (-qb - disc) / (2 * qa))
    rts <- rts[rts >= L - 1e-9 * N & rts <= U + 1e-9 * N]
    if (!length(rts)) return(stats::uniroot(function(k) g(k, psi), c(L, U), tol = 1e-14)$root)
    min(max(rts[1], L), U)
  }
  lo <- NULL; hi <- NULL
  if (P$integer) {
    gq <- function(k, psi) { k <- gmp::as.bigq(k); kq <- k * (gmp::as.bigq(N - a - b) + k) -
      psi * (gmp::as.bigq(a) - k) * (gmp::as.bigq(b) - k); kq }
    if (psiL > 0) {
      ps <- gmp::as.bigq(psiL); k <- max(L, min(U, floor(root(psiL)) - 1))
      while (k <= U && gq(k, ps) < 0) k <- k + 1
      while (k - 1 >= L && gq(k - 1, ps) >= 0) k <- k - 1
      lo <- k                                    # may be U + 1: empty
    }
    if (is.finite(psiU)) {
      ps <- gmp::as.bigq(psiU); k <- min(U, max(L, ceiling(root(psiU)) + 1))
      while (k >= L && gq(k, ps) > 0) k <- k - 1
      while (k + 1 <= U && gq(k + 1, ps) <= 0) k <- k + 1
      hi <- k
    }
  } else {
    if (psiL > 0) lo <- root(psiL)
    if (is.finite(psiU)) hi <- root(psiU)
  }
  list(lo = lo, hi = hi)
}

# ---------------------------------------------------------------------------
# unit sets
# ---------------------------------------------------------------------------
# returns indices into the kept units
.ec_units <- function(P, units = NULL, where = NULL, env = parent.frame(), what = "units") {
  if (!is.null(where)) {
    if (inherits(where, "formula")) where <- where[[length(where)]]
    if (is.language(where)) {
      if (is.null(P$data)) stop(what, ": an expression in 'where' needs 'data'")
      where <- eval(where, P$data, env)
    }
    units <- where
  }
  if (is.null(units)) return(seq_len(P$p))
  if (is.logical(units)) {
    if (length(units) != P$n_orig) stop(what, ": a logical unit selection must have one value per row of the data")
    sel <- which(units & !is.na(units))
    return(sort(unique(stats::na.omit(match(sel, P$orig)))))
  }
  if (is.character(units) || is.factor(units)) {
    m <- match(as.character(units), P$units)
    if (any(is.na(m))) {
      # dropped (N = 0) units are silently out; unknown labels are an error
      stop(what, ": unknown unit label(s): ", paste(units[is.na(m)], collapse = ", "))
    }
    return(sort(unique(m)))
  }
  if (is.numeric(units)) {
    m <- match(as.integer(units), P$orig)
    return(sort(unique(stats::na.omit(m))))
  }
  stop(what, ": units must be labels, row numbers, a logical vector or an expression")
}

# ---------------------------------------------------------------------------
# restrictions -> rows (common rows, and CR unions)
# ---------------------------------------------------------------------------
.ec_types <- c("known", "zero", "box", "total", "order_unit", "odds", "sign", "cr")
.ec_default_role <- c(known = "identity", total = "identity", zero = "belief", box = "belief",
                      order_unit = "belief", odds = "belief", sign = "belief", cr = "belief")

.ec_label_of <- function(P, s, type) {
  rr <- function(r) paste(P$rows[r], collapse = "+"); cc <- function(c) paste(P$cols[c], collapse = "+")
  us <- function(S) if (length(S) == P$p) "all units" else if (length(S) <= 3)
    paste0("units {", paste(P$units[S], collapse = ", "), "}") else sprintf("%d units", length(S))
  switch(type,
         known = sprintf("known cell %s -> %s in %s = %s", rr(s$row), cc(s$col), us(s$units),
                         paste(format(s$value), collapse = ",")),
         zero = sprintf("structural zero %s -> %s in %s", rr(s$row), cc(s$col), us(s$units)),
         box = sprintf("box %s -> %s in %s: [%s, %s]", rr(s$row), cc(s$col), us(s$units),
                       format(s$lower), format(s$upper)),
         total = sprintf("total %s -> %s over %s in [%s, %s]", rr(s$row), cc(s$col), us(s$units),
                         format(s$lower), format(s$upper)),
         order_unit = sprintf("realised ordering in each unit (%s): %s %s %s on %s%s", us(s$units), rr(s$row),
                              if (s$sign > 0) ">=" else "<=",
                              if (is.null(s$over)) "rest" else rr(s$over), cc(s$col),
                              if (s$slack > 0) sprintf(" (slack %g)", s$slack) else ""),
         odds = sprintf("odds-ratio dial [%g, %g] in %s", s$psi[1], s$psi[2], us(s$units)),
         sign = sprintf("realised %s %s %s on %s (%s%s, %s)%s",
                        if (s$type == "dB") "delta_B" else "delta_W",
                        if (s$version == "pair") paste0("(", rr(s$row), " vs ", rr(s$other), ")") else rr(s$row),
                        if (s$sign > 0) ">= 0" else "<= 0", cc(s$col), s$version,
                        if (s$version == "W1" && s$type == "dB") " vs rest" else "", us(s$units),
                        if (s$slack > 0) sprintf(" (slack %g)", s$slack) else ""),
         cr = sprintf("realised CR %s on %s (%s%s, %s)%s", rr(s$row), cc(s$col), s$version,
                      if (s$version == "pair") paste0(" vs ", rr(s$other)) else "", us(s$units),
                      if (s$slack > 0) sprintf(" (slack %g)", s$slack) else ""))
}

.ec_idx <- function(v, labels, what) {
  if (is.null(v)) stop(what, " is required")
  if (is.numeric(v)) { if (any(v < 1 | v > length(labels))) stop(what, " out of range"); return(as.integer(v)) }
  i <- match(as.character(v), labels)
  if (any(is.na(i))) stop("unknown label(s) in ", what, ": ", paste(v[is.na(i)], collapse = ", "))
  i
}

# parse one user spec (labels) into an internal spec (indices), checking the object tag
.ec_parse_spec <- function(P, type, s, env) {
  if (!type %in% .ec_types) stop("unknown count restriction type '", type, "'")
  obj <- if (is.null(s$object)) "actual" else as.character(s$object)
  if (length(obj) != 1 || !obj %in% c("actual", "model"))
    stop("restriction '", type, "': object must be one of \"actual\" or \"model\"; the same restriction ",
         "may not be imposed under one label in both layers (realised and model beliefs are two assumptions: ",
         "give them as two restrictions)")
  if (obj == "model")
    stop("restriction '", type, "' has object = \"model\": rows on the model's coefficients belong to ",
         "eiBounds(restrictions = ...), never to the count layer (6a Remark 4.1)")
  role <- if (is.null(s$role)) .ec_default_role[[type]] else match.arg(s$role, c("identity", "belief"))
  ri <- function(v, w = "row") .ec_idx(v, P$rows, w); ci <- function(v, w = "col") .ec_idx(v, P$cols, w)
  units <- .ec_units(P, if (!is.null(s$units)) s$units else s$unit, s$where, env, type)
  sl <- if (is.null(s$slack)) 0 else as.numeric(s$slack)
  if (length(sl) != 1 || !is.finite(sl) || sl < 0) stop(type, ": 'slack' must be a nonnegative number")
  out <- list(type = type, object = obj, role = role, units = units, slack = sl)
  if (type %in% c("known", "zero", "box", "total", "order_unit", "sign", "cr")) {
    out$row <- ri(s$row); out$col <- ci(if (is.null(s$col)) 1 else s$col)
    if (length(out$row) != 1 || length(out$col) != 1) stop(type, ": 'row' and 'col' must name one cell")
  }
  if (type == "known") {
    if (is.null(s$value)) stop("known: 'value' is required")
    if (is.null(s$units) && is.null(s$unit) && is.null(s$where)) stop("known: 'unit' is required")
    out$value <- rep_len(as.numeric(s$value), length(units))
  } else if (type == "box") {
    out$lower <- if (is.null(s$lower)) -Inf else as.numeric(s$lower)
    out$upper <- if (is.null(s$upper)) Inf else as.numeric(s$upper)
  } else if (type == "total") {
    if (!is.null(s$value)) { out$lower <- out$upper <- as.numeric(s$value) } else {
      out$lower <- if (is.null(s$lower)) -Inf else as.numeric(s$lower)
      out$upper <- if (is.null(s$upper)) Inf else as.numeric(s$upper)
    }
  } else if (type == "order_unit") {
    ov <- if (is.null(s$over)) "rest" else s$over
    out$over <- if (identical(ov, "rest")) NULL else ri(ov, "over")
    if (!is.null(out$over) && out$over == out$row) stop("order_unit: 'over' must differ from 'row'")
    if (!is.null(out$over) && P$R == 2) out$over <- NULL           # R = 2: the other row is the rest
    out$sign <- sign(if (is.null(s$sign)) 1 else as.numeric(s$sign))
  } else if (type == "odds") {
    if (P$R != 2 || P$C != 2) stop("odds: the odds-ratio dial is defined for 2 x 2 tables only")
    ps <- as.numeric(s$psi)
    if (length(ps) != 2 || any(is.na(ps)) || ps[1] < 0 || ps[2] < ps[1]) stop("odds: 'psi' must be c(lower, upper) with 0 <= lower <= upper")
    out$psi <- ps
  } else if (type %in% c("sign", "cr")) {
    if (type == "sign") {
      out$type_a <- as.character(if (is.null(s$type)) "dB" else s$type)
      if (!out$type_a %in% c("dB", "dW")) stop("sign: type must be \"dB\" or \"dW\"")
      out$sign <- sign(if (is.null(s$sign)) 1 else as.numeric(s$sign))
      if (out$sign == 0) stop("sign: sign must be +1 or -1")
      allowed <- if (out$type_a == "dB") c("W1", "pair") else c("W1", "W2")
    } else allowed <- c("W1", "W2", "pair")
    out$version <- if (is.null(s$version)) "W1" else as.character(s$version)
    if (!out$version %in% allowed) stop(type, ": version must be one of ", paste(allowed, collapse = ", "))
    if (out$version == "pair") {
      out$other <- if (is.null(s$other)) { if (P$R == 2) setdiff(1:2, out$row) else stop(type, ": version \"pair\" needs 'other'") } else ri(s$other, "other")
      if (out$other == out$row) stop(type, ": 'other' must differ from 'row'")
    }
    if (length(units) < 2) stop(type, ": the unit set needs at least two units")
  }
  out
}

# rows of one parsed spec: list(rows = list(...), union = NULL or list(plus = rows, minus = rows))
.ec_spec_rows <- function(P, s, id) {
  rows <- list(); nm <- P$num
  add <- function(row) rows[[length(rows) + 1]] <<- row
  bnd <- function(i, r, c, lo, hi) add(.ec_row(P, id, "bound", .ec_vix(P, i, r, c), 1, lo, hi,
                                                cell = c(r, c), units = i))
  if (s$type == "known") {
    for (j in seq_along(s$units)) bnd(s$units[j], s$row, s$col, s$value[j], s$value[j])
  } else if (s$type == "zero") {
    for (i in s$units) bnd(i, s$row, s$col, 0, 0)
  } else if (s$type == "box") {
    for (i in s$units) bnd(i, s$row, s$col, if (is.finite(s$lower)) s$lower else NULL,
                            if (is.finite(s$upper)) s$upper else NULL)
  } else if (s$type == "total") {
    lo <- if (is.finite(s$lower)) s$lower else NULL; hi <- if (is.finite(s$upper)) s$upper else NULL
    if (length(s$units) == 1) bnd(s$units, s$row, s$col, lo, hi) else
      add(.ec_row(P, id, "total", .ec_vix(P, s$units, s$row, s$col), rep(1, length(s$units)), lo, hi,
                  cell = c(s$row, s$col), units = s$units))
  } else if (s$type == "order_unit") {
    r <- s$row; c <- s$col
    for (i in s$units) {
      a <- nm(P$A[i, r]); Ni <- nm(P$N[i]); b <- nm(P$B[i, c])
      if (is.null(s$over)) {
        if (P$A[i, r] == 0 || P$A[i, r] == P$N[i]) next          # no rest group, or no group r
        # n/a >= (b - n)/(N - a) - slack  <=>  n >= a b / N - slack a (N - a) / N   (a box; 6a §1.3)
        th <- a * b / Ni - nm(s$slack) * a * (Ni - a) / Ni
        if (s$sign > 0) bnd(i, r, c, th, NULL) else bnd(i, r, c, NULL, a * b / Ni + nm(s$slack) * a * (Ni - a) / Ni)
      } else {
        q <- s$over; aq <- nm(P$A[i, q])
        if (P$A[i, r] == 0 || P$A[i, q] == 0) next
        # n_r / a_r - n_q / a_q >= -slack  <=>  a_q n_r - a_r n_q >= -slack a_r a_q
        cf <- if (P$integer) c(aq, -a) else c(P$A[i, q], -P$A[i, r])
        if (s$sign < 0) cf <- -cf
        add(.ec_row(P, id, "gen", c(.ec_vix(P, i, r, c), .ec_vix(P, i, q, c)), cf,
                    -nm(s$slack) * a * aq, NULL))
      }
    }
  } else if (s$type == "odds") {
    for (i in s$units) {
      bx <- .ec_odds_box(P, i, s$psi[1], s$psi[2])
      if (is.null(bx) || (is.null(bx$lo) && is.null(bx$hi))) next
      bnd(i, 1, 1, bx$lo, bx$hi)
    }
  } else if (s$type == "sign") {
    S <- s$units; r <- s$row; c <- s$col
    if (s$version == "W1") {
      th <- .ec_ovr_thresholds(P, S, r, c, slack_B = s$slack, slack_W = s$slack)
      if (s$type_a == "dW" && is.null(th$ER))
        stop("sign (dW): Var_p(X_r) > 0 and 0 < A_r < N are needed over the sign's unit set (6a Lemma C1(d))")
      if (s$type_a == "dB") {
        if (s$sign > 0) lohi <- list(th$YNM_lo, NULL) else lohi <- list(NULL, th$YNM_hi)
      } else {
        if (s$sign > 0) lohi <- list(NULL, th$YER_hi) else lohi <- list(th$YER_lo, NULL)
      }
      add(.ec_row(P, id, "total", .ec_vix(P, S, r, c), rep(1, length(S)), lohi[[1]], lohi[[2]],
                  cell = c(r, c), units = S))
    } else {
      g <- .ec_assoc_gen(P, S, s$type_a, s$version, r, c, s$other)
      cf <- if (s$sign > 0) g$coef else -g$coef
      add(.ec_row(P, id, "gen", g$idx, cf, -nm(s$slack) * g$scale, NULL))
    }
  } else if (s$type == "cr") {
    S <- s$units; r <- s$row; c <- s$col
    mk <- function(tp, ver, sg) {
      if (ver == "W1") {
        th <- .ec_ovr_thresholds(P, S, r, c, slack_B = s$slack, slack_W = s$slack)
        if (tp == "dW" && is.null(th$ER))
          stop("cr: Var_p(X_r) > 0 and 0 < A_r < N are needed over the unit set (6a Lemma C1(d))")
        lohi <- if (tp == "dB") { if (sg > 0) list(th$YNM_lo, NULL) else list(NULL, th$YNM_hi) } else
          { if (sg > 0) list(NULL, th$YER_hi) else list(th$YER_lo, NULL) }
        .ec_row(P, id, "total", .ec_vix(P, S, r, c), rep(1, length(S)), lohi[[1]], lohi[[2]],
                cell = c(r, c), units = S)
      } else {
        g <- .ec_assoc_gen(P, S, tp, ver, r, c, s$other)
        .ec_row(P, id, "gen", g$idx, if (sg > 0) g$coef else -g$coef, -nm(s$slack) * g$scale, NULL)
      }
    }
    vB <- if (s$version == "pair") "pair" else "W1"; vW <- if (s$version == "W1") "W1" else "W2"
    return(list(rows = list(), union = list(plus = list(mk("dB", vB, 1), mk("dW", vW, 1)),
                                           minus = list(mk("dB", vB, -1), mk("dW", vW, -1)))))
  }
  list(rows = rows, union = NULL)
}

# ---------------------------------------------------------------------------
# TU classes of R19(b) for the whole row set of one system (contract C4)
# ---------------------------------------------------------------------------
.ec_cross <- function(S1, S2) {
  it <- length(intersect(S1, S2))
  it > 0 && it < length(S1) && it < length(S2)
}
.ec_laminar <- function(sets) {
  n <- length(sets)
  if (n < 2) return(TRUE)
  for (i in 1:(n - 1)) for (j in (i + 1):n) if (.ec_cross(sets[[i]], sets[[j]])) return(FALSE)
  TRUE
}
# two laminar families iff the crossing graph is bipartite (2-colouring by BFS)
.ec_two_laminar <- function(sets) {
  n <- length(sets)
  if (n < 2) return(TRUE)
  col <- rep(NA_integer_, n)
  for (s0 in seq_len(n)) {
    if (!is.na(col[s0])) next
    col[s0] <- 0L; queue <- s0
    while (length(queue)) {
      u <- queue[1]; queue <- queue[-1]
      for (v in seq_len(n)) {
        if (v == u || !.ec_cross(sets[[u]], sets[[v]])) next
        if (is.na(col[v])) { col[v] <- 1L - col[u]; queue <- c(queue, v) }
        else if (col[v] == col[u]) return(FALSE)
      }
    }
  }
  TRUE
}

# returns list(class = "(i)" etc. or NA, reason)
.ec_classify <- function(P, rows) {
  kinds <- vapply(rows, function(r) r$kind, "")
  if (any(kinds == "gen"))
    return(list(class = NA_character_, reason = "rows with unequal coefficients across cells or units (pairwise ordering, pairwise or W2 sign)"))
  tots <- rows[kinds == "total"]
  if (!length(tots)) return(list(class = "(i)", reason = "margins and single-cell bounds only"))
  cells <- t(vapply(tots, function(r) as.numeric(r$cell), c(0, 0)))
  sets <- lapply(tots, function(r) sort(r$units))
  key <- paste(cells[, 1], cells[, 2])
  if (length(unique(key)) == 1 && .ec_two_laminar(sets)) return(list(class = "(i)", reason = "one cell, two laminar families"))
  per_cell_lam <- function(keys) all(vapply(split(sets, keys), .ec_laminar, TRUE))
  if (length(unique(cells[, 2])) == 1 && per_cell_lam(key)) return(list(class = "(ii)", reason = "one column, one laminar family per cell"))
  if (length(unique(cells[, 1])) == 1 && per_cell_lam(key)) return(list(class = "(ii, transposed)", reason = "one row, one laminar family per cell"))
  if (P$R == 2 && per_cell_lam(cells[, 2])) return(list(class = "(iii)", reason = "R = 2, one laminar family per column"))
  if (P$C == 2 && per_cell_lam(cells[, 1])) return(list(class = "(iii, transposed)", reason = "C = 2, one laminar family per row"))
  list(class = NA_character_, reason = "totals outside the classes of R19(b) (mixture of classes, crossing families, or cells in different rows and columns)")
}

# ---------------------------------------------------------------------------
# LP / MILP (lpSolve). Basic optimal solutions are vertices (R19(d)).
# ---------------------------------------------------------------------------
.ec_lp_call <- function(nv, obj, trip, dir, rhs, sense, int) {
  sol <- lpSolve::lp(if (sense > 0) "max" else "min", objective.in = obj, const.dir = dir,
                     const.rhs = rhs, dense.const = trip, all.int = int)
  if (sol$status == 0) list(ok = TRUE, x = sol$solution) else list(ok = FALSE, x = NULL, status = sol$status)
}

# rows -> triplets; 'rd' uses the rounded rows
.ec_rows_trip <- function(rows, start, rd, map = NULL) {
  trip <- list(); dir <- character(0); rhs <- numeric(0); j <- start
  for (rw in rows) {
    cf <- if (rd) rw$coef_rd else rw$coef
    lo <- if (rd) rw$lo_rd else rw$lo; hi <- if (rd) rw$hi_rd else rw$hi
    idx <- if (is.null(map)) rw$idx else map[rw$idx]
    if (is.finite(lo) && is.finite(hi) && lo == hi) {
      j <- j + 1; trip[[length(trip) + 1]] <- cbind(j, idx, cf); dir <- c(dir, "="); rhs <- c(rhs, lo)
      next
    }
    if (is.finite(lo)) { j <- j + 1; trip[[length(trip) + 1]] <- cbind(j, idx, cf); dir <- c(dir, ">="); rhs <- c(rhs, lo) }
    if (is.finite(hi)) { j <- j + 1; trip[[length(trip) + 1]] <- cbind(j, idx, cf); dir <- c(dir, "<="); rhs <- c(rhs, hi) }
  }
  list(trip = if (length(trip)) do.call(rbind, trip) else NULL, dir = dir, rhs = rhs, n = j)
}

# north-west corner table for one unit (feasible; integer for integer margins)
.ec_nw <- function(a, b) {
  R <- length(a); C <- length(b); n <- matrix(0, R, C); r <- 1; c <- 1
  a <- a; b <- b
  while (r <= R && c <= C) {
    v <- min(a[r], b[c]); n[r, c] <- v; a[r] <- a[r] - v; b[c] <- b[c] - v
    if (a[r] <= 1e-12 * max(1, v) && r < R) r <- r + 1 else if (b[c] <= 1e-12 * max(1, v)) c <- c + 1 else r <- r + 1
  }
  n
}

# optimise obj' n over the margins and rows; rd = rounded rows; int = MILP.
# Separable problems (every row within one unit) are solved unit by unit.
.ec_optimise <- function(P, rows, obj, sense, rd = FALSE, int = FALSE) {
  RC <- P$R * P$C
  ru <- vapply(rows, function(rw) {
    u <- unique(.ec_unit_of(P, rw$idx)); if (length(u) == 1) u else 0L }, 0L)
  if (all(ru > 0)) {
    x <- numeric(P$nv)
    byunit <- split(seq_along(rows), factor(ru, levels = seq_len(P$p)))
    for (i in seq_len(P$p)) {
      ks <- (i - 1L) * RC + seq_len(RC)
      rr <- rows[byunit[[i]]]
      oi <- obj[ks]
      if (!length(rr) && all(oi == 0)) { x[ks] <- as.vector(t(.ec_nw(P$A[i, ], P$B[i, ]))); next }
      m <- .ec_margin_trip(P, i, local = TRUE)
      map <- integer(P$nv); map[ks] <- seq_len(RC)
      ex <- .ec_rows_trip(rr, m$n, rd, map)
      sol <- .ec_lp_call(RC, oi, rbind(m$trip, ex$trip), c(rep("=", m$n), ex$dir), c(m$rhs, ex$rhs), sense, int)
      if (!sol$ok) return(list(ok = FALSE))
      x[ks] <- sol$x
    }
    return(list(ok = TRUE, x = x, value = sum(obj * x)))
  }
  m <- P$margin_all
  ex <- .ec_rows_trip(rows, m$n, rd)
  sol <- .ec_lp_call(P$nv, obj, rbind(m$trip, ex$trip), c(rep("=", m$n), ex$dir), c(m$rhs, ex$rhs), sense, int)
  if (!sol$ok) return(list(ok = FALSE))
  list(ok = TRUE, x = sol$x, value = sum(obj * sol$x))
}

# check a solution against margins and every row (rounded rows for integer modes, and the
# original rows exactly); returns TRUE/FALSE
.ec_check <- function(P, rows, x, integer) {
  if (integer) {
    if (any(x != round(x))) return(FALSE)
    if (!.ec_margins_ok(P, x, TRUE)) return(FALSE)
    return(all(vapply(rows, function(rw) .ec_row_ok(rw, x, P, TRUE), TRUE)))
  }
  .ec_margins_ok(P, x, FALSE) && all(vapply(rows, function(rw) .ec_row_ok(rw, x, P, FALSE), TRUE))
}
