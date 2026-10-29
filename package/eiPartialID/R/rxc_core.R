# Core engine for R x C bounds (internal).
#
# Notation follows docs/extension/PHASE2a_rxc.md and PHASE2b_districts.md.
# For each column c, the homogeneous contextual-effect matrix H_c = Q_c + K_c:
# Q_c (symmetric) is identified from the regression of T_c on the quadratic
# monomials X_r X_s; K_c (skew-symmetric, sum over columns zero) is not.
# The unknown K is parametrised by z: K_c[r, s] = z_{c, (r, s)} for c < C and
# r < s, K_c[s, r] = -K_c[r, s], K_C = -sum_{c < C} K_c.
# Every constraint is linear: s(beta) + A z >= 0 with s(beta) = S beta + s0, and
# every target is c(beta) + gamma' z with c(beta) = cvec' beta + c0.

.rxc_pairs <- function(R) {
  p <- which(upper.tri(matrix(0, R, R)), arr.ind = TRUE)
  p <- p[order(p[, 1], p[, 2]), , drop = FALSE]
  colnames(p) <- c("r", "s")
  p
}

.rxc_mono <- function(R) {
  m <- which(upper.tri(matrix(0, R, R), diag = TRUE), arr.ind = TRUE)
  m <- m[order(m[, 1], m[, 2]), , drop = FALSE]
  colnames(m) <- c("r", "s")
  m
}

# coefficient vector a (length d0) with (K_c(z) v)_r = a' z (base block)
.rxc_Krow <- function(R, C, c, r, v) {
  P <- .rxc_pairs(R)
  np <- nrow(P)
  d0 <- (C - 1) * np
  a <- numeric(d0)
  if (np == 0) return(a)
  if (c < C) {
    for (p in seq_len(np)) {
      i <- P[p, 1]; j <- P[p, 2]
      if (i == r) a[(c - 1) * np + p] <- a[(c - 1) * np + p] + v[j]
      if (j == r) a[(c - 1) * np + p] <- a[(c - 1) * np + p] - v[i]
    }
  } else {
    for (cc in seq_len(C - 1)) a <- a - .rxc_Krow(R, C, cc, r, v)
  }
  a
}

# ---------------------------------------------------------------------------
# Regression design and covariance
# ---------------------------------------------------------------------------

# type: "pooled" (M0, L4 without groups), "shift" (M1, L4 with groups),
# "block" (M2, M1delta). Returns the design matrix and bookkeeping.
.rxc_design <- function(X, g, G, type) {
  R <- ncol(X)
  M <- .rxc_mono(R)
  mono <- sapply(seq_len(nrow(M)), function(j) X[, M[j, 1]] * X[, M[j, 2]])
  mono <- matrix(mono, nrow = nrow(X))
  m <- ncol(mono)
  if (type == "pooled" || G == 1) {
    Z <- mono
  } else if (type == "shift") {
    sh <- matrix(0, nrow(X), R * (G - 1))
    for (gg in 2:G) {
      sel <- g == gg
      sh[sel, (gg - 2) * R + seq_len(R)] <- X[sel, , drop = FALSE]
    }
    Z <- cbind(mono, sh)
  } else {
    Z <- matrix(0, nrow(X), m * G)
    for (gg in seq_len(G)) {
      sel <- g == gg
      Z[sel, (gg - 1) * m + seq_len(m)] <- mono[sel, , drop = FALSE]
    }
  }
  list(Z = Z, m = m, p = ncol(Z), type = if (G == 1) "pooled" else type)
}

# Heteroskedasticity- (or cluster-) robust covariance of stacked OLS
# coefficients for several responses sharing one design. Falls back
# HC3 -> HC1 -> HC0 -> classical when an estimator cannot be computed.
.rxc_vcov <- function(Z, E, type = "HC3", cluster = NULL) {
  n <- nrow(Z); k <- ncol(Z); C1 <- ncol(E)
  ZtZi <- solve(crossprod(Z))
  h <- rowSums((Z %*% ZtZi) * Z)
  build <- function(U) {
    V <- matrix(0, C1 * k, C1 * k)
    for (a in seq_len(C1)) for (b in a:C1) {
      if (is.null(cluster)) {
        Sab <- crossprod(Z, Z * (U[, a] * U[, b]))
      } else {
        Sa <- rowsum(Z * U[, a], cluster); Sb <- rowsum(Z * U[, b], cluster)
        Sab <- crossprod(Sa, Sb)
      }
      blk <- ZtZi %*% Sab %*% ZtZi
      ia <- (a - 1) * k + seq_len(k); ib <- (b - 1) * k + seq_len(k)
      V[ia, ib] <- blk; V[ib, ia] <- t(blk)
    }
    V
  }
  ok_psd <- function(V) {
    if (any(!is.finite(V))) return(FALSE)
    ev <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    min(ev) > -1e-10 * max(1, max(abs(ev)))
  }
  order <- if (!is.null(cluster)) c("CR1", "classical") else
    unique(c(type, "HC1", "HC0", "classical"))
  for (tp in order) {
    V <- NULL
    if (tp == "HC3") {
      if (any(h > 1 - 1e-8)) next
      V <- build(E / (1 - h))
    } else if (tp == "HC1") {
      if (n <= k) next
      V <- build(E) * n / (n - k)
    } else if (tp == "HC0") {
      V <- build(E)
    } else if (tp == "CR1") {
      Gc <- length(unique(cluster))
      if (Gc < 2 || n <= k) next
      V <- build(E) * (Gc / (Gc - 1)) * ((n - 1) / (n - k))
    } else if (tp == "classical") {
      df <- max(n - k, 1)
      Sig <- crossprod(E) / df
      V <- kronecker(Sig, ZtZi)
    }
    if (!is.null(V) && ok_psd(V)) return(list(V = (V + t(V)) / 2, type = tp))
  }
  stop("could not compute a covariance matrix for the regression coefficients")
}

# square root of a PSD matrix (for Gaussian draws)
.rxc_sqrtm <- function(V) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}

# ---------------------------------------------------------------------------
# Composition range (vertices where the linear model is imposed)
# ---------------------------------------------------------------------------

.rxc_vertices <- function(X, offset = 0) {
  R <- ncol(X)
  if (nrow(X) == 0) return(X)
  if (R == 2) {
    lo <- min(X[, 1]) + offset; hi <- max(X[, 1]) - offset
    V <- rbind(c(lo, 1 - lo), c(hi, 1 - hi))
  } else if (R == 3) {
    idx <- grDevices::chull(X[, 1], X[, 2])
    V <- X[idx, , drop = FALSE]
  } else {
    V <- unique(X)
  }
  V
}

# ---------------------------------------------------------------------------
# Linear pieces: identified part and unknown part of an implied rate
# ---------------------------------------------------------------------------

# identified part of e_r' H^(g)_c v as (coef on beta, constant)
.rxc_idlin <- function(st, g, c, r, v) {
  R <- st$R; C <- st$C; m <- st$m; p <- st$p
  M <- st$mono
  if (c == C) {
    coef <- numeric(st$nbeta); const <- sum(v)
    for (cc in seq_len(C - 1)) {
      pc <- .rxc_idlin(st, g, cc, r, v)
      coef <- coef - pc$coef; const <- const - pc$const
    }
    return(list(coef = coef, const = const))
  }
  coef <- numeric(st$nbeta)
  off <- (c - 1) * p
  if (st$type %in% c("goodman", "goodman_block")) {
    # H^(g)_c = alpha_{c,g} 1': rate = alpha_{c,g,r} * sum(v)
    goff <- if (st$type == "goodman_block") (g - 1) * R else 0
    coef[off + goff + r] <- sum(v)
    return(list(coef = coef, const = 0))
  }
  qoff <- if (st$type == "block") (g - 1) * m else 0
  for (j in seq_len(nrow(M))) {
    a <- M[j, 1]; b <- M[j, 2]
    if (a == r && b == r) coef[off + qoff + j] <- coef[off + qoff + j] + v[r]
    else if (a == r) coef[off + qoff + j] <- coef[off + qoff + j] + v[b] / 2
    else if (b == r) coef[off + qoff + j] <- coef[off + qoff + j] + v[a] / 2
  }
  if (st$type == "shift" && g >= 2) {
    idx <- off + m + (g - 2) * R + r
    coef[idx] <- coef[idx] + sum(v)
  }
  list(coef = coef, const = 0)
}

# unknown part (coefficients on z) of e_r' H^(g)_c v
.rxc_unklin <- function(st, g, c, r, v) {
  a <- .rxc_Krow(st$R, st$C, c, r, v)
  z <- numeric(st$nz)
  if (st$d0 == 0) return(z)
  if (st$zlayout == "common") {
    z[seq_len(st$d0)] <- a
  } else {
    off <- st$zoff + (g - 1) * st$d0
    z[off + seq_len(st$d0)] <- a
  }
  z
}

.rxc_rate <- function(st, g, c, r, v) {
  id <- .rxc_idlin(st, g, c, r, v)
  list(coef = id$coef, const = id$const, a = .rxc_unklin(st, g, c, r, v))
}

# ---------------------------------------------------------------------------
# Constraint system
# ---------------------------------------------------------------------------

.rxc_add_row <- function(acc, coef, const, a, kind) {
  acc$S[[length(acc$S) + 1]] <- coef
  acc$s0 <- c(acc$s0, const)
  acc$A[[length(acc$A) + 1]] <- a
  acc$kind <- c(acc$kind, kind)
  acc
}

.rxc_constraints <- function(st, eta = 0) {
  acc <- list(S = list(), s0 = numeric(0), A = list(), kind = character(0))
  R <- st$R; C <- st$C
  if (st$type %in% c("goodman", "goodman_block")) {
    return(list(S = matrix(0, 0, st$nbeta), s0 = numeric(0), A = matrix(0, 0, 0),
                kind = character(0)))
  }
  # (1) implied rates are valid probabilities at every vertex of each group's range
  for (g in seq_len(st$G)) {
    Vg <- st$vertices[[g]]
    for (k in seq_len(nrow(Vg))) for (c in seq_len(C)) for (r in seq_len(R)) {
      rt <- .rxc_rate(st, g, c, r, Vg[k, ])
      acc <- .rxc_add_row(acc, rt$coef, rt$const + eta, rt$a, "validity")
    }
  }
  # (2) M1(delta): district skew parts within delta of the M1 prediction
  if (st$model == "M1delta" && st$d0 > 0) {
    P <- .rxc_pairs(R); np <- nrow(P); m <- st$m; p <- st$p
    dg <- which(st$mono[, 1] == st$mono[, 2])            # positions of Q_rr in a block
    for (g in seq_len(st$G)) for (c in seq_len(C - 1)) for (q in seq_len(np)) {
      j <- (c - 1) * np + q
      r <- P[q, 1]; s <- P[q, 2]
      cg <- numeric(st$nbeta)
      if (g >= 2) {
        bo <- (c - 1) * p
        cg[bo + (g - 1) * m + dg[r]] <- cg[bo + (g - 1) * m + dg[r]] + 0.5
        cg[bo + dg[r]] <- cg[bo + dg[r]] - 0.5
        cg[bo + (g - 1) * m + dg[s]] <- cg[bo + (g - 1) * m + dg[s]] - 0.5
        cg[bo + dg[s]] <- cg[bo + dg[s]] + 0.5
      }
      a <- numeric(st$nz)
      a[st$zoff + (g - 1) * st$d0 + j] <- 1
      a[j] <- -1
      acc <- .rxc_add_row(acc, -cg, st$delta, a, "link")
      acc <- .rxc_add_row(acc, cg, st$delta, -a, "link")
    }
  }
  # (3) L4: designated row has no contextual effect: H_c[r0, s] = H_c[r0, r0]
  if (!is.null(st$L4row)) {
    r0 <- st$L4row
    e <- diag(R)
    for (c in seq_len(C - 1)) for (s in setdiff(seq_len(R), r0)) {
      a1 <- .rxc_rate(st, 1, c, r0, e[s, ]); a0 <- .rxc_rate(st, 1, c, r0, e[r0, ])
      coef <- a1$coef - a0$coef; const <- a1$const - a0$const; a <- a1$a - a0$a
      acc <- .rxc_add_row(acc, coef, const, a, "L4")
      acc <- .rxc_add_row(acc, -coef, -const, -a, "L4")
    }
  }
  # (4) optional restrictions
  acc <- .rxc_extra_rows(st, acc)
  S <- do.call(rbind, acc$S); A <- do.call(rbind, acc$A)
  if (is.null(dim(A))) A <- matrix(A, ncol = st$nz)
  list(S = S, s0 = acc$s0, A = A, kind = acc$kind)
}

# optional restrictions; each is applied in every group (shifts cancel in differences)
.rxc_extra_rows <- function(st, acc) {
  R <- st$R; C <- st$C; e <- diag(R)
  ex <- st$extra
  if (is.null(ex)) return(acc)
  groups <- seq_len(st$G)
  rate_diff <- function(g, c, r1, v1, r2, v2) {
    a1 <- .rxc_rate(st, g, c, r1, v1); a2 <- .rxc_rate(st, g, c, r2, v2)
    list(coef = a1$coef - a2$coef, const = a1$const - a2$const, a = a1$a - a2$a)
  }
  # structural zeros: rate of row r to column c is 0 throughout the range
  for (zr in ex$zeros) {
    gs <- if (is.null(zr$group)) groups else zr$group
    for (g in gs) {
      Vg <- st$vertices[[g]]
      for (k in seq_len(nrow(Vg))) {
        rt <- .rxc_rate(st, g, zr$col, zr$row, Vg[k, ])
        acc <- .rxc_add_row(acc, -rt$coef, -rt$const, -rt$a, "zero")
      }
    }
  }
  # ordering: mean rate of row 'high' to column c >= that of row 'low', at every vertex
  for (od in ex$order) {
    for (g in groups) {
      Vg <- st$vertices[[g]]
      for (k in seq_len(nrow(Vg))) {
        d <- rate_diff(g, od$col, od$high, Vg[k, ], od$low, Vg[k, ])
        acc <- .rxc_add_row(acc, d$coef, d$const, d$a, "order")
      }
    }
  }
  # sign of a contextual effect: moving composition from row 'from' to row 'to'
  # changes the mean rate of row r to column c in the direction 'sign'
  for (sg in ex$sign) {
    for (g in groups) {
      d <- rate_diff(g, sg$col, sg$row, e[sg$to, ], sg$row, e[sg$from, ])
      acc <- .rxc_add_row(acc, sg$sign * d$coef, sg$sign * d$const, sg$sign * d$a, "sign")
    }
  }
  # bounded slopes: |H_c[r, s] - H_c[r, s']| <= lambda for all r, c, s, s'
  if (!is.null(ex$slope)) {
    lam <- ex$slope
    for (g in groups) for (c in seq_len(C)) for (r in seq_len(R)) {
      for (s1 in seq_len(R - 1)) for (s2 in (s1 + 1):R) {
        d <- rate_diff(g, c, r, e[s2, ], r, e[s1, ])
        acc <- .rxc_add_row(acc, d$coef, d$const + lam, d$a, "slope")
        acc <- .rxc_add_row(acc, -d$coef, -d$const + lam, -d$a, "slope")
      }
    }
  }
  # generic linear restrictions on H entries: sum value * H_c[r, s] >= rhs
  for (lc in ex$linear) {
    gs <- if (is.null(lc$group)) groups else lc$group
    for (g in gs) {
      coef <- numeric(st$nbeta); const <- 0; a <- numeric(st$nz)
      for (q in seq_len(nrow(lc$terms))) {
        tq <- lc$terms[q, ]
        rt <- .rxc_rate(st, g, tq[["col"]], tq[["row"]], e[tq[["context"]], ])
        coef <- coef + tq[["value"]] * rt$coef; const <- const + tq[["value"]] * rt$const
        a <- a + tq[["value"]] * rt$a
      }
      acc <- .rxc_add_row(acc, coef, const - lc$rhs, a, "linear")
    }
  }
  # external benchmarks on a target functional: lower <= phi <= upper
  for (eb in ex$external) {
    tg <- .rxc_target(st, eb$weights)
    acc <- .rxc_add_row(acc, tg$cvec, tg$c0 - eb$lower, tg$gamma, "external")
    acc <- .rxc_add_row(acc, -tg$cvec, eb$upper - tg$c0, -tg$gamma, "external")
  }
  acc
}

# Target sum_{r,c} a_rc * Btilde_rc with sample weights w_ir = N_i X_ir / sum N X_r
# (or user-supplied row weights w, each column summing to one)
.rxc_target <- function(st, a, w = st$w) {
  R <- st$R; C <- st$C
  cvec <- numeric(st$nbeta); c0 <- 0; gam <- numeric(st$nz)
  for (r in seq_len(R)) for (c in seq_len(C)) {
    if (a[r, c] == 0) next
    wr <- w[, r]
    for (g in seq_len(st$G)) {
      sel <- st$g == g
      if (!any(sel)) next
      vbar <- colSums(st$X[sel, , drop = FALSE] * wr[sel])
      rt <- .rxc_rate(st, g, c, r, vbar)
      cvec <- cvec + a[r, c] * rt$coef; c0 <- c0 + a[r, c] * rt$const
      gam <- gam + a[r, c] * rt$a
    }
  }
  list(cvec = cvec, c0 = c0, gamma = gam)
}

# ---------------------------------------------------------------------------
# LP: dual form min s' lambda s.t. A' lambda = -g, lambda >= 0
# max {g' z : s + A z >= 0} = min {s' lambda : A' lambda = -g, lambda >= 0}
# ---------------------------------------------------------------------------

.rxc_dual <- function(s, A, g) {
  nz <- ncol(A)
  if (nz == 0) {
    feas <- all(s >= -1e-10)
    return(list(status = if (feas) 0 else 3, value = 0, lambda = numeric(length(s))))
  }
  sol <- lpSolve::lp("min", objective.in = s, const.mat = t(A),
                     const.dir = rep("=", nz), const.rhs = -g)
  # lpSolve: 0 optimal, 2 infeasible (primal unbounded), 3 unbounded (primal infeasible)
  list(status = sol$status, value = sol$objval, lambda = sol$solution)
}

# plug-in bounds for a target at coefficient vector beta
.rxc_plugin <- function(cons, tg, beta) {
  s <- as.vector(cons$S %*% beta) + cons$s0
  cb <- sum(tg$cvec * beta) + tg$c0
  up <- .rxc_dual(s, cons$A, tg$gamma)
  lo <- .rxc_dual(s, cons$A, -tg$gamma)
  if (up$status == 3 || lo$status == 3) return(list(empty = TRUE))
  if (up$status != 0 || lo$status != 0) return(list(empty = FALSE, unbounded = TRUE))
  list(empty = FALSE, lower = cb - lo$value, upper = cb + up$value,
       lam_up = up$lambda, lam_lo = lo$lambda)
}
