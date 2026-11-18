# Setup, plug-in bounds and inference for R x C bounds (internal).
# See docs/extension/PHASE3_inference.md.

.rxc_models <- c("DD", "M0", "M2", "M1delta", "M1", "L4", "L5", "L6")

.rxc_setup <- function(X, T, N, group = NULL, model = "M0", delta = 0, L4row = NULL,
                       range = c("hull", "trim", "extrapolate"), trim = NULL,
                       extrapolate = NULL, range_offset = 0, weights = c("none", "n"),
                       vcov = "HC3", cluster = NULL, extra = NULL) {
  range <- match.arg(range); weights <- match.arg(weights)
  model <- match.arg(model, setdiff(.rxc_models, "DD"))
  n <- nrow(X); R <- ncol(X); C <- ncol(T)
  if (is.null(group)) group <- rep(1L, n)
  group <- as.factor(group)
  # ---- trimming (polls kept out of the model are bounded by Duncan-Davis) ----
  keep <- rep(TRUE, n)
  if (range == "trim") {
    if (is.null(trim)) trim <- list(quantile = 0.01)
    if (!is.null(trim$quantile) && trim$quantile > 0) {
      q <- trim$quantile
      for (r in seq_len(R)) {
        lo <- stats::quantile(X[, r], q); hi <- stats::quantile(X[, r], 1 - q)
        keep <- keep & X[, r] >= lo & X[, r] <= hi
      }
    }
    if (!is.null(trim$min_n)) keep <- keep & N >= trim$min_n
    if (!is.null(trim$drop)) keep[trim$drop] <- FALSE
    if (sum(keep) < 2 * R * (R + 1)) stop("trimming leaves too few units")
  }
  g_all <- as.integer(droplevels(group))
  lev <- levels(droplevels(group))
  Xk <- X[keep, , drop = FALSE]; Tk <- T[keep, , drop = FALSE]; Nk <- N[keep]
  gk_f <- droplevels(group[keep]); gk <- as.integer(gk_f); G <- nlevels(gk_f)
  type <- switch(model,
                 M0 = "pooled", M1 = "shift", L4 = if (G > 1) "shift" else "pooled",
                 M2 = "block", M1delta = "block",
                 L5 = "goodman_block", L6 = "goodman")
  if (model %in% c("M0", "L6")) { gk <- rep(1L, length(gk)); G <- 1; type <- if (model == "M0") "pooled" else "goodman" }
  if (model %in% c("M1", "M2", "M1delta", "L5") && G < 2)
    stop("model ", model, " needs a grouping variable with at least two groups")
  # ---- regression ----
  if (type %in% c("goodman", "goodman_block")) {
    if (type == "goodman") Z <- Xk else {
      Z <- matrix(0, nrow(Xk), R * G)
      for (gg in seq_len(G)) Z[gk == gg, (gg - 1) * R + seq_len(R)] <- Xk[gk == gg, ]
    }
    des <- list(Z = Z, m = R, p = ncol(Z), type = type)
  } else {
    des <- .rxc_design(Xk, gk, G, type)
    type <- des$type
  }
  Z <- des$Z
  sw <- if (weights == "n") sqrt(Nk) else rep(1, length(Nk))
  Zw <- Z * sw; Yw <- Tk[, seq_len(C - 1), drop = FALSE] * sw
  qrZ <- qr(Zw)
  if (qrZ$rank < ncol(Zw)) {
    bad <- if (type %in% c("block", "shift", "goodman_block"))
      paste0(" (check that every group has enough units with varying composition; groups: ",
             paste(lev, collapse = ", "), ")") else ""
    stop("the regression design is rank deficient", bad, call. = FALSE)
  }
  B <- qr.coef(qrZ, Yw)
  B <- matrix(B, ncol = C - 1)
  E <- Yw - Zw %*% B
  vc <- .rxc_vcov(Zw, E, type = vcov, cluster = if (is.null(cluster)) NULL else cluster[keep])
  beta <- as.vector(B)
  # ---- composition range ----
  vert <- vector("list", G)
  for (gg in seq_len(G)) {
    Xg <- Xk[gk == gg, , drop = FALSE]
    vert[[gg]] <- if (range == "extrapolate") {
      if (is.null(extrapolate) || identical(extrapolate, "simplex")) diag(R) else as.matrix(extrapolate)
    } else .rxc_vertices(Xg, offset = range_offset)
  }
  # ---- unknowns ----
  d0 <- if (type %in% c("goodman", "goodman_block")) 0 else (C - 1) * R * (R - 1) / 2
  zlayout <- if (model %in% c("M2", "M1delta")) "group" else "common"
  zoff <- if (model == "M1delta") d0 else 0
  nz <- if (model == "M2") d0 * G else if (model == "M1delta") d0 * (G + 1) else d0
  w <- sweep(Nk * Xk, 2, colSums(Nk * Xk), "/")
  st <- list(X = Xk, T = Tk, N = Nk, g = gk, G = G, groups = levels(gk_f), R = R, C = C,
             model = model, type = type, delta = delta, L4row = L4row,
             mono = .rxc_mono(R), m = des$m, p = des$p, nbeta = length(beta),
             beta = beta, V = vc$V, vcov_type = vc$type, Z = Z,
             vertices = vert, d0 = d0, nz = nz, zlayout = zlayout, zoff = zoff, w = w,
             keep = keep, X_all = X, T_all = T, N_all = N, extra = extra, range = range,
             resid = E)
  st
}

# ---------------------------------------------------------------------------
# Duncan-Davis bounds for the actual rate (over a set of units)
# ---------------------------------------------------------------------------
.rxc_dd_cells <- function(X, T, N) {
  R <- ncol(X); C <- ncol(T)
  lo <- hi <- matrix(0, R, C)
  for (r in seq_len(R)) for (c in seq_len(C)) {
    den <- sum(N * X[, r])
    lo[r, c] <- sum(N * pmax(0, X[, r] + T[, c] - 1)) / den
    hi[r, c] <- sum(N * pmin(X[, r], T[, c])) / den
  }
  list(lower = lo, upper = hi)
}

.rxc_interval_lin <- function(a, lo, hi) {
  c(sum(ifelse(a >= 0, a * lo, a * hi)), sum(ifelse(a >= 0, a * hi, a * lo)))
}

# ---------------------------------------------------------------------------
# Candidate endpoints (dual vertices) and CLR-type intersection bounds
# ---------------------------------------------------------------------------
# 'relax = TRUE': when the system is infeasible at every draw, each draw's rows are shifted
# by the smallest uniform slack that makes them feasible. Any dual-feasible vertex is a
# valid CLR candidate (weak duality), so this only supplies candidates; it lets the CLR
# ends be computed, and cross, for (near-)empty sets.
.rxc_collect <- function(cons, tg, beta, Vsq, ndraw, side, relax = FALSE, draws = NULL) {
  g <- if (side == "upper") tg$gamma else -tg$gamma
  if (is.null(draws))
    draws <- c(list(beta), lapply(seq_len(ndraw), function(i) beta + as.vector(Vsq %*% stats::rnorm(ncol(Vsq)))))
  found <- list(); keys <- character(0)
  for (b in draws) {
    s <- as.vector(cons$S %*% b) + cons$s0
    if (relax) {
      tt <- .rxc_meet(numeric(0), matrix(0, 0, ncol(cons$A)), s, cons$A)   # max_z min_j (s + A z)_j
      if (is.finite(tt) && tt < 0) s <- s - tt + 1e-9
    }
    d <- .rxc_dual(s, cons$A, g)
    if (d$status != 0) next
    lam <- d$lambda
    key <- paste(round(lam, 8), collapse = ",")
    if (!(key %in% keys)) { keys <- c(keys, key); found[[length(found) + 1]] <- lam }
  }
  if (length(found) == 0) return(NULL)
  L <- do.call(rbind, found)
  if (side == "upper") {
    Gr <- sweep(L %*% cons$S, 2, tg$cvec, "+"); ic <- tg$c0 + as.vector(L %*% cons$s0)
  } else {
    Gr <- sweep(-L %*% cons$S, 2, tg$cvec, "+"); ic <- tg$c0 - as.vector(L %*% cons$s0)
  }
  list(G = Gr, ic = ic)
}

# upper confidence bound for theta_U = min_j theta_j (CLR 2013, analytic, with selection)
.rxc_clr_min <- function(theta, Gr, Zs, alpha, n) {
  se <- sqrt(pmax(rowSums((Gr %*% Zs$V) * Gr), 0))
  if (length(theta) == 1) {
    return(theta + stats::qnorm(1 - alpha) * se)
  }
  W <- (Zs$draws %*% t(Gr)) / matrix(pmax(se, 1e-300), nrow(Zs$draws), length(se), byrow = TRUE)
  W[, se == 0] <- 0
  gam_n <- 1 - 0.1 / log(max(n, 3))
  kp <- stats::quantile(apply(W, 1, max), gam_n, names = FALSE)
  upb <- min(theta + kp * se)
  sel <- theta <= upb + 2 * kp * se
  k <- stats::quantile(apply(W[, sel, drop = FALSE], 1, max), 1 - alpha, names = FALSE)
  min(theta[sel] + k * se[sel])
}

.rxc_zdraws <- function(V, nz_draw) {
  Vsq <- .rxc_sqrtm(V)
  list(V = V, Vsq = Vsq, draws = matrix(stats::rnorm(nz_draw * ncol(Vsq)), nz_draw) %*% t(Vsq))
}

# CLR confidence interval for the model-implied target (two sides at alpha_lo, alpha_up)
# Emptiness rule for this version (R10 §3.1): the CI is empty iff the ends cross. If no draw
# is feasible, candidates come from relaxed rows (see .rxc_collect); the CLR rule is kept and
# is not mixed with the conservative rule. The caller reports a crossing as "rejected".
.rxc_ci_clr <- function(st, cons, tg, alpha_lo, alpha_up, Zs, ndraw) {
  if (length(cons$disj)) stop("internal: CLR is not valid with disjunctive (slab) rows; use branches")
  draws <- c(list(st$beta), lapply(seq_len(ndraw), function(i)
    st$beta + as.vector(Zs$Vsq %*% stats::rnorm(ncol(Zs$Vsq)))))
  cu <- .rxc_collect(cons, tg, st$beta, Zs$Vsq, ndraw, "upper", draws = draws)
  cl <- .rxc_collect(cons, tg, st$beta, Zs$Vsq, ndraw, "lower", draws = draws)
  if (is.null(cu)) cu <- .rxc_collect(cons, tg, st$beta, Zs$Vsq, ndraw, "upper", relax = TRUE, draws = draws)
  if (is.null(cl)) cl <- .rxc_collect(cons, tg, st$beta, Zs$Vsq, ndraw, "lower", relax = TRUE, draws = draws)
  if (is.null(cu) || is.null(cl)) return(c(NA_real_, NA_real_))
  n <- nrow(st$X)
  U <- .rxc_clr_min(as.vector(cu$G %*% st$beta) + cu$ic, cu$G, Zs, alpha_up, n)
  Ln <- .rxc_clr_min(-(as.vector(cl$G %*% st$beta) + cl$ic), -cl$G, Zs, alpha_lo, n)
  c(-Ln, U)
}

# conservative version (Theorem P3.5(a)): inflate every constraint simultaneously
# Disjunctive (slab) rows, R10(b): every piece joins the critical-value maximum, each piece
# is inflated by its own k * se, and the max over pieces is taken afterwards. Emptiness rule
# for this version: the CI is empty iff the inflated set is empty.
.rxc_ci_conservative <- function(st, cons, tg, alpha, Zs) {
  DS <- if (length(cons$disj)) do.call(rbind, lapply(cons$disj, function(d) d$S)) else NULL
  Gall <- rbind(cons$S, DS, tg$cvec)
  se <- sqrt(pmax(rowSums((Gall %*% Zs$V) * Gall), 0))
  W <- (Zs$draws %*% t(Gall)) / matrix(pmax(se, 1e-300), nrow(Zs$draws), length(se), byrow = TRUE)
  W[, se == 0] <- 0
  last <- ncol(W); W[, last] <- abs(W[, last])
  k <- stats::quantile(apply(W, 1, max), 1 - alpha, names = FALSE)
  ns <- nrow(cons$S)
  cons2 <- cons; cons2$s0 <- cons$s0 + k * se[seq_len(ns)]
  if (length(cons$disj)) {
    np <- vapply(cons$disj, function(d) nrow(d$S), 0L)
    off <- ns + c(0, cumsum(np))
    infl <- lapply(seq_along(np), function(j) k * se[off[j] + seq_len(np[j])])
    cons2 <- .rxc_fix_disj(cons2, st$beta, infl)
  }
  pl <- .rxc_plugin(cons2, tg, st$beta)
  if (isTRUE(pl$empty) || isTRUE(pl$unbounded)) return(list(ci = c(NA_real_, NA_real_), empty = TRUE, k = k))
  list(ci = c(pl$lower - k * se[last], pl$upper + k * se[last]), empty = FALSE, k = k)
}

# finite-population allowance: sd bound for E = sum_i w_i e_i (one cell)
.rxc_fpc_sd <- function(st, cons, r, c, method, alpha_v, Zs, ndraw) {
  wr <- st$w[, r]
  sw2 <- sum(wr^2)
  if (method == "popoviciu") return(0.5 * sqrt(sw2))
  # Bhatia-Davis: Var(e_i) <= m_i (1 - m_i); by concavity sum pi_i m_i(1-m_i) <= mbar(1-mbar)
  # with pi_i = w_i^2 / sum w^2 and mbar = sum pi_i m_i, bounded by its own CI.
  wpi <- st$w; wpi[, r] <- wr^2 / sw2
  a <- matrix(0, st$R, st$C); a[r, c] <- 1
  tg <- .rxc_target(st, a, w = wpi)
  ci <- .rxc_ci_clr(st, cons, tg, alpha_v / 2, alpha_v / 2, Zs, ndraw)
  if (any(is.na(ci))) return(0.5 * sqrt(sw2))
  lo <- max(0, min(1, ci[1])); hi <- max(0, min(1, ci[2]))
  vmax <- if (lo <= 0.5 && hi >= 0.5) 0.25 else max(lo * (1 - lo), hi * (1 - hi))
  sqrt(sw2 * vmax)
}

# JKST CI_x for a target: plug-in candidates at beta-hat, S = S1 (Popoviciu) + se
.rxc_jkst <- function(st, cons, tg, pl, dd, xs, a) {
  V <- st$V
  gU <- tg$cvec + as.vector(t(cons$S) %*% pl$lam_up)
  gL <- tg$cvec - as.vector(t(cons$S) %*% pl$lam_lo)
  seU <- sqrt(max(0, sum(gU * (V %*% gU)))); seL <- sqrt(max(0, sum(gL * (V %*% gL))))
  s1 <- sum(sapply(seq_len(st$R), function(r) sum(abs(a[r, ])) * 0.5 * sqrt(sum(st$w[, r]^2))))
  SU <- s1 + seU; SL <- s1 + seL
  sel <- !(min(pl$upper, dd[2]) < max(pl$lower, dd[1]))
  out <- t(sapply(xs, function(x) {
    if (!sel) return(c(x, stats::pnorm(x), dd[1], dd[2], FALSE))
    c(x, stats::pnorm(x), max(pl$lower - x * SL, dd[1]), min(pl$upper + x * SU, dd[2]), TRUE)
  }))
  colnames(out) <- c("x", "nominal", "lower", "upper", "selected")
  out <- as.data.frame(out); out$selected <- out$selected == 1
  out
}

# ---------------------------------------------------------------------------
# Full estimation for a list of targets
# ---------------------------------------------------------------------------
.rxc_candidates_point <- function(tg) list(G = matrix(tg$cvec, 1), ic = tg$c0)

.rxc_estimate <- function(st, targets, alpha = 0.05, ci = c("clr", "conservative", "none"),
                          fpc = c("bhatia-davis", "popoviciu"), jkst = FALSE,
                          jkst_x = seq(0, 2, by = 0.25), eta = 0, ndraw = 100,
                          nz_draw = 4000, alpha_split = c(est = 0.04, E = 0.005, v = 0.005),
                          union = c("pieces", "hull")) {
  ci <- match.arg(ci); fpc <- match.arg(fpc); union <- match.arg(union)
  cons <- .rxc_constraints(st, eta = eta)
  # CR restrictions make the restricted set a union of systems (R3(c)). Point bounds and the
  # conservative CI use the W1/2x2 slab form (R10(b)); CLR uses the branch union (R10(c)).
  sys_pt <- .rxc_systems(cons, "slab")
  sys_ci <- .rxc_systems(cons, if (ci == "clr") "branch" else "slab")
  base <- sys_pt[[1]]; base$disj <- NULL                 # rung + linear rows (no CR)
  has_cr <- length(cons$cr) > 0
  Zs <- .rxc_zdraws(st$V, nz_draw)
  R <- st$R; C <- st$C
  share <- colSums(st$N * st$X) / colSums(st$N_all * st$X_all)          # kept share by row
  trimmed <- !st$keep
  dd_all <- .rxc_dd_cells(st$X_all, st$T_all, st$N_all)
  dd_trim <- if (any(trimmed)) .rxc_dd_cells(st$X_all[trimmed, , drop = FALSE],
                                            st$T_all[trimmed, , drop = FALSE],
                                            st$N_all[trimmed]) else NULL
  a_est <- alpha * alpha_split[["est"]] / sum(alpha_split)
  a_E <- alpha * alpha_split[["E"]] / sum(alpha_split)
  a_v <- alpha * alpha_split[["v"]] / sum(alpha_split)
  if (fpc == "popoviciu") { a_est <- a_est + a_v; a_v <- 0 }      # no variance-bound CI needed
  zE <- stats::qnorm(1 - a_E / 2)
  point <- st$nz == 0
  fpc_cache <- list()
  out <- list(); jk <- list(); pieces <- list()
  usum <- function(u) {
    sm <- .rxc_union_summary(u)
    if (union == "hull" && !is.null(u)) { sm$pieces <- 1L; sm$gap_lower <- sm$gap_upper <- NA_real_ }
    sm
  }
  for (nm in names(targets)) {
    a <- targets[[nm]]
    tg <- .rxc_target(st, a)
    wid_eta <- eta * sum(abs(a))
    pt <- .rxc_sys_interval(sys_pt, tg, st$beta)
    est_u <- pt$union
    if (!is.null(est_u)) { est_u[, 1] <- est_u[, 1] - wid_eta; est_u[, 2] <- est_u[, 2] + wid_eta }
    unb <- pt$unbounded
    if (unb) est_u <- NULL
    se_ <- usum(est_u)
    dd <- .rxc_interval_lin(a, dd_all$lower, dd_all$upper)
    row <- list(target = nm, estimate_lower = se_$lower, estimate_upper = se_$upper,
                ci_model_lower = NA, ci_model_upper = NA,
                ci_actual_lower = NA, ci_actual_upper = NA,
                dd_lower = dd[1], dd_upper = dd[2], empty = is.null(pt$union) && !unb,
                ci_rejected = FALSE, n_pieces = se_$pieces,
                estimate_gap_lower = se_$gap_lower, estimate_gap_upper = se_$gap_upper,
                ci_model_gap_lower = NA_real_, ci_model_gap_upper = NA_real_,
                ci_actual_gap_lower = NA_real_, ci_actual_gap_upper = NA_real_)
    cm_u <- ca_u <- NULL
    if (ci != "none") {
      rng <- .rxc_interval_lin(a, matrix(0, R, C), matrix(1, R, C))
      a2 <- sweep(a, 1, share, "*")
      tg2 <- .rxc_target(st, a2)
      # unit-level allowance and trimmed part do not depend on the system
      sdE <- 0
      for (r in seq_len(R)) for (c in seq_len(C)) {
        if (a2[r, c] == 0) next
        key <- paste(r, c)
        if (is.null(fpc_cache[[key]])) {
          fpc_cache[[key]] <- if (point) {
            # point rungs: Bhatia-Davis at the CI of the cell's own rate (or Popoviciu)
            if (fpc == "popoviciu") 0.5 * sqrt(sum(st$w[, r]^2)) else {
              wpi <- st$w; wpi[, r] <- st$w[, r]^2 / sum(st$w[, r]^2)
              e1 <- matrix(0, R, C); e1[r, c] <- 1
              tgp <- .rxc_target(st, e1, w = wpi)
              m <- sum(tgp$cvec * st$beta) + tgp$c0
              sem <- sqrt(max(0, sum(tgp$cvec * (st$V %*% tgp$cvec))))
              lo <- max(0, min(1, m - stats::qnorm(1 - a_v / 2) * sem))
              hi <- max(0, min(1, m + stats::qnorm(1 - a_v / 2) * sem))
              vmax <- if (lo <= 0.5 && hi >= 0.5) 0.25 else max(lo * (1 - lo), hi * (1 - hi))
              sqrt(sum(st$w[, r]^2) * vmax)
            }
          } else .rxc_fpc_sd(st, base, r, c, fpc, a_v, Zs, ndraw)   # rung + linear rows: valid
        }
        sdE <- sdE + abs(a2[r, c]) * fpc_cache[[key]]
      }
      tr <- if (!is.null(dd_trim)) {
        a3 <- sweep(a, 1, 1 - share, "*")
        .rxc_interval_lin(a3, dd_trim$lower, dd_trim$upper)
      } else c(0, 0)
      cms <- list(); cas <- list()
      for (sy in sys_ci) {
        # model-implied rate
        if (point) {
          pl <- .rxc_plugin(sy, tg, st$beta)
          se <- sqrt(max(0, sum(tg$cvec * (st$V %*% tg$cvec))))
          cm <- c(pl$lower, pl$upper) + c(-1, 1) * stats::qnorm(1 - alpha / 2) * se
          pl2 <- sum(tg2$cvec * st$beta) + tg2$c0
          se2 <- sqrt(max(0, sum(tg2$cvec * (st$V %*% tg2$cvec))))
          ce <- pl2 + c(-1, 1) * stats::qnorm(1 - a_est / 2) * se2
        } else if (ci == "clr") {
          cm <- .rxc_ci_clr(st, sy, tg, alpha / 2, alpha / 2, Zs, ndraw)
          ce <- .rxc_ci_clr(st, sy, tg2, a_est / 2, a_est / 2, Zs, ndraw)
        } else {
          cm <- .rxc_ci_conservative(st, sy, tg, alpha, Zs)$ci
          ce <- .rxc_ci_conservative(st, sy, tg2, a_est, Zs)$ci
        }
        # the model-implied rate is an average of valid rates: intersect with its logical range
        cms[[length(cms) + 1]] <- c(max(cm[1] - wid_eta, rng[1]), min(cm[2] + wid_eta, rng[2]))
        # actual rate: kept part (scaled by kept share) + trimmed part (Duncan-Davis)
        lo <- ce[1] - zE * sdE + tr[1] - eta * sum(abs(a2))
        hi <- ce[2] + zE * sdE + tr[2] + eta * sum(abs(a2))
        cas[[length(cas) + 1]] <- c(max(lo, dd[1]), min(hi, dd[2]))
      }
      # a branch whose interval is empty (inflated set empty, CLR ends crossed, or wholly
      # outside the logical range / Duncan-Davis) drops out; all empty = rejected (R10(c))
      cm_u <- .rxc_union(cms, tol = 0); ca_u <- .rxc_union(cas, tol = 0)
      sm <- usum(cm_u); sa <- usum(ca_u)
      row$ci_model_lower <- sm$lower; row$ci_model_upper <- sm$upper
      row$ci_model_gap_lower <- sm$gap_lower; row$ci_model_gap_upper <- sm$gap_upper
      row$ci_actual_lower <- sa$lower; row$ci_actual_upper <- sa$upper
      row$ci_actual_gap_lower <- sa$gap_lower; row$ci_actual_gap_upper <- sa$gap_upper
      row$ci_rejected <- is.null(cm_u) || is.null(ca_u)
    }
    if (jkst && !row$empty && !point && !has_cr) {
      pl <- .rxc_plugin(base, tg, st$beta)
      jk[[nm]] <- .rxc_jkst(st, base, tg, pl, dd, jkst_x, a)
    }
    pieces[[nm]] <- list(estimate = est_u, ci_model = cm_u, ci_actual = ca_u)
    out[[nm]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  list(table = res, jkst = if (jkst) jk else NULL, cons = cons, pieces = pieces)
}

# Emptiness test (Proposition P3.6(a), conservative version): reject iff the inflated set is
# empty. With CR, every system (slab form for W1/2x2, branches otherwise) is tested and the
# union is rejected only if all are empty (R10(c)).
.rxc_empty_test <- function(st, alpha = 0.05, nz_draw = 4000, detail = FALSE) {
  if (st$nz == 0) return(NA)
  cons <- .rxc_constraints(st)
  a <- matrix(0, st$R, st$C); a[1, 1] <- 1
  tg <- .rxc_target(st, a)
  Zs <- .rxc_zdraws(st$V, nz_draw)
  sys <- .rxc_systems(cons, "slab")
  emp <- vapply(sys, function(sy) isTRUE(.rxc_ci_conservative(st, sy, tg, alpha, Zs)$empty), logical(1))
  # LCE alone (rung without the optional restrictions), to name the rejected combination
  base_only <- NA
  if (all(emp) && !is.null(st$extra)) {
    st0 <- st; st0$extra <- NULL
    c0 <- .rxc_constraints(st0)
    base_only <- isTRUE(.rxc_ci_conservative(st0, c0, tg, alpha, Zs)$empty)
  }
  if (!detail) return(all(emp))
  list(rejected = all(emp), systems = emp, branches = vapply(sys, function(s) s$branch, ""),
       rung_alone_rejected = base_only)
}

# Wald test of M1's testable part (overidentifying restrictions, 2b Lemma D1(b))
.rxc_m1_test <- function(X, T, N, group, weights = "none", vcov = "HC3") {
  st <- tryCatch(.rxc_setup(X, T, N, group, model = "M2", weights = weights, vcov = vcov),
                 error = function(e) NULL)
  if (is.null(st)) return(list(statistic = NA, df = NA, p.value = NA,
                               note = "M2 regression not estimable (too few units per group)"))
  R <- st$R; C <- st$C; m <- st$m; p <- st$p; M <- st$mono
  jj <- function(r, s) which(M[, 1] == min(r, s) & M[, 2] == max(r, s))
  rows <- list()
  P <- .rxc_pairs(R)
  for (c in seq_len(C - 1)) for (g in 2:st$G) for (q in seq_len(nrow(P))) {
    r <- P[q, 1]; s <- P[q, 2]
    v <- numeric(st$nbeta)
    bo <- (c - 1) * p
    add <- function(gg, j, val) v[bo + (gg - 1) * m + j] <<- v[bo + (gg - 1) * m + j] + val
    # (Q^g - Q^1)_rs - ((Q^g - Q^1)_rr + (Q^g - Q^1)_ss)/2 = 0, with Q_rs = b_rs / 2
    add(g, jj(r, s), 0.5); add(1, jj(r, s), -0.5)
    add(g, jj(r, r), -0.5); add(1, jj(r, r), 0.5)
    add(g, jj(s, s), -0.5); add(1, jj(s, s), 0.5)
    rows[[length(rows) + 1]] <- v
  }
  Rm <- do.call(rbind, rows)
  d <- as.vector(Rm %*% st$beta)
  Vd <- Rm %*% st$V %*% t(Rm)
  e <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  Wd <- sum((t(e$vectors[, keep, drop = FALSE]) %*% d)^2 / e$values[keep])
  df <- sum(keep)
  list(statistic = Wd, df = df, p.value = stats::pchisq(Wd, df, lower.tail = FALSE),
       vcov_type = st$vcov_type)
}
