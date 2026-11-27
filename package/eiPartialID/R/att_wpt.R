# Changes and DiD, causal layer in expectation: within-precinct (WPT) and pooled (PPT) parallel trends
# in expected rates. Theory: docs/extension/PHASE2e_did_general.md (2e) Sec. 5.8-5.10:
#   Prop. E28 (identity, rate limits, sharp set under a known law; (d) no comparison voters;
#   (e) FD-LCE projection; (f) spillover sensitivity), Prop. E29 (realised target D_W and the Hoeffding
#   allowance), Prop. E30 (baseline / post / change LCE in one LP, baseline-curvature allowance epsilon),
#   Prop. E33(a) (local trend dial |b_it - b_i,t-1| <= M).
# Deferred (not built): E31, E32, E33(b)-(c), E34.
#
# Mean-model LP (E28(c), E30, E33(a)). Variables, for every precinct i and election e in {base, posts}:
#   mV[i,e], mN[i,e]  expected factual rates of the treated (V) and comparison (N) groups (only where
#                     the group exists), with x_i mV + (1 - x_i) mN = E[T_ie] and mV within the averaged
#                     Duncan-Davis ends of the law (the Aumann expectation, 2 x 2);
#   WPT: b[i,k], S[i,k] on supported precincts (a_iN >= a_min), theta[i,k] on unsupported ones;
#   PPT: b[k], S[k] (one aggregate gap and spillover per post-period);
#   LCE coefficients (alpha, beta) per group and profile; eps and M (fixed, or free for a breakdown LP).
# Target (Lemma E14 in means): ATT_k = sum_S wV (dV - dN + S - b) + sum_U wV (mV_t - theta)  (WPT),
#                              ATT_k = sum wV dV - sum wN dN + S - b                          (PPT).

#' Parallel trends in expectation: within each precinct, or pooled
#'
#' Trend assumptions for \code{\link{eiATT}}. Neither is a default: the user chooses (Def. 5.1 item 2).
#'
#' \strong{Realised reading (default of \code{within_pt()}, Prop. E35).} The target is the realised
#' \eqn{ATT_t = \sum_i w_{Vi}(B_{iVt} - \theta_i)} on the actual treated persons. In every supported
#' precinct the treated group's untreated change equals the comparison group's realised change up to
#' \eqn{\tau_i}: \eqn{|(\theta_i - B_{iV0}) - (B_{iNt} - S_i - B_{iN0})| \le \tau_i}, with default
#' \eqn{\tau_i = M\sqrt{1/a_{iV} + 1/a_{iN}}} (or a user vector \code{tau}), one spillover \eqn{|S_i| \le \Lambda}
#' per supported precinct, and \eqn{\theta_i \in [0,1]} on unsupported precincts (share \eqn{\pi_U} reported).
#' The set is an identification set of the count layer (continuous tables), with the breakdown value
#' \eqn{M^*} of its sign from one LP. \eqn{M} must hold in every precinct at once: it is a tolerance, not a
#' confidence level. With \code{lce = TRUE} in \code{eiATT()} the fit's model layer enters only through
#' E13 (a row \eqn{D_W^S(n) \in CI}), and the result is a confidence set.
#'
#' \strong{Expected-rate reading} (\code{within_pt(target = "expected")}, \code{pooled_pt()}; 2e Sec. 5.8):
#' the estimand is the expected programme effect on the same treated persons,
#' \eqn{ATT^{mean}_t = \sum_i w_{Vi}(m_{iVt} - m^{(0)}_{iVt})}.
#' \itemize{
#'   \item \code{within_pt(a_min, M)}: within-precinct parallel trends (WPT) on the supported precincts
#'     \eqn{\{i: a_{iN} \ge a_{min}\}}; a precinct's weight in the comparison trend is its share of the
#'     treated group. Unsupported precincts are bounded with \eqn{\theta_i \in [0,1]} and their share
#'     \eqn{\pi_U} is always reported. \code{M > 0} is the local trend dial of Prop. E33(a),
#'     \eqn{|b_{i,t} - b_{i,t-1}| \le M} in every supported precinct (exact WPT at \code{M = 0}).
#'   \item \code{pooled_pt(M)}: pooled parallel trends (PPT) in expectation; the weights are each
#'     group's own shares; \code{M} bounds the aggregate gap's increments.
#' }
#' The realised-table assumption \code{\link{parallel_trends}()} (Thm E16) is a different statement and
#' stays available under its own label.
#' @param a_min support threshold: a precinct is supported when it has at least \code{a_min}
#'   comparison persons (\eqn{a_{min} \ge 1}).
#' @param M nonnegative tolerance scale: realised WPT, \eqn{\tau_i = M\sqrt{1/a_{iV}+1/a_{iN}}} (E35);
#'   expected WPT/PPT, the local or aggregate dial on the increments of the expected trend gap (E33(a)).
#' @param target \code{"realised"} (default, Prop. E35) or \code{"expected"} (Props. E28-E30).
#' @param tau optional nonnegative tolerance per unit (or one number) for the realised reading; replaces
#'   the default \eqn{\tau_i}, and the breakdown value then scales the whole vector.
#' @return an object of class \code{"edTrend"}.
#' @export
within_pt <- function(a_min = 1, M = 0, target = c("realised", "expected"), tau = NULL) {
  target <- match.arg(target)
  if (length(a_min) != 1 || !is.numeric(a_min) || !is.finite(a_min) || a_min < 1) stop("within_pt(a_min): one number >= 1")
  if (!is.null(tau) && (target != "realised" || !is.numeric(tau))) stop("within_pt(tau): a numeric tolerance vector for target = \"realised\" (Prop. E35)")
  .att_trend("wpt", a_min = as.numeric(a_min), M = .att_nonneg(M, "within_pt(M)"), target = target, tau = tau)
}
#' @rdname within_pt
#' @export
pooled_pt <- function(M = 0) .att_trend("ppt", M = .att_nonneg(M, "pooled_pt(M)"))

#' LCE restrictions on expected rates for the mean model of eiATT()
#'
#' Linear contextual effects on the expected factual rates used by \code{within_pt()} and
#' \code{pooled_pt()} (2e Prop. E30): group \eqn{g}'s expected rate at the base election
#' (\code{base}), at each post-election (\code{post}) and its change (\code{change}) affine in the
#' composition \eqn{x_i}. Any two imply the third (E30(a)). \code{epsilon} is the baseline-curvature
#' allowance of E30(b): \eqn{|m_{g0}(x_i) - (a_g + b_g x_i)| \le \epsilon} (with \code{base = TRUE});
#' \eqn{\epsilon \ge 1/2} gives the change-only set.
#' @param base,post,change logical: which profiles are affine.
#' @param epsilon nonnegative baseline-curvature allowance.
#' @return an object of class \code{"edLCEmean"}.
#' @export
lce_mean <- function(base = FALSE, post = FALSE, change = TRUE, epsilon = 0) {
  structure(list(base = isTRUE(base), post = isTRUE(post), change = isTRUE(change),
                 epsilon = .att_nonneg(epsilon, "lce_mean(epsilon)")), class = "edLCEmean")
}

.wpt_lce_label <- function(l) {
  p <- c(if (l$base) sprintf("base%s", if (l$epsilon > 0) sprintf(" (curvature allowance epsilon = %g)", l$epsilon) else ""),
         if (l$post) "post", if (l$change) "change")
  if (!length(p)) "none" else paste(p, collapse = " + ")
}

# ---------------------------------------------------------------------------
# design: precinct quantities for the mean model and D_W
# ---------------------------------------------------------------------------
.wpt_design <- function(P, cs, tr, law) {
  r <- cs$r; cc <- cs$cc
  if (!is.null(cs$s) && P$R > 2)
    stop("within_pt()/pooled_pt(): Prop. E28 is stated for two groups; with R > 2 use comparison = \"rest\"")
  num <- function(v) as.double(v)
  N <- num(P$N0); aV <- num(P$A0[, r]); aN <- N - aV
  els <- c(cs$base, cs$posts)
  p <- P$p0; E <- length(els)
  Tm <- lo <- hi <- matrix(NA_real_, p, E)
  for (j in seq_len(E)) {
    e <- els[j]
    st <- law[[P$elections[e]]]
    if (is.null(st)) st <- matrix(num(P$Bt[[e]][, cc]), p, 1)
    st <- as.matrix(st)
    if (nrow(st) != p || any(!is.finite(st)) || any(st < -1e-9) || any(st > N + 1e-9))
      stop("law[[\"", P$elections[e], "\"]]: one row per unit (", p, ") of outcome counts between 0 and N_i")
    Tm[, j] <- rowMeans(st) / N
    lo[, j] <- rowMeans(matrix(pmax(0, st - aN), p)) / pmax(aV, 1e-300)     # per-state Duncan-Davis ends, averaged
    hi[, j] <- rowMeans(matrix(pmin(aV, st), p)) / pmax(aV, 1e-300)
  }
  sup <- aV > 0 & aN >= tr$a_min
  wV <- aV / sum(aV); wN <- if (sum(aN) > 0) aN / sum(aN) else rep(0, p)
  list(p = p, E = E, K = E - 1, els = els, N = N, aV = aV, aN = aN, x = aV / N, wV = wV, wN = wN,
       sup = if (tr$kind == "wpt") sup else rep(TRUE, p), Tm = Tm, lo = lo, hi = hi, units = P$units0,
       pi_U = if (tr$kind == "wpt") sum(wV[!sup]) else 0, n_eff = 1 / sum(wV^2))
}

# ---------------------------------------------------------------------------
# the mean-model LP
# ---------------------------------------------------------------------------
# returns a closure: solve(k, sense, mode = "range"/"eps"/"M", side) for post-period k
.wpt_lp <- function(des, kind, lce, M, Lambda) {
  p <- des$p; E <- des$E; K <- des$K
  nv <- 0; new <- function(n = 1) { i <- nv + seq_len(n); nv <<- nv + n; i }
  lb <- numeric(0); ub <- numeric(0)
  setb <- function(i, l, u) { lb[i] <<- l; ub[i] <<- u }
  mV <- matrix(NA_integer_, p, E); mN <- matrix(NA_integer_, p, E)
  for (i in seq_len(p)) for (j in seq_len(E)) {
    if (des$aV[i] > 0) { mV[i, j] <- new(); setb(mV[i, j], max(0, des$lo[i, j]), min(1, des$hi[i, j])) }
    if (des$aN[i] > 0) { mN[i, j] <- new(); setb(mN[i, j], 0, 1) }
  }
  iE <- new(); iM <- new()          # epsilon and M (fixed unless a breakdown LP frees one)
  setb(iE, lce$epsilon, lce$epsilon); setb(iM, M, M)
  rows <- list(); add <- function(idx, cf, dir, rhs) rows[[length(rows) + 1]] <<- list(idx = idx, cf = cf, dir = dir, rhs = rhs)
  # adding-up in each precinct and election
  for (i in seq_len(p)) for (j in seq_len(E)) {
    idx <- c(mV[i, j], mN[i, j]); cf <- c(des$x[i], 1 - des$x[i])[!is.na(idx)]; idx <- idx[!is.na(idx)]
    add(idx, cf, "==", des$Tm[i, j])
  }
  # LCE profiles (E30): explicit affine coefficients over the precincts where the group exists
  lcerow <- function(g, f, tol) {
    have <- if (g == "V") des$aV > 0 else des$aN > 0
    ab <- new(2); setb(ab, -Inf, Inf)
    for (i in which(have)) {
      fi <- f(i)
      idx <- c(fi$idx, ab); cf <- c(fi$cf, -1, -des$x[i])
      if (tol) { add(c(idx, iE), c(cf, -1), "<=", 0); add(c(idx, iE), c(cf, 1), ">=", 0) } else add(idx, cf, "==", 0)
    }
  }
  for (g in c("V", "N")) {
    m <- if (g == "V") mV else mN
    if (lce$base) lcerow(g, function(i) list(idx = m[i, 1], cf = 1), TRUE)
    for (k in seq_len(K)) {
      if (lce$post) lcerow(g, function(i) list(idx = m[i, k + 1], cf = 1), FALSE)
      if (lce$change) lcerow(g, function(i) list(idx = c(m[i, k + 1], m[i, 1]), cf = c(1, -1)), FALSE)
    }
  }
  obj <- vector("list", K)
  sv <- which(des$aV > 0 & des$sup); uv <- which(des$aV > 0 & !des$sup)
  if (kind == "wpt") {
    bprev <- rep(NA_integer_, p)
    for (k in seq_len(K)) {
      oi <- integer(0); oc <- numeric(0)
      for (i in sv) {
        b <- new(); S <- new(); setb(b, -Inf, Inf); setb(S, -Lambda, Lambda)
        # local dial (E33(a)): |b_ik - b_i,k-1 | <= M, b_i0 = 0
        if (is.na(bprev[i])) { add(c(b, iM), c(1, -1), "<=", 0); add(c(b, iM), c(1, 1), ">=", 0) }
        else { add(c(b, bprev[i], iM), c(1, -1, -1), "<=", 0); add(c(b, bprev[i], iM), c(1, -1, 1), ">=", 0) }
        bprev[i] <- b
        # rate limits (E28(b), (f)): m_iNt - S in [0,1]; m_iV0 + d_iN - S + b in [0,1]
        add(c(mN[i, k + 1], S), c(1, -1), ">=", 0); add(c(mN[i, k + 1], S), c(1, -1), "<=", 1)
        rl <- c(mV[i, 1], mN[i, k + 1], mN[i, 1], S, b); rc <- c(1, 1, -1, -1, 1)
        add(rl, rc, ">=", 0); add(rl, rc, "<=", 1)
        w <- des$wV[i]
        oi <- c(oi, mV[i, k + 1], mV[i, 1], mN[i, k + 1], mN[i, 1], S, b); oc <- c(oc, w, -w, -w, w, w, -w)
      }
      for (i in uv) {
        th <- new(); setb(th, 0, 1)
        oi <- c(oi, mV[i, k + 1], th); oc <- c(oc, des$wV[i], -des$wV[i])
      }
      obj[[k]] <- list(idx = oi, cf = oc)
    }
  } else {
    bprev <- NA_integer_
    vv <- which(des$aV > 0); nn <- which(des$aN > 0)
    for (k in seq_len(K)) {
      b <- new(); S <- new(); setb(b, -Inf, Inf); setb(S, -Lambda, Lambda)
      if (is.na(bprev)) { add(c(b, iM), c(1, -1), "<=", 0); add(c(b, iM), c(1, 1), ">=", 0) }
      else { add(c(b, bprev, iM), c(1, -1, -1), "<=", 0); add(c(b, bprev, iM), c(1, -1, 1), ">=", 0) }
      bprev <- b
      # aggregate rate limits only (E28(b)); the witness is one common untreated rate p0 (E28(c) step 3)
      add(c(mN[nn, k + 1], S), c(des$wN[nn], -1), ">=", 0); add(c(mN[nn, k + 1], S), c(des$wN[nn], -1), "<=", 1)
      rl <- c(mV[vv, 1], mN[nn, k + 1], mN[nn, 1], S, b); rc <- c(des$wV[vv], des$wN[nn], -des$wN[nn], -1, 1)
      add(rl, rc, ">=", 0); add(rl, rc, "<=", 1)
      obj[[k]] <- list(idx = c(mV[vv, k + 1], mV[vv, 1], mN[nn, k + 1], mN[nn, 1], S, b),
                       cf = c(des$wV[vv], -des$wV[vv], -des$wN[nn], des$wN[nn], 1, -1),
                       p0 = list(idx = rl[seq_len(length(rl) - 2)], cf = rc[seq_len(length(rc) - 2)], S = S, b = b))
    }
  }
  dense <- function(o) { v <- numeric(nv); for (q in seq_along(o$idx)) v[o$idx[q]] <- v[o$idx[q]] + o$cf[q]; v }
  solve <- function(k, sense, mode = "range", side = NULL) {
    rr <- rows; l <- lb; u <- ub; ob <- dense(obj[[k]])
    if (mode != "range") {
      free <- if (mode == "eps") iE else iM
      l[free] <- 0; u[free] <- Inf
      rr[[length(rr) + 1]] <- list(idx = obj[[k]]$idx, cf = obj[[k]]$cf, dir = if (side == "le") "<=" else ">=", rhs = 0)
      ob <- numeric(nv); ob[free] <- 1; sense <- -1
    }
    tri <- do.call(rbind, lapply(seq_along(rr), function(q) cbind(q, rr[[q]]$idx, rr[[q]]$cf)))
    mm <- slam::simple_triplet_matrix(as.integer(tri[, 1]), as.integer(tri[, 2]), as.double(tri[, 3]), nrow = length(rr), ncol = nv)
    fin_l <- is.finite(l); fin_u <- is.finite(u)
    bn <- list(lower = list(ind = seq_len(nv), val = ifelse(fin_l, l, -Inf)),
               upper = list(ind = which(fin_u), val = u[fin_u]))
    s <- Rglpk::Rglpk_solve_LP(ob, mm, vapply(rr, `[[`, "", "dir"), vapply(rr, `[[`, 0, "rhs"), bounds = bn, max = sense > 0,
                               control = list(canonicalize_status = FALSE, tm_limit = 60000L))
    if (s$status == 5) {
      out <- list(value = s$optimum, x = s$solution)
      if (kind == "ppt" && mode == "range") {
        pz <- obj[[k]]$p0; out$p0 <- sum(pz$cf * s$solution[pz$idx]) - s$solution[pz$S] + s$solution[pz$b]
      }
      return(out)
    }
    if (s$status %in% c(3, 4)) return(NULL)
    if (s$status == 6) return(list(value = if (sense > 0) Inf else -Inf))
    stop("the mean-model LP solver failed (GLPK status ", s$status, "); no result is reported rather than a false 'empty' set")
  }
  list(solve = solve, K = K, mV = mV, mN = mN)
}

# one-LP breakdown value of the sign conclusion for setting 'mode' ("eps" or "M"), the other held fixed:
# the smallest setting at which the interval reaches 0 (Inf if never; NA if it already contains 0 or is empty)
.wpt_breakdown <- function(lp, k, iv, mode) {
  if (anyNA(iv) || (iv[1] <= 0 && iv[2] >= 0)) return(NA_real_)
  s <- lp$solve(k, -1, mode, side = if (iv[1] > 0) "le" else "ge")
  if (is.null(s)) Inf else s$value
}

# E28(e): FD-LCE projection (two groups): supported part = I_W + pi_S dkappa (WPT), I_P + s dkappa (PPT)
.wpt_e28e <- function(des, k) {
  x <- des$x; dT <- des$Tm[, k + 1] - des$Tm[, 1]
  Z <- cbind(x^2, 2 * x * (1 - x), (1 - x)^2)
  if (qr(Z)$rank < 3) return(NULL)
  Q <- qr.solve(Z, dT)               # (dQvv, dQvn, dQnn); exact when change LCE holds for the law
  qd <- (Q[1] - Q[2]) * x + (Q[2] - Q[3]) * (1 - x)
  sv <- des$aV > 0 & des$sup
  list(Q = Q, I_W = sum(des$wV[sv] * qd[sv]), pi_S = 1 - des$pi_U,
       I_P = sum(des$wV * (Q[1] * x + Q[2] * (1 - x))) - sum(des$wN * (Q[2] * x + Q[3] * (1 - x))),
       s = sum(des$wV * (1 - x)) + sum(des$wN * x), resid = max(abs(Z %*% Q - dT)))
}

# E29: realised target D_W (count layer) and the Hoeffding allowance
.wpt_dw <- function(P, fb, cs, des, k, alpha, independence, roles, M = 0, Lambda = 0) {
  t <- cs$posts[k]; b0 <- cs$base; cc <- cs$cc; r <- cs$r
  nrows <- if (is.null(cs$s)) setdiff(seq_len(P$R), r) else cs$s
  terms <- list()
  for (i in which(des$aV > 0)) {
    w <- des$wV[i]
    if (des$sup[i]) {
      terms <- c(terms, list(.ed_term(P, t, i, r, cc, w / des$aV[i]), .ed_term(P, b0, i, r, cc, -w / des$aV[i]),
                             .ed_term(P, t, i, nrows, cc, -w / des$aN[i]), .ed_term(P, b0, i, nrows, cc, w / des$aN[i])))
    } else terms <- c(terms, list(.ed_term(P, t, i, r, cc, w / des$aV[i])))
  }
  LU <- .att_nrange(P, fb, .ec_obj(P, terms))
  sV <- des$aV > 0 & des$sup; uV <- des$aV > 0 & !des$sup
  rho2 <- if (independence == "precincts") sum((4 * des$wV[sV])^2) + sum(des$wV[uV]^2) else
    sum(4 * des$wV[sV]^2 * (1 / des$aV[sV] + 1 / des$aN[sV])) + sum(des$wV[uV]^2 / des$aV[uV])
  ra <- sqrt(0.5 * log(2 / alpha) * rho2)
  # E29(d): the local dial |b_it - b_i,t-1| <= M (so |b_ik| <= k M) and |S| <= Lambda restrict means, not
  # the realised [L, U]; they widen both ends by pi_S (h M + Lambda) at horizon h = k (code review R2-1)
  wd <- (1 - des$pi_U) * (k * M + Lambda)
  ci <- if (anyNA(LU)) c(NA_real_, NA_real_) else c(max(-1, LU[1] - des$pi_U - ra - wd), min(1, LU[2] + ra + wd))
  list(LU = LU, r_alpha = ra, ci = ci, belief = any(roles == "belief"), widen = wd, h = k)
}

.wpt_run <- function(fit, P, fb, cs, trend, Lambda, lce, law, alpha, independence, tgt) {
  if (!cs$has_comp) stop("internal: comparison required")
  des <- .wpt_design(P, cs, trend, law)
  lp <- .wpt_lp(des, trend$kind, lce, trend$M, Lambda)
  roles <- if (is.null(fit$restrictions)) character(0) else fit$restrictions$role
  res <- list()
  for (t in tgt) {
    k <- match(t, cs$posts)
    lo <- lp$solve(k, -1); hi <- if (is.null(lo)) NULL else lp$solve(k, 1)
    iv <- if (is.null(lo) || is.null(hi)) c(NA_real_, NA_real_) else c(lo$value, hi$value)
    e <- list(post = P$elections[t], lower = iv[1], upper = iv[2], empty = anyNA(iv),
              pieces = if (anyNA(iv)) NULL else matrix(iv, 1, dimnames = list(NULL, c("lower", "upper"))))
    e$breakdown_M <- .wpt_breakdown(lp, k, iv, "M")
    e$breakdown_eps <- if (lce$base) .wpt_breakdown(lp, k, iv, "eps") else NA_real_
    if (trend$kind == "ppt" && !is.null(lo)) e$p0 <- c(lo$p0, hi$p0)
    if (lce$change) e$e28e <- .wpt_e28e(des, k)
    if (trend$kind == "wpt") e$dw <- .wpt_dw(P, fb, cs, des, k, alpha, independence, roles, trend$M, Lambda)
    else e$D <- .att_nrange(P, fb, .att_rate(P, t, cs$r) - .att_rate_s(P, t, cs$r, cs$s) -
                                     (.att_rate(P, cs$base, cs$r) - .att_rate_s(P, cs$base, cs$r, cs$s)))
    res[[length(res) + 1]] <- e
  }
  list(results = res, design = des, lce = lce, law_given = !is.null(law), independence = independence, alpha = alpha,
       support = data.frame(unit = des$units, a_treated = des$aV, a_comparison = des$aN, w_treated = des$wV,
                            supported = des$sup, stringsAsFactors = FALSE))
}

.wpt_print <- function(x, digits = 3) {
  w <- x$wpt; des <- w$design
  f <- function(v) if (is.null(v) || anyNA(v)) "EMPTY" else sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), v[1], v[2])
  wpt <- x$trend$kind == "wpt"
  cat("  estimand: ATT^mean_t = sum_i w_Vi (m_iVt - m^(0)_iVt), the expected programme effect on the same treated persons (2e Sec. 5.8)\n")
  cat(sprintf("  mean model: expected rates %s; LCE on expected rates: %s (Prop. E30)\n",
              if (w$law_given) "from the supplied law (known law, E28(c))" else
                "taken to be the observed returns (known-law reading; with one realised election pair the law is not known: E29 is the confidence route)",
              .wpt_lce_label(w$lce)))
  if (wpt) {
    cat(sprintf("  support: a_iN >= a_min = %g; pi_U = %.4g (treated share in unsupported precincts, bounded with theta_i in [0,1]); n_eff = 1/sum w_Vi^2 = %.4g\n",
                x$trend$a_min, des$pi_U, des$n_eff))
    sp <- w$support[w$support$a_treated > 0, ]
    cat("  local comparison counts a_iN (treated precincts): ",
        paste(sprintf("%s %s%s", sp$unit, format(sp$a_comparison, trim = TRUE), ifelse(sp$supported, "", " (unsupported)")), collapse = ", "), "\n", sep = "")
    cat(sprintf("  E29 independence unit: %s\n", if (w$independence == "precincts")
      "precincts (independent across precincts given the design; any dependence within a precinct)" else
        "persons (ballot pairs independent across persons; precincts with few comparison voters enter through 1/a_iN)"))
  }
  cat("  count rows are not used by the mean model (they restrict realised counts); they enter the realised target D_W\n")
  for (e in w$results) {
    cat(sprintf("\nPost-election %s\n", e$post))
    cat(sprintf("  [IDENTIFICATION BOUND, mean model (E28/E30)] ATT^mean: %s%s\n", f(c(e$lower, e$upper)),
                if (e$empty) sprintf(" -- refutes the mean model %s, with %s and its rate limits", if (w$law_given) "under the supplied law" else "at the observed returns taken as the law", if (wpt) "WPT" else "PPT") else ""))
    if (!e$empty) {
      bd <- function(v, nm) if (is.na(v)) "none (the interval contains 0)" else if (!is.finite(v)) "infinite (no setting reaches 0)" else
        sprintf("%s* = %.4g (the sign holds strictly below %s*; at %s* zero first becomes feasible)", nm, v, nm, nm)
      cat(sprintf("    breakdown of the sign, %s dial M (other settings fixed; one LP): %s\n", if (wpt) "local" else "aggregate", bd(e$breakdown_M, "M")))
      if (w$lce$base) cat(sprintf("    breakdown of the sign, baseline curvature epsilon (M fixed; one LP): %s\n", bd(e$breakdown_eps, "epsilon")))
    }
    if (!is.null(e$e28e)) {
      z <- e$e28e
      if (wpt) cat(sprintf("    E28(e) FD-LCE projection: supported part = I_W + pi_S dkappa, I_W = %.4g, pi_S = %.4g%s\n",
                           z$I_W, z$pi_S, if (z$resid > 1e-8) sprintf(" (the law's changes are not exactly quadratic: max residual %.2g)", z$resid) else ""))
      else cat(sprintf("    E28(e) FD-LCE projection: ATT = I + s dkappa, I = %.4g, s = %.4g\n", z$I_P, z$s))
    }
    if (!is.null(e$p0)) cat(sprintf("    PPT witness: one common untreated rate p0 in every precinct with treated persons (%.4g at the lower end, %.4g at the upper end)\n", e$p0[1], e$p0[2]))
    if (!is.null(e$dw)) {
      d <- e$dw
      cat(sprintf("  [DESCRIPTIVE, realised counts] D_W (E29), count layer with the fit's rows: %s\n", f(d$LU)))
      cat(sprintf("  [CONFIDENCE INTERVAL, one observed election pair (E29(c))] ATT^mean: %s at level %g (r_alpha = %.4g; [L - pi_U - r - w, U + r + w], clipped to [-1, 1])\n",
                  f(d$ci), 1 - w$alpha, d$r_alpha))
      cat(sprintf("    E29(d) widening w = pi_S (h M + Lambda) = %.4g x (%d x %g + %g) = %.4g (local dial and spillover restrict means, not the realised [L, U])\n",
                  1 - des$pi_U, d$h, x$trend$M, x$Lambda, d$widen))
      cat(sprintf("    coverage >= %s; E29 covers the realised-count route only: it does not make the identification bound above a confidence interval\n",
                  if (d$belief) sprintf("1 - %g - P(the belief rows are false) (belief rows in the count layer, E13)", w$alpha) else sprintf("1 - %g (identity rows only)", w$alpha)))
    }
    if (!is.null(e$D)) cat(sprintf("  [DESCRIPTIVE, realised counts] pooled DiD D_t: %s (no confidence interval for PPT is implemented; E29 is for WPT)\n", f(e$D)))
  }
  cat(sprintf("\n  spillover: %s\n", if (x$Lambda == 0) "zero mean spillover" else sprintf("|S| <= %g, optimised jointly (E28(f))", x$Lambda)))
}

# ---------------------------------------------------------------------------
# Prop. E35: realised within-precinct parallel trends with tolerance tau_i (RWPT-tau), the default
# reading of within_pt(). Rows in Theorem E16's joint LP over (n, theta_i, S_i, M), per post-period k:
#   |(theta_ik - B_iV,base(n)) - (B_iN,t(n) - S_ik - B_iN,base(n))| <= c_i M   (i supported),
#   B_iN,t(n) - S_ik in [0, 1], |S_ik| <= Lambda (supported), theta_ik in [0, 1] (every precinct with
#   treated persons; free on unsupported ones). c_i = sqrt(1/a_iV + 1/a_iN), or the user's tau_i with M = 1.
# Target: ATT_t = sum_i w_Vi (B_iV,t(n) - theta_i). LCE only through E13: a row D_W^S(n) in CI.
# ---------------------------------------------------------------------------
.rwpt_lp <- function(P, cs, des, tau) {
  p <- des$p; K <- length(cs$posts); r <- cs$r; cc <- cs$cc; b0 <- cs$base
  nrows <- if (is.null(cs$s)) setdiff(seq_len(P$R), r) else cs$s
  nx <- 2 * p * K + 1; iM <- nx
  ith <- function(k, i) (k - 1) * 2 * p + i; iS <- function(k, i) (k - 1) * 2 * p + p + i
  BV <- function(t, i) .ec_obj(P, list(.ed_term(P, t, i, r, cc, 1 / des$aV[i])))
  BN <- function(t, i) .ec_obj(P, list(.ed_term(P, t, i, nrows, cc, 1 / des$aN[i])))
  cf <- if (is.null(tau)) sqrt(1 / des$aV + 1 / pmax(des$aN, 1e-300)) else tau
  xv <- function(...) { v <- numeric(nx); a <- list(...); for (q in a) v[q[1]] <- v[q[1]] + q[2]; v }
  lower <- numeric(nx); upper <- numeric(nx)
  rows <- list()
  for (k in seq_len(K)) {
    t <- cs$posts[k]
    for (i in seq_len(p)) {
      if (des$aV[i] > 0) upper[ith(k, i)] <- 1
      if (!(des$sup[i] && des$aV[i] > 0)) next
      lower[iS(k, i)] <- -des$Lambda; upper[iS(k, i)] <- des$Lambda
      gn <- -BV(b0, i) - BN(t, i) + BN(b0, i)                        # gap = theta + S + gn(n)
      rows[[length(rows) + 1]] <- list(n = gn, x = xv(c(ith(k, i), 1), c(iS(k, i), 1), c(iM, -cf[i])), lo = -Inf, hi = 0)
      rows[[length(rows) + 1]] <- list(n = gn, x = xv(c(ith(k, i), 1), c(iS(k, i), 1), c(iM, cf[i])), lo = 0, hi = Inf)
      rows[[length(rows) + 1]] <- list(n = BN(t, i), x = xv(c(iS(k, i), -1)), lo = 0, hi = 1)
    }
  }
  att <- function(k) {
    t <- cs$posts[k]; on <- numeric(P$nv); ox <- numeric(nx)
    for (i in which(des$aV > 0)) { on <- on + des$wV[i] * BV(t, i); ox[ith(k, i)] <- -des$wV[i] }
    list(n = on, x = ox)
  }
  list(nx = nx, iM = iM, rows = rows, lower = lower, upper = upper, att = att, cf = cf)
}

# FD-LCE interval for the realised supported target D_W^S (E13 with E12's allowance, rho_i = 4 w_Vi)
.rwpt_lce_ci <- function(fit, P, cs, des, k, alpha) {
  m <- fit$model
  if (is.null(m)) stop("within_pt(): lce = TRUE needs an eiChange fit with a model layer (model = fd() or lce_periods()); LCE enters only by E13")
  t <- cs$posts[k]; b0 <- cs$base; r <- cs$r; cc <- cs$cc
  nrows <- if (is.null(cs$s)) setdiff(seq_len(P$R), r) else cs$s
  sv <- which(des$sup & des$aV > 0)
  W <- matrix(0, P$p0, P$R); terms <- list()
  for (i in sv) {
    w <- des$wV[i]; W[i, r] <- w
    for (q in nrows) W[i, q] <- -w * as.double(P$A0[i, q]) / des$aN[i]
    terms <- c(terms, list(.ed_term(P, t, i, r, cc, w / des$aV[i]), .ed_term(P, b0, i, r, cc, -w / des$aV[i]),
                           .ed_term(P, t, i, nrows, cc, -w / des$aN[i]), .ed_term(P, b0, i, nrows, cc, w / des$aN[i])))
  }
  mb <- m$mb; spec <- m$spec
  f <- .ed_model_change(mb, b0, t, cc, W)
  tgm <- list(cvec = f$coef, c0 = f$const, gamma = f$a)
  a_est <- 0.9 * alpha; a_E <- 0.1 * alpha
  al <- .ed_allowance(P, list(list(terms = terms)), spec$allowance, a_E)
  Zs <- .ed_with_seed(spec$seed, .rxc_zdraws(mb$V, 4000))
  ivs <- lapply(.ed_mb_systems(m$base, m$unions), function(sy)
    .ed_with_seed(spec$seed, .ed_proj(mb, sy$cons, tgm, a_est, spec$ci, Zs, spec$ndraw)) + c(-1, 1) * al$value)
  list(pieces = .rxc_union(ivs, tol = 0), allowance = al$value, rho = al$rho, obj = .ec_obj(P, terms), rung = mb$rung)
}

.rwpt_run <- function(fit, P, fb, cs, trend, Lambda, lce, alpha, tgt) {
  if (!cs$has_comp) stop("internal: comparison required")
  des <- .wpt_design(P, cs, trend, NULL); des$Lambda <- Lambda
  tau <- trend$tau
  if (!is.null(tau)) {
    if (length(tau) == 1) tau <- rep(tau, des$p)
    if (length(tau) != des$p || any(!is.finite(tau)) || any(tau < 0)) stop("within_pt(tau): one nonnegative number per unit (", des$p, ") or one number")
  }
  lp <- .rwpt_lp(P, cs, des, tau)
  Mfix <- if (is.null(tau)) trend$M else 1
  roles <- if (is.null(fit$restrictions)) character(0) else fit$restrictions$role
  res <- list()
  for (t in tgt) {
    k <- match(t, cs$posts); a <- lp$att(k)
    lc <- if (isTRUE(lce)) .rwpt_lce_ci(fit, P, cs, des, k, alpha) else NULL
    lrows <- if (is.null(lc)) list(list()) else if (is.null(lc$pieces)) list() else
      lapply(seq_len(nrow(lc$pieces)), function(j) list(list(n = lc$obj, x = numeric(lp$nx), lo = lc$pieces[j, 1], hi = lc$pieces[j, 2])))
    solve_all <- function(Mlo, Mhi, extra, obj_n, obj_x, sense) {
      best <- NULL
      for (lr in lrows) for (b in fb) {
        lo <- lp$lower; hi <- lp$upper; lo[lp$iM] <- Mlo; hi[lp$iM] <- Mhi
        s <- .att_lp(P, b$rows, c(lp$rows, lr, extra), lp$nx, list(lower = lo, upper = hi), obj_n, obj_x, sense)
        if (!is.null(s) && (is.null(best) || sense * s$value > sense * best)) best <- s$value
      }
      best
    }
    ivs <- list()
    for (lr in lrows) for (b in fb) {
      lo <- lp$lower; hi <- lp$upper; lo[lp$iM] <- Mfix; hi[lp$iM] <- Mfix
      l <- .att_lp(P, b$rows, c(lp$rows, lr), lp$nx, list(lower = lo, upper = hi), a$n, a$x, -1)
      h <- if (is.null(l)) NULL else .att_lp(P, b$rows, c(lp$rows, lr), lp$nx, list(lower = lo, upper = hi), a$n, a$x, 1)
      ivs[[length(ivs) + 1]] <- if (is.null(l) || is.null(h)) c(NA_real_, NA_real_) else c(l$value, h$value)
    }
    pcs <- .rxc_union(ivs, tol = 1e-10)
    e <- list(post = P$elections[t], pieces = pcs, lower = if (is.null(pcs)) NA_real_ else min(pcs[, 1]),
              upper = if (is.null(pcs)) NA_real_ else max(pcs[, 2]), empty = is.null(pcs), lce = lc)
    # E35(b): breakdown of the sign, one LP with M a variable (the ATT row on the other side of 0)
    e$breakdown_M <- NA_real_
    if (!e$empty && (e$lower > 0 || e$upper < 0)) {
      side <- if (e$lower > 0) list(n = a$n, x = a$x, lo = -Inf, hi = 0) else list(n = a$n, x = a$x, lo = 0, hi = Inf)
      ox <- numeric(lp$nx); ox[lp$iM] <- 1
      v <- solve_all(0, Inf, list(side), numeric(P$nv), ox, -1)
      e$breakdown_M <- if (is.null(v)) Inf else v
    }
    e$dw <- .wpt_dw(P, fb, cs, des, k, alpha, "precincts", roles)$LU
    res[[length(res) + 1]] <- e
  }
  list(results = res, design = des, tau_given = !is.null(tau), M = trend$M, lce = isTRUE(lce), alpha = alpha,
       belief = any(roles == "belief"),
       support = data.frame(unit = des$units, a_treated = des$aV, a_comparison = des$aN, w_treated = des$wV,
                            supported = des$sup, stringsAsFactors = FALSE))
}

.rwpt_print <- function(x, digits = 3) {
  w <- x$rwpt; des <- w$design
  f <- function(v) if (is.null(v) || anyNA(v)) "EMPTY" else sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), v[1], v[2])
  cat(sprintf("  tolerance: %s\n", if (w$tau_given) "user tau_i (M* below scales the whole vector)" else
    if (w$M == 0) "tau_i = 0 (M = 0: exact realised WPT, Prop. E35(c))" else sprintf("tau_i = M sqrt(1/a_iV + 1/a_iN) with M = %g", w$M)))
  cat("  note (R4-1): M is a tolerance that must hold in every supported precinct at once, not a confidence level\n")
  cat(sprintf("  support: a_iN >= a_min = %g; pi_U = %.4g (treated share in unsupported precincts, theta_i in [0,1], no assumption); n_eff = %.4g\n",
              x$trend$a_min, des$pi_U, des$n_eff))
  sp <- w$support[w$support$a_treated > 0, ]
  cat("  local comparison counts a_iN (treated precincts): ",
      paste(sprintf("%s %s%s", sp$unit, format(sp$a_comparison, trim = TRUE), ifelse(sp$supported, "", " (unsupported)")), collapse = ", "), "\n", sep = "")
  cat("  continuous tables (Prop. E35(c): the integer check is optional and not applied; whole-person ends may be narrower, Remark 7.1)\n")
  for (e in w$results) {
    cat(sprintf("\nPost-election %s\n", e$post))
    lab <- if (w$lce) "CONFIDENCE SET, realised ATT (E35 + E13)" else "IDENTIFICATION SET, realised ATT, count layer (E35)"
    cat(sprintf("  [%s] ATT: %s%s\n", lab, if (e$empty) "EMPTY" else paste(apply(e$pieces, 1, f), collapse = " U "),
                if (e$empty) paste0(" -- refutes (FP) + NA + count rows + RWPT-tau + the spillover limit",
                                     if (w$lce) " + the E13 LCE row D_W^S(n) in CI (with probability at most alpha when the other assumptions hold)" else "") else ""))
    if (!e$empty) cat(sprintf("    breakdown of the sign (one LP, M a variable): %s\n",
      if (is.na(e$breakdown_M)) "none (the set contains 0)" else if (!is.finite(e$breakdown_M)) "infinite (no table puts the ATT on the other side of 0)" else
        sprintf("M* = %.4g", e$breakdown_M)))
    cat(sprintf("  [DESCRIPTIVE, realised counts] D_W (E29), count layer: %s\n", f(e$dw)))
    if (w$lce) {
      l <- e$lce
      cat(sprintf("    E13 row: D_W^S(n) in the %s interval %s at level %g (E12 allowance %.4g, rho_i = 4 w_Vi)\n", l$rung,
                  if (is.null(l$pieces)) "EMPTY" else paste(apply(l$pieces, 1, f), collapse = " U "), 1 - w$alpha, l$allowance))
      cat(sprintf("    coverage >= %s under FP, NA, RWPT-tau, |S_i| <= Lambda%s; the confidence level comes only from the LCE part\n",
                  if (w$belief) sprintf("1 - %g - P(the belief rows are false)", w$alpha) else sprintf("1 - %g", w$alpha),
                  if (w$belief) " and the belief rows" else " (identity rows only)"))
    }
  }
  cat(sprintf("\n  spillover: |S_i| <= %g per supported precinct%s\n", x$Lambda, if (x$Lambda == 0) " (no spillover, stated explicitly)" else ""))
}
