# Changes and DiD, causal layer: eiATT() on an eiChange() fit.
# Theory: docs/extension/PHASE2e_did_general.md (2e) Sec. 5 (Lemma E14, Def. 5.1, Props E15-E22,
# Thm E16, Cor. E17, Defs 5.2-5.3) and Sec. 6 (Prop. E25 for model-based pre-trends).
#
# Joint causal LP (Thm E16): variables the stacked realised table n (count layer) and, for each
# post-period t, theta_t (the treated group's untreated rate), S_t (spillover onto the comparison
# group's rate) and b_t (counterfactual trend gap), linked by Lemma E14:
#   ATT_t = B_{r*ct}(n) - theta_t,   b_t = theta_t - thetaPT_t(n) + S_t,
#   thetaPT_t(n) = B_{r*c0} + B_{s*ct} - B_{s*c0}     (so ATT_t = D_t - b_t + S_t).
# Rate limits: theta_t in [0, 1], B_{s*ct}(n) - S_t in [0, 1]; |S_t| <= Lambda.

# ---------------------------------------------------------------------------
# trend assumptions (Def. 5.1 item 2)
# ---------------------------------------------------------------------------

#' Trend assumptions for eiATT()
#'
#' Constructors for the trend assumption of \code{\link{eiATT}} (2e Def. 5.1 item 2). Each
#' restricts the counterfactual trend gap \eqn{b_t} (the difference between the treated and
#' comparison groups' untreated changes from the base election to \eqn{t}).
#' \itemize{
#'   \item \code{no_trend()}: no trend assumption: Prop. E15, the set is
#'     \eqn{[B^L_{r^*ct}-1, B^U_{r^*ct}]}, nothing causal follows.
#'   \item \code{trend_dial(M)}: \eqn{|b_t| \le M} (Cor. E17); \code{parallel_trends()} is
#'     \code{trend_dial(0)}.
#'   \item \code{rel_mag(Mbar, pre)}: Rambachan-Roth relative magnitudes (Def. 5.2):
#'     \eqn{|\delta_t-\delta_{t-1}| \le \bar M A(n)}, \eqn{A(n)} the largest consecutive pre-period
#'     DiD in absolute value, over the \emph{joint} set of latent tables (Prop. E18).
#'   \item \code{smoothness(M)}: second differences of \eqn{\delta} at most \eqn{M} (Def. 5.3).
#'   \item \code{mtr_nd_ps(mtr, nd, ps)}: monotone treatment response, no decline of the
#'     untreated outcome, positive selection on the untreated comparison rate (Prop. E21).
#' }
#' @param M,Mbar nonnegative numbers.
#' @param pre pre-period elections used by \code{rel_mag} (labels or positions; contiguous,
#'   ending just before the base); \code{NULL} for all.
#' @param mtr,nd,ps logical: which of the three assumptions to impose.
#' @return an object of class \code{"edTrend"}.
#' @name trend_assumptions
NULL

.att_trend <- function(kind, ...) structure(list(kind = kind, ...), class = "edTrend")
.att_nonneg <- function(x, what) {
  if (length(x) != 1 || !is.numeric(x) || !is.finite(x) || x < 0) stop(what, " must be one nonnegative number")
  as.numeric(x)
}
#' @rdname trend_assumptions
#' @export
no_trend <- function() .att_trend("none")
#' @rdname trend_assumptions
#' @export
trend_dial <- function(M) .att_trend("dial", M = .att_nonneg(M, "trend_dial(M)"))
#' @rdname trend_assumptions
#' @export
parallel_trends <- function() .att_trend("dial", M = 0)
#' @rdname trend_assumptions
#' @export
rel_mag <- function(Mbar, pre = NULL) .att_trend("rm", Mbar = .att_nonneg(Mbar, "rel_mag(Mbar)"), pre = pre)
#' @rdname trend_assumptions
#' @export
smoothness <- function(M) .att_trend("sd", M = .att_nonneg(M, "smoothness(M)"))
#' @rdname trend_assumptions
#' @export
mtr_nd_ps <- function(mtr = TRUE, nd = TRUE, ps = TRUE) {
  if (!any(c(mtr, nd, ps))) stop("mtr_nd_ps(): impose at least one of mtr, nd, ps")
  .att_trend("mtr", mtr = isTRUE(mtr), nd = isTRUE(nd), ps = isTRUE(ps))
}

.att_trend_label <- function(tr) switch(tr$kind,
  none = "none (Prop. E15)",
  dial = if (tr$M == 0) "parallel trends, |b_t| = 0 (Cor. E17 with M = 0)" else sprintf("trend dial |b_t| <= %g (Cor. E17)", tr$M),
  rm = sprintf("relative magnitudes, Mbar = %g (Def. 5.2, Prop. E18)", tr$Mbar),
  sd = sprintf("smoothness, M = %g (Def. 5.3)", tr$M),
  mtr = paste0(paste(c(if (tr$mtr) "MTR", if (tr$nd) "ND", if (tr$ps) "PS"), collapse = " + "), " (Prop. E21)"),
  wpt = sprintf("within-precinct parallel trends in expectation (WPT, Prop. E28)%s, support a_iN >= %g",
                if (tr$M > 0) sprintf(" relaxed by the local dial |b_it - b_i,t-1| <= %g (Prop. E33(a))", tr$M) else "", tr$a_min),
  ppt = sprintf("pooled parallel trends in expectation (PPT, Prop. E28)%s",
                if (tr$M > 0) sprintf(" relaxed by the aggregate dial |b_t - b_t-1| <= %g", tr$M) else ""))
.att_uses_comparison <- function(tr) tr$kind %in% c("dial", "rm", "sd", "wpt", "ppt") || (tr$kind == "mtr" && tr$ps)

# ---------------------------------------------------------------------------
# linear forms on the stacked table
# ---------------------------------------------------------------------------
.att_rate <- function(P, t, r) .ec_obj(P, list(.ed_term(P, t, seq_len(P$p0), r, P$cc, 1 / sum(P$A0[, r]))))
.att_rate_s <- function(P, t, r, s) {
  if (!is.null(s)) return(.att_rate(P, t, s))
  rest <- setdiff(seq_len(P$R), r)
  .ec_obj(P, list(.ed_term(P, t, seq_len(P$p0), rest, P$cc, 1 / (sum(P$N0) - sum(P$A0[, r])))))
}

# one LP: optimise obj_n' n + obj_x' x over margins, count rows and extra rows.
# extra row: list(n = numeric(nv) or NULL, x = numeric(nx), lo, hi). xb: list(lower, upper).
.att_lp <- function(P, rows_n, extra, nx, xb, obj_n, obj_x, sense) {
  m <- P$margin_all
  ex <- .ec_rows_trip(rows_n, m$n, rd = FALSE)
  trip <- rbind(m$trip, ex$trip); dir <- c(rep("==", m$n), sub("^=$", "==", ex$dir)); rhs <- c(m$rhs, ex$rhs)
  j <- if (is.null(ex$trip)) m$n else ex$n
  tl <- list(trip)
  for (e in extra) {
    idx <- c(if (!is.null(e$n)) which(e$n != 0), P$nv + which(e$x != 0))
    cf <- c(if (!is.null(e$n)) e$n[e$n != 0], e$x[e$x != 0])
    if (!length(idx)) next
    add <- function(d, v) { j <<- j + 1; tl[[length(tl) + 1]] <<- cbind(j, idx, cf); dir <<- c(dir, d); rhs <<- c(rhs, v) }
    if (is.finite(e$lo) && is.finite(e$hi) && e$lo == e$hi) add("==", e$lo) else {
      if (is.finite(e$lo)) add(">=", e$lo)
      if (is.finite(e$hi)) add("<=", e$hi)
    }
  }
  trip <- do.call(rbind, tl)
  nvar <- P$nv + nx
  mm <- slam::simple_triplet_matrix(as.integer(trip[, 1]), as.integer(trip[, 2]), as.double(trip[, 3]),
                                    nrow = j, ncol = nvar)
  bnds <- NULL
  if (nx > 0) {
    fin <- is.finite(xb$upper)
    bnds <- list(lower = list(ind = P$nv + seq_len(nx), val = xb$lower))
    if (any(fin)) bnds$upper <- list(ind = P$nv + which(fin), val = xb$upper[fin])
  }
  # GLPK statuses as in .ec_optimise (code review issue 1): 5 optimal; 3/4 infeasible (an empty
  # branch); 6 unbounded (the end is +/-Inf); anything else is a solver failure, never "empty"
  s <- Rglpk::Rglpk_solve_LP(c(obj_n, obj_x), mm, dir, rhs, bounds = bnds, max = sense > 0,
                             control = list(canonicalize_status = FALSE, tm_limit = 60000L))
  if (s$status == 5) return(list(value = s$optimum, x = s$solution))
  if (s$status %in% c(3, 4)) return(NULL)
  if (s$status == 6) return(list(value = if (sense > 0) Inf else -Inf, x = NULL))
  stop("the causal LP solver failed (GLPK status ", s$status, "); no result is reported rather than a false 'empty' set")
}

# F branches of the fit: common rows plus one side of each CR union
.att_fbranches <- function(fit) {
  b <- fit$built; U <- length(b$unions)
  combos <- if (U) as.matrix(expand.grid(rep(list(1:2), U))) else matrix(0L, 1, 0)
  lapply(seq_len(nrow(combos)), function(k) {
    rows <- b$common; lab <- character(0)
    for (u in seq_len(U)) {
      side <- if (combos[k, u] == 1) "plus" else "minus"
      rows <- c(rows, b$unions[[u]][[side]]); lab <- c(lab, paste0(b$unions[[u]]$id, if (side == "plus") "+" else "-"))
    }
    list(rows = rows, label = paste(lab, collapse = ", "))
  })
}

# range of a linear form in n over F (union over branches): c(lower, upper) or NA
.att_nrange <- function(P, fb, obj, extra_rows = list()) {
  lo <- Inf; hi <- -Inf
  for (b in fb) {
    l <- .att_lp(P, b$rows, extra_rows, 0, list(lower = numeric(0), upper = numeric(0)), obj, numeric(0), -1)
    if (is.null(l)) next
    h <- .att_lp(P, b$rows, extra_rows, 0, list(lower = numeric(0), upper = numeric(0)), obj, numeric(0), 1)
    lo <- min(lo, l$value); hi <- max(hi, h$value)
  }
  if (!is.finite(lo)) c(NA_real_, NA_real_) else c(lo, hi)
}

# ---------------------------------------------------------------------------
# the causal system for one trend assumption: Gamma branches (each a list of extra rows)
# ---------------------------------------------------------------------------
# x layout: for post-period k (k = 1..K): theta = 3k - 2, S = 3k - 1, b = 3k
.att_system <- function(P, cs, tr, Lambda, rate_limits = TRUE) {
  K <- length(cs$posts); nx <- 3 * K
  ix <- function(k, w) 3 * (k - 1) + switch(w, theta = 1, S = 2, b = 3)
  ex <- function(...) { v <- numeric(nx); a <- list(...); for (nm in names(a)) for (k in seq_along(a[[nm]])) NULL; v }
  xv <- function(pairs) { v <- numeric(nx); for (p in pairs) v[p[1]] <- v[p[1]] + p[2]; v }
  r <- cs$r; s <- cs$s; b0 <- cs$base
  Br <- function(t) .att_rate(P, t, r); Bs <- function(t) .att_rate_s(P, t, r, s)
  g <- function(t) Br(t) - Bs(t)
  thPT <- function(t) Br(b0) + Bs(t) - Bs(b0)
  lower <- rep(-Inf, nx); upper <- rep(Inf, nx)
  common <- list()
  for (k in seq_len(K)) {
    t <- cs$posts[k]
    lower[ix(k, "theta")] <- if (rate_limits) 0 else -Inf; upper[ix(k, "theta")] <- if (rate_limits) 1 else Inf
    lam <- if (cs$has_comp) Lambda else 0
    lower[ix(k, "S")] <- -lam; upper[ix(k, "S")] <- lam
    # b_t - theta_t - S_t + thetaPT_t(n) = 0  (Lemma E14)
    if (cs$has_comp)
      common[[length(common) + 1]] <- list(n = thPT(t), x = xv(list(c(ix(k, "b"), 1), c(ix(k, "theta"), -1), c(ix(k, "S"), -1))), lo = 0, hi = 0)
    if (cs$has_comp && rate_limits)
      common[[length(common) + 1]] <- list(n = Bs(t), x = xv(list(c(ix(k, "S"), -1))), lo = 0, hi = 1)
  }
  delta <- function(k) if (k == 0) list(n = NULL, x = numeric(nx)) else list(n = NULL, x = xv(list(c(ix(k, "b"), 1))))
  dif <- function(a, b) list(n = if (is.null(a$n) && is.null(b$n)) NULL else (if (is.null(a$n)) 0 else a$n) - (if (is.null(b$n)) 0 else b$n),
                             x = a$x - b$x)
  branches <- list()
  if (tr$kind == "none") branches <- list(list(rows = list(), label = ""))
  if (tr$kind == "dial") {
    for (k in seq_len(K)) { lower[ix(k, "b")] <- -tr$M; upper[ix(k, "b")] <- tr$M }
    branches <- list(list(rows = list(), label = ""))
  }
  if (tr$kind == "mtr") {
    rr <- list()
    for (k in seq_len(K)) {
      t <- cs$posts[k]; th <- xv(list(c(ix(k, "theta"), 1)))
      if (tr$mtr) rr[[length(rr) + 1]] <- list(n = -Br(t), x = th, lo = -Inf, hi = 0)          # theta <= B_r,t
      if (tr$nd) rr[[length(rr) + 1]] <- list(n = -Br(b0), x = th, lo = 0, hi = Inf)           # theta >= B_r,0
      if (tr$ps) rr[[length(rr) + 1]] <- list(n = -Bs(t), x = th + xv(list(c(ix(k, "S"), 1))), lo = 0, hi = Inf)  # theta >= B_s,t - S
    }
    branches <- list(list(rows = rr, label = ""))
  }
  if (tr$kind == "rm") {
    pp <- cs$pre_pairs
    if (!length(pp)) stop("rel_mag() needs at least one pre-period before the base election (2e Sec. 1.1: L >= 1)")
    for (q in seq_along(pp)) for (sg in c(1, -1)) {
      d <- sg * tr$Mbar * (g(pp[[q]][2]) - g(pp[[q]][1]))
      rr <- list()
      for (k in seq_len(K)) {
        dd <- dif(delta(k), delta(k - 1))
        rr[[length(rr) + 1]] <- list(n = d, x = -dd$x, lo = 0, hi = Inf)    # Mbar sg Dpre - (delta_t - delta_t-1) >= 0
        rr[[length(rr) + 1]] <- list(n = d, x = dd$x, lo = 0, hi = Inf)
      }
      branches[[length(branches) + 1]] <- list(rows = rr, label = sprintf("A(n) at %s -> %s, sign %s",
                                                                          P$elections[pp[[q]][1]], P$elections[pp[[q]][2]], if (sg > 0) "+" else "-"))
    }
  }
  if (tr$kind == "sd") {
    # delta over all elections in order: pre: D_s(n) = g_s - g_0; base: 0; post: b_t
    els <- seq_len(P$T)
    dl <- lapply(els, function(e) {
      if (e < b0) list(n = g(e) - g(b0), x = numeric(nx)) else if (e == b0) list(n = NULL, x = numeric(nx))
      else { k <- match(e, cs$posts); if (is.na(k)) NULL else list(n = NULL, x = xv(list(c(ix(k, "b"), 1)))) }
    })
    rr <- list()
    for (j in seq_len(max(0, length(els) - 2)) + 1) {     # code review issue 8: no triple with T = 2
      if (is.null(dl[[j - 1]]) || is.null(dl[[j]]) || is.null(dl[[j + 1]])) next
      sdif <- dif(dif(dl[[j + 1]], dl[[j]]), dif(dl[[j]], dl[[j - 1]]))
      rr[[length(rr) + 1]] <- list(n = sdif$n, x = sdif$x, lo = -tr$M, hi = tr$M)
    }
    branches <- list(list(rows = rr, label = ""))
  }
  list(nx = nx, ix = ix, common = common, branches = branches, lower = lower, upper = upper,
       Br = Br, Bs = Bs, g = g, thPT = thPT)
}

# ATT range for post-period k: union over Gamma branches x F branches (Thm E16(c)).
# Branch intervals are merged only where they overlap; holes are kept (Prop. E18(a), finding 12).
.att_range <- function(P, fb, sys, k, t, cs = NULL) {
  obj_x <- numeric(sys$nx); obj_x[sys$ix(k, "theta")] <- -1
  obj_n <- sys$Br(t)
  ivs <- list(); tab <- list(); sols <- list()
  for (gb in sys$branches) for (b in fb) {
    xb <- list(lower = sys$lower, upper = sys$upper)
    ext <- c(sys$common, gb$rows)
    lo <- .att_lp(P, b$rows, ext, sys$nx, xb, obj_n, obj_x, -1)
    hi <- if (is.null(lo)) NULL else .att_lp(P, b$rows, ext, sys$nx, xb, obj_n, obj_x, 1)
    iv <- if (is.null(lo) || is.null(hi)) c(NA_real_, NA_real_) else c(lo$value, hi$value)
    ivs[[length(ivs) + 1]] <- iv
    sols[[length(sols) + 1]] <- list(lo = lo$x, hi = hi$x)
    tab[[length(tab) + 1]] <- data.frame(trend_branch = gb$label, count_branch = b$label, lower = iv[1], upper = iv[2],
                                         stringsAsFactors = FALSE)
  }
  pcs <- .rxc_union(ivs, tol = 1e-10)
  out <- list(pieces = pcs, lower = if (is.null(pcs)) NA_real_ else min(pcs[, 1]), upper = if (is.null(pcs)) NA_real_ else max(pcs[, 2]),
              branches = do.call(rbind, tab))
  if (!is.null(cs) && !is.null(pcs) && isTRUE(P$integer)) out$witness <- .att_int_witness(P, sys, cs, ivs, sols, out)
  out
}

# Thm E16(b), finite persons: the LP solution at an end gives a whole-person witness when its table n is
# integer, A_r* theta_t is an integer and (B_s*,t(n) - S_t) A_s* is an integer at every post-period.
# Returns c(lower = TRUE/FALSE, upper = TRUE/FALSE); FALSE means "not verified" (the integer end is open).
.att_int_witness <- function(P, sys, cs, ivs, sols, out) {
  Ar <- sum(P$A0[, cs$r]); As <- if (is.null(cs$s)) sum(P$N0) - Ar else sum(P$A0[, cs$s])
  isint <- function(v) all(abs(v - round(v)) <= 1e-7 * pmax(1, abs(v)))
  ok <- function(x) {
    if (is.null(x)) return(FALSE)
    n <- x[seq_len(P$nv)]; xx <- x[P$nv + seq_len(sys$nx)]
    if (!isint(n)) return(FALSE)
    for (j in seq_along(cs$posts)) {
      th <- xx[sys$ix(j, "theta")]
      if (!isint(Ar * th)) return(FALSE)
      if (cs$has_comp && !isint((sum(sys$Bs(cs$posts[j]) * n) - xx[sys$ix(j, "S")]) * As)) return(FALSE)
    }
    TRUE
  }
  at <- function(side) {
    v <- if (side == "lo") out$lower else out$upper
    for (i in seq_along(ivs)) {
      iv <- ivs[[i]]; if (anyNA(iv)) next
      if (abs((if (side == "lo") iv[1] else iv[2]) - v) <= 1e-9 * max(1, abs(v)) && ok(sols[[i]][[side]])) return(TRUE)
    }
    FALSE
  }
  c(lower = at("lo"), upper = at("hi"))
}

# ---------------------------------------------------------------------------
# Def. 5.1: the labelling rule
# ---------------------------------------------------------------------------
.att_label_check <- function(fit, given) {
  miss <- character(0)
  if (!given$population) miss <- c(miss, "the population (population = \"fixed\" or \"stayers\"; Def. 5.1 item 0)")
  if (!given$treated) miss <- c(miss, "the treated group (treated = ; item 1)")
  if (!given$trend) miss <- c(miss, "a trend assumption (trend = no_trend(), trend_dial(M), parallel_trends(), rel_mag(Mbar), smoothness(M) or mtr_nd_ps(); item 2)")
  if (given$trend && given$uses_comp && !given$comparison)
    miss <- c(miss, "the comparison group (comparison = a group or \"rest\"; item 1: this trend assumption uses one)")
  # item 3: only when the trend assumption uses a comparison group (MTR and ND alone need none)
  if (!given$spillover && (!given$trend || given$uses_comp))
    miss <- c(miss, "a spillover assumption (spillover = Lambda, |S_t| <= Lambda; spillover = 0 must be stated explicitly; item 3: this trend assumption uses a comparison group)")
  if (!given$na) miss <- c(miss, "no anticipation (no_anticipation = TRUE; item 4)")
  if (length(miss))
    stop("eiATT refuses a causal label (2e Def. 5.1). Missing:\n  - ", paste(miss, collapse = "\n  - "),
         "\nA statement without these is descriptive: report the eiChange() fit instead.", call. = FALSE)
  invisible(TRUE)
}

#' Causal effects on the treated (ATT) on top of an eiChange() fit
#'
#' \code{eiATT()} is the causal layer of DiD with ecological data (2e Sec. 5). The estimand is
#' \eqn{ATT_t = B_{r^*ct} - B^{(0)}_{r^*ct}}, the effect of the treatment programme on the treated
#' group's rate at post-election \eqn{t}. By Lemma E14, \eqn{ATT_t = D_t - b_t + S_t}: the
#' descriptive DiD minus the counterfactual trend gap plus the spillover onto the comparison
#' group. Every causal set is the projection of a joint linear program over the realised tables
#' of all elections and the potential-outcome quantities (Thm E16): valid, sharp for continuous
#' tables, and refutable (an empty set refutes the conjunction).
#'
#' \strong{Labelling rule (Def. 5.1).} A causal label needs: the population (\code{"fixed"},
#' or \code{"stayers"} with fixed group membership); the treated group; the comparison group
#' when the trend assumption uses one; a trend assumption (there is no default, also not between
#' electorate-wide and within-precinct parallel trends); a spillover assumption
#' \eqn{|S_t| \le \Lambda} (\code{spillover = 0} stated explicitly if assumed) whenever the trend
#' assumption uses a comparison group (MTR and ND alone need none); and no anticipation. If any is missing \code{eiATT()} stops. Descriptive information is printed
#' separately and labelled DESCRIPTIVE.
#'
#' \strong{Forms.} \code{no_trend()}: Prop. E15 (always contains 0). \code{trend_dial(M)}:
#' Cor. E17, with the outer form \eqn{D \oplus [-(M+\Lambda), M+\Lambda]}, the exact criterion
#' of E17(b) for when the sharp ends equal the outer ends, and the breakdown value
#' \eqn{M^* = \min D_t - \Lambda} for "ATT > 0" (E17(c)) with the sharp breakdown value (one LP with
#' M a variable). \code{rel_mag(Mbar)}: Prop. E18, a union of 2L LPs over the joint set, bounding each
#' post-election's increment of the trend gap, so violations accumulate (the marginal outer set at
#' horizon h is \eqn{D_h \oplus [-h\bar M\sup A - \Lambda, h\bar M\sup A + \Lambda]}); branch
#' intervals are merged only where they overlap and holes are reported; the print always shows the implied
#' \eqn{M = \bar M \sup_F A(n)} and the degeneration report of Prop. E19: whether the set
#' equals the trend dial with that M, whether it equals the no-assumption set (no information),
#' the vacuity criterion E19(b) (these two with one post-period only), and the margins-only artefact
#' \eqn{\bar M w_{-1}/2} (E19(c2)).
#' \code{smoothness(M)}: Def. 5.3, one LP. \code{mtr_nd_ps()}: Prop. E21 (with the
#' spillover dial; closed form E21(c) shown under margins only).
#'
#' \strong{Parallel trends in expectation} (\code{within_pt()}, \code{pooled_pt()}; 2e Sec. 5.8-5.9).
#' The estimand is the expected effect \eqn{ATT^{mean}_t}. Three kinds of interval are printed, each
#' labelled: an IDENTIFICATION BOUND under the mean model (Props E28, E30: the LP over expected rates,
#' with LCE rows from \code{lce}, the local dial of E33(a), the rate limits and, for WPT, unsupported
#' precincts bounded with \eqn{\theta_i \in [0,1]}); the DESCRIPTIVE realised-count target \eqn{D_W}
#' (E29, count layer with the fit's rows); and, for WPT, the CONFIDENCE INTERVAL of E29(c),
#' \eqn{[L - \pi_U - r_\alpha, U + r_\alpha]} from one observed election pair (Hoeffding over the
#' independence unit). E29 covers the realised-count route only; it never makes the identification
#' bound a confidence interval. Without \code{law} the observed returns are taken as the expected
#' rates (the known-law reading). Each sign conclusion has a one-LP breakdown value for the local dial
#' \eqn{M} and, with baseline LCE, for the curvature allowance \eqn{\epsilon}, the other held fixed.
#' \eqn{\pi_U}, \eqn{n_{eff} = 1/\sum_i w_{Vi}^2}, the local comparison counts and the independence
#' unit are always printed. Deferred: E31, E32, E33(b)-(c), E34.
#'
#' \strong{Labels.} Every printed interval is labelled IDENTIFICATION SET/BOUND (no sampling
#' statement), DESCRIPTIVE (realised counts) or CONFIDENCE INTERVAL/SET.
#'
#' \strong{Spillovers (Prop. E22).} Pre-periods carry no information on spillover onto the
#' comparison group; \eqn{\Lambda} is its own dial and is always printed.
#'
#' \strong{Model-based pre-trends} (\code{pretrends = "model"}, Prop. E25): with an eiChange fit
#' that has a model layer (\code{fd()} on consecutive pairs, or \code{lce_periods()}), one joint
#' projection region for the stacked regression parameter, box variables for the residual of
#' every actual quantity in the constraints and the target, and the union bound. Under FD alone
#' the rate limits are dropped and the set is labelled OUTER ((H2)).
#'
#' \strong{Inference status.} Count-layer causal sets are exact given the margins, rows and
#' assumptions (no sampling inference); model-based sets inherit the model interval's status.
#' Integer counts: the causal LP is solved for continuous tables (the integer ends are open,
#' Remark 7.1).
#'
#' @param fit an \code{\link{eiChange}} fit.
#' @param treated the treated group (row label or index).
#' @param comparison the comparison group (a row, or \code{"rest"}).
#' @param col the outcome (column); defaults to the fit's default column.
#' @param post post-election(s) to report (default: every election after the base).
#' @param population \code{"fixed"} or \code{"stayers"}; must agree with the fit.
#' @param trend a trend assumption (\code{\link{trend_assumptions}}).
#' @param spillover \eqn{\Lambda \ge 0}.
#' @param no_anticipation must be \code{TRUE}.
#' @param pretrends \code{"counts"} (count layer, default) or \code{"model"} (Prop. E25).
#' @param alpha level for \code{pretrends = "model"}.
#' @param rate_limits \code{TRUE} (default) imposes \eqn{\theta_t \in [0,1]} and
#'   \eqn{B_{s^*ct} - S_t \in [0,1]} (the sharp set). \code{FALSE} drops them: an OUTER set,
#'   printed as such (used to reproduce NS N21's numbers; 2e N21 correction).
#' @param lce LCE rows on expected rates for \code{within_pt()}/\code{pooled_pt()}
#'   (\code{\link{lce_mean}}); default none.
#' @param law optional known aggregate law for the mean model: a list named by election, each a matrix
#'   with one row per unit of equally likely outcome counts (E28(c)); default: the observed returns.
#' @param independence the independence unit of E29(b): \code{"precincts"} (any dependence within a
#'   precinct; default) or \code{"persons"}.
#' @return an object of class \code{"eiATT"}.
#' @examples
#' d0 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(30, 45), no = c(70, 55))
#' d1 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(32, 48), no = c(68, 52))
#' d2 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(30, 60), no = c(70, 40))
#' fit <- eiChange(list(early = d0, base = d1, post = d2), cbind(yes, no) ~ cbind(vet, oth),
#'                 unit = "p", base = "base", row = "vet", col = "yes")
#' eiATT(fit, treated = "vet", comparison = "oth", population = "fixed",
#'       trend = rel_mag(1), spillover = 0, no_anticipation = TRUE)
#' eiATT(fit, treated = "vet", comparison = "oth", population = "fixed",
#'       trend = within_pt(a_min = 1), spillover = 0, no_anticipation = TRUE)
#' @export
eiATT <- function(fit, treated, comparison, col = NULL, post = NULL, population, trend, spillover,
                  no_anticipation, pretrends = c("counts", "model"), alpha = 0.05, rate_limits = TRUE,
                  lce = NULL, law = NULL, independence = c("precincts", "persons")) {
  if (!inherits(fit, "eiChange")) stop("'fit' must be an eiChange() fit")
  pretrends <- match.arg(pretrends)
  if (!missing(trend) && is.character(trend) && identical(trend, "none")) trend <- no_trend()
  if (!missing(trend) && !inherits(trend, "edTrend")) stop("'trend' must be built with no_trend(), trend_dial(), parallel_trends(), rel_mag(), smoothness(), mtr_nd_ps(), within_pt() or pooled_pt()")
  given <- list(population = !missing(population), treated = !missing(treated), trend = !missing(trend),
                comparison = !missing(comparison), spillover = !missing(spillover),
                na = !missing(no_anticipation) && isTRUE(no_anticipation),
                uses_comp = !missing(trend) && .att_uses_comparison(trend))
  .att_label_check(fit, given)
  population <- match.arg(population, c("fixed", "stayers"))
  if (population != fit$population)
    stop("population = \"", population, "\" but the eiChange fit was made for population = \"", fit$population,
         "\": the ATT is defined for the population whose margins were supplied (Def. 5.1 item 0; G1-G4)")
  Lambda <- if (missing(spillover)) 0 else .att_nonneg(spillover, "spillover")
  P <- fit$P
  r <- .ed_rowidx(P, treated, "treated")
  if (length(r) != 1) stop("one treated group")
  s <- if (given$comparison) .ed_versus(P, comparison, r) else NULL
  has_comp <- given$comparison
  cc <- .ed_colidx(P, if (is.null(col)) fit$defaults$col else col, "col")
  P$cc <- cc
  b0 <- fit$base
  if (b0 == P$T) stop("the base election is the last election: there is no post-period")
  posts_all <- seq(b0 + 1, P$T)
  tgt <- if (is.null(post)) posts_all else vapply(post, function(e) .ed_elec(P, e, "post"), 0L)
  if (any(tgt <= b0)) stop("post-periods must come after the base election")
  posts <- posts_all          # code review issue 3: every post-period's constraints, always (Thm E16, vector form)
  pre_el <- if (b0 > 1) seq_len(b0 - 1) else integer(0)
  if (trend$kind == "rm" && !is.null(trend$pre)) {
    pe <- sort(vapply(trend$pre, function(e) .ed_elec(P, e, "rel_mag(pre)"), 0L))
    if (any(pe >= b0) || !identical(as.integer(pe), as.integer(seq(min(pe), b0 - 1)))) stop("rel_mag(pre): contiguous pre-period elections ending just before the base")
    pre_el <- pe
  }
  pre_pairs <- if (length(pre_el)) lapply(pre_el, function(e) c(e, e + 1L)) else list()
  cs <- list(r = r, s = s, base = b0, posts = posts, pre_pairs = pre_pairs, has_comp = has_comp, cc = cc)
  out <- list(call = match.call(), fit = fit, treated = P$rows[r], comparison = if (!has_comp) NA_character_ else if (is.null(s)) "rest" else P$rows[s],
              col = P$cols[cc], base = P$elections[b0], posts = P$elections[tgt], population = population,
              trend = trend, Lambda = Lambda, pretrends = pretrends, alpha = alpha)
  mean_kind <- trend$kind %in% c("wpt", "ppt")
  if (!mean_kind && (!is.null(lce) || !is.null(law)))
    stop("'lce' and 'law' belong to the mean model: use them with within_pt() or pooled_pt()")
  if (mean_kind) {
    if (pretrends == "model") stop("within_pt()/pooled_pt() have their own mean model (lce = lce_mean(...)); pretrends = \"model\" is for Prop. E25")
    if (!isTRUE(rate_limits)) stop("within_pt()/pooled_pt(): the rate limits of Prop. E28(b) are part of the sharp set and are always imposed")
    if (is.null(lce)) lce <- lce_mean(change = FALSE)
    if (!inherits(lce, "edLCEmean")) stop("'lce' must be built with lce_mean()")
    if (!is.null(law) && (!is.list(law) || is.null(names(law)) || !all(names(law) %in% P$elections)))
      stop("'law' must be a list named by election labels")
    independence <- match.arg(independence)
    if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
    out$wpt <- .wpt_run(fit, P, .att_fbranches(fit), cs, trend, Lambda, lce, law, alpha, independence, tgt)
    out$results <- out$wpt$results
    return(structure(out, class = "eiATT"))
  }
  if (pretrends == "model") {
    out$model <- .att_model(fit, cs, trend, Lambda, alpha, tgt)
    out$results <- out$model$results
    return(structure(out, class = "eiATT"))
  }
  fb <- .att_fbranches(fit)
  rate_limits <- isTRUE(rate_limits)
  out$rate_limits <- rate_limits
  sys <- .att_system(P, cs, trend, Lambda, rate_limits)
  none <- .att_system(P, cs, no_trend(), Lambda, TRUE)
  if (trend$kind == "sd" && !length(sys$branches[[1]]$rows))
    out$sd_note <- "smoothness imposes nothing here: no three consecutive elections (Def. 5.3), so the set is the no-assumption set of Prop. E15 (with the spillover limits)"
  res <- list()
  for (t in tgt) {
    k <- match(t, posts)
    rg <- .att_range(P, fb, sys, k, t, cs)
    e <- list(post = P$elections[t], pieces = rg$pieces, lower = rg$lower, upper = rg$upper, branches = rg$branches, witness = rg$witness)
    # descriptive basis
    e$B_r <- .att_nrange(P, fb, sys$Br(t))
    if (has_comp) e$D <- .att_nrange(P, fb, sys$g(t) - sys$g(b0))
    e$none <- .att_range(P, fb, none, k, t)[c("lower", "upper")]
    e$empty <- is.na(rg$lower)
    if (trend$kind == "dial") e$dial <- .att_dial_report(P, fb, sys, cs, k, t, b0, trend$M, Lambda, e)
    if (trend$kind == "rm") e$rm <- .att_rm_report(P, fb, cs, trend, Lambda, k, t, e, rate_limits)
    if (trend$kind == "mtr" && !length(fit$built$common) && !length(fit$built$unions))
      e$e21 <- .att_e21(P, cs, trend, Lambda, t)
    if (trend$kind == "sd" && length(pre_el) == 1 && P$T == b0 + 1 && !length(fit$built$common) && !length(fit$built$unions))
      e$e20 <- .att_e20(P, cs, trend$M, Lambda, t)
    res[[length(res) + 1]] <- e
  }
  out$results <- res
  structure(out, class = "eiATT")
}

# Cor. E17: outer form, exact criterion (b), breakdown (c)
.att_dial_report <- function(P, fb, sys, cs, k, t, b0, M, Lambda, e) {
  D <- e$D; ML <- M + Lambda
  if (is.null(D) || anyNA(D)) return(NULL)        # code review issue 2: refuted count layer -> main set EMPTY (E16(d))
  Dobj <- sys$g(t) - sys$g(b0)
  crit <- function(side) {
    for (b in fb) {
      dr <- .att_nrange(P, list(b), Dobj)
      if (is.na(dr[1])) next
      v <- if (side == "upper") D[2] else D[1]
      if (abs((if (side == "upper") dr[2] else dr[1]) - v) > 1e-9 * max(1, abs(v))) next
      tol <- 1e-9
      rows <- list(if (side == "upper") list(n = Dobj, x = numeric(0), lo = v - tol, hi = Inf) else list(n = Dobj, x = numeric(0), lo = -Inf, hi = v + tol),
                   if (side == "upper") list(n = sys$thPT(t), x = numeric(0), lo = ML - tol, hi = 1 + ML + tol) else list(n = sys$thPT(t), x = numeric(0), lo = -ML - tol, hi = 1 - ML + tol),
                   if (side == "upper") list(n = sys$Bs(t), x = numeric(0), lo = Lambda - tol, hi = Inf) else list(n = sys$Bs(t), x = numeric(0), lo = -Inf, hi = 1 - Lambda + tol))
      if (!is.null(.att_lp(P, b$rows, rows, 0, list(lower = numeric(0), upper = numeric(0)), numeric(P$nv), numeric(0), 1))) return(TRUE)
    }
    FALSE
  }
  list(outer = c(D[1] - ML, D[2] + ML), criterion_upper = crit("upper"), criterion_lower = crit("lower"),
       breakdown = if (D[1] - Lambda > 0) D[1] - Lambda else if (D[2] + Lambda < 0) -(D[2] + Lambda) else NA_real_,
       sharp_breakdown = .att_dial_breakdown(P, fb, cs, Lambda, k, t, e))
}

# E17(c): sharp breakdown value, one LP per count branch with M a variable:
# min M s.t. (n, theta, S, b) jointly feasible, |b_j| <= M for every post-period j, and ATT_t on the other side of 0.
# NA when the set already contains 0 or is empty; Inf when no such witness exists.
.att_dial_breakdown <- function(P, fb, cs, Lambda, k, t, e) {
  if (is.na(e$lower) || (e$lower <= 0 && e$upper >= 0)) return(NA_real_)
  pos <- e$lower > 0
  sys <- .att_system(P, cs, no_trend(), Lambda, TRUE)
  nx <- sys$nx + 1; iM <- nx
  K <- length(cs$posts)
  rows <- sys$common
  rows <- lapply(rows, function(r) { r$x <- c(r$x, 0); r })
  for (j in seq_len(K)) {
    v <- numeric(nx); v[sys$ix(j, "b")] <- 1; v[iM] <- -1
    rows[[length(rows) + 1]] <- list(n = NULL, x = v, lo = -Inf, hi = 0)       # b_j - M <= 0
    v[iM] <- 1
    rows[[length(rows) + 1]] <- list(n = NULL, x = v, lo = 0, hi = Inf)        # b_j + M >= 0
  }
  xa <- numeric(nx); xa[sys$ix(k, "theta")] <- -1                             # ATT = B_r(n) - theta
  rows[[length(rows) + 1]] <- if (pos) list(n = sys$Br(t), x = xa, lo = -Inf, hi = 0) else list(n = sys$Br(t), x = xa, lo = 0, hi = Inf)
  xb <- list(lower = c(sys$lower, 0), upper = c(sys$upper, Inf))
  obj <- numeric(nx); obj[iM] <- 1
  best <- Inf
  for (b in fb) {
    s <- .att_lp(P, b$rows, rows, nx, xb, numeric(P$nv), obj, -1)
    if (!is.null(s)) best <- min(best, s$value)
  }
  best
}

# Prop. E19: degeneration report for relative magnitudes
.att_rm_report <- function(P, fb, cs, tr, Lambda, k, t, e, rate_limits = TRUE) {
  g <- function(tt) .att_rate(P, tt, cs$r) - .att_rate_s(P, tt, cs$r, cs$s)
  pre <- lapply(cs$pre_pairs, function(pp) .att_nrange(P, fb, g(pp[2]) - g(pp[1])))
  if (any(vapply(pre, anyNA, NA)) || is.null(e$D) || anyNA(e$D)) return(NULL)   # code review issue 2
  supA <- max(vapply(pre, function(v) max(abs(v)), 0))
  Mimp <- tr$Mbar * supA
  K <- length(cs$posts); h <- k                    # horizon of post-period k (k-th election after the base)
  one_post <- K == 1                               # E19 is stated for one post-period (code review issue 4)
  dial <- if (one_post) .att_range(P, fb, .att_system(P, cs, trend_dial(Mimp), Lambda, rate_limits), k, t) else list(lower = NA_real_, upper = NA_real_)
  tol <- 1e-7
  eqv <- function(a, b) !is.na(a[1]) && !is.na(b[1]) && abs(a[1] - b[1]) < tol && abs(a[2] - b[2]) < tol
  rmset <- c(e$lower, e$upper)
  # E19(b) vacuity (one post-period, Lambda = 0): upper end = B^U iff exists n with B_r1 = B^U and |thetaPT| <= Mbar A(n)
  vac <- c(NA, NA)
  if (one_post && Lambda == 0 && length(cs$pre_pairs)) {
    sys <- .att_system(P, cs, no_trend(), 0)
    one <- function(side) {
      BU <- e$B_r
      for (b in fb) for (pp in cs$pre_pairs) for (sg in c(1, -1)) {
        d <- sg * tr$Mbar * (g(pp[2]) - g(pp[1]))
        tolr <- 1e-9
        th <- sys$thPT(t)
        # upper: |thetaPT| <= d ; lower: |1 - thetaPT| <= d, i.e. d - thetaPT >= -1 and d + thetaPT >= 1
        r0 <- list(if (side == "upper") list(n = sys$Br(t), x = numeric(0), lo = BU[2] - tolr, hi = Inf) else list(n = sys$Br(t), x = numeric(0), lo = -Inf, hi = BU[1] + tolr),
                   list(n = d - th, x = numeric(0), lo = if (side == "upper") -tolr else -1 - tolr, hi = Inf),
                   list(n = d + th, x = numeric(0), lo = if (side == "upper") -tolr else 1 - tolr, hi = Inf))
        if (!is.null(.att_lp(P, b$rows, r0, 0, list(lower = numeric(0), upper = numeric(0)), numeric(P$nv), numeric(0), 1))) return(TRUE)
      }
      FALSE
    }
    vac <- c(one("lower"), one("upper"))
  }
  b0m <- cs$base - 1
  margins_only <- length(fb) == 1 && !length(fb[[1]]$rows)      # E19(c2) is stated under margins only
  w1 <- if (b0m >= 1 && margins_only) { G <- .att_nrange(P, list(list(rows = list())), g(b0m)); G[2] - G[1] } else NA_real_
  # E18(b): violations accumulate over the horizon, |b_h| <= h Mbar sup A
  ew <- c(e$D[1] - h * Mimp - Lambda, e$D[2] + h * Mimp + Lambda)
  list(pre_ranges = pre, supA = supA, M_implied = Mimp, horizon = h, endpointwise = ew, dial_set = c(dial$lower, dial$upper),
       one_post = one_post, collapses_to_dial = one_post && !is.null(e$pieces) && nrow(e$pieces) == 1 && eqv(rmset, c(dial$lower, dial$upper)),
       no_information = !is.null(e$pieces) && nrow(e$pieces) == 1 && eqv(rmset, c(e$none$lower, e$none$upper)),
       vacuous_lower = vac[1], vacuous_upper = vac[2], artefact = if (is.na(w1)) NA_real_ else tr$Mbar * w1 / 2)
}

# Prop. E21(c): closed form under margins only
.att_e21 <- function(P, cs, tr, Lambda, t) {
  if (!(tr$nd && tr$ps && tr$mtr)) return(NULL)
  e0 <- list(list(rows = list()))
  BU <- .att_nrange(P, e0, .att_rate(P, t, cs$r))[2]
  BL0 <- .att_nrange(P, e0, .att_rate(P, cs$base, cs$r))[1]
  Bs_low <- if (!is.null(cs$s)) .att_nrange(P, e0, .att_rate(P, t, cs$s))[1] else {
    YU <- BU * sum(P$A0[, cs$r]); (sum(P$Bt[[t]][, P$cc]) - YU) / (sum(P$N0) - sum(P$A0[, cs$r])) }
  c(0, BU - max(BL0, Bs_low - Lambda))
}

# Prop. E20: smoothness closed form (margins only, L = 1, one post-period; without rate limits)
.att_e20 <- function(P, cs, M, Lambda, t) {
  e0 <- list(list(rows = list()))
  G <- function(tt) .att_nrange(P, e0, .att_rate(P, tt, cs$r) - .att_rate_s(P, tt, cs$r, cs$s))
  G1 <- G(t); G0 <- G(cs$base); Gm <- G(cs$base - 1)
  c(G1[1] - 2 * G0[2] + Gm[1] - M - Lambda, G1[2] - 2 * G0[1] + Gm[2] + M + Lambda)
}

#' @export
print.eiATT <- function(x, digits = 3, ...) {
  f <- function(v) if (is.null(v) || is.na(v[1])) "empty" else sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), v[1], v[2])
  if (is.null(x$wpt)) cat("CAUSAL (2e Def. 5.1): ATT_t = B_{r*,t} - B^(0)_{r*,t}, the effect on the treated group's realised rate\n")
  else cat("CAUSAL IN EXPECTATION (2e Def. 5.1, Sec. 5.8): the expected effect on the treated group's rate\n")
  cat("  labels: IDENTIFICATION SET/BOUND (no sampling statement), DESCRIPTIVE (realised counts), CONFIDENCE INTERVAL/SET (sampling)\n")
  cat(sprintf("  population: %s\n  treated: %s; comparison: %s; outcome: %s; base election: %s\n",
              if (x$population == "fixed") "fixed population (2e (FP))" else "stayers with fixed group membership (G1-G4: the trend assumption is stated for the stayers)",
              x$treated, if (is.na(x$comparison)) "none (not used by this trend assumption)" else x$comparison, x$col, x$base))
  cat(sprintf("  trend assumption: %s\n  spillover onto the comparison group: |S_t| <= %g%s\n  no anticipation: assumed\n",
              .att_trend_label(x$trend), x$Lambda, if (x$Lambda == 0) " (no spillover, stated explicitly)" else ""))
  fit <- x$fit
  cat(sprintf("  descriptive basis: %s%s\n", if (is.null(fit$restrictions)) "count layer, margins only (Delta-0)" else
    paste("count rows", paste(fit$restrictions$id, collapse = ", ")),
    if (x$pretrends == "model") sprintf("; model rung %s%s", fit$model$rung,
      if (length(fit$model$row_labels)) paste0(" with ", paste(fit$model$row_labels, collapse = "; ")) else "") else ""))
  if (!is.null(x$wpt)) { .wpt_print(x, digits); return(invisible(x)) }
  if (x$pretrends == "model") { .att_print_model(x, digits); return(invisible(x)) }
  if (isFALSE(x$rate_limits))
    cat("  OUTER: the rate limits theta_t in [0,1], B_s,t - S_t in [0,1] are dropped (the sharp set imposes them; Thm E16)\n")
  cat("  inference: count layer, exact given the margins, rows and assumptions (no sampling inference);\n",
      "  continuous tables (integer causal ends open, Remark 7.1)\n", sep = "")
  if (!is.null(x$sd_note)) cat("  ", x$sd_note, "\n", sep = "")
  for (e in x$results) {
    cat(sprintf("\nPost-election %s\n", e$post))
    pcs <- e$pieces
    cat(sprintf("  [IDENTIFICATION SET, realised ATT, count layer (Thm E16)%s] ATT: %s%s\n", if (isFALSE(x$rate_limits)) ", OUTER" else "",
                if (is.null(pcs)) "EMPTY" else paste(apply(pcs, 1, f), collapse = " U "),
                if (is.null(pcs)) " -- refutes (FP) + no anticipation + count rows + the trend and spillover assumptions (Thm E16(d))" else ""))
    cat(sprintf("  [IDENTIFICATION SET, realised ATT, no trend assumption (Prop. E15)] ATT: %s\n", f(c(e$none$lower, e$none$upper))))
    cat(sprintf("  [DESCRIPTIVE, realised counts] B_r* at %s: %s%s\n", e$post, f(e$B_r),
                if (!is.null(e$D)) sprintf("; DiD D_t (base -> %s): %s", e$post, f(e$D)) else ""))
    if (!is.null(e$dial)) {
      d <- e$dial
      cat(sprintf("  [IDENTIFICATION SET, outer] E17(a) outer form D (+) [-(M+L), M+L]: %s\n", f(d$outer)))
      cat(sprintf("  E17(b) exact criterion: upper end %s max D + M + Lambda; lower end %s min D - M - Lambda\n",
                  if (d$criterion_upper) "equals" else "is below", if (d$criterion_lower) "equals" else "is above"))
      cat(sprintf("  E17(c) breakdown value for the sign (outer form, M* = min D - Lambda): %s\n",
                  if (is.na(d$breakdown)) "none (the outer form contains 0 at M = 0)" else sprintf("%.4g", d$breakdown)))
      cat(sprintf("  E17(c) sharp breakdown value (one LP, M a variable): %s\n",
                  if (is.na(d$sharp_breakdown)) "none (the set contains 0 or is empty)" else if (!is.finite(d$sharp_breakdown))
                    "infinite (no table and trend gap put the ATT on the other side of 0)" else sprintf("%.4g", d$sharp_breakdown)))
    }
    if (!is.null(e$rm)) {
      m <- e$rm
      cat(sprintf("  E19 implied M = Mbar * sup|consecutive pre DiD| = %.4g (sup A = %.4g; DESCRIPTIVE pre DiD ranges %s)\n",
                  m$M_implied, m$supA, paste(vapply(m$pre_ranges, f, ""), collapse = ", ")))
      cat(sprintf("  [IDENTIFICATION SET, outer] E18(b) marginal outer set at horizon %d, D (+) [-h Mbar sup A - Lambda, +]: %s (not sharp; joint statements keep the path constraints)\n",
                  m$horizon, f(m$endpointwise)))
      if (m$one_post) cat(sprintf("  [IDENTIFICATION SET] E19(a) trend dial with that M: %s; %s\n", f(m$dial_set),
                  if (m$collapses_to_dial) "the RM set COLLAPSES to this trend dial" else "the RM set is strictly inside it"))
      else cat("  E19(a)-(b) not computed: E19 is stated for one post-period, and later post-periods' constraints are imposed here\n")
      if (m$no_information) cat("  the RM set equals the no-assumption set: NO INFORMATION (Prop. E19(b))\n")
      if (!is.na(m$vacuous_upper))
        cat(sprintf("  E19(b) vacuity: lower end %s; upper end %s\n",
                    if (m$vacuous_lower) "uninformative (= B^L - 1)" else "informative", if (m$vacuous_upper) "uninformative (= B^U)" else "informative"))
      if (!is.na(m$artefact)) cat(sprintf("  E19(c2) under margins only the allowance is at least Mbar w_-1 / 2 = %.4g (an artefact of the pre-period's Duncan-Davis width)\n", m$artefact))
      cat("  E19(c1) every row that narrows the pre-period DiDs narrows the RM set (signs/CR on pre-period changes, contrasts, known cells, a model rung)\n")
    }
    if (!is.null(e$e21)) cat(sprintf("  [IDENTIFICATION SET] E21(c) closed form (margins only): %s\n", f(e$e21)))
    if (!is.null(e$e20)) cat(sprintf("  [IDENTIFICATION SET, outer] E20 closed form (margins only, no rate limits): %s\n", f(e$e20)))
    if (!is.null(e$witness))
      cat(sprintf("  E16(b) whole-person witness at the LP ends: lower %s; upper %s\n",
                  if (e$witness[["lower"]]) "yes" else "not verified (integer end open, Remark 7.1)",
                  if (e$witness[["upper"]]) "yes" else "not verified (integer end open, Remark 7.1)"))
    br <- e$branches
    if (!is.null(pcs) && nrow(pcs) > 1)
      cat(sprintf("  the set has %d pieces: branch intervals merged only where they overlap; the holes are part of the answer (Prop. E18(a))\n", nrow(pcs)))
    if (!is.null(br) && nrow(br) > 1)
      cat(sprintf("  union of %d branch LPs (no multiplicity correction; never a selected branch, Prop. E24)\n", nrow(br)))
  }
  cat("\n  Prop. E22: pre-periods carry no information on spillover onto the comparison group; Lambda is its own dial.\n")
  invisible(x)
}
