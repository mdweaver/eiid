# Changes and DiD, model layer with a baseline-curvature allowance (2e Prop. E36): two groups, two
# elections (base u, post t). Change LCE exact (each group's change mean d_g(x) = c_g + e_g x) and each
# group's base-level mean within epsilon of an affine function at every observed composition:
# m_g0,i = a_g + b_g x_i + r_g,i, |r_g,i| <= epsilon (a fixed constant). beta = (gamma, delta): OLS of
# T_u and of Delta T on (1, x, x^2), one joint sandwich clustered by unit. Rows in P3's form s(beta) + A k >= 0,
# k = (m_V0,i, m_N0,i, a_V, b_V, a_N, b_N, c_V, e_V, c_N, e_N):
#   base-level link (an identity): gamma = P (x m_V0 + (1 - x) m_N0), P = (Z'WZ)^-1 Z'W;
#   curvature: +-(m_g0,i - a_g - b_g x_i) + epsilon >= 0 (epsilon in the intercept);
#   change: delta = (c_N, c_V - c_N + e_N, e_V - e_N);
#   rate limits: m_g0,i in [0, 1], m_g0,i + d_g(x_i) in [0, 1].
# epsilon = 0 is the two-election lce_periods() system; epsilon >= 1/2 leaves change LCE plus the link.

#' Two-election LCE with a baseline-curvature allowance (Prop. E36)
#'
#' A model layer for \code{\link{eiChange}} with two groups: the change of each group's expected
#' rate is affine in the composition (change LCE, exact), and each group's base-level expected
#' rate is within \code{epsilon} of an affine function at every observed composition. The
#' base-level regression enters through its defining identity (the link), so \code{epsilon = 0}
#' is the exact two-election LCE system of \code{lce_periods()} and \code{epsilon >= 1/2} leaves
#' change LCE plus the link (inside \code{fd()}). Intervals are projection intervals for the
#' model-implied DiD plus Lemma E12's allowance. Used by \code{eiATT(trend = within_pt(), lce = TRUE)}
#' (Prop. E35 with E13), which also reports the breakdown value \eqn{\epsilon^*} of the sign.
#' @param epsilon nonnegative baseline-curvature allowance, a fixed constant.
#' @param post the post election (default: the election after the fit's base).
#' @inheritParams change_models
#' @return a model-layer specification.
#' @export
lce_levels_change <- function(epsilon = 0, post = NULL, vcov = "HC3", weights = c("none", "n"),
                              ci = c("projection", "conservative"), allowance = c("popoviciu", "hoeffding"), ndraw = 200, seed = 1)
  structure(list(rung = "LCEeps", epsilon = .att_nonneg(epsilon, "lce_levels_change(epsilon)"), post = post, rows = list(),
                 vcov = vcov, weights = match.arg(weights), ci = match.arg(ci), allowance = match.arg(allowance),
                 ndraw = ndraw, seed = seed), class = "edModel")

# the regression and the row system for groups (r = V, other = N), column cc, elections u -> t
.lce_eps_build <- function(P, spec, r, cc, u, t) {
  if (P$R != 2) stop("lce_levels_change(): Prop. E36 is stated for two groups (R = 2)")
  q <- setdiff(1:2, r); p <- P$p0
  N <- as.double(P$N0); x <- as.double(P$A0[, r]) / N
  Z <- cbind(1, x, x^2)
  if (qr(Z)$rank < 3) stop("lce_levels_change(): needs at least three distinct compositions")
  sw <- if (spec$weights == "n") sqrt(N) else rep(1, p)
  Zw <- Z * sw
  T0 <- as.double(P$Bt[[u]][, cc]) / N; dT <- as.double(P$Bt[[t]][, cc]) / N - T0
  ZtZi <- solve(crossprod(Zw)); Pm <- ZtZi %*% t(Zw * sw)          # gamma = Pm %*% tau0
  Y <- cbind(T0, dT) * sw
  B <- ZtZi %*% crossprod(Zw, Y)
  vc <- .rxc_vcov(Zw, Y - Zw %*% B, type = spec$vcov)
  nz <- 2 * p + 8
  iV <- seq_len(p); iN <- p + seq_len(p)
  ia <- function(g) 2 * p + 2 * (g - 1) + 1; ib <- function(g) ia(g) + 1     # g = 1: V, 2: N
  ic <- function(g) 2 * p + 4 + 2 * (g - 1) + 1; ie <- function(g) ic(g) + 1
  mb <- list(beta = as.vector(B), X = cbind(x, 1 - x), V = vc$V, vcov_type = vc$type, nbeta = 6L, nz = nz, rung = "LCEeps",
             grouped = FALSE, spec = spec)
  zf <- function(pairs, const = 0, coef = numeric(6)) { a <- numeric(nz); for (pp in pairs) a[pp[1]] <- a[pp[1]] + pp[2]
    list(coef = coef, const = const, a = a) }
  eb <- function(j) { v <- numeric(6); v[j] <- 1; v }
  cons_at <- function(eps) {
    cons <- .ed_cons_new(mb)
    add <- function(f, kind) cons <<- .ed_cons_add(cons, f, kind)
    eq <- function(f, kind) { add(f, kind); add(.ed_fscale(f, -1), kind) }
    for (j in 1:3) {                                                   # link: gamma_j - sum_i P_ji tau0_i = 0
      pr <- c(lapply(iV, function(i) c(i, -Pm[j, i] * x[i])), lapply(seq_len(p), function(i) c(iN[i], -Pm[j, i] * (1 - x[i]))))
      eq(zf(pr, coef = eb(j)), "link")
    }
    eq(zf(list(c(ic(2), -1)), coef = eb(4)), "change")                                     # delta1 = c_N
    eq(zf(list(c(ic(1), -1), c(ic(2), 1), c(ie(2), -1)), coef = eb(5)), "change")          # delta2 = c_V - c_N + e_N
    eq(zf(list(c(ie(1), -1), c(ie(2), 1)), coef = eb(6)), "change")                        # delta3 = e_V - e_N
    for (g in 1:2) for (i in seq_len(p)) {
      m <- if (g == 1) iV[i] else iN[i]
      dev <- list(c(m, 1), c(ia(g), -1), c(ib(g), -x[i]))
      add(.ed_fsum(zf(dev), zf(list(), eps)), "curvature")
      add(.ed_fsum(.ed_fscale(zf(dev), -1), zf(list(), eps)), "curvature")
      add(zf(list(c(m, 1))), "rate"); add(zf(list(c(m, -1)), 1), "rate")
      lv <- list(c(m, 1), c(ic(g), 1), c(ie(g), x[i]))
      add(zf(lv), "rate"); add(.ed_fscale(zf(lv), -1), "rate"); cons$s0[length(cons$s0)] <- cons$s0[length(cons$s0)] + 1
    }
    cons
  }
  # model-implied sum_i (WV_i d_V(x_i) + WN_i d_N(x_i))
  change_form <- function(WV, WN) zf(list(c(ic(1), sum(WV)), c(ie(1), sum(WV * x)), c(ic(2), sum(WN)), c(ie(2), sum(WN * x))))
  list(mb = mb, cons_at = cons_at, change_form = change_form, r = r, q = q, u = u, t = t, cc = cc)
}

.lce_eps_post <- function(P, spec, b0) {
  t <- if (is.null(spec$post)) b0 + 1L else .ed_elec(P, spec$post, "lce_levels_change(post)")
  if (t > P$T || t == b0) stop("lce_levels_change(): the post election must differ from the base")
  t
}

# model interval for a form at epsilon: projection (or conservative) union over systems, plus allowance.
# spec$.plugin / spec$.t_E are internal (tests reproduce the theory's population example: plug-in, no allowance)
.lce_eps_interval <- function(bd, eps, f, alpha, reps_terms, P) {
  spec <- bd$mb$spec
  tgm <- list(cvec = f$coef, c0 = f$const, gamma = f$a)
  a_est <- 0.9 * alpha; a_E <- 0.1 * alpha
  al <- if (!is.null(spec$.t_E)) list(value = spec$.t_E, S = NA_real_, lindeberg = NA_real_, rho = .ed_rho(P, reps_terms)) else
    .ed_allowance(P, list(list(terms = reps_terms)), spec$allowance, a_E)
  cons <- bd$cons_at(eps)
  pl <- .ed_plugin(bd$mb, cons, tgm)
  iv <- if (isTRUE(spec$.plugin)) pl else {
    Zs <- .ed_with_seed(spec$seed, .rxc_zdraws(bd$mb$V, 4000))
    .ed_with_seed(spec$seed, .ed_proj(bd$mb, cons, tgm, a_est, spec$ci, Zs, spec$ndraw))
  }
  list(model = if (anyNA(iv)) NULL else matrix(iv + c(-1, 1) * al$value, 1, dimnames = list(NULL, c("lower", "upper"))),
       plugin = pl, allowance = al$value, S = al$S, lindeberg = al$lindeberg, rho = al$rho)
}

# eiChange model block for lce_levels_change(): DiD targets over (base, post) of the default groups
.lce_eps_block <- function(fit, spec, alpha) {
  P <- fit$P; r <- .ed_rowidx(P, fit$defaults$row, "row"); cc <- .ed_colidx(P, fit$defaults$col, "col")
  b0 <- fit$base; t <- .lce_eps_post(P, spec, b0)
  bd <- .lce_eps_build(P, spec, r, cc, b0, t)
  belief <- !is.null(fit$restrictions) && any(fit$restrictions$role == "belief")
  out <- list()
  for (j in seq_along(fit$targets)) {
    tg <- fit$targets[[j]]; tb <- fit$table[j, ]
    row <- data.frame(target = tg$label, count_lower = tb$lower, count_upper = tb$upper, plugin_lower = NA_real_, plugin_upper = NA_real_,
                      model_lower = NA_real_, model_upper = NA_real_, lower = NA_real_, upper = NA_real_, disjoint = NA,
                      allowance = NA_real_, S_phi = NA_real_, lindeberg = NA_real_, status = NA_character_, note = "", stringsAsFactors = FALSE)
    ok <- identical(tg$kind, "did") && identical(as.integer(tg$c), as.integer(cc)) && setequal(c(tg$u, tg$t), c(b0, t)) && isTRUE(tg$design)
    if (!ok) { row$note <- sprintf("lce_levels_change() models the change %s -> %s of column %s only (Prop. E36)",
                                        P$elections[b0], P$elections[t], P$cols[cc]); out[[j]] <- row; next }
    W <- .ed_gapW(P, tg$S, tg$r, tg$s); sg <- if (tg$t == t) 1 else -1
    f <- .ed_fscale(bd$change_form(W[, r], W[, bd$q]), sg)
    iv <- .lce_eps_interval(bd, spec$epsilon, f, alpha, tg$reps[[1]]$terms, P)
    row$plugin_lower <- iv$plugin[1]; row$plugin_upper <- iv$plugin[2]
    row$allowance <- iv$allowance; row$S_phi <- iv$S; row$lindeberg <- iv$lindeberg
    row$status <- "asymptotic guarantee for each fixed epsilon (Prop. E36(a): P3.5 with Lemma E12's allowance)"
    if (is.null(iv$model)) { row$note <- "the model interval is empty: the model layer rejects the model by itself"; row$disjoint <- TRUE }
    else {
      row$model_lower <- iv$model[1]; row$model_upper <- iv$model[2]
      I <- .ec_intersect(fit$solutions[[j]]$pieces, iv$model)
      row$disjoint <- is.null(I) && !is.null(fit$solutions[[j]]$pieces)
      if (!is.null(I)) { row$lower <- min(I[, 1]); row$upper <- max(I[, 2]) }
    }
    out[[j]] <- row
  }
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(table = tab, extra = list(), rung = "LCEeps", spec = spec, row_labels = character(0), alpha = alpha, belief = belief,
       mb = bd$mb, base = bd$cons_at(spec$epsilon), unions = list(), grouped = FALSE, epsilon = spec$epsilon, post = t)
}
