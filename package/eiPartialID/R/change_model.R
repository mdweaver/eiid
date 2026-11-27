# Changes and DiD, model layer: first-difference LCE (FD) and per-period pooled LCE with model rows
# on changes; actual-target intervals with the two-election allowance; the cross-layer rule; and the
# joint projection region for model-based pre-trends in eiATT().
# Theory: docs/extension/PHASE2e_did_general.md (2e) Lemma E7, Props E8, E9, E11, Lemma E12,
# Thm E13, Props E23, E24, E25; inference machinery of PHASE3 (P3.5(a), (b2)(ii)) in rxc_infer.R.
#
# A model layer is a stack of blocks sharing one regression design (same X under (FP)):
#   FD: block k is the regression of Delta T over the consecutive pair (t_{k-1}, t_k) on the
#       homogeneous quadratic (beta_k = Delta Q, z_k = Delta k); st$change = TRUE (sum_c Delta H_c = 0).
#   LCE: block t is pooled M0 (or a grouped rung) for election t's levels.
# The stacked beta has one joint sandwich covariance (.rxc_vcov with all responses), as P3.5(b2)(ii)
# (ii-a) needs for a stacked parameter (2e E25 (H1)); consecutive non-overlapping pairs and
# per-election levels are not singular by construction.

# ---------------------------------------------------------------------------
# constructors
# ---------------------------------------------------------------------------

#' Model layers on changes for eiChange()
#'
#' \code{fd()}: the first-difference linear contextual-effects model (2e Sec. 3): in each pair of
#' consecutive elections \eqn{E[\Delta\beta \mid X] = \Delta H X} with levels unrestricted.
#' \eqn{\Delta Q} is identified by the regression of \eqn{\Delta T} on the homogeneous quadratic
#' (Lemma E7(a)); \eqn{\Delta K} is not, so the model DiD is \eqn{I + \ell(\Delta K)}
#' (Lemma E7(b)). Rows (\code{rows}, a list of specifications, each optionally with
#' \code{pair = c(from, to)}, default every pair):
#' \itemize{
#'   \item \code{list(type = "R1", col)}, \code{"R2"}, \code{"R4"} (R = 2; 2d Prop. FD2): no
#'     contextual change for the comparison group, for the treated group, or on average;
#'     \code{list(type = "delta", col, delta)}: \eqn{|\Delta w_1| \le \delta} (2d FD3). The first
#'     group is the treated ("v") group, the second the comparison.
#'   \item \code{list(type = "didS", rows = c(r, s), col, at = "hull" or a matrix of
#'     compositions)}: DiD-S, the ordering of changes and the composition-weighted spillover sign
#'     (Prop. E9); for R = 2 the closed form \eqn{[I - s q_{min}, I + s q_{min}]} is printed.
#'   \item \code{list(type = "sign", sign_type = "dW" or "dB", row, col, sign, version)}: a model
#'     sign on changes (Prop. E23(b)).
#' }
#' Validity rows \eqn{|\Delta m_{rc}(x)| \le 1} at the vertices of the composition hull always
#' hold (\eqn{\Delta\beta \in [-1, 1]}). \code{sharp = TRUE} (2 x 2, one pair) adds the sharp
#' identified set of Prop. E8(d) evaluated at the empirical law, with FD4's outer set beside it:
#' a plug-in with no sampling inference (inference for the sharp set is open, E23(c)).
#'
#' \code{lce_periods()}: per-period LCE (pooled M0 by default) in every election, with model rows
#' on changes (Prop. E11): \code{list(type = "sign", sign_type, row, col, from, to, sign,
#' version)} or \code{list(type = "cr", row, col, from, to, version)} (a union of two branches,
#' Prop. E24). A grouped rung (\code{model = "M2"} with \code{group}) is computed but flagged:
#' no coverage theorem for model rows on changes under grouped rungs (R16 open).
#'
#' Intervals are projection intervals (P3.5(b2)(ii); \code{ci = "conservative"}: P3.5(a)) for the
#' model-implied target plus the residual allowance of Lemma E12 (\eqn{\rho_i = 2(w_{ir} + w_{is})}
#' for a DiD; Popoviciu, or Hoeffding with \code{allowance = "hoeffding"}), so they cover the
#' actual target; the count set is then intersected with them (Thm E13).
#' @param pairs FD: list of \code{c(from, to)} consecutive pairs (default: every consecutive pair).
#' @param rows list of row specifications (see Details).
#' @param sharp FD: also compute the plug-in sharp set of Prop. E8(d) (2 x 2).
#' @param model,group \code{lce_periods}: the levels rung per election and its grouping variable
#'   (a vector with one value per unit, in the order of the first election's data).
#' @param vcov,weights regression covariance (\code{"HC3"} etc.) and weights (\code{"none"}, \code{"n"}).
#' @param ci \code{"projection"} or \code{"conservative"}.
#' @param allowance \code{"popoviciu"} or \code{"hoeffding"}.
#' @param ndraw,seed draws for candidate collection and the random seed.
#' @return a model-layer specification.
#' @name change_models
NULL

#' @rdname change_models
#' @export
fd <- function(pairs = NULL, rows = list(), sharp = FALSE, vcov = "HC3", weights = c("none", "n"),
               ci = c("projection", "conservative"), allowance = c("popoviciu", "hoeffding"), ndraw = 200, seed = 1)
  structure(list(rung = "FD", pairs = pairs, rows = rows, sharp = sharp, vcov = vcov, weights = match.arg(weights),
                 ci = match.arg(ci), allowance = match.arg(allowance), ndraw = ndraw, seed = seed), class = "edModel")

#' @rdname change_models
#' @export
lce_periods <- function(rows = list(), model = c("M0", "M2", "M1"), group = NULL, vcov = "HC3", weights = c("none", "n"),
                        ci = c("projection", "conservative"), allowance = c("popoviciu", "hoeffding"), ndraw = 200, seed = 1)
  structure(list(rung = "LCE", rows = rows, model = match.arg(model), group = group, vcov = vcov, weights = match.arg(weights),
                 ci = match.arg(ci), allowance = match.arg(allowance), ndraw = ndraw, seed = seed), class = "edModel")

# ---------------------------------------------------------------------------
# the stacked model layer
# ---------------------------------------------------------------------------
.ed_mb <- function(P, spec) {
  X <- P$A0 / P$N0; N <- P$N0
  Tt <- lapply(P$Bt, function(B) B / N)
  if (spec$rung == "FD") {
    pairs <- if (is.null(spec$pairs)) lapply(seq_len(P$T - 1), function(k) c(k, k + 1L)) else
      lapply(spec$pairs, function(pr) c(.ed_elec(P, pr[1], "fd(pairs)"), .ed_elec(P, pr[2], "fd(pairs)")))
    for (pr in pairs) if (pr[2] != pr[1] + 1)
      stop("fd(pairs): use consecutive elections. Overlapping or non-consecutive pairs make the stacked covariance ",
           "singular by an exact identity (2e Prop. E25 (H1)); a DiD over a longer span is the sum of consecutive blocks")
    if (anyDuplicated(vapply(pairs, function(p) p[1], 0L))) stop("fd(pairs): each pair once")
    blocks <- lapply(pairs, function(pr) {
      st <- .rxc_setup(X, Tt[[pr[2]]] - Tt[[pr[1]]], N, model = "M0", vcov = spec$vcov, weights = spec$weights)
      st$change <- TRUE; st$pair <- pr; st
    })
  } else {
    pairs <- NULL
    blocks <- lapply(seq_len(P$T), function(t) {
      st <- .rxc_setup(X, Tt[[t]], N, group = spec$group, model = spec$model, vcov = spec$vcov, weights = spec$weights)
      st$election <- t; st
    })
  }
  nb <- vapply(blocks, function(b) as.integer(b$nbeta), 0L); nzb <- vapply(blocks, function(b) as.integer(b$nz), 0L)
  boff <- c(0L, cumsum(nb))[seq_along(blocks)]; zoff <- c(0L, cumsum(nzb))[seq_along(blocks)]
  # one joint sandwich over all responses (same design Z and weights in every block)
  sw <- if (spec$weights == "n") sqrt(N) else rep(1, length(N))
  Zw <- blocks[[1]]$Z * sw
  for (b in blocks) if (!isTRUE(all.equal(b$Z, blocks[[1]]$Z))) stop("internal: blocks with different designs")
  vc <- .rxc_vcov(Zw, do.call(cbind, lapply(blocks, function(b) b$resid)), type = spec$vcov)
  list(rung = spec$rung, spec = spec, blocks = blocks, pairs = pairs, boff = boff, zoff = zoff,
       nbeta = sum(nb), nz = sum(nzb), beta = unlist(lapply(blocks, function(b) b$beta)), V = vc$V, vcov_type = vc$type,
       X = X, N = N, grouped = spec$rung == "LCE" && !identical(spec$model, "M0"))
}

# a block form (coef on beta_k, const, a on z_k) in joint coordinates
.ed_pad <- function(mb, k, f) {
  cf <- numeric(mb$nbeta); a <- numeric(mb$nz)
  cf[mb$boff[k] + seq_along(f$coef)] <- f$coef
  if (length(f$a)) a[mb$zoff[k] + seq_along(f$a)] <- f$a
  list(coef = cf, const = f$const, a = a)
}
.ed_fsum <- function(...) {
  fs <- Filter(Negate(is.null), list(...))
  Reduce(function(x, y) list(coef = x$coef + y$coef, const = x$const + y$const, a = x$a + y$a), fs)
}
.ed_fscale <- function(f, s) list(coef = s * f$coef, const = s * f$const, a = s * f$a)
.ed_fzero <- function(mb) list(coef = numeric(mb$nbeta), const = 0, a = numeric(mb$nz))

# weights W[i, s] of the gap g_{rs,c}(S) (s = NULL: rest), units in the model's order
.ed_gapW <- function(P, S, r, s) {
  W <- matrix(0, P$p0, P$R)
  W[S, r] <- P$A0[S, r] / sum(P$A0[S, r])
  if (is.null(s)) { rest <- setdiff(seq_len(P$R), r); den <- sum(P$N0[S]) - sum(P$A0[S, r])
    for (q in rest) W[S, q] <- -P$A0[S, q] / den } else W[S, s] <- -P$A0[S, s] / sum(P$A0[S, s])
  W
}
.ed_rateW <- function(P, S, r) { W <- matrix(0, P$p0, P$R); W[S, r] <- P$A0[S, r] / sum(P$A0[S, r]); W }

# model-implied sum_i sum_s W[i, s] m_{s,c}(X_i) at election t (LCE) as a joint form
.ed_lce_lin <- function(mb, t, c, W) .ed_pad(mb, t, .rxc_mlin(mb$blocks[[t]], c, W))
# model-implied change of sum W m over (u, t) under FD: the sum of the consecutive blocks
.ed_fd_change <- function(mb, u, t, c, W) {
  sg <- if (t > u) 1 else -1; lo <- min(u, t); hi <- max(u, t)
  out <- .ed_fzero(mb)
  for (e in seq(lo + 1, hi)) {
    k <- which(vapply(mb$pairs, function(p) p[2] == e, TRUE))
    if (length(k) != 1) stop("the FD layer has no block for the pair (", e - 1, ", ", e, "): fit fd() on every consecutive pair the target spans")
    out <- .ed_fsum(out, .ed_pad(mb, k, .rxc_mlin(mb$blocks[[k]], c, W)))
  }
  .ed_fscale(out, sg)
}

# model-implied change of a quantity with weights W between u and t (either rung)
.ed_model_change <- function(mb, u, t, c, W) {
  if (mb$rung == "FD") return(.ed_fd_change(mb, u, t, c, W))
  .ed_fsum(.ed_lce_lin(mb, t, c, W), .ed_fscale(.ed_lce_lin(mb, u, c, W), -1))
}

# model-implied form of an eiChange target; NULL with a reason when undefined
.ed_model_target <- function(P, mb, tg) {
  if (identical(tg$kind, "did")) return(.ed_model_change(mb, tg$u, tg$t, tg$c, .ed_gapW(P, tg$S, tg$r, tg$s)))
  if (identical(tg$kind, "bin_did"))       # code review issue 10
    return(structure(list(), reason = paste("no model-implied version of a binned DiD is implemented: bin aggregates are",
                                            "evaluated at row-specific weighted compositions, which E10's remark excludes")))
  if (mb$rung == "FD") return(structure(list(), reason = "FD leaves levels free: a level target has no model-implied value (2e Lemma E12(i))"))
  if (identical(tg$kind, "gap_at")) return(.ed_lce_lin(mb, tg$t, tg$c, .ed_gapW(P, tg$S, tg$r, tg$s)))
  if (identical(tg$kind, "rate_at")) return(.ed_lce_lin(mb, tg$t, tg$c, .ed_rateW(P, tg$S, tg$r)))
  structure(list(), reason = "no model-implied version of this target is implemented")
}

# ---------------------------------------------------------------------------
# constraint system in P3's form s(beta) + A z >= 0
# ---------------------------------------------------------------------------
.ed_cons_new <- function(mb) list(S = matrix(0, 0, mb$nbeta), s0 = numeric(0), A = matrix(0, 0, mb$nz), kind = character(0), module = character(0))
.ed_cons_add <- function(cons, f, kind, module = kind) {        # row f >= 0
  cons$S <- rbind(cons$S, f$coef); cons$s0 <- c(cons$s0, f$const); cons$A <- rbind(cons$A, f$a)
  cons$kind <- c(cons$kind, kind); cons$module <- c(cons$module, module); cons
}
.ed_cons_bind <- function(a, b) list(S = rbind(a$S, b$S), s0 = c(a$s0, b$s0), A = rbind(a$A, b$A),
                                     kind = c(a$kind, b$kind), module = c(a$module, b$module))

.ed_mb_validity <- function(mb) {
  cons <- .ed_cons_new(mb)
  for (k in seq_along(mb$blocks)) {
    st <- mb$blocks[[k]]
    if (mb$rung == "LCE") {
      c0 <- .rxc_constraints(st)
      if (length(c0$s0)) for (j in seq_along(c0$s0))
        cons <- .ed_cons_add(cons, .ed_pad(mb, k, list(coef = c0$S[j, ], const = c0$s0[j], a = c0$A[j, ])), "validity")
      next
    }
    Vg <- st$vertices[[1]]
    for (v in seq_len(nrow(Vg))) for (c in seq_len(st$C)) for (r in seq_len(st$R)) {
      f <- .ed_pad(mb, k, .rxc_rate(st, 1, c, r, Vg[v, ]))
      cons <- .ed_cons_add(cons, list(coef = f$coef, const = f$const + 1, a = f$a), "validity")    # Delta m >= -1
      cons <- .ed_cons_add(cons, list(coef = -f$coef, const = 1 - f$const, a = -f$a), "validity") # Delta m <= 1
    }
  }
  cons
}

.ed_blocks_for <- function(P, mb, s) {
  if (is.null(s$pair)) return(seq_along(mb$blocks))
  u <- .ed_elec(P, s$pair[1], "pair"); t <- .ed_elec(P, s$pair[2], "pair")
  k <- which(vapply(mb$pairs, function(p) p[1] == u && p[2] == t, TRUE))
  if (!length(k)) stop("no FD block for pair ", P$elections[u], " -> ", P$elections[t])
  k
}

# FD rows (2d FD2/FD3 for R = 2; DiD-S, E9; model signs on changes, E23(b))
.ed_fd_rows <- function(P, mb, rows, meta) {
  cons <- .ed_cons_new(mb); labs <- character(0); unions <- list()
  e <- function(j) { v <- numeric(P$R); v[j] <- 1; v }
  for (s in rows) {
    tp <- s$type
    if (is.null(tp)) stop("fd rows: each row needs a 'type'")
    c <- .ed_colidx(P, if (is.null(s$col)) meta$col else s$col, "col")
    for (k in .ed_blocks_for(P, mb, s)) {
      st <- mb$blocks[[k]]; pl <- sprintf("%s -> %s", P$elections[mb$pairs[[k]][1]], P$elections[mb$pairs[[k]][2]])
      rate <- function(r, v) .ed_pad(mb, k, .rxc_rate(st, 1, c, r, v))
      if (tp %in% c("R1", "R2", "R4", "delta")) {
        if (P$R != 2) stop("fd rows R1, R2, R4, delta are stated for R = 2 (2d Prop. FD2, FD3); for R >= 3 use didS or sign rows")
        v <- meta$r; n <- setdiff(1:2, v)
        dw1 <- .ed_fsum(rate(n, e(v)), .ed_fscale(rate(n, e(n)), -1))      # slope of Delta m_n in x_v
        db1 <- .ed_fsum(rate(v, e(v)), .ed_fscale(rate(v, e(n)), -1))      # slope of Delta m_v in x_v
        if (tp == "R1") f <- dw1
        if (tp == "R2") f <- db1
        if (tp == "R4") {
          p <- P$N0 / sum(P$N0); x <- P$A0[, v] / P$N0
          wv <- function(w) { w <- w / sum(w); sum(w * x^2) - sum(w * x)^2 }
          a <- sum(p * x) * wv(p * x); b <- sum(p * (1 - x)) * wv(p * (1 - x))
          f <- .ed_fsum(.ed_fscale(db1, a), .ed_fscale(dw1, b))
        }
        if (tp == "delta") {
          d <- .att_nonneg(s$delta, "delta row: delta")
          cons <- .ed_cons_add(cons, list(coef = dw1$coef, const = dw1$const + d, a = dw1$a), "delta")
          cons <- .ed_cons_add(cons, list(coef = -dw1$coef, const = d - dw1$const, a = -dw1$a), "delta")
          labs <- c(labs, sprintf("|Delta w_1| <= %g on %s, %s (2d FD3)", d, P$cols[c], pl))
        } else {
          cons <- .ed_cons_add(cons, f, tp); cons <- .ed_cons_add(cons, .ed_fscale(f, -1), tp)
          labs <- c(labs, sprintf("%s on %s, %s (2d Prop. FD2)", tp, P$cols[c], pl))
        }
      } else if (tp == "didS") {
        rs <- if (is.null(s$rows)) c(meta$r, if (is.null(meta$s)) setdiff(1:2, meta$r)[1] else meta$s) else
          vapply(s$rows, function(z) .ed_rowidx(P, z, "didS rows"), 0L)
        if (length(rs) != 2 || anyNA(rs)) stop("didS: 'rows' must name two groups")
        A <- if (is.null(s$at) || identical(s$at, "hull")) st$vertices[[1]] else as.matrix(s$at)
        if (ncol(A) != P$R) stop("didS: compositions must have one column per group")
        for (j in seq_len(nrow(A))) {
          x <- A[j, ]
          O <- .ed_fsum(rate(rs[1], x), .ed_fscale(rate(rs[2], x), -1))                    # (Delta H x)_r - (Delta H x)_s
          L <- .ed_fsum(Reduce(.ed_fsum, lapply(seq_len(P$R), function(q) .ed_fscale(rate(q, e(rs[1])), x[q]))),
                        .ed_fscale(Reduce(.ed_fsum, lapply(seq_len(P$R), function(q) .ed_fscale(rate(q, e(rs[2])), x[q]))), -1))
          cons <- .ed_cons_add(cons, O, "didS"); cons <- .ed_cons_add(cons, L, "didS")
        }
        labs <- c(labs, sprintf("DiD-S (ordering of changes and composition-weighted spillover sign) %s vs %s at %d compositions on %s, %s (Prop. E9)",
                                P$rows[rs[1]], P$rows[rs[2]], nrow(A), P$cols[c], pl))
        meta$didS[[length(meta$didS) + 1]] <- list(k = k, c = c, rs = rs, A = A)
      } else if (tp == "sign") {
        sty <- if (is.null(s$sign_type)) "dB" else s$sign_type
        ver <- if (is.null(s$version)) "W1" else s$version
        r <- .ed_rowidx(P, if (is.null(s$row)) meta$r else s$row, "row")
        oth <- if (!is.null(s$other)) .ed_rowidx(P, s$other, "other") else NULL
        sg <- .ed_sign_arg(s[["sign"]], "fd sign row")
        f <- .ed_pad(mb, k, .rxc_assoc_form(st, sty, ver, r, c, oth))
        cons <- .ed_cons_add(cons, .ed_fscale(f, sg), "sign")
        labs <- c(labs, sprintf("model %s on changes %s %s on %s, %s (%s; E23(b))", sty, P$rows[r], if (sg > 0) ">= 0" else "<= 0", P$cols[c], pl, ver))
      } else stop("unknown fd row type '", tp, "'")
    }
  }
  list(cons = cons, labels = labs, unions = unions, meta = meta)
}

# LCE rows on changes (Prop. E11): model signs on changes, CR (a union)
.ed_lce_rows <- function(P, mb, rows, meta) {
  cons <- .ed_cons_new(mb); labs <- character(0); unions <- list()
  for (s in rows) {
    tp <- s$type
    c <- .ed_colidx(P, if (is.null(s$col)) meta$col else s$col, "col")
    r <- .ed_rowidx(P, if (is.null(s$row)) meta$r else s$row, "row")
    u <- .ed_elec(P, s$from, "lce row (from)"); t <- .ed_elec(P, s$to, "lce row (to)")
    ver <- if (is.null(s$version)) "W1" else s$version
    oth <- if (!is.null(s$other)) .ed_rowidx(P, s$other, "other") else NULL
    af <- function(sty, v) .ed_fsum(.ed_pad(mb, t, .rxc_assoc_form(mb$blocks[[t]], sty, v, r, c, oth)),
                                   .ed_fscale(.ed_pad(mb, u, .rxc_assoc_form(mb$blocks[[u]], sty, v, r, c, oth)), -1))
    pl <- sprintf("%s -> %s", P$elections[u], P$elections[t])
    if (tp == "sign") {
      sty <- if (is.null(s$sign_type)) "dB" else s$sign_type
      sg <- .ed_sign_arg(s[["sign"]], "lce sign row")
      cons <- .ed_cons_add(cons, .ed_fscale(af(sty, ver), sg), "sign")
      labs <- c(labs, sprintf("model %s* on changes %s %s on %s, %s (%s; Prop. E11)", sty, P$rows[r], if (sg > 0) ">= 0" else "<= 0", P$cols[c], pl, ver))
    } else if (tp == "cr") {
      vB <- if (ver == "pair") "pair" else "W1"; vW <- if (ver == "W1") "W1" else "W2"
      fB <- af("dB", vB); fW <- af("dW", vW)
      plus <- .ed_cons_add(.ed_cons_add(.ed_cons_new(mb), fB, "cr"), fW, "cr")
      minus <- .ed_cons_add(.ed_cons_add(.ed_cons_new(mb), .ed_fscale(fB, -1), "cr"), .ed_fscale(fW, -1), "cr")
      unions[[length(unions) + 1]] <- list(plus = plus, minus = minus, id = sprintf("model CR %s", pl))
      labs <- c(labs, sprintf("model CR on changes %s on %s, %s (%s; Prop. E11: a union of two branches)", P$rows[r], P$cols[c], pl, ver))
    } else stop("unknown lce_periods row type '", tp, "'")
  }
  list(cons = cons, labels = labs, unions = unions, meta = meta)
}

# all systems: validity + rows + one side of each union
.ed_mb_systems <- function(base, unions) {
  U <- length(unions)
  if (!U) return(list(list(cons = base, label = "")))
  combos <- as.matrix(expand.grid(rep(list(1:2), U)))
  lapply(seq_len(nrow(combos)), function(k) {
    cons <- base; lab <- character(0)
    for (u in seq_len(U)) { sd <- if (combos[k, u] == 1) "plus" else "minus"
      cons <- .ed_cons_bind(cons, unions[[u]][[sd]]); lab <- c(lab, paste0(unions[[u]]$id, if (sd == "plus") "+" else "-")) }
    list(cons = cons, label = paste(lab, collapse = ", "))
  })
}

# ---------------------------------------------------------------------------
# residual allowance (Lemma E12): rho_i = sum_r a_ir sum_t (max_c v_irct - min_c v_irct)
# ---------------------------------------------------------------------------
.ed_rho <- function(P, terms) {
  v <- array(0, c(P$p0, P$T, P$R, P$C))
  for (tm in terms) {
    st <- tm$units; t <- (st - 1L) %/% P$p0 + 1L; i <- (st - 1L) %% P$p0 + 1L
    for (j in seq_along(st)) v[i[j], t[j], tm$rows, tm$cols] <- v[i[j], t[j], tm$rows, tm$cols] + tm$wd
  }
  rho <- numeric(P$p0)
  for (r in seq_len(P$R)) for (t in seq_len(P$T)) {
    vv <- matrix(v[, t, r, ], P$p0, P$C)
    rho <- rho + P$A0[, r] * (apply(vv, 1, max) - apply(vv, 1, min))
  }
  rho
}
.ed_allowance <- function(P, reps, type, a_E) {
  # code review issue 5: Lemma E12 bounds phi(n) - phi~ for the representation whose phi~ is targeted.
  # The model-implied forms (.ed_model_target, .ed_gapW, E25) use representation 1, so the allowance does too.
  reps <- reps[1]
  rhos <- lapply(reps, function(rp) .ed_rho(P, rp$terms))
  S <- vapply(rhos, function(r) 0.5 * sqrt(sum(r^2)), 0)
  j <- which.min(S); rho <- rhos[[j]]
  list(value = if (type == "popoviciu") stats::qnorm(1 - a_E / 2) * S[j] else sqrt(0.5 * log(2 / a_E) * sum(rho^2)),
       S = S[j], lindeberg = if (sum(rho^2) > 0) max(rho^2) / sum(rho^2) else NA_real_, rho = rho)
}

# ---------------------------------------------------------------------------
# intervals
# ---------------------------------------------------------------------------
.ed_stJ <- function(mb) list(beta = mb$beta, X = mb$X, V = mb$V, nbeta = mb$nbeta, nz = mb$nz)

.ed_proj <- function(mb, cons, tgm, alpha, ci, Zs, ndraw) {
  stJ <- .ed_stJ(mb)
  if (ci == "conservative") return(.rxc_ci_conservative(stJ, cons, tgm, alpha, Zs)$ci)
  as.vector(attr(.rxc_ci_clr(stJ, cons, tgm, alpha / 2, alpha / 2, Zs, ndraw, proj_alpha = alpha), "projection"))
}

.ed_plugin <- function(mb, cons, tgm) {
  pl <- .rxc_plugin(cons, tgm, mb$beta)
  if (isTRUE(pl$empty)) return(c(NA_real_, NA_real_))
  if (isTRUE(pl$unbounded)) return(c(-Inf, Inf))
  c(pl$lower, pl$upper)
}

.ed_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

.ed_model_status <- function(mb, ci, has_rows) {
  thm <- if (ci == "conservative") "P3.5(a)" else "P3.5(b2)(ii)"
  if (mb$rung == "FD") return(sprintf("asymptotic guarantee (%s for the model DiD, Prop. E23(a)-(b); Lemma E12 allowance), pooled FD", thm))
  if (mb$grouped && has_rows) return("no coverage theorem for model rows on changes under grouped rungs (R16 open)")
  if (mb$grouped) return(sprintf("%s: asymptotic guarantee only with few, large groups (P3 regime R1)", thm))
  sprintf("asymptotic guarantee (%s; Prop. E11 per-period pooled LCE; Lemma E12 allowance)", thm)
}

# E9(a) closed form for R = 2: |Delta kappa| <= q_min, model DiD in [I - |s| q_min, I + |s| q_min] (plug-in)
.ed_e9_closed <- function(P, mb, tgm, d) {
  if (P$R != 2) return(NULL)
  st <- mb$blocks[[d$k]]; e <- function(j) { v <- numeric(2); v[j] <- 1; v }
  rate <- function(r, v) .ed_pad(mb, d$k, .rxc_rate(st, 1, d$c, r, v))
  q <- vapply(seq_len(nrow(d$A)), function(j) {
    x <- d$A[j, ]
    O <- .ed_fsum(rate(d$rs[1], x), .ed_fscale(rate(d$rs[2], x), -1))
    L <- .ed_fsum(.ed_fscale(rate(1, e(d$rs[1])), x[1]), .ed_fscale(rate(2, e(d$rs[1])), x[2]),
                  .ed_fscale(rate(1, e(d$rs[2])), -x[1]), .ed_fscale(rate(2, e(d$rs[2])), -x[2]))
    (sum((O$coef + L$coef) * mb$beta) + O$const + L$const) / 2 }, 0)
  zk <- mb$zoff[d$k] + seq_len(mb$blocks[[d$k]]$nz)
  zc <- zk[d$c]
  if (d$c > length(zk)) return(NULL)                          # the last column's Delta k is minus the sum: not shown
  sl <- tgm$gamma[zc]
  if (any(tgm$gamma[-zc] != 0)) return(NULL)                  # target spans other blocks/columns
  I <- sum(tgm$cvec * mb$beta) + tgm$c0
  qm <- min(q)
  list(q_min = qm, s = sl, interval = if (qm < 0) c(NA_real_, NA_real_) else c(I - abs(sl) * qm, I + abs(sl) * qm), I = I)
}

# E8(d): sharp FD set (2 x 2, one pair) at the empirical law, and FD4's outer set
.ed_e8 <- function(P, mb, tgm) {
  if (P$R != 2 || P$C != 2 || length(mb$blocks) != 1) return(NULL)
  pr <- mb$pairs[[1]]; st <- mb$blocks[[1]]
  x <- P$A0[, 1] / P$N0
  Tu <- P$Bt[[pr[1]]][, 1] / P$N0; Tt <- P$Bt[[pr[2]]][, 1] / P$N0
  Lf <- function(t, x) pmax(0, (t - 1 + x) / x); Uf <- function(t, x) pmin(1, t / x)
  key <- format(x, digits = 15); ux <- unique(key)
  one <- function(outer) {
    cons <- .ed_cons_new(mb)
    for (k in ux) {
      j <- key == k; xx <- x[j][1]
      if (xx <= 0 || xx >= 1) next
      if (outer) { tu <- mean(Tu[j]); tt <- mean(Tt[j]); lo <- Lf(tt, xx) - Uf(tu, xx); hi <- Uf(tt, xx) - Lf(tu, xx) }
      else { lo <- mean(Lf(Tt[j], xx)) - mean(Uf(Tu[j], xx)); hi <- mean(Uf(Tt[j], xx)) - mean(Lf(Tu[j], xx)) }
      f <- .ed_pad(mb, 1, .rxc_rate(st, 1, 1, 1, c(xx, 1 - xx)))          # Delta m_v(x)
      cons <- .ed_cons_add(cons, list(coef = f$coef, const = f$const - lo, a = f$a), "E8")
      cons <- .ed_cons_add(cons, list(coef = -f$coef, const = hi - f$const, a = -f$a), "E8")
    }
    .ed_plugin(mb, cons, tgm)
  }
  list(sharp = one(FALSE), outer = one(TRUE), n_support = length(ux))
}

.ed_model_block <- function(fit, spec, alpha) {
  if (!inherits(spec, "edModel")) stop("'model' must be fd(), lce_periods() or lce_levels_change()")
  if (identical(spec$rung, "LCEeps")) return(.lce_eps_block(fit, spec, alpha))
  if (!is.null(fit$P$margins) && fit$P$margins == "estimated")
    warning("model layer with margins stated as estimated (G7): R26 on changes is deferred; compositions are taken as exact", call. = FALSE)
  P <- fit$P; mb <- .ed_mb(P, spec)
  meta <- list(r = .ed_rowidx(P, fit$defaults$row, "row"), s = .ed_versus(P, fit$defaults$versus, .ed_rowidx(P, fit$defaults$row, "row")),
               col = fit$defaults$col, didS = list())
  rr <- if (mb$rung == "FD") .ed_fd_rows(P, mb, spec$rows, meta) else .ed_lce_rows(P, mb, spec$rows, meta)
  base <- .ed_cons_bind(.ed_mb_validity(mb), rr$cons)
  systems <- .ed_mb_systems(base, rr$unions)
  a_est <- 0.9 * alpha; a_E <- 0.1 * alpha
  Zs <- .ed_with_seed(spec$seed, .rxc_zdraws(mb$V, 4000))
  belief <- !is.null(fit$restrictions) && any(fit$restrictions$role == "belief")
  out <- list(); extra <- list()
  for (j in seq_along(fit$targets)) {
    tg <- fit$targets[[j]]; tb <- fit$table[j, ]
    row <- data.frame(target = tg$label, count_lower = tb$lower, count_upper = tb$upper, plugin_lower = NA_real_, plugin_upper = NA_real_,
                      model_lower = NA_real_, model_upper = NA_real_, lower = NA_real_, upper = NA_real_, disjoint = NA,
                      allowance = NA_real_, S_phi = NA_real_, lindeberg = NA_real_, status = NA_character_, note = "", stringsAsFactors = FALSE)
    f <- .ed_model_target(P, mb, tg)
    if (!length(f)) { row$note <- attr(f, "reason"); out[[j]] <- row; next }
    if (!tg$design) { row$note <- "count layer only: the unit set uses the outcomes (Lemma E12 needs design coefficients)"; out[[j]] <- row; next }
    tgm <- list(cvec = f$coef, c0 = f$const, gamma = f$a)
    al <- .ed_allowance(P, tg$reps, spec$allowance, a_E)
    ivs <- list(); pls <- list()
    for (sy in systems) {
      iv <- .ed_with_seed(spec$seed, .ed_proj(mb, sy$cons, tgm, a_est, spec$ci, Zs, spec$ndraw))
      ivs[[length(ivs) + 1]] <- iv; pls[[length(pls) + 1]] <- .ed_plugin(mb, sy$cons, tgm)
    }
    mpcs <- .rxc_union(lapply(ivs, function(v) v + c(-1, 1) * al$value), tol = 0)
    ppcs <- .rxc_union(pls, tol = 0)
    if (!is.null(ppcs)) { row$plugin_lower <- min(ppcs[, 1]); row$plugin_upper <- max(ppcs[, 2]) }
    row$allowance <- al$value; row$S_phi <- al$S; row$lindeberg <- al$lindeberg
    row$status <- .ed_model_status(mb, spec$ci, length(spec$rows) > 0)
    if (is.null(mpcs)) { row$note <- "the model interval is empty: the model layer rejects the model (+ its rows) by itself"; row$disjoint <- TRUE }
    else {
      row$model_lower <- min(mpcs[, 1]); row$model_upper <- max(mpcs[, 2])
      I <- .ec_intersect(fit$solutions[[j]]$pieces, mpcs)
      row$disjoint <- is.null(I) && !is.null(fit$solutions[[j]]$pieces)
      if (!is.null(I)) { row$lower <- min(I[, 1]); row$upper <- max(I[, 2]) }
    }
    if (isTRUE(tg$given_set)) row$note <- "the unit set must be fixed from the design (Lemma E12)"
    ex <- list(model_pieces = mpcs, branches = data.frame(branch = vapply(systems, function(s) s$label, ""),
                                                          lower = vapply(ivs, `[`, 0, 1), upper = vapply(ivs, `[`, 0, 2)))
    for (d in rr$meta$didS) { e9 <- .ed_e9_closed(P, mb, tgm, d); if (!is.null(e9)) ex$e9 <- e9 }
    if (isTRUE(spec$sharp) && mb$rung == "FD") ex$e8 <- .ed_e8(P, mb, tgm)
    out[[j]] <- row; extra[[j]] <- ex
  }
  tab <- do.call(rbind, out); rownames(tab) <- NULL
  list(table = tab, extra = extra, rung = mb$rung, spec = spec, row_labels = rr$labels, alpha = alpha, belief = belief,
       mb = mb, base = base, unions = rr$unions, vcov_type = mb$vcov_type, grouped = mb$grouped,
       identity_rows = c("margins", if (!is.null(fit$restrictions)) fit$restrictions$restriction[fit$restrictions$role == "identity"]),
       belief_rows = if (!is.null(fit$restrictions)) fit$restrictions$restriction[fit$restrictions$role == "belief"] else character(0))
}

.ed_print_model <- function(m, digits = 3) {
  tb <- m$table
  cat(sprintf("\nModel layer (%s; DESCRIPTIVE about model-implied and, through the allowance, actual DiDs): level %g%%\n",
              if (m$rung == "FD") "first-difference LCE, 2e Sec. 3" else if (m$rung == "LCEeps")
                sprintf("two-election LCE, change exact, base level within epsilon = %g of affine, Prop. E36", m$epsilon) else
                sprintf("per-period LCE (%s), Prop. E11", m$spec$model), 100 * (1 - m$alpha)))
  if (identical(m$rung, "LCEeps") && m$epsilon >= 0.5)
    cat("  epsilon >= 1/2: the curvature rows are slack; this is change LCE plus the base-level link (inside fd(), Prop. E36(c))\n")
  if (length(m$row_labels)) cat("  model rows (object = \"model\"):", paste(m$row_labels, collapse = "; "), "\n")
  cat("  validity rows: |Delta m(x)| <= 1 at the hull vertices (FD) or rates in [0, 1] (LCE)\n")
  cat("  labels: plugin = the model set at the estimate (an estimate, not an interval); actual_target_interval = CONFIDENCE INTERVAL\n",
      "  for the actual DiD (model interval + allowance); count_set = DESCRIPTIVE, realised counts; intersection = both combined (Thm E13)\n", sep = "")
  out <- data.frame(target = tb$target, plugin = .ec_fmt(tb$plugin_lower, tb$plugin_upper, digits),
                    actual_target_interval = .ec_fmt(tb$model_lower, tb$model_upper, digits),
                    count_set = .ec_fmt(tb$count_lower, tb$count_upper, digits),
                    intersection = ifelse(tb$disjoint %in% TRUE, "disjoint", .ec_fmt(tb$lower, tb$upper, digits)), stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  for (j in seq_len(nrow(tb))) {
    if (nzchar(tb$note[j])) cat(sprintf("  %s: %s\n", tb$target[j], tb$note[j]))
    if (!is.na(tb$status[j])) cat(sprintf("  %s: %s interval + %s allowance %.4g (Lemma E12; S = %.4g); %s\n", tb$target[j],
                                          m$spec$ci, m$spec$allowance, tb$allowance[j], tb$S_phi[j], tb$status[j]))
    ex <- if (j <= length(m$extra)) m$extra[[j]] else NULL
    if (!is.null(ex$e9)) cat(sprintf("  E9(a) closed form (plug-in): q_min = %.4g, model DiD in %s\n", ex$e9$q_min, .ec_fmt(ex$e9$interval[1], ex$e9$interval[2], digits)))
    if (isTRUE(m$spec$sharp) && is.null(ex$e8) && m$rung == "FD")
      cat("  E8 not computed: the plug-in sharp set is for one pair in 2 x 2 (pair-by-pair stacking over several pairs would be OUTER, E8(e))\n")
    if (!is.null(ex$e8)) cat(sprintf("  E8 sharp FD set at the empirical law (plug-in, %d support points; inference open, E23(c)): %s; FD4 outer: %s\n",
                                     ex$e8$n_support, .ec_fmt(ex$e8$sharp[1], ex$e8$sharp[2], digits), .ec_fmt(ex$e8$outer[1], ex$e8$outer[2], digits)))
    if (!is.null(ex$branches) && nrow(ex$branches) > 1) cat(sprintf("  %s: union of %d branch intervals (Prop. E24; no multiplicity correction)\n", tb$target[j], nrow(ex$branches)))
  }
  cat("  plug-in: the model set at the estimate (not an interval). Cross-layer rule (Thm E13, 6a R27): count rows and model rows never share an LP.\n")
  if (!m$belief) cat("  ", .ec_txt("cov_id"), "\n", sep = "") else cat("  ", .ec_txt("cov_bel"), "\n", sep = "")
  if (any(tb$disjoint %in% TRUE)) cat("  ", if (m$belief) .ec_txt("rej_bel") else .ec_txt("rej_id"), "\n", sep = "")
}

# ---------------------------------------------------------------------------
# eiATT with model-based pre-trends: Prop. E25
# ---------------------------------------------------------------------------
.att_model <- function(fit, cs, trend, Lambda, alpha, tgt) {
  m <- fit$model
  if (identical(m$rung, "LCEeps")) stop("lce_levels_change() is for eiATT(trend = within_pt(), lce = TRUE) (Props. E35-E36); pretrends = \"model\" uses fd() or lce_periods()")
  if (is.null(m)) stop("pretrends = \"model\" needs an eiChange fit with a model layer (model = fd() or lce_periods())")
  if (!trend$kind %in% c("dial", "rm", "sd"))
    stop("pretrends = \"model\" is for the trend dial, relative magnitudes and smoothness (Prop. E25)")
  if (!cs$has_comp) stop("internal: comparison required")
  P <- fit$P; P$cc <- cs$cc; mb <- m$mb; spec <- m$spec
  levels_ok <- mb$rung == "LCE"                   # (H2): rate limits only with a rung that models levels
  b0 <- cs$base; posts <- cs$posts; K <- length(posts)
  S_all <- seq_len(P$p0)
  gW <- .ed_gapW(P, S_all, cs$r, cs$s)
  # actual quantities q: list(form = model-implied joint form, reps = count terms for the allowance)
  gap_reps <- function(u, t) {
    gt <- .ed_gap_terms(P, t, S_all, cs$r, cs$s, P$cc, 1); gu <- .ed_gap_terms(P, u, S_all, cs$r, cs$s, P$cc, -1)
    list(list(terms = c(gt$terms, gu$terms)), list(terms = c(gt$terms2, gu$terms2)))
  }
  qs <- list()
  addq <- function(name, form, reps) { qs[[name]] <<- list(form = form, reps = reps); length(qs) }
  for (k in seq_len(K)) addq(paste0("D", k), .ed_model_change(mb, b0, posts[k], P$cc, gW), gap_reps(b0, posts[k]))
  if (trend$kind == "rm") {
    if (!length(cs$pre_pairs)) stop("rel_mag() needs a pre-period")
    for (q in seq_along(cs$pre_pairs)) { pp <- cs$pre_pairs[[q]]
      addq(paste0("pre", q), .ed_model_change(mb, pp[1], pp[2], P$cc, gW), gap_reps(pp[1], pp[2])) }
  }
  if (trend$kind == "sd") for (e in seq_len(b0 - 1))
    addq(paste0("sdpre", e), .ed_model_change(mb, b0, e, P$cc, gW), gap_reps(b0, e))
  if (levels_ok) for (k in seq_len(K)) {
    t <- posts[k]
    WBs <- -gW; WBs[, cs$r] <- 0                                   # comparison rate weights (positive)
    WBr <- .ed_rateW(P, S_all, cs$r)
    th <- .ed_fsum(.ed_lce_lin(mb, b0, P$cc, WBr), .ed_lce_lin(mb, t, P$cc, WBs), .ed_fscale(.ed_lce_lin(mb, b0, P$cc, WBs), -1))
    nm <- P$num
    Bs_terms <- function(e, sg = 1) if (is.null(cs$s))
      list(.ed_term(P, e, S_all, setdiff(seq_len(P$R), cs$r), P$cc, sg / (sum(nm(P$N0)) - sum(nm(P$A0[, cs$r]))))) else
      list(.ed_term(P, e, S_all, cs$s, P$cc, sg / sum(nm(P$A0[, cs$s]))))
    Br_terms <- list(.ed_term(P, b0, S_all, cs$r, P$cc, 1 / sum(nm(P$A0[, cs$r]))))
    addq(paste0("thPT", k), th, list(list(terms = c(Br_terms, Bs_terms(t), Bs_terms(b0, -1)))))
    addq(paste0("Bs", k), .ed_lce_lin(mb, t, P$cc, WBs), list(list(terms = Bs_terms(t))))
  }
  nq <- length(qs)
  a_proj <- 0.9 * alpha; a_q <- 0.1 * alpha / nq
  tq <- vapply(qs, function(q) .ed_allowance(P, q$reps, spec$allowance, a_q)$value, 0)
  # LP variables: z = (k [mb$nz], theta_k, S_k, b_k [3K], e_q [nq])
  nzx <- mb$nz + 3 * K + nq
  ix <- function(k, w) mb$nz + 3 * (k - 1) + switch(w, theta = 1, S = 2, b = 3)
  iq <- function(name) mb$nz + 3 * K + match(name, names(qs))
  ext <- function(f) list(coef = f$coef, const = f$const, a = c(f$a, numeric(3 * K + nq)))
  zv <- function(pairs, const = 0) { a <- numeric(nzx); for (p in pairs) a[p[1]] <- a[p[1]] + p[2]
    list(coef = numeric(mb$nbeta), const = const, a = a) }
  actual <- function(name) .ed_fsum(ext(qs[[name]]$form), zv(list(c(iq(name), 1))))     # q~ + e_q
  base <- m$base; base$A <- cbind(base$A, matrix(0, nrow(base$A), 3 * K + nq))
  add <- function(cons, f, kind) .ed_cons_add(cons, f, kind)
  cons0 <- base
  for (nmq in names(qs)) { cons0 <- add(cons0, zv(list(c(iq(nmq), -1)), tq[[nmq]]), "box"); cons0 <- add(cons0, zv(list(c(iq(nmq), 1)), tq[[nmq]]), "box") }
  for (k in seq_len(K)) {
    cons0 <- add(cons0, zv(list(c(ix(k, "S"), -1)), Lambda), "spill"); cons0 <- add(cons0, zv(list(c(ix(k, "S"), 1)), Lambda), "spill")
    if (levels_ok) {
      th <- .ed_fsum(actual(paste0("thPT", k)), zv(list(c(ix(k, "b"), 1), c(ix(k, "S"), -1))))   # theta = thetaPT + b - S
      cons0 <- add(cons0, th, "rate"); cons0 <- add(cons0, .ed_fsum(.ed_fscale(th, -1), zv(list(), 1)), "rate")
      bs <- .ed_fsum(actual(paste0("Bs", k)), zv(list(c(ix(k, "S"), -1))))
      cons0 <- add(cons0, bs, "rate"); cons0 <- add(cons0, .ed_fsum(.ed_fscale(bs, -1), zv(list(), 1)), "rate")
    }
  }
  delta <- function(k) if (k == 0) zv(list()) else zv(list(c(ix(k, "b"), 1)))
  gbr <- list()
  if (trend$kind == "dial") {
    cc <- cons0
    for (k in seq_len(K)) { cc <- add(cc, zv(list(c(ix(k, "b"), -1)), trend$M), "dial"); cc <- add(cc, zv(list(c(ix(k, "b"), 1)), trend$M), "dial") }
    gbr <- list(list(cons = cc, label = ""))
  } else if (trend$kind == "rm") {
    for (q in seq_along(cs$pre_pairs)) for (sg in c(1, -1)) {
      cc <- cons0; d <- .ed_fscale(actual(paste0("pre", q)), sg * trend$Mbar)
      for (k in seq_len(K)) { dd <- .ed_fsum(delta(k), .ed_fscale(delta(k - 1), -1))
        cc <- add(cc, .ed_fsum(d, .ed_fscale(dd, -1)), "rm"); cc <- add(cc, .ed_fsum(d, dd), "rm") }
      pp <- cs$pre_pairs[[q]]
      gbr[[length(gbr) + 1]] <- list(cons = cc, label = sprintf("A at %s -> %s, sign %s", P$elections[pp[1]], P$elections[pp[2]], if (sg > 0) "+" else "-"))
    }
  } else {
    dl <- lapply(seq_len(P$T), function(e) if (e < b0) actual(paste0("sdpre", e)) else if (e == b0) zv(list()) else {
      k <- match(e, posts); if (is.na(k)) NULL else delta(k) })
    cc <- cons0
    for (j in seq_len(max(0, P$T - 2)) + 1) {          # code review issue 8
      if (is.null(dl[[j - 1]]) || is.null(dl[[j]]) || is.null(dl[[j + 1]])) next
      s2 <- .ed_fsum(dl[[j + 1]], .ed_fscale(dl[[j]], -2), dl[[j - 1]])
      cc <- add(cc, .ed_fsum(.ed_fscale(s2, -1), zv(list(), trend$M)), "sd"); cc <- add(cc, .ed_fsum(s2, zv(list(), trend$M)), "sd")
    }
    gbr <- list(list(cons = cc, label = ""))
  }
  stJ <- list(beta = mb$beta, X = mb$X, V = mb$V)
  Zs <- .ed_with_seed(spec$seed, .rxc_zdraws(mb$V, 4000))
  unions_ext <- lapply(m$unions, function(u) { f <- function(cs0) { cs0$A <- cbind(cs0$A, matrix(0, nrow(cs0$A), 3 * K + nq)); cs0 }
    list(plus = f(u$plus), minus = f(u$minus), id = u$id) })
  res <- list()
  for (t in tgt) {
    k <- match(t, posts)
    tf <- .ed_fsum(actual(paste0("D", k)), zv(list(c(ix(k, "b"), -1), c(ix(k, "S"), 1))))   # ATT = D - b + S
    tgm <- list(cvec = tf$coef, c0 = tf$const, gamma = tf$a)
    ivs <- list(); labs <- character(0)
    for (g in gbr) for (sy in .ed_mb_systems(g$cons, unions_ext)) {
      iv <- .ed_with_seed(spec$seed, as.vector(attr(.rxc_ci_clr(stJ, sy$cons, tgm, a_proj / 2, a_proj / 2, Zs, spec$ndraw, proj_alpha = a_proj), "projection")))
      ivs[[length(ivs) + 1]] <- iv; labs <- c(labs, paste(c(g$label, sy$label)[nzchar(c(g$label, sy$label))], collapse = "; "))
    }
    # an ATT of a rate lies in [-1, 1]: clipping the projection interval there is valid
    # an interval wholly outside [-1, 1] is a rejection at the stated level: EMPTY (code review issue 13(a))
    ivc <- lapply(ivs, function(v) if (any(is.na(v)) || v[1] > 1 || v[2] < -1) c(NA_real_, NA_real_) else pmin(pmax(v, -1), 1))
    pcs <- .rxc_union(ivc, tol = 0)
    res[[length(res) + 1]] <- list(post = P$elections[t], pieces = pcs, lower = if (is.null(pcs)) NA_real_ else min(pcs[, 1]),
                                   upper = if (is.null(pcs)) NA_real_ else max(pcs[, 2]),
                                   branches = data.frame(branch = labs, lower = vapply(ivs, `[`, 0, 1), upper = vapply(ivs, `[`, 0, 2)))
  }
  list(results = res, outer = !levels_ok, rung = mb$rung, radii = tq, alpha_proj = a_proj, alpha_q = a_q,
       status = if (mb$grouped && length(m$row_labels)) "no coverage theorem for model rows under grouped rungs (R16 open)" else
         if (mb$grouped) "grouped rung without model rows: asymptotic guarantee only with few, large groups (P3 regime R1)" else
         sprintf("asymptotic: coverage >= 1 - %.3g - %d x %.3g (Prop. E25(a): one projection event for the stacked beta, union bound over the boxes)", a_proj, nq, a_q))
}

.att_print_model <- function(x, digits = 3) {
  m <- x$model
  f <- function(v) if (is.null(v) || is.na(v[1])) "empty" else sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), v[1], v[2])
  cat(sprintf("  model-based pre-trends (Prop. E25): rung %s; one joint projection region (P3.5(b2)(ii)) for the stacked beta;\n", m$rung))
  cat(sprintf("  box variables for the residual of each actual quantity (radii %s; Lemma E12, %s)\n",
              paste(sprintf("%s %.3g", names(m$radii), m$radii), collapse = ", "), x$fit$model$spec$allowance))
  cat("  inference:", m$status, "\n  (intervals clipped to [-1, 1], the logical range of an ATT of a rate)\n")
  if (m$outer) cat("  OUTER: under FD alone the levels are free, so the rate limits are dropped (Prop. E25 (H2))\n")
  for (e in m$results) {
    cat(sprintf("\nPost-election %s\n  [CONFIDENCE SET, realised ATT (Prop. E25)%s] ATT: %s\n", e$post, if (m$outer) ", OUTER" else "",
                if (is.null(e$pieces)) "EMPTY (rejects the conjunction at the stated level)" else paste(apply(e$pieces, 1, f), collapse = " U ")))
    if (nrow(e$branches) > 1) cat(sprintf("  union of %d branch intervals (Prop. E24; never a selected branch)\n", nrow(e$branches)))
  }
  cat("\n  Prop. E22: pre-periods carry no information on spillover onto the comparison group; Lambda is its own dial.\n")
}
