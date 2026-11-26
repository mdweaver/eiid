# Changes and DiD, descriptive layer: eiChange(), closed forms, the joint (event-study) set,
# print and witness. Theory: docs/extension/PHASE2e_did_general.md (2e) Secs 2-4, 7.

# ---------------------------------------------------------------------------
# closed forms for one-vs-rest DiDs with W1 rows on the same (r, c, u, t, S): 2e Props E5, E26
# ---------------------------------------------------------------------------
# Every W1 sign row on changes is a half-line in Z = Y_rct(S) - Y_rcu(S) (Lemma E4(b)), CR the
# union of two such systems, and D^r_c is affine increasing in Z. Returns NULL when the target is
# not a one-vs-rest DiD (R = 2: the named pair is the rest) or some row is not of this kind.
.ed_closed_ovr <- function(P, tg, specs) {
  if (!identical(tg$kind, "did") || !(is.null(tg$s) || P$R == 2)) return(NULL)
  th <- .ed_ovr_change(P, tg$S, tg$r, tg$c, tg$u, tg$t)
  if (is.null(th)) return(NULL)
  same <- function(tp) !is.null(tp) && tp$r == tg$r && tp$c == tg$c && tp$u == tg$u && tp$t == tg$t && setequal(tp$S, tg$S)
  for (s in specs) if (!same(s$tuple)) return(NULL)
  # Z-intervals per branch (exact); each W1 row has one finite side
  side <- function(row) list(lo = row$loq, hi = row$hiq, lod = row$lo, hid = row$hi)
  base <- list(lo = th$Zlo, hi = th$Zhi)
  cut <- function(iv, row) {
    sd <- side(row)
    if (!is.null(sd$lo)) iv$lo <- if (iv$lo >= sd$lo) iv$lo else sd$lo
    if (!is.null(sd$hi)) iv$hi <- if (iv$hi <= sd$hi) iv$hi else sd$hi
    if (P$integer) return(iv)
    if (is.finite(sd$lod)) iv$lo <- max(iv$lo, sd$lod)
    if (is.finite(sd$hid)) iv$hi <- min(iv$hi, sd$hid)
    iv
  }
  com <- base; crs <- list()
  for (s in specs) {
    if (length(s$rows)) for (rw in s$rows) com <- cut(com, rw)
    if (!is.null(s$union)) crs[[length(crs) + 1]] <- s$union
  }
  combos <- if (length(crs)) as.matrix(expand.grid(rep(list(1:2), length(crs)))) else matrix(0L, 1, 0)
  br <- list()
  for (k in seq_len(nrow(combos))) {
    iv <- com
    for (j in seq_along(crs)) for (rw in crs[[j]][[if (combos[k, j] == 1) "plus" else "minus"]]) iv <- cut(iv, rw)
    br[[k]] <- iv
  }
  toD <- function(iv) c(as.double(th$Dof(iv$lo)), as.double(th$Dof(iv$hi)))
  cont <- lapply(br, function(iv) if (as.double(iv$lo) <= as.double(iv$hi) + 1e-12) toD(iv) else c(NA_real_, NA_real_))
  int <- NULL; Zint <- NULL
  if (P$integer) {
    Zint <- lapply(br, function(iv) { lo <- .ec_ceil(iv$lo); hi <- floor(iv$hi); if (lo <= hi) c(lo, hi) else NULL })
    int <- lapply(Zint, function(z) if (is.null(z)) c(NA_real_, NA_real_) else c(as.double(th$Dof(z[1])), as.double(th$Dof(z[2]))))
  }
  list(th = th, cont = cont, int = int, Zint = Zint,
       ER = as.double(th$ER), NM = as.double(th$NM), gamma = as.double(th$gamma),
       Dlo = as.double(th$Dof(th$Zlo)), Dhi = as.double(th$Dof(th$Zhi)), step = as.double(th$kappa),
       labels = if (length(crs)) apply(combos, 1, function(z) paste(vapply(seq_along(z), function(j)
         paste0(crs[[j]]$id, if (z[j] == 1) "+" else "-"), ""), collapse = ", ")) else "")
}

# integer witness for a value Z of Y_rct(S) - Y_rcu(S): one total row per election (one cell, two
# disjoint unit sets: TU class (i) of 6a R19(b)), LP vertex checked exactly; MILP otherwise
.ed_witness_Z <- function(P, tg, th, Z) {
  nm <- P$num; S <- tg$S; r <- tg$r; c <- tg$c
  ddi <- function(e) { a <- P$A0[S, r]; b <- P$Bt[[e]][S, c]; Ni <- P$N0[S]
    c(sum(nm(pmax(0, a + b - Ni))), sum(nm(pmin(a, b)))) }
  lt <- ddi(tg$t); lu <- ddi(tg$u)
  Yt <- Z + lu[1]; if (Yt < lt[1]) Yt <- lt[1]; if (Yt > lt[2]) Yt <- lt[2]
  Yu <- Yt - Z
  rows <- list(.ec_row(P, "w", "total", .ed_vix(P, S, tg$t, r, c), rep(1, length(S)), Yt, Yt, cell = c(r, c), units = .ed_su(P, tg$t, S)),
               .ec_row(P, "w", "total", .ed_vix(P, S, tg$u, r, c), rep(1, length(S)), Yu, Yu, cell = c(r, c), units = .ed_su(P, tg$u, S)))
  rows <- lapply(rows, .ec_round_row)
  s <- .ec_optimise(P, rows, numeric(P$nv), 1, rd = TRUE)
  if (isTRUE(s$ok) && max(abs(s$x - round(s$x))) <= 1e-6 && .ec_check(P, rows, round(s$x), TRUE)) return(round(s$x))
  m <- .ec_milp(P, rows, tg, 1, P$time_limit)
  if (m$status == "optimal") m$x else NULL
}

# the integer range of 2e Prop. E26 (lattice of step kappa_r), in the format of .ec_range
.ed_e26_range <- function(P, tg, cf) {
  ok <- !vapply(cf$Zint, is.null, TRUE)
  cont_ok <- !vapply(cf$cont, function(v) is.na(v[1]), TRUE)
  out <- list(cont = any(cont_ok), milp = FALSE, class = "E26", unresolved = NULL, fallback = FALSE,
              lp = if (any(cont_ok)) c(min(vapply(cf$cont[cont_ok], `[`, 0, 1)), max(vapply(cf$cont[cont_ok], `[`, 0, 2))) else c(NA, NA),
              status = "integer-exact (E26 lattice)")
  out$int_ok <- any(ok)
  if (nzchar(cf$labels[1])) out$branches <- data.frame(branch = cf$labels, continuous = cont_ok,
                                                      lower = vapply(cf$int, `[`, 0, 1), upper = vapply(cf$int, `[`, 0, 2),
                                                      stringsAsFactors = FALSE)
  if (!any(ok)) { out$lower <- out$upper <- NA_real_; out$pieces <- NULL; out$x <- list(NULL, NULL); return(out) }
  out$pieces <- .rxc_union(cf$int[ok], tol = 1e-12)
  jl <- which(ok)[which.min(vapply(cf$int[ok], `[`, 0, 1))]; ju <- which(ok)[which.max(vapply(cf$int[ok], `[`, 0, 2))]
  out$lower <- cf$int[[jl]][1]; out$upper <- cf$int[[ju]][2]
  out$x <- list(.ed_witness_Z(P, tg, cf$th, cf$Zint[[jl]][1]), .ed_witness_Z(P, tg, cf$th, cf$Zint[[ju]][2]))
  out$sys <- c(cf$labels[jl], cf$labels[ju])
  out
}

# closed form of 2e Prop. E6(b) for the bin-averaged DiD with the contrast beliefs on the same bins
.ed_closed_contrast <- function(P, tg, specs, dd) {
  if (!identical(tg$kind, "bin_did") || !length(specs)) return(NULL)
  same <- function(s) s$type == "contrast_change" && s$bins$c == tg$c && s$bins$u == tg$u && s$bins$t == tg$t &&
    setequal(s$bins$L, tg$L) && setequal(s$bins$H, tg$H)
  if (!all(vapply(specs, same, TRUE))) return(NULL)
  bel <- unique(unlist(lapply(specs, function(s) s$beliefs)))
  cap <- as.double(tg$cap)
  lo <- if ("mix" %in% bel) max(0, dd$lower) else dd$lower
  hi <- if ("rate" %in% bel) min(dd$upper, cap) else dd$upper
  why <- character(0)
  if (lo > hi) {
    if (all(c("mix", "rate") %in% bel) && as.double(tg$dtau) < 0) why <- c(why, "Dtau_H < Dtau_L: the two beliefs are jointly refuted")
    if ("mix" %in% bel && dd$upper < 0) why <- c(why, "Dbar+ < 0: the mix belief is refuted by the margins")
    if ("rate" %in% bel && cap < dd$lower) why <- c(why, "cap < Dbar-: the rate belief is refuted")
  }
  list(lower = lo, upper = hi, empty = lo > hi, why = why, cap = cap, lambda = as.double(tg$lambda),
       dtau = as.double(tg$dtau), beliefs = bel)
}

# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------
.ed_fit <- function(P, tl, specs, integer, binding) {
  built <- .ed_build_rows(P, specs)
  empty <- list(common = list(), unions = list())
  rows_out <- list(); sol <- list(); cfs <- list()
  integer_mode <- P$integer && integer != "lp"
  for (j in seq_along(tl)) {
    tg <- tl[[j]]
    dd <- .ec_range(P, empty, tg, integer, keep_x = FALSE)
    cf <- .ed_closed_ovr(P, tg, specs)
    cc <- .ed_closed_contrast(P, tg, specs, dd)
    if (!is.null(cf) && integer_mode) rg <- .ed_e26_range(P, tg, cf)
    else rg <- if (length(specs)) .ec_range(P, built, tg, integer) else dd
    bd <- if (isTRUE(binding) && length(specs)) .ec_binding(P, specs, built, tg, rg, integer) else c(NA, NA)
    st <- rg$status
    if (!isTRUE(rg$cont)) st <- "infeasible: the rows contradict the margins (continuous tables)"
    else if (identical(rg$int_ok, FALSE)) st <- "integer-infeasible: no whole-number table satisfies these rows"
    # closed form (continuous) as a cross-check of the LP: E5 (one-vs-rest rows), E6 (contrast)
    cl <- c(NA_real_, NA_real_); cform <- NA_character_
    if (!is.null(cf)) {
      pcs <- .rxc_union(Filter(function(v) !is.na(v[1]), cf$cont), tol = 1e-12)
      if (!is.null(pcs)) cl <- c(min(pcs[, 1]), max(pcs[, 2]))
      cform <- if (length(specs)) "E5 (W1 rows on changes)" else "E3(b)/E5(a) (margins only)"
    } else if (!is.null(cc)) {
      if (!cc$empty) cl <- c(cc$lower, cc$upper)
      cform <- "E6(b) (count contrast on changes)"
    }
    rows_out[[j]] <- data.frame(
      target = tg$label, denominator = tg$den, dd_lower = dd$lower, dd_upper = dd$upper,
      lower = rg$lower, upper = rg$upper, width = rg$upper - rg$lower,
      integer_status = st, tu_class = if (is.null(rg$class)) NA_character_ else rg$class,
      lp_lower = rg$lp[1], lp_upper = rg$lp[2], closed_lower = cl[1], closed_upper = cl[2], closed_form = cform,
      continuous_feasible = isTRUE(rg$cont), integer_feasible = rg$int_ok,
      n_pieces = if (is.null(rg$pieces)) 0L else nrow(rg$pieces),
      binding_lower = bd[1], binding_upper = bd[2], design = tg$design, change_target = isTRUE(tg$change),
      stringsAsFactors = FALSE)
    sol[[j]] <- list(x = rg$x, sys = rg$sys, pieces = rg$pieces, branches = rg$branches, fallback = isTRUE(rg$fallback))
    cfs[[j]] <- list(ovr = if (is.null(cf)) NULL else cf[c("ER", "NM", "gamma", "Dlo", "Dhi", "step")], contrast = cc)
  }
  tab <- do.call(rbind, rows_out); rownames(tab) <- NULL
  # cross-check: closed form = LP (continuous) wherever both exist
  dis <- !is.na(tab$closed_lower) & !is.na(tab$lp_lower) &
    (abs(tab$closed_lower - tab$lp_lower) > 1e-7 | abs(tab$closed_upper - tab$lp_upper) > 1e-7)
  if (any(dis)) warning("closed form and LP disagree for: ", paste(tab$target[dis], collapse = "; "),
                        " (please report; the LP is used)")
  rt <- if (length(specs)) data.frame(id = vapply(specs, function(s) s$id, ""), restriction = vapply(specs, function(s) s$label, ""),
                                      scope = vapply(specs, function(s) s$scope, ""), object = "actual",
                                      role = vapply(specs, function(s) s$role, ""),
                                      experimental = vapply(specs, function(s) isTRUE(s$experimental), TRUE),
                                      stringsAsFactors = FALSE) else NULL
  list(table = tab, targets = tl, specs = specs, built = built, restrictions = rt, solutions = sol, closed = cfs)
}

# ---------------------------------------------------------------------------
# the joint (event-study) set: 2e Prop. E3(c)
# ---------------------------------------------------------------------------
.ed_joint_feasible <- function(P, built, tgs, values, tol = 1e-9) {
  eq <- lapply(seq_along(tgs), function(j) {
    tg <- tgs[[j]]; idx <- which(tg$obj != 0)
    v <- values[j] - as.double(tg$const); sc <- tol * max(1, abs(v))
    rw <- .ec_row(P, "joint", "gen", idx, tg$obj[idx], v - sc, v + sc)
    rw$coefq <- NULL; rw$loq <- NULL; rw$hiq <- NULL
    rw
  })
  U <- length(built$unions)
  combos <- if (U) as.matrix(expand.grid(rep(list(1:2), U))) else matrix(0L, 1, 0)
  for (k in seq_len(nrow(combos))) {
    rows <- c(built$common, eq)
    for (u in seq_len(U)) rows <- c(rows, built$unions[[u]][[if (combos[k, u] == 1) "plus" else "minus"]])
    if (isTRUE(.ec_optimise(P, rows, numeric(P$nv), 1)$ok)) return(TRUE)
  }
  FALSE
}

#' Is a vector of target values jointly attainable?
#'
#' \code{eiJoint()} tests whether given values of several targets of an \code{\link{eiChange}}
#' fit are attained by \emph{one} family of realised tables (continuous tables) consistent with
#' the margins and every count row (any branch of a CR union). The joint set of event-study
#' DiDs is not the box of their separate ranges (2e Prop. E3(c)): a corner of the box can be
#' infeasible because the DiDs share the base election's table.
#' @param fit an \code{eiChange} fit.
#' @param values numeric vector of target values (named by target label, or in target order).
#' @param targets optional target labels or indices that \code{values} refer to.
#' @param rows \code{TRUE} (default) to impose the fit's count rows, \code{FALSE} for margins only.
#' @return \code{TRUE} or \code{FALSE}.
#' @export
eiJoint <- function(fit, values, targets = NULL, rows = TRUE) {
  if (!inherits(fit, "eiChange")) stop("'fit' must be an eiChange fit")
  labs <- fit$table$target
  js <- if (!is.null(targets)) { if (is.numeric(targets)) targets else match(targets, labs) } else
    if (!is.null(names(values))) match(names(values), labs) else seq_along(values)
  if (any(is.na(js)) || length(js) != length(values)) stop("unknown targets in eiJoint()")
  b <- if (rows) fit$built else list(common = list(), unions = list())
  .ed_joint_feasible(fit$P, b, fit$targets[js], as.numeric(values))
}

.ed_event_study <- function(P, tl, built, base) {
  es <- which(vapply(tl, function(tg) identical(tg$kind, "did") && tg$u == base, TRUE))
  if (!length(es)) return(NULL)
  tg1 <- tl[[es[1]]]
  # per-election gap ranges G_t (margins only; 2e Prop. E3(b))
  G <- t(vapply(seq_len(P$T), function(t) {
    g <- .ed_gap_terms(P, t, tg1$S, tg1$r, tg1$s, tg1$c, 1)
    tgg <- list(terms = g$terms, const = P$num(0)); tgg$obj <- .ec_obj(P, g$terms)
    r <- .ec_range(P, list(common = list(), unions = list()), tgg, "lp", keep_x = FALSE)
    c(r$lower, r$upper) }, c(0, 0)))
  out <- list(targets = es, G = data.frame(election = P$elections, G_lower = G[, 1], G_upper = G[, 2], stringsAsFactors = FALSE))
  if (length(es) >= 2) {
    a <- es[1]; b <- es[2]
    same <- function(k) tl[[k]]$r == tg1$r && identical(tl[[k]]$s, tg1$s) && tl[[k]]$c == tg1$c && setequal(tl[[k]]$S, tg1$S)
    if (same(b)) {
      rng <- function(k, bl) { r <- .ec_range(P, bl, tl[[k]], "lp", keep_x = FALSE); c(r$lower, r$upper) }
      e0 <- list(common = list(), unions = list())
      ra <- rng(a, e0); rb <- rng(b, e0)
      out$corner <- list(targets = c(tl[[a]]$label, tl[[b]]$label), values = c(ra[2], rb[1]),
                         margins_only = .ed_joint_feasible(P, e0, tl[c(a, b)], c(ra[2], rb[1])),
                         base_width = G[base, 2] - G[base, 1])
    }
  }
  out
}

# ---------------------------------------------------------------------------
# user-facing
# ---------------------------------------------------------------------------

#' Descriptive changes and DiD for T >= 2 elections on the same units
#'
#' \code{eiChange()} is the descriptive layer of difference-in-differences with ecological
#' data (theory: \code{docs/extension/PHASE2e_did_general.md}, cited as 2e, Secs 2-4 and 7).
#' Under a fixed population (the same people, with the same group, at every election; 2e (FP))
#' the T-election problem is a count-layer problem whose units are the unit-elections
#' (Lemma E1), so every result of \code{\link{eiCounts}} applies, including rows that couple
#' elections. Everything here is \strong{descriptive}: a statement about the realised tables
#' of the actual groups. Causal statements need \code{\link{eiATT}}.
#'
#' \strong{Targets.} By default, the event-study DiDs \eqn{D_t = g_t - g_{base}} of group
#' \code{row} against \code{versus} on outcome \code{col}, for every election \eqn{t} other
#' than \code{base} (2e Sec. 1.3). Other targets: \code{\link{did}} (any pair of elections),
#' \code{\link{gap_at}}, \code{\link{rate_at}}, \code{\link{bin_did}}.
#'
#' \strong{Margins only (Delta-0, Prop. E3).} The range of each DiD is the Minkowski sum of
#' the elections' gap ranges (exact LP). The \emph{joint} set of the event-study DiDs is not
#' the box of their separate ranges: they share the base election's table (E3(c)). The print
#' says whether the corner of the box is attainable; \code{\link{eiJoint}} tests any vector.
#'
#' \strong{Count rows} (\code{restrictions}, a named list as in \code{eiCounts}; every row is
#' \code{object = "actual"}: a belief or fact about the realised tables).
#' \itemize{
#'   \item \code{sign_change = list(type = "dW" or "dB", row, col, from, to, sign, version,
#'     other, where, slack)}: realised association signs on changes (Lemma E4, Prop. E5).
#'     \code{"W1"} (one-vs-rest, default) is a half-line in \eqn{Z = Y_{rct}(S) - Y_{rcu}(S)}:
#'     for the one-vs-rest DiD on the same (row, col, from, to, units) the set has the closed
#'     form of E5(b); \code{"pair"} (\eqn{\delta_B}) and \code{"W2"} (\eqn{\delta_W}) are LP rows.
#'   \item \code{cr_change = list(row, col, from, to, version, other, where, slack)}: realised
#'     contextual reinforcement on changes, a union of two row systems (E5(d), (f)); the result
#'     is the union of the branch sets.
#'   \item \code{order_change = list(row, over = "rest" or a row, col, from, to, units, sign)}:
#'     the ordering of changes in every unit, \eqn{\Delta\beta_{irc} \ge \Delta\beta_{isc}}.
#'   \item \code{contrast_change = list(col, from, to, low, high, beliefs = c("mix", "rate"))}:
#'     the count contrast on changes of Prop. E6 (R = 2), for prespecified disjoint bins.
#'   \item within-election rows \code{known, zero, box, total, linear, order_unit, sign, cr,
#'     odds}, as in \code{\link{eiCounts}}, each with \code{election =}.
#' }
#' Several period pairs are rows on the same product: one LP; CR in K pairs is a union of at
#' most \eqn{2^K} branch LPs (E5(f)). Sign and CR rows print EXPERIMENTAL, as in eiCounts.
#'
#' \strong{Integer tables} (\code{integer}, as in \code{eiCounts}). With exactly integer
#' margins: for a one-vs-rest DiD whose only rows are W1 signs or CR on the same changes the
#' integer set is the lattice of Prop. E26 (step \eqn{\kappa_r}), computed in exact
#' arithmetic ("integer-exact (E26 lattice)"); with R = 2, change rows on one column whose
#' period pairs form paths and whose unit sets are laminar per pair (and single-cell boxes) are
#' in the totally unimodular class of Theorem E27 ("integer-exact (TU class)", class "E27");
#' rows of 6a's classes on the unit-elections are handled as in eiCounts; anything else is a
#' certified MILP per endpoint (GLPK, \code{time_limit}). Non-integer margins switch the
#' check off ("exact only for the stated margins"; gap-register item G7).
#'
#' \strong{Model layer} (\code{model}): \code{\link{fd}()} (first-difference LCE, 2e Sec. 3:
#' Lemma E7, Props E8, E9, E23) or \code{\link{lce_periods}()} (per-period pooled LCE with
#' model rows on changes, Prop. E11). Each change target gets an actual-target model interval
#' (projection interval + the two-election residual allowance of Lemma E12), intersected with
#' the count set (Theorem E13, 6a R27). Count rows and model rows never share one LP.
#'
#' \strong{Inference.} Count-layer sets are exact given the margins and rows: no sampling
#' inference (2e Sec. 6). Unions (CR branches) are reported whole, never the branch that holds
#' at an estimate (Prop. E24).
#'
#' \strong{Gap register} (2e Sec. 8; deferred, not implemented). Different group counts across
#' elections (entry, exit, aging-in, interval bases, time-varying membership: G1-G3, G5, G9),
#' outcome totals that differ from the group totals (a voters-only base), or different unit
#' sets across elections (G6) are refused with a message. \code{margins = "estimated"} (G7)
#' switches the integer check off and warns that results are exact only for the stated margins.
#'
#' @param data a list of data frames, one per election in time order (names = election
#'   labels), or one data frame with an election variable named in \code{election}.
#' @param formula \code{cbind(outcome counts) ~ cbind(group counts)}; or give \code{rows} and
#'   \code{cols}.
#' @param rows,cols names of the group-count and outcome-count variables.
#' @param unit name of the unit-label variable (units are matched across elections by label).
#' @param election name of the election variable when \code{data} is one data frame.
#' @param elections optional election labels in time order.
#' @param base the base election (the last before treatment); required when there are three or
#'   more elections, defaults to the first of two.
#' @param row,versus,col group, comparison (\code{"rest"} or a group) and outcome of the
#'   default event-study targets.
#' @param targets optional list of targets (see \code{\link{change_targets}}).
#' @param restrictions named list of count rows (see Details).
#' @param model optional model layer: \code{\link{fd}()} or \code{\link{lce_periods}()}.
#' @param integer \code{"auto"}, \code{"lp"} or \code{"exact"} (as in \code{eiCounts}).
#' @param margins \code{"exact"} (default) or \code{"estimated"} (G7: integer check off, warning).
#' @param population \code{"fixed"} (2e (FP), default) or \code{"stayers"} (the margins supplied are
#'   those of the stayers with fixed group membership; 2e Def. 5.1 item 0).
#' @param binding compute binding rows (drop-one test).
#' @param time_limit seconds per integer program.
#' @param alpha level of the model-layer intervals.
#' @return an object of class \code{"eiChange"}.
#' @examples
#' d0 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(30, 45), no = c(70, 55))
#' d1 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(32, 48), no = c(68, 52))
#' d2 <- data.frame(p = c("A", "B"), vet = c(20, 60), oth = c(80, 40), yes = c(30, 60), no = c(70, 40))
#' fit <- eiChange(list(early = d0, base = d1, post = d2), cbind(yes, no) ~ cbind(vet, oth),
#'                 unit = "p", base = "base", row = "vet", col = "yes")
#' fit
#' @export
eiChange <- function(data, formula = NULL, rows = NULL, cols = NULL, unit = NULL, election = NULL, elections = NULL,
                     base = NULL, row = 1, versus = "rest", col = 1, targets = NULL, restrictions = NULL,
                     model = NULL, integer = c("auto", "lp", "exact"), margins = c("exact", "estimated"),
                     population = c("fixed", "stayers"), binding = FALSE, time_limit = 10, alpha = 0.05) {
  integer <- match.arg(integer); margins <- match.arg(margins); population <- match.arg(population)
  inp <- .ed_input(data, rows, cols, unit, election, elections, formula)
  P <- .ed_problem(inp, margins, population)
  P$time_limit <- time_limit
  notes <- character(0)
  if (margins == "estimated") {
    warning("margins = \"estimated\": uncertain margins and compositions are gap-register item G7, deferred and not ",
            "implemented (2e Sec. 8). Results are exact only for the stated margins; no coverage for the true margins ",
            "is claimed. The integer check is off", call. = FALSE)
    notes <- c(notes, "G7: margins stated as estimated; results exact only for the stated margins (integer check off)")
  }
  if (integer == "exact" && !P$integer)
    stop("integer = \"exact\" needs exactly integer margins; these are not (or margins = \"estimated\"): use \"auto\" or \"lp\"")
  if (is.null(base)) {
    if (P$T > 2) stop("with three or more elections give 'base' (the last election before treatment; 2e Sec. 1.1)")
    base <- 1L
  }
  b0 <- .ed_elec(P, base, "base")
  tl <- .ed_targets(P, targets, row, versus, col, b0)
  specs <- .ed_parse_restrictions(P, restrictions, parent.frame())
  f <- .ed_fit(P, tl, specs, integer, binding)
  out <- structure(c(list(call = match.call(), P = P, integer = integer, margins_integer = P$integer, base = b0,
                          population = population, margins = margins, notes = notes,
                          defaults = list(row = row, versus = versus, col = col)), f), class = "eiChange")
  out$event_study <- .ed_event_study(P, tl, f$built, b0)
  if (!is.null(model)) out$model <- .ed_model_block(out, model, alpha)
  out
}

# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------
#' @export
print.eiChange <- function(x, digits = 3, ...) {
  P <- x$P; tb <- x$table
  cat("DESCRIPTIVE: changes and DiD on the realised tables (2e Secs 2-4; count layer, exact given the margins and rows)\n")
  cat(sprintf("Elections: %s (base %s); units: %d%s; population: %s; margins: %s; integer mode: %s\n",
              paste(P$elections, collapse = ", "), P$elections[x$base], P$p0,
              if (P$n_dropped0) sprintf(" (%d with N_i = 0 dropped)", P$n_dropped0) else "",
              if (x$population == "fixed") "fixed (2e (FP))" else "stayers with fixed group membership (supplied margins)",
              if (P$integer) "integer counts" else if (x$margins == "estimated") "stated as estimated (G7)" else "not integers",
              x$integer))
  out <- data.frame(target = tb$target, margins_only = .ec_fmt(tb$dd_lower, tb$dd_upper, digits), stringsAsFactors = FALSE)
  has_rows <- !is.null(x$restrictions)
  if (has_rows) {
    out$with_rows <- .ec_fmt(tb$lower, tb$upper, digits)
    for (j in seq_len(nrow(tb))) {
      pc <- x$solutions[[j]]$pieces
      if (!is.null(pc) && nrow(pc) > 1) out$with_rows[j] <- paste(.ec_fmt(pc[, 1], pc[, 2], digits), collapse = " U ")
    }
  }
  out$integer_status <- tb$integer_status
  cat("\nTarget  |  margins only (Delta-0, Prop. E3; exact LP)", if (has_rows) "  |  with the count rows", "\n", sep = "")
  print(out, row.names = FALSE, right = FALSE)
  for (j in seq_len(nrow(tb))) {
    co <- x$closed[[j]]$ovr
    if (!is.null(co))
      cat(sprintf("  %s: one-vs-rest on changes: ER_D = %.4g, NM_D = %.4g (gamma = %.3g); closed form %s; lattice step %.4g\n",
                  tb$target[j], co$ER, co$NM, co$gamma, .ec_fmt(tb$closed_lower[j], tb$closed_upper[j], digits), co$step))
    cc <- x$closed[[j]]$contrast
    if (!is.null(cc))
      cat(sprintf("  %s: contrast on changes (E6): lambda = %.4g, Dtau_H - Dtau_L = %.4g, cap = %.4g; set %s%s\n",
                  tb$target[j], cc$lambda, cc$dtau, cc$cap, if (cc$empty) "empty" else .ec_fmt(cc$lower, cc$upper, digits),
                  if (length(cc$why)) paste0(" (", paste(cc$why, collapse = "; "), ")") else ""))
  }
  es <- x$event_study
  if (!is.null(es) && !is.null(es$corner)) {
    cat(sprintf("\nJoint (event-study) set (Prop. E3(c)): not the box of the separate ranges. The corner (%s at its maximum, %s at its minimum) is %s under the margins",
                es$corner$targets[1], es$corner$targets[2], if (es$corner$margins_only) "attainable" else "NOT attainable"))
    cat(if (!es$corner$margins_only) sprintf(" (the base election's gap range has width %.3g; both DiDs use it). Use eiJoint() to test a vector.\n", es$corner$base_width) else ".\n")
  }
  if (has_rows) {
    cat("\nCount rows (object = \"actual\": beliefs or facts about the realised tables):\n")
    rt <- x$restrictions
    for (k in seq_len(nrow(rt)))
      cat(sprintf("  %-18s %s [%s row]%s\n", rt$id[k], rt$restriction[k], rt$role[k],
                  if (rt$experimental[k]) " EXPERIMENTAL (2c M2; review batch B comment 1 open)" else ""))
    for (j in seq_len(nrow(tb))) {
      br <- x$solutions[[j]]$branches
      if (is.null(br)) next
      cat(sprintf("  branches for %s (union reported; never a selected branch, Prop. E24): %s\n", tb$target[j],
                  paste(sprintf("%s %s", br$branch, ifelse(is.na(br$lower), "empty", .ec_fmt(br$lower, br$upper, digits))), collapse = "; ")))
    }
    if (any(!is.na(tb$binding_lower)))
      for (j in seq_len(nrow(tb)))
        cat(sprintf("  binding rows for %s: lower %s; upper %s\n", tb$target[j],
                    if (nzchar(tb$binding_lower[j])) tb$binding_lower[j] else "none",
                    if (nzchar(tb$binding_upper[j])) tb$binding_upper[j] else "none"))
  }
  cat("\n")
  if (any(grepl("E26", tb$integer_status)))
    cat("  integer-exact (E26 lattice): the integer values of the one-vs-rest DiD are D(Z) for integer Z (Prop. E26), exact arithmetic\n")
  if (any(tb$tu_class %in% "E27"))
    cat("  integer-exact (TU class E27): change rows on one column, paths of period pairs, laminar unit sets (Theorem E27)\n")
  if (any(grepl("MILP", tb$integer_status)))
    cat("  integer-exact (MILP): rows outside the integer classes; each endpoint a certified integer program\n")
  if (any(grepl("unresolved", tb$integer_status)))
    cat("  integer status unresolved: the range shown is the LP over the rounded rows, an outer bound (never integer-exact)\n")
  if (any(grepl("integer-infeasible", tb$integer_status)))
    cat("  integer-infeasible: continuous tables satisfy the rows, but no whole-number table does (e.g. CR on changes, Prop. E26(c))\n")
  if (!P$integer) cat("  Note:", .ec_txt("nonint"), "(gap-register item G7)\n")
  if (length(P$log$rat))
    cat(sprintf("  Note: values given as doubles were rationalised: %s\n", paste(utils::head(unique(P$log$rat), 5), collapse = "; ")))
  for (nt in unique(c(P$log$notes, x$notes))) cat("  Note:", nt, "\n")
  cat("  No sampling inference: count-layer sets are exact given the margins and rows (2e Sec. 6).\n")
  if (!is.null(x$model)) .ed_print_model(x$model, digits)
  cat("  Causal statements: eiATT(fit, ...) with a stated population, comparison group, trend assumption,\n",
      "  spillover allowance and no anticipation (2e Def. 5.1).\n", sep = "")
  invisible(x)
}

#' @export
witness.eiChange <- function(object, target = NULL, ...) {
  w <- witness.eiCounts(object, target = target)
  P <- object$P
  for (k in seq_along(w)) {
    arr <- w[[k]]$table                                  # stacked units x rows x cols
    a4 <- array(NA_real_, c(P$p0, P$R, P$C, P$T), dimnames = list(P$units0, P$rows, P$cols, P$elections))
    for (t in seq_len(P$T)) a4[, , , t] <- arr[.ed_su(P, t, seq_len(P$p0)), , , drop = FALSE]
    w[[k]]$table4 <- a4
  }
  w
}
