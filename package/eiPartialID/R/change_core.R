# Changes and DiD, descriptive layer (internal core): the stacked T-election count problem.
# Theory: docs/extension/PHASE2e_did_general.md (cited as 2e): Lemma E1 (stacking), Props E3,
# E5, E6, E26 and Theorem E27; built on the count layer of PHASE6a (counts_core.R, counts_api.R).
#
# Lemma E1: under the fixed population (FP) the T-election problem is the count-layer problem
# whose units are the unit-elections (i, t), with row margins A[i, ] (the same at every t) and
# column margins B_t[i, ]. Stacked unit index k = (t - 1) p0 + i, so cell (i, t, r, c) is the
# count-layer variable .ec_vix(P, k, r, c). Every count-layer routine (margins, exact rounding,
# LP, certified MILP, witnesses, unions) then applies verbatim, including rows that couple
# elections ("change rows").
#
# Row kinds used here, in addition to those of counts_core.R:
#   "change"  coefficient +1 on cell (r, c) at election t and -1 at election u, over a unit set
#             U (2e Theorem E27: sum_{i in U} (n_{i r c t} - n_{i r c u}) in [lo, hi]); row$pair =
#             c(u, t), row$units = U (indices of the p0 units), row$cell = c(r, c).

# stacked unit indices of units S (indices 1..p0) at election t
.ed_su <- function(P, t, S) (t - 1L) * P$p0 + S
# stacked variable index of cell (i, t, r, c)
.ed_vix <- function(P, i, t, r, c) .ec_vix(P, .ed_su(P, t, i), r, c)

# ---------------------------------------------------------------------------
# input: a list of per-election data frames, or a long data frame
# ---------------------------------------------------------------------------
.ed_input <- function(data, rows, cols, unit, election, elections, formula) {
  if (!is.null(formula)) { cols <- all.vars(formula[[2]]); rows <- all.vars(formula[[3]]) }
  if (is.null(rows) || is.null(cols)) stop("give the group (row) and outcome (column) count variables: 'formula' or 'rows' and 'cols'")
  if (is.data.frame(data)) {
    if (is.null(election)) stop("with one data frame, name its election variable in 'election' (or pass a list of data frames, one per election)")
    ev <- data[[election]]
    if (is.null(ev)) stop("election variable '", election, "' not found in data")
    lev <- if (!is.null(elections)) as.character(elections) else if (is.factor(ev)) levels(droplevels(ev)) else
      if (is.numeric(ev)) as.character(sort(unique(ev))) else unique(as.character(ev))
    if (any(!lev %in% as.character(ev))) stop("elections not found in the election variable: ", paste(setdiff(lev, as.character(ev)), collapse = ", "))
    dl <- lapply(lev, function(l) data[as.character(ev) == l, , drop = FALSE]); names(dl) <- lev
  } else if (is.list(data) && all(vapply(data, is.data.frame, TRUE))) {
    dl <- data
    if (is.null(names(dl)) || any(!nzchar(names(dl)))) names(dl) <- paste0("e", seq_along(dl))
    if (!is.null(elections)) {
      if (!all(as.character(elections) %in% names(dl))) stop("'elections' must name elements of the data list")
      dl <- dl[as.character(elections)]
    }
  } else stop("'data' must be a list of data frames (one per election, in time order) or one data frame with an election variable")
  if (length(dl) < 2) stop("eiChange needs at least two elections (2e Sec. 1.1: T >= 2)")
  if (anyDuplicated(names(dl))) stop("election labels must be distinct")
  lab_of <- function(d) {
    if (is.null(unit)) return(if (!is.null(rownames(d)) && !identical(rownames(d), as.character(seq_len(nrow(d))))) rownames(d) else NULL)
    if (is.character(unit) && length(unit) == 1 && unit %in% names(d)) as.character(d[[unit]]) else
      stop("'unit' must name a unit-label variable present in every election's data")
  }
  labs <- lapply(dl, lab_of)
  for (e in seq_along(dl)) {
    miss <- setdiff(c(rows, cols), names(dl[[e]]))
    if (length(miss)) stop("variables not found in election ", names(dl)[e], ": ", paste(miss, collapse = ", "))
  }
  if (any(vapply(labs, is.null, TRUE))) {
    if (length(unique(vapply(dl, nrow, 0L))) != 1)
      stop("the elections have different numbers of units and no 'unit' labels are given: changing unit sets ",
           "(boundaries) are gap-register item G6, not implemented (2e Sec. 8; Lemma G6: coarsen to common units)")
    labs <- rep(list(as.character(seq_len(nrow(dl[[1]])))), length(dl))
  }
  for (e in seq_along(dl)) if (anyDuplicated(labs[[e]])) stop("duplicated unit labels in election ", names(dl)[e])
  u1 <- labs[[1]]
  for (e in seq_along(dl)[-1]) if (!setequal(labs[[e]], u1))
    stop("the unit sets differ across elections (", names(dl)[1], " vs ", names(dl)[e], "): changing unit ",
         "boundaries are gap-register item G6, deferred and not implemented (2e Sec. 8). Lemma G6: aggregate to ",
         "units that are stable across elections (valid, weakly wider) and call eiChange on those")
  # align every election to the first election's unit order
  dl <- lapply(seq_along(dl), function(e) dl[[e]][match(u1, labs[[e]]), , drop = FALSE]); names(dl) <- names(labs)
  A <- lapply(dl, function(d) { m <- as.matrix(d[, rows, drop = FALSE]); storage.mode(m) <- "double"; m })
  B <- lapply(dl, function(d) { m <- as.matrix(d[, cols, drop = FALSE]); storage.mode(m) <- "double"; m })
  list(data = dl, A = A, B = B, rows = rows, cols = cols, units = u1, elections = names(dl))
}

# fixed population (FP, 2e Sec. 1.2) and margin checks; refusals for the gap register (2e Sec. 8)
.ed_check_fp <- function(inp, margins) {
  A <- inp$A; B <- inp$B; el <- inp$elections
  for (e in seq_along(A)) {
    if (any(!is.finite(A[[e]])) || any(!is.finite(B[[e]]))) stop("missing or non-finite counts in election ", el[e])
    if (any(A[[e]] < 0) || any(B[[e]] < 0)) stop("negative counts in election ", el[e])
  }
  for (e in seq_along(A)[-1]) {
    d <- abs(A[[e]] - A[[1]])
    if (any(d > 1e-9 * pmax(1, abs(A[[1]])))) {
      stop("the group counts differ across elections (", el[1], " vs ", el[e], "; first at unit ",
           inp$units[which(rowSums(d > 1e-9 * pmax(1, abs(A[[1]]))) > 0)[1]], "). eiChange needs a fixed population (2e (FP)): ",
           "the same people, with the same group, at every election. A changing population is gap-register item G1 (entry), ",
           "G2 (exit), G3 (aging-in), G5 (interval bases) or G9 (time-varying group membership); these are deferred and ",
           "not implemented (2e Sec. 8). Supply the stayers' margins (population = \"stayers\") if they are known", call. = FALSE)
    }
  }
  N <- rowSums(A[[1]])
  for (e in seq_along(B)) {
    NB <- rowSums(B[[e]])
    bad <- abs(NB - N) > 1e-9 * pmax(1, N)
    if (any(bad))
      stop("in election ", el[e], " the outcome counts do not add up to the group counts (first at unit ",
           inp$units[which(bad)[1]], "). The outcome must be on the common base of everyone in the fixed population ",
           "(e.g. include abstention: 2e Sec. 1.1). A voters-only base or a population that changes between elections ",
           "is gap-register item G1-G5, deferred and not implemented (2e Sec. 8)", call. = FALSE)
  }
  int <- all(vapply(c(A, B), function(m) all(m == round(m)), TRUE)) && sum(N) < 2^53
  list(N = N, integer = int && margins == "exact")
}

# ---------------------------------------------------------------------------
# the stacked problem
# ---------------------------------------------------------------------------
.ed_problem <- function(inp, margins, population) {
  fp <- .ed_check_fp(inp, margins)
  keep <- fp$N > 0
  if (!any(keep)) stop("no units with N_i > 0")
  T <- length(inp$A); p0 <- sum(keep)
  As <- do.call(rbind, lapply(inp$A, function(m) m[keep, , drop = FALSE]))
  Bs <- do.call(rbind, lapply(inp$B, function(m) m[keep, , drop = FALSE]))
  lab <- as.vector(outer(inp$units[keep], inp$elections, function(u, e) paste0(u, "@", e)))
  P <- .ec_problem(As, Bs, inp$rows, inp$cols, lab, integer_input = fp$integer, data = NULL)
  if (P$p != T * p0) stop("internal: stacked problem has the wrong number of units")
  P$T <- T; P$p0 <- p0; P$elections <- inp$elections; P$units0 <- inp$units[keep]
  P$n_dropped0 <- sum(!keep); P$classify <- .ed_classify; P$margins <- margins; P$population <- population
  P$A0 <- inp$A[[1]][keep, , drop = FALSE]; P$N0 <- fp$N[keep]
  P$Bt <- lapply(inp$B, function(m) m[keep, , drop = FALSE])
  # one-election problems (same margins), sharing the stacked problem's number type and log:
  # used to parse within-election rows and to compute one-election thresholds exactly
  P$Pt <- lapply(seq_len(T), function(t) {
    Q <- .ec_problem(P$A0, P$Bt[[t]], inp$rows, inp$cols, P$units0, integer_input = fp$integer,
                     data = inp$data[[t]][keep, , drop = FALSE])
    Q$integer <- P$integer; Q$log <- P$log; Q$num <- P$num
    Q
  })
  P
}

.ed_elec <- function(P, e, what) {
  if (is.null(e)) stop(what, ": an election is required")
  if (length(e) != 1) stop(what, ": one election")
  if (is.numeric(e)) { if (e < 1 || e > P$T || e != round(e)) stop(what, ": election index out of range"); return(as.integer(e)) }
  k <- match(as.character(e), P$elections)
  if (is.na(k)) stop(what, ": unknown election '", e, "' (elections: ", paste(P$elections, collapse = ", "), ")")
  k
}

# ---------------------------------------------------------------------------
# targets: did(), rate_at(), gap_at(), bin_did()
# ---------------------------------------------------------------------------

#' Targets for eiChange()
#'
#' Constructors for the targets of \code{\link{eiChange}}: linear functionals of the realised
#' tables of several elections on the same units (2e Sec. 1.3), each over a known denominator
#' fixed by the design.
#' \itemize{
#'   \item \code{did(row, versus, col, from, to, where)}: the DiD
#'     \eqn{D_c(r,s;u,t) = g_{rs,ct} - g_{rs,cu}}, the change in the gap between group
#'     \code{row} and group \code{versus} (another row, or \code{"rest"} for everyone else)
#'     between election \code{from} (u) and election \code{to} (t).
#'   \item \code{rate_at(row, col, election, where)}: the group rate \eqn{B_{rct}}.
#'   \item \code{gap_at(row, versus, col, election, where)}: the gap \eqn{g_{rs,ct}}.
#'   \item \code{bin_did(col, from, to, low, high)}: the bin-averaged DiD
#'     \eqn{\bar D} of 2e Prop. E6 for two prespecified disjoint unit sets (R = 2 only).
#' }
#' Elections are given by label or position. \code{where}, \code{low} and \code{high} select
#' units as in \code{\link{count_targets}} (expressions are evaluated in the first election's
#' data; they should use design variables only).
#'
#' @param row,versus group labels or indices (\code{versus = "rest"}: all other groups).
#' @param col outcome label or index.
#' @param from,to,election election labels or positions.
#' @param where,low,high unit selections.
#' @param label optional label.
#' @return an object of class \code{"edTarget"}.
#' @name change_targets
NULL

.ed_tspec <- function(kind, args, label, env) structure(list(kind = kind, args = args, label = label, env = env),
                                                       class = "edTarget")

#' @rdname change_targets
#' @export
did <- function(row, versus = "rest", col, from, to, where = NULL, label = NULL)
  .ed_tspec("did", list(row = row, versus = versus, col = col, from = from, to = to, where = substitute(where)),
            label, parent.frame())

#' @rdname change_targets
#' @export
rate_at <- function(row, col, election, where = NULL, label = NULL)
  .ed_tspec("rate_at", list(row = row, col = col, election = election, where = substitute(where)), label, parent.frame())

#' @rdname change_targets
#' @export
gap_at <- function(row, versus = "rest", col, election, where = NULL, label = NULL)
  .ed_tspec("gap_at", list(row = row, versus = versus, col = col, election = election, where = substitute(where)),
            label, parent.frame())

#' @rdname change_targets
#' @export
bin_did <- function(col, from, to, low, high, label = NULL)
  .ed_tspec("bin_did", list(col = col, from = from, to = to, low = substitute(low), high = substitute(high)),
            label, parent.frame())

# unit selection in the first election's problem: indices 1..p0 with attribute wkind
.ed_where <- function(P, w, env, what) .ec_where(P$Pt[[1]], w, env, what)

.ed_term <- function(P, t, S, rws, cls, w) list(units = .ed_su(P, t, S), rows = rws, cols = cls, w = w, wd = as.double(w))

# gap g_{rs,ct}(S) as terms; s = NULL means the rest. Returns terms, plus the constant-form
# representation for the rest (kappa Y_rc - Y_c / (N - A_r)), used for the allowance.
.ed_gap_terms <- function(P, t, S, r, s, c, sign = 1) {
  nm <- P$num
  Ar <- sum(nm(P$A0[S, r]))
  if (as.double(Ar) <= 0) stop("group ", P$rows[r], " has no members in the unit set")
  if (is.null(s)) {
    N <- sum(nm(P$N0[S])); rest <- setdiff(seq_len(P$R), r)
    if (as.double(Ar) >= as.double(N)) stop("'rest' is empty in the unit set (A_r = N)")
    list(terms = list(.ed_term(P, t, S, r, c, sign / Ar), .ed_term(P, t, S, rest, c, -sign / (N - Ar))),
         terms2 = list(.ed_term(P, t, S, r, c, sign * N / (Ar * (N - Ar)))),
         const2 = -sign * sum(nm(P$Bt[[t]][S, c])) / (N - Ar))
  } else {
    As <- sum(nm(P$A0[S, s]))
    if (as.double(As) <= 0) stop("group ", P$rows[s], " has no members in the unit set")
    tm <- list(.ed_term(P, t, S, r, c, sign / Ar), .ed_term(P, t, S, s, c, -sign / As))
    list(terms = tm, terms2 = tm, const2 = nm(0))
  }
}

.ed_rowidx <- function(P, v, what) .ec_idx(v, P$rows, what)
.ed_colidx <- function(P, v, what) .ec_idx(v, P$cols, what)
.ed_versus <- function(P, versus, r) {
  if (identical(versus, "rest")) return(if (P$R == 2) setdiff(1:2, r) else NULL)
  s <- .ed_rowidx(P, versus, "versus")
  if (length(s) != 1 || s == r) stop("'versus' must be one other group, or \"rest\"")
  s
}

.ed_resolve <- function(P, ts, name = NULL) {
  nm <- P$num; a <- ts$args
  us <- function(S) if (length(S) == P$p0) "" else if (length(S) <= 3) sprintf(" [%s]", paste(P$units0[S], collapse = ", ")) else
    sprintf(" [%d units]", length(S))
  rl <- function(r) paste(P$rows[r], collapse = "+")
  out <- list(kind = ts$kind)
  if (ts$kind %in% c("did", "gap_at", "rate_at")) {
    S <- .ed_where(P, a$where, ts$env, ts$kind)
    r <- .ed_rowidx(P, a$row, "row"); c <- .ed_colidx(P, a$col, "col")
    if (length(r) != 1 || length(c) != 1) stop(ts$kind, "(): one row and one column")
    out$S <- as.vector(S); out$r <- r; out$c <- c; out$wkinds <- attr(S, "wkind")
  }
  if (ts$kind == "did") {
    u <- .ed_elec(P, a$from, "did(from)"); t <- .ed_elec(P, a$to, "did(to)")
    if (u == t) stop("did(): 'from' and 'to' must differ")
    s <- .ed_versus(P, a$versus, out$r)
    gt <- .ed_gap_terms(P, t, out$S, out$r, s, out$c, 1); gu <- .ed_gap_terms(P, u, out$S, out$r, s, out$c, -1)
    out$terms <- c(gt$terms, gu$terms); out$const <- nm(0)
    out$reps <- list(list(terms = out$terms, const = out$const),
                     list(terms = c(gt$terms2, gu$terms2), const = gt$const2 + gu$const2))
    out$s <- s; out$u <- u; out$t <- t; out$change <- TRUE
    out$label <- sprintf("DiD %s vs %s on %s, %s -> %s%s", rl(out$r), if (is.null(s)) "rest" else rl(s), P$cols[out$c],
                         P$elections[u], P$elections[t], us(out$S))
    out$den <- sprintf("A_%s = %s", rl(out$r), format(as.double(sum(nm(P$A0[out$S, out$r])))))
  } else if (ts$kind == "gap_at") {
    t <- .ed_elec(P, a$election, "gap_at(election)"); s <- .ed_versus(P, a$versus, out$r)
    g <- .ed_gap_terms(P, t, out$S, out$r, s, out$c, 1)
    out$terms <- g$terms; out$const <- nm(0)
    out$reps <- list(list(terms = g$terms, const = nm(0)), list(terms = g$terms2, const = g$const2))
    out$s <- s; out$t <- t; out$change <- FALSE
    out$label <- sprintf("gap %s vs %s on %s at %s%s", rl(out$r), if (is.null(s)) "rest" else rl(s), P$cols[out$c],
                         P$elections[t], us(out$S))
    out$den <- "A_r(S)"
  } else if (ts$kind == "rate_at") {
    t <- .ed_elec(P, a$election, "rate_at(election)")
    Ar <- sum(nm(P$A0[out$S, out$r]))
    if (as.double(Ar) <= 0) stop("rate_at(): the denominator A_r(S) is zero")
    out$terms <- list(.ed_term(P, t, out$S, out$r, out$c, 1 / Ar)); out$const <- nm(0)
    out$reps <- list(list(terms = out$terms, const = out$const))
    out$t <- t; out$change <- FALSE
    out$label <- sprintf("rate %s -> %s at %s%s", rl(out$r), P$cols[out$c], P$elections[t], us(out$S))
    out$den <- sprintf("A_%s = %s", rl(out$r), format(as.double(Ar)))
  } else if (ts$kind == "bin_did") {
    bd <- .ed_bins(P, a, ts$env, "bin_did()")
    out <- c(out, bd)
    out$terms <- bd$terms; out$const <- nm(0)
    out$reps <- list(list(terms = out$terms, const = out$const))
    out$change <- TRUE
    out$label <- sprintf("bin-averaged DiD on %s, %s -> %s (low %s, high %s)", P$cols[bd$c], P$elections[bd$u],
                         P$elections[bd$t], us(bd$L), us(bd$H))
    out$den <- "A_r(bin)"
  } else stop("unknown eiChange target kind '", ts$kind, "'")
  out$design <- !("outcome" %in% out$wkinds)
  if (!out$design) out$design_reason <- "outcome_set"
  out$given_set <- "given" %in% out$wkinds
  if (!is.null(ts$label)) out$label <- ts$label else if (!is.null(name) && nzchar(name)) out$label <- name
  out$obj <- .ec_obj(P, out$terms)
  out
}

# bins of 2e Prop. E6 (R = 2): D-bar = (1/2) sum_b (Dm_{1,b} - Dm_{2,b}) as terms; lambda, Dtau
.ed_bins <- function(P, a, env, what) {
  if (P$R != 2) stop(what, ": the contrast on changes of 2e Prop. E6 is stated for R = 2 (with R >= 3 use named pairs; not covered by E6)")
  nm <- P$num
  c <- .ed_colidx(P, a$col, "col"); u <- .ed_elec(P, a$from, paste0(what, " from")); t <- .ed_elec(P, a$to, paste0(what, " to"))
  if (u == t) stop(what, ": 'from' and 'to' must differ")
  L <- .ed_where(P, a$low, env, paste0(what, " low")); H <- .ed_where(P, a$high, env, paste0(what, " high"))
  if (length(intersect(L, H))) stop(what, ": the bins 'low' and 'high' must be disjoint")
  if (!length(L) || !length(H)) stop(what, ": empty bin")
  terms <- list()
  for (b in list(L, H)) {
    A1 <- sum(nm(P$A0[b, 1])); A2 <- sum(nm(P$A0[b, 2]))
    if (as.double(A1) <= 0 || as.double(A2) <= 0) stop(what, ": each bin needs members of both groups")
    terms <- c(terms, list(.ed_term(P, t, b, 1, c, 1 / (2 * A1)), .ed_term(P, u, b, 1, c, -1 / (2 * A1)),
                           .ed_term(P, t, b, 2, c, -1 / (2 * A2)), .ed_term(P, u, b, 2, c, 1 / (2 * A2))))
  }
  xb <- function(b) sum(nm(P$A0[b, 1])) / sum(nm(P$N0[b]))
  dtau <- function(b) (sum(nm(P$Bt[[t]][b, c])) - sum(nm(P$Bt[[u]][b, c]))) / sum(nm(P$N0[b]))
  lambda <- xb(H) - xb(L)
  if (as.double(lambda) <= 0) stop(what, ": the 'high' bin must have the larger pooled share of ", P$rows[1], " (lambda > 0)")
  list(c = c, u = u, t = t, L = as.vector(L), H = as.vector(H), terms = terms, lambda = lambda,
       dtau = dtau(H) - dtau(L), cap = (dtau(H) - dtau(L)) / lambda,
       wkinds = unique(c(attr(L, "wkind"), attr(H, "wkind"))))
}

# default targets: the event-study DiDs relative to the base election (2e Sec. 1.3)
.ed_default_targets <- function(P, row, versus, col, base) {
  lapply(setdiff(seq_len(P$T), base), function(t)
    structure(list(kind = "did", args = list(row = row, versus = versus, col = col, from = base, to = t, where = NULL),
                   label = NULL, env = globalenv()), class = "edTarget"))
}

.ed_targets <- function(P, targets, row, versus, col, base) {
  if (is.null(targets)) targets <- .ed_default_targets(P, row, versus, col, base)
  if (inherits(targets, "edTarget")) targets <- list(targets)
  if (!all(vapply(targets, inherits, TRUE, "edTarget")))
    stop("eiChange targets must be built with did(), rate_at(), gap_at() or bin_did()")
  nms <- names(targets); if (is.null(nms)) nms <- rep("", length(targets))
  out <- lapply(seq_along(targets), function(j) .ed_resolve(P, targets[[j]], nms[j]))
  labs <- make.unique(vapply(out, function(t) t$label, ""), sep = " #")
  for (j in seq_along(out)) out[[j]]$label <- labs[j]
  out
}

# ---------------------------------------------------------------------------
# restrictions
# ---------------------------------------------------------------------------
.ed_change_types <- c("sign_change", "cr_change", "order_change", "contrast_change")
.ed_within_types <- c("known", "zero", "box", "total", "linear", "order_unit", "sign", "cr", "odds")
.ed_aliases <- c(assoc = "sign", order = "order_unit", known_cell = "known", zeros = "zero",
                 sign_on_changes = "sign_change", cr_on_changes = "cr_change")

.ed_parse_restrictions <- function(P, restrictions, env) {
  if (is.null(restrictions) || !length(restrictions)) return(list())
  if (is.null(names(restrictions)) || any(names(restrictions) == "")) stop("'restrictions' must be a named list")
  out <- list()
  for (j in seq_along(restrictions)) {
    type <- names(restrictions)[j]
    if (type %in% names(.ed_aliases)) type <- .ed_aliases[[type]]
    if (!type %in% c(.ed_change_types, .ed_within_types)) stop("unknown eiChange restriction type '", type, "'")
    v <- restrictions[[j]]
    if (is.null(v) || !length(v)) next
    specs <- if (is.list(v) && length(v) && all(vapply(v, is.list, TRUE))) v else list(v)
    for (s in specs) {
      k <- sum(vapply(out, function(o) o$type == type, TRUE)) + 1
      id <- paste0(type, "[", k, "]")
      ps <- if (type %in% .ed_change_types) .ed_parse_change(P, type, s, env, id) else
        if (type == "linear" && !is.null(s$cells) && "election" %in% names(s$cells)) .ed_parse_linear_multi(P, s, env, id) else
          .ed_parse_within(P, type, s, env, id)
      ps$id <- id
      out[[length(out) + 1]] <- ps
    }
  }
  out
}

.ed_object_check <- function(s, type) {
  obj <- if (is.null(s$object)) "actual" else as.character(s$object)
  if (length(obj) != 1 || !obj %in% c("actual", "model"))
    stop("restriction '", type, "': object must be \"actual\" (count layer) or \"model\"; one belief cannot be imposed in both ",
         "layers under one label (realised and model beliefs are two assumptions; 2e Sec. 4.2, NS Part 3A fn. 2)")
  if (obj == "model")
    stop("restriction '", type, "' has object = \"model\": model rows on changes belong to the model layer ",
         "(eiChange(model = fd(...)) or lce_periods(...)), never to the count layer (2e Sec. 4.1)")
  if (!is.null(s$on) && !identical(s$on, "changes") && type %in% .ed_change_types)
    stop("restriction '", type, "': 'on' must be \"changes\"")
  invisible(TRUE)
}

# within-election rows: parsed and built on that election's one-election problem, then mapped to
# the stacked unit-elections (Lemma E1)
.ed_parse_within <- function(P, type, s, env, id) {
  .ed_object_check(s, type)
  t <- .ed_elec(P, s$election, paste0(type, " (election)"))
  Q <- P$Pt[[t]]
  s2 <- s; s2$election <- NULL; s2$on <- NULL
  ps <- .ec_parse_spec(Q, type, s2, env)
  ps$id <- id
  ps2 <- ps; if (type == "sign") ps2$type <- ps$type_a
  lab <- paste0(.ec_label_of(Q, ps2, type), " at ", P$elections[t])
  rs <- .ec_spec_rows(Q, ps, id)
  off <- (t - 1L) * P$p0 * P$R * P$C
  mp <- function(row) {
    if (is.null(row)) return(NULL)
    row$idx <- row$idx + off
    if (!is.null(row$units)) row$units <- .ed_su(P, t, row$units)
    if (!P$integer) row$coefq <- NULL
    row
  }
  rows <- lapply(rs$rows, mp)
  un <- if (is.null(rs$union)) NULL else list(plus = lapply(rs$union$plus, mp), minus = lapply(rs$union$minus, mp))
  list(type = type, scope = "within", election = t, role = ps$role, label = lab, rows = rows, union = un,
       experimental = type %in% c("sign", "cr"), tuple = NULL)
}

# a general linear row on cells of several elections (a row on the product, Lemma E1 / E5(e)):
# cells = data.frame(unit, row, col, election, coef); lower / upper bounds
.ed_parse_linear_multi <- function(P, s, env, id) {
  .ed_object_check(s, "linear")
  cl <- s$cells
  if (!all(c("unit", "row", "col", "election", "coef") %in% names(cl)))
    stop("linear: 'cells' must have columns unit, row, col, election, coef")
  ui <- .ec_units(P$Pt[[1]], units = cl$unit, env = env, what = "linear", ordered = FALSE)
  um <- if (is.character(cl$unit) || is.factor(cl$unit)) match(as.character(cl$unit), P$units0) else match(as.integer(cl$unit), P$Pt[[1]]$orig)
  if (any(is.na(um))) stop("linear: unknown units")
  el <- vapply(seq_len(nrow(cl)), function(j) .ed_elec(P, cl$election[j], "linear (election)"), 0L)
  idx <- .ed_vix(P, um, el, .ed_rowidx(P, cl$row, "row"), .ed_colidx(P, cl$col, "col"))
  u <- sort(unique(idx))
  cfq <- .ec_rat(if (is.character(cl$coef)) cl$coef else as.numeric(cl$coef), P$log)
  cf <- gmp::as.bigq(rep(0, length(u)))
  for (j in seq_along(u)) cf[j] <- sum(cfq[idx == u[j]])
  if (!P$integer) cf <- as.double(cf)
  lo <- if (is.null(s$lower) || !is.finite(as.numeric(s$lower))) NULL else s$lower
  hi <- if (is.null(s$upper) || !is.finite(as.numeric(s$upper))) NULL else s$upper
  if (!is.null(s$value)) lo <- hi <- s$value
  role <- if (is.null(s$role)) "belief" else match.arg(s$role, c("identity", "belief"))
  list(type = "linear", scope = "change", role = role,
       label = sprintf("linear row on %d cells across elections in [%s, %s]", length(u), format(if (is.null(lo)) -Inf else lo), format(if (is.null(hi)) Inf else hi)),
       rows = list(.ec_row(P, id, "gen", u, cf, lo, hi)), union = NULL, experimental = FALSE, tuple = NULL)
}

.ed_sign_arg <- function(sg, what) {
  if (is.null(sg)) sg <- 1
  if (is.character(sg)) sg <- switch(gsub(" ", "", sg), ">=0" = 1, ">0" = 1, "+" = 1, "<=0" = -1, "<0" = -1, "-" = -1,
                                     stop(what, ": sign must be +1 or -1 (or \">=0\" / \"<=0\")"))
  sg <- sign(as.numeric(sg))
  if (length(sg) != 1 || is.na(sg) || sg == 0) stop(what, ": sign must be +1 or -1")
  sg
}

# exact one-vs-rest thresholds on changes (2e Lemma E4(b), Prop. E5): differences of the
# one-election thresholds of counts_core.R (each affine in that election's outcome margins)
.ed_ovr_change <- function(P, S, r, c, u, t) {
  th_t <- .ec_ovr_thresholds(P$Pt[[t]], S, r, c); th_u <- .ec_ovr_thresholds(P$Pt[[u]], S, r, c)
  if (is.null(th_t$ER) || is.null(th_u$ER)) return(NULL)        # Var_p(X_r) = 0 or A_r in {0, N}
  nm <- P$num
  Yt <- sum(nm(P$Bt[[t]][S, c])); Yu <- sum(nm(P$Bt[[u]][S, c]))
  out <- list(N = th_t$N, Ar = th_t$Ar, kappa = th_t$kappa, varX = th_t$varX, mu = th_t$mu,
              scale_B = th_t$scale_B, scale_W = th_t$scale_W,
              ER = th_t$ER - th_u$ER, dYc = Yt - Yu,
              ZNM = th_t$YNM - th_u$YNM,             # delta_B^Delta >= 0  <=>  Z >= ZNM
              ZER = th_t$YER_hi - th_u$YER_hi)       # delta_W^Delta >= 0  <=>  Z <= ZER
  out$gamma <- out$varX / (out$mu * (1 - out$mu))
  out$NM <- out$ER * out$gamma
  out$Dof <- function(Z) out$kappa * Z - out$dYc / (out$N - out$Ar)   # D^r_c as a function of Z (E4(b))
  # Duncan-Davis extremes of each election's aggregate Y_rce(S) (6a Lemma 2.1(a))
  dd <- function(e) {
    a <- P$A0[S, r]; b <- P$Bt[[e]][S, c]; Ni <- P$N0[S]
    c(sum(nm(pmax(0, a + b - Ni))), sum(nm(pmin(a, b))))
  }
  ddt <- dd(t); ddu <- dd(u)
  out$Zlo <- ddt[1] - ddu[2]; out$Zhi <- ddt[2] - ddu[1]
  out
}

.ed_parse_change <- function(P, type, s, env, id) {
  .ed_object_check(s, type)
  nm <- P$num
  sl <- if (is.null(s$slack)) 0 else as.numeric(s$slack)
  if (length(sl) != 1 || !is.finite(sl) || sl < 0) stop(type, ": 'slack' must be a nonnegative number")
  role <- if (is.null(s$role)) "belief" else match.arg(s$role, c("identity", "belief"))
  rl <- function(r) paste(P$rows[r], collapse = "+")
  if (type == "contrast_change") {
    bd <- .ed_bins(P, list(col = s$col, from = s$from, to = s$to, low = s$low, high = s$high), env, "contrast_change")
    beliefs <- if (is.null(s$beliefs)) c("mix", "rate") else match.arg(s$beliefs, c("mix", "rate"), several.ok = TRUE)
    obj <- .ec_obj(P, bd$terms); idx <- which(obj != 0)
    cq <- if (P$integer) .ec_objq(P, list(terms = bd$terms), idx) else obj[idx]
    rows <- list()
    if ("mix" %in% beliefs) rows[[length(rows) + 1]] <- .ec_row(P, id, "gen", idx, cq, 0, NULL)       # lambda D-bar >= 0
    if ("rate" %in% beliefs) rows[[length(rows) + 1]] <- .ec_row(P, id, "gen", idx, cq, NULL, bd$cap) # rho >= 0
    lab <- sprintf("count contrast on changes %s on %s, %s -> %s: %s (2e Prop. E6)", paste(P$rows, collapse = " vs "),
                   P$cols[bd$c], P$elections[bd$u], P$elections[bd$t], paste(beliefs, ">= 0", collapse = " and "))
    return(list(type = type, scope = "change", role = role, label = lab, rows = rows, union = NULL, experimental = FALSE,
                bins = bd, beliefs = beliefs))
  }
  r <- .ed_rowidx(P, s$row, "row"); c <- .ed_colidx(P, if (is.null(s$col)) 1 else s$col, "col")
  if (length(r) != 1 || length(c) != 1) stop(type, ": 'row' and 'col' must name one cell")
  u <- .ed_elec(P, s$from, paste0(type, " (from)")); t <- .ed_elec(P, s$to, paste0(type, " (to)"))
  if (u == t) stop(type, ": 'from' and 'to' must differ")
  S <- as.vector(.ec_units(P$Pt[[1]], if (!is.null(s$units)) s$units else s$unit, s$where, env, type))
  pr <- sprintf("%s -> %s", P$elections[u], P$elections[t])
  us <- if (length(S) == P$p0) "all units" else sprintf("%d units", length(S))
  chg <- function(Srows, lo, hi) {
    row <- .ec_row(P, id, "change", c(.ed_vix(P, Srows, t, r, c), .ed_vix(P, Srows, u, r, c)),
                   c(rep(1, length(Srows)), rep(-1, length(Srows))), lo, hi, cell = c(r, c), units = Srows)
    row$pair <- c(u, t)
    row
  }
  gen_change <- function(g_t, g_u) {
    off <- function(e) (e - 1L) * P$p0 * P$R * P$C
    idx <- c(g_t$idx + off(t), g_u$idx + off(u))
    cf <- if (P$integer) c(gmp::as.bigq(g_t$coef), -gmp::as.bigq(g_u$coef)) else c(as.double(g_t$coef), -as.double(g_u$coef))
    list(idx = idx, coef = cf)
  }
  if (type == "order_change") {
    ov <- if (is.null(s$over)) "rest" else s$over
    q <- if (identical(ov, "rest")) NULL else .ed_rowidx(P, ov, "over")
    if (!is.null(q) && q == r) stop("order_change: 'over' must differ from 'row'")
    if (!is.null(q) && P$R == 2) q <- NULL
    sg <- .ed_sign_arg(s$sign, "order_change")
    rows <- list()
    for (i in S) {
      a <- nm(P$A0[i, r]); Ni <- nm(P$N0[i])
      if (is.null(q)) {
        if (P$A0[i, r] == 0 || P$A0[i, r] == P$N0[i]) next
        db <- nm(P$Bt[[t]][i, c]) - nm(P$Bt[[u]][i, c])
        # Dbeta_r - Dbeta_rest = (n_t - n_u) kappa_i - Db / (N_i - a) >= -slack  (sign +), kappa_i = N_i / (a (N_i - a))
        kap <- Ni / (a * (Ni - a))
        th <- (db / (Ni - a) - sg * nm(sl)) / kap
        rows[[length(rows) + 1]] <- if (sg > 0) chg(i, th, NULL) else chg(i, NULL, (db / (Ni - a) + nm(sl)) / kap)
      } else {
        aq <- nm(P$A0[i, q])
        if (P$A0[i, r] == 0 || P$A0[i, q] == 0) next
        # aq (n_rt - n_ru) - a (n_qt - n_qu) >= -slack a aq   (sign +)
        idx <- c(.ed_vix(P, i, t, r, c), .ed_vix(P, i, u, r, c), .ed_vix(P, i, t, q, c), .ed_vix(P, i, u, q, c))
        cf <- if (P$integer) c(aq, -aq, -a, a) else as.double(c(P$A0[i, q], -P$A0[i, q], -P$A0[i, r], P$A0[i, r]))
        if (sg < 0) cf <- -cf
        rows[[length(rows) + 1]] <- .ec_row(P, id, "gen", idx, cf, -nm(sl) * a * aq, NULL)
      }
    }
    lab <- sprintf("realised ordering of changes in each unit (%s): Delta %s %s Delta %s on %s, %s%s", us, rl(r),
                   if (sg > 0) ">=" else "<=", if (is.null(q)) "rest" else rl(q), P$cols[c], pr,
                   if (sl > 0) sprintf(" (slack %g)", sl) else "")
    return(list(type = type, scope = "change", role = role, label = lab, rows = rows, union = NULL, experimental = FALSE,
                tuple = NULL))
  }
  # signs and CR on changes (2e Lemma E4, Prop. E5)
  if (length(S) < 2) stop(type, ": the unit set needs at least two units")
  version <- if (is.null(s$version)) "W1" else as.character(s$version)
  other <- NULL
  if (version == "pair") {
    other <- if (is.null(s$other)) { if (P$R == 2) setdiff(1:2, r) else stop(type, ": version \"pair\" needs 'other'") } else .ed_rowidx(P, s$other, "other")
    if (other == r) stop(type, ": 'other' must differ from 'row'")
  }
  w1 <- function(tp, sg) {
    th <- .ed_ovr_change(P, S, r, c, u, t)
    if (is.null(th)) stop(type, ": Var_p(X_r) > 0 and 0 < A_r < N are needed over the unit set (2e Prop. E5)")
    if (tp == "dB") {
      eps <- nm(sl) * th$scale_B * th$N
      if (sg > 0) chg(S, th$ZNM - eps, NULL) else chg(S, NULL, th$ZNM + eps)
    } else {
      eps <- nm(sl) * th$scale_W / th$varX / th$kappa
      if (sg > 0) chg(S, NULL, th$ZER + eps) else chg(S, th$ZER - eps, NULL)
    }
  }
  gen <- function(tp, ver, sg) {
    gt <- .ec_assoc_gen(P$Pt[[t]], S, tp, ver, r, c, other); gu <- .ec_assoc_gen(P$Pt[[u]], S, tp, ver, r, c, other)
    g <- gen_change(gt, gu)
    if (!length(g$coef) || all(as.double(g$coef) == 0)) {
      P$log$notes <- c(P$log$notes, paste0(id, ": the association on changes has all-zero coefficients; the row is vacuous and was dropped"))
      return(NULL)
    }
    .ec_row(P, id, "gen", g$idx, if (sg > 0) g$coef else -g$coef, -nm(sl) * gt$scale, NULL)
  }
  mk <- function(tp, ver, sg) if (ver == "W1") w1(tp, sg) else gen(tp, ver, sg)
  sl_lab <- if (sl > 0) sprintf(" (slack %g)", sl) else ""
  if (type == "sign_change") {
    tp <- as.character(if (is.null(s$type)) "dB" else s$type)
    if (!tp %in% c("dB", "dW")) stop("sign_change: type must be \"dB\" or \"dW\"")
    allowed <- if (tp == "dB") c("W1", "pair") else c("W1", "W2")
    if (!version %in% allowed) stop("sign_change: version must be one of ", paste(allowed, collapse = ", "))
    sg <- .ed_sign_arg(s$sign, "sign_change")
    row <- mk(tp, version, sg)
    lab <- sprintf("realised %s on changes %s %s on %s, %s (%s%s, %s)%s",
                   if (tp == "dB") "delta_B" else "delta_W",
                   if (version == "pair") paste0("(", rl(r), " vs ", rl(other), ")") else rl(r),
                   if (sg > 0) ">= 0" else "<= 0", P$cols[c], pr, version,
                   if (version == "W1" && tp == "dB" && P$R > 2) " vs rest" else "", us, sl_lab)
    tup <- if (version == "W1") list(r = r, c = c, u = u, t = t, S = S, kind = paste0("W1", tp), sign = sg, slack = sl) else NULL
    return(list(type = type, scope = "change", role = role, label = lab, rows = Filter(Negate(is.null), list(row)),
                union = NULL, experimental = TRUE, tuple = tup))
  }
  # cr_change: union of the two sign systems (2e Prop. E5(d), (f))
  if (!version %in% c("W1", "W2", "pair")) stop("cr_change: version must be W1, W2 or pair")
  vB <- if (version == "pair") "pair" else "W1"; vW <- if (version == "W1") "W1" else "W2"
  un <- list(plus = Filter(Negate(is.null), list(mk("dB", vB, 1), mk("dW", vW, 1))),
             minus = Filter(Negate(is.null), list(mk("dB", vB, -1), mk("dW", vW, -1))))
  lab <- sprintf("realised CR on changes %s on %s, %s (%s%s, %s)%s", rl(r), P$cols[c], pr, version,
                 if (version == "pair") paste0(" vs ", rl(other)) else "", us, sl_lab)
  tup <- if (version == "W1") list(r = r, c = c, u = u, t = t, S = S, kind = "CR", sign = 0, slack = sl) else NULL
  list(type = type, scope = "change", role = role, label = lab, rows = list(), union = un, experimental = TRUE, tuple = tup)
}

.ed_build_rows <- function(P, specs) {
  common <- list(); unions <- list()
  for (s in specs) {
    common <- c(common, s$rows)
    if (!is.null(s$union)) unions[[length(unions) + 1]] <- c(s$union, list(id = s$id))
  }
  if (P$integer) {
    common <- lapply(common, .ec_round_row)
    unions <- lapply(unions, function(u) { u$plus <- lapply(u$plus, .ec_round_row); u$minus <- lapply(u$minus, .ec_round_row); u })
  }
  list(common = common, unions = unions)
}

# ---------------------------------------------------------------------------
# integer classes on the stacked problem: 6a R19(b) (Lemma E1) and 2e Theorem E27
# ---------------------------------------------------------------------------
.ed_classify <- function(P, rows) {
  kinds <- vapply(rows, function(r) r$kind, "")
  if (any(kinds == "gen"))
    return(list(class = NA_character_, reason = "rows with unequal coefficients (pairwise or W2 signs, pairwise orderings, contrasts, linear rows)"))
  if (!any(kinds == "change")) return(.ec_classify(P, rows))          # 6a classes on the unit-elections (E1)
  # Theorem E27: R = 2; boxes on single cells; change rows on one column; period pairs form a
  # graph of maximum degree 2 without cycles; laminar unit sets per pair
  if (P$R != 2) return(list(class = NA_character_, reason = "change rows with R > 2 (outside 2e Theorem E27)"))
  if (any(kinds == "total")) return(list(class = NA_character_, reason = "change rows mixed with totals within an election (outside 2e Theorem E27)"))
  ch <- rows[kinds == "change"]
  cols <- unique(vapply(ch, function(r) r$cell[2], 0))
  if (length(cols) != 1) return(list(class = NA_character_, reason = "change rows on several columns (outside 2e Theorem E27)"))
  pairs <- unique(t(vapply(ch, function(r) as.numeric(sort(r$pair)), c(0, 0))))
  if (!is.matrix(pairs)) pairs <- matrix(pairs, ncol = 2)
  deg <- table(factor(as.vector(pairs), levels = seq_len(P$T)))
  if (any(deg > 2)) return(list(class = NA_character_, reason = "an election lies in three or more period pairs (degree > 2; not TU, 2e Thm E27 checks)"))
  # cycle check (union-find over elections)
  par <- seq_len(P$T); fnd <- function(i) { while (par[i] != i) i <- par[i]; i }
  for (k in seq_len(nrow(pairs))) {
    a <- fnd(pairs[k, 1]); b <- fnd(pairs[k, 2])
    if (a == b) return(list(class = NA_character_, reason = "the period pairs form a cycle (outside 2e Theorem E27 as stated)"))
    par[a] <- b
  }
  keys <- vapply(ch, function(r) paste(sort(r$pair), collapse = "-"), "")
  for (k in unique(keys)) {
    sets <- lapply(ch[keys == k], function(r) sort(r$units))
    if (!.ec_laminar(sets)) return(list(class = NA_character_, reason = "crossing unit sets within one period pair (outside 2e Theorem E27)"))
  }
  list(class = "E27", reason = "R = 2, change rows on one column, paths of period pairs, laminar unit sets (2e Theorem E27)")
}
