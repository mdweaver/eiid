# The count layer (M-COUNT): eiCounts(), count targets, witness() and the cross-layer
# (R27) block. Theory: docs/extension/PHASE6a_count_layer.md; contract Sec. 6 (C1-C10).

# ---------------------------------------------------------------------------
# targets
# ---------------------------------------------------------------------------

#' Count targets for eiCounts()
#'
#' Constructors for the targets of \code{\link{eiCounts}}: linear functionals of the
#' realised table over a known denominator (PHASE6a Sec. 1.2, contract C2).
#' \itemize{
#'   \item \code{rate(row, col, where)}: \eqn{Y_{rc}(S)/A_r(S)}. \code{row} may name several
#'     rows (a merged group) and \code{col} several columns (a sum over columns).
#'   \item \code{gap(row1, row2, col, where, where2)}: \code{rate(row1, col, where) -
#'     rate(row2, col, where2)}; with \code{row1 == row2} and different unit sets it is a
#'     difference over unit sets.
#'   \item \code{ovr_gap(row, col, where)}: the one-vs-rest gap \eqn{D^r_c = B_{rc} -
#'     B_{\neg r, c}} (6a Lemma C1(a)).
#'   \item \code{share(row, col, where)}: \eqn{Y_{rc}(S)/Y_c(S)}, a share among the
#'     \code{col}-voters. Its denominator uses outcome margins, so it gets the count layer only
#'     (no cross-layer interval; 6a Lemma 4.1).
#'   \item \code{combo(..., weights)}: a linear combination of targets.
#' }
#' \code{where} selects units: an expression evaluated in the data of the
#' \code{eiCounts()} call (e.g. \code{where = county == "Dane"}), a logical vector with one
#' value per row of the data, unit labels, or row numbers. \code{NULL} means all units.
#' For the cross-layer interval the unit set must be fixed from the design (6a Lemma 4.1):
#' an expression that uses an outcome (column) variable makes the target count-layer only,
#' and for a set given as labels, numbers or a vector the printout says that a set chosen
#' using the outcomes voids the cross-layer interval.
#'
#' @param row,row1,row2 row (group) label(s) or index(es).
#' @param col column (outcome) label(s) or index(es).
#' @param where,where2 unit selection (see Details).
#' @param label optional label for printing.
#' @param ... targets to combine.
#' @param weights weights of the combination.
#' @return an object of class \code{"ecTarget"}.
#' @name count_targets
NULL

.ec_tspec <- function(kind, args, label, env) structure(list(kind = kind, args = args, label = label, env = env),
                                                       class = "ecTarget")

#' @rdname count_targets
#' @export
rate <- function(row, col, where = NULL, label = NULL)
  .ec_tspec("rate", list(row = row, col = col, where = substitute(where)), label, parent.frame())

#' @rdname count_targets
#' @export
gap <- function(row1, row2, col, where = NULL, where2 = NULL, label = NULL) {
  w2 <- substitute(where2); w1 <- substitute(where)
  .ec_tspec("gap", list(row1 = row1, row2 = row2, col = col, where = w1,
                        where2 = if (is.null(w2)) w1 else w2), label, parent.frame())
}

#' @rdname count_targets
#' @export
ovr_gap <- function(row, col, where = NULL, label = NULL)
  .ec_tspec("ovr_gap", list(row = row, col = col, where = substitute(where)), label, parent.frame())

#' @rdname count_targets
#' @export
share <- function(row, col, where = NULL, label = NULL)
  .ec_tspec("share", list(row = row, col = col, where = substitute(where)), label, parent.frame())

#' @rdname count_targets
#' @export
combo <- function(..., weights, label = NULL) {
  ts <- list(...)
  if (!all(vapply(ts, inherits, TRUE, "ecTarget"))) stop("combo(): every argument must be a count target")
  if (length(weights) != length(ts)) stop("combo(): one weight per target")
  .ec_tspec("combo", list(targets = ts, weights = as.numeric(weights)), label, parent.frame())
}

# where-expressions are evaluated in the data with the constructor's environment
# where-expressions are evaluated in the data with the constructor's environment. The result
# carries attr "wkind": "all", "design" (an expression using no outcome variable), "outcome" (an
# expression using an outcome (column) variable: not design-measurable, 6a Lemma 4.1, C2) or
# "given" (a vector, labels or row numbers: the package cannot tell how it was chosen).
.ec_where <- function(P, w, env, what) {
  tagged <- function(u, k) structure(u, wkind = k)
  if (is.null(w)) return(tagged(seq_len(P$p), "all"))
  if (inherits(w, "formula")) w <- w[[length(w)]]
  if (is.language(w) && !is.symbol(w)) {
    if (any(all.vars(w) %in% P$cols)) return(tagged(.ec_units(P, where = w, env = env, what = what), "outcome"))
    v <- if (is.null(P$data)) eval(w, env) else eval(w, P$data, env)
    # a predicate on design variables is a design set; an expression that yields labels or row
    # numbers is a given set (how it was chosen is unknown)
    return(tagged(.ec_units(P, units = v, env = env, what = what), if (is.logical(v)) "design" else "given"))
  }
  v <- if (is.symbol(w)) {
    if (!is.null(P$data) && as.character(w) %in% names(P$data)) eval(w, P$data, env) else eval(w, env)
  } else w
  if (inherits(v, "formula")) v <- v[[length(v)]]
  if (is.language(v)) return(.ec_where(P, v, env, what))
  tagged(.ec_units(P, units = v, env = env, what = what), "given")
}

# resolve a target: terms (units, rows, cols, weight), constant, design-measurability and the
# representations available for the residual allowance (6a Remark 4.1a)
.ec_resolve <- function(P, ts, name = NULL) {
  nm <- P$num
  ri <- function(v) .ec_idx(v, P$rows, "row"); ci <- function(v) .ec_idx(v, P$cols, "col")
  rr <- function(r) paste(P$rows[r], collapse = "+"); cc <- function(c) paste(P$cols[c], collapse = "+")
  us <- function(S) if (length(S) == P$p) "" else if (length(S) <= 3) sprintf(" [%s]", paste(P$units[S], collapse = ", ")) else
    sprintf(" [%d units: %s, ...]", length(S), paste(P$units[S][1:2], collapse = ", "))
  A_of <- function(S, rws) sum(nm(P$A[S, rws, drop = FALSE]))
  term <- function(S, rws, cls, w) list(units = S, rows = rws, cols = cls, w = w, wd = as.double(w))
  a <- ts$args; out <- NULL
  if (ts$kind == "rate") {
    S <- .ec_where(P, a$where, ts$env, "rate"); r <- ri(a$row); c <- ci(a$col)
    den <- A_of(S, r)
    if (as.double(den) <= 0) stop("rate(): the denominator A_r(S) is zero")
    out <- list(terms = list(term(S, r, c, 1 / den)), const = nm(0), design = TRUE,
                label = paste0(rr(r), " -> ", cc(c), us(S)), den = sprintf("A_%s = %s", rr(r), format(as.double(den))))
    out$reps <- list(list(terms = out$terms, const = out$const))
  } else if (ts$kind == "gap") {
    S1 <- .ec_where(P, a$where, ts$env, "gap"); S2 <- .ec_where(P, a$where2, ts$env, "gap")
    r1 <- ri(a$row1); r2 <- ri(a$row2); c <- ci(a$col)
    d1 <- A_of(S1, r1); d2 <- A_of(S2, r2)
    if (as.double(d1) <= 0 || as.double(d2) <= 0) stop("gap(): a denominator is zero")
    out <- list(terms = list(term(S1, r1, c, 1 / d1), term(S2, r2, c, -1 / d2)), const = nm(0), design = TRUE,
                label = paste0(rr(r1), us(S1), " - ", rr(r2), us(S2), " on ", cc(c)),
                den = sprintf("A_%s = %s, A_%s = %s", rr(r1), format(as.double(d1)), rr(r2), format(as.double(d2))))
    out$reps <- list(list(terms = out$terms, const = out$const))
  } else if (ts$kind == "ovr_gap") {
    S <- .ec_where(P, a$where, ts$env, "ovr_gap"); r <- ri(a$row); c <- ci(a$col)
    if (length(r) != 1 || length(c) != 1) stop("ovr_gap(): one row and one column")
    Ar <- A_of(S, r); N <- sum(nm(P$N[S])); Yc <- sum(nm(P$B[S, c]))
    if (as.double(Ar) <= 0 || as.double(Ar) >= as.double(N)) stop("ovr_gap(): needs 0 < A_r < N over the unit set")
    rest <- setdiff(seq_len(P$R), r)
    # representation over cells (no constant), and kappa Y_rc - Y_c / (N - A_r) (Lemma C1(a))
    rep1 <- list(terms = list(term(S, r, c, 1 / Ar), term(S, rest, c, -1 / (N - Ar))), const = nm(0))
    rep2 <- list(terms = list(term(S, r, c, N / (Ar * (N - Ar)))), const = -Yc / (N - Ar))
    out <- list(terms = rep1$terms, const = rep1$const, design = TRUE, reps = list(rep1, rep2),
                label = paste0(rr(r), " vs rest on ", cc(c), us(S)),
                den = sprintf("A_%s = %s, N - A_%s = %s", rr(r), format(as.double(Ar)), rr(r), format(as.double(N - Ar))))
  } else if (ts$kind == "share") {
    S <- .ec_where(P, a$where, ts$env, "share"); r <- ri(a$row); c <- ci(a$col)
    Yc <- sum(nm(P$B[S, c, drop = FALSE]))
    if (as.double(Yc) <= 0) stop("share(): the denominator Y_c(S) is zero")
    out <- list(terms = list(term(S, r, c, 1 / Yc)), const = nm(0), design = FALSE,
                label = paste0(rr(r), " among ", cc(c), us(S)), den = sprintf("Y_%s = %s (outcome margin)", cc(c), format(as.double(Yc))))
    out$reps <- list(list(terms = out$terms, const = out$const))
  } else if (ts$kind == "combo") {
    subs <- lapply(a$targets, function(t) .ec_resolve(P, t))
    terms <- list(); const <- nm(0)
    for (j in seq_along(subs)) {
      wj <- nm(a$weights[j])
      for (tm in subs[[j]]$terms) terms[[length(terms) + 1]] <- term(tm$units, tm$rows, tm$cols, tm$w * wj)
      const <- const + subs[[j]]$const * wj
    }
    out <- list(terms = terms, const = const, design = all(vapply(subs, function(s) s$design, TRUE)),
                label = paste(sprintf("%g*(%s)", a$weights, vapply(subs, function(s) s$label, "")), collapse = " + "),
                den = "combination")
    out$reps <- list(list(terms = terms, const = const))
  } else stop("unknown target kind")
  # design-measurability of the unit sets (code review issue 7)
  kinds <- if (ts$kind == "combo") unlist(lapply(subs, function(s) s$wkinds)) else
    unlist(lapply(out$terms, function(tm) attr(tm$units, "wkind")))
  out$wkinds <- unique(kinds)
  if ("outcome" %in% out$wkinds) { out$design <- FALSE; out$design_reason <- "outcome_set" }
  else if (!out$design && is.null(out$design_reason)) out$design_reason <- "outcome_margin"
  out$given_set <- "given" %in% out$wkinds
  if (!is.null(ts$label)) out$label <- ts$label else if (!is.null(name) && nzchar(name)) out$label <- name
  out$obj <- .ec_obj(P, out$terms)
  out
}

.ec_obj <- function(P, terms) {
  obj <- numeric(P$nv)
  for (tm in terms) {
    ks <- as.vector(outer(outer((tm$units - 1L) * P$R * P$C, (tm$rows - 1L) * P$C, "+"), tm$cols, "+"))
    obj[ks] <- obj[ks] + tm$wd
  }
  obj
}

# value of a target at table x; exact rational arithmetic for integer x with integer margins
.ec_value <- function(P, tg, x) {
  sums <- vapply(tg$terms, function(tm) {
    ks <- as.vector(outer(outer((tm$units - 1L) * P$R * P$C, (tm$rows - 1L) * P$C, "+"), tm$cols, "+"))
    sum(x[ks]) }, 0)
  if (P$integer && all(x == round(x))) {
    v <- tg$const
    for (j in seq_along(tg$terms)) v <- v + tg$terms[[j]]$w * gmp::as.bigq(sums[j])
    return(as.double(v))
  }
  sum(vapply(tg$terms, function(tm) tm$wd, 0) * sums) + as.double(tg$const)
}

# eiBounds-style targets (named c(row, col) pairs or R x C weight matrices over B_rc)
.ec_from_matrix <- function(P, a, name) {
  terms <- list()
  for (r in seq_len(P$R)) for (c in seq_len(P$C)) if (a[r, c] != 0) {
    den <- sum(P$num(P$A[, r]))
    if (as.double(den) <= 0) stop("target '", name, "': row ", P$rows[r], " has no members")
    w <- P$num(a[r, c]) / den
    terms[[length(terms) + 1]] <- list(units = seq_len(P$p), rows = r, cols = c, w = w, wd = as.double(w))
  }
  out <- list(terms = terms, const = P$num(0), design = TRUE, label = name, den = "A_r (all units)")
  out$reps <- list(list(terms = terms, const = out$const))
  out$obj <- .ec_obj(P, terms)
  out
}

.ec_targets <- function(P, targets) {
  if (is.null(targets)) {
    out <- list()
    for (r in seq_len(P$R)) for (c in seq_len(P$C)) {
      if (sum(P$A[, r]) == 0) next                  # a row with no members has no rate
      a <- matrix(0, P$R, P$C); a[r, c] <- 1
      out[[length(out) + 1]] <- .ec_from_matrix(P, a, paste0(P$rows[r], " -> ", P$cols[c]))
    }
    return(out)
  }
  if (inherits(targets, "ecTarget")) targets <- list(targets)
  nms <- names(targets); if (is.null(nms)) nms <- rep("", length(targets))
  out <- lapply(seq_along(targets), function(j) {
    tg <- targets[[j]]
    if (inherits(tg, "ecTarget")) return(.ec_resolve(P, tg, nms[j]))
    if (is.matrix(tg)) {
      if (!all(dim(tg) == c(P$R, P$C))) stop("target matrices must be ", P$R, " x ", P$C)
      return(.ec_from_matrix(P, tg, if (nzchar(nms[j])) nms[j] else paste0("target", j)))
    }
    a <- matrix(0, P$R, P$C); a[.ec_idx(tg[1], P$rows, "row"), .ec_idx(tg[2], P$cols, "col")] <- 1
    .ec_from_matrix(P, a, if (nzchar(nms[j])) nms[j] else paste0(P$rows[a != 0][1], " -> ", tg[2]))
  })
  # distinct labels (witness(target = label) must be able to tell targets apart; issue 14)
  labs <- make.unique(vapply(out, function(t) t$label, ""), sep = " #")
  for (j in seq_along(out)) out[[j]]$label <- labs[j]
  out
}

# ---------------------------------------------------------------------------
# restrictions
# ---------------------------------------------------------------------------
.ec_aliases <- c(assoc = "sign", order = "order_unit", known_cell = "known", zeros = "zero")

# a named list of specs (one spec, or a list of specs, per type) -> parsed specs with ids
.ec_parse_restrictions <- function(P, restrictions, env) {
  if (is.null(restrictions) || !length(restrictions)) return(list())
  if (is.null(names(restrictions)) || any(names(restrictions) == "")) stop("'restrictions' must be a named list")
  out <- list()
  for (j in seq_along(restrictions)) {
    type <- names(restrictions)[j]
    if (type %in% names(.ec_aliases)) type <- .ec_aliases[[type]]
    v <- restrictions[[j]]
    if (is.null(v) || !length(v)) next
    specs <- if (is.list(v) && length(v) && all(vapply(v, is.list, TRUE))) v else list(v)
    for (s in specs) {
      ps <- .ec_parse_spec(P, type, s, env)
      k <- sum(vapply(out, function(o) o$type == type, TRUE)) + 1
      ps$id <- paste0(type, "[", k, "]")
      ps2 <- ps; if (type == "sign") ps2$type <- ps$type_a       # label needs dB / dW (issue 13)
      ps$label <- .ec_label_of(P, ps2, type)
      ps$experimental <- type %in% c("sign", "cr")
      out[[length(out) + 1]] <- ps
    }
  }
  out
}

.ec_build_rows <- function(P, specs) {
  common <- list(); unions <- list()
  for (s in specs) {
    rs <- .ec_spec_rows(P, s, s$id)
    common <- c(common, rs$rows)
    if (!is.null(rs$union)) unions[[length(unions) + 1]] <- c(rs$union, list(id = s$id))
  }
  if (P$integer) {
    common <- lapply(common, .ec_round_row)
    unions <- lapply(unions, function(u) { u$plus <- lapply(u$plus, .ec_round_row); u$minus <- lapply(u$minus, .ec_round_row); u })
  }
  list(common = common, unions = unions)
}

# ---------------------------------------------------------------------------
# ranges
# ---------------------------------------------------------------------------
# margins only: per-unit closed forms (Lemma 2.1(a), (b)) when each unit's objective has one
# nonzero cell, or two cells of one column with opposite signs; NULL otherwise
.ec_fast <- function(P, tg) {
  ob <- array(tg$obj, c(P$C, P$R, P$p))
  xlo <- xhi <- numeric(P$nv)
  for (i in seq_len(P$p)) {
    v <- matrix(ob[, , i], P$R, P$C, byrow = TRUE)
    nz <- which(v != 0, arr.ind = TRUE)
    if (!nrow(nz)) next
    a <- P$A[i, ]; b <- P$B[i, ]; N <- P$N[i]
    lo1 <- function(r, c) max(0, a[r] + b[c] - N); hi1 <- function(r, c) min(a[r], b[c])
    k <- function(r, c) (i - 1L) * P$R * P$C + (r - 1L) * P$C + c
    if (nrow(nz) == 1) {
      r <- nz[1, 1]; c <- nz[1, 2]
      if (v[r, c] > 0) { xhi[k(r, c)] <- hi1(r, c); xlo[k(r, c)] <- lo1(r, c) } else
        { xhi[k(r, c)] <- lo1(r, c); xlo[k(r, c)] <- hi1(r, c) }
    } else if (nrow(nz) == 2 && nz[1, 2] == nz[2, 2] && sign(v[nz[1, 1], nz[1, 2]]) != sign(v[nz[2, 1], nz[2, 2]])) {
      c <- nz[1, 2]; rp <- nz[v[nz] > 0, 1]; rm <- nz[v[nz] < 0, 1]
      # Lemma 2.1(b): both separate optima are attained simultaneously
      xhi[k(rp, c)] <- hi1(rp, c); xhi[k(rm, c)] <- lo1(rm, c)
      xlo[k(rp, c)] <- lo1(rp, c); xlo[k(rm, c)] <- hi1(rm, c)
    } else return(NULL)
  }
  list(lower = .ec_value(P, tg, xlo), upper = .ec_value(P, tg, xhi))
}

# one system (no unions): continuous and integer endpoints with witnesses
.ec_system <- function(P, rows, tg, mode) {
  out <- list(cont = FALSE, int_ok = NA, lp = c(NA_real_, NA_real_), iv = c(NA_real_, NA_real_),
              x = list(NULL, NULL), class = NA_character_, milp = FALSE, fallback = FALSE)
  integer_mode <- P$integer && mode != "lp"
  if (mode != "exact" || !P$integer) {
    clo <- .ec_optimise(P, rows, tg$obj, -1)
    if (!clo$ok) return(out)
    chi <- .ec_optimise(P, rows, tg$obj, 1)
    out$cont <- TRUE
    out$lp <- c(.ec_value(P, tg, clo$x), .ec_value(P, tg, chi$x))
    if (!integer_mode) { out$iv <- out$lp; out$x <- list(clo$x, chi$x); return(out) }
  }
  cls <- .ec_classify(P, rows)
  out$class <- cls$class; out$class_reason <- cls$reason
  use_milp <- mode == "exact" || is.na(cls$class)
  if (mode == "exact" && !out$cont) {
    clo <- .ec_optimise(P, rows, tg$obj, -1)
    if (clo$ok) { chi <- .ec_optimise(P, rows, tg$obj, 1); out$cont <- TRUE
      out$lp <- c(.ec_value(P, tg, clo$x), .ec_value(P, tg, chi$x)) }
  }
  # the LP over the rounded rows: an outer bound for the integer values (R19(c) item 1), exact
  # in the TU classes (R19(b)); infeasible => no integer table (F^Z is contained in F_rd)
  same <- out$cont && all(vapply(rows, function(rw) isTRUE(rw$rounded) && identical(rw$coef_rd, rw$coef) &&
                                    identical(rw$lo_rd, rw$lo) && identical(rw$hi_rd, rw$hi), TRUE))
  rl <- if (same) list(clo, chi) else                  # rounding changed nothing: reuse the LP
    list(.ec_optimise(P, rows, tg$obj, -1, rd = TRUE), .ec_optimise(P, rows, tg$obj, 1, rd = TRUE))
  if (!isTRUE(rl[[1]]$ok) || !isTRUE(rl[[2]]$ok)) { out$int_ok <- FALSE; return(out) }
  out$rd <- c(.ec_value(P, tg, rl[[1]]$x), .ec_value(P, tg, rl[[2]]$x))
  sol <- list(NULL, NULL); unres <- character(0)
  for (j in 1:2) {
    sense <- c(-1, 1)[j]; x <- NULL
    if (!use_milp) {
      # accept the vertex only if it is integral (before rounding), passes the exact check and
      # has the LP optimum's value (code review issue 6); otherwise fall back to the MILP
      sx <- rl[[j]]$x
      if (max(abs(sx - round(sx))) <= 1e-6) {
        xr <- round(sx)
        if (.ec_check(P, rows, xr, TRUE) &&
            abs(sum(tg$obj * xr) - rl[[j]]$value) <= 1e-9 * max(1, abs(rl[[j]]$value))) x <- xr
      }
      if (is.null(x)) out$fallback <- TRUE
    }
    if (is.null(x)) {
      out$milp <- TRUE
      m <- .ec_milp(P, rows, tg, sense, P$time_limit)
      if (m$status == "infeasible") { out$int_ok <- FALSE; return(out) }
      if (m$status == "optimal") x <- m$x else { unres <- c(unres, m$reason); x <- NULL }
    }
    sol[j] <- list(x)
  }
  if (length(unres)) {
    # integer status unresolved: report the rounded-LP outer bound at the unresolved end(s)
    out$unresolved <- paste(unique(unres), collapse = "; ")
    out$int_ok <- NA
    out$x <- lapply(1:2, function(j) if (is.null(sol[[j]])) rl[[j]]$x else sol[[j]])
    out$iv <- vapply(1:2, function(j) .ec_value(P, tg, out$x[[j]]), 0)
    return(out)
  }
  out$int_ok <- TRUE
  out$iv <- c(.ec_value(P, tg, sol[[1]]), .ec_value(P, tg, sol[[2]]))
  out$x <- sol
  out
}

# all systems (CR branches), their union, and the integer status
.ec_range <- function(P, built, tg, mode, keep_x = TRUE) {
  U <- length(built$unions)
  combos <- if (U) as.matrix(expand.grid(rep(list(1:2), U))) else matrix(0L, 1, 0)
  if (U == 0 && !length(built$common)) {
    fs <- .ec_fast(P, tg)
    if (!is.null(fs)) {
      cls <- if (P$integer && mode != "lp") "(i)" else NA_character_
      return(list(lower = fs$lower, upper = fs$upper, lp = c(fs$lower, fs$upper), cont = TRUE,
                  int_ok = if (P$integer && mode != "lp") TRUE else NA, pieces = matrix(c(fs$lower, fs$upper), 1),
                  status = .ec_status(P, mode, FALSE, cls), class = cls, x = list(NULL, NULL),
                  branches = NULL, milp = FALSE))
    }
  }
  res <- list(); labs <- character(0)
  for (s in seq_len(nrow(combos))) {
    rows <- built$common
    lab <- character(0)
    for (u in seq_len(U)) {
      br <- if (combos[s, u] == 1) "plus" else "minus"
      rows <- c(rows, built$unions[[u]][[br]])
      lab <- c(lab, paste0(built$unions[[u]]$id, if (br == "plus") "+" else "-"))
    }
    res[[s]] <- .ec_system(P, rows, tg, mode)
    labs[s] <- paste(lab, collapse = ", ")
  }
  feas_c <- vapply(res, function(r) r$cont, TRUE)
  integer_mode <- P$integer && mode != "lp"
  unres <- unique(unlist(lapply(res, function(r) r$unresolved)))
  ok <- if (integer_mode) vapply(res, function(r) isTRUE(r$int_ok) || !is.null(r$unresolved), TRUE) else feas_c
  milp <- any(vapply(res, function(r) r$milp, TRUE))
  cls <- if (integer_mode && !milp) unique(stats::na.omit(vapply(res, function(r) r$class, "")))[1] else NA_character_
  int_ok <- if (!integer_mode) NA else if (any(vapply(res, function(r) isTRUE(r$int_ok), TRUE))) TRUE else
    if (length(unres)) NA else FALSE
  out <- list(cont = any(feas_c), int_ok = int_ok, milp = milp, class = cls, unresolved = if (length(unres)) paste(unres, collapse = "; "),
              lp = if (any(feas_c)) c(min(vapply(res[feas_c], function(r) r$lp[1], 0)),
                                       max(vapply(res[feas_c], function(r) r$lp[2], 0))) else c(NA_real_, NA_real_),
              fallback = any(vapply(res, function(r) isTRUE(r$fallback), TRUE)))
  out$status <- .ec_status(P, mode, milp, cls, out$unresolved)
  if (U) out$branches <- data.frame(branch = labs, continuous = feas_c,
                                    lower = vapply(res, function(r) r$iv[1], 0),
                                    upper = vapply(res, function(r) r$iv[2], 0), stringsAsFactors = FALSE)
  if (!any(ok)) { out$lower <- out$upper <- NA_real_; out$pieces <- NULL; out$x <- list(NULL, NULL); return(out) }
  ivs <- lapply(res[ok], function(r) r$iv)
  out$pieces <- .rxc_union(ivs, tol = 1e-12)
  jl <- which(ok)[which.min(vapply(res[ok], function(r) r$iv[1], 0))]
  ju <- which(ok)[which.max(vapply(res[ok], function(r) r$iv[2], 0))]
  out$lower <- res[[jl]]$iv[1]; out$upper <- res[[ju]]$iv[2]
  out$x <- if (keep_x) list(res[[jl]]$x[[1]], res[[ju]]$x[[2]]) else list(NULL, NULL)
  out$sys <- c(labs[jl], labs[ju])
  out
}

.ec_status <- function(P, mode, milp, cls, unresolved = NULL) {
  if (!P$integer && mode != "lp") return("not applicable: non-integer margins")
  if (mode == "lp") return("continuous only")
  if (!is.null(unresolved))
    return(paste0("integer status unresolved (", unresolved, "): LP outer bound (R19(c) item 1)"))
  if (milp) "integer-exact (MILP)" else "integer-exact (TU class)"
}

# drop-one (added-last) test: a restriction binds an endpoint iff removing it moves the
# endpoint outward; equivalently it removes the whole optimal face (PD R2 face test)
.ec_binding <- function(P, specs, built_all, tg, rg, mode, tol = 1e-9) {
  if (!length(specs) || is.na(rg$lower)) return(c("", ""))
  bl <- character(0); bu <- character(0)
  # LP duality: a restriction whose rows are all slack at the optimal LP vertex has zero duals,
  # so removing it leaves that vertex optimal; only active restrictions are refitted. (Not
  # used for MILP or unresolved results, nor for CR unions.)
  lp_path <- !isTRUE(rg$milp) && is.null(rg$unresolved) && !length(built_all$unions)
  rd <- P$integer && mode != "lp"
  active <- function(s, x) {
    if (is.null(x)) return(TRUE)
    for (rw in Filter(function(r) r$id == s$id, built_all$common)) {
      cf <- if (rd && !is.null(rw$coef_rd)) rw$coef_rd else rw$coef
      lo <- if (rd && !is.null(rw$lo_rd)) rw$lo_rd else rw$lo; hi <- if (rd && !is.null(rw$hi_rd)) rw$hi_rd else rw$hi
      v <- sum(cf * x[rw$idx]); sc <- 1e-7 * max(1, sum(abs(cf * x[rw$idx])))
      if ((is.finite(lo) && v <= lo + sc) || (is.finite(hi) && v >= hi - sc)) return(TRUE)
    }
    FALSE
  }
  for (s in specs) {
    if (lp_path && !active(s, rg$x[[1]]) && !active(s, rg$x[[2]])) next
    b <- list(common = Filter(function(r) r$id != s$id, built_all$common),
              unions = Filter(function(u) u$id != s$id, built_all$unions))
    r0 <- .ec_range(P, b, tg, mode, keep_x = FALSE)
    sc <- tol * max(1, abs(rg$lower), abs(rg$upper))
    if (!is.na(r0$lower) && r0$lower < rg$lower - sc) bl <- c(bl, s$id)
    if (!is.na(r0$upper) && r0$upper > rg$upper + sc) bu <- c(bu, s$id)
  }
  c(paste(bl, collapse = ", "), paste(bu, collapse = ", "))
}

# ---------------------------------------------------------------------------
# data input
# ---------------------------------------------------------------------------
.ec_input <- function(formula, data, rows, cols, x, t, n, unit) {
  if (!is.null(formula)) {
    if (is.null(data)) stop("'data' is required with a formula")
    cols <- all.vars(formula[[2]]); rows <- all.vars(formula[[3]])
  }
  labs <- NULL
  if (!is.null(unit)) labs <- if (is.character(unit) && length(unit) == 1 && !is.null(data) && unit %in% names(data))
    data[[unit]] else unit
  if (!is.null(rows) && !is.null(cols)) {
    if (is.null(data)) stop("'data' is required with 'rows' and 'cols'")
    A <- as.matrix(data[, rows, drop = FALSE]); B <- as.matrix(data[, cols, drop = FALSE])
    int_in <- TRUE
  } else if (!is.null(x) && !is.null(t) && !is.null(n)) {
    X <- if (is.null(dim(x))) { m <- cbind(x, 1 - x); colnames(m) <- c("X", "notX"); m } else as.matrix(x)
    T <- if (is.null(dim(t))) { m <- cbind(t, 1 - t); colnames(m) <- c("T", "notT"); m } else as.matrix(t)
    rows <- if (is.null(colnames(X))) { if (is.null(dim(x))) c("X", "notX") else paste0("row", seq_len(ncol(X))) } else colnames(X)
    cols <- if (is.null(colnames(T))) { if (is.null(dim(t))) c("T", "notT") else paste0("col", seq_len(ncol(T))) } else colnames(T)
    A <- X * as.numeric(n); B <- T * as.numeric(n)
    int_in <- FALSE                     # shares: N * X recomputed in floating point is not a count
  } else stop("supply a formula and data, 'data' with 'rows' and 'cols', or 'x', 't' and 'n'")
  if (is.null(labs)) labs <- if (!is.null(data) && !is.null(rownames(data)) &&
                                 !identical(rownames(data), as.character(seq_len(nrow(A))))) rownames(data) else seq_len(nrow(A))
  if (length(labs) != nrow(A)) stop("'unit' has the wrong length")
  list(A = A, B = B, rows = rows, cols = cols, units = labs, integer_input = int_in)
}

# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------
.ec_fit <- function(P, targets, restrictions, integer, binding, env, time_limit = 10) {
  P$time_limit <- time_limit
  specs <- .ec_parse_restrictions(P, restrictions, env)
  if (integer == "exact" && !P$integer)
    stop("integer = \"exact\" needs exactly integer margins; these margins are not integers (or were given ",
         "as shares). The integer check does not apply (6a R19(e)): use integer = \"auto\" or \"lp\"")
  built <- .ec_build_rows(P, specs)
  tl <- .ec_targets(P, targets)
  empty <- list(common = list(), unions = list())
  rows_out <- list(); wit <- list()
  for (j in seq_along(tl)) {
    tg <- tl[[j]]
    dd <- .ec_range(P, empty, tg, integer, keep_x = FALSE)
    rg <- if (length(specs)) .ec_range(P, built, tg, integer) else dd
    bd <- if (isTRUE(binding)) .ec_binding(P, specs, built, tg, rg, integer) else c(NA, NA)
    st <- rg$status
    if (!isTRUE(rg$cont)) st <- "infeasible: the rows contradict the margins (continuous tables)"
    else if (identical(rg$int_ok, FALSE)) st <- "integer-infeasible: no whole-number table satisfies these rows"
    rows_out[[j]] <- data.frame(
      target = tg$label, denominator = tg$den,
      dd_lower = dd$lower, dd_upper = dd$upper,
      lower = rg$lower, upper = rg$upper, width = rg$upper - rg$lower,
      integer_status = st, tu_class = if (is.null(rg$class)) NA_character_ else rg$class,
      lp_lower = rg$lp[1], lp_upper = rg$lp[2],
      continuous_feasible = isTRUE(rg$cont), integer_feasible = rg$int_ok,
      n_pieces = if (is.null(rg$pieces)) 0L else nrow(rg$pieces),
      binding_lower = bd[1], binding_upper = bd[2], design = tg$design,
      stringsAsFactors = FALSE)
    wit[[j]] <- list(x = rg$x, sys = rg$sys, pieces = rg$pieces, branches = rg$branches,
                     dd_status = dd$status, fallback = isTRUE(rg$fallback))
  }
  tab <- do.call(rbind, rows_out); rownames(tab) <- NULL
  rt <- if (length(specs)) data.frame(
    id = vapply(specs, function(s) s$id, ""), restriction = vapply(specs, function(s) s$label, ""),
    object = "actual", role = vapply(specs, function(s) s$role, ""),
    experimental = vapply(specs, function(s) s$experimental, TRUE), stringsAsFactors = FALSE) else NULL
  list(table = tab, targets = tl, specs = specs, built = built, restrictions = rt, solutions = wit)
}

#' Exact bounds on actual counts and rates (the count layer, M-COUNT)
#'
#' \code{eiCounts()} bounds linear functionals of the \emph{realised} unit tables (actual
#' counts, rates, gaps, one-vs-rest gaps, rates over subsets of units) using only the unit
#' margins and stated rows on the realised table: known cells, structural zeros, cell boxes,
#' published totals, realised orderings, the odds-ratio dial, realised association signs and
#' contextual reinforcement (CR). No model, regression or sampling theory is involved. Theory:
#' \code{docs/extension/PHASE6a_count_layer.md} (R19, R27) in the source repository.
#'
#' \strong{What is computed.} The range of each target over all tables that reproduce the
#' margins and satisfy every row (R19(a)): a linear program (LP), sharp for continuous tables.
#' The first line of each target is the margins-only range (Duncan-Davis, computed as the
#' exact LP; for multi-cell targets it is narrower than cellwise interval arithmetic,
#' 6a Lemma 2.1(c)).
#'
#' \strong{Integer tables} (\code{integer}). With exactly integer margins (counts in the
#' data; shares are not counts) \code{"auto"} rounds every row inward in exact rational
#' arithmetic (package \pkg{gmp}; 6a Sec. 1.4) and checks whether the \emph{whole} row set lies in
#' one total-unimodularity class of R19(b): (i) totals on one cell whose unit sets form two
#' laminar families (decided by 2-colouring the crossing graph); (ii) totals on cells of one
#' column (or row), one laminar family per cell; (iii) \eqn{R = 2} (or \eqn{C = 2}), one
#' laminar family per column (row). Single-cell rows (known cells, zeros, boxes, one-vs-rest
#' orderings, the odds-ratio dial) and one-vs-rest (W1) sign and CR rows, which are totals
#' (6a Lemma C1), fit the classes; then the LP over the rounded rows returns integer
#' endpoints with integer witnesses ("integer-exact (TU class)"); a vertex is accepted only if
#' it is integral before rounding and passes the exact check, otherwise the MILP is used.
#' Otherwise (pairwise orderings or signs, W2 signs, \code{linear} rows, mixtures of classes,
#' crossing families) each endpoint is a mixed-integer program ("integer-exact (MILP)"),
#' solved by GLPK (package \pkg{Rglpk}) with a time limit per program
#' (\code{time_limit}), the objective scaled to coprime integers from the exact weights, an
#' exact check of the solution, and a certificate of optimality: a second program with
#' \eqn{c'n \ge v^* + 1} must be infeasible. The certificate is solved in the shifted variables
#' \eqn{d = n - n^*}, so the cut reads \eqn{c'd \ge 1} and GLPK's relative tolerances stay far
#' below its unit gap; a table it returns must improve on \eqn{v^*} exactly. The scaled integer
#' objective coefficients must not exceed \eqn{10^7} in absolute value (and the objective range
#' \eqn{10^{12}}, so values are exact doubles); above the cap the endpoint is reported as
#' unresolved. On timeout, above the cap, or when the certificate cannot be completed, the endpoint is the LP over the rounded rows, an outer
#' bound, with the status "integer status unresolved (...): LP outer bound"; it is never
#' labelled integer-exact. The rows are rounded, never the endpoints (check T5). Doubles typed
#' by the user (bounds, coefficients, slack, \code{psi}) are rationalised: a fraction with
#' denominator at most \eqn{10^6} within \eqn{10^{-12}} relative is used when one exists (so
#' \code{1/3} means 1/3), and the printout lists the values changed; character fractions
#' (\code{"1/3"}) in \code{linear} rows are exact. With non-integer margins (e.g. estimated counts) the integer check
#' does not apply and the range is exact only for the stated margins (R19(e)).
#' \code{"lp"} gives the continuous LP only; \code{"exact"} always solves the MILP.
#'
#' \strong{Restrictions} are a named list; each element is one specification or a list of
#' them. Every specification may carry \code{object = "actual"} (the default; \code{"model"}
#' is refused: model rows belong to \code{eiBounds()}) and \code{role = "identity"} or
#' \code{"belief"}. Defaults: identity for \code{known} and for a \code{total} given by
#' \code{value} (its right-hand side is the realised value); belief for a \code{total} given
#' as a floor or interval (\code{lower}/\code{upper}: a stated floor is a belief, 6a Sec. 4.1
#' and C8; set \code{role = "identity"} for a published bracket) and for every other row. Units are given by \code{units} (labels or row numbers),
#' \code{unit} (alias), or \code{where} (a logical vector, or a formula / quoted expression
#' evaluated in \code{data}, e.g. \code{where = ~ county == "Dane"}).
#' \itemize{
#'   \item \code{known = list(unit, row, col, value)}: audited cell counts.
#'   \item \code{zero = list(row, col, units)}: structural zeros.
#'   \item \code{box = list(row, col, units, lower, upper)}: cell bounds (counts).
#'   \item \code{total = list(units or where, row, col, value)} or \code{lower}/\code{upper}:
#'     a published total over a unit set.
#'   \item \code{order_unit = list(row, over = "rest" or a row, col, units, sign = 1, slack = 0)}:
#'     the realised ordering in every unit (row at least as supportive as \code{over}); the
#'     one-vs-rest version is a box, the pairwise version (\eqn{R \ge 3}) a two-cell row.
#'   \item \code{odds = list(psi = c(lower, upper), units)}: the 2 x 2 odds-ratio dial (a box;
#'     ignored in units with a zero margin). A dial containing 1 is labelled untestable: with
#'     continuous tables and no other rows the equal-rates table (odds ratio 1) satisfies it,
#'     so the margins alone do not contradict it (with whole people, or with other rows, they
#'     can).
#'   \item \code{sign = list(type = "dB" or "dW", row, col, sign, version, other, where, slack)}:
#'     realised association signs (2c Definition 1 with realised rates); \code{"W1"}
#'     (one-vs-rest, default) is a total on \eqn{Y_{rc}(S)}; \code{"pair"} and \code{"W2"}
#'     have unequal coefficients. \code{assoc} is accepted as an alias. Slack in the same
#'     probability units as \code{eiBounds()}, but the scale (and every design quantity) is
#'     computed over the sign's own unit set S, whereas the model layer uses all units.
#'     \code{sign} may be +1/-1 or \code{">=0"}/\code{"<=0"}. An association whose coefficients
#'     are all zero over S (e.g. X_r constant there) is vacuous and dropped with a note.
#'   \item \code{cr = list(row, col, version, other, where, slack)}: realised contextual
#'     reinforcement, a union of two row systems; the result may be two intervals.
#' }
#' Sign and CR rows print EXPERIMENTAL (their theory rests on 2c M2, review batch B).
#'
#' \strong{Cross-layer block} (R27). With \code{model = } an \code{eiBounds()} fit on the
#' same units, each design-measurable target gets an actual-target model interval: the
#' projection interval (P3.5(b2)(ii); conservative interval P3.5(a) with \code{ci =
#' "conservative"}) for the model-implied value plus the residual allowance of 6a Lemma 4.1:
#' Popoviciu \eqn{S_\varphi} (default) or Hoeffding \eqn{t_E(\varphi)}. The count set and the
#' model interval are intersected; count rows and model rows never share one LP. With
#' identity rows only (margins, published totals, audited known cells) the intersection
#' covers at the model interval's level and disjointness rejects the model and its rows at
#' level alpha; with belief rows the printout says what holds instead (contract C8).
#'
#' @param formula \code{cbind(column counts) ~ cbind(row counts)}.
#' @param data data frame with the counts (and variables used in \code{where}).
#' @param rows,cols names of the row-count and column-count variables (alternative to
#'   \code{formula}).
#' @param x,t,n shares and sizes (alternative input; the integer check is then off).
#' @param unit name of a unit-label variable in \code{data}, or a vector of labels.
#' @param targets a list of count targets (\code{\link{rate}}, \code{\link{gap}},
#'   \code{\link{ovr_gap}}, \code{\link{share}}, \code{\link{combo}}), or as in
#'   \code{eiBounds()} (named \code{c(row, col)} pairs or R x C weight matrices over the
#'   rates); \code{NULL} for every cell rate.
#' @param restrictions named list of count restrictions (see Details).
#' @param integer \code{"auto"}, \code{"lp"} or \code{"exact"} (see Details).
#' @param model optional \code{eiBounds()} fit for the cross-layer block.
#' @param allowance residual allowance for the cross-layer block: \code{"popoviciu"}
#'   (\eqn{S_\varphi}, Lemma 4.1(a)) or \code{"hoeffding"} (\eqn{t_E(\varphi)}, Lemma 4.1(b)).
#' @param binding compute the binding rows (drop-one test) for each endpoint: \code{TRUE},
#'   \code{FALSE}, or \code{"auto"} (default: on for at most 200 units, off above, with a note).
#'   Only restrictions active at the optimal LP vertex are refitted (an inactive one has zero
#'   duals and cannot bind); rows that bind only jointly are not shown.
#' @param time_limit seconds allowed per integer program (MILP; solver GLPK via package
#'   \pkg{Rglpk}). On timeout the endpoint is the LP over the rounded rows, an outer bound,
#'   with the status "integer status unresolved".
#' @return an object of class \code{"eiCounts"}; \code{$table} has one row per target with the
#'   margins-only range (\code{dd_*}), the range under the rows, the integer status, the
#'   continuous LP range (\code{lp_*}) and the binding rows; \code{$cross} the cross-layer
#'   block. Use \code{\link{witness}} for tables attaining each endpoint.
#' @examples
#' d <- data.frame(precinct = c("P1", "P2"), vet = c(20, 60), other = c(80, 40),
#'                 yes = c(30, 50), no = c(70, 50))
#' fit <- eiCounts(cbind(yes, no) ~ cbind(vet, other), data = d, unit = "precinct",
#'                 targets = list(rate("vet", "yes"), gap("vet", "other", "yes")),
#'                 restrictions = list(order_unit = list(row = "vet", col = "yes")))
#' fit
#' witness(fit)
#' @export
eiCounts <- function(formula = NULL, data = NULL, rows = NULL, cols = NULL, x = NULL, t = NULL, n = NULL,
                     unit = NULL, targets = NULL, restrictions = NULL,
                     integer = c("auto", "lp", "exact"), model = NULL,
                     allowance = c("popoviciu", "hoeffding"), binding = "auto", time_limit = 10) {
  integer <- match.arg(integer); allowance <- match.arg(allowance)
  inp <- .ec_input(formula, data, rows, cols, x, t, n, unit)
  P <- .ec_problem(inp$A, inp$B, inp$rows, inp$cols, inp$units, inp$integer_input, data = data)
  if (P$p == 0) stop("no units with N_i > 0")
  if (identical(binding, "auto")) binding <- if (P$p <= 200) TRUE else structure(FALSE, auto = TRUE)
  f <- .ec_fit(P, targets, restrictions, integer, binding, parent.frame(), time_limit)
  if (isTRUE(attr(binding, "auto")) && length(f$specs))
    f$notes <- c(f$notes, "binding rows not computed (more than 200 units; set binding = TRUE)")
  out <- structure(c(list(call = match.call(), P = P, integer = integer, margins_integer = P$integer,
                          shares_input = !inp$integer_input), f), class = "eiCounts")
  if (!is.null(model)) out$cross <- .ec_cross_layer(out, model, allowance)
  out
}

# ---------------------------------------------------------------------------
# cross-layer block (R27)
# ---------------------------------------------------------------------------
.ec_utf8 <- function() isTRUE(l10n_info()[["UTF-8"]])
# the printed wording of contract C8 (Unicode, with an ASCII fallback in non-UTF-8 locales)
.ec_txt <- function(key) {
  u <- .ec_utf8()
  al <- if (u) "\u03b1" else "alpha"; mi <- if (u) "\u2212" else "-"; le <- if (u) "\u2264" else "<="
  switch(key,
         rej_id = paste0("disjoint: rejects (model + its rows) at level ", al),
         cov_bel = paste0("coverage at least 1 ", mi, " ", al, " ", mi, " P(belief rows false), which is unknown; ",
                          "level 1 ", mi, " ", al, " only if the belief rows are taken as certain"),
         rej_bel = paste0("disjoint: the model (+ its rows) or a belief row is wrong; P(reject while the belief rows are true) ",
                          le, " ", al),
         cov_id = "the intersection covers the actual value at the model interval's level (identity rows only)",
         nonint = paste("exact only for the stated margins, which are not integers; if they are estimates,",
                        "see estimated margins (Step 1b)"))
}

# status of the model interval (6a Sec. 4.3, with the rows added for L4, L5, L6 and M1(delta) on
# claude/mcount-theory, bbf8cf46; strings as given by the coordinating session, in the format of
# .rxc_ci_status). {what} names the interval used here; {thm} its theorem.
.ec_model_status <- function(fit) {
  st <- fit$setup
  msign <- length(st$extra$assoc) > 0 || length(st$extra$cr) > 0
  pooled <- identical(st$type, "pooled") || st$G == 1
  thm <- if (identical(fit$ci, "conservative")) "P3.5(a)" else "P3.5(b2)(ii)"
  what <- if (st$nz == 0) "normal interval" else if (identical(fit$ci, "clr")) "proj_*" else "ci_* (conservative)"
  r1 <- paste0(what, ": asymptotic guarantee only with few, large groups (P3 regime R1); with many small ",
               "groups a simulation finding")
  if (st$nz == 0) {
    pt <- "P3.5(a) with no free k: point-identified rung, one linear functional"
    if (fit$model == "L6" || st$G == 1) return(paste0(what, ": asymptotic guarantee (", pt, "), pooled model"))
    return(paste0(what, ": asymptotic guarantee (", pt, ") only with few, large groups (P3 regime R1); ",
                  "with many small groups a simulation finding"))
  }
  if (!pooled && msign) return(paste0(what, ": no coverage theorem for sign/CR restrictions under grouped rungs (R16 open)"))
  if (pooled && fit$model == "L4")
    return(paste0(what, ": asymptotic guarantee (", thm, "; L4 equalities as rows, as in PHASE5a R10), pooled model"))
  if (pooled && fit$model == "M1delta")
    return(paste0(what, ": asymptotic guarantee (", thm, "; M1(delta) with one group is M0), pooled model"))
  if (pooled) return(paste0(what, ": asymptotic guarantee (", thm, if (msign) "; with sign/CR rows PHASE5a R10, pooled M0" else "",
                            "), pooled model"))
  r1
}

.ec_check_model <- function(cnt, fit) {
  if (!inherits(fit, "eiRxC"))
    stop("'model' must be an eiBounds() fit. The count set is intersected only with an actual-target model ",
         "interval (with the residual allowance), never with a model estimate or an interval for the ",
         "model-implied rate (6a R27, check X3)")
  if (identical(fit$ci, "none"))
    stop("the model fit has ci = \"none\": a model estimate is not an actual-target confidence interval; ",
         "refit with ci = \"clr\" or \"conservative\" (6a R27, check X3)")
  st <- fit$setup; P <- cnt$P
  if (any(!st$keep)) stop("the cross-layer block needs the model fitted on all units (range = \"trim\" is not supported)")
  if (isTRUE(fit$eta > 0)) stop("the cross-layer block is not available with eta > 0")
  # align on the units with N_i > 0, in order (the x/t/n path of eiBounds keeps N_i = 0 units)
  pos <- which(st$N_all > 0)
  if (length(pos) != P$p || max(abs(st$X_all[pos, , drop = FALSE] - P$A / P$N)) > 1e-8 ||
      max(abs(st$N_all[pos] - P$N)) > 1e-8 * max(P$N))
    stop("the model fit and the count layer are not on the same units and margins")
  invisible(pos)
}

# actual-target model interval for one count target: model interval for the model-implied
# value (projection or conservative, at level a_est) plus the residual allowance (Lemma 4.1)
.ec_model_interval <- function(P, fit, tg, allowance, pos = seq_len(P$p)) {
  st <- fit$setup; alpha <- fit$alpha
  a_est <- alpha * 0.9; a_E <- alpha * 0.1          # as eiBounds with fpc = "popoviciu" (0.045 / 0.005)
  # representation fixed from the design alone: the one with the smallest S_phi (Remark 4.1a)
  rho_of <- function(rep) {
    v <- array(0, c(P$p, P$R, P$C))
    for (tm in rep$terms) v[tm$units, tm$rows, tm$cols] <- v[tm$units, tm$rows, tm$cols] + tm$wd
    rho <- numeric(P$p)
    for (r in seq_len(P$R)) rho <- rho + P$A[, r] * (apply(v[, r, , drop = FALSE], 1, max) - apply(v[, r, , drop = FALSE], 1, min))
    rho
  }
  rhos <- lapply(tg$reps, rho_of)
  S <- vapply(rhos, function(r) 0.5 * sqrt(sum(r^2)), 0)
  j <- which.min(S); rep <- tg$reps[[j]]; rho <- rhos[[j]]
  allow <- if (allowance == "popoviciu") stats::qnorm(1 - a_E / 2) * S[j] else sqrt(0.5 * log(2 / a_E) * sum(rho^2))
  lind <- if (sum(rho^2) > 0) max(rho^2) / sum(rho^2) else NA_real_
  # model-implied value as a linear functional of (beta, z)
  tgm <- list(cvec = numeric(st$nbeta), c0 = as.double(rep$const), gamma = numeric(st$nz))
  for (tm in rep$terms) {
    a <- matrix(0, st$R, st$C); a[tm$rows, tm$cols] <- 1
    W <- matrix(0, nrow(st$X), st$R)
    W[pos[tm$units], tm$rows] <- tm$wd * (st$N * st$X)[pos[tm$units], tm$rows]
    t1 <- .rxc_target(st, a, w = W)
    tgm$cvec <- tgm$cvec + t1$cvec; tgm$c0 <- tgm$c0 + t1$c0; tgm$gamma <- tgm$gamma + t1$gamma
  }
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()))
  cons <- .rxc_constraints(st)
  ivs <- list()
  if (st$nz == 0) {
    pl <- sum(tgm$cvec * st$beta) + tgm$c0
    se <- sqrt(max(0, sum(tgm$cvec * (st$V %*% tgm$cvec))))
    ivs[[1]] <- pl + c(-1, 1) * stats::qnorm(1 - a_est / 2) * se
    type <- "normal interval (point-identified rung)"
  } else {
    Zs <- if (!is.null(fit$zdraws)) fit$zdraws else .rxc_zdraws(st$V, fit$settings$nz_draw)
    if (identical(fit$ci, "clr")) {
      for (sy in .rxc_systems(cons, "branch"))
        ivs[[length(ivs) + 1]] <- attr(.rxc_ci_clr(st, sy, tgm, a_est / 2, a_est / 2, Zs, fit$settings$ndraw,
                                                  proj_alpha = a_est), "projection")
      type <- "projection interval (P3.5(b2)(ii))"
    } else {
      for (sy in .rxc_systems(cons, "slab")) ivs[[length(ivs) + 1]] <- .rxc_ci_conservative(st, sy, tgm, a_est, Zs)$ci
      type <- "conservative interval (P3.5(a))"
    }
  }
  ivs <- Filter(function(v) !any(is.na(v)) && v[1] <= v[2], ivs)
  ivs <- lapply(ivs, function(v) v + c(-1, 1) * allow)
  list(pieces = if (length(ivs)) .rxc_union(ivs, tol = 0) else NULL, type = type, S_phi = S[j],
       allowance = allow, allowance_type = if (allowance == "popoviciu") "Popoviciu S_phi (6a Lemma 4.1(a))" else
         "Hoeffding t_E(phi) (6a Lemma 4.1(b))", representation = j, lindeberg = lind)
}

.ec_intersect <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NULL)
  out <- list()
  for (i in seq_len(nrow(A))) for (j in seq_len(nrow(B))) {
    lo <- max(A[i, 1], B[j, 1]); hi <- min(A[i, 2], B[j, 2])
    if (lo <= hi) out[[length(out) + 1]] <- c(lo, hi)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}

.ec_cross_layer <- function(cnt, fit, allowance) {
  pos <- .ec_check_model(cnt, fit)
  P <- cnt$P
  belief <- !is.null(cnt$restrictions) && any(cnt$restrictions$role == "belief")
  rows <- list()
  for (j in seq_along(cnt$targets)) {
    tg <- cnt$targets[[j]]; tb <- cnt$table[j, ]
    base <- data.frame(target = tg$label, count_lower = tb$lower, count_upper = tb$upper,
                       model_lower = NA_real_, model_upper = NA_real_, lower = NA_real_, upper = NA_real_,
                       n_pieces = NA_integer_, disjoint = NA, disjoint_margins_only = NA, model_type = NA_character_,
                       status = NA_character_, allowance = NA_character_, S_phi = NA_real_, lindeberg = NA_real_,
                       note = "", stringsAsFactors = FALSE)
    if (!tg$design) {
      base$note <- if (identical(tg$design_reason, "outcome_set"))
        "count layer only: the unit set is selected using the outcomes, so the coefficients are not functions of the design (6a Lemma 4.1, C2)" else
        "count layer only: the coefficients use outcome margins (6a Lemma 4.1 does not apply)"
      rows[[j]] <- base; next
    }
    if (isTRUE(tg$given_set))
      base$note <- paste("the unit set must be fixed from the design (6a Lemma 4.1); a set chosen using",
                         "the outcomes voids the cross-layer interval")
    mi <- .ec_model_interval(P, fit, tg, allowance, pos)
    base$model_type <- mi$type; base$status <- .ec_model_status(fit); base$allowance <- mi$allowance_type
    base$S_phi <- mi$S_phi; base$lindeberg <- mi$lindeberg
    if (is.null(mi$pieces)) {
      base$note <- "the model interval is empty: the model layer rejects the model (+ its rows) by itself"
      base$disjoint <- TRUE; base$disjoint_margins_only <- TRUE
      rows[[j]] <- base; next
    }
    base$model_lower <- min(mi$pieces[, 1]); base$model_upper <- max(mi$pieces[, 2])
    Cs <- cnt$solutions[[j]]$pieces
    I <- .ec_intersect(Cs, mi$pieces)
    Idd <- .ec_intersect(matrix(c(tb$dd_lower, tb$dd_upper), 1), mi$pieces)
    base$disjoint <- is.null(I) && !is.null(Cs)
    base$disjoint_margins_only <- is.null(Idd)
    if (!is.null(I)) { base$lower <- min(I[, 1]); base$upper <- max(I[, 2]); base$n_pieces <- nrow(I) }
    else base$n_pieces <- 0L
    if (is.null(Cs)) base$note <- "the count set is empty"
    rows[[j]] <- base
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  list(table = tab, belief = belief, alpha = fit$alpha, model = fit$model,
       model_rows = if (is.null(fit$modules)) character(0) else fit$modules$restriction,
       identity_rows = c("margins", if (!is.null(cnt$restrictions)) cnt$restrictions$restriction[cnt$restrictions$role == "identity"]),
       belief_rows = if (!is.null(cnt$restrictions)) cnt$restrictions$restriction[cnt$restrictions$role == "belief"] else character(0))
}

# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------
.ec_fmt <- function(l, u, digits) ifelse(is.na(l), "      --      ", sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), l, u))

.ec_print_block <- function(x, digits = 3, header = TRUE) {
  P <- x$P; tb <- x$table
  if (header) {
    cat("Count layer (realised table): exact bounds on actual counts and rates\n")
    cat(sprintf("Units: %d%s; margins: %s; integer mode: %s\n", P$p,
                if (P$n_dropped) sprintf(" (%d with N_i = 0 dropped)", P$n_dropped) else "",
                if (P$integer) "integer counts" else if (isTRUE(x$shares_input)) "given as shares (not counts)" else "not integers",
                x$integer))
  }
  out <- data.frame(target = tb$target, margins_only = .ec_fmt(tb$dd_lower, tb$dd_upper, digits),
                    stringsAsFactors = FALSE)
  has_rows <- !is.null(x$restrictions)
  if (has_rows) {
    out$with_rows <- .ec_fmt(tb$lower, tb$upper, digits)
    for (j in seq_len(nrow(tb))) {
      pc <- x$solutions[[j]]$pieces
      if (!is.null(pc) && nrow(pc) > 1) out$with_rows[j] <- paste(.ec_fmt(pc[, 1], pc[, 2], digits), collapse = " U ")
    }
  }
  out$width <- ifelse(is.na(tb$lower), "--", sprintf(paste0("%.", digits, "f"), tb$width))
  out$integer_status <- tb$integer_status
  cat("\nTarget (denominator)  |  margins only = Duncan-Davis (exact LP)", if (has_rows) "  |  with the count rows", "\n", sep = "")
  print(out, row.names = FALSE, right = FALSE)
  for (j in seq_len(nrow(tb))) cat(sprintf("  %s: denominator %s\n", tb$target[j], tb$denominator[j]))
  if (has_rows) {
    cat("\nCount rows (object = \"actual\": beliefs or facts about the realised table):\n")
    rt <- x$restrictions
    for (k in seq_len(nrow(rt)))
      cat(sprintf("  %-14s %s [%s row]%s\n", rt$id[k], rt$restriction[k], rt$role[k],
                  if (rt$experimental[k]) " EXPERIMENTAL (2c M2; review batch B comment 1 open)" else ""))
    if (any(!is.na(tb$binding_lower))) {
      cat("Binding rows (drop-one face test: removing the row moves the endpoint; rows that bind only together are not shown):\n")
      for (j in seq_len(nrow(tb)))
        cat(sprintf("  %s: lower %s; upper %s\n", tb$target[j],
                    if (nzchar(tb$binding_lower[j])) tb$binding_lower[j] else "none",
                    if (nzchar(tb$binding_upper[j])) tb$binding_upper[j] else "none"))
    }
    for (j in seq_len(nrow(tb))) {
      br <- x$solutions[[j]]$branches
      if (is.null(br)) next
      cat(sprintf("  CR branches for %s: %s; the union is %s\n", tb$target[j],
                  paste(sprintf("%s %s", br$branch, ifelse(is.na(br$lower), if (P$integer && x$integer != "lp") "empty (integer)" else "empty",
                                                           .ec_fmt(br$lower, br$upper, digits))), collapse = "; "),
                  if (tb$n_pieces[j] == 0) "empty" else if (tb$n_pieces[j] == 1) "one interval" else "two intervals"))
    }
    if (any(vapply(x$specs, function(s) s$type == "odds" && s$psi[1] <= 1 && s$psi[2] >= 1, TRUE)))
      cat("  the odds-ratio dial contains 1: untestable (with continuous tables and no other rows the margins alone\n",
          "  do not contradict it: the equal-rates table has odds ratio 1)\n", sep = "")
  }
  cat("\n")
  if (any(grepl("TU class", tb$integer_status)))
    cat("  integer-exact (TU class): rows rounded inward in exact rational arithmetic, LP over the rounded rows (R19(b));",
        "\n  the endpoints and witnesses are whole-number tables\n")
  if (any(grepl("MILP", tb$integer_status)))
    cat("  integer-exact (MILP): rows outside the TU classes of R19(b); each endpoint solved as an integer program (R19(c))\n")
  if (any(tb$integer_status == "continuous only")) cat("  continuous only: sharp for continuous tables (R19(a))\n")
  if (any(grepl("integer-infeasible", tb$integer_status)))
    cat("  integer-infeasible: continuous tables satisfy the rows, but no whole-number table does (R19(c))\n")
  if (!P$integer) cat("  Note:", .ec_txt("nonint"), "\n")
  if (length(P$log$rat))
    cat(sprintf("  Note: values given as doubles were rationalised (tolerance %g, denominators <= 1e6): %s\n",
                .ec_rat_tol, paste(utils::head(unique(P$log$rat), 5), collapse = "; ")))
  for (nt in unique(P$log$notes)) cat("  Note:", nt, "\n")
  if (any(grepl("unresolved", tb$integer_status)))
    cat("  integer status unresolved: the MILP hit the time limit (or no MILP solver is installed); the range shown is\n",
        "  the LP over the rounded rows, an outer bound for whole-number tables (R19(c) item 1); never integer-exact\n", sep = "")
  for (nt in x$notes) cat("  Note:", nt, "\n")
  if (any(!is.na(tb$lp_lower) & !is.na(tb$lower) & (abs(tb$lp_lower - tb$lower) > 1e-12 | abs(tb$lp_upper - tb$upper) > 1e-12)))
    cat("  Continuous LP range (lp_lower/lp_upper in $table) differs from the integer range for some targets\n")
}

.ec_print_cross <- function(cr, digits = 3) {
  if (is.null(cr)) return(invisible())
  tb <- cr$table
  cat("\nCross-layer block (R27): count set intersected with an actual-target model interval\n")
  cat(sprintf("Model: %s; level %g%%. Count rows and model rows are never put into one LP.\n",
              cr$model, 100 * (1 - cr$alpha)))
  out <- data.frame(target = tb$target, count_set = .ec_fmt(tb$count_lower, tb$count_upper, digits),
                    model_interval = .ec_fmt(tb$model_lower, tb$model_upper, digits),
                    intersection = ifelse(tb$disjoint %in% TRUE, "disjoint", .ec_fmt(tb$lower, tb$upper, digits)),
                    stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  for (j in seq_len(nrow(tb))) {
    if (nzchar(tb$note[j])) { cat(sprintf("  %s: %s\n", tb$target[j], tb$note[j])); if (is.na(tb$model_type[j])) next }
    cat(sprintf("  %s: model interval = %s + %s; status %s\n", tb$target[j], tb$model_type[j], tb$allowance[j], tb$status[j]))
  }
  cat("  Identity rows (hold in every replication):", paste(cr$identity_rows, collapse = "; "), "\n")
  if (length(cr$belief_rows)) cat("  Belief rows (statements about the realised table):", paste(cr$belief_rows, collapse = "; "), "\n")
  if (!cr$belief) {
    cat("  ", .ec_txt("cov_id"), "\n", sep = "")
    if (any(tb$disjoint %in% TRUE)) cat("  ", .ec_txt("rej_id"), sprintf(" (alpha = %g): %s", cr$alpha,
                                                                    paste(tb$target[tb$disjoint %in% TRUE], collapse = ", ")), "\n", sep = "")
  } else {
    cat("  ", .ec_txt("cov_bel"), "\n", sep = "")
    if (any(tb$disjoint %in% TRUE)) {
      cat("  ", .ec_txt("rej_bel"), sprintf(" (alpha = %g): %s", cr$alpha, paste(tb$target[tb$disjoint %in% TRUE], collapse = ", ")), "\n", sep = "")
      dm <- tb$disjoint_margins_only %in% TRUE
      if (any(dm)) cat("  ", .ec_txt("rej_id"), sprintf(" with the margins-only (identity) count set too: %s", paste(tb$target[dm], collapse = ", ")), "\n", sep = "")
    }
  }
  hoef <- grepl("Hoeffding", tb$allowance)
  if (any(hoef)) cat("  Hoeffding allowance: the model part is asymptotic; the allowance holds in finite samples\n")
  big <- !hoef & !is.na(tb$lindeberg) & tb$lindeberg > 0.1
  if (any(big))
    cat(sprintf("  Note (heuristic): Lemma 4.1(a) needs no unit to dominate; max rho_i^2 / sum rho_i^2 = %s (%s)\n",
                paste(sprintf("%.3f", tb$lindeberg[big]), collapse = ", "), paste(tb$target[big], collapse = ", ")))
  cat("  Realised beliefs (object = \"actual\", count layer):", if (length(cr$belief_rows) + length(cr$identity_rows) > 1)
    paste(c(cr$identity_rows[-1], cr$belief_rows), collapse = "; ") else "none", "\n")
  cat("  Model beliefs (object = \"model\", model layer):", if (length(cr$model_rows)) paste(cr$model_rows, collapse = "; ") else "none", "\n")
}

#' @export
print.eiCounts <- function(x, digits = 3, ...) {
  .ec_print_block(x, digits)
  .ec_print_cross(x[["cross"]], digits)
  cat("  witness(fit) gives a table attaining each endpoint (a certificate of attainability, not an estimate)\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# witnesses
# ---------------------------------------------------------------------------

#' Tables attaining each endpoint
#'
#' \code{witness()} returns, for each target of an \code{\link{eiCounts}} fit (or the count
#' block of an \code{eiBounds()} fit) and each endpoint, a unit-by-row-by-column table that
#' reproduces every margin, satisfies every count row and attains the endpoint (6a R19(d)):
#' a basic optimal LP solution, or the MILP solution; a whole-number table when the integer
#' status is integer-exact. Each is checked against the margins and every row (exactly, in
#' rational arithmetic, for integer tables). A witness certifies attainability; it is one
#' of many tables with the same value and is \emph{not} an estimate of the realised table.
#'
#' @param object an \code{eiCounts} fit, or an \code{eiBounds()} fit with a count block.
#' @param target optional target label(s) or index(es).
#' @param ... unused.
#' @return an object of class \code{"ecWitness"}: a list with one element per target and
#'   endpoint, each with \code{table} (array units x rows x columns), \code{value},
#'   \code{margins_ok}, \code{rows_ok} and \code{integer}.
#' @export
witness <- function(object, target = NULL, ...) UseMethod("witness")

#' @export
witness.eiRxC <- function(object, target = NULL, ...) {
  if (is.null(object$counts)) stop("this fit has no count block")
  witness(object$counts, target = target, ...)
}

#' @export
witness.eiCounts <- function(object, target = NULL, ...) {
  P <- object$P; tb <- object$table
  js <- if (is.null(target)) seq_len(nrow(tb)) else if (is.numeric(target)) target else match(target, tb$target)
  if (any(is.na(js))) stop("unknown target")
  out <- list()
  for (j in js) {
    tg <- object$targets[[j]]; so <- object$solutions[[j]]
    for (e in 1:2) {
      x <- so$x[[e]]; sys <- if (is.null(so$sys)) "" else so$sys[e]
      int <- P$integer && object$integer != "lp"
      if (is.null(x) && !is.na(tb$lower[j])) {
        # margins-only targets computed in closed form: solve now (separable, unit by unit)
        s <- .ec_optimise(P, list(), tg$obj, c(-1, 1)[e], rd = int, int = FALSE)
        x <- if (int) round(s$x) else s$x
      }
      if (is.null(x)) next
      rows <- object$built$common
      if (nzchar(sys)) for (b in strsplit(sys, ", ")[[1]]) {
        id <- sub("[+-]$", "", b); br <- if (grepl("\\+$", b)) "plus" else "minus"
        u <- Filter(function(u) u$id == id, object$built$unions)[[1]]
        rows <- c(rows, u[[br]])
      }
      exact <- int && all(x == round(x))
      arr <- aperm(array(x, c(P$C, P$R, P$p), dimnames = list(P$cols, P$rows, P$units)), c(3, 2, 1))
      out[[length(out) + 1]] <- list(
        target = tb$target[j], end = c("lower", "upper")[e], value = .ec_value(P, tg, x),
        reported = c(tb$lower[j], tb$upper[j])[e], branch = sys, table = arr, integer = exact,
        margins_ok = .ec_margins_ok(P, x, exact),
        rows_ok = all(vapply(rows, function(rw) .ec_row_ok(rw, x, P, exact), TRUE)))
    }
  }
  structure(out, class = "ecWitness")
}

#' @export
print.ecWitness <- function(x, max_units = 10, ...) {
  cat("Witness tables: each reproduces every margin, satisfies every count row and attains the\n",
      "endpoint. A certificate of attainability, not an estimate of the realised table.\n", sep = "")
  for (w in x) {
    cat(sprintf("\n%s, %s end = %.6g%s: margins %s, rows %s%s\n", w$target, w$end, w$value,
                if (nzchar(w$branch)) paste0(" (branch ", w$branch, ")") else "",
                if (w$margins_ok) "reproduced" else "NOT reproduced", if (w$rows_ok) "all satisfied" else "VIOLATED",
                if (w$integer) "; whole numbers" else ""))
    arr <- w$table; p <- dim(arr)[1]
    for (i in seq_len(min(p, max_units))) {
      cat("  unit", dimnames(arr)[[1]][i], "\n")
      m <- arr[i, , , drop = TRUE]; if (!is.matrix(m)) m <- matrix(m, dim(arr)[2], dim(arr)[3], dimnames = dimnames(arr)[2:3])
      print(round(m, 6))
    }
    if (p > max_units) cat(sprintf("  ... %d more units\n", p - max_units))
  }
  invisible(x)
}
