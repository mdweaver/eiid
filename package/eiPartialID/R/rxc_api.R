# User-facing functions for R x C bounds.

# ---------------------------------------------------------------------------
# Input parsing: formula, data + column names, or shares (x, t, n)
# ---------------------------------------------------------------------------
.rxc_parse <- function(formula, data, rows, cols, x, t, n, group, mismatch, tol = 1e-6) {
  grp <- NULL
  if (!is.null(group)) {
    if (is.character(group) && length(group) == 1 && !is.null(data) && group %in% names(data)) {
      grp <- data[[group]]
    } else grp <- group
  }
  if (!is.null(formula)) {
    if (is.null(data)) stop("'data' is required with a formula")
    cols <- all.vars(formula[[2]]); rows <- all.vars(formula[[3]])
  }
  if (!is.null(rows) && !is.null(cols)) {
    if (is.null(data)) stop("'data' is required with 'rows' and 'cols'")
    Rc <- as.matrix(data[, rows, drop = FALSE]); Cc <- as.matrix(data[, cols, drop = FALSE])
    storage.mode(Rc) <- "double"; storage.mode(Cc) <- "double"
    if (any(Rc < 0, na.rm = TRUE) || any(Cc < 0, na.rm = TRUE)) stop("counts must be nonnegative")
    nr <- rowSums(Rc); nc <- rowSums(Cc)
    bad <- abs(nr - nc) > tol * pmax(nr, nc, 1)
    if (any(bad)) {
      msg <- sprintf("row and column totals differ in %d of %d units (same population in both margins is required)",
                     sum(bad), length(bad))
      if (mismatch == "error") stop(msg, "; use mismatch = \"rescale\" to rescale each margin to shares")
      warning(msg, "; each margin rescaled to shares")
    }
    ok <- nr > 0 & nc > 0
    if (any(!ok)) warning(sum(!ok), " units with zero totals dropped")
    X <- Rc[ok, , drop = FALSE] / nr[ok]; T <- Cc[ok, , drop = FALSE] / nc[ok]; N <- nr[ok]
    if (!is.null(grp)) grp <- grp[ok]
    colnames(X) <- rows; colnames(T) <- cols
  } else if (!is.null(x) && !is.null(t) && !is.null(n)) {
    if (is.null(dim(x))) { X <- cbind(x, 1 - x); colnames(X) <- c("X", "notX") } else X <- as.matrix(x)
    if (is.null(dim(t))) { T <- cbind(t, 1 - t); colnames(T) <- c("T", "notT") } else T <- as.matrix(t)
    if (is.null(colnames(X))) colnames(X) <- paste0("row", seq_len(ncol(X)))
    if (is.null(colnames(T))) colnames(T) <- paste0("col", seq_len(ncol(T)))
    N <- as.numeric(n)
    if (any(abs(rowSums(X) - 1) > 1e-6) || any(abs(rowSums(T) - 1) > 1e-6))
      stop("rows of 'x' and 't' must be shares summing to one")
  } else stop("supply a formula and data, 'data' with 'rows' and 'cols', or 'x', 't' and 'n'")
  if (any(!is.finite(X)) || any(!is.finite(T)) || any(!is.finite(N))) stop("missing or non-finite values")
  if (!is.null(grp) && length(grp) != nrow(X)) stop("'group' has the wrong length")
  list(X = X, T = T, N = N, group = grp)
}

.rxc_index <- function(v, labels) {
  if (is.numeric(v)) return(as.integer(v))
  i <- match(v, labels)
  if (any(is.na(i))) stop("unknown label(s): ", paste(v[is.na(i)], collapse = ", "))
  i
}

# targets: NULL (all cells), or a named list of c(row, col) pairs or R x C weight matrices
.rxc_targets <- function(targets, rl, cl) {
  R <- length(rl); C <- length(cl)
  if (is.null(targets)) {
    out <- list()
    for (r in seq_len(R)) for (c in seq_len(C)) {
      a <- matrix(0, R, C); a[r, c] <- 1
      out[[paste0(rl[r], " -> ", cl[c])]] <- a
    }
    return(out)
  }
  if (is.null(names(targets)) || any(names(targets) == "")) stop("'targets' must be a named list")
  lapply(targets, function(tg) {
    if (is.matrix(tg)) {
      if (!all(dim(tg) == c(R, C))) stop("target matrices must be ", R, " x ", C)
      return(tg)
    }
    a <- matrix(0, R, C)
    a[.rxc_index(tg[1], rl), .rxc_index(tg[2], cl)] <- 1
    a
  })
}

.rxc_restrictions <- function(res, rl, cl, groups) {
  if (is.null(res)) return(NULL)
  ri <- function(v) .rxc_index(v, rl); ci <- function(v) .rxc_index(v, cl)
  gi <- function(v) if (is.null(v)) NULL else .rxc_index(as.character(v), groups)
  out <- list()
  if (!is.null(res$zeros)) out$zeros <- lapply(res$zeros, function(z)
    list(row = ri(z[["row"]]), col = ci(z[["col"]]),
         group = if ("group" %in% names(z)) gi(z[["group"]]) else NULL))
  if (!is.null(res$order)) out$order <- lapply(res$order, function(o)
    list(high = ri(o[["high"]]), low = ri(o[["low"]]), col = ci(o[["col"]])))
  if (!is.null(res$sign)) out$sign <- lapply(res$sign, function(s)
    list(row = ri(s[["row"]]), col = ci(s[["col"]]), from = ri(s[["from"]]), to = ri(s[["to"]]),
         sign = sign(as.numeric(s[["sign"]]))))
  if (!is.null(res$slope)) out$slope <- as.numeric(res$slope)
  if (!is.null(res$linear)) out$linear <- lapply(res$linear, function(l) {
    tm <- l$terms
    m <- cbind(row = ri(tm$row), col = ci(tm$col), context = ri(tm$context), value = as.numeric(tm$value))
    list(terms = m, rhs = as.numeric(l$rhs), group = gi(l$group))
  })
  if (!is.null(res$external)) out$external <- lapply(res$external, function(e) {
    a <- if (is.matrix(e$target)) e$target else {
      a <- matrix(0, length(rl), length(cl)); a[ri(e$target[1]), ci(e$target[2])] <- 1; a }
    list(weights = a, lower = e$lower, upper = e$upper)
  })
  R <- length(rl)
  ver <- function(v, allowed, what) {
    v <- if (is.null(v)) "W1" else as.character(v)
    if (!v %in% allowed) stop(what, ": version must be one of ", paste(allowed, collapse = ", "))
    v
  }
  oth <- function(o, r, v, what) {
    if (v != "pair") return(NULL)
    if (is.null(o)) { if (R == 2) return(setdiff(1:2, r)) else stop(what, ": version \"pair\" needs 'other'") }
    o <- ri(o); if (o == r) stop(what, ": 'other' must differ from 'row'")
    o
  }
  if (!is.null(res$assoc)) out$assoc <- lapply(res$assoc, function(s) {
    type <- as.character(s[["type"]])
    if (length(type) != 1 || !type %in% c("dB", "dW")) stop("assoc: type must be \"dB\" or \"dW\"")
    v <- ver(s[["version"]], if (type == "dB") c("W1", "pair") else c("W1", "W2"),
             paste0("assoc (", type, ")"))
    r <- ri(s[["row"]])
    sg <- sign(as.numeric(if (is.null(s[["sign"]])) 1 else s[["sign"]]))
    if (sg == 0) stop("assoc: sign must be +1 or -1")
    list(type = type, version = v, row = r, col = ci(if (is.null(s[["col"]])) 1 else s[["col"]]),
         sign = sg, other = oth(s[["other"]], r, v, "assoc"))
  })
  if (!is.null(res$cr)) out$cr <- lapply(res$cr, function(s) {
    v <- ver(s[["version"]], c("W1", "W2", "pair"), "cr")
    r <- ri(s[["row"]])
    list(row = r, col = ci(if (is.null(s[["col"]])) 1 else s[["col"]]), version = v,
         other = oth(s[["other"]], r, v, "cr"))
  })
  for (nm in setdiff(names(res), c("zeros", "order", "sign", "slope", "linear", "external", "assoc", "cr")))
    warning("unknown restriction type '", nm, "' ignored")
  out
}

#' Bounds for R x C ecological inference
#'
#' \code{eiBounds()} computes partial-identification bounds and confidence
#' intervals for group-by-outcome rates in an R x C ecological inference
#' problem, under a chosen rung of the assumption ladder. It extends the 2 x 2
#' method of Jiang et al. (2019) to R x C tables, to grouped units (e.g. polls
#' within ridings), and to optional additional restrictions.
#' See \code{docs/extension/PLAIN_LANGUAGE_GUIDE.md} in the source repository.
#'
#' Data can be given in three ways: a formula \code{cbind(c1, c2, ...) ~
#' cbind(r1, r2, ...)} with counts in \code{data}; \code{data} with the names
#' of the row-count and column-count variables in \code{rows} and \code{cols};
#' or shares \code{x} (vector for two rows, or matrix), \code{t} (vector or
#' matrix) and unit sizes \code{n}.
#'
#' Models (rungs): \code{"M0"} pooled linear contextual effects;
#' \code{"M2"} group-specific linear contextual effects; \code{"M1delta"}
#' shared contextual effects with group deviations up to \code{delta};
#' \code{"M1"} shared contextual effects with group baselines (group dummies);
#' \code{"L4"} as M1 (or M0 without groups) with the row
#' \code{row_constant} having no contextual effect; \code{"L5"} constancy within
#' groups (Goodman within each group); \code{"L6"} Goodman regression.
#'
#' Reported for each target: \code{estimate} (the estimated identified set for
#' the model-implied rate; not a confidence interval), \code{ci_model} (CI for
#' the model-implied rate at the observed compositions), \code{ci_actual} (CI
#' for the actual rate in these units, including an allowance for unit-level
#' deviations bounded by Bhatia-Davis or Popoviciu) and the Duncan-Davis bounds.
#'
#' @param formula optional formula \code{cbind(column counts) ~ cbind(row counts)}.
#' @param data data frame with counts (and the grouping variable).
#' @param rows,cols names of the row-count and column-count variables.
#' @param x,t,n shares and unit sizes (alternative input).
#' @param group grouping variable (name in \code{data} or a vector).
#' @param model rung of the assumption ladder (see Details).
#' @param delta deviation bound for \code{"M1delta"} (probability units).
#' @param row_constant row (label or index) without contextual effect, for \code{"L4"}.
#' @param targets \code{NULL} for all cells, or a named list of \code{c(row, col)}
#'   pairs or R x C weight matrices (linear combinations of cell rates).
#' @param range composition range: \code{"hull"} (observed units, default),
#'   \code{"trim"} (drop extreme units; their contribution is bounded by
#'   Duncan-Davis) or \code{"extrapolate"} (impose the model beyond the data).
#' @param trim list with \code{quantile}, \code{min_n} and/or \code{drop} (indices).
#' @param extrapolate \code{"simplex"} (pure units) or a matrix of compositions.
#' @param weights regression weights: \code{"none"} (as Jiang et al.) or \code{"n"}.
#' @param vcov \code{"HC3"} (default; falls back to HC1, HC0, classical), \code{"HC1"} or \code{"HC0"}.
#' @param cluster optional cluster variable (pooled models), giving CR1 variances.
#' @param restrictions optional list with elements \code{zeros}, \code{order},
#'   \code{sign}, \code{slope}, \code{linear}, \code{external}, \code{assoc}
#'   and \code{cr}; see the vignette. \code{assoc} is a list of signs on mean
#'   associations, each \code{list(type = "dB" or "dW", row, col, sign = +1 or -1,
#'   version, other)}: \code{"dB"} is the between-group association
#'   \eqn{\delta_B} (version \code{"W1"}, row vs the rest, default; or
#'   \code{"pair"}, row vs \code{other}), \code{"dW"} the within-group
#'   (contextual) association \eqn{\delta_W} (version \code{"W1"}, one-vs-rest,
#'   default; or \code{"W2"}, \eqn{\Delta_r}). \code{cr} is a list of
#'   contextual-reinforcement restrictions (\eqn{\delta_B} and \eqn{\delta_W}
#'   have the same sign, which is not assumed known), each \code{list(row, col,
#'   version = "W1", "W2" or "pair", other)}. W1 (and every 2 x 2 CR) is one LP
#'   (a slab); W2 and pair are a union of two branches and may give two disjoint
#'   intervals. All associations are model-implied (LCE) means with design
#'   weights \eqn{N_i}. Signs on \eqn{\delta_W} can be refuted by the data;
#'   with two rows a sign on \eqn{\delta_B} and CR cannot be refuted by
#'   aggregates on their own (refutation is always joint with the rung).
#' @param alpha 1 - confidence level.
#' @param ci \code{"clr"} (default), \code{"conservative"} or \code{"none"}.
#' @param fpc unit-level allowance for \code{ci_actual}: \code{"bhatia-davis"} (default) or \code{"popoviciu"}.
#' @param jkst also report Jiang et al. CI_x intervals.
#' @param jkst_x values of x for CI_x.
#' @param eta sensitivity: allowed undetectable nonlinearity (probability units).
#' @param ndraw draws used to collect candidate endpoints.
#' @param nz_draw Gaussian draws for critical values.
#' @param seed optional random seed (for reproducible CIs).
#' @param tests run the emptiness test and (for M1) the test of its testable part.
#' @param range_offset shrink the observed range by this amount (two-row tables
#'   only; Jiang et al. use 1e-5).
#' @param mismatch \code{"error"} or \code{"rescale"} when row and column totals differ.
#' @param union with CR restrictions the estimate and CIs can be two disjoint
#'   intervals: \code{"pieces"} (default) reports the hull in
#'   \code{*_lower}/\code{*_upper} and the excluded gap in \code{*_gap_lower}/
#'   \code{*_gap_upper}; \code{"hull"} reports the hull only (valid, wider or equal).
#' @param diagnostics compute \code{eiRedundancy()} diagnostics and store them in
#'   \code{$sign} (default: when \code{assoc} or \code{cr} restrictions are given).
#' @return an object of class \code{"eiRxC"}.
#' @references Jiang, W., King, G., Schmaltz, A. and Tanner, M. A. (2019).
#'   Ecological Regression with Partial Identification. Political Analysis.
#' @examples
#' library(eco); data(census)
#' fit <- eiBounds(x = census$X, t = census$Y, n = census$N, model = "M0", seed = 1)
#' fit
#' @export
eiBounds <- function(formula = NULL, data = NULL, rows = NULL, cols = NULL,
                     x = NULL, t = NULL, n = NULL, group = NULL,
                     model = c("M0", "M1", "M1delta", "M2", "L4", "L5", "L6"),
                     delta = 0, row_constant = NULL, targets = NULL,
                     range = c("hull", "trim", "extrapolate"), trim = NULL, extrapolate = NULL,
                     weights = c("none", "n"), vcov = c("HC3", "HC1", "HC0"), cluster = NULL,
                     restrictions = NULL, alpha = 0.05, ci = c("clr", "conservative", "none"),
                     fpc = c("bhatia-davis", "popoviciu"), jkst = FALSE,
                     jkst_x = seq(0, 2, by = 0.25), eta = 0, ndraw = 100, nz_draw = 4000,
                     seed = NULL, tests = TRUE, range_offset = 0,
                     mismatch = c("error", "rescale"), union = c("pieces", "hull"),
                     diagnostics = NULL) {
  cl <- match.call()
  union <- match.arg(union)
  model <- match.arg(model); range <- match.arg(range); weights <- match.arg(weights)
  vcov <- match.arg(vcov); ci <- match.arg(ci); fpc <- match.arg(fpc); mismatch <- match.arg(mismatch)
  if (!is.null(seed)) set.seed(seed)
  dat <- .rxc_parse(formula, data, rows, cols, x, t, n, group, mismatch)
  rl <- colnames(dat$X); cll <- colnames(dat$T)
  if (!is.null(cluster) && is.character(cluster) && length(cluster) == 1 && !is.null(data))
    cluster <- data[[cluster]]
  grp_levels <- if (is.null(dat$group)) "all" else levels(as.factor(dat$group))
  L4row <- NULL
  if (model == "L4") {
    if (is.null(row_constant)) stop("model \"L4\" needs 'row_constant'")
    L4row <- .rxc_index(row_constant, rl)
  }
  extra <- .rxc_restrictions(restrictions, rl, cll, grp_levels)
  st <- .rxc_setup(dat$X, dat$T, dat$N, dat$group, model = model, delta = delta, L4row = L4row,
                   range = range, trim = trim, extrapolate = extrapolate,
                   range_offset = range_offset, weights = weights, vcov = vcov,
                   cluster = cluster, extra = extra)
  tl <- .rxc_targets(targets, rl, cll)
  est <- .rxc_estimate(st, tl, alpha = alpha, ci = ci, fpc = fpc, jkst = jkst, jkst_x = jkst_x,
                       eta = eta, ndraw = ndraw, nz_draw = nz_draw, union = union)
  tst <- NULL
  if (tests) {
    ed <- .rxc_empty_test(st, alpha = alpha, nz_draw = nz_draw, detail = TRUE)
    tst <- list(empty_rejected = if (is.list(ed)) ed$rejected else ed,
                empty_detail = if (is.list(ed)) ed else NULL)
    if (model %in% c("M1", "L4") && st$G > 1)
      tst$m1 <- .rxc_m1_test(dat$X[st$keep, , drop = FALSE], dat$T[st$keep, , drop = FALSE],
                             dat$N[st$keep], dat$group[st$keep], weights = weights, vcov = vcov)
  }
  if (range == "extrapolate")
    warning("range = \"extrapolate\": the linear model is assumed where no units are observed; this cannot be checked")
  notes <- character(0)
  mlabels <- .rxc_module_labels(extra, rl, cll)
  combo <- if (length(mlabels)) paste0("the rung (", model, ") + ", paste(mlabels, collapse = " + ")) else
    paste0("the rung (", model, ")")
  if (any(est$table$empty))
    notes <- c(notes, paste0("empty bound: ", combo, " is contradicted by the data (or by sampling noise; ",
                             "see the emptiness test). Refutation is joint: the combination is rejected, ",
                             "never a restriction on its own"))
  if (isTRUE(tst$empty_rejected)) {
    ra <- tst$empty_detail$rung_alone_rejected
    notes <- c(notes, paste0("emptiness test rejects ", combo,
                             if (isTRUE(ra)) "; the rung alone is also rejected" else
                               if (identical(ra, FALSE)) "; the rung alone is not rejected" else ""))
  }
  sign_like <- c(vapply(extra$assoc, function(a) a$type == "dB", TRUE), rep(TRUE, length(extra$cr)))
  if (length(sign_like) && all(sign_like) && length(rl) == 2 &&
      (any(est$table$empty) || isTRUE(tst$empty_rejected)))
    notes <- c(notes, "with two rows a sign on delta_B and CR can never be refuted by aggregate data on their own; the rejection concerns their combination with the linear-contextual-effects rung")
  if (length(extra$assoc) || length(extra$cr)) {
    if (any(est$table$ci_rejected) || any(est$table$empty) || isTRUE(tst$empty_rejected))
      notes <- c(notes, "a narrow interval after a (near-)rejection is not reassurance: coverage holds only if the restrictions are true")
    if (st$nz == 0) notes <- c(notes, "sign/CR restrictions are not imposed on point-identified rungs (L5, L6)")
    if (jkst && length(extra$cr)) notes <- c(notes, "JKST CI_x intervals are not reported with CR restrictions")
  }
  if (st$vcov_type != vcov) notes <- c(notes, paste0("variance estimator fell back from ", vcov, " to ", st$vcov_type))
  if (any(!st$keep)) notes <- c(notes, sprintf("%d units trimmed; their contribution to the actual rate is bounded by Duncan-Davis, and the model-implied rate refers to the kept units", sum(!st$keep)))
  if (any(est$table$ci_rejected))
    notes <- c(notes, "a confidence interval lies wholly outside the logically possible range: the rung is contradicted by the data (reported as --)")
  if (model %in% c("L5", "L6") && any(est$table$estimate_lower < 0 | est$table$estimate_upper > 1, na.rm = TRUE))
    notes <- c(notes, "point estimate outside [0, 1]: constancy is contradicted by the data")
  out <- structure(list(call = cl, model = model, delta = delta, rows = rl, cols = cll,
                 n = nrow(dat$X), n_kept = sum(st$keep), groups = st$groups, G = st$G,
                 vcov_type = st$vcov_type, alpha = alpha, ci = ci, fpc = fpc,
                 table = est$table, jkst = est$jkst, tests = tst, notes = notes,
                 n_vertices = sum(sapply(st$vertices, nrow)), eta = eta,
                 setup = st, targets = tl, pieces = est$pieces,
                 settings = list(ndraw = ndraw, nz_draw = nz_draw, range = range,
                                 union = union, weights = weights)),
            class = "eiRxC")
  if (is.null(diagnostics)) diagnostics <- length(extra$assoc) + length(extra$cr) > 0
  if (diagnostics) out$sign <- tryCatch(eiRedundancy(out, ci = ci != "none"),
                                        error = function(e) { warning("diagnostics failed: ", conditionMessage(e)); NULL })
  out
}

.rxc_model_label <- function(model, delta = 0, row = NULL) {
  switch(model,
         DD = "Duncan-Davis (no model)",
         M0 = "M0: pooled linear contextual effects",
         M2 = "M2: group-specific linear contextual effects",
         M1delta = sprintf("M1(delta = %g): shared contextual effects, group deviations up to delta", delta),
         M1 = "M1: shared contextual effects with group baselines",
         L4 = sprintf("L4: row '%s' has no contextual effect", row),
         L5 = "L5: constancy within groups (Goodman within each group)",
         L6 = "L6: Goodman constancy")
}

#' @export
print.eiRxC <- function(x, digits = 3, ...) {
  cat("R x C ecological inference bounds\n")
  cat("Model: ", .rxc_model_label(x$model, x$delta, if (!is.null(x$setup$L4row)) x$rows[x$setup$L4row]), "\n", sep = "")
  cat(sprintf("Units: %d (%d used in the model)%s; variance: %s; %g%% intervals (%s)\n",
              x$n, x$n_kept, if (x$G > 1) sprintf(", %d groups", x$G) else "",
              x$vcov_type, 100 * (1 - x$alpha), x$ci))
  if (x$eta > 0) cat(sprintf("Sensitivity: eta = %g (undetectable nonlinearity allowed)\n", x$eta))
  tb <- x$table
  f <- .rxc_fmt_iv(digits)
  g <- function(p) if (paste0(p, "_gap_lower") %in% names(tb))
    f(tb[[paste0(p, "_lower")]], tb[[paste0(p, "_upper")]], tb[[paste0(p, "_gap_lower")]], tb[[paste0(p, "_gap_upper")]]) else
    f(tb[[paste0(p, "_lower")]], tb[[paste0(p, "_upper")]])
  out <- data.frame(target = tb$target,
                    estimate = ifelse(tb$empty, "empty", g("estimate")),
                    ci_model = g("ci_model"), ci_actual = g("ci_actual"),
                    duncan_davis = f(tb$dd_lower, tb$dd_upper), stringsAsFactors = FALSE)
  print(out, row.names = FALSE, right = FALSE)
  ex <- x$setup$extra
  if (length(ex$assoc) || length(ex$cr)) {
    cat("\nSign / CR restrictions (model-implied mean associations, design-weighted):\n")
    for (lb in .rxc_module_labels(list(assoc = ex$assoc, cr = ex$cr), x$rows, x$cols)) cat("  ", lb, "\n")
    if (any(!is.na(tb$estimate_gap_lower)) || any(!is.na(tb$ci_model_gap_lower)))
      cat("  'U' marks a union of two disjoint intervals (CR branches; use union = \"hull\" for the hull)\n")
    cat("  estimate is sharp for the regression (PHASE5a R3(b))\n")
  }
  cat("\n  estimate:  estimated identified set for the model-implied rate (not a confidence interval)\n")
  cat("  ci_model:  CI for the rate implied by the model at the observed compositions\n")
  cat(sprintf("  ci_actual: CI for the actual rate in these units (unit-level allowance: %s)\n", x$fpc))
  cat("  duncan_davis: model-free bounds for the actual rate\n")
  if (!is.null(x$tests)) {
    if (isTRUE(x$tests$empty_rejected)) cat("\n  Emptiness test: the model is REJECTED at level", x$alpha, "\n")
    if (!is.null(x$tests$m1) && !is.na(x$tests$m1$p.value))
      cat(sprintf("  Test of M1's testable part: Wald = %.2f, df = %d, p = %.3f\n",
                  x$tests$m1$statistic, x$tests$m1$df, x$tests$m1$p.value))
  }
  for (nt in x$notes) cat("  Note:", nt, "\n")
  invisible(x)
}

.rxc_fmt_iv <- function(digits) {
  fmt <- paste0("[%.", digits, "f, %.", digits, "f]")
  function(l, u, gl = NULL, gu = NULL) {
    out <- ifelse(is.na(l), "      --      ", sprintf(fmt, l, u))
    if (!is.null(gl)) {
      j <- !is.na(l) & !is.na(gl)
      out[j] <- paste(sprintf(fmt, l[j], gl[j]), "U", sprintf(fmt, gu[j], u[j]))
    }
    out
  }
}

#' @export
summary.eiRxC <- function(object, ...) {
  print(object, ...)
  st <- object$setup
  cat(sprintf("\nComposition range: %s; %d constraint points in total\n",
              object$settings$range, object$n_vertices))
  cat(sprintf("Unknown parameters (contextual-effect splits): %d\n", st$nz))
  if (!is.null(object$sign)) { cat("\n"); print(object$sign) }
  if (!is.null(object$jkst)) {
    cat("\nJiang et al. (2019) CI_x intervals (valid with one active constraint per side):\n")
    for (nm in names(object$jkst)) {
      cat(" ", nm, "\n")
      jk <- object$jkst[[nm]]
      jk[, 1:4] <- round(jk[, 1:4], 4)
      print(jk, row.names = FALSE)
    }
  }
  invisible(object)
}

#' @export
plot.eiRxC <- function(x, ...) {
  tb <- x$table
  k <- nrow(tb)
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, k + 0.5), yaxt = "n", xlab = "rate", ylab = "",
                 main = .rxc_model_label(x$model, x$delta), ...)
  graphics::axis(2, at = k:1, labels = tb$target, las = 1, cex.axis = 0.7)
  for (i in seq_len(k)) {
    y <- k - i + 1
    graphics::segments(tb$dd_lower[i], y, tb$dd_upper[i], y, col = "grey80", lwd = 8)
    if (!is.na(tb$ci_actual_lower[i])) graphics::segments(tb$ci_actual_lower[i], y + 0.15, tb$ci_actual_upper[i], y + 0.15, col = "grey40", lwd = 2)
    if (!is.na(tb$ci_model_lower[i])) graphics::segments(tb$ci_model_lower[i], y - 0.15, tb$ci_model_upper[i], y - 0.15, col = "black", lwd = 2)
    if (!is.na(tb$estimate_lower[i])) graphics::segments(tb$estimate_lower[i], y - 0.15, tb$estimate_upper[i], y - 0.15, col = "black", lwd = 5)
  }
  graphics::legend("bottomright", c("Duncan-Davis", "CI actual", "CI model", "estimate"),
                   col = c("grey80", "grey40", "black", "black"), lwd = c(8, 2, 2, 5), cex = 0.7, bty = "n")
  invisible(x)
}

#' Run the assumption ladder
#'
#' \code{ladder()} computes Duncan-Davis bounds and the bounds and confidence
#' intervals under each rung, from the most to the least plausible, so that
#' readers can see how much each additional assumption buys.
#'
#' @param ... data arguments passed to \code{eiBounds()} (formula/data, rows/cols, or x/t/n, and group).
#' @param rungs rungs to run; rungs needing groups are skipped without groups.
#' @param delta values of delta for M1delta.
#' @param row_constant row for L4 (skipped if \code{NULL}).
#' @param seed random seed.
#' @param sign_rungs optional named list of extra rungs with sign or CR
#'   restrictions, run after the standard rungs; each element is
#'   \code{list(model = "M0", restrictions = list(cr = ..., assoc = ...))}
#'   (model defaults to \code{"M0"}). Example:
#'   \code{list("M0 + CR" = list(restrictions = list(cr = list(list(row = 1, col = 1)))))}.
#' @return a data frame with one row per rung and target.
#' @export
ladder <- function(..., rungs = c("DD", "M0", "M2", "M1delta", "M1", "L4", "L5", "L6"),
                   delta = c(0.10, 0.05, 0.02), row_constant = NULL, seed = 1,
                   sign_rungs = NULL) {
  args <- list(...)
  has_group <- !is.null(args$group)
  out <- list()
  run <- function(model, d = 0, restr = NULL, label = NULL) {
    a <- args; a$model <- model; a$delta <- d; a$seed <- seed; a$row_constant <- row_constant
    a$tests <- FALSE; a$diagnostics <- FALSE
    if (!is.null(restr)) a$restrictions <- c(a$restrictions, restr)
    fit <- tryCatch(do.call(eiBounds, a), error = function(e) e)
    if (inherits(fit, "error")) {
      warning("rung ", if (is.null(label)) model else label, " skipped: ", conditionMessage(fit)); return(NULL)
    }
    tb <- fit$table
    tb$rung <- if (!is.null(label)) label else if (model == "M1delta") sprintf("M1(delta=%g)", d) else model
    tb
  }
  base <- NULL
  for (rg in rungs) {
    if (rg %in% c("M2", "M1delta", "M1", "L5") && !has_group) next
    if (rg == "L4" && is.null(row_constant)) next
    if (rg == "DD") {
      a <- args; a$model <- "M0"; a$ci <- "none"; a$tests <- FALSE; a$seed <- seed
      fit <- do.call(eiBounds, a)
      tb <- fit$table
      tb[, c("estimate_lower", "estimate_upper", "ci_model_lower", "ci_model_upper")] <- NA
      tb[, grep("_gap_", names(tb))] <- NA; tb$n_pieces <- 1L
      tb$ci_actual_lower <- tb$dd_lower; tb$ci_actual_upper <- tb$dd_upper
      tb$rung <- "DD"; tb$empty <- FALSE
      out[[length(out) + 1]] <- tb
    } else if (rg == "M1delta") {
      for (d in sort(delta, decreasing = TRUE)) out[[length(out) + 1]] <- run(rg, d)
    } else out[[length(out) + 1]] <- run(rg)
  }
  for (nm in names(sign_rungs)) {
    sr <- sign_rungs[[nm]]
    m <- if (is.null(sr$model)) "M0" else sr$model
    out[[length(out) + 1]] <- run(m, if (is.null(sr$delta)) 0 else sr$delta, sr$restrictions, nm)
  }
  res <- do.call(rbind, out)
  res <- res[, c("rung", setdiff(names(res), "rung"))]
  rownames(res) <- NULL
  class(res) <- c("eiLadder", "data.frame")
  res
}

#' @export
print.eiLadder <- function(x, digits = 3, ...) {
  f <- .rxc_fmt_iv(digits)
  for (tg in unique(x$target)) {
    s <- x[x$target == tg, ]
    cat("Target:", tg, "\n")
    print(data.frame(rung = s$rung, estimate = ifelse(s$empty, "empty", f(s$estimate_lower, s$estimate_upper, s$estimate_gap_lower, s$estimate_gap_upper)),
                     ci_model = f(s$ci_model_lower, s$ci_model_upper, s$ci_model_gap_lower, s$ci_model_gap_upper),
                     ci_actual = f(s$ci_actual_lower, s$ci_actual_upper, s$ci_actual_gap_lower, s$ci_actual_gap_upper)),
          row.names = FALSE, right = FALSE)
    cat("\n")
  }
  cat("Rungs run from the most plausible (top) to the most restrictive (bottom).\n")
  invisible(x)
}

#' Sensitivity to undetectable nonlinearity
#'
#' Recomputes the bounds of a fitted \code{eiRxC} object allowing average
#' rates to depart from the linear model by up to \code{eta} (in probability
#' units) in directions the data cannot detect (Phase 3, Proposition S1).
#'
#' @param object an \code{eiRxC} object.
#' @param eta grid of values.
#' @param seed random seed.
#' @return a data frame of bounds by \code{eta} and target (class \code{"eiSensitivity"}).
#' @export
sensitivity <- function(object, eta = seq(0, 0.1, by = 0.01), seed = 1) {
  st <- object$setup
  res <- lapply(eta, function(e) {
    set.seed(seed)
    tb <- .rxc_estimate(st, object$targets, alpha = object$alpha, ci = object$ci, fpc = object$fpc,
                        eta = e, ndraw = object$settings$ndraw, nz_draw = object$settings$nz_draw,
                        union = .rxc_union_setting(object))$table
    tb$eta <- e
    tb
  })
  res <- do.call(rbind, res)
  class(res) <- c("eiSensitivity", "data.frame")
  res
}

.rxc_union_setting <- function(object) if (is.null(object$settings$union)) "pieces" else object$settings$union

#' @export
plot.eiSensitivity <- function(x, target = unique(x$target)[1], ...) {
  s <- x[x$target == target, ]
  graphics::plot(s$eta, s$ci_actual_upper, type = "n", ylim = c(0, 1), xlab = "eta", ylab = "rate",
                 main = target, ...)
  graphics::polygon(c(s$eta, rev(s$eta)), c(s$dd_lower, rev(s$dd_upper)), col = "grey90", border = NA)
  graphics::lines(s$eta, s$ci_actual_lower, lty = 2); graphics::lines(s$eta, s$ci_actual_upper, lty = 2)
  graphics::lines(s$eta, s$estimate_lower, lwd = 2); graphics::lines(s$eta, s$estimate_upper, lwd = 2)
  graphics::legend("topright", c("estimate", "CI actual", "Duncan-Davis"), lwd = c(2, 1, 8),
                   lty = c(1, 2, 1), col = c("black", "black", "grey90"), bty = "n", cex = 0.8)
  invisible(x)
}

#' Breakdown value of a conclusion
#'
#' The largest \code{eta} (undetectable departure from linear contextual
#' effects) for which a conclusion of the form "target > threshold" (or
#' "target < threshold") still holds for every value in the chosen interval.
#'
#' @param object an \code{eiRxC} object.
#' @param target name of a target in \code{object$table}.
#' @param threshold the threshold in the conclusion.
#' @param direction \code{"above"} (target > threshold) or \code{"below"}.
#' @param which interval used: \code{"ci_actual"}, \code{"ci_model"} or \code{"estimate"}.
#' @param max_eta upper end of the search.
#' @param tol tolerance of the bisection.
#' @param seed random seed.
#' @return the breakdown value (0 if the conclusion fails already at eta = 0).
#' @export
breakdown <- function(object, target, threshold, direction = c("above", "below"),
                      which = c("ci_actual", "ci_model", "estimate"), max_eta = 0.5,
                      tol = 1e-3, seed = 1) {
  direction <- match.arg(direction); which <- match.arg(which)
  st <- object$setup
  tl <- object$targets[target]
  holds <- function(e) {
    set.seed(seed)
    tb <- .rxc_estimate(st, tl, alpha = object$alpha, ci = if (which == "estimate") "none" else object$ci,
                        fpc = object$fpc, eta = e, ndraw = object$settings$ndraw,
                        nz_draw = object$settings$nz_draw, union = .rxc_union_setting(object))$table
    lo <- tb[[paste0(which, "_lower")]]; hi <- tb[[paste0(which, "_upper")]]
    if (is.na(lo) || isTRUE(tb$empty)) return(FALSE)
    # with a CR union the hull is used: a conclusion holds only if every piece satisfies it
    if (direction == "above") lo > threshold else hi < threshold
  }
  if (!holds(0)) return(0)
  lo <- 0; hi <- max_eta
  if (holds(hi)) return(hi)
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (holds(mid)) lo <- mid else hi <- mid
  }
  lo
}

#' Benchmark rungs on data with known cell values
#'
#' For each dataset with known cell counts, runs the requested rungs and reports
#' whether the true aggregate rates fall inside each interval, and the widths.
#'
#' @param datasets a named list; each element a list with \code{data},
#'   \code{rows}, \code{cols}, optional \code{group}, \code{truth} (an
#'   R x C matrix of true aggregate cell counts, rows x columns) and optional
#'   \code{meta} (a named list of descriptors copied to the output).
#' @param rungs either a character vector of models, or a named list in which
#'   each element is a list of arguments for \code{eiBounds()} (e.g.
#'   \code{list(M0 = list(model = "M0"), "M0+order" = list(model = "M0",
#'   restrictions = ...))}). Restrictions given as a function of the dataset
#'   (\code{function(d) list(...)}) are evaluated per dataset.
#' @param cores number of cores (uses \code{parallel::mclapply} where available).
#' @param ... further arguments passed to \code{eiBounds()} for every rung.
#' @return a data frame with one row per dataset, rung and cell: the truth, the
#'   bounds, and coverage indicators.
#' @examples
#' \dontrun{
#' library(eiPack); data(senc)
#' tr <- matrix(colSums(senc[, c("whdem", "whrep", "whnon", "bldem", "blrep", "blnon",
#'   "natamdem", "natamrep", "natamnon")]), 3, 3, byrow = TRUE)
#' eiBenchmark(list(senc = list(data = senc, rows = c("white", "black", "natam"),
#'   cols = c("dem", "rep", "non"), group = "county", truth = tr)))
#' }
#' @export
eiBenchmark <- function(datasets, rungs = c("M0", "M2", "M1", "L5", "L6"), cores = 1, ...) {
  if (is.character(rungs)) rungs <- stats::setNames(lapply(rungs, function(m) list(model = m)), rungs)
  common <- list(...)
  one <- function(nm) {
    d <- datasets[[nm]]
    truth <- d$truth / rowSums(d$truth)
    tv <- as.vector(t(truth))
    res <- list()
    for (rg in names(rungs)) {
      args <- utils::modifyList(common, rungs[[rg]])
      if (is.function(args$restrictions)) args$restrictions <- args$restrictions(d)
      if (args$model %in% c("M2", "M1delta", "M1", "L5") && is.null(d$group)) next
      args <- c(list(data = d$data, rows = d$rows, cols = d$cols, group = d$group, tests = FALSE), args)
      fit <- tryCatch(do.call(eiBounds, args), error = function(e) e)
      if (inherits(fit, "error")) { warning(nm, "/", rg, ": ", conditionMessage(fit)); next }
      tb <- fit$table
      out <- data.frame(dataset = nm, rung = rg, target = tb$target, truth = tv,
                        estimate_lower = tb$estimate_lower, estimate_upper = tb$estimate_upper, empty = tb$empty,
                        ci_model_lower = tb$ci_model_lower, ci_model_upper = tb$ci_model_upper,
                        ci_actual_lower = tb$ci_actual_lower, ci_actual_upper = tb$ci_actual_upper,
                        dd_lower = tb$dd_lower, dd_upper = tb$dd_upper, stringsAsFactors = FALSE)
      tol <- 1e-9
      out$in_estimate <- !out$empty & tv >= out$estimate_lower - tol & tv <= out$estimate_upper + tol
      out$ci_rejected <- tb$ci_rejected
      out$in_ci_actual <- !is.na(out$ci_actual_lower) & tv >= out$ci_actual_lower - tol & tv <= out$ci_actual_upper + tol
      out$in_dd <- tv >= out$dd_lower - tol & tv <= out$dd_upper + tol
      out$width_estimate <- ifelse(out$empty, NA, out$estimate_upper - out$estimate_lower)
      out$width_ci_actual <- out$ci_actual_upper - out$ci_actual_lower
      out$width_dd <- out$dd_upper - out$dd_lower
      if (!is.null(d$meta)) for (k in names(d$meta)) out[[k]] <- d$meta[[k]]
      res[[length(res) + 1]] <- out
    }
    do.call(rbind, res)
  }
  nms <- names(datasets)
  parts <- if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(nms, one, mc.cores = cores, mc.preschedule = FALSE)
  } else lapply(nms, one)
  bad <- vapply(parts, function(p) inherits(p, "try-error"), logical(1))
  if (any(bad)) warning(sum(bad), " datasets failed")
  out <- do.call(rbind, parts[!bad])
  rownames(out) <- NULL
  out
}

#' Benchmark tables from the ei.Datasets package
#'
#' Builds ecological-inference tables with known truth from the New Zealand
#' general elections in \pkg{ei.Datasets} (Pavia 2023): in each electorate,
#' polling places are the units, the party vote gives the rows, the electorate
#' (candidate) vote gives the columns, and the official electorate-level
#' cross-tabulation of party by candidate vote gives the true table. The
#' chosen parties form their own rows and columns; all other parties,
#' candidates and informal votes are pooled into \code{"OTH"}. Columns with no
#' votes in an electorate (e.g. no National candidate) are dropped.
#'
#' @param elections NZ election years (2002, 2005, ..., 2020).
#' @param parties party names as used in \pkg{ei.Datasets}.
#' @param labels short labels for the parties.
#' @return a named list of benchmark entries for \code{eiBenchmark()}.
#' @references Pavia, J. M. (2023). ei.Datasets: Real Datasets for Assessing
#'   Ecological Inference Algorithms. Social Science Computer Review.
#' @export
eiDatasetsTables <- function(elections = c(2002, 2005, 2008, 2011, 2014, 2017, 2020),
                             parties = c("Labour Party", "National Party"),
                             labels = c("LAB", "NAT")) {
  if (!requireNamespace("ei.Datasets", quietly = TRUE)) stop("package 'ei.Datasets' is required")
  out <- list()
  for (yr in elections) {
    nm <- paste0("ei_NZ_", yr)
    env <- new.env()
    utils::data(list = nm, package = "ei.Datasets", envir = env)
    x <- get(nm, envir = env)
    for (i in seq_len(nrow(x))) {
      P <- as.data.frame(x$Votes_to_parties[[i]]); Cn <- as.data.frame(x$Votes_to_candidates[[i]])
      cv <- as.data.frame(x$District_cross_votes[[i]])
      Pm <- as.matrix(P[, -(1:2)]); Cm <- as.matrix(Cn[, -(1:2)]); Xm <- as.matrix(cv[, -1])
      pc <- colnames(Pm); cc <- colnames(cv)[-1]
      # party (row) and candidate (column) groups; cross-table rows/columns follow the same order
      rg <- rep("OTH", length(pc)); cg <- rep("OTH", length(cc))
      for (k in seq_along(parties)) {
        rg[pc == parties[k]] <- labels[k]
        cg[grepl(paste0("(", parties[k], ")"), cc, fixed = TRUE)] <- labels[k]
      }
      lv <- c(labels, "OTH")
      agg_c <- function(M, g) sapply(lv, function(l) rowSums(M[, g == l, drop = FALSE]))
      R <- agg_c(Pm, rg); C <- agg_c(Cm, cg)
      tr <- sapply(lv, function(l) rowSums(Xm[, cg == l, drop = FALSE]))       # rows: cross rows
      tr <- t(sapply(lv, function(l) colSums(tr[rg == l, , drop = FALSE])))  # aggregate rows by party
      keepc <- colSums(C) > 0; keepr <- colSums(R) > 0
      R <- R[, keepr, drop = FALSE]; C <- C[, keepc, drop = FALSE]; tr <- tr[keepr, keepc, drop = FALSE]
      rn <- paste0("p_", lv[keepr]); cn <- paste0("c_", lv[keepc])
      colnames(R) <- rn; colnames(C) <- cn; dimnames(tr) <- list(rn, cn)
      dat <- as.data.frame(cbind(R, C))
      dat <- dat[rowSums(R) > 0, , drop = FALSE]
      key <- sprintf("NZ%d_%02d_%s", yr, i, gsub("[^A-Za-z]", "", x$District[i]))
      out[[key]] <- list(data = dat, rows = rn, cols = cn, truth = tr,
                         meta = list(election = yr, district = x$District[i], units = nrow(dat)))
    }
  }
  out
}
