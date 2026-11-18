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
#' Reported for each target: \code{estimate} (the estimated outer identified
#' set for the model-implied rate, sharp only relative to the regression
#' function; not a confidence interval), \code{ci_model} (CI for the
#' finite-design conditional-mean rate, i.e. the rate the model implies at the
#' observed compositions), \code{ci_actual} (CI for the realised rate in these
#' units, including an allowance for unit-level deviations bounded by
#' Bhatia-Davis or Popoviciu) and the Duncan-Davis bounds. With \code{ci = "clr"}
#' the table also has \code{proj_model_*} and \code{proj_actual_*}: projection
#' intervals built from the same candidate endpoints.
#'
#' Status of the intervals (PHASE3_inference.md, sections 0.5 and 5): the CLR
#' interval as implemented collects candidate endpoints by simulation and is
#' \emph{algorithmic} (its coverage is a conjecture supported by simulation);
#' the projection interval and the conservative interval have an asymptotic
#' guarantee for a fixed-dimensional regression with a nonsingular covariance
#' matrix \eqn{\hat V} (M0 in any full-rank design; grouped models with few,
#' large groups). If \eqn{\hat V} is numerically singular, candidate endpoints
#' whose gradients use its null directions are dropped from the projection and
#' CLR intervals (a valid, wider interval) and a note says so. The
#' \code{proj_actual} allowance for unit-level deviations uses a variance bound
#' taken from a projection interval, so the projection columns rest only on
#' P3.5(b2). Inference is conditional on the design and on the
#' group effects; units are assumed independent given them, and no clustered
#' variance is used (group-level shocks belong in the mean: use M1 or M2).
#'
#' Sign and contextual-reinforcement restrictions (\code{assoc}, \code{cr};
#' the M-SIGN module) are \strong{experimental} in this release: both referee
#' reports of review batch B are being answered. Under grouped rungs (M1,
#' M1delta, M2, L4 with groups) each sign is one pooled restriction across
#' groups, and sharpness (PHASE5a R3(b)) is proved only for the pooled model M0.
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
#' @param restrictions optional list with elements \code{zeros}, \code{order},
#'   \code{sign}, \code{slope}, \code{linear}, \code{external}, \code{assoc}
#'   and \code{cr}; see the vignette. \code{assoc} is a list of signs on mean
#'   associations, each \code{list(type = "dB" or "dW", row, col, sign = +1 or -1,
#'   version, other)}: \code{"dB"} is the between-group association
#'   \eqn{\delta_B} (version \code{"W1"}, row vs the rest, default; or
#'   \code{"pair"}, row vs \code{other}), \code{"dW"} the within-group
#'   (contextual) association \eqn{\delta_W} (version \code{"W1"}, one-vs-rest,
#'   default; or \code{"W2"}, \eqn{\Delta_r}, the within-group associations
#'   of every original group with the row's share, \emph{pooled} into one sum:
#'   its sign does not sign each group's association). \code{cr} is a list of
#'   contextual-reinforcement restrictions (\eqn{\delta_B} and \eqn{\delta_W}
#'   have the same sign, which is not assumed known), each \code{list(row, col,
#'   version = "W1", "W2" or "pair", other)}. W1 (and every 2 x 2 CR) is one LP
#'   (a slab); W2 and pair are a union of two branches and may give two disjoint
#'   intervals. All associations are model-implied (LCE) means with design
#'   weights \eqn{N_i}. Signs on \eqn{\delta_W} can be refuted by the data.
#'   A sign on \eqn{\delta_B} and CR are always feasible only in the
#'   model-free problem with continuous unit tables (Duncan-Davis margins, no
#'   linear rung, no other restrictions); with the linear rung they can be
#'   rejected jointly with it, and with integer counts or other restrictions
#'   CR can be refuted outright. An empty set rejects the combination; it does
#'   not say which part fails (the fit's notes say whether the rung alone is
#'   rejected).
#' @param alpha 1 - confidence level.
#' @param ci \code{"clr"} (default; reported with the projection interval
#'   alongside), \code{"conservative"} or \code{"none"}.
#' @param fpc unit-level allowance for \code{ci_actual}: \code{"bhatia-davis"} (default) or \code{"popoviciu"}.
#' @param jkst also report Jiang et al. CI_x intervals.
#' @param jkst_x values of x for CI_x.
#' @param eta sensitivity: allowed undetectable nonlinearity (probability units).
#' @param ndraw draws used to collect candidate endpoints.
#' @param nz_draw Gaussian draws for critical values.
#' @param seed optional random seed (for reproducible CIs).
#' @param tests run the emptiness test and, for grouped M1/L4, the test of M1's
#'   testable part (against M2) and the pooling test (M0 against M1: all group
#'   shifts zero).
#' @param range_offset shrink the observed range by this amount (two-row tables
#'   only; Jiang et al. use 1e-5).
#' @param mismatch \code{"error"} or \code{"rescale"} when row and column totals differ.
#' @param union with CR restrictions the estimate and CIs can be a union of
#'   disjoint intervals: \code{"pieces"} (default) reports the hull in
#'   \code{*_lower}/\code{*_upper}, the number of pieces in \code{n_pieces} and
#'   the \emph{widest} excluded gap in \code{*_gap_lower}/\code{*_gap_upper};
#'   \code{print()} shows every piece, and \code{fit$pieces} holds them all;
#'   \code{"hull"} reports the hull only (valid, wider or equal).
#' @param diagnostics compute \code{eiRedundancy()} diagnostics and store them in
#'   \code{$sign} (default: when \code{assoc} or \code{cr} restrictions are given).
#' @param base also fit the rung alone (without the optional restrictions) and
#'   report its estimate and CIs in the \code{base_*} columns, printed beside
#'   the restricted result (default: when restrictions are given).
#' @return an object of class \code{"eiRxC"}.
#' @references Jiang, W., King, G., Schmaltz, A. and Tanner, M. A. (2019).
#'   Ecological Regression with Partial Identification. Political Analysis.
#' @examples
#' if (requireNamespace("eco", quietly = TRUE)) {
#'   data("census", package = "eco")
#'   fit <- eiBounds(x = census$X, t = census$Y, n = census$N, model = "M0", seed = 1)
#'   fit
#' }
#' @export
eiBounds <- function(formula = NULL, data = NULL, rows = NULL, cols = NULL,
                     x = NULL, t = NULL, n = NULL, group = NULL,
                     model = c("M0", "M1", "M1delta", "M2", "L4", "L5", "L6"),
                     delta = 0, row_constant = NULL, targets = NULL,
                     range = c("hull", "trim", "extrapolate"), trim = NULL, extrapolate = NULL,
                     weights = c("none", "n"), vcov = c("HC3", "HC1", "HC0"),
                     restrictions = NULL, alpha = 0.05, ci = c("clr", "conservative", "none"),
                     fpc = c("bhatia-davis", "popoviciu"), jkst = FALSE,
                     jkst_x = seq(0, 2, by = 0.25), eta = 0, ndraw = 100, nz_draw = 4000,
                     seed = NULL, tests = TRUE, range_offset = 0,
                     mismatch = c("error", "rescale"), union = c("pieces", "hull"),
                     diagnostics = NULL, base = NULL) {
  cl <- match.call()
  union <- match.arg(union)
  model <- match.arg(model); range <- match.arg(range); weights <- match.arg(weights)
  vcov <- match.arg(vcov); ci <- match.arg(ci); fpc <- match.arg(fpc); mismatch <- match.arg(mismatch)
  if (!is.null(seed)) set.seed(seed)
  dat <- .rxc_parse(formula, data, rows, cols, x, t, n, group, mismatch)
  rl <- colnames(dat$X); cll <- colnames(dat$T)
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
                   extra = extra)
  tl <- .rxc_targets(targets, rl, cll)
  est <- .rxc_estimate(st, tl, alpha = alpha, ci = ci, fpc = fpc, jkst = jkst, jkst_x = jkst_x,
                       eta = eta, ndraw = ndraw, nz_draw = nz_draw, union = union)
  tst <- NULL
  if (tests) {
    ed <- .rxc_empty_test(st, alpha = alpha, nz_draw = nz_draw, detail = TRUE)
    tst <- list(empty_rejected = if (is.list(ed)) ed$rejected else ed,
                empty_detail = if (is.list(ed)) ed else NULL)
    if (model %in% c("M1", "L4") && st$G > 1)
    {
      tst$m1 <- .rxc_m1_test(dat$X[st$keep, , drop = FALSE], dat$T[st$keep, , drop = FALSE],
                             dat$N[st$keep], dat$group[st$keep], weights = weights, vcov = vcov)
      tst$pool <- .rxc_pool_test(dat$X[st$keep, , drop = FALSE], dat$T[st$keep, , drop = FALSE],
                                 dat$N[st$keep], dat$group[st$keep], weights = weights, vcov = vcov)
    }
  }
  if (range == "extrapolate")
    warning("range = \"extrapolate\": the linear model is assumed where no units are observed; this cannot be checked")
  notes <- character(0)
  mlabels <- .rxc_module_labels(extra, rl, cll)
  combo <- if (length(mlabels)) paste0("the rung (", model, ") + ", paste(mlabels, collapse = " + ")) else
    paste0("the rung (", model, ")")
  if (any(est$table$empty))
    notes <- c(notes, paste0("empty estimate: at the point estimate ", combo, " has no solution. ",
                             "This is a sample diagnostic, not a test (see the emptiness test). An empty ",
                             "set rejects the combination; it does not say which part fails"))
  if (isTRUE(tst$empty_rejected)) {
    ra <- tst$empty_detail$rung_alone_rejected
    notes <- c(notes, paste0("emptiness test rejects ", combo,
                             if (isTRUE(ra)) "; the rung alone is also rejected" else
                               if (identical(ra, FALSE)) "; the rung alone is not rejected" else ""))
  }
  sign_like <- c(vapply(extra$assoc, function(a) a$type == "dB", TRUE), rep(TRUE, length(extra$cr)))
  if (length(sign_like) && all(sign_like) &&
      (any(est$table$empty) || isTRUE(tst$empty_rejected)))
    notes <- c(notes, paste0("a sign on delta_B and CR are always feasible only in the model-free problem with ",
                             "continuous tables (no linear rung, no other restrictions). With the rung (", model,
                             ") they can be rejected jointly with it, and with integer counts CR can be ",
                             "refuted outright. The rejection is of the combination; it is not evidence that the belief holds"))
  if (length(extra$assoc) || length(extra$cr)) {
    if (st$G > 1 && model %in% c("M1", "M2", "M1delta", "L4"))
      notes <- c(notes, paste0("sign/CR restrictions under a grouped rung (", model, ") are experimental: ",
                               "each sign is one pooled restriction across groups, not one per group, and ",
                               "sharpness (PHASE5a R3(b)) is proved only for the pooled model M0"))
    if (any(vapply(c(extra$assoc, extra$cr), function(a) identical(a$version, "W2"), TRUE)))
      notes <- c(notes, paste0("W2 signs one pooled sum over the original groups; it does not sign each ",
                               "group's association (a per-group version is planned)"))
    if (any(est$table$ci_rejected) || any(est$table$empty) || isTRUE(tst$empty_rejected))
      notes <- c(notes, "a narrow interval after a (near-)rejection is not reassurance: coverage holds only if the restrictions are true")
    if (st$nz == 0) notes <- c(notes, "sign/CR restrictions are not imposed on point-identified rungs (L5, L6)")
    if (jkst && length(extra$cr)) notes <- c(notes, "JKST CI_x intervals are not reported with CR restrictions")
  }
  if (isTRUE(est$slack < 0.005) && !any(est$table$empty))
    notes <- c(notes, sprintf("near the boundary: for every K in the estimated set some implied rate is within %.3f of 0 at a range vertex. The CIs remain valid, but the rate result for the plug-in (P3 Prop. P3.2) assumes strict interiority", max(est$slack, 0)))
  if (st$G > 1 && model %in% c("M1", "M2", "M1delta", "L4")) {
    nmed <- stats::median(tabulate(st$g, st$G))
    if (nmed < 30) notes <- c(notes, sprintf("few units per group (median %g): asymptotic coverage is established only for few, large groups; here the CLR interval's coverage rests on simulation (P3 section 0.5)", nmed))
  }
  if (any(is.na(est$table$estimate_lower) & !est$table$empty))
    notes <- c(notes, "unbounded estimate: some composition range is not full-dimensional in the simplex (2a Lemma 3(b)); add units with more varied compositions or restrict the targets")
  if (ci == "clr")
    notes <- c(notes, "CLR intervals are algorithmic (candidate endpoints collected by simulation; coverage supported by simulation, not proved). proj_* columns: projection intervals from the same candidates, with an asymptotic guarantee for a nonsingular V-hat in a fixed-dimension design (P3 Prop. P3.5(b2))")
  if (isTRUE(est$proj_dropped > 0))
    notes <- c(notes, sprintf("the estimated covariance matrix of the regression is numerically singular in directions that %d candidate endpoint(s) use; they were dropped from the CLR and projection intervals (valid, wider; P3.5(b2)(ii) covers only the range of V-hat)", est$proj_dropped))
  if (st$vcov_type != vcov) notes <- c(notes, paste0("variance estimator fell back from ", vcov, " to ", st$vcov_type))
  if (any(!st$keep)) notes <- c(notes, sprintf("%d units trimmed; their contribution to the actual rate is bounded by Duncan-Davis, and the model-implied rate refers to the kept units", sum(!st$keep)))
  if (any(est$table$ci_rejected))
    notes <- c(notes, paste0("a confidence interval is empty or lies wholly outside the logically possible range: ",
                             combo, " is contradicted by the data (reported as --)"))
  if (model %in% c("L5", "L6") && any(est$table$estimate_lower < 0 | est$table$estimate_upper > 1, na.rm = TRUE))
    notes <- c(notes, "point estimate outside [0, 1]: constancy is contradicted by the data")
  out <- structure(list(call = cl, model = model, delta = delta, rows = rl, cols = cll,
                 n = nrow(dat$X), n_kept = sum(st$keep), groups = st$groups, G = st$G,
                 vcov_type = st$vcov_type, alpha = alpha, ci = ci, fpc = fpc,
                 table = est$table, jkst = est$jkst, tests = tst, notes = notes,
                 n_vertices = sum(sapply(st$vertices, nrow)), eta = eta, slack = est$slack,
                 setup = st, targets = tl, pieces = est$pieces,
                 modules = .rxc_module_table(extra, rl, cll),
                 settings = list(ndraw = ndraw, nz_draw = nz_draw, range = range,
                                 union = union, weights = weights)),
            class = "eiRxC")
  if (is.null(diagnostics)) diagnostics <- length(extra$assoc) + length(extra$cr) > 0
  if (diagnostics) out$sign <- tryCatch(eiRedundancy(out, ci = ci != "none"),
                                        error = function(e) { warning("diagnostics failed: ", conditionMessage(e)); NULL })
  # the rung alone, beside the restricted result (computed last: it does not change the draws
  # used above)
  if (is.null(base)) base <- !is.null(extra) && length(extra) > 0
  if (base && !is.null(extra)) {
    st0 <- st; st0$extra <- NULL
    if (!is.null(seed)) set.seed(seed)      # the same draws as a separate fit of the rung alone
    b <- .rxc_estimate(st0, tl, alpha = alpha, ci = ci, fpc = fpc, eta = eta, ndraw = ndraw,
                       nz_draw = nz_draw, union = union)$table
    for (v in c("estimate_lower", "estimate_upper", "ci_model_lower", "ci_model_upper",
                "ci_actual_lower", "ci_actual_upper"))
      out$table[[paste0("base_", v)]] <- b[[v]]
    out$table$base_empty <- b$empty
    if (any(est$table$empty))
      out$notes <- c(out$notes, if (any(b$empty)) "the rung alone also has an empty estimate" else
        paste0("the rung alone (", model, ") has a non-empty estimate (base_* columns)"))
  }
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
  hull <- identical(x$settings$union, "hull")
  # every piece of a union (P-E15): from fit$pieces when available, else hull + widest gap
  g <- function(p) {
    out <- if (paste0(p, "_gap_lower") %in% names(tb))
      f(tb[[paste0(p, "_lower")]], tb[[paste0(p, "_upper")]], tb[[paste0(p, "_gap_lower")]], tb[[paste0(p, "_gap_upper")]]) else
      f(tb[[paste0(p, "_lower")]], tb[[paste0(p, "_upper")]])
    if (!hull && !is.null(x$pieces)) for (i in seq_len(nrow(tb))) {
      u <- x$pieces[[tb$target[i]]][[p]]
      if (!is.null(u) && nrow(u) > 2 && !is.na(tb[[paste0(p, "_lower")]][i]))
        out[i] <- paste(f(u[, 1], u[, 2]), collapse = " U ")
    }
    out
  }
  out <- data.frame(target = tb$target,
                    estimate = ifelse(tb$empty, "empty", g("estimate")),
                    ci_model = g("ci_model"), ci_actual = g("ci_actual"),
                    duncan_davis = f(tb$dd_lower, tb$dd_upper), stringsAsFactors = FALSE)
  has_proj <- "proj_model_lower" %in% names(tb) && any(!is.na(tb$proj_model_lower)) && identical(x$ci, "clr")
  if (has_proj) {
    out$proj_model <- f(tb$proj_model_lower, tb$proj_model_upper)
    out$proj_actual <- f(tb$proj_actual_lower, tb$proj_actual_upper)
    out <- out[, c("target", "estimate", "ci_model", "proj_model", "ci_actual", "proj_actual", "duncan_davis")]
  }
  print(out, row.names = FALSE, right = FALSE)
  has_base <- "base_estimate_lower" %in% names(tb)
  if (has_base) {
    cat("\nThe rung alone (", x$model, ", without the optional restrictions), for comparison:\n", sep = "")
    print(data.frame(target = tb$target,
                     estimate = ifelse(tb$base_empty, "empty", f(tb$base_estimate_lower, tb$base_estimate_upper)),
                     ci_model = f(tb$base_ci_model_lower, tb$base_ci_model_upper),
                     ci_actual = f(tb$base_ci_actual_lower, tb$base_ci_actual_upper), stringsAsFactors = FALSE),
          row.names = FALSE, right = FALSE)
  }
  ex <- x$setup$extra
  has_cr <- length(ex$cr) > 0
  if (length(ex$assoc) || has_cr) {
    cat("\nSign / CR restrictions (model-implied mean associations, design-weighted):\n")
    for (lb in .rxc_module_labels(list(assoc = ex$assoc, cr = ex$cr), x$rows, x$cols)) cat("  ", lb, "\n")
    if (any(tb$n_pieces > 1, na.rm = TRUE) || any(!is.na(tb$ci_model_gap_lower)))
      cat("  'U' marks a union of disjoint intervals (CR branches; every piece is shown; use union = \"hull\" for the hull)\n")
    if (identical(x$setup$type, "pooled"))
      cat("  estimate is sharp for the regression (PHASE5a R3(b), pooled model)\n") else
      cat("  estimate is an outer set for the regression; sharpness under grouped rungs is not yet proved\n")
    if (any(vapply(c(ex$assoc, ex$cr), function(a) identical(a$version, "W2"), TRUE)))
      cat("  W2 is one pooled sum over the original groups, not a sign for each group\n")
  }
  if (!is.null(x$modules)) {
    exp <- x$modules[!grepl("^proved", x$modules$proof_status), , drop = FALSE]
    for (k in seq_len(nrow(exp))) cat(sprintf("  EXPERIMENTAL module %s: %s\n", exp$restriction[k], exp$proof_status[k]))
  }
  clr <- identical(x$ci, "clr")
  cat("\n  estimate:  estimated outer identified set for the model-implied rate (not a confidence interval)\n")
  cat("  ci_model:  CI for the rate implied by the model at the observed compositions",
      if (clr) "\n             (CLR: algorithmic, coverage supported by simulation, not proved)" else "", "\n", sep = "")
  cat(sprintf("  ci_actual: CI for the actual rate in these units (unit-level allowance: %s)%s\n", x$fpc,
              if (clr) "; CLR, algorithmic" else ""))
  if (has_cr && clr) cat("  with CR, ci_* are the union of the CLR branch intervals (algorithmic); proj_* are the\n",
                         "  union of the projection branch intervals\n", sep = "")
  if (has_proj) cat("  proj_*:    projection intervals (asymptotic guarantee for a nonsingular V-hat in a\n",
                    "             fixed-dimension design: M0, or few large groups; P3.5(b2))\n", sep = "")
  cat("  duncan_davis: model-free bounds for the actual rate\n")
  if (!is.null(x$tests)) {
    if (isTRUE(x$tests$empty_rejected)) cat("\n  Emptiness test: the model is REJECTED at level", x$alpha, "\n")
    if (!is.null(x$tests$m1) && !is.na(x$tests$m1$p.value))
      cat(sprintf("  Test of M1's testable part: Wald = %.2f, df = %d, p = %.3f\n",
                  x$tests$m1$statistic, x$tests$m1$df, x$tests$m1$p.value))
    if (!is.null(x$tests$pool) && !is.na(x$tests$pool$p.value)) {
      pm <- x$tests$pool$method
      cat(sprintf("  Pooling test (M0 vs M1, group shifts = 0): Wald = %.4g, df = %d, p = %.3f (%s)%s\n",
                  x$tests$pool$statistic, x$tests$pool$df, x$tests$pool$p.value,
                  if (identical(pm, "wild")) sprintf("wild bootstrap, B = %d", x$tests$pool$B) else
                    if (is.null(pm)) "chi-squared" else pm,
                  if (isTRUE(x$tests$pool$spanning)) "" else " (some group's compositions are not affinely spanning)"))
    }
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
#' @param delta values of delta for M1delta. \code{delta = 0} is a sensitivity
#'   benchmark rather than a substantive rung: it keeps the shared contextual
#'   effects of M1 but drops M1's testable restrictions on the symmetric parts,
#'   so it can be wider than M1.
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
    a$tests <- FALSE; a$diagnostics <- FALSE; a$base <- FALSE
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

#' Sensitivity to undetectable nonlinearity, or to the slack on sign beliefs
#'
#' Recomputes the bounds of a fitted \code{eiRxC} object along a dial.
#' \code{dial = "eta"} (default) allows average rates to depart from the linear
#' model by up to \code{eta} (in probability units) in directions the data
#' cannot detect (Phase 3, Proposition S1; the result is a valid outer bound,
#' not the sharp set of the relaxed model). \code{dial = "slack"} relaxes the
#' sign and contextual-reinforcement restrictions (\code{assoc}, \code{cr}):
#' each association may have the "wrong" sign by up to \code{slack}, measured
#' as an average within-unit gap (\eqn{\delta_B}) or an average slope in own
#' share (\eqn{\delta_W}), in probability units. Bounds widen monotonically in
#' either dial.
#'
#' @param object an \code{eiRxC} object.
#' @param eta grid of values of the dial.
#' @param seed random seed.
#' @param dial \code{"eta"} or \code{"slack"}.
#' @param modules for \code{dial = "slack"}: module names (e.g. \code{"cr[1]"},
#'   \code{"assoc[2]"}; see \code{fit$modules}) whose slack is set; default all
#'   \code{assoc} and \code{cr} modules.
#' @return a data frame of bounds by dial value and target (class
#'   \code{"eiSensitivity"}); the dial value is in column \code{eta} or \code{slack}.
#' @export
sensitivity <- function(object, eta = seq(0, 0.1, by = 0.01), seed = 1, dial = c("eta", "slack"),
                        modules = NULL) {
  dial <- match.arg(dial)
  res <- lapply(eta, function(e) {
    set.seed(seed)
    tb <- .rxc_dial_estimate(object, object$targets, dial, e, modules, ci = object$ci)$table
    tb[[dial]] <- e
    tb
  })
  res <- do.call(rbind, res)
  class(res) <- c("eiSensitivity", "data.frame")
  attr(res, "dial") <- dial
  res
}

.rxc_union_setting <- function(object) if (is.null(object$settings$union)) "pieces" else object$settings$union

# set the slack of the sign / CR modules (all, or those named) to a common value
.rxc_set_slack <- function(ex, value, modules = NULL) {
  for (ty in c("assoc", "cr")) for (i in seq_along(ex[[ty]])) {
    if (is.null(modules) || paste0(ty, "[", i, "]") %in% modules) ex[[ty]][[i]]$slack <- value
  }
  ex
}

# re-estimate a fit at one value of a dial
.rxc_dial_estimate <- function(object, targets, dial, value, modules = NULL, ci = object$ci) {
  st <- object$setup
  eta <- object$eta
  if (dial == "eta") eta <- value else {
    if (!length(st$extra$assoc) && !length(st$extra$cr)) stop("dial = \"slack\" needs assoc or cr restrictions")
    st$extra <- .rxc_set_slack(st$extra, value, modules)
  }
  .rxc_estimate(st, targets, alpha = object$alpha, ci = ci, fpc = object$fpc, eta = eta,
                ndraw = object$settings$ndraw, nz_draw = object$settings$nz_draw,
                union = .rxc_union_setting(object))
}

#' @export
plot.eiSensitivity <- function(x, target = unique(x$target)[1], ...) {
  s <- x[x$target == target, ]
  dial <- if (is.null(attr(x, "dial"))) "eta" else attr(x, "dial")
  dv <- s[[dial]]
  graphics::plot(dv, s$ci_actual_upper, type = "n", ylim = c(0, 1), xlab = dial, ylab = "rate",
                 main = target, ...)
  graphics::polygon(c(dv, rev(dv)), c(s$dd_lower, rev(s$dd_upper)), col = "grey90", border = NA)
  graphics::lines(dv, s$ci_actual_lower, lty = 2); graphics::lines(dv, s$ci_actual_upper, lty = 2)
  graphics::lines(dv, s$estimate_lower, lwd = 2); graphics::lines(dv, s$estimate_upper, lwd = 2)
  graphics::legend("topright", c("estimate", "CI actual", "Duncan-Davis"), lwd = c(2, 1, 8),
                   lty = c(1, 2, 1), col = c("black", "black", "grey90"), bty = "n", cex = 0.8)
  invisible(x)
}

#' Breakdown value of a conclusion
#'
#' The largest value of a dial for which a conclusion of the form "target >
#' threshold" (or "target < threshold") still holds for every value in the
#' chosen interval. With \code{dial = "eta"} the dial is the undetectable
#' departure from linear contextual effects; with \code{dial = "slack"} it is
#' the slack on the sign / contextual-reinforcement beliefs (see
#' \code{sensitivity()}), the natural dial for beliefs the data cannot refute.
#'
#' @param object an \code{eiRxC} object.
#' @param target name of a target in \code{object$table}.
#' @param threshold the threshold in the conclusion.
#' @param direction \code{"above"} (target > threshold) or \code{"below"}.
#' @param which interval used: \code{"ci_actual"}, \code{"ci_model"} or \code{"estimate"}.
#' @param max_eta upper end of the search.
#' @param tol tolerance of the bisection.
#' @param seed random seed.
#' @param dial \code{"eta"} or \code{"slack"}.
#' @param modules for \code{dial = "slack"}: modules whose slack is varied (default all sign / CR modules).
#' @return the breakdown value (0 if the conclusion fails already at 0).
#' @export
breakdown <- function(object, target, threshold, direction = c("above", "below"),
                      which = c("ci_actual", "ci_model", "estimate"), max_eta = 0.5,
                      tol = 1e-3, seed = 1, dial = c("eta", "slack"), modules = NULL) {
  direction <- match.arg(direction); which <- match.arg(which); dial <- match.arg(dial)
  tl <- object$targets[target]
  holds <- function(e) {
    set.seed(seed)
    tb <- .rxc_dial_estimate(object, tl, dial, e, modules,
                             ci = if (which == "estimate") "none" else object$ci)$table
    lo <- tb[[paste0(which, "_lower")]]; hi <- tb[[paste0(which, "_upper")]]
    if (is.na(lo) || isTRUE(tb$empty)) return(FALSE)
    # with a CR union the hull is used: a conclusion holds only if every piece satisfies it
    if (direction == "above") lo > threshold else hi < threshold
  }
  base <- 0
  if (!holds(base)) return(0)
  lo <- base; hi <- max_eta
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
      if (is.null(args$base)) args$base <- FALSE
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
