# Module registry for optional restrictions (PACKAGE_DESIGN.md §1.1).
#
# Each restriction type is one registry entry:
#   parse(spec, ctx)            user spec (labels allowed) -> internal spec (indices)
#   rows(st, acc, spec, md)     append linear rows s(beta) + a'z >= 0 (form "rows"), or
#   union                       TRUE for either/or modules (built by .rxc_cr_parts)
#   label(spec, rl, cl)         one-line description used in messages and tables
#   meta                        source, object, target, refutable, form, sharp, dial,
#                               proof_status (PACKAGE_DESIGN.md §1.1-1.2, §4)
# 'scalar' entries take a single value rather than a list of specs.
# eiBounds() parses restrictions through the registry, .rxc_extra_rows() builds the rows,
# and fit$modules / eiModules() report the meta.

.rxc_meta <- function(source, object, target, refutable, form, sharp, dial, proof_status)
  list(source = source, object = object, target = target, refutable = refutable, form = form,
       sharp = sharp, dial = dial, proof_status = proof_status)

.rxc_rate_diff <- function(st, g, c, r1, v1, r2, v2) {
  a1 <- .rxc_rate(st, g, c, r1, v1); a2 <- .rxc_rate(st, g, c, r2, v2)
  list(coef = a1$coef - a2$coef, const = a1$const - a2$const, a = a1$a - a2$a)
}

.rxc_slack_arg <- function(s, what) {
  sl <- if (is.null(s[["slack"]])) 0 else as.numeric(s[["slack"]])
  if (length(sl) != 1 || !is.finite(sl) || sl < 0) stop(what, ": 'slack' must be a nonnegative number")
  sl
}

.rxc_version_arg <- function(v, allowed, what) {
  v <- if (is.null(v)) "W1" else as.character(v)
  if (!v %in% allowed) stop(what, ": version must be one of ", paste(allowed, collapse = ", "))
  v
}

.rxc_other_arg <- function(o, r, v, what, ctx) {
  if (v != "pair") return(NULL)
  if (is.null(o)) {
    if (length(ctx$rl) == 2) return(setdiff(1:2, r))
    stop(what, ": version \"pair\" needs 'other'")
  }
  o <- ctx$ri(o)
  if (o == r) stop(what, ": 'other' must differ from 'row'")
  o
}

.rxc_registry <- list(
  zeros = list(
    parse = function(z, ctx) list(row = ctx$ri(z[["row"]]), col = ctx$ci(z[["col"]]),
                                  group = if ("group" %in% names(z)) ctx$gi(z[["group"]]) else NULL),
    rows = function(st, acc, zr, md) {
      # structural zero: rate of row r to column c is 0 throughout the range
      for (g in if (is.null(zr$group)) seq_len(st$G) else zr$group) {
        Vg <- st$vertices[[g]]
        for (k in seq_len(nrow(Vg))) {
          rt <- .rxc_rate(st, g, zr$col, zr$row, Vg[k, ])
          acc <- .rxc_add_row(acc, -rt$coef, -rt$const, -rt$a, "zero", md)
        }
      }
      acc
    },
    label = function(s, rl, cl) sprintf("structural zero %s -> %s", rl[s$row], cl[s$col]),
    meta = .rxc_meta("known cells (external)", "mean rate on the composition range", "model",
                     "yes (empty set)", "rows", "sharp for the regression", "none", "proved (PHASE4)")),
  order = list(
    parse = function(o, ctx) list(high = ctx$ri(o[["high"]]), low = ctx$ri(o[["low"]]), col = ctx$ci(o[["col"]])),
    rows = function(st, acc, od, md) {
      # mean rate of row 'high' to column c >= that of row 'low', at every vertex
      for (g in seq_len(st$G)) {
        Vg <- st$vertices[[g]]
        for (k in seq_len(nrow(Vg))) {
          d <- .rxc_rate_diff(st, g, od$col, od$high, Vg[k, ], od$low, Vg[k, ])
          acc <- .rxc_add_row(acc, d$coef, d$const, d$a, "order", md)
        }
      }
      acc
    },
    label = function(s, rl, cl) sprintf("ordering %s >= %s for '%s'", rl[s$high], rl[s$low], cl[s$col]),
    meta = .rxc_meta("belief (M-ORD)", "mean rate at every composition", "model", "yes (empty set)",
                     "rows", "sharp for the regression", "none", "proved (2c M5(b), M2)")),
  sign = list(
    parse = function(s, ctx) list(row = ctx$ri(s[["row"]]), col = ctx$ci(s[["col"]]), from = ctx$ri(s[["from"]]),
                                  to = ctx$ri(s[["to"]]), sign = sign(as.numeric(s[["sign"]]))),
    rows = function(st, acc, sg, md) {
      # moving composition from 'from' to 'to' changes row r's mean rate in direction 'sign'
      e <- diag(st$R)
      for (g in seq_len(st$G)) {
        d <- .rxc_rate_diff(st, g, sg$col, sg$row, e[sg$to, ], sg$row, e[sg$from, ])
        acc <- .rxc_add_row(acc, sg$sign * d$coef, sg$sign * d$const, sg$sign * d$a, "sign", md)
      }
      acc
    },
    label = function(s, rl, cl) sprintf("contextual-effect sign (row %s, col %s)", rl[s$row], cl[s$col]),
    meta = .rxc_meta("belief", "contextual-effect matrix entry", "model", "yes (empty set)", "rows",
                     "sharp for the regression", "none", "proved (PHASE4)")),
  slope = list(
    scalar = TRUE,
    parse = function(s, ctx) as.numeric(s),
    rows = function(st, acc, lam, md) {
      # bounded slopes: |H_c[r, s] - H_c[r, s']| <= lambda for all r, c, s, s'
      e <- diag(st$R)
      for (g in seq_len(st$G)) for (c in seq_len(st$C)) for (r in seq_len(st$R)) {
        for (s1 in seq_len(st$R - 1)) for (s2 in (s1 + 1):st$R) {
          d <- .rxc_rate_diff(st, g, c, r, e[s2, ], r, e[s1, ])
          acc <- .rxc_add_row(acc, d$coef, d$const + lam, d$a, "slope", md)
          acc <- .rxc_add_row(acc, -d$coef, -d$const + lam, -d$a, "slope", md)
        }
      }
      acc
    },
    label = function(s, rl, cl) sprintf("bounded slopes (%g)", s),
    meta = .rxc_meta("belief", "contextual-effect matrix", "model", "yes (empty set)", "rows",
                     "sharp for the regression", "the bound itself", "proved (PHASE4)")),
  linear = list(
    parse = function(l, ctx) {
      tm <- l$terms
      m <- cbind(row = ctx$ri(tm$row), col = ctx$ci(tm$col), context = ctx$ri(tm$context), value = as.numeric(tm$value))
      list(terms = m, rhs = as.numeric(l$rhs), group = ctx$gi(l$group))
    },
    rows = function(st, acc, lc, md) {
      # sum value * H_c[r, s] >= rhs
      e <- diag(st$R)
      for (g in if (is.null(lc$group)) seq_len(st$G) else lc$group) {
        coef <- numeric(st$nbeta); const <- 0; a <- numeric(st$nz)
        for (q in seq_len(nrow(lc$terms))) {
          tq <- lc$terms[q, ]
          rt <- .rxc_rate(st, g, tq[["col"]], tq[["row"]], e[tq[["context"]], ])
          coef <- coef + tq[["value"]] * rt$coef; const <- const + tq[["value"]] * rt$const
          a <- a + tq[["value"]] * rt$a
        }
        acc <- .rxc_add_row(acc, coef, const - lc$rhs, a, "linear", md)
      }
      acc
    },
    label = function(s, rl, cl) "linear restriction",
    meta = .rxc_meta("belief", "contextual-effect matrix", "model", "yes (empty set)", "rows",
                     "sharp for the regression", "none", "proved (PHASE4)")),
  external = list(
    parse = function(e, ctx) {
      a <- if (is.matrix(e$target)) e$target else {
        a <- matrix(0, length(ctx$rl), length(ctx$cl)); a[ctx$ri(e$target[1]), ctx$ci(e$target[2])] <- 1; a }
      list(weights = a, lower = e$lower, upper = e$upper)
    },
    rows = function(st, acc, eb, md) {
      # lower <= phi <= upper for a target functional
      tg <- .rxc_target(st, eb$weights)
      acc <- .rxc_add_row(acc, tg$cvec, tg$c0 - eb$lower, tg$gamma, "external", md)
      .rxc_add_row(acc, -tg$cvec, eb$upper - tg$c0, -tg$gamma, "external", md)
    },
    label = function(s, rl, cl) sprintf("external benchmark [%g, %g]", s$lower, s$upper),
    meta = .rxc_meta("external (survey, census)", "model-implied target", "model", "yes (empty set)", "rows",
                     "sharp given the benchmark", "width of [lower, upper]", "proved (PHASE4)")),
  assoc = list(
    parse = function(s, ctx) {
      type <- as.character(s[["type"]])
      if (length(type) != 1 || !type %in% c("dB", "dW")) stop("assoc: type must be \"dB\" or \"dW\"")
      v <- .rxc_version_arg(s[["version"]], if (type == "dB") c("W1", "pair") else c("W1", "W2"),
                            paste0("assoc (", type, ")"))
      r <- ctx$ri(s[["row"]])
      sg <- sign(as.numeric(if (is.null(s[["sign"]])) 1 else s[["sign"]]))
      if (sg == 0) stop("assoc: sign must be +1 or -1")
      list(type = type, version = v, row = r, col = ctx$ci(if (is.null(s[["col"]])) 1 else s[["col"]]),
           sign = sg, other = .rxc_other_arg(s[["other"]], r, v, "assoc", ctx),
           slack = .rxc_slack_arg(s, "assoc"))
    },
    rows = function(st, acc, s, md) {
      # sign * delta + slack * scale >= 0: one design-weighted row, linear in beta (PHASE5a R10(a));
      # 'scale' turns delta into an average gap / slope in probability units (dial S9)
      f <- .rxc_assoc_form(st, s$type, s$version, s$row, s$col, s$other)
      sc <- .rxc_assoc_scale(st, s$type, s$version, s$row, s$other)
      .rxc_add_row(acc, s$sign * f$coef, s$sign * f$const + s$slack * sc, s$sign * f$a, "assoc", md)
    },
    label = function(s, rl, cl) .rxc_assoc_label(s, rl, cl),
    meta = .rxc_meta("belief (M-SIGN)", "mean association (model-implied, N-weighted)", "model",
                     "delta_W: yes; delta_B: only jointly with the rung (never alone with 2 rows)",
                     "rows", "sharp for the regression (PHASE5a R3(b))", "slack (probability units)",
                     "proved (PHASE5a R3, R10(a); 2c M2-M5)")),
  cr = list(
    union = TRUE,
    parse = function(s, ctx) {
      v <- .rxc_version_arg(s[["version"]], c("W1", "W2", "pair"), "cr")
      r <- ctx$ri(s[["row"]])
      list(row = r, col = ctx$ci(if (is.null(s[["col"]])) 1 else s[["col"]]), version = v,
           other = .rxc_other_arg(s[["other"]], r, v, "cr", ctx), slack = .rxc_slack_arg(s, "cr"))
    },
    label = function(s, rl, cl) .rxc_cr_label(s, rl, cl),
    meta = .rxc_meta("belief (M-SIGN, contextual reinforcement)", "mean associations (model-implied, N-weighted)",
                     "model", "only jointly with the rung (never alone with 2 rows)",
                     "W1 / 2x2: slab (one LP); W2 / pair: union of two LPs",
                     "sharp for the regression (PHASE5a R3(b, c))", "slack (probability units)",
                     "proved (PHASE5a R3(c, d), R10(b, c))"))
)

# parse user restrictions through the registry
.rxc_restrictions <- function(res, rl, cl, groups) {
  if (is.null(res)) return(NULL)
  ctx <- list(rl = rl, cl = cl,
              ri = function(v) .rxc_index(v, rl), ci = function(v) .rxc_index(v, cl),
              gi = function(v) if (is.null(v)) NULL else .rxc_index(as.character(v), groups))
  out <- list()
  for (nm in names(res)) {
    reg <- .rxc_registry[[nm]]
    if (is.null(reg)) { warning("unknown restriction type '", nm, "' ignored"); next }
    if (is.null(res[[nm]])) next
    out[[nm]] <- if (isTRUE(reg$scalar)) reg$parse(res[[nm]], ctx) else lapply(res[[nm]], reg$parse, ctx = ctx)
  }
  out
}

.rxc_module_ids <- function(ex) {
  ids <- character(0); types <- character(0); idx <- integer(0)
  for (nm in names(.rxc_registry)) {
    if (is.null(ex[[nm]])) next
    if (isTRUE(.rxc_registry[[nm]]$scalar)) { ids <- c(ids, nm); types <- c(types, nm); idx <- c(idx, NA) }
    else for (i in seq_along(ex[[nm]])) { ids <- c(ids, paste0(nm, "[", i, "]")); types <- c(types, nm); idx <- c(idx, i) }
  }
  data.frame(module = ids, type = types, index = idx, stringsAsFactors = FALSE)
}

# linear rows of every "rows" module (unions are handled by .rxc_cr_parts)
.rxc_extra_rows <- function(st, acc) {
  ex <- st$extra
  if (is.null(ex)) return(acc)
  ids <- .rxc_module_ids(ex)
  for (j in seq_len(nrow(ids))) {
    reg <- .rxc_registry[[ids$type[j]]]
    if (isTRUE(reg$union)) next
    spec <- if (is.na(ids$index[j])) ex[[ids$type[j]]] else ex[[ids$type[j]]][[ids$index[j]]]
    acc <- reg$rows(st, acc, spec, ids$module[j])
  }
  acc
}

.rxc_module_labels <- function(ex, rl, cl) {
  ids <- .rxc_module_ids(ex)
  out <- character(0)
  for (j in seq_len(nrow(ids))) {
    spec <- if (is.na(ids$index[j])) ex[[ids$type[j]]] else ex[[ids$type[j]]][[ids$index[j]]]
    out[ids$module[j]] <- .rxc_registry[[ids$type[j]]]$label(spec, rl, cl)
  }
  out
}

# active modules of a fit with their meta
.rxc_module_table <- function(ex, rl, cl) {
  ids <- .rxc_module_ids(ex)
  if (!nrow(ids)) return(NULL)
  lab <- .rxc_module_labels(ex, rl, cl)
  rows <- lapply(seq_len(nrow(ids)), function(j) {
    m <- .rxc_registry[[ids$type[j]]]$meta
    spec <- if (is.na(ids$index[j])) ex[[ids$type[j]]] else ex[[ids$type[j]]][[ids$index[j]]]
    sl <- if (is.list(spec) && !is.null(spec$slack)) spec$slack else NA_real_
    data.frame(module = ids$module[j], restriction = unname(lab[ids$module[j]]), slack = sl,
               as.data.frame(m, stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Registry of optional restrictions
#'
#' Lists the restriction types accepted by \code{eiBounds(restrictions = ...)}
#' with their properties: information source, the object restricted, the
#' target, whether aggregate data can refute them, the form (linear rows or a
#' union of LPs), sharpness, the sensitivity dial and the proof status
#' (see \code{docs/extension/PACKAGE_DESIGN.md}). The modules active in a fit
#' are in \code{fit$modules}.
#'
#' @return a data frame with one row per restriction type.
#' @export
eiModules <- function() {
  do.call(rbind, lapply(names(.rxc_registry), function(nm)
    data.frame(type = nm, as.data.frame(.rxc_registry[[nm]]$meta, stringsAsFactors = FALSE),
               stringsAsFactors = FALSE)))
}
