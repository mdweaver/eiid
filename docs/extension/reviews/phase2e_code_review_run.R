# R side of the oracle comparison for the code review of eiChange() / eiATT() (PHASE2e).
# Usage: Rscript phase2e_code_review_run.R INSTANCES.json RESULTS.json [lib]
# Reads the instances written by phase2e_code_review_oracle.py, runs the package on each (integer
# "auto" and "lp"; eiATT for causal instances) and writes endpoints, statuses and reports.
args <- commandArgs(TRUE)
if (length(args) >= 3) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({ library(eiPartialID); library(jsonlite) })
inst <- fromJSON(args[1], simplifyVector = FALSE)

rl <- function(k) paste0("g", k + 1); cl <- function(k) paste0("c", k + 1); el <- function(k) paste0("e", k + 1)
ul <- function(k) paste0("U", unlist(k) + 1)
mkdata <- function(d) {
  A <- do.call(rbind, lapply(d$A, unlist)); out <- list()
  for (t in seq_len(d$T)) {
    B <- do.call(rbind, lapply(d$B[[t]], unlist))
    df <- data.frame(u = paste0("U", seq_len(d$p)))
    for (r in seq_len(d$R)) df[[rl(r - 1)]] <- A[, r]
    for (c in seq_len(d$C)) df[[cl(c - 1)]] <- B[, c]
    out[[el(t - 1)]] <- df
  }
  out
}
fo <- function(d) as.formula(sprintf("cbind(%s) ~ cbind(%s)", paste(cl(seq_len(d$C) - 1), collapse = ","),
                                     paste(rl(seq_len(d$R) - 1), collapse = ",")))
vs <- function(v) if (identical(v, "rest")) "rest" else rl(v)
mkrow <- function(s) {
  tp <- s$type
  if (tp == "sign_change") return(list(type = s$stype, row = rl(s$row), col = cl(s$col), from = el(s$from), to = el(s$to),
                                       sign = s$sign, version = if (is.null(s$version)) "W1" else s$version,
                                       other = if (is.null(s$other)) NULL else rl(s$other), slack = if (is.null(s$slack)) 0 else s$slack))
  if (tp == "cr_change") return(list(row = rl(s$row), col = cl(s$col), from = el(s$from), to = el(s$to),
                                     version = if (is.null(s$version)) "W1" else s$version,
                                     other = if (is.null(s$other)) NULL else rl(s$other), slack = if (is.null(s$slack)) 0 else s$slack))
  if (tp == "order_change") return(list(row = rl(s$row), col = cl(s$col), from = el(s$from), to = el(s$to),
                                        over = if (is.null(s$over) || identical(s$over, "rest")) "rest" else rl(s$over),
                                        sign = if (is.null(s$sign)) 1 else s$sign, slack = if (is.null(s$slack)) 0 else s$slack))
  if (tp == "contrast_change") return(list(col = cl(s$col), from = el(s$from), to = el(s$to), low = ul(s$low), high = ul(s$high),
                                           beliefs = if (is.null(s$beliefs)) c("mix", "rate") else unlist(s$beliefs)))
  if (tp == "known") return(list(row = rl(s$row), col = cl(s$col), election = el(s$election), unit = ul(s$unit), value = s$value))
  if (tp == "total") return(list(row = rl(s$row), col = cl(s$col), election = el(s$election), units = ul(s$units),
                                 lower = s$lower, upper = s$upper))
  if (tp == "linear") {
    cells <- do.call(rbind, lapply(s$cells, function(z) data.frame(unit = paste0("U", z[[1]] + 1), row = rl(z[[2]]), col = cl(z[[3]]),
                                                                  election = el(z[[4]]), coef = z[[5]])))
    return(list(cells = cells, lower = s$lower, upper = s$upper))
  }
  stop("row type ", tp)
}
mkrows <- function(d) {
  if (!length(d$rows)) return(NULL)
  out <- list()
  for (s in d$rows) out[[length(out) + 1]] <- setNames(list(mkrow(s)), s$type)
  # group by type (eiChange takes a named list; several specs of one type as a list of lists)
  types <- unique(vapply(d$rows, function(s) s$type, ""))
  setNames(lapply(types, function(tp) lapply(Filter(function(s) s$type == tp, d$rows), mkrow)), types)
}
mktg <- function(tg) {
  switch(tg$kind,
         did = did(rl(tg$row), vs(tg$versus), cl(tg$col), el(tg$from), el(tg$to)),
         rate_at = rate_at(rl(tg$row), cl(tg$col), el(tg$election)),
         gap_at = gap_at(rl(tg$row), vs(tg$versus), cl(tg$col), el(tg$election)),
         bin_did = bin_did(cl(tg$col), el(tg$from), el(tg$to), low = ul(tg$low), high = ul(tg$high)))
}
mktrend <- function(tr) switch(tr$kind, none = no_trend(), dial = trend_dial(tr$M), rm = rel_mag(tr$Mbar), sd = smoothness(tr$M),
                               mtr = mtr_nd_ps(tr$mtr, tr$nd, tr$ps))
iv <- function(l, u) if (is.na(l)) NULL else c(l, u)

out <- list()
for (d in inst) {
  res <- tryCatch({
    dl <- mkdata(d); rows <- mkrows(d)
    r <- list(id = d$id)
    tgl <- if (length(d$targets)) lapply(d$targets, mktg) else NULL
    base <- el(d$base)
    if (length(d$targets)) {
      fa <- eiChange(dl, fo(d), unit = "u", base = base, targets = tgl, restrictions = rows, integer = "auto")
      fl <- eiChange(dl, fo(d), unit = "u", base = base, targets = tgl, restrictions = rows, integer = "lp")
      r$int <- lapply(seq_len(nrow(fa$table)), function(j) iv(fa$table$lower[j], fa$table$upper[j]))
      r$status <- as.list(fa$table$integer_status)
      r$lp <- lapply(seq_len(nrow(fl$table)), function(j) iv(fl$table$lower[j], fl$table$upper[j]))
      r$dd <- lapply(seq_len(nrow(fl$table)), function(j) iv(fl$table$dd_lower[j], fl$table$dd_upper[j]))
    }
    if (!is.null(d$causal)) {
      cs <- d$causal
      fc <- eiChange(dl, fo(d), unit = "u", base = base, row = rl(cs$treated), col = cl(cs$col), restrictions = rows, integer = "lp")
      args2 <- list(fc, treated = rl(cs$treated), population = "fixed", trend = mktrend(cs$trend), spillover = cs$Lambda,
                    no_anticipation = TRUE, rate_limits = isTRUE(cs$rate_limits))
      if (!is.null(cs$comparison)) args2$comparison <- vs(cs$comparison)
      a <- do.call(eiATT, args2)
      r$causal <- list()
      for (e in a$results) {
        t <- match(e$post, names(dl)) - 1
        r$causal[as.character(t)] <- list(iv(e$lower, e$upper))
      }
      e1 <- a$results[[1]]
      r$none <- iv(e1$none$lower, e1$none$upper); r$B_r <- iv(e1$B_r[1], e1$B_r[2])
      if (!is.null(e1$dial)) r$dial <- list(crit_upper = e1$dial$criterion_upper, crit_lower = e1$dial$criterion_lower)
      if (!is.null(e1$rm)) r$rm <- list(supA = e1$rm$supA, dial_set = e1$rm$dial_set,
                                       vac_upper = if (is.na(e1$rm$vacuous_upper)) NULL else e1$rm$vacuous_upper,
                                       vac_lower = if (is.na(e1$rm$vacuous_lower)) NULL else e1$rm$vacuous_lower,
                                       no_information = e1$rm$no_information, collapses = e1$rm$collapses_to_dial)
    }
    r
  }, error = function(e) list(id = d$id, error = conditionMessage(e)))
  out[[length(out) + 1]] <- res
}
writeLines(toJSON(out, auto_unbox = TRUE, digits = NA, null = "null"), args[2])
cat("ran", length(out), "instances;", sum(vapply(out, function(r) !is.null(r$error), TRUE)), "errors\n")
