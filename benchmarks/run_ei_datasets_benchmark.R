# Benchmark of eiPartialID (>= 0.2.0) on real data with known truth: ei.Datasets (Pavia 2023).
#
# New Zealand general elections 2002-2020: in each electorate, polling places are the units,
# the party vote gives the rows (Labour, National, other incl. informal), the electorate
# (candidate) vote gives the columns (Labour candidate, National candidate, other), and the
# official electorate-level party x candidate cross-tabulation is the truth. 492 tables
# (450 3x3, 42 3x2 where National fielded no candidate).
#
# Two parts:
#   1. one table per electorate (as before): rungs DD, M0, M0+order, M0 trim, L6;
#   2. one table per election, pooling the 3x3 electorates (polling places of all electorates
#      where National stood), fitted by M0 (pooled) and by M1 with electorates as groups; the
#      truth is the sum of the electorates' cross-tabulations. Also reported: the pooling test
#      (M0 against M1, Wald test of all electorate shifts = 0, P3 Prop. P3.6(c)) and the
#      emptiness test.
# Both parts report the CLR interval (algorithmic, Conjecture P3.5(b3)) and the projection
# interval from the same candidates (P3.5(b2)(ii)), each for the actual rate and intersected
# with Duncan-Davis.
#
# Usage (from the repository root, with eiPartialID and ei.Datasets installed):
#   Rscript benchmarks/run_ei_datasets_benchmark.R [cores] [ndraw]
# Writes benchmarks/results/ei_datasets_nz.csv.gz (one row per table x rung x cell)
# and benchmarks/results/ei_datasets_nz_summary.md.
suppressMessages(library(eiPartialID))
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else max(1, parallel::detectCores() - 1)
ndraw <- if (length(args) >= 2) as.integer(args[2]) else 100
out_dir <- Sys.getenv("BENCH_OUT", file.path("benchmarks", "results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tables <- eiDatasetsTables()
if (nzchar(Sys.getenv("BENCH_QUICK"))) tables <- tables[grepl("^NZ2002_0[1-5]", names(tables))]   # smoke test only
cat(sprintf("%d tables; cores = %d; ndraw = %d\n", length(tables), cores, ndraw))

# ordering restriction: each party's voters are at least as likely as the other party's voters
# to vote for that party's candidate (only where the candidate stood)
order_restr <- function(d) {
  od <- list()
  if ("c_LAB" %in% d$cols) od[[length(od) + 1]] <- c(high = "p_LAB", low = "p_NAT", col = "c_LAB")
  if ("c_NAT" %in% d$cols) od[[length(od) + 1]] <- c(high = "p_NAT", low = "p_LAB", col = "c_NAT")
  list(order = od)
}
rungs <- list(
  DD = list(model = "M0", ci = "none"),               # Duncan-Davis columns are read from this run
  M0 = list(model = "M0"),
  `M0+order` = list(model = "M0", restrictions = order_restr),
  `M0 trim` = list(model = "M0", range = "trim", trim = list(min_n = 100)),   # polling places < 100 voters
  L6 = list(model = "L6"))

# one fit -> one row per cell (as eiBenchmark, plus the projection and model-implied columns)
fit_rows <- function(d, rg, args, tests = FALSE) {
  truth <- d$truth / rowSums(d$truth)
  tv <- as.vector(t(truth))
  if (is.function(args$restrictions)) args$restrictions <- args$restrictions(d)
  args <- c(list(data = d$data, rows = d$rows, cols = d$cols, group = d$group, tests = tests,
                 seed = 2026, ndraw = ndraw), args)
  fit <- tryCatch(suppressWarnings(do.call(eiBounds, args)), error = function(e) e)
  if (inherits(fit, "error")) { warning(rg, ": ", conditionMessage(fit)); return(NULL) }
  tb <- fit$table
  tol <- 1e-9
  inside <- function(lo, hi) !is.na(lo) & tv >= lo - tol & tv <= hi + tol
  out <- data.frame(rung = rg, target = tb$target, truth = tv,
                    estimate_lower = tb$estimate_lower, estimate_upper = tb$estimate_upper, empty = tb$empty,
                    ci_model_lower = tb$ci_model_lower, ci_model_upper = tb$ci_model_upper,
                    ci_actual_lower = tb$ci_actual_lower, ci_actual_upper = tb$ci_actual_upper,
                    proj_model_lower = tb$proj_model_lower, proj_model_upper = tb$proj_model_upper,
                    proj_actual_lower = tb$proj_actual_lower, proj_actual_upper = tb$proj_actual_upper,
                    dd_lower = tb$dd_lower, dd_upper = tb$dd_upper, ci_rejected = tb$ci_rejected,
                    stringsAsFactors = FALSE)
  out$in_estimate <- !out$empty & inside(out$estimate_lower, out$estimate_upper)
  out$in_ci_actual <- inside(out$ci_actual_lower, out$ci_actual_upper)
  out$in_proj_actual <- inside(out$proj_actual_lower, out$proj_actual_upper)
  out$in_dd <- inside(out$dd_lower, out$dd_upper)
  out$width_estimate <- ifelse(out$empty, NA, out$estimate_upper - out$estimate_lower)
  out$width_ci_actual <- out$ci_actual_upper - out$ci_actual_lower
  out$width_proj_actual <- out$proj_actual_upper - out$proj_actual_lower
  out$width_dd <- out$dd_upper - out$dd_lower
  if (tests) {
    out$empty_test_rejects <- isTRUE(fit$tests$empty_rejected)
    out$pool_stat <- if (is.null(fit$tests$pool)) NA else fit$tests$pool$statistic
    out$pool_df <- if (is.null(fit$tests$pool)) NA else fit$tests$pool$df
    out$pool_p <- if (is.null(fit$tests$pool)) NA else fit$tests$pool$p.value
    out$m1_test_p <- if (is.null(fit$tests$m1)) NA else fit$tests$m1$p.value
  }
  out
}

# Set BENCH_SUMMARY_ONLY=1 to rebuild the summary from an existing CSV without refitting.
csv_path <- file.path(out_dir, "ei_datasets_nz.csv.gz")
if (!nzchar(Sys.getenv("BENCH_SUMMARY_ONLY"))) {
# ---------------- part 1: one table per electorate ----------------
t0 <- Sys.time()
one <- function(nm) {
  d <- tables[[nm]]
  res <- lapply(names(rungs), function(rg) fit_rows(d, rg, rungs[[rg]]))
  res <- do.call(rbind, res)
  if (is.null(res)) return(NULL)
  res <- cbind(dataset = nm, res)
  for (k in names(d$meta)) res[[k]] <- d$meta[[k]]
  res
}
parts <- parallel::mclapply(names(tables), one, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(parts, function(p) inherits(p, "try-error") || is.null(p), logical(1))
if (any(bad)) warning(sum(bad), " tables failed")
res <- do.call(rbind, parts[!bad])
elapsed1 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
# the DD rung reports Duncan-Davis bounds as its interval
dd <- res$rung == "DD"
res$ci_actual_lower[dd] <- res$proj_actual_lower[dd] <- res$dd_lower[dd]
res$ci_actual_upper[dd] <- res$proj_actual_upper[dd] <- res$dd_upper[dd]
res$in_ci_actual[dd] <- res$in_proj_actual[dd] <- res$in_dd[dd]
res$width_ci_actual[dd] <- res$width_proj_actual[dd] <- res$width_dd[dd]
res$estimate_lower[dd] <- res$estimate_upper[dd] <- NA; res$in_estimate[dd] <- NA; res$width_estimate[dd] <- NA
res$empty[dd] <- NA
res$diagonal <- sub("p_", "", sub(" ->.*", "", res$target)) == sub(".*-> c_", "", res$target)
res$part <- "electorate"

# ---------------- part 2: one table per election, electorates as groups ----------------
t0 <- Sys.time()
years <- sort(unique(vapply(tables, function(e) e$meta$election, 0)))
pooled_one <- function(yr) {
  tb3 <- tables[vapply(tables, function(e) e$meta$election == yr && length(e$cols) == 3, TRUE)]
  dat <- do.call(rbind, lapply(names(tb3), function(k) cbind(tb3[[k]]$data, electorate = k)))
  d <- list(data = dat, rows = tb3[[1]]$rows, cols = tb3[[1]]$cols, group = NULL,
            truth = Reduce(`+`, lapply(tb3, `[[`, "truth")))
  r0 <- fit_rows(d, "pooled M0", list(model = "M0"), tests = TRUE)
  d$group <- "electorate"
  r1 <- fit_rows(d, "M1 by electorate", list(model = "M1"), tests = TRUE)
  # the pooling test is computed with the M1 fit; copy it to the M0 rows for reference
  r0$pool_stat <- r1$pool_stat[1]; r0$pool_df <- r1$pool_df[1]; r0$pool_p <- r1$pool_p[1]
  out <- rbind(r0, r1)
  cbind(dataset = sprintf("NZ%d_pooled", yr), out, election = yr, district = "all (3x3 electorates)",
        units = nrow(dat), electorates = length(tb3))
}
pooled <- do.call(rbind, parallel::mclapply(years, pooled_one, mc.cores = min(cores, length(years)),
                                            mc.preschedule = FALSE))
elapsed2 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
pooled$diagonal <- sub("p_", "", sub(" ->.*", "", pooled$target)) == sub(".*-> c_", "", pooled$target)
pooled$part <- "election"

common <- intersect(names(res), names(pooled))
allres <- rbind(cbind(res[, common], electorates = 1, empty_test_rejects = NA, pool_stat = NA, pool_df = NA,
                      pool_p = NA, m1_test_p = NA),
                pooled[, c(common, "electorates", "empty_test_rejects", "pool_stat", "pool_df", "pool_p", "m1_test_p")])
allres$run_minutes <- ifelse(allres$part == "electorate", elapsed1, elapsed2)
allres$run_cores <- cores; allres$run_ndraw <- ndraw; allres$run_date <- format(Sys.Date())
utils::write.csv(allres, gzfile(csv_path), row.names = FALSE)
}
allres <- utils::read.csv(gzfile(csv_path), stringsAsFactors = FALSE)
res <- allres[allres$part == "electorate", ]
pooled <- allres[allres$part == "election", ]
elapsed1 <- res$run_minutes[1]; elapsed2 <- pooled$run_minutes[1]
cores <- allres$run_cores[1]; ndraw <- allres$run_ndraw[1]

# ---------------- summary ----------------
rung_order <- names(rungs)
summ <- function(df) {
  data.frame(cells = nrow(df),
             ci_covers = mean(df$in_ci_actual, na.rm = TRUE),
             ci_rejected = mean(df$ci_rejected | is.na(df$ci_actual_lower), na.rm = TRUE),
             ci_misses = mean(!is.na(df$ci_actual_lower) & !df$in_ci_actual, na.rm = TRUE),
             median_ci_width = stats::median(df$width_ci_actual, na.rm = TRUE),
             proj_covers = mean(df$in_proj_actual, na.rm = TRUE),
             proj_misses = mean(!is.na(df$proj_actual_lower) & !df$in_proj_actual, na.rm = TRUE),
             median_proj_width = stats::median(df$width_proj_actual, na.rm = TRUE),
             ci_width_vs_dd = stats::median(df$width_ci_actual / pmax(df$width_dd, 1e-12), na.rm = TRUE),
             proj_width_vs_dd = stats::median(df$width_proj_actual / pmax(df$width_dd, 1e-12), na.rm = TRUE),
             estimate_empty = if (all(is.na(df$empty))) NA else mean(df$empty, na.rm = TRUE),
             estimate_covers = mean(df$in_estimate, na.rm = TRUE))
}
tab <- do.call(rbind, lapply(rung_order, function(r) cbind(rung = r, summ(res[res$rung == r, ]))))
tab_diag <- do.call(rbind, lapply(rung_order, function(r) cbind(rung = r, summ(res[res$rung == r & res$diagonal, ]))))
tab_year <- do.call(rbind, lapply(sort(unique(res$election)), function(y)
  do.call(rbind, lapply(c("M0", "M0 trim"), function(r) cbind(election = y, rung = r, summ(res[res$rung == r & res$election == y, ]))))))
p_rungs <- c("pooled M0", "M1 by electorate")
tab_pool <- do.call(rbind, lapply(p_rungs, function(r) cbind(rung = r, summ(pooled[pooled$rung == r, ]))))
mark <- function(inside, lo) if (is.na(lo)) " rejected" else if (isTRUE(inside)) "" else " **miss**"
f3 <- function(v) ifelse(is.na(v), "--", sprintf("%.3f", v))
cell_lines <- c("| election | electorates | units | rung | cell | truth | DD | CLR CI (actual) | projection CI (actual) |",
                "|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(pooled))) {
  p <- pooled[i, ]
  cell_lines <- c(cell_lines, sprintf("| %d | %d | %d | %s | %s | %.3f | [%s, %s] | [%s, %s]%s | [%s, %s]%s |",
    p$election, p$electorates, p$units, p$rung, p$target, p$truth, f3(p$dd_lower), f3(p$dd_upper),
    f3(p$ci_actual_lower), f3(p$ci_actual_upper), mark(p$in_ci_actual, p$ci_actual_lower),
    f3(p$proj_actual_lower), f3(p$proj_actual_upper), mark(p$in_proj_actual, p$proj_actual_lower)))
}
tests_tab <- unique(pooled[pooled$rung == "M1 by electorate",
                           c("election", "electorates", "units", "pool_stat", "pool_df", "pool_p", "m1_test_p")])
e0 <- unique(pooled[pooled$rung == "pooled M0", c("election", "empty_test_rejects")])
e1 <- unique(pooled[pooled$rung == "M1 by electorate", c("election", "empty_test_rejects")])
tests_tab$M0_estimate_empty <- vapply(tests_tab$election, function(y) all(pooled$empty[pooled$rung == "pooled M0" & pooled$election == y]), TRUE)
tests_tab$M0_emptiness_test_rejects <- e0$empty_test_rejects[match(tests_tab$election, e0$election)]
tests_tab$M1_estimate_empty <- vapply(tests_tab$election, function(y) all(pooled$empty[pooled$rung == "M1 by electorate" & pooled$election == y]), TRUE)
tests_tab$M1_emptiness_test_rejects <- e1$empty_test_rejects[match(tests_tab$election, e1$election)]

fmt <- function(df) {
  num <- vapply(df, is.numeric, logical(1))
  # counts (whole numbers, some above 1) as integers; everything else to 3 decimals
  df[num] <- lapply(df[num], function(v) {
    cnt <- all(is.na(v) | abs(v - round(v)) < 1e-12) && any(abs(v) > 1, na.rm = TRUE)
    ifelse(is.na(v), "--", if (cnt) sprintf("%.0f", v) else ifelse(abs(v) >= 1000, sprintf("%.0f", v), sprintf("%.3f", v)))
  })
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}
md <- c(
  "# Benchmark: ei.Datasets, New Zealand general elections 2002-2020",
  "",
  sprintf("Generated by `benchmarks/run_ei_datasets_benchmark.R` with eiPartialID %s and ei.Datasets %s on %s (part 1: %.1f minutes, part 2: %.1f minutes; %d cores, ndraw = %d).",
          utils::packageVersion("eiPartialID"), utils::packageVersion("ei.Datasets"), allres$run_date[1], elapsed1, elapsed2, cores, ndraw),
  "",
  sprintf("%d electorate tables: units are polling places (median %d per table); rows are the party vote (Labour, National, other incl. informal); columns are the electorate vote (Labour candidate, National candidate, other); the truth is the official party x candidate cross-tabulation, which aggregates exactly the same voters as the units. Targets are all row-by-column rates (transition proportions).",
          length(unique(res$dataset)), as.integer(stats::median(res$units[!duplicated(res$dataset)]))),
  "",
  "Intervals are 95% confidence intervals for the actual rate, intersected with Duncan-Davis: `ci` is the CLR interval (default; algorithmic, P3 Conjecture P3.5(b3)), `proj` the projection interval from the same candidates (P3 Prop. P3.5(b2)(ii)). Both use the Bhatia-Davis allowance. `estimate` is the plug-in estimate of the identified set (not a confidence interval). `*_width_vs_dd` is the median ratio of the CI width to the Duncan-Davis width. `ci_rejected`: the CLR ends cross or lie outside the logical range (the rung is contradicted); `ci_misses`/`proj_misses`: share of cells where a non-empty interval excludes the truth.",
  "",
  "Rungs: `DD` Duncan-Davis (no model); `M0` pooled linear contextual effects; `M0+order` M0 plus ordering (each party's voters are at least as likely as the other party's voters to vote for that party's candidate); `M0 trim` M0 with polling places of fewer than 100 voters removed from the composition range (their votes enter through Duncan-Davis); `L6` Goodman constancy.",
  "",
  "## Part 1: one table per electorate", "",
  "### All cells", "", fmt(tab), "",
  "### Diagonal cells (same party in both votes)", "", fmt(tab_diag), "",
  "### By election (M0 and M0 trim, all cells)", "", fmt(tab_year), "",
  "## Part 2: one table per election, electorates as groups", "",
  "Each election's polling places in the electorates where National stood (3x3 tables) form one table. `pooled M0` ignores the electorates; `M1 by electorate` adds electorate shifts (groups = electorates). The truth is the sum of the electorates' official cross-tabulations. The pooling test is the Wald test of all electorate shifts = 0 in the M1 regression (P3 Prop. P3.6(c)); the M1 test is M1 against M2 (P3.6(b)).", "",
  "### Tests", "", fmt(tests_tab), "",
  "### Coverage summary (7 elections x 9 cells)", "", fmt(tab_pool), "",
  "### Every cell", "", cell_lines, "")
writeLines(md, file.path(out_dir, "ei_datasets_nz_summary.md"))
cat(md, sep = "\n")
