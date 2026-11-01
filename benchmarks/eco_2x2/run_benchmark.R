# Benchmark: eiDPM (blocked Gibbs in C++, and Stan HMC) against eco / ecoNP on
# the two eco example datasets with known unit truth: reg (n = 275) and
# census (n = 1040). Models: CAR and NCAR, parametric (K = 1) and DP mixture.
#
#   Rscript benchmarks/eco_2x2/run_benchmark.R [cores] [quick]
#
# Writes benchmarks/eco_2x2/results/{summary.csv, units.csv, report.md}.
# Pre-registered pass criteria (docs/bayes/PHASE1_proposal.md, sec. 5):
# aggregate posterior means within max(3 combined MCSE, 0.005) of eco, and unit
# posterior means correlated > 0.98 with eco's.

suppressMessages({ library(eco); library(eiDPM); library(posterior); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4
quick <- length(args) >= 2 && args[2] == "quick"
out_dir <- file.path("benchmarks", "eco_2x2", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data(reg, package = "eco"); data(census, package = "eco")
datasets <- list(reg = reg, census = census)

K_DP <- 30
sz <- if (quick) {
  list(eco = c(1500, 500, 1), gib = c(1500, 500, 1), stan = c(200, 200))
} else {
  list(eco = c(5000, 1000, 4), gib = c(5000, 1000, 4), stan = c(1000, 1000))
}
n_eco_chains <- 2; n_gib_chains <- 4; n_stan_chains <- 4
specs <- expand.grid(context = c(FALSE, TRUE), dp = c(FALSE, TRUE), data = names(datasets),
                     stringsAsFactors = FALSE)

# ---- one fit = list(W1 = draws x n, W2, chain, seconds = per-chain CPU seconds, extra)
fit_eco <- function(dd, context, dp, seed) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  fn <- if (dp) ecoNP else eco
  f <- fn(Y ~ X, N = N, data = dd, context = context,
          n.draws = sum(sz$eco[1:2]), burnin = sz$eco[2], thin = sz$eco[3] - 1)
  list(W1 = f$W[, 1, ], W2 = f$W[, 2, ], seconds = proc.time()[["elapsed"]] - t0)
}
run_eco <- function(dd, context, dp) {
  ch <- mclapply(seq_len(n_eco_chains), function(s) fit_eco(dd, context, dp, 100 + s),
                 mc.cores = min(cores, n_eco_chains))
  list(W1 = do.call(rbind, lapply(ch, `[[`, "W1")), W2 = do.call(rbind, lapply(ch, `[[`, "W2")),
       chain = rep(seq_along(ch), sapply(ch, function(c) nrow(c$W1))),
       seconds = sapply(ch, `[[`, "seconds"), extra = list())
}
run_dpm <- function(dd, context, dp, method) {
  K <- if (dp) K_DP else 1
  if (method == "gibbs") {
    ch <- mclapply(seq_len(n_gib_chains), function(s)
      eiDPM(dd$X, dd$Y, dd$N, context = context, K = K, method = "gibbs",
            n.draws = sz$gib[1], burnin = sz$gib[2], thin = sz$gib[3], seed = 200 + s),
      mc.cores = min(cores, n_gib_chains))
    f <- list(W1 = do.call(rbind, lapply(ch, `[[`, "W1")), W2 = do.call(rbind, lapply(ch, `[[`, "W2")),
              chain = rep(seq_along(ch), sapply(ch, function(c) nrow(c$W1))),
              seconds = sapply(ch, function(c) c$diagnostics$chain_seconds),
              extra = list(pi_last_q99 = max(sapply(ch, function(c) quantile(c$diagnostics$pi_last, 0.99))),
                           alpha_med = median(unlist(lapply(ch, function(c) c$diagnostics$alpha)))))
  } else {
    g <- eiDPM(dd$X, dd$Y, dd$N, context = context, K = K, method = "stan",
               n.draws = sz$stan[1], burnin = sz$stan[2], chains = n_stan_chains,
               cores = min(cores, n_stan_chains), seed = 300)
    f <- list(W1 = g$W1, W2 = g$W2, chain = g$chain, seconds = g$diagnostics$chain_seconds,
              extra = list(pi_last_q99 = unname(quantile(g$diagnostics$pi_last, 0.99)),
                           alpha_med = median(g$diagnostics$alpha),
                           divergences = g$diagnostics$divergences,
                           treedepth_hits = g$diagnostics$treedepth_max))
  }
  f
}

aggregates <- function(f, dd) {
  a1 <- dd$N * dd$X; a2 <- dd$N * (1 - dd$X)
  cbind(mean_W1 = rowMeans(f$W1), mean_W2 = rowMeans(f$W2),
        agg_W1 = as.vector(f$W1 %*% a1) / sum(a1), agg_W2 = as.vector(f$W2 %*% a2) / sum(a2))
}
truth_agg <- function(dd) {
  a1 <- dd$N * dd$X; a2 <- dd$N * (1 - dd$X)
  c(mean_W1 = mean(dd$W1), mean_W2 = mean(dd$W2),
    agg_W1 = sum(a1 * dd$W1) / sum(a1), agg_W2 = sum(a2 * dd$W2) / sum(a2))
}
dd_bounds <- function(dd) {
  L1 <- pmax(0, (dd$X + dd$Y - 1) / dd$X); U1 <- pmin(1, dd$Y / dd$X)
  L2 <- pmax(0, (dd$Y - dd$X) / (1 - dd$X)); U2 <- pmin(1, dd$Y / (1 - dd$X))
  a1 <- dd$N * dd$X; a2 <- dd$N * (1 - dd$X)
  b <- rbind(lower = c(mean(L1), mean(L2), sum(a1 * L1) / sum(a1), sum(a2 * L2) / sum(a2)),
             upper = c(mean(U1), mean(U2), sum(a1 * U1) / sum(a1), sum(a2 * U2) / sum(a2)))
  colnames(b) <- c("mean_W1", "mean_W2", "agg_W1", "agg_W2")
  b
}
chain_array <- function(v, chain) {
  m <- min(table(chain))
  sapply(sort(unique(chain)), function(ch) v[chain == ch][seq_len(m)])
}

summ <- list(); units <- list()
for (i in seq_len(nrow(specs))) {
  sp <- specs[i, ]; dd <- datasets[[sp$data]]
  label <- sprintf("%s %s %s", sp$data, if (sp$context) "NCAR" else "CAR", if (sp$dp) "DP" else "parametric")
  methods <- if (sp$dp) c("eco", "gibbs", "stan") else c("eco", "gibbs", "stan")
  fits <- list()
  for (m in methods) {
    cat(format(Sys.time()), label, m, "\n")
    fits[[m]] <- if (m == "eco") run_eco(dd, sp$context, sp$dp) else run_dpm(dd, sp$context, sp$dp, m)
  }
  tr <- truth_agg(dd); bd <- dd_bounds(dd)
  eco_agg <- aggregates(fits$eco, dd)
  eco_unit <- colMeans(fits$eco$W1)
  eco_mcse <- apply(eco_agg, 2, function(v) mcse_mean(chain_array(v, fits$eco$chain)))
  for (m in methods) {
    f <- fits[[m]]; ag <- aggregates(f, dd)
    for (j in colnames(ag)) {
      arr <- chain_array(ag[, j], f$chain)
      ess <- ess_bulk(arr); mc <- mcse_mean(arr)
      d_eco <- mean(ag[, j]) - mean(eco_agg[, j])
      summ[[length(summ) + 1]] <- data.frame(
        data = sp$data, assumption = if (sp$context) "NCAR" else "CAR",
        model = if (sp$dp) "DP" else "parametric", method = m, quantity = j,
        truth = tr[[j]], dd_lower = bd["lower", j], dd_upper = bd["upper", j],
        mean = mean(ag[, j]), sd = sd(ag[, j]),
        q05 = unname(quantile(ag[, j], 0.05)), q95 = unname(quantile(ag[, j], 0.95)),
        chain_mean_sd = sd(colMeans(arr)), rhat = rhat(arr), ess_bulk = ess, mcse = mc,
        cpu_seconds = sum(f$seconds), ess_per_cpu_s = ess / sum(f$seconds),
        diff_vs_eco = d_eco,
        z_vs_eco = if (m == "eco") NA else d_eco / sqrt(mc^2 + eco_mcse[[j]]^2),
        pi_last_q99 = if (is.null(f$extra$pi_last_q99)) NA else f$extra$pi_last_q99,
        alpha_med = if (is.null(f$extra$alpha_med)) NA else f$extra$alpha_med,
        divergences = if (is.null(f$extra$divergences)) NA else f$extra$divergences)
    }
    pm1 <- colMeans(f$W1); pm2 <- colMeans(f$W2)
    q1 <- apply(f$W1, 2, quantile, c(0.05, 0.95)); q2 <- apply(f$W2, 2, quantile, c(0.05, 0.95))
    units[[length(units) + 1]] <- data.frame(
      data = sp$data, assumption = if (sp$context) "NCAR" else "CAR",
      model = if (sp$dp) "DP" else "parametric", method = m,
      rmse_W1 = sqrt(mean((pm1 - dd$W1)^2)), rmse_W2 = sqrt(mean((pm2 - dd$W2)^2)),
      cover90_W1 = mean(dd$W1 >= q1[1, ] & dd$W1 <= q1[2, ]),
      cover90_W2 = mean(dd$W2 >= q2[1, ] & dd$W2 <= q2[2, ]),
      cor_eco_W1 = cor(pm1, eco_unit))
  }
  rm(fits); gc(verbose = FALSE)
  write.csv(do.call(rbind, summ), file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, units), file.path(out_dir, "units.csv"), row.names = FALSE)
}

S <- do.call(rbind, summ); U <- do.call(rbind, units)
S$pass_mean <- with(S, ifelse(method == "eco", NA, abs(z_vs_eco) <= 3 | abs(diff_vs_eco) <= 0.005))
U$pass_cor <- with(U, ifelse(method == "eco", NA, cor_eco_W1 > 0.98))
write.csv(S, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(U, file.path(out_dir, "units.csv"), row.names = FALSE)
cat("done\n")
