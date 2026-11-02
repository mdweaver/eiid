# Calibration of the eiDPM exact-DP sampler: simulation-based calibration with
# equivalence criteria, as pre-registered in docs/bayes/calibration_plan.md.
#
#   Rscript benchmarks/eco_2x2/calibration.R [cores] [scale]
#   scale (default 1) multiplies the replication counts (use 0.01 for a smoke test).
#   Rscript benchmarks/eco_2x2/calibration.R analyse      # re-run only the analysis
#
# Outputs in benchmarks/eco_2x2/results/:
#   calibration_raw.csv     one row per replication x quantity (rank, posterior summaries, ESS)
#   calibration_metrics.csv one row per design x quantity x metric (estimate, CI, region, verdict)
#   calibration_report.md   summary tables
#   calibration_ecdf_<design>.png  ECDF-difference plots with simultaneous 95% bands
suppressMessages({ library(eiDPM); library(parallel); library(posterior) })
args <- commandArgs(trailingOnly = TRUE)
analyse_only <- length(args) >= 1 && args[1] == "analyse"
cores <- if (!analyse_only && length(args) >= 1) as.integer(args[1]) else 4
scale <- if (!analyse_only && length(args) >= 2) as.numeric(args[2]) else 1
out_dir <- Sys.getenv("CALIB_OUT", file.path("benchmarks", "eco_2x2", "results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_units <- 100; L <- 199; thin <- 25; burn <- 1000
prior <- list(tau0 = 2, nu0 = 4, S0 = 10, a0 = 1, b0 = 0.1)
designs <- list(
  D1  = list(gen = "dp",   p = 2, K = Inf, R = 4000, seed0 = 1e5),
  D1b = list(gen = "dp",   p = 2, K = 50,  R = 2000, seed0 = 1e5),   # same simulations as D1
  D2  = list(gen = "par",  p = 2, K = 1,   R = 2000, seed0 = 2e5),
  D3  = list(gen = "dp",   p = 3, K = Inf, R = 2000, seed0 = 3e5))

# ---------------------------------------------------------------- simulation
simulate <- function(seed, gen, p) {
  set.seed(seed)
  S0 <- diag(prior$S0, p); attempts <- 0
  repeat {
    attempts <- attempts + 1
    alpha <- rgamma(1, prior$a0, prior$b0)
    if (gen == "dp") {                          # Chinese restaurant process
      z <- integer(n_units); z[1] <- 1; k <- 1
      for (i in 2:n_units) {
        pr <- c(tabulate(z[1:(i - 1)], k), alpha)
        z[i] <- sample.int(k + 1, 1, prob = pr)
        if (z[i] > k) k <- k + 1
      }
    } else { z <- rep(1L, n_units); k <- 1 }
    Ws <- matrix(NA, n_units, p)
    for (j in seq_len(k)) {
      Sig <- solve(rWishart(1, prior$nu0, solve(S0))[, , 1])
      mu <- MASS::mvrnorm(1, rep(0, p), Sig / prior$tau0)
      Ws[z == j, ] <- MASS::mvrnorm(sum(z == j), mu, Sig)
    }
    W <- plogis(Ws[, 1:2, drop = FALSE])
    X <- if (p == 3) plogis(Ws[, 3]) else runif(n_units, 0.05, 0.95)
    Y <- X * W[, 1] + (1 - X) * W[, 2]
    ok <- all(Y > 1e-6 & Y < 1 - 1e-6) && (p == 2 || all(X > 1e-6 & X < 1 - 1e-6))
    if (ok) break                               # data-only rejection
  }
  list(X = X, Y = Y, W = W, alpha = alpha, clusters = length(unique(z)), attempts = attempts)
}

quantities <- function(W1, W2) cbind(
  W1_u1 = W1[, 1], W1_u2 = W1[, 2], W1_u3 = W1[, 3],
  mean_W1 = rowMeans(W1), mean_W2 = rowMeans(W2),
  sd_W1 = apply(W1, 1, sd), cor_W1W2 = sapply(seq_len(nrow(W1)), function(s) cor(W1[s, ], W2[s, ])))

one_fit <- function(r, dname) {
  ds <- designs[[dname]]
  s <- simulate(ds$seed0 + r, ds$gen, ds$p)
  t0 <- proc.time()[["elapsed"]]
  f <- eiDPM(s$X, s$Y, context = ds$p == 3, K = ds$K, n.draws = L * thin, burnin = burn,
             thin = thin, prior = prior, seed = r)
  secs <- proc.time()[["elapsed"]] - t0
  dr <- quantities(f$W1, f$W2)
  tv <- c(W1_u1 = s$W[1, 1], W1_u2 = s$W[2, 1], W1_u3 = s$W[3, 1],
          mean_W1 = mean(s$W[, 1]), mean_W2 = mean(s$W[, 2]),
          sd_W1 = sd(s$W[, 1]), cor_W1W2 = cor(s$W[, 1], s$W[, 2]))
  if (ds$K != 1) {
    dr <- cbind(dr, alpha = f$diagnostics$alpha, clusters = f$diagnostics$n_clusters)
    tv <- c(tv, alpha = s$alpha, clusters = s$clusters)
  }
  ddw <- c(pmin(1, s$Y[1:3] / s$X[1:3]) - pmax(0, (s$X[1:3] + s$Y[1:3] - 1) / s$X[1:3]),
           rep(NA, length(tv) - 3))
  do.call(rbind, lapply(seq_along(tv), function(j) {
    d <- dr[, j]; t <- tv[[j]]
    lo <- sum(d < t); eq <- sum(d == t)
    rank <- lo + if (eq > 0) sample.int(eq + 1, 1) - 1 else 0
    data.frame(design = dname, rep = r, quantity = names(tv)[j], rank = rank, truth = t,
               post_mean = mean(d), post_sd = sd(d),
               width90 = unname(diff(quantile(d, c(0.05, 0.95)))), dd_width = ddw[j],
               ess = suppressWarnings(ess_bulk(matrix(d, ncol = 1))),
               seconds = secs, attempts = s$attempts)
  }))
}

raw_file <- file.path(out_dir, "calibration_raw.csv")
if (!analyse_only) {
  raw <- NULL
  for (dname in names(designs)) {
    R <- max(4, round(designs[[dname]]$R * scale))
    cat(format(Sys.time()), dname, R, "replications\n")
    rs <- mclapply(seq_len(R), function(r) tryCatch(one_fit(r, dname), error = function(e) {
      message(dname, " rep ", r, ": ", conditionMessage(e)); NULL }),
      mc.cores = cores, mc.preschedule = FALSE)
    cat("  failed:", sum(sapply(rs, is.null)), "\n")
    raw <- rbind(raw, do.call(rbind, rs))
    write.csv(raw, raw_file, row.names = FALSE)
  }
}
raw <- read.csv(raw_file)

# ---------------------------------------------------------------- analysis
set.seed(20261102)
nrank <- L + 1
var_unif <- (nrank^2 - 1) / (12 * nrank^2)
cp <- function(k, n) binom.test(k, n, conf.level = 0.90)$conf.int
verdict <- function(lo, hi, a, b) if (lo >= a && hi <= b) "equivalent" else if (hi < a || lo > b) "miscalibrated" else "inconclusive"
ksD <- function(u) { u <- sort(u); n <- length(u); max(pmax(seq_len(n) / n - u, u - (seq_len(n) - 1) / n)) }
met <- NULL
for (dname in unique(raw$design)) for (q in unique(raw$quantity[raw$design == dname])) {
  x <- raw[raw$design == dname & raw$quantity == q, ]
  rk <- x$rank; R <- length(rk); u <- (rk + 0.5) / nrank
  add <- function(metric, est, lo, hi, a, b, v = verdict(lo, hi, a, b))
    met <<- rbind(met, data.frame(design = dname, quantity = q, R = R, metric = metric,
                                  estimate = est, ci_lo = lo, ci_hi = hi, region_lo = a, region_hi = b, verdict = v))
  for (lev in c(0.5, 0.8, 0.9, 0.95)) {
    lo_r <- round((1 - lev) / 2 * nrank); k <- sum(rk >= lo_r & rk <= nrank - lo_r - 1)
    ci <- cp(k, R); add(sprintf("coverage_%02d", round(100 * lev)), k / R, ci[1], ci[2], lev - 0.03, lev + 0.03)
  }
  kl <- sum(rk < 10); ci <- cp(kl, R); add("left_miss_90", kl / R, ci[1], ci[2], 0.03, 0.07)
  kr <- sum(rk >= nrank - 10); ci <- cp(kr, R); add("right_miss_90", kr / R, ci[1], ci[2], 0.03, 0.07)
  m <- mean(u); se <- sd(u) / sqrt(R); add("mean_rank", m, m - 1.645 * se, m + 1.645 * se, 0.47, 0.53)
  bv <- replicate(2000, var(sample(u, replace = TRUE))) / var_unif
  add("var_ratio", var(u) / var_unif, quantile(bv, 0.05), quantile(bv, 0.95), 0.90, 1.10)
  bd <- replicate(2000, ksD(sample(u, replace = TRUE)))
  D <- ksD(u); up <- unname(quantile(bd, 0.95)); low <- unname(quantile(bd, 0.05))
  add("ks_D", D, low, up, 0, 0.05, if (up < 0.05) "equivalent" else if (low > 0.05) "miscalibrated" else "inconclusive")
}
write.csv(met, file.path(out_dir, "calibration_metrics.csv"), row.names = FALSE)

# per quantity and design verdicts
qv <- aggregate(verdict ~ design + quantity, met, function(v)
  if (any(v == "miscalibrated")) "miscalibrated" else if (all(v == "equivalent")) "equivalent" else "inconclusive")
names(qv)[3] <- "overall"
cov90 <- met[met$metric == "coverage_90", c("design", "quantity", "verdict")]; names(cov90)[3] <- "cov90"
qv <- merge(qv, cov90)
dv <- do.call(rbind, lapply(split(qv, qv$design), function(d) data.frame(
  design = d$design[1],
  calibrated = all(d$cov90 == "equivalent") && !any(d$overall == "miscalibrated"))))

# diagnostics
diagn <- do.call(rbind, lapply(split(raw, list(raw$design, raw$quantity), drop = TRUE), function(x) {
  z <- (x$post_mean - x$truth) / x$post_sd
  data.frame(design = x$design[1], quantity = x$quantity[1],
             ess_median = median(x$ess), ess_share_below_100 = mean(x$ess < 100),
             z_mean = mean(z[is.finite(z)]), z_sd = sd(z[is.finite(z)]),
             contraction = 1 - mean(x$post_sd^2) / var(x$truth),
             width90_over_dd = if (all(is.na(x$dd_width))) NA else median(x$width90 / x$dd_width, na.rm = TRUE))
}))
timing <- aggregate(seconds ~ design, unique(raw[, c("design", "rep", "seconds")]), median)
attempts <- aggregate(attempts ~ design, unique(raw[, c("design", "rep", "attempts")]), mean)

# ECDF-difference plots with simultaneous bands (Sailynoja, Burkner & Vehtari 2022)
band_gamma <- function(R, zg, M = 3000) {
  mins <- replicate(M, {
    u <- (sample.int(nrank, R, replace = TRUE) - 0.5) / nrank
    N <- sapply(zg, function(z) sum(u <= z))
    min(2 * pmin(pbinom(N, R, zg), 1 - pbinom(N - 1, R, zg)))
  })
  unname(quantile(mins, 0.05))
}
zg <- (1:(nrank - 1)) / nrank
for (dname in unique(raw$design)) {
  qs <- unique(raw$quantity[raw$design == dname]); R <- sum(raw$design == dname & raw$quantity == qs[1])
  g <- band_gamma(R, zg)
  lo <- qbinom(g / 2, R, zg) / R - zg; hi <- qbinom(1 - g / 2, R, zg) / R - zg
  png(file.path(out_dir, sprintf("calibration_ecdf_%s.png", dname)), width = 1500, height = 900, res = 130)
  op <- par(mfrow = c(ceiling(length(qs) / 3), 3), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
  for (q in qs) {
    u <- (raw$rank[raw$design == dname & raw$quantity == q] + 0.5) / nrank
    dd <- ecdf(u)(zg) - zg
    plot(zg, dd, type = "n", ylim = range(c(lo, hi, dd)), xlab = "fractional rank", ylab = "ECDF - uniform",
         main = sprintf("%s: %s", dname, q))
    polygon(c(zg, rev(zg)), c(lo, rev(hi)), col = "grey85", border = NA)
    abline(h = 0, col = "grey50"); lines(zg, dd, lwd = 2, col = ifelse(any(dd < lo | dd > hi), "firebrick", "black"))
  }
  par(op); dev.off()
}

# report
f3 <- function(x) formatC(x, format = "f", digits = 3)
lines <- c("# Calibration of the eiDPM exact-DP sampler", "",
           "Pre-registered in `docs/bayes/calibration_plan.md`. Generated by `calibration.R` from `calibration_raw.csv`.", "",
           "## Design verdicts (primary criterion)", "",
           "| design | calibrated | replications | median seconds per fit | mean data draws per replication |", "|---|---|---|---|---|")
for (d in dv$design) lines <- c(lines, sprintf("| %s | %s | %d | %s | %s |", d, dv$calibrated[dv$design == d],
  max(met$R[met$design == d]), f3(timing$seconds[timing$design == d]), f3(attempts$attempts[attempts$design == d])))
lines <- c(lines, "", "## Verdicts by quantity", "", "| design | quantity | 90% coverage | overall |", "|---|---|---|---|")
for (k in seq_len(nrow(qv))) lines <- c(lines, sprintf("| %s | %s | %s | %s |", qv$design[k], qv$quantity[k], qv$cov90[k], qv$overall[k]))
lines <- c(lines, "", "## Metrics", "", "| design | quantity | metric | estimate | 90% CI | region | verdict |", "|---|---|---|---|---|---|---|")
for (k in seq_len(nrow(met))) lines <- c(lines, sprintf("| %s | %s | %s | %s | [%s, %s] | [%s, %s] | %s |",
  met$design[k], met$quantity[k], met$metric[k], f3(met$estimate[k]), f3(met$ci_lo[k]), f3(met$ci_hi[k]),
  f3(met$region_lo[k]), f3(met$region_hi[k]), met$verdict[k]))
lines <- c(lines, "", "## Diagnostics (no pass/fail)", "",
           "| design | quantity | median ESS (of 199) | share ESS < 100 | mean z | sd z | contraction | 90% width / DD width |",
           "|---|---|---|---|---|---|---|---|")
for (k in seq_len(nrow(diagn))) lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
  diagn$design[k], diagn$quantity[k], f3(diagn$ess_median[k]), f3(diagn$ess_share_below_100[k]), f3(diagn$z_mean[k]),
  f3(diagn$z_sd[k]), f3(diagn$contraction[k]), ifelse(is.na(diagn$width90_over_dd[k]), "", f3(diagn$width90_over_dd[k]))))
writeLines(lines, file.path(out_dir, "calibration_report.md"))
print(dv); print(qv)
