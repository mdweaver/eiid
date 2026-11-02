# Per-sweep time of the corrected parametric eco (eco_param_cpp) on reg and census,
# against eco() as distributed (2 chains x 6000 sweeps in run_benchmark.R).
#   Rscript benchmarks/eco_2x2/eco_fixed/param_timing.R
suppressMessages({ library(eiDPM); library(eco) })
Rcpp::sourceCpp("benchmarks/eco_2x2/eco_fixed/eco_alg1.cpp")
data(reg); data(census)
old <- read.csv("benchmarks/eco_2x2/results/summary.csv")
for (dn in c("reg", "census")) for (ctx in c(FALSE, TRUE)) {
  dd <- get(dn); d <- eiDPM:::dpmData(dd$X, dd$Y, ctx); pr <- eiDPM:::.dpm_prior(list(), d$p)
  for (ws in 0:1) {
    set.seed(1); t0 <- proc.time()[["elapsed"]]
    r <- eco_param_cpp(d$p, d$X, d$Y, d$L, d$U, d$Wobs, d$w1obs, d$w2obs, pr$mu0, pr$tau0, pr$nu0, pr$S0,
                       22000L, 2000L, 10L, (d$L + d$U) / 2, as.integer(ws))
    s <- proc.time()[["elapsed"]] - t0
    a1 <- dd$N * dd$X                                  # aggregate W1: weights N * X, all units
    agg <- mean((r$W1 %*% a1[d$lat] + sum(a1[d$yb] * d$Wyb)) / sum(a1[c(d$lat, d$yb)]))
    o <- old[old$data == dn & old$assumption == (if (ctx) "NCAR" else "CAR") & old$model == "parametric" &
             old$method == "eco" & old$quantity == "agg_W1", ]
    cat(sprintf("%-7s %-5s %-6s corrected %.4f ms/sweep (eco original %.4f) | agg_W1 %.4f (eco original %.4f)\n",
                dn, if (ctx) "NCAR" else "CAR", if (ws) "slice" else "MH", 1000 * s / 22000,
                1000 * o$cpu_seconds / 12000, agg, o$mean))
  }
}
