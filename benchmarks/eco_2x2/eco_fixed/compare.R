# Corrected eco (efficient code) versus the eiDPM collapsed sampler, CAR and NCAR
# DP mixtures on reg and census. Separates the sources of the speed difference:
#   eco_units   : eco's algorithm as eco codes it (loop over units, n^2), MH rate step
#   eco_grouped : the same with the configuration sum taken over clusters (n k)
#   eco_slice   : eco_grouped with the slice step for the rates instead of MH
#   collapsed   : eiDPM, Neal (2000) Algorithm 3, slice step
# The original eco timings come from results/summary.csv (run_benchmark.R).
#   Rscript benchmarks/eco_2x2/eco_fixed/compare.R [cores]
suppressMessages({ library(eiDPM); library(parallel); library(posterior) })
Rcpp::sourceCpp("benchmarks/eco_2x2/eco_fixed/eco_alg1.cpp")
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[1]) else 4
out <- "benchmarks/eco_2x2/eco_fixed"
data(reg, package = "eco"); data(census, package = "eco")
datasets <- list(reg = reg, census = census)
burn <- 2000; keep <- 2000; thin <- 10; chains <- 2
samplers <- c("eco_units", "eco_grouped", "eco_slice", "collapsed")

fit <- function(dn, context, sampler, chain) {
  dd <- datasets[[dn]]
  d <- eiDPM:::dpmData(dd$X, dd$Y, context); pr <- eiDPM:::.dpm_prior(list(), d$p)
  set.seed(1000 * chain + match(sampler, samplers) + 10 * context)
  args <- list(d$p, d$X, d$Y, d$L, d$U, d$Wobs, d$w1obs, d$w2obs, pr$mu0, pr$tau0, pr$nu0, pr$S0,
               pr$a0, pr$b0, as.integer(burn + keep * thin), as.integer(burn), as.integer(thin),
               (d$L + d$U) / 2)
  t0 <- proc.time()[["elapsed"]]
  r <- switch(sampler,
    eco_units = do.call(eco_alg1_cpp, c(args, list(0L, 0L, 1.0))),
    eco_grouped = do.call(eco_alg1_cpp, c(args, list(0L, 1L, 1.0))),
    eco_slice = do.call(eco_alg1_cpp, c(args, list(1L, 1L, 1.0))),
    collapsed = do.call(eiDPM:::dpm_collapsed_cpp, c(args, list(1L, 1.0))))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$N <- dd$N[d$lat]
  r
}

jobs <- expand.grid(chain = seq_len(chains), sampler = samplers, context = c(FALSE, TRUE),
                    data = names(datasets), stringsAsFactors = FALSE)
jobs <- jobs[order(jobs$data == "reg", jobs$sampler != "eco_units"), ]   # long jobs first
cat(format(Sys.time()), nrow(jobs), "fits\n")
res <- mclapply(seq_len(nrow(jobs)), function(j) with(jobs[j, ], fit(data, context, sampler, chain)),
                mc.cores = cores, mc.preschedule = FALSE)
cat(format(Sys.time()), "done\n")

ess2 <- function(x) ess_bulk(matrix(x, ncol = chains))   # draws x chains
rows <- list(); unitmeans <- list()
for (dn in names(datasets)) for (ctx in c(FALSE, TRUE)) for (s in samplers) {
  ix <- which(jobs$data == dn & jobs$context == ctx & jobs$sampler == s)
  rr <- res[ix]
  W1 <- do.call(rbind, lapply(rr, `[[`, "W1")); W2 <- do.call(rbind, lapply(rr, `[[`, "W2"))
  N <- rr[[1]]$N
  agg <- as.vector(W1 %*% N / sum(N))
  alpha <- unlist(lapply(rr, `[[`, "alpha")); k <- unlist(lapply(rr, `[[`, "n_clusters"))
  unit_ess <- apply(W1, 2, ess2)
  secs <- sum(sapply(rr, `[[`, "seconds")); sweeps <- chains * (burn + keep * thin)
  e <- c(agg_W1 = ess2(agg), alpha = ess2(alpha), clusters = ess2(k),
         unit_median = median(unit_ess), unit_p10 = unname(quantile(unit_ess, 0.1)))
  rows[[length(rows) + 1]] <- data.frame(
    data = dn, model = if (ctx) "NCAR" else "CAR", sampler = s,
    ms_per_sweep = 1000 * secs / sweeps, seconds = secs,
    mean_agg_W1 = mean(agg), mean_alpha = mean(alpha), mean_clusters = mean(k),
    rhat_alpha = rhat(matrix(alpha, ncol = chains)), rhat_agg = rhat(matrix(agg, ncol = chains)),
    t(e), t(setNames(e / secs, paste0(names(e), "_per_s"))))
  unitmeans[[paste(dn, ctx, s)]] <- colMeans(W1)
}
tab <- do.call(rbind, rows)
# agreement of unit posterior means with the collapsed sampler
tab$unit_cor_vs_collapsed <- NA; tab$unit_rmsd_vs_collapsed <- NA
for (i in seq_len(nrow(tab))) {
  ref <- unitmeans[[paste(tab$data[i], tab$model[i] == "NCAR", "collapsed")]]
  u <- unitmeans[[paste(tab$data[i], tab$model[i] == "NCAR", tab$sampler[i])]]
  tab$unit_cor_vs_collapsed[i] <- cor(u, ref); tab$unit_rmsd_vs_collapsed[i] <- sqrt(mean((u - ref)^2))
}
# original eco, from the earlier benchmark (2 chains x 6000 sweeps; ESS of the unweighted mean W1)
old <- read.csv("benchmarks/eco_2x2/results/summary.csv")
old <- old[old$model == "DP" & old$method == "eco" & old$quantity == "agg_W1", ]
tab_old <- data.frame(data = old$data, model = old$assumption, sampler = "eco_original",
                      ms_per_sweep = 1000 * old$cpu_seconds / 12000, seconds = old$cpu_seconds,
                      mean_agg_W1 = old$mean, agg_W1 = old$ess_bulk, agg_W1_per_s = old$ess_per_cpu_s)
write.csv(tab, file.path(out, "compare_results.csv"), row.names = FALSE)
write.csv(tab_old, file.path(out, "compare_eco_original.csv"), row.names = FALSE)
options(width = 200)
print(tab[, c("data", "model", "sampler", "ms_per_sweep", "mean_agg_W1", "mean_alpha", "mean_clusters",
              "rhat_alpha", "agg_W1_per_s", "alpha_per_s", "clusters_per_s", "unit_median_per_s",
              "unit_p10_per_s", "unit_cor_vs_collapsed", "unit_rmsd_vs_collapsed")], digits = 3)
print(tab_old, digits = 3)
