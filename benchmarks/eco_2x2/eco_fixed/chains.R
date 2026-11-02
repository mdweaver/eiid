# Follow-up to compare.R: per-chain behaviour with many chains.
#  (a) census CAR/NCAR: 8 chains each of eco_slice (Algorithm 1 + remixing, slice step) and collapsed,
#      22000 sweeps: does Algorithm 1 have a metastable extra cluster, and does the collapsed sampler find it?
#  (b) reg CAR: eco's MH rate step, 4 chains of 220000 sweeps: does it reach the slice / collapsed answer?
#   Rscript benchmarks/eco_2x2/eco_fixed/chains.R [cores]
suppressMessages({ library(eiDPM); library(parallel) })
Rcpp::sourceCpp("benchmarks/eco_2x2/eco_fixed/eco_alg1.cpp")
cores <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 4 }
data(reg, package = "eco"); data(census, package = "eco")
ds <- list(reg = reg, census = census)
one <- function(dn, ctx, sampler, chain, sweeps, init_k_alpha = 1.0) {
  dd <- ds[[dn]]; d <- eiDPM:::dpmData(dd$X, dd$Y, ctx); pr <- eiDPM:::.dpm_prior(list(), d$p)
  set.seed(5000 + chain)
  burn <- sweeps %/% 11; thin <- 10; keep <- (sweeps - burn) %/% thin
  a <- list(d$p, d$X, d$Y, d$L, d$U, d$Wobs, d$w1obs, d$w2obs, pr$mu0, pr$tau0, pr$nu0, pr$S0, pr$a0, pr$b0,
            as.integer(burn + keep * thin), as.integer(burn), as.integer(thin), (d$L + d$U) / 2)
  r <- switch(sampler,
    eco_mh = do.call(eco_alg1_cpp, c(a, list(0L, 1L, 1.0))),
    eco_slice = do.call(eco_alg1_cpp, c(a, list(1L, 1L, 1.0))),
    collapsed = do.call(eiDPM:::dpm_collapsed_cpp, c(a, list(1L, 1.0))))
  N <- dd$N[d$lat]; agg <- as.vector(r$W1 %*% N / sum(N))
  half <- seq_len(length(agg)) > length(agg) / 2
  data.frame(data = dn, model = if (ctx) "NCAR" else "CAR", sampler, chain, sweeps,
             agg_W1 = mean(agg), agg_W1_2nd_half = mean(agg[half]), alpha = mean(r$alpha),
             clusters = mean(r$n_clusters), share_k1 = mean(r$n_clusters == 1),
             unit_W1 = I(list(colMeans(r$W1))))
}
jobs <- rbind(expand.grid(dn = "census", ctx = c(FALSE, TRUE), sampler = c("eco_slice", "collapsed"),
                          chain = 1:8, sweeps = 22000, stringsAsFactors = FALSE),
              expand.grid(dn = "reg", ctx = FALSE, sampler = c("eco_mh", "eco_slice", "collapsed"),
                          chain = 1:4, sweeps = 220000, stringsAsFactors = FALSE))
res <- mclapply(seq_len(nrow(jobs)), function(j) with(jobs[j, ], one(dn, ctx, sampler, chain, sweeps)),
                mc.cores = cores, mc.preschedule = FALSE)
tab <- do.call(rbind, res)
saveRDS(tab, "benchmarks/eco_2x2/eco_fixed/chains.rds")
options(width = 160)
print(tab[, c("data", "model", "sampler", "chain", "sweeps", "agg_W1", "agg_W1_2nd_half", "alpha", "clusters", "share_k1")], digits = 3)
# unit agreement of each sampler's pooled unit means with pooled collapsed means
for (key in unique(paste(tab$data, tab$model))) {
  t <- tab[paste(tab$data, tab$model) == key, ]
  ref <- Reduce(`+`, t$unit_W1[t$sampler == "collapsed"]) / sum(t$sampler == "collapsed")
  for (s in unique(t$sampler)) {
    u <- Reduce(`+`, t$unit_W1[t$sampler == s]) / sum(t$sampler == s)
    cat(sprintf("%-12s %-10s pooled unit means vs collapsed: cor %.4f  RMSD %.4f\n", key, s, cor(u, ref), sqrt(mean((u - ref)^2))))
  }
}
