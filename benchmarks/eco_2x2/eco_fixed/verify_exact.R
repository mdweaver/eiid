# Brute-force exactness check of the corrected eco sampler (eco_alg1.cpp),
# reusing docs/bayes/verify_collapsed_exact.R with the sampler swapped in.
#   Rscript benchmarks/eco_2x2/eco_fixed/verify_exact.R <variant>
# variant: mh_grouped (eco's MH step, cluster-form sum), mh_units (eco's MH step,
# eco's loop over units), slice_grouped.
v <- commandArgs(trailingOnly = TRUE)[1]
Rcpp::sourceCpp("benchmarks/eco_2x2/eco_fixed/eco_alg1.cpp")
opt <- switch(v, mh_grouped = c(0, 1), mh_units = c(0, 0), slice_grouped = c(1, 1))
cat("variant:", v, "\n")
RUN_SAMPLER <- function(p, X, Y, lo, hi, Wobs, w1obs, pr, iters)
  eco_alg1_cpp(p, X, Y, lo, hi, Wobs, w1obs, numeric(0), pr$mu0, pr$tau0, pr$nu0, pr$S0,
               pr$a0, pr$b0, iters, 1000, 1, (lo + hi) / 2, opt[1], opt[2], 1.0)
source("docs/bayes/verify_collapsed_exact.R")
