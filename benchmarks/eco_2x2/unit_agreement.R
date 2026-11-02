# Supplement to run_benchmark.R: is unit-level disagreement between eiDPM's
# Gibbs sampler and ecoNP larger than the disagreement between independent
# runs of the same sampler? For each DP model, compares unit posterior means of W1:
#   eco run A vs eco run B (2 chains each), Gibbs run A vs Gibbs run B
#   (2 chains each), and pooled Gibbs vs pooled eco.
#   Rscript benchmarks/eco_2x2/unit_agreement.R [cores]
suppressMessages({ library(eco); library(eiDPM); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4
data(reg, package = "eco"); data(census, package = "eco")
datasets <- list(reg = reg, census = census)
out_dir <- file.path("benchmarks", "eco_2x2", "results")
n_draws <- 5000; burnin <- 1000; thin <- 4

one_eco <- function(dd, context, seed) {
  set.seed(seed)
  f <- ecoNP(Y ~ X, N = N, data = dd, context = context, n.draws = n_draws + burnin,
             burnin = burnin, thin = thin - 1)
  colMeans(f$W[, 1, ])
}
one_gibbs <- function(dd, context, seed) {
  f <- eiDPM(dd$X, dd$Y, dd$N, context = context, K = 30, n.draws = n_draws, burnin = burnin,
             thin = thin, seed = seed)
  colMeans(f$W1)
}
res <- NULL
for (dn in names(datasets)) for (context in c(FALSE, TRUE)) {
  dd <- datasets[[dn]]
  cat(format(Sys.time()), dn, if (context) "NCAR" else "CAR", "\n")
  jobs <- c(lapply(1:4, function(s) list("eco", 500 + s)), lapply(1:4, function(s) list("gibbs", 600 + s)))
  pm <- mclapply(jobs, function(j) if (j[[1]] == "eco") one_eco(dd, context, j[[2]]) else
                   one_gibbs(dd, context, j[[2]]), mc.cores = cores, mc.preschedule = FALSE)
  eA <- (pm[[1]] + pm[[2]]) / 2; eB <- (pm[[3]] + pm[[4]]) / 2
  gA <- (pm[[5]] + pm[[6]]) / 2; gB <- (pm[[7]] + pm[[8]]) / 2
  res <- rbind(res, data.frame(
    data = dn, assumption = if (context) "NCAR" else "CAR",
    cor_eco_A_vs_eco_B = cor(eA, eB), cor_gibbs_A_vs_gibbs_B = cor(gA, gB),
    cor_gibbs_vs_eco = cor((gA + gB) / 2, (eA + eB) / 2),
    rmsd_eco_A_vs_B = sqrt(mean((eA - eB)^2)), rmsd_gibbs_A_vs_B = sqrt(mean((gA - gB)^2)),
    rmsd_gibbs_vs_eco = sqrt(mean(((gA + gB) / 2 - (eA + eB) / 2)^2))))
  write.csv(res, file.path(out_dir, "unit_agreement.csv"), row.names = FALSE)
}
print(res, digits = 3)
