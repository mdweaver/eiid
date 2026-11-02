# Follow-up to verify_collapsed_exact.R test 2 (NCAR): batch-means SEs from one chain gave |z| up to 2.8,
# so rerun with 16 independent chains of 4e5 sweeps and between-chain SEs. Exact values from verify_collapsed_exact_output.txt.
#   Rscript docs/bayes/verify_collapsed_exact_rep.R
suppressMessages(library(eiDPM)); library(parallel)
pr3 <- list(mu0 = c(0, 0, 0), tau0 = 1, nu0 = 5, S0 = diag(3), a0 = 2, b0 = 1)
X <- 0.3; Y <- 0.6; lo <- max(0, (Y - (1 - X)) / X); hi <- min(1, Y / X)
Wobs3 <- cbind(c(-1, 0.5, -0.8), c(-1.2, 0.3, -1), c(1.5, -1, 0.7), c(1.2, -0.7, 1))
res <- mclapply(1:16, function(s) { set.seed(100 + s)
  r <- eiDPM:::dpm_collapsed_cpp(3, X, Y, lo, hi, Wobs3, numeric(0), numeric(0), pr3$mu0, pr3$tau0, pr3$nu0, pr3$S0, pr3$a0, pr3$b0, 4e5, 1000, 1, (lo+hi)/2, 1, 1.0)
  c(sapply(1:5, function(k) mean(r$n_clusters == k)), alpha = mean(r$alpha), W1 = mean(r$W1[,1])) }, mc.cores = 4)
m <- do.call(rbind, res); ex <- c(0.230254, 0.454067, 0.240701, 0.066879, 0.008099, 1.6383, 0.4053)
se <- apply(m, 2, sd) / sqrt(nrow(m))
print(round(rbind(exact = ex, sampler = colMeans(m), se = se, z = (colMeans(m) - ex) / se), 5))
