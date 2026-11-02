suppressMessages({library(eco); library(eiDPM)}); data(reg, package="eco")
# 1. eco's rWish, transcribed from rand.c, versus the Wishart mean df * S
rWish_eco <- function(S, df) { p <- nrow(S); V <- sapply(0:(p-1), function(i) rchisq(1, df - i - 1))
  Tm <- matrix(0, p, p); diag(Tm) <- sqrt(V); for (i in 1:p) for (j in 1:p) if (i < j) Tm[i, j] <- rnorm(1)
  C <- t(chol(S)); C %*% (t(Tm) %*% Tm) %*% t(C) }
set.seed(1); S <- matrix(c(2, .5, .5, 1), 2)
m <- Reduce(`+`, replicate(40000, rWish_eco(S, 6), simplify = FALSE)) / 40000
cat("eco rWish(S, df=6): mean / S =", round(m / S, 3), " (Wishart(6) gives 6; Wishart(5) gives 5)\n")
# 2. reg, CAR, DP: Gibbs with nu0 = 4 (stated prior) and nu0 = 3 (eco's effective draws)
ref <- NULL
for (s in 1:2) { set.seed(10 + s); f <- ecoNP(Y ~ X, N = N, data = reg, n.draws = 22000, burnin = 2000, thin = 3, parameter = TRUE)
  ref <- cbind(ref, colMeans(f$W[, 1, ])); cat(sprintf("eco run %d: alpha %.2f clusters %.1f\n", s, median(f$alpha), mean(f$nstar))) }
ref <- rowMeans(ref)
for (nu in c(4, 3)) for (s in 1:2) {
  g <- eiDPM(reg$X, reg$Y, reg$N, K = 30, n.draws = 20000, burnin = 2000, thin = 4, prior = list(nu0 = nu), seed = 20 + s)
  pm <- colMeans(g$W1)
  cat(sprintf("Gibbs nu0=%d seed %d: alpha %.2f clusters %.1f  cor(eco) %.3f  rmsd(eco) %.4f\n", nu, s,
              median(g$diagnostics$alpha), mean(g$diagnostics$n_clusters), cor(pm, ref), sqrt(mean((pm - ref)^2)))) }
