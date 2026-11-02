# Exactness check of the collapsed sampler (src/collapsed.cpp) against the
# posterior computed by brute force.
#
# With five units, every partition of the units (Bell(5) = 52) can be
# enumerated. The exact joint posterior of (partition, alpha, latent rates) is
#   p ∝ Gamma(alpha; a0, b0) * alpha^k Gamma(alpha) / Gamma(alpha + n) * prod_c Gamma(n_c)
#       * prod_c ev(D_c) * Jacobian,
# where ev(D_c) is the normal-inverse-Wishart marginal likelihood of the
# logit-scale points in cluster c, computed here from the closed form
# (not from the sampler's sequential predictive t), and alpha is integrated
# numerically. A partially observed point is integrated over its missing
# coordinate numerically. The sampler's output must match:
#   P(k clusters), E[alpha], and the posterior of the latent unit's W1.
#
#   Rscript docs/bayes/verify_collapsed_exact.R
suppressMessages(library(eiDPM))
cpp <- eiDPM:::dpm_collapsed_cpp
set.seed(1)

lmvgamma <- function(a, p) p * (p - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(p)) / 2))
# log NIW evidence of points (columns of D); mu | Sigma ~ N(mu0, Sigma / tau0), Sigma ~ IW(nu0, S0)
lev <- function(D, pr) {
  p <- nrow(D); m <- ncol(D)
  xb <- rowMeans(D); tm <- pr$tau0 + m; nm <- pr$nu0 + m
  Sm <- pr$S0 + tcrossprod(D - xb) + pr$tau0 * m / tm * tcrossprod(xb - pr$mu0)
  -m * p / 2 * log(pi) + lmvgamma(nm / 2, p) - lmvgamma(pr$nu0 / 2, p) +
    pr$nu0 / 2 * determinant(pr$S0)$modulus - nm / 2 * determinant(Sm)$modulus +
    p / 2 * log(pr$tau0 / tm)
}
partitions <- function(n) {   # restricted growth strings
  out <- list(1L)
  for (i in 2:n) out <- unlist(lapply(out, function(r) lapply(1:(max(r) + 1), function(v) c(r, v))), recursive = FALSE)
  out
}
logA <- function(k, n, pr) {  # log of the integral over alpha, and of E[alpha] times it
  f <- function(a, e) exp(dgamma(a, pr$a0, pr$b0, log = TRUE) + k * log(a) + lgamma(a) - lgamma(a + n) + e * log(a))
  c(log(integrate(f, 0, Inf, e = 0, rel.tol = 1e-10)$value),
    log(integrate(f, 0, Inf, e = 1, rel.tol = 1e-10)$value))
}

# pts: p x n matrix of logit points with NA in the latent unit (column 1) and
# in the missing coordinate of a partial unit; latent unit lies on X*w1 + (1-X)*w2 = Y.
exact <- function(pts, X, Y, lo, hi, pr, partial = NULL, ngrid = 400) {
  n <- ncol(pts); P <- partitions(n)
  g <- lo + (hi - lo) * (seq_len(ngrid) - 0.5) / ngrid
  lA <- t(sapply(1:n, logA, n = n, pr = pr))
  logpost <- matrix(-Inf, length(P), ngrid); lpk <- integer(length(P))
  for (j in seq_along(g)) {
    w1 <- g[j]; w2 <- (Y - X * w1) / (1 - X)
    D <- pts; D[1:2, 1] <- qlogis(c(w1, w2))
    jac <- -sum(log(c(w1, 1 - w1, w2, 1 - w2)))
    for (q in seq_along(P)) {
      r <- P[[q]]; k <- max(r); lpk[q] <- k
      s <- lA[k, 1] + sum(lgamma(tabulate(r))) + jac
      for (c in 1:k) {
        Dc <- D[, r == c, drop = FALSE]
        if (!is.null(partial) && r[partial$col] == c) {
          col <- which(which(r == c) == partial$col); ref <- NULL
          f <- function(v) sapply(v, function(y) { Dc[partial$miss, col] <- y; lev(Dc, pr) })
          ref <- max(f(seq(-15, 15, 1)))
          s <- s + ref + log(integrate(function(v) exp(f(v) - ref), -60, 60, subdivisions = 500, rel.tol = 1e-8)$value)
        } else s <- s + lev(Dc, pr)
      }
      logpost[q, j] <- s
    }
  }
  mx <- max(logpost); w <- exp(logpost - mx)
  dens <- colSums(w); pk <- tapply(rowSums(w), lpk, sum)
  Ealpha <- sum(rowSums(w) * exp(lA[lpk, 2] - lA[lpk, 1])) / sum(w)
  list(grid = g, dens = dens / sum(dens), pk = pk / sum(pk), Ealpha = Ealpha)
}

run <- function(p, X, Y, lo, hi, Wobs, w1obs, pr, iters = 4e5) {
  if (exists("RUN_SAMPLER")) return(RUN_SAMPLER(p, X, Y, lo, hi, Wobs, w1obs, pr, iters))
  r <- cpp(p, X, Y, lo, hi, Wobs, w1obs, numeric(0), pr$mu0, pr$tau0, pr$nu0, pr$S0, pr$a0, pr$b0,
           iters, 1000, 1, (lo + hi) / 2, 1, 1.0)
  r
}
# batch-means standard error
bse <- function(x, b = 100) { m <- colMeans(matrix(x[1:(b * floor(length(x) / b))], ncol = b)); sd(m) / sqrt(b) }

compare <- function(name, ex, r) {
  ks <- max(r$n_clusters)
  cat("\n==", name, "==\n")
  tab <- data.frame(k = 1:5, exact = as.numeric(ex$pk[as.character(1:5)]))
  tab$exact[is.na(tab$exact)] <- 0
  tab$sampler <- sapply(1:5, function(k) mean(r$n_clusters == k))
  tab$se <- sapply(1:5, function(k) bse(as.numeric(r$n_clusters == k)))
  tab$z <- (tab$sampler - tab$exact) / tab$se
  print(format(tab, digits = 4), row.names = FALSE)
  za <- (mean(r$alpha) - ex$Ealpha) / bse(r$alpha)
  cat(sprintf("E[alpha]: exact %.4f  sampler %.4f  z = %.2f\n", ex$Ealpha, mean(r$alpha), za))
  w1 <- r$W1[, 1]; m_ex <- sum(ex$grid * ex$dens)
  zm <- (mean(w1) - m_ex) / bse(w1)
  cdf_ex <- cumsum(ex$dens); cdf_s <- ecdf(w1)(ex$grid + diff(ex$grid[1:2]) / 2)
  cat(sprintf("E[W1]:    exact %.4f  sampler %.4f  z = %.2f\n", m_ex, mean(w1), zm))
  qs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  cat("W1 quantiles exact:  ", sprintf("%.4f", sapply(qs, function(q) ex$grid[which(cdf_ex >= q)[1]])), "\n")
  cat("W1 quantiles sampler:", sprintf("%.4f", quantile(w1, qs)), "\n")
  cat(sprintf("max |CDF difference| = %.4f\n", max(abs(cdf_ex - cdf_s))))
  ok <- all(abs(c(tab$z[tab$se > 0], za, zm)) < 4) && max(abs(cdf_ex - cdf_s)) < 0.01
  cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
  ok
}

if (!exists("VERIFY_LIB_ONLY")) {  # sourcing with VERIFY_LIB_ONLY defined loads the functions only
# --- Test 1: CAR (p = 2). 1 latent unit, 3 observed units, 1 unit with only W1 observed.
pr2 <- list(mu0 = c(0, 0), tau0 = 2, nu0 = 4, S0 = diag(10, 2), a0 = 1, b0 = 0.1)
X <- 0.6; Y <- 0.3; lo <- 0; hi <- 0.5           # W1 in [0, 0.5]
Wobs <- cbind(c(-2, -1.5), c(-1.8, -1.2), c(2, 2.5))
w1obs <- 1.0
pts <- cbind(NA, Wobs, c(w1obs, NA)); pts <- rbind(pts[1, ], pts[2, ])
ex1 <- exact(pts, X, Y, lo, hi, pr2, partial = list(col = 5, miss = 2))
r1 <- run(2, X, Y, lo, hi, Wobs, w1obs, pr2)
ok1 <- compare("CAR, p = 2, with a partially observed unit", ex1, r1)

# --- Test 2: NCAR (p = 3). 1 latent unit, 4 observed units; tighter prior.
pr3 <- list(mu0 = c(0, 0, 0), tau0 = 1, nu0 = 5, S0 = diag(3), a0 = 2, b0 = 1)
X <- 0.3; Y <- 0.6; lo <- 0; hi <- 1              # W1 in [0, 1], W2 in [0.43, 0.86]
lo <- max(0, (Y - (1 - X)) / X); hi <- min(1, Y / X)
Wobs3 <- cbind(c(-1, 0.5, -0.8), c(-1.2, 0.3, -1), c(1.5, -1, 0.7), c(1.2, -0.7, 1))
pts3 <- cbind(c(NA, NA, qlogis(X)), Wobs3)
ex2 <- exact(pts3, X, Y, lo, hi, pr3)
r2 <- run(3, X, Y, lo, hi, Wobs3, numeric(0), pr3)
ok2 <- compare("NCAR, p = 3", ex2, r2)

cat("\nOVERALL:", if (ok1 && ok2) "PASS" else "FAIL", "\n")
}
