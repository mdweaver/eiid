# Exactness check of the corrected parametric eco sampler (eco_param_cpp) against
# the posterior computed on a grid. With one normal-inverse-Wishart component,
#   p(W_latent | data) ∝ ev(all logit points) * Jacobian,
# with ev the closed-form NIW marginal likelihood (docs/bayes/verify_collapsed_exact.R),
# and E[mu | data], E[Sigma | data] follow from the NIW posterior given the points.
# Both rate steps (eco's uniform MH, wstep = 0; slice, wstep = 1) are checked.
#   Rscript benchmarks/eco_2x2/eco_fixed/verify_param_exact.R
VERIFY_LIB_ONLY <- TRUE
source("docs/bayes/verify_collapsed_exact.R")   # lev(), bse()
Rcpp::sourceCpp("benchmarks/eco_2x2/eco_fixed/eco_alg1.cpp")
set.seed(2)

post_moments <- function(D, pr) {   # E[mu], E[Sigma] given complete points D (p x m)
  p <- nrow(D); m <- ncol(D); xb <- rowMeans(D); tm <- pr$tau0 + m; nm <- pr$nu0 + m
  Sm <- pr$S0 + tcrossprod(D - xb) + pr$tau0 * m / tm * tcrossprod(xb - pr$mu0)
  list(mu = (pr$tau0 * pr$mu0 + m * xb) / tm, Sig = Sm / (nm - p - 1))
}
line_pt <- function(w1, X, Y) { w2 <- (Y - X * w1) / (1 - X); list(z = qlogis(c(w1, w2)), jac = -sum(log(c(w1, 1 - w1, w2, 1 - w2)))) }

check <- function(name, ex, r) {
  z <- c(W1 = (mean(r$W1[, 1]) - ex$W1) / bse(r$W1[, 1]),
         mu1 = (mean(r$mu[, 1]) - ex$mu1) / bse(r$mu[, 1]),
         Sig11 = (mean(r$Sigma[, 1]) - ex$Sig11) / bse(r$Sigma[, 1]))
  est <- c(mean(r$W1[, 1]), mean(r$mu[, 1]), mean(r$Sigma[, 1]))
  cat(sprintf("%-40s exact W1 %.4f mu1 %.4f Sig11 %.4f | sampler %.4f %.4f %.4f | z %s\n", name,
              ex$W1, ex$mu1, ex$Sig11, est[1], est[2], est[3], paste(sprintf("%.2f", z), collapse = " ")))
  if (!is.null(ex$cdf)) {
    d <- max(abs(ex$cdf - ecdf(r$W1[, 1])(ex$grid)))   # grid = right edges (test 1: grid points)
    cat(sprintf("%-40s max |CDF difference| for W1 = %.4f\n", "", d))
  } else d <- 0
  all(abs(z) < 4) && d < 0.01
}
run <- function(p, X, Y, lo, hi, Wobs, w1obs, pr, wstep, iters = 4e5)
  eco_param_cpp(p, X, Y, lo, hi, Wobs, w1obs, numeric(0), pr$mu0, pr$tau0, pr$nu0, pr$S0,
                iters, 1000, 1, (lo + hi) / 2, wstep)
ok <- TRUE

# Test 1: CAR, two latent units (2-D grid), three observed units.
pr2 <- list(mu0 = c(0, 0), tau0 = 2, nu0 = 4, S0 = diag(10, 2))
Xl <- c(0.6, 0.3); Yl <- c(0.3, 0.6)
lo <- pmax(0, (Xl + Yl - 1) / Xl); hi <- pmin(1, Yl / Xl)
Wobs <- cbind(c(-2, -1.5), c(-1.8, -1.2), c(2, 2.5))
# Grid on a logistic scale along each line, v in [-25, 25], so the mass near the
# ends of the lines (where E[Sigma | W] is large) is resolved; a midpoint grid in
# W1 is not accurate enough for E[Sigma] (it gave 4.280 instead of 4.328).
G <- 300; v <- seq(-25, 25, length.out = G)
pts <- function(k) { s <- plogis(v); w1 <- lo[k] + (hi[k] - lo[k]) * s; w2 <- (Yl[k] - Xl[k] * w1) / (1 - Xl[k])
  list(z = cbind(qlogis(w1), qlogis(w2)), w1 = w1,
       lj = -log(w1) - log1p(-w1) - log(w2) - log1p(-w2) + log(hi[k] - lo[k]) + log(s) + log1p(-s)) }
A <- pts(1); B <- pts(2)
lw <- M1 <- S11 <- matrix(0, G, G)
for (a in 1:G) for (b in 1:G) {
  D <- cbind(A$z[a, ], B$z[b, ], Wobs)
  lw[a, b] <- lev(D, pr2) + A$lj[a] + B$lj[b]
  pm <- post_moments(D, pr2); M1[a, b] <- pm$mu[1]; S11[a, b] <- pm$Sig[1, 1]
}
w <- exp(lw - max(lw)); w <- w / sum(w)
ex1 <- list(W1 = sum(rowSums(w) * A$w1), mu1 = sum(w * M1), Sig11 = sum(w * S11),
            grid = A$w1, cdf = cumsum(rowSums(w)) - rowSums(w) / 2)   # CDF at the grid points
for (ws in 0:1) ok <- check(paste("CAR, 2 latent, wstep", ws), ex1, run(2, Xl, Yl, lo, hi, Wobs, numeric(0), pr2, ws)) && ok

# Test 2: CAR, one latent unit, three observed, one with only W1 observed (W1 posterior only).
X <- 0.6; Y <- 0.3; lo <- 0; hi <- 0.5; w1obs <- 1.0
g <- lo + (hi - lo) * (seq_len(400) - 0.5) / 400
lwt <- sapply(g, function(w1) {
  pt <- line_pt(w1, X, Y); D <- cbind(pt$z, Wobs, c(w1obs, NA))
  f <- function(v) sapply(v, function(y) { D[2, 5] <- y; lev(D, pr2) })
  ref <- max(f(seq(-15, 15, 1)))
  ref + log(integrate(function(v) exp(f(v) - ref), -60, 60, subdivisions = 500, rel.tol = 1e-8)$value) + pt$jac
})
wt <- exp(lwt - max(lwt)); wt <- wt / sum(wt)
ex2 <- list(W1 = sum(wt * g), mu1 = NA, Sig11 = NA, grid = g + (hi - lo) / 800, cdf = cumsum(wt))
for (ws in 0:1) {
  r <- run(2, X, Y, lo, hi, Wobs, w1obs, pr2, ws)
  z <- (mean(r$W1[, 1]) - ex2$W1) / bse(r$W1[, 1]); d <- max(abs(ex2$cdf - ecdf(r$W1[, 1])(ex2$grid)))
  cat(sprintf("%-40s exact W1 %.4f | sampler %.4f | z %.2f | max |CDF difference| %.4f\n",
              paste("CAR, partial unit, wstep", ws), ex2$W1, mean(r$W1[, 1]), z, d))
  ok <- ok && abs(z) < 4 && d < 0.01
}

# Test 3: NCAR (p = 3), one latent unit, four observed units.
pr3 <- list(mu0 = c(0, 0, 0), tau0 = 1, nu0 = 5, S0 = diag(3))
X <- 0.3; Y <- 0.6; lo <- max(0, (Y - (1 - X)) / X); hi <- min(1, Y / X)
Wobs3 <- cbind(c(-1, 0.5, -0.8), c(-1.2, 0.3, -1), c(1.5, -1, 0.7), c(1.2, -0.7, 1))
g <- lo + (hi - lo) * (seq_len(600) - 0.5) / 600
tmp <- lapply(g, function(w1) { pt <- line_pt(w1, X, Y); D <- cbind(c(pt$z, qlogis(X)), Wobs3)
  pm <- post_moments(D, pr3); c(lev(D, pr3) + pt$jac, pm$mu[1], pm$Sig[1, 1]) })
tmp <- do.call(rbind, tmp); w3 <- exp(tmp[, 1] - max(tmp[, 1])); w3 <- w3 / sum(w3)
ex3 <- list(W1 = sum(w3 * g), mu1 = sum(w3 * tmp[, 2]), Sig11 = sum(w3 * tmp[, 3]),
            grid = g + (hi - lo) / 1200, cdf = cumsum(w3))
for (ws in 0:1) ok <- check(paste("NCAR, wstep", ws), ex3, run(3, X, Y, lo, hi, Wobs3, numeric(0), pr3, ws)) && ok

cat("\nOVERALL:", if (ok) "PASS" else "FAIL", "\n")
