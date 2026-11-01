# Verifies the density of Y given X under the logit-normal model (PHASE1 sec. 4.1):
#   p(y | x) = (1 / (1 - x)) * int_L^U f_W(w1, w2(w1)) dw1,
#   f_W(w) = phi(logit w) / prod_j w_j (1 - w_j),
# against Monte Carlo, and contrasts it with the Euclidean arc-length measure
# in logit space used by ecoML's E-step (fintegrate.c, `pfact`).
# Run: Rscript docs/bayes/verify_line_density.R > docs/bayes/verify_line_density_output.txt
set.seed(1)
mu <- c(0.8, -0.6); s <- c(1.2, 0.9); rho <- 0.4
S <- diag(s) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(s); Si <- solve(S)
phi <- function(z) exp(-0.5 * rowSums((z %*% Si) * z)) / (2 * pi * sqrt(det(S)))
dens <- function(y, x, measure = c("correct", "arc")) {
  measure <- match.arg(measure)
  L <- max(0, (x + y - 1) / x); U <- min(1, y / x)
  f <- function(w1) {
    w2 <- (y - x * w1) / (1 - x)
    ok <- w1 > 0 & w1 < 1 & w2 > 0 & w2 < 1; out <- numeric(length(w1))
    a <- w1[ok]; b <- w2[ok]
    z <- cbind(qlogis(a) - mu[1], qlogis(b) - mu[2])
    if (measure == "correct") out[ok] <- phi(z) / (a * (1 - a) * b * (1 - b)) / (1 - x)
    else {                         # arc length in logit space, per unit of w1
      d1 <- 1 / (a * (1 - a)); d2 <- x / (1 - x) / (b * (1 - b))
      out[ok] <- phi(z) * sqrt(d1^2 + d2^2)
    }
    out[!is.finite(out)] <- 0; out }
  integrate(f, L, U, rel.tol = 1e-9, subdivisions = 500)$value
}
nmc <- 2e6
cat("Check 1: p(y|x) by the line integral vs Monte Carlo (bin probabilities, width 0.02)\n")
res <- NULL
for (x in c(0.15, 0.5, 0.85)) {
  W <- plogis(MASS::mvrnorm(nmc, mu, S)); Y <- x * W[, 1] + (1 - x) * W[, 2]
  grid <- seq(0.01, 0.99, by = 0.02)
  mc <- sapply(grid, function(g) mean(abs(Y - g) < 0.01))
  lc <- sapply(grid, function(g) integrate(Vectorize(function(y) dens(y, x)), g - 0.01, g + 0.01)$value)
  la <- sapply(grid, function(g) integrate(Vectorize(function(y) dens(y, x, "arc")), g - 0.01, g + 0.01)$value)
  tot_a <- integrate(Vectorize(function(y) dens(y, x, "arc")), 0, 1)$value
  la <- la / tot_a  # normalise the arc-length version to a density in y
  se <- sqrt(pmax(lc, 1 / nmc) * (1 - lc) / nmc)  # binomial SE under the predicted probability
  res <- rbind(res, data.frame(x = x, max_abs_z_correct = max(abs(lc - mc) / se),
                               total_correct = sum(lc) + 0,
                               max_abs_diff_correct = max(abs(lc - mc)),
                               max_abs_diff_arc_normalised = max(abs(la - mc))))
}
print(res, digits = 4)
pass1 <- all(res$max_abs_z_correct < 5)
cat(sprintf("PASS(line integral matches Monte Carlo, all |z| < 5): %s\n", pass1))
cat("\nCheck 2: is the arc-length factor proportional to the correct Jacobian along a line?\n")
x <- 0.4; y <- 0.55; L <- max(0, (x + y - 1) / x); U <- min(1, y / x)
w1 <- seq(L + 0.01, U - 0.01, length.out = 7); w2 <- (y - x * w1) / (1 - x)
ratio <- sqrt((1 / (w1 * (1 - w1)))^2 + (x / (1 - x) / (w2 * (1 - w2)))^2) /
         (1 / (w1 * (1 - w1) * w2 * (1 - w2)))
print(data.frame(w1 = round(w1, 3), arc_over_correct = signif(ratio, 4)))
pass2 <- diff(range(ratio)) / mean(ratio) > 0.05
cat(sprintf("PASS(ratio not constant along the line, so the two measures differ): %s\n", pass2))
cat(sprintf("\nOVERALL: %s\n", if (pass1 && pass2) "PASS" else "FAIL"))
