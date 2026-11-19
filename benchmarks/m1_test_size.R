# Size of the test of M1's testable part (M1 against M2; P3 Prop. P3.6(b)) in the WS6 design S0 of
# docs/extension/PHASE3_inference.md section 8.3: 2 x 2, 20 districts x 50 polls, N = 400,
# no district shifts, Beta(60 m, 60 (1 - m)) unit rates. (G - 1)(C - 1)R(R - 1)/2 = 19 restrictions; M1 holds (no shifts).
# Reports the rejection rate at nominal 5% for the HC3 Wald statistic against chi^2, the F form
# (W / df against F(df, n - k)) and the restricted wild bootstrap (the default from 0.3.0, decision 3(b)).
# Usage: Rscript benchmarks/m1_test_size.R [reps] [B] [cores]
args <- as.numeric(commandArgs(TRUE))
reps <- if (length(args) >= 1) args[1] else 1000
B <- if (length(args) >= 2) args[2] else 199
cores <- if (length(args) >= 3) args[3] else 1
library(eiPartialID)
simS0 <- function(seed, D = 20, PD = 50) {
  set.seed(seed)
  m <- stats::runif(D, 0.1, 0.9)
  d <- rep(seq_len(D), each = PD)
  x <- pmin(pmax(m[d] + stats::runif(D * PD, -0.15, 0.15), 0.01), 0.99)
  mw <- 0.55 - 0.10 * x; mb <- 0.35 + 0.20 * x
  bw <- stats::rbeta(length(x), 60 * mw, 60 * (1 - mw)); bb <- stats::rbeta(length(x), 60 * mb, 60 * (1 - mb))
  t <- x * bb + (1 - x) * bw
  list(X = cbind(x, 1 - x), T = cbind(t, 1 - t), N = rep(400, length(x)), g = d)
}
one <- function(s) {
  d <- simS0(70000 + s)
  set.seed(90000 + s)
  r <- eiPartialID:::.rxc_m1_test(d$X, d$T, d$N, d$g, B = B)
  c(chisq = r$p.value_chisq, F = r$p.value_F, wild = r$p.value, W = r$statistic, df = r$df)
}
res <- if (cores > 1) parallel::mclapply(seq_len(reps), one, mc.cores = cores) else lapply(seq_len(reps), one)
res <- do.call(rbind, res)
rej <- colMeans(res[, c("chisq", "F", "wild")] <= 0.05)
se <- sqrt(rej * (1 - rej) / reps)
cat(sprintf("design S0, %d replications, B = %d, df = %d\n", reps, B, res[1, "df"]))
cat(sprintf("Wald statistic: mean %.1f, variance %.1f (chi^2_%d: %d, %d)\n", mean(res[, "W"]), stats::var(res[, "W"]),
            res[1, "df"], res[1, "df"], 2 * res[1, "df"]))
for (k in names(rej)) cat(sprintf("rejection at 5%%, %-5s: %.3f (MC SE %.3f)\n", k, rej[k], se[k]))
