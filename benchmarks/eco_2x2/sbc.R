# Simulation-based calibration (Talts et al. 2018) of the eiDPM Gibbs sampler
# and of eco::ecoNP for the 2x2 CAR DP-mixture model of Imai, Lu & Strauss (2008).
#
# For each replication: draw (alpha, sticks, components) from the prior, unit
# rates W, compositions X ~ U(0.05, 0.95) (independent of W: CAR), and Y. Fit,
# and record the rank of each true quantity among L thinned posterior draws.
# If the sampler targets the right posterior and mixes, ranks are uniform.
# Replications with any Y within 1e-6 of 0 or 1 are discarded: this conditions on
# the observed data only, which leaves the posterior, and hence SBC, unchanged.
#
#   A: generate from the truncated prior with K = 30, fit Gibbs with K = 30 (exactness)
#   B: generate with K = 300 (effectively the DP), fit Gibbs with K = 50 (truncation)
#   C: ecoNP on the B simulations (first n_eco of them)
#
#   Rscript benchmarks/eco_2x2/sbc.R [cores] [reps] [n_eco]
suppressMessages({ library(eco); library(eiDPM); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4
reps <- if (length(args) >= 2) as.integer(args[2]) else 400
n_eco <- if (length(args) >= 3) as.integer(args[3]) else 200
n_units <- 100
L <- 100                       # posterior draws kept per fit
burn <- 1000; thin <- 50       # Gibbs: 1000 + 100 * 50 sweeps
out_dir <- file.path("benchmarks", "eco_2x2", "results")
prior <- list(mu0 = c(0, 0), tau0 = 2, nu0 = 4, S0 = diag(10, 2), a0 = 1, b0 = 0.1)

simulate <- function(seed, Kgen) {
  set.seed(seed)
  repeat {
    alpha <- rgamma(1, prior$a0, prior$b0)
    v <- c(rbeta(Kgen - 1, 1, alpha), 1)
    pi <- v * cumprod(c(1, 1 - v[-Kgen]))
    z <- sample.int(Kgen, n_units, replace = TRUE, prob = pi)
    Ws <- matrix(NA, n_units, 2)
    for (k in unique(z)) {
      Sig <- solve(rWishart(1, prior$nu0, solve(prior$S0))[, , 1])
      mu <- MASS::mvrnorm(1, prior$mu0, Sig / prior$tau0)
      Ws[z == k, ] <- MASS::mvrnorm(sum(z == k), mu, Sig)
    }
    W <- plogis(Ws); X <- runif(n_units, 0.05, 0.95)
    Y <- X * W[, 1] + (1 - X) * W[, 2]
    if (all(Y > 1e-6 & Y < 1 - 1e-6)) break   # data-only rejection
  }
  list(X = X, Y = Y, W = W, alpha = alpha, clusters = length(unique(z)))
}

truth_vec <- function(s) c(W1_u1 = s$W[1, 1], W1_u2 = s$W[2, 1], W1_u3 = s$W[3, 1],
                           mean_W1 = mean(s$W[, 1]), mean_W2 = mean(s$W[, 2]),
                           alpha = s$alpha, clusters = s$clusters)
draws_mat <- function(W1, W2, alpha, clusters)
  cbind(W1_u1 = W1[, 1], W1_u2 = W1[, 2], W1_u3 = W1[, 3],
        mean_W1 = rowMeans(W1), mean_W2 = rowMeans(W2), alpha = alpha, clusters = clusters)
rank_of <- function(d, t) {   # ties (discrete cluster counts) broken at random
  lo <- sum(d < t); eq <- sum(d == t); lo + if (eq > 0) sample.int(eq + 1, 1) - 1 else 0
}

fit_gibbs <- function(s, K, seed) {
  f <- eiDPM(s$X, s$Y, K = K, n.draws = L * thin, burnin = burn, thin = thin, prior = prior, seed = seed)
  draws_mat(f$W1, f$W2, f$diagnostics$alpha, f$diagnostics$n_clusters)
}
fit_eco <- function(s, seed) {
  set.seed(seed)
  f <- ecoNP(Y ~ X, data = data.frame(X = s$X, Y = s$Y), n.draws = burn + L * thin,
             burnin = burn, thin = thin - 1, parameter = TRUE)
  draws_mat(f$W[, 1, ], f$W[, 2, ], as.vector(f$alpha), as.vector(f$nstar))
}

one <- function(r, design) {
  Kgen <- if (design == "A") 30 else 300
  s <- simulate(10000 + r + if (design == "A") 0 else 50000, Kgen)
  t0 <- proc.time()[["elapsed"]]
  d <- switch(design, A = fit_gibbs(s, 30, r), B = fit_gibbs(s, 50, r), C = fit_eco(s, r))
  tv <- truth_vec(s)
  c(rep = r, sapply(names(tv), function(q) rank_of(d[, q], tv[[q]])),
    seconds = proc.time()[["elapsed"]] - t0)
}

res <- list()
for (design in c("A", "B", "C")) {
  R <- if (design == "C") n_eco else reps
  cat(format(Sys.time()), "design", design, R, "replications\n")
  rk <- mclapply(seq_len(R), function(r) tryCatch(one(r, design), error = function(e) NULL),
                 mc.cores = cores, mc.preschedule = FALSE)
  ok <- !sapply(rk, is.null)
  m <- as.data.frame(do.call(rbind, rk[ok])); m$design <- design
  res[[design]] <- m
  cat("  failed fits:", sum(!ok), "\n")
  write.csv(do.call(rbind, res), file.path(out_dir, "sbc_ranks.csv"), row.names = FALSE)
}

# Uniformity: chi-square over 10 equal bins of the rank 0..L
qs <- c("W1_u1", "W1_u2", "W1_u3", "mean_W1", "mean_W2", "alpha", "clusters")
tab <- NULL
for (design in names(res)) for (q in qs) {
  x <- res[[design]][[q]]
  b <- table(cut(x, breaks = seq(-0.5, L + 0.5, length.out = 11)))
  ch <- suppressWarnings(chisq.test(b))
  tab <- rbind(tab, data.frame(design = design, quantity = q, reps = length(x),
                               chisq_p = ch$p.value, bins = paste(as.vector(b), collapse = " "),
                               mean_rank_frac = mean(x) / L))
}
tab$median_seconds <- sapply(tab$design, function(d) median(res[[d]]$seconds))
write.csv(tab, file.path(out_dir, "sbc_summary.csv"), row.names = FALSE)
print(tab, digits = 3)
