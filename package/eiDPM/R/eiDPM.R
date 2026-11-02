#' Prepare 2x2 data following eco's boundary rule
#'
#' Units with \code{0 < X < 1} and \code{0 < Y < 1} have latent rates on their
#' tomography line. Units with \code{Y} equal to 0 or 1 have both rates known.
#' Units with \code{X = 1} (\code{X = 0}) have \code{W1} (\code{W2}) known; these
#' are allowed only without context. Exact 0 and 1 are replaced by
#' \code{eps} and \code{1 - eps} before the logit, as in \pkg{eco}.
#'
#' @param X,Y shares in [0, 1].
#' @param context logical; \code{TRUE} for NCAR (the joint model of rates and
#'   composition).
#' @param eps replacement for exact 0 / 1.
#' @param shrink relative shrinkage of each unit's bounds, keeping implied rates
#'   strictly inside (0, 1).
#' @return a list used by \code{eiDPM()}.
#' @export
dpmData <- function(X, Y, context = FALSE, eps = 1e-4, shrink = 1e-8) {
  if (length(X) != length(Y)) stop("X and Y must have the same length")
  if (any(!is.finite(X) | !is.finite(Y)) || any(X < 0 | X > 1 | Y < 0 | Y > 1))
    stop("X and Y must lie in [0, 1]")
  cl <- function(w) pmin(pmax(w, eps), 1 - eps)
  lg <- function(w) stats::qlogis(cl(w))
  idx <- seq_along(X)
  x1 <- idx[X == 1]; x0 <- idx[X == 0]
  yb <- idx[X > 0 & X < 1 & (Y == 0 | Y == 1)]
  lat <- idx[X > 0 & X < 1 & Y > 0 & Y < 1]
  if (context && length(c(x1, x0)) > 0)
    stop("context = TRUE needs 0 < X < 1 for every unit")
  p <- if (context) 3L else 2L
  Xl <- X[lat]; Yl <- Y[lat]
  L <- pmax(0, (Xl + Yl - 1) / Xl); U <- pmin(1, Yl / Xl)
  w <- U - L
  L <- L + shrink * w; U <- U - shrink * w
  Wobs <- matrix(0, p, length(yb))
  if (length(yb)) {
    Wobs[1, ] <- lg(Y[yb]); Wobs[2, ] <- lg(Y[yb])
    if (context) Wobs[3, ] <- stats::qlogis(X[yb])
  }
  list(p = p, n_all = length(X), lat = lat, yb = yb, x1 = x1, x0 = x0,
       X = Xl, Y = Yl, L = L, U = U, Wobs = Wobs,
       w1obs = lg(Y[x1]), w2obs = lg(Y[x0]),
       Wyb = cl(Y[yb]), Wx1 = cl(Y[x1]), Wx0 = cl(Y[x0]))
}

.dpm_prior <- function(prior, p) {
  pr <- list(mu0 = 0, tau0 = 2, nu0 = 4, S0 = 10, a0 = 1, b0 = 0.1)
  pr[names(prior)] <- prior
  if (length(pr$mu0) == 1) pr$mu0 <- rep(pr$mu0, p)
  if (!is.matrix(pr$S0)) pr$S0 <- diag(pr$S0, p)
  pr
}

.stan_model_cache <- new.env(parent = emptyenv())

.dpm_stan_model <- function(cache = getOption("eiDPM.cache", TRUE)) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("method = \"stan\" needs the 'rstan' package")
  if (is.null(.stan_model_cache$m)) {
    f <- system.file("stan", "dpm2x2.stan", package = "eiDPM")
    key <- paste0("dpm2x2-", unname(tools::md5sum(f)), "-rstan",
                  utils::packageVersion("rstan"), ".rds")
    dir <- tools::R_user_dir("eiDPM", "cache")
    path <- file.path(dir, key)
    m <- NULL
    if (isTRUE(cache) && file.exists(path)) m <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(m)) {
      m <- rstan::stan_model(file = f, model_name = "dpm2x2", auto_write = FALSE)
      if (isTRUE(cache)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        try(saveRDS(m, path), silent = TRUE)
      }
    }
    .stan_model_cache$m <- m
  }
  .stan_model_cache$m
}

#' Bayesian mixture models for 2x2 ecological inference
#'
#' Fits the parametric (\code{K = 1}) or Dirichlet-process mixture models of
#' Imai, Lu and Strauss (2008) under CAR
#' (\code{context = FALSE}) or NCAR (\code{context = TRUE}), with the same
#' default priors as \pkg{eco}.
#'
#' The DP mixture is fitted exactly (\code{K = Inf}, the default) by a collapsed
#' Gibbs sampler (Neal 2000, Algorithm 3) with the conjugate normal-inverse-Wishart
#' base measure integrated out and the Escobar-West (1995) update for
#' \code{alpha}. It has no truncation and no label ordering. A finite
#' \code{K > 1} instead fits the truncated stick-breaking approximation with a
#' blocked Gibbs sampler (Ishwaran and James 2001), or with Stan; that sampler
#' can under-represent the number of clusters when \code{alpha} is large
#' (see \code{docs/bayes/eco_sampler_review.md}).
#'
#' The Stan model is compiled at first use and cached under
#' \code{tools::R_user_dir("eiDPM", "cache")}; set
#' \code{options(eiDPM.cache = FALSE)} to compile without caching.
#' Both Gibbs samplers start with every unit in one cluster and
#' \code{alpha = 1}. For the truncated sampler, check \code{diagnostics$pi_last}
#' (the mass of the last component).
#'
#' @param X,Y shares (composition and outcome) in [0, 1].
#' @param N optional unit sizes, used for weighted aggregates.
#' @param context logical; NCAR if \code{TRUE}.
#' @param K \code{Inf} (default) for the exact DP mixture; 1 for the parametric
#'   model; a finite value > 1 for the truncated DP (blocked Gibbs or Stan).
#' @param method \code{"gibbs"} (C++) or \code{"stan"} (HMC; finite \code{K} only).
#' @param n.draws,burnin,thin iterations (Gibbs) or per-chain sampling /
#'   warm-up iterations (Stan: \code{n.draws} samples after \code{burnin} warm-up).
#' @param chains number of chains.
#' @param cores cores for Stan chains.
#' @param prior list overriding \code{mu0, tau0, nu0, S0, a0, b0}.
#' @param seed random seed.
#' @param slice_steps slice-sampling updates per latent rate per sweep (Gibbs).
#' @param ... passed to \code{rstan::sampling()}.
#' @return an object of class \code{"eiDPM"} with draws of the unit rates
#'   (\code{W1}, \code{W2}: draws x units, all units in input order).
#' @export
eiDPM <- function(X, Y, N = NULL, context = FALSE, K = Inf,
                  method = c("gibbs", "stan"), n.draws = 5000, burnin = 1000,
                  thin = 1, chains = 1, cores = 1, prior = list(), seed = NULL,
                  slice_steps = 1, ...) {
  method <- match.arg(method)
  d <- dpmData(X, Y, context)
  pr <- .dpm_prior(prior, d$p)
  if (!is.null(seed)) set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  diag <- list()
  if (!(K == Inf || (K >= 1 && K == round(K)))) stop("K must be Inf or a positive integer")
  if (method == "stan" && is.infinite(K)) stop("method = \"stan\" needs a finite K")
  if (method == "gibbs" && is.infinite(K)) {
    res <- lapply(seq_len(chains), function(ch) {
      tc <- proc.time()[["elapsed"]]
      r <- dpm_collapsed_cpp(d$p, d$X, d$Y, d$L, d$U, d$Wobs, d$w1obs, d$w2obs,
                             pr$mu0, pr$tau0, pr$nu0, pr$S0, pr$a0, pr$b0,
                             as.integer(n.draws + burnin), as.integer(burnin), as.integer(thin),
                             (d$L + d$U) / 2, as.integer(slice_steps), 1)
      r$seconds <- proc.time()[["elapsed"]] - tc
      r
    })
    W1l <- do.call(rbind, lapply(res, `[[`, "W1"))
    W2l <- do.call(rbind, lapply(res, `[[`, "W2"))
    chain <- rep(seq_len(chains), sapply(res, function(r) nrow(r$W1)))
    diag$alpha <- unlist(lapply(res, `[[`, "alpha"))
    diag$n_clusters <- unlist(lapply(res, `[[`, "n_clusters"))
    diag$largest_cluster <- unlist(lapply(res, `[[`, "largest_cluster"))
    diag$chain_seconds <- sapply(res, `[[`, "seconds")
  } else if (method == "gibbs") {
    res <- lapply(seq_len(chains), function(ch) {
      tc <- proc.time()[["elapsed"]]
      r <- dpm_gibbs_cpp(d$p, as.integer(K), d$X, d$Y, d$L, d$U, d$Wobs, d$w1obs, d$w2obs,
                    pr$mu0, pr$tau0, pr$nu0, pr$S0, pr$a0, pr$b0,
                    as.integer(n.draws + burnin), as.integer(burnin), as.integer(thin),
                    (d$L + d$U) / 2, as.integer(slice_steps))
      r$seconds <- proc.time()[["elapsed"]] - tc
      r
    })
    W1l <- do.call(rbind, lapply(res, `[[`, "W1"))
    W2l <- do.call(rbind, lapply(res, `[[`, "W2"))
    chain <- rep(seq_len(chains), sapply(res, function(r) nrow(r$W1)))
    diag$alpha <- unlist(lapply(res, `[[`, "alpha"))
    diag$pi_last <- unlist(lapply(res, `[[`, "pi_last"))
    diag$n_clusters <- unlist(lapply(res, `[[`, "n_clusters"))
    diag$chain_seconds <- sapply(res, `[[`, "seconds")
  } else {
    m <- .dpm_stan_model()
    sd <- list(p = d$p, K = as.integer(K), n = length(d$X), X = d$X, Y = d$Y,
               L = d$L, U = d$U, n_obs = ncol(d$Wobs), Wobs = d$Wobs,
               n_x1 = length(d$w1obs), w1obs = as.array(d$w1obs),
               n_x0 = length(d$w2obs), w2obs = as.array(d$w2obs),
               mu0 = pr$mu0, tau0 = pr$tau0, nu0 = pr$nu0, S0 = pr$S0,
               a0 = pr$a0, b0 = pr$b0)
    init <- function() list(W1 = d$L + (d$U - d$L) * stats::runif(length(d$L), 0.25, 0.75),
                            alpha = 1, v = as.array(rep(0.5, K - 1)),
                            csq = lapply(seq_len(K), function(k) as.array(pr$nu0 - seq_len(d$p) + 1)),
                            zoff = lapply(seq_len(K), function(k) as.array(rep(0, d$p * (d$p - 1) / 2))),
                            zmu = lapply(seq_len(K), function(k) as.array(stats::rnorm(d$p, 0, 0.5))))
    fit <- rstan::sampling(m, data = sd, chains = chains, cores = cores,
                           iter = burnin + n.draws, warmup = burnin, thin = thin,
                           init = init, seed = if (is.null(seed)) sample.int(1e6, 1) else seed,
                           refresh = 0, ...)
    W1a <- as.array(fit, pars = "W1"); W2a <- as.array(fit, pars = "W2")
    flat <- function(a) do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
    W1l <- flat(W1a); W2l <- flat(W2a)
    if (is.null(dim(W1l))) { W1l <- matrix(W1l, ncol = 1); W2l <- matrix(W2l, ncol = 1) }
    chain <- rep(seq_len(dim(W1a)[2]), each = dim(W1a)[1])
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    diag$divergences <- sum(sapply(sp, function(s) sum(s[, "divergent__"])))
    diag$treedepth_max <- sum(sapply(sp, function(s) sum(s[, "treedepth__"] >= 10)))
    diag$pi_last <- as.vector(as.array(fit, pars = "pi_last"))
    diag$alpha <- as.vector(as.array(fit, pars = "alpha"))
    diag$chain_seconds <- rowSums(rstan::get_elapsed_time(fit))
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  nd <- nrow(W1l)
  W1 <- W2 <- matrix(NA_real_, nd, d$n_all)
  W1[, d$lat] <- W1l; W2[, d$lat] <- W2l
  if (length(d$yb)) { W1[, d$yb] <- rep(d$Wyb, each = nd); W2[, d$yb] <- rep(d$Wyb, each = nd) }
  if (length(d$x1)) W1[, d$x1] <- rep(d$Wx1, each = nd)
  if (length(d$x0)) W2[, d$x0] <- rep(d$Wx0, each = nd)
  structure(list(W1 = W1, W2 = W2, chain = chain, X = X, Y = Y, N = N,
                 context = context, K = K, method = method, elapsed = elapsed,
                 diagnostics = diag, prior = pr, call = match.call()),
            class = "eiDPM")
}

#' Aggregate summaries of an eiDPM fit
#'
#' Draws of the in-sample mean rates (unweighted, over units with the rate
#' defined) and, if \code{N} was given, the size-weighted aggregates
#' \eqn{\sum N X W_1 / \sum N X} and \eqn{\sum N (1-X) W_2 / \sum N (1-X)}.
#' @param object an \code{"eiDPM"} fit.
#' @return a matrix of draws, one column per aggregate.
#' @export
dpmAggregates <- function(object) {
  W1 <- object$W1; W2 <- object$W2
  out <- cbind(mean_W1 = rowMeans(W1, na.rm = TRUE), mean_W2 = rowMeans(W2, na.rm = TRUE))
  if (!is.null(object$N)) {
    a1 <- object$N * object$X; a2 <- object$N * (1 - object$X)
    W1z <- W1; W1z[is.na(W1z)] <- 0; W2z <- W2; W2z[is.na(W2z)] <- 0
    out <- cbind(out, agg_W1 = as.vector(W1z %*% a1) / sum(a1),
                 agg_W2 = as.vector(W2z %*% a2) / sum(a2))
  }
  out
}

#' Unit-level posterior summaries
#' @param object an \code{"eiDPM"} fit.
#' @param prob interval probability.
#' @return data frame with posterior means and equal-tailed intervals of W1, W2.
#' @export
dpmUnits <- function(object, prob = 0.9) {
  q <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
  f <- function(W) t(apply(W, 2, function(w) c(mean(w), stats::quantile(w, q, na.rm = TRUE))))
  a <- f(object$W1); b <- f(object$W2)
  data.frame(W1 = a[, 1], W1_lo = a[, 2], W1_hi = a[, 3],
             W2 = b[, 1], W2_lo = b[, 2], W2_hi = b[, 3])
}

#' @export
print.eiDPM <- function(x, ...) {
  cat(sprintf("eiDPM: %s, %s, method = %s, %d draws, %.1f s\n",
              if (x$context) "NCAR" else "CAR",
              if (x$K == 1) "parametric" else if (is.infinite(x$K)) "DP mixture (exact, collapsed)"
              else sprintf("DP mixture (truncated, K = %d)", x$K), x$method, nrow(x$W1), x$elapsed))
  print(summary(x))
  invisible(x)
}

#' @export
summary.eiDPM <- function(object, ...) {
  a <- dpmAggregates(object)
  t(apply(a, 2, function(v) c(mean = mean(v), sd = stats::sd(v),
                              q05 = unname(stats::quantile(v, 0.05)),
                              q95 = unname(stats::quantile(v, 0.95)))))
}
