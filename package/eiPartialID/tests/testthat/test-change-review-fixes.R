# Regression tests for the PHASE2e code review (docs/extension/reviews/PHASE2e_code_review.md, issues
# 1-13 and 15), the Def. 5.1 item-3 scope, and an R = 3 brute force of rows on changes (oracle set C).

rw_data <- function() {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10))
  yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  list(early = mk(1), base = mk(2), post = mk(3))
}
rw_fit <- function(d = rw_data(), ...) eiChange(d, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet",
                                                col = "yes", integer = "lp", ...)
rw_att <- function(f, tr, ...) eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = 0,
                                     no_anticipation = TRUE, ...)

test_that("issue 1: unbounded causal LPs are infinite, never EMPTY", {
  f <- rw_fit()
  a <- rw_att(f, no_trend(), rate_limits = FALSE)
  expect_equal(c(a$results[[1]]$lower, a$results[[1]]$upper), c(-Inf, Inf))
  expect_false(a$results[[1]]$empty)
  b <- eiATT(f, treated = "vet", population = "fixed", trend = mtr_nd_ps(nd = FALSE, ps = FALSE), no_anticipation = TRUE,
             rate_limits = FALSE)
  expect_equal(b$results[[1]]$upper, Inf)
  expect_false(any(grepl("EMPTY", capture.output(print(b)))))
})

test_that("issue 2: an infeasible count layer gives EMPTY for the dial and RM (no crash)", {
  f6 <- rw_fit(restrictions = list(known = list(row = "vet", col = "yes", election = "post", unit = "P1", value = 11)))
  for (tr in list(parallel_trends(), trend_dial(0.1), rel_mag(1), no_trend(), mtr_nd_ps())) {
    a <- rw_att(f6, tr)
    expect_true(a$results[[1]]$empty)
    expect_output(print(a), "EMPTY")
  }
})

test_that("issue 3: post = does not drop later post-periods' constraints", {
  d <- c(rw_data(), list(post2 = transform(rw_data()$post, yes = c(90, 100, 160, 190), no = 200 - c(90, 100, 160, 190))))
  f <- rw_fit(d)
  for (tr in list(smoothness(0.01), rel_mag(0.2))) {
    a <- rw_att(f, tr); b <- rw_att(f, tr, post = "post")
    expect_equal(c(b$results[[1]]$lower, b$results[[1]]$upper), c(a$results[[1]]$lower, a$results[[1]]$upper), tolerance = 1e-9)
  }
})

test_that("issue 4: E19(a)-(b) are computed only with one post-period", {
  d <- c(rw_data(), list(post2 = rw_data()$post))
  a <- rw_att(rw_fit(d), rel_mag(0.5))
  expect_true(is.na(a$results[[1]]$rm$vacuous_upper))
  expect_false(a$results[[1]]$rm$one_post)
  expect_output(print(a), "E19\\(a\\)-\\(b\\) not computed")
  a1 <- rw_att(rw_fit(), rel_mag(1))
  expect_false(is.na(a1$results[[1]]$rm$vacuous_upper))
})

test_that("issue 5: the allowance uses the representation of the model-implied form (R = 3, one-vs-rest)", {
  set.seed(11)
  p <- 60; X <- matrix(runif(p * 3), p); X <- X / rowSums(X); N <- sample(200:400, p, TRUE)
  A <- round(X * N); A[, 3] <- N - A[, 1] - A[, 2]
  yy <- function(m) rbinom(p, A[, 1], m[1]) + rbinom(p, A[, 2], m[2]) + rbinom(p, A[, 3], m[3])
  mk <- function(y) data.frame(u = 1:p, r1 = A[, 1], r2 = A[, 2], r3 = A[, 3], yes = y, no = N - y)
  f <- eiChange(list(e0 = mk(yy(c(.3, .5, .6))), e1 = mk(yy(c(.45, .5, .62)))), cbind(yes, no) ~ cbind(r1, r2, r3), unit = "u",
                row = "r1", versus = "rest", col = "yes", integer = "lp", model = fd())
  tg <- f$targets[[1]]
  rho1 <- eiPartialID:::.ed_rho(f$P, tg$reps[[1]]$terms)
  expect_equal(f$model$table$S_phi, 0.5 * sqrt(sum(rho1^2)), tolerance = 1e-12)
  expect_equal(eiPartialID:::.ed_model_target(f$P, f$model$mb, tg)$const, 0)
})

test_that("issue 6: numeric unit selections refer to data rows (an N = 0 unit first)", {
  dz <- lapply(rw_data(), function(d) rbind(data.frame(p = "Z", vet = 0, oth = 0, yes = 0, no = 0), d))
  fz <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                 targets = list(rate_at("vet", "yes", "post", where = 2)))
  fl <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                 targets = list(rate_at("vet", "yes", "post", where = "P1")))
  expect_equal(fz$table[, c("lower", "upper")], fl$table[, c("lower", "upper")])
  fk <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                 restrictions = list(known = list(row = "vet", col = "yes", election = "post", unit = 2, value = 10)))
  fkl <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                  restrictions = list(known = list(row = "vet", col = "yes", election = "post", unit = "P1", value = 10)))
  expect_equal(fk$table[, c("lower", "upper")], fkl$table[, c("lower", "upper")])
  expect_error(eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                        targets = list(rate_at("vet", "yes", "post", where = 1))), "N_i = 0|dropped")
})

test_that("issue 7: a character election variable needs 'elections'", {
  long <- do.call(rbind, lapply(names(rw_data()), function(e) cbind(rw_data()[[e]], el = e)))
  set.seed(2); long <- long[sample(nrow(long)), ]
  expect_error(eiChange(long, cbind(yes, no) ~ cbind(vet, oth), unit = "p", election = "el", base = "base", row = "vet",
                        col = "yes", integer = "lp"), "elections")
  f <- eiChange(long, cbind(yes, no) ~ cbind(vet, oth), unit = "p", election = "el", elections = c("early", "base", "post"),
                base = "base", row = "vet", col = "yes", integer = "lp")
  expect_equal(f$table[, c("lower", "upper")], rw_fit()$table[, c("lower", "upper")])
  long$el <- factor(long$el, levels = c("early", "base", "post"))
  expect_silent(eiChange(long, cbind(yes, no) ~ cbind(vet, oth), unit = "p", election = "el", base = "base", row = "vet",
                         col = "yes", integer = "lp"))
})

test_that("issue 8: smoothness with T = 2 imposes nothing and says so", {
  f2 <- eiChange(rw_data()[2:3], cbind(yes, no) ~ cbind(vet, oth), unit = "p", row = "vet", col = "yes", integer = "lp")
  a <- rw_att(f2, smoothness(0.05)); n <- rw_att(f2, no_trend())
  expect_equal(c(a$results[[1]]$lower, a$results[[1]]$upper), c(n$results[[1]]$lower, n$results[[1]]$upper))
  expect_output(print(a), "smoothness imposes nothing")
})

test_that("issue 9: model sign rows without 'sign' default to sign >= 0", {
  set.seed(3); p <- 80; x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
  mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = round(N * t), no = N - round(N * t))
  t0 <- pmin(pmax(.3 + .2 * x + rnorm(p, 0, .03), 0), 1); t1 <- pmin(pmax(t0 + .05 + rnorm(p, 0, .02), 0), 1)
  dl <- list(a = mk(t0), b = mk(t1))
  f1 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
                 model = fd(rows = list(list(type = "sign", sign_type = "dW"))))
  f2 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
                 model = fd(rows = list(list(type = "sign", sign_type = "dW", sign = 1))))
  expect_equal(f1$model$table$model_lower, f2$model$table$model_lower)
  expect_silent(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
                         model = lce_periods(rows = list(list(type = "sign", sign_type = "dB", from = "a", to = "b")))))
})

test_that("issues 10, 11, 13: bin_did note, LP cross-check on the E26 path, E25 status for grouped rungs", {
  f10 <- eiChange(rw_data()[2:3], cbind(yes, no) ~ cbind(vet, oth), unit = "p", integer = "lp",
                  targets = list(bin_did("yes", "base", "post", low = c("P1", "P2"), high = c("P3", "P4"))), model = fd())
  expect_match(f10$model$table$note, "binned DiD")
  f9 <- eiChange(rw_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes",
                 integer = "exact", restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
  fl <- rw_fit(restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
  expect_true(any(grepl("E26", f9$table$integer_status)))
  expect_equal(f9$table$lp_lower, fl$table$lower, tolerance = 1e-9)   # a real continuous LP, not the closed form
  expect_equal(f9$table$lp_upper, fl$table$upper, tolerance = 1e-9)
})

test_that("Def. 5.1 item 3: spillover is required only when the trend assumption uses a comparison group", {
  f <- rw_fit()
  a <- eiATT(f, treated = "vet", population = "fixed", trend = mtr_nd_ps(ps = FALSE), no_anticipation = TRUE)
  expect_equal(a$Lambda, 0)
  expect_equal(a$results[[1]]$lower, 0)
  expect_error(eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = mtr_nd_ps(), no_anticipation = TRUE),
               "spillover")
  expect_error(eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = parallel_trends(), no_anticipation = TRUE),
               "spillover")
})

test_that("E16(b): whole-person witness flag at the LP ends (integer margins)", {
  a <- rw_att(rw_fit(), parallel_trends())
  expect_true(is.logical(a$results[[1]]$witness) && length(a$results[[1]]$witness) == 2)
  # E16 [F]: one treated and one comparison person, both voting yes; S = 1/2 needs half a person
  d <- data.frame(u = c("A", "B"), tr = c(1, 0), co = c(0, 1), yes = c(1, 1), no = c(0, 0))
  f <- eiChange(list(e0 = d, e1 = d), cbind(yes, no) ~ cbind(tr, co), unit = "u", row = "tr", col = "yes", integer = "lp")
  b <- eiATT(f, treated = "tr", comparison = "co", population = "fixed", trend = parallel_trends(), spillover = 0.5, no_anticipation = TRUE)
  expect_equal(c(b$results[[1]]$lower, b$results[[1]]$upper), c(0, 0.5))
  expect_equal(b$results[[1]]$witness, c(lower = TRUE, upper = FALSE))     # the upper end needs theta = 1/2 of one person
  expect_output(print(b), "E16\\(b\\) whole-person witness")
})

# ---- R = 3 brute force of rows on changes (2c Def. 1, written here without the package's builders) ----
bf_R3 <- function(A, Y, rows) {
  # A: p x 3 group counts; Y: 2 x p yes totals (C = 2). Enumerate every integer table of both elections.
  p <- nrow(A); N <- rowSums(A); pw <- N / sum(N); X <- A / N
  tabs <- function(i, y) { g <- expand.grid(0:A[i, 1], 0:A[i, 2], 0:A[i, 3]); as.matrix(g[rowSums(g) == y, , drop = FALSE]) }
  per <- list(); for (t in 1:2) for (i in 1:p) per[[length(per) + 1]] <- tabs(i, Y[t, i])
  idx <- as.matrix(expand.grid(lapply(per, function(m) seq_len(nrow(m)))))
  beta <- function(k, t) t(vapply(1:p, function(i) per[[(t - 1) * p + i]][idx[k, (t - 1) * p + i], ] / pmax(A[i, ], 1), numeric(3)))
  mu <- colSums(pw * X); xr <- X[, 1]
  dB_pair <- function(b, s) sum(pw * xr * X[, s] * (b[, 1] - b[, s]))
  dB_W1 <- function(b) dB_pair(b, 2) + dB_pair(b, 3)
  dW_W2 <- function(b) sum(vapply(1:3, function(g) { xb <- sum(pw * X[, g] * xr) / mu[g]; sum(pw * X[, g] * (xr - xb) * b[, g]) }, 0))
  did <- function(k) {
    v <- 0
    for (t in 1:2) { n <- t(vapply(1:p, function(i) per[[(t - 1) * p + i]][idx[k, (t - 1) * p + i], ], numeric(3)))
      g <- sum(n[, 1]) / sum(A[, 1]) - sum(n[, 2:3]) / sum(A[, 2:3]); v <- v + if (t == 2) g else -g }
    v
  }
  vals <- c()
  for (k in seq_len(nrow(idx))) {
    b <- beta(k, 2) - beta(k, 1)
    ok <- switch(rows, pair = dB_pair(b, 3) >= -1e-12, W2 = dW_W2(b) >= -1e-12,
                 crW2 = (dB_W1(b) >= -1e-12 && dW_W2(b) >= -1e-12) || (dB_W1(b) <= 1e-12 && dW_W2(b) <= 1e-12))
    if (ok) vals <- c(vals, did(k))
  }
  if (!length(vals)) c(NA, NA) else range(vals)
}

test_that("R = 3 rows on changes: integer ends equal a brute force (pair dB, W2 dW, CR W2)", {
  set.seed(7)
  specs <- list(pair = list(sign_change = list(type = "dB", row = "g1", col = "yes", from = "e0", to = "e1", version = "pair", other = "g3", sign = 1)),
                W2 = list(sign_change = list(type = "dW", row = "g1", col = "yes", from = "e0", to = "e1", version = "W2", sign = 1)),
                crW2 = list(cr_change = list(row = "g1", col = "yes", from = "e0", to = "e1", version = "W2")))
  n_ok <- 0
  for (k in 1:6) {
    repeat { A <- matrix(sample(1:2, 6, TRUE), 2, 3); N <- rowSums(A); if (A[1, 1] / N[1] != A[2, 1] / N[2]) break }
    Y <- rbind(vapply(N, function(n) sample(0:n, 1), 0), vapply(N, function(n) sample(0:n, 1), 0))
    mk <- function(t) data.frame(u = c("U1", "U2"), g1 = A[, 1], g2 = A[, 2], g3 = A[, 3], yes = Y[t, ], no = N - Y[t, ])
    for (nm in names(specs)) {
      bf <- bf_R3(A, Y, nm)
      f <- eiChange(list(e0 = mk(1), e1 = mk(2)), cbind(yes, no) ~ cbind(g1, g2, g3), unit = "u", row = "g1", versus = "rest",
                    col = "yes", integer = "exact", restrictions = specs[[nm]])
      if (anyNA(bf)) { expect_true(is.na(f$table$lower[1]) || identical(f$table$integer_feasible[1], FALSE)); next }
      expect_equal(c(f$table$lower[1], f$table$upper[1]), bf, tolerance = 1e-9, info = paste(nm, k))
      n_ok <- n_ok + 1
    }
  }
  expect_gt(n_ok, 10)
})

test_that("code review round 2: R2-1 E29 CI widened by pi_S (h M + Lambda); N2-N3 wording", {
  set.seed(11); p <- 12
  aV <- sample(5:60, p, TRUE); aN <- sample(5:60, p, TRUE); N <- aV + aN
  y0 <- round(N * runif(p, .3, .6)); y1 <- pmin(N, y0 + round(N * runif(p, 0, .1)))
  mk <- function(y) data.frame(u = seq_len(p), v = aV, o = aN, yes = y, no = N - y)
  f <- eiChange(list(b = mk(y0), p = mk(y1)), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  A <- function(M, L) eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(M = M, target = "expected"),
                            spillover = L, no_anticipation = TRUE, alpha = 0.1)
  d0 <- A(0, 0)$results[[1]]$dw; a1 <- A(0.03, 0.02); d1 <- a1$results[[1]]$dw
  piS <- 1 - a1$wpt$design$pi_U
  expect_equal(d1$widen, piS * (1 * 0.03 + 0.02), tolerance = 1e-12)
  expect_equal(d1$ci, c(max(-1, d0$ci[1] - d1$widen), min(1, d0$ci[2] + d1$widen)), tolerance = 1e-12)
  expect_equal(d0$widen, 0)
  out <- capture.output(print(a1))
  expect_true(any(grepl("E29(d) widening w = pi_S (h M + Lambda)", out, fixed = TRUE)))
  # N2: no "tau_i = 0 sqrt" for the realised reading at M = 0
  o2 <- capture.output(print(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(),
                                   spillover = 0, no_anticipation = TRUE)))
  expect_false(any(grepl("0 sqrt", o2, fixed = TRUE)))
  expect_true(any(grepl("tau_i = 0", o2, fixed = TRUE)))
})
