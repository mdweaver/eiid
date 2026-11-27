# Probes for the code review of eiChange() / eiATT() (PHASE2e): edge cases, alignment, the model
# layer, E25, old outputs. Usage: Rscript phase2e_code_review_probes.R [lib_new] [lib_old]
# Each section prints one line per finding; see PHASE2e_code_review.md for the numbered issues.
args <- commandArgs(TRUE)
LIB_NEW <- if (length(args) >= 1) args[1] else .libPaths()[1]
LIB_OLD <- if (length(args) >= 2) args[2] else NA
.libPaths(c(LIB_NEW, .libPaths()))
suppressPackageStartupMessages(library(eiPartialID))
section <- function(s) cat("\n=====", s, "=====\n")

section("probe1")
local({
w_data <- function() {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10))
  yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  list(early = mk(1), base = mk(2), post = mk(3))
}
f0 <- eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp")
A <- function(f, tr, ...) tryCatch({ a <- eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE, ...)
  c(a$results[[1]]$lower, a$results[[1]]$upper) }, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("P1 no_trend, rate_limits=FALSE:", A(f0, no_trend(), rate_limits = FALSE), "\n")
a <- eiATT(f0, treated = "vet", population = "fixed", trend = mtr_nd_ps(nd = FALSE, ps = FALSE), spillover = 0, no_anticipation = TRUE, rate_limits = FALSE)
cat("P1b MTR only, rate_limits=FALSE:", a$results[[1]]$lower, a$results[[1]]$upper, "\n")
print(a)
# P2: smoothness with T = 2
d2 <- w_data()[2:3]
f2 <- eiChange(d2, cbind(yes, no) ~ cbind(vet, oth), unit = "p", row = "vet", col = "yes", integer = "lp")
cat("P2 SD with T=2:", A(f2, smoothness(0.05)), "\n")
cat("P2b RM with T=2:", A(f2, rel_mag(1)), "\n")
# P3: numeric unit selection with an N = 0 unit before the selected unit
dz <- lapply(w_data(), function(d) rbind(data.frame(p = "Z", vet = 0, oth = 0, yes = 0, no = 0), d))
fz <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
               targets = list(rate_at("vet", "yes", "post", where = 2)))       # row 2 of the data = P1
fr <- eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
               targets = list(rate_at("vet", "yes", "post", where = 1)))       # P1 without the N = 0 unit
cat("P3 rate of P1 via row number with an N=0 unit first:", fz$table$target, fz$table$lower, fz$table$upper, " vs direct:", fr$table$target, fr$table$lower, fr$table$upper, "\n")
fz2 <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                targets = list(rate_at("vet", "yes", "post", where = c("P1"))))
cat("P3 by label:", fz2$table$lower, fz2$table$upper, "\n")
r <- tryCatch(eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                targets = list(rate_at("vet", "yes", "post", where = vet > 100))), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("P3c where-expression with N=0 unit:", if (is.character(r)) r else c(r$table$lower, r$table$upper), "\n")
fk <- tryCatch(eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                restrictions = list(known = list(row = "vet", col = "yes", election = "post", unit = 2, value = 10))), error = function(e) conditionMessage(e))
if (!is.character(fk)) cat("P3d known with unit = 2 (data row 2 = P1):", fk$restrictions$restriction, "\n") else cat("P3d", fk, "\n")
# P4: row order differs across elections
d4 <- w_data(); d4$post <- d4$post[c(3, 1, 4, 2), ]
f4 <- eiChange(d4, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
               restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
f4r <- eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
                restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
cat("P4 shuffled order equal:", isTRUE(all.equal(f4$table[, c("lower", "upper")], f4r$table[, c("lower", "upper")])), "\n")
# P4b: long data frame with an election variable, unordered rows
long <- do.call(rbind, lapply(names(w_data()), function(e) cbind(w_data()[[e]], el = e)))
long <- long[sample(nrow(long)), ]
f4l <- eiChange(long, cbind(yes, no) ~ cbind(vet, oth), unit = "p", election = "el", elections = c("early", "base", "post"), base = "base", row = "vet", col = "yes", integer = "lp")
cat("P4b long format equal:", isTRUE(all.equal(f4l$table[, c("lower", "upper")], f0$table[, c("lower", "upper")])), "\n")
f4l2 <- tryCatch(eiChange(long, cbind(yes, no) ~ cbind(vet, oth), unit = "p", election = "el", base = "base", row = "vet", col = "yes", integer = "lp"),
                 error = function(e) conditionMessage(e))
cat("P4c long format without 'elections' (character election variable):", if (is.character(f4l2)) f4l2 else paste(f4l2$P$elections, collapse = ","), "\n")
# P5: N = 0 units dropped: same answer
f5 <- eiChange(dz, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp")
cat("P5 N=0 dropped equal:", isTRUE(all.equal(f5$table[, c("lower", "upper")], f0$table[, c("lower", "upper")])), "\n")
# P6: infeasible rows
f6 <- eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp",
               restrictions = list(known = list(row = "vet", col = "yes", election = "post", unit = "P1", value = 11)))
cat("P6 status:", f6$table$integer_status, "\n")
for (tr in list(parallel_trends(), rel_mag(1), no_trend(), mtr_nd_ps())) { a6 <- tryCatch(eiATT(f6, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE), error = function(e) conditionMessage(e)); cat("P6", tr$kind, if (is.character(a6)) paste("ERROR", a6) else a6$results[[1]]$lower, "\n") }

# P7: the event-study corner when the first two default targets are on the same side of the base
d7 <- c(w_data(), list(post2 = w_data()$post))
f7 <- eiChange(d7, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp")
cat("P7 corner targets:", f7$event_study$corner$targets, "attainable:", f7$event_study$corner$margins_only, "\n")
# P8: SD with posts restricted (post = first only) vs all posts
d8 <- c(w_data(), list(post2 = transform(w_data()$post, yes = c(90, 100, 160, 190), no = 200 - c(90, 100, 160, 190))))
f8 <- eiChange(d8, cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "lp")
for (tr in list(smoothness(0.01), rel_mag(0.2))) {
  a8 <- eiATT(f8, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE)
  a8p <- eiATT(f8, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE, post = "post")
  cat("P8", tr$kind, "ATT_post with all posts:", a8$results[[1]]$lower, a8$results[[1]]$upper, " with post = 'post':", a8p$results[[1]]$lower, a8p$results[[1]]$upper, "\n")
}
# P9: integer = "exact" with W1-only rows uses the E26 lattice
f9 <- eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes", integer = "exact",
               restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
cat("P9 exact:", f9$table$integer_status, "\n")
# P10: FD model, bin_did target note
f10 <- eiChange(w_data()[2:3], cbind(yes, no) ~ cbind(vet, oth), unit = "p", integer = "lp",
                targets = list(bin_did("yes", "base", "post", low = c("P1", "P2"), high = c("P3", "P4"))), model = fd())
cat("P10 bin_did under FD:", f10$model$table$note, "\n")
})

section("probe2")
local({
set.seed(4); found <- 0
for (k in 1:60) {
  p <- 4; N <- sample(30:80, p, TRUE); a1 <- round(N * runif(p, .1, .9))
  mk <- function(y) data.frame(u = paste0("U", 1:p), v = a1, o = N - a1, yes = y, no = N - y)
  ys <- lapply(1:4, function(t) round(N * runif(p, if (t == 4) .6 else 0, 1)))
  dl <- setNames(lapply(ys, mk), c("pre", "base", "p1", "p2"))
  f <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes", integer = "lp",
                restrictions = if (k %% 2) list(cr_change = list(row = "v", col = "yes", from = "pre", to = "base")) else NULL)
  for (tr in list(smoothness(0.02), rel_mag(0.5))) {
    a <- eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE)
    b <- eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE, post = "p1")
    x <- c(a$results[[1]]$lower, a$results[[1]]$upper); y <- c(b$results[[1]]$lower, b$results[[1]]$upper)
    if (!isTRUE(all.equal(x, y, tolerance = 1e-8))) { found <- found + 1; if (found <= 3) cat(k, tr$kind, "all posts:", format(x, digits = 5), " post='p1':", format(y, digits = 5), "\n") }
  }
}
cat("designs where restricting post changes ATT_p1:", found, "\n")
})

section("probe3")
local({
set.seed(11)
p <- 300; X <- matrix(runif(p * 3), p); X <- X / rowSums(X); N <- sample(200:400, p, TRUE)
A <- round(X * N); A[, 3] <- N - A[, 1] - A[, 2]; X <- A / N
m0 <- c(.3, .5, .6); m1 <- c(.45, .5, .62)                 # true group rates per election (no contextual effect)
yy <- function(m) { y <- rbinom(p, A[, 1], m[1]) + rbinom(p, A[, 2], m[2]) + rbinom(p, A[, 3], m[3]); y }
mk <- function(y) data.frame(u = 1:p, r1 = A[, 1], r2 = A[, 2], r3 = A[, 3], yes = y, no = N - y)
dl <- list(e0 = mk(yy(m0)), e1 = mk(yy(m1)))
f <- eiChange(dl, cbind(yes, no) ~ cbind(r1, r2, r3), unit = "u", row = "r1", versus = "rest", col = "yes", integer = "lp", model = fd())
tg <- f$targets[[1]]
P <- f$P
rho <- lapply(tg$reps, function(rp) eiPartialID:::.ed_rho(P, rp$terms))
S <- vapply(rho, function(r) 0.5 * sqrt(sum(r^2)), 0)
cat("Q1 representations S_phi:", S, " used S_phi:", f$model$table$S_phi, "\n")
mb <- f$model$mb
fm <- eiPartialID:::.ed_model_target(P, mb, tg)
cat("Q1 model target form constant (rep 1 has 0):", fm$const, "  rep-2 constant (outcome-based):", as.double(tg$reps[[2]]$const), "\n")
# the model-implied value under each representation at the true rates (no contextual effects: Delta m = m1 - m0)
dm <- m1 - m0; Ar <- sum(A[, 1]); Nn <- sum(N)
phi1 <- dm[1] - sum(A[, 2] * dm[2] + A[, 3] * dm[3]) / (Nn - Ar)
dY <- sum(dl$e1$yes) - sum(dl$e0$yes)
phi2 <- Nn / (Ar * (Nn - Ar)) * sum(A[, 1] * dm[1]) - dY / (Nn - Ar)
cat(sprintf("Q1 phi~(rep 1) = %.5f, phi~(rep 2) = %.5f, difference %.5f; allowance used %.5f\n", phi1, phi2, phi2 - phi1, f$model$table$allowance))
# Monte Carlo: actual DiD in the actual-target interval?
cov <- c(0, 0); nrep <- 150
for (b in seq_len(nrep)) {
  n0 <- cbind(rbinom(p, A[, 1], m0[1]), rbinom(p, A[, 2], m0[2]), rbinom(p, A[, 3], m0[3]))
  n1 <- cbind(rbinom(p, A[, 1], m1[1]), rbinom(p, A[, 2], m1[2]), rbinom(p, A[, 3], m1[3]))
  d <- list(e0 = mk(rowSums(n0)), e1 = mk(rowSums(n1)))
  ff <- eiChange(d, cbind(yes, no) ~ cbind(r1, r2, r3), unit = "u", row = "r1", versus = "rest", col = "yes", integer = "lp", model = fd(ndraw = 50))
  act <- (sum(n1[, 1]) - sum(n0[, 1])) / Ar - (sum(n1[, 2:3]) - sum(n0[, 2:3])) / (Nn - Ar)
  m <- ff$model$table
  cov[1] <- cov[1] + (act >= m$model_lower & act <= m$model_upper)
}
cat("Q1 MC coverage of the actual DiD by the actual-target interval (FD, R=3 rest):", cov[1] / nrep, "\n")
})

section("probe4")
local({
sim3 <- function(p = 200, seed = 3) {
  set.seed(seed)
  x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
  t0 <- .3 + .2 * x + rnorm(p, 0, .03); t1 <- t0 + .02 + .1 * x + rnorm(p, 0, .02); t2 <- t1 + .05 + .15 * x + rnorm(p, 0, .02)
  mkT <- function(l) round(N * pmin(pmax(l, 0), 1))
  mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = mkT(t), no = N - mkT(t))
  list(m1 = mk(t0), b = mk(t1), p = mk(t2))
}
dl <- sim3()
fF <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
               model = fd(rows = list(list(type = "R1"))))
mb <- fF$model$mb
for (k in seq_along(mb$blocks)) {
  idx <- mb$boff[k] + seq_len(mb$blocks[[k]]$nbeta)
  cat("H1 block", k, "joint V block == own V:", isTRUE(all.equal(mb$V[idx, idx], mb$blocks[[k]]$V, tolerance = 1e-10)),
      " beta:", isTRUE(all.equal(mb$beta[idx], mb$blocks[[k]]$beta)), "\n")
}
cat("H1 stacked V positive definite: min eigen", min(eigen(mb$V, only.values = TRUE)$values), "\n")
# overlapping / non-consecutive pairs refused
r <- tryCatch(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
              model = fd(pairs = list(c("m1", "b"), c("b", "p"), c("m1", "p")))), error = function(e) conditionMessage(e))
cat("H1 overlapping pairs:", substr(if (is.character(r)) r else "accepted", 1, 90), "\n")
for (tr in list(trend_dial(0.02), rel_mag(1), smoothness(0.02))) {
  a <- eiATT(fF, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = 0.01, no_anticipation = TRUE, pretrends = "model")
  cat("H2/H3 FD", tr$kind, ": outer", a$model$outer, " nq", length(a$model$radii), " alpha_proj", a$model$alpha_proj, " alpha_q", a$model$alpha_q,
      " sum", a$model$alpha_proj + length(a$model$radii) * a$model$alpha_q, " radii", paste(names(a$model$radii), collapse = ","),
      " set", format(c(a$results[[1]]$lower, a$results[[1]]$upper), digits = 4), "\n")
}
fL <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp", model = lce_periods())
for (tr in list(trend_dial(0.02), rel_mag(1))) {
  a <- eiATT(fL, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE, pretrends = "model")
  cat("H2/H3 LCE", tr$kind, ": outer", a$model$outer, " radii", paste(names(a$model$radii), collapse = ","), " set", format(c(a$results[[1]]$lower, a$results[[1]]$upper), digits = 4), "\n")
}
r <- tryCatch(eiATT(fF, treated = "v", comparison = "o", population = "fixed", trend = mtr_nd_ps(), spillover = 0, no_anticipation = TRUE, pretrends = "model"), error = function(e) conditionMessage(e))
cat("MTR model:", substr(r, 1, 80), "\n")
r <- tryCatch(eiATT(fF, treated = "v", comparison = "o", population = "fixed", trend = no_trend(), spillover = 0, no_anticipation = TRUE, pretrends = "model"), error = function(e) conditionMessage(e))
cat("none model:", substr(r, 1, 80), "\n")
# grouped LCE rung
g <- rep(1:4, length.out = nrow(dl$b))
fG <- tryCatch(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
               model = lce_periods(model = "M2", group = g, rows = list(list(type = "sign", sign_type = "dW", from = "b", to = "p")))), error = function(e) conditionMessage(e))
if (is.character(fG)) cat("grouped:", fG, "\n") else cat("grouped status:", unique(fG$model$table$status), "\n")
# FD with a model sign row (deviation 11) and R4/delta rows; status
fS <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
               model = fd(rows = list(list(type = "sign", sign_type = "dW", sign = 1), list(type = "delta", delta = 0.05))))
cat("FD + sign + delta status:", unique(fS$model$table$status), "\n")
print(fS$model$table[, c("target", "plugin_lower", "plugin_upper", "model_lower", "model_upper", "lower", "upper")])
# E8 sharp, 2x2 one pair
f8 <- eiChange(dl[2:3], cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp", model = fd(sharp = TRUE))
print(f8$model$extra[[1]]$e8)
# coarse design with repeated compositions: E8 sharp vs outer
set.seed(9); p <- 200; x <- sample(c(.2, .4, .6, .8), p, TRUE); N <- rep(100, p); a1 <- round(N * x)
t0 <- pmin(pmax(.3 + .2 * x + rnorm(p, 0, .15), 0), 1); t1 <- pmin(pmax(t0 + .05 + rnorm(p, 0, .1), 0), 1)
mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = round(N * t), no = N - round(N * t))
f9 <- eiChange(list(a = mk(t0), b = mk(t1)), cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp", model = fd(sharp = TRUE))
print(f9$model$extra[[1]]$e8)
})

section("probe7")
local({
set.seed(3); p <- 120; x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = round(N * t), no = N - round(N * t))
t0 <- pmin(pmax(.3 + .2 * x + rnorm(p, 0, .03), 0), 1); t1 <- pmin(pmax(t0 + .05 + rnorm(p, 0, .02), 0), 1)
dl <- list(a = mk(t0), b = mk(t1))
r <- tryCatch(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
              model = fd(rows = list(list(type = "sign", sign_type = "dW")))), error = function(e) conditionMessage(e))
cat("FD sign row without 'sign':", if (is.character(r)) r else "ok", "\n")
r <- tryCatch(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
              model = lce_periods(rows = list(list(type = "sign", sign_type = "dB", from = "a", to = "b")))), error = function(e) conditionMessage(e))
cat("LCE sign row without 'sign':", if (is.character(r)) r else "ok", "\n")
r <- tryCatch(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
              model = lce_periods(rows = list(list(type = "sign", from = "a", to = "b")))), error = function(e) conditionMessage(e))
cat("LCE sign row without sign_type and sign:", if (is.character(r)) r else "ok", "\n")
# Def. 5.1 corner cases
f <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp")
r <- tryCatch(eiATT(f, treated = "v", population = "fixed", trend = mtr_nd_ps(ps = FALSE), no_anticipation = TRUE), error = function(e) conditionMessage(e))
cat("MTR+ND, no comparison, no spillover stated:", substr(if (is.character(r)) r else "accepted", 1, 120), "\n")
r <- tryCatch(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = "none", spillover = 0, no_anticipation = TRUE), error = function(e) conditionMessage(e))
cat("trend = 'none' string:", if (is.character(r)) r else "accepted", "\n")
r <- tryCatch(eiATT(f, treated = "v", comparison = "v", population = "fixed", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE), error = function(e) conditionMessage(e))
cat("comparison = treated:", if (is.character(r)) r else "accepted", "\n")
r <- tryCatch(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = parallel_trends(), spillover = -0.1, no_anticipation = TRUE), error = function(e) conditionMessage(e))
cat("negative spillover:", if (is.character(r)) r else "accepted", "\n")
fs <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp", population = "stayers")
a <- eiATT(fs, treated = "v", comparison = "o", population = "stayers", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE)
cat(grep("population", capture.output(print(a)), value = TRUE), sep = "\n")
})

section("probe8")
local({
set.seed(21); bad <- 0; nj <- 0; badj <- 0
for (k in 1:15) {
  p <- 5; N <- sample(20:60, p, TRUE); a1 <- pmax(1, pmin(N - 1, round(N * runif(p)))) 
  ys <- lapply(1:3, function(t) sample(0:max(N), p, TRUE) %% (N + 1))
  mk <- function(y, s = 1) data.frame(u = paste0("U", 1:p), v = a1 * s, o = (N - a1) * s, yes = y * s, no = (N - y) * s)
  dl <- setNames(lapply(ys, mk), c("a", "b", "c")); dh <- setNames(lapply(ys, mk, s = 0.37), c("a", "b", "c"))
  rows <- list(cr_change = list(row = "v", col = "yes", from = "b", to = "c", slack = 0.03),
               order_change = list(row = "v", col = "yes", from = "a", to = "b", sign = -1, slack = 0.2))
  f1 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp", restrictions = rows)
  f2 <- eiChange(dh, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", restrictions = rows)
  if (!isTRUE(all.equal(f1$table[, c("lower", "upper")], f2$table[, c("lower", "upper")], tolerance = 1e-7))) bad <- bad + 1
  # eiJoint vs the interval criterion of E3(c) (margins only)
  f0 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  G <- f0$event_study$G
  for (j in 1:20) {
    v <- c(runif(1, f0$table$dd_lower[1], f0$table$dd_upper[1]), runif(1, f0$table$dd_lower[2], f0$table$dd_upper[2]))
    lo <- max(G$G_lower[2], G$G_lower[1] - v[1], G$G_lower[3] - v[2]); hi <- min(G$G_upper[2], G$G_upper[1] - v[1], G$G_upper[3] - v[2])
    crit <- lo <= hi
    nj <- nj + 1; if (crit != eiJoint(f0, v, rows = FALSE)) badj <- badj + 1
  }
}
cat("non-integer margins (x0.37) vs integer, LP sets equal: failures", bad, "of 15\n")
cat("eiJoint vs interval criterion of E3(c): failures", badj, "of", nj, "\n")
})

section("probe9")
local({
set.seed(8); mis <- 0; n <- 0; nna <- 0
for (k in 1:40) {
  p <- 3; N <- sample(30:60, p, TRUE); a1 <- round(N * runif(p, .2, .8))
  mk <- function(y) data.frame(u = paste0("U", 1:p), v = a1, o = N - a1, yes = y, no = N - y)
  ys <- list(round(N * runif(p, .3, .5)), round(N * runif(p, .3, .5)), round(N * runif(p, 0, .3)), round(N * runif(p, .8, 1)))
  dl <- setNames(lapply(ys, mk), c("pre", "base", "p1", "p2"))
  rows <- if (k %% 2) list(cr_change = list(row = "v", col = "yes", from = "pre", to = "base")) else NULL
  f <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes", integer = "lp", restrictions = rows)
  a <- eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(runif(1, .05, .6)), spillover = 0, no_anticipation = TRUE)
  e <- a$results[[1]]
  if (is.na(e$lower)) next
  n <- n + 1
  vu <- abs(e$upper - e$none$upper) < 1e-7; vl <- abs(e$lower - e$none$lower) < 1e-7
  if (is.na(e$rm$vacuous_upper)) { nna <- nna + 1; next }   # round 2: not computed with several posts (issue 4 fix)
  if (vu != e$rm$vacuous_upper || vl != e$rm$vacuous_lower) { mis <- mis + 1; if (mis <= 3) cat("k", k, "set", e$lower, e$upper, "none", e$none$lower, e$none$upper, "criterion", e$rm$vacuous_lower, e$rm$vacuous_upper, "\n") }
}
cat("E19(b) criterion with two post-periods: mismatches", mis, "of", n, "; not computed (NA):", nna, "\n")
})

section("probe10")
local({
set.seed(3); p <- 150; x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
t0 <- .3 + .2 * x + rnorm(p, 0, .03); t1 <- t0 + .02 + .1 * x + rnorm(p, 0, .02)
mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = round(N * t), no = N - round(N * t))
dl <- list(a = mk(t0), b = mk(t1))
f <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
              restrictions = list(total = list(row = "v", col = "yes", election = "a", value = sum(round(a1 * t0)))),
              model = fd(rows = list(list(type = "didS"))))
print(f)
f2 <- eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", row = "v", col = "yes", integer = "lp",
              restrictions = list(cr_change = list(row = "v", col = "yes", from = "a", to = "b")), model = fd(rows = list(list(type = "didS"))))
cat(tail(capture.output(print(f2)), 8), sep = "\n")
})

section("probe5: performance, 500 units x 3 elections")
if (!identical(Sys.getenv("SKIP_PERF"), "1")) local({
set.seed(1); p <- 500
x <- runif(p, .05, .95); N <- sample(200:900, p, TRUE); a1 <- round(N * x)
t0 <- plogis(qlogis(.4) + .5 * x + rnorm(p, 0, .2)); t1 <- plogis(qlogis(t0) + .05 + rnorm(p, 0, .1)); t2 <- plogis(qlogis(t1) + .2 * x + rnorm(p, 0, .1))
mk <- function(t) { y <- round(N * t); data.frame(u = seq_len(p), v = a1, o = N - a1, yes = y, no = N - y) }
dl <- list(pre = mk(t0), base = mk(t1), post = mk(t2))
tm <- function(expr) { s <- proc.time()[3]; v <- expr; cat(sprintf("  %.1fs\n", proc.time()[3] - s)); v }
cat("margins only, integer auto:"); f0 <- tm(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes"))
print(f0$table[, c("target", "lower", "upper", "integer_status")])
cat("CR on both changes, integer lp:"); f1 <- tm(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes", integer = "lp",
  restrictions = list(cr_change = list(list(row = "v", col = "yes", from = "pre", to = "base"), list(row = "v", col = "yes", from = "base", to = "post")))))
cat("CR on both changes, integer auto:"); f1a <- tm(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes",
  restrictions = list(cr_change = list(list(row = "v", col = "yes", from = "pre", to = "base"), list(row = "v", col = "yes", from = "base", to = "post")))))
print(f1a$table[, c("target", "lower", "upper", "integer_status")])
cat("eiATT RM on CR fit:"); a <- tm(eiATT(f1, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(1), spillover = 0, no_anticipation = TRUE))
cat("eiATT SD:"); a <- tm(eiATT(f1, treated = "v", comparison = "o", population = "fixed", trend = smoothness(.02), spillover = 0, no_anticipation = TRUE))
cat("order_change in every unit, integer auto (TU E27?):"); f2 <- tm(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes",
  restrictions = list(order_change = list(row = "v", col = "yes", from = "base", to = "post"))))
print(f2$table[, c("target", "lower", "upper", "integer_status", "tu_class")])
cat("FD + didS model layer:"); f3 <- tm(eiChange(dl, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "base", row = "v", col = "yes", integer = "lp", model = fd(rows = list(list(type = "didS")))))
cat("E25 RM model:"); a3 <- tm(eiATT(f3, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(1), spillover = 0, no_anticipation = TRUE, pretrends = "model"))
})

