# eiATT in expectation (2e Sec. 5.8-5.9): within-precinct and pooled parallel trends (E28), the realised
# target D_W and the Hoeffding allowance (E29), baseline / post / change LCE with the curvature allowance
# (E30), the local dial (E33(a)); plus the RM hole and horizon examples (E18) and the reporting labels.
# Reference values: docs/extension/reviews/phase2e_independent_checks_pt.py (published to 2 decimals
# of a percentage point, hence tolerance 5.1e-5 in rate units) and 2e Sec. 5.9(c) (exact fractions).

mk_mean <- function(size, vet, tot) data.frame(u = paste0("P", seq_along(vet)), v = vet, o = size - vet, yes = tot, no = size - tot)
fit_mean <- function(size, vet, pre, post_states) {
  post_states <- as.matrix(post_states)
  eiChange(list(b = mk_mean(size, vet, pre), p = mk_mean(size, vet, rowMeans(post_states))), cbind(yes, no) ~ cbind(v, o),
           unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
}
att_mean <- function(f, tr, lce = NULL, law = NULL, L = 0, ...)
  eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = L, no_anticipation = TRUE,
        lce = lce, law = law, ...)
iv <- function(a, k = 1) c(a$results[[k]]$lower, a$results[[k]]$upper)

rare <- list(size = 20000, vet = c(5000, 8000, 11000, 14000, 17000), pre = c(40, 12000, 12000, 12000, 12000),
             post = rbind(c(50, 30), c(9964, 9956), c(8188, 8192), c(6682, 6698), c(5446, 5474)))
ordinary <- list(size = 1000, vet = c(100, 300, 500, 700, 900), pre = rep(500, 5),
                 post = rbind(c(682, 650), c(606, 590), c(530, 530), c(454, 470), c(378, 410)))
three <- list(size = 800, vet = c(200, 400, 600), pre = c(4, 164, 324), post = matrix(c(204, 124, 4), 3, 1))

test_that("E28: the five-precinct WPT/PPT table (known aggregate law, FD-LCE on changes)", {
  for (ex in list(list(d = rare, w = c(10.39, 12.25), p = c(-3.96, -2.22)),
                  list(d = ordinary, w = c(-67.33, 43.78), p = c(-80.56, 52.13)))) {
    f <- fit_mean(ex$d$size, ex$d$vet, ex$d$pre, ex$d$post)
    law <- list(p = ex$d$post)
    expect_true(all(abs(100 * iv(att_mean(f, within_pt(), lce_mean(change = TRUE), law)) - ex$w) <= 0.0051))
    expect_true(all(abs(100 * iv(att_mean(f, pooled_pt(), lce_mean(change = TRUE), law)) - ex$p) <= 0.0051))
  }
})

test_that("E30: the three-precinct table, exact fractions, and the curvature breakdown 1/600", {
  f <- fit_mean(three$size, three$vet, three$pre, three$post)
  W <- function(l) iv(att_mean(f, within_pt(), l))
  expect_equal(W(lce_mean(change = FALSE)), c(-29 / 50, 41 / 150), tolerance = 1e-9)
  expect_equal(W(lce_mean(change = TRUE)), c(-37 / 75, 11 / 150), tolerance = 1e-9)
  for (l in list(lce_mean(base = TRUE), lce_mean(base = TRUE, post = TRUE, change = FALSE), lce_mean(post = TRUE),
                 lce_mean(base = TRUE, post = TRUE)))
    expect_equal(W(l), c(2 / 75, 11 / 150), tolerance = 1e-9)
  a <- att_mean(f, within_pt(), lce_mean(base = TRUE))
  expect_equal(a$results[[1]]$breakdown_eps, 1 / 600, tolerance = 1e-9)
  expect_true(W(lce_mean(base = TRUE, epsilon = 1 / 600 - 1e-6))[1] > 0)                  # strictly positive below eps*
  expect_equal(W(lce_mean(base = TRUE, epsilon = 1 / 600))[1], 0, tolerance = 1e-9)       # zero first feasible at eps*
  eps <- c(0, 0.001, 0.01, 0.1, 0.3)
  ivs <- lapply(eps, function(e) W(lce_mean(base = TRUE, epsilon = e)))
  for (j in 2:length(eps)) expect_true(ivs[[j]][1] <= ivs[[j - 1]][1] + 1e-12 && ivs[[j]][2] >= ivs[[j - 1]][2] - 1e-12)
  expect_equal(W(lce_mean(base = TRUE, epsilon = 0.5)), W(lce_mean(change = TRUE)), tolerance = 1e-9)  # eps >= 1/2: change only
  expect_output(print(a), "epsilon\\* = 0.001667 \\(the sign holds strictly below")
})

test_that("E33(a): local dial M nested, exact WPT at 0, one-LP breakdown; either setting varied with the other fixed", {
  f <- fit_mean(three$size, three$vet, three$pre, three$post)
  W <- function(M, e = 0) iv(att_mean(f, within_pt(M = M), lce_mean(base = TRUE, epsilon = e)))
  expect_equal(W(0), c(2 / 75, 11 / 150), tolerance = 1e-9)
  Ms <- c(0, 0.01, 0.05, 0.2)
  for (j in 2:4) { a <- W(Ms[j]); b <- W(Ms[j - 1]); expect_true(a[1] <= b[1] + 1e-12 && a[2] >= b[2] - 1e-12) }
  a <- att_mean(f, within_pt(), lce_mean(base = TRUE))
  Ms <- a$results[[1]]$breakdown_M
  expect_true(is.finite(Ms) && Ms > 0)
  expect_equal(W(Ms)[1], 0, tolerance = 1e-9)
  expect_true(W(Ms * 0.99)[1] > 0)
  # the curvature breakdown with M held at a positive value is smaller than with M = 0
  a2 <- att_mean(f, within_pt(M = Ms / 2), lce_mean(base = TRUE))
  expect_true(a2$results[[1]]$breakdown_eps < 1 / 600)
  # M* with epsilon held fixed > 0 is smaller than at epsilon = 0
  a3 <- att_mean(f, within_pt(), lce_mean(base = TRUE, epsilon = 1 / 1200))
  expect_true(a3$results[[1]]$breakdown_M < Ms)
})

test_that("E28(b)-(c): PPT witness uses one common untreated rate where the local imputation is not a rate", {
  # veteran baselines (.1, .9), comparison changes (-.3, +.3) pinned by the bounds (2e E28 [F], round 2)
  des <- list(p = 2, E = 2, K = 1, aV = c(50, 50), aN = c(50, 50), x = c(.5, .5), wV = c(.5, .5), wN = c(.5, .5),
              sup = c(TRUE, TRUE), Tm = cbind(c(.3, .7), c(.25, .75)),
              lo = cbind(c(.1, .9), c(.3, .7)), hi = cbind(c(.1, .9), c(.3, .7)))
  lce0 <- lce_mean(change = FALSE)
  w <- eiPartialID:::.wpt_lp(des, "wpt", lce0, 0, 0)
  expect_null(w$solve(1, -1))                      # WPT refuted: .1 - .3 = -.2 is not a rate
  p <- eiPartialID:::.wpt_lp(des, "ppt", lce0, 0, 0)
  s <- p$solve(1, -1)
  expect_equal(s$value, 0, tolerance = 1e-12)
  expect_equal(s$p0, 0.5, tolerance = 1e-12)       # the common untreated rate p0 lies in [0, 1]
})

test_that("E28(d): unsupported precincts are bounded with theta in [0, 1]; pi_U reported; a_min threshold", {
  vet <- c(20, 50, 80, 100); N <- 100
  d0 <- data.frame(u = paste0("P", 1:4), v = vet, o = N - vet, yes = c(30, 45, 60, 70), no = N - c(30, 45, 60, 70))
  d1 <- data.frame(u = paste0("P", 1:4), v = vet, o = N - vet, yes = c(35, 55, 62, 80), no = N - c(35, 55, 62, 80))
  f <- eiChange(list(b = d0, p = d1), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  a <- att_mean(f, within_pt())
  w <- a$wpt$design
  expect_equal(w$pi_U, 100 / 250)
  expect_equal(a$wpt$support$supported, c(TRUE, TRUE, TRUE, FALSE))
  fs <- eiChange(list(b = d0[1:3, ], p = d1[1:3, ]), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  sub <- iv(att_mean(fs, within_pt()))
  piS <- 150 / 250; shift <- (100 / 250) * 0.80                 # E[T] of the all-treated precinct at the post election
  expect_equal(iv(a), c(piS * sub[1] + shift - 0.4, piS * sub[2] + shift), tolerance = 1e-9)
  # a_min: precinct 3 (20 comparison persons) becomes unsupported at a_min = 25
  a25 <- att_mean(f, within_pt(a_min = 25))
  expect_equal(a25$wpt$design$pi_U, 180 / 250)
  expect_output(print(a25), "P3 20 \\(unsupported\\)")
  expect_output(print(a25), "pi_U = 0.72")
})

test_that("E29: D_W from the count layer, the Hoeffding allowance and its coverage statement", {
  f <- fit_mean(three$size, three$vet, three$pre, three$post)
  a <- att_mean(f, within_pt(), alpha = 0.1)
  d <- a$results[[1]]$dw; w <- a$wpt$design
  rho <- 4 * w$wV
  expect_equal(d$r_alpha, sqrt(0.5 * log(2 / 0.1) * sum(rho^2)), tolerance = 1e-12)
  expect_equal(d$ci, c(max(-1, d$LU[1] - w$pi_U - d$r_alpha), min(1, d$LU[2] + d$r_alpha)))
  ap <- att_mean(f, within_pt(), alpha = 0.1, independence = "persons")
  rp <- sqrt(0.5 * log(2 / 0.1) * sum(4 * w$wV^2 * (1 / w$aV + 1 / w$aN)))
  expect_equal(ap$results[[1]]$dw$r_alpha, rp, tolerance = 1e-12)
  expect_true(rp < d$r_alpha)
  expect_equal(w$n_eff, 1 / sum(w$wV^2))
  # D_W by brute force: the count layer's range of the realised target (margins only: per-precinct boxes)
  aV <- three$vet; aN <- 800 - aV
  lo <- function(T) pmax(0, T - aN); hi <- function(T) pmin(aV, T)
  wv <- aV / sum(aV)
  # term_i = (yV1 - yV0)/aV - ((T1 - yV1) - (T0 - yV0))/aN, increasing in yV1, decreasing in yV0
  tl <- (lo(three$post[, 1]) - hi(three$pre)) * (1 / aV + 1 / aN) - (three$post[, 1] - three$pre) / aN
  th <- (hi(three$post[, 1]) - lo(three$pre)) * (1 / aV + 1 / aN) - (three$post[, 1] - three$pre) / aN
  expect_equal(d$LU, c(sum(wv * tl), sum(wv * th)), tolerance = 1e-9)
  # belief rows: the coverage statement carries the penalty
  fc <- eiChange(list(b = mk_mean(800, three$vet, three$pre), p = mk_mean(800, three$vet, three$post[, 1])),
                 cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp",
                 restrictions = list(cr_change = list(row = "v", col = "yes", from = "b", to = "p")))
  expect_output(print(att_mean(fc, within_pt())), "P\\(the belief rows are false\\)")
})

test_that("reporting rule: every printed interval says what it is", {
  f <- fit_mean(three$size, three$vet, three$pre, three$post)
  a <- att_mean(f, within_pt(), lce_mean(base = TRUE))
  out <- paste(capture.output(print(a)), collapse = "\n")
  expect_match(out, "[IDENTIFICATION BOUND, mean model (E28/E30)]", fixed = TRUE)
  expect_match(out, "[DESCRIPTIVE, realised counts] D_W", fixed = TRUE)
  expect_match(out, "[CONFIDENCE INTERVAL, one observed election pair (E29(c))]", fixed = TRUE)
  expect_match(out, "does not make the identification bound above a confidence interval", fixed = TRUE)
  expect_match(out, "observed returns (known-law reading", fixed = TRUE)
  expect_match(out, "n_eff", fixed = TRUE); expect_match(out, "pi_U", fixed = TRUE)
  expect_match(out, "local comparison counts", fixed = TRUE); expect_match(out, "independence unit: precincts", fixed = TRUE)
  ap <- att_mean(f, pooled_pt(), lce_mean(change = TRUE))
  op <- paste(capture.output(print(ap)), collapse = "\n")
  expect_match(op, "[IDENTIFICATION BOUND, mean model (E28/E30)]", fixed = TRUE)
  expect_match(op, "[DESCRIPTIVE, realised counts] pooled DiD", fixed = TRUE)
  expect_false(grepl("[CONFIDENCE INTERVAL", op, fixed = TRUE))
  ar <- eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = parallel_trends(), spillover = 0, no_anticipation = TRUE)
  orr <- paste(capture.output(print(ar)), collapse = "\n")
  expect_match(orr, "[IDENTIFICATION SET, realised ATT, count layer (Thm E16)]", fixed = TRUE)
  expect_match(orr, "[DESCRIPTIVE, realised counts]", fixed = TRUE)
  # every interval line carries a bracketed label
  ivlines <- grep("\\[-?[0-9.]+, -?[0-9.]+\\]", strsplit(out, "\n")[[1]], value = TRUE)
  expect_true(all(grepl("^\\s*\\[(IDENTIFICATION|DESCRIPTIVE|CONFIDENCE)", ivlines)))
  ivr <- grep("\\[-?[0-9.]+, -?[0-9.]+\\]", strsplit(orr, "\n")[[1]], value = TRUE)
  expect_true(all(grepl("\\[(IDENTIFICATION|DESCRIPTIVE|CONFIDENCE)", ivr)))
})

test_that("E28(f): joint spillover S_i widens the set monotonically; Def. 5.1 needs comparison and spillover", {
  f <- fit_mean(three$size, three$vet, three$pre, three$post)
  a0 <- iv(att_mean(f, within_pt(), lce_mean(change = TRUE))); a1 <- iv(att_mean(f, within_pt(), lce_mean(change = TRUE), L = 0.02))
  expect_true(a1[1] <= a0[1] + 1e-12 && a1[2] >= a0[2] - 1e-12)
  expect_true(a1[2] <= a0[2] + 0.02 + 1e-9)                      # at most pi_S * Lambda per fixed fit, here not exceeded
  expect_error(eiATT(f, treated = "v", population = "fixed", trend = within_pt(), spillover = 0, no_anticipation = TRUE),
               "comparison group")
  expect_error(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(), no_anticipation = TRUE),
               "spillover")
  expect_error(att_mean(f, within_pt(), rate_limits = FALSE), "always imposed")
  expect_error(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = parallel_trends(), spillover = 0,
                     no_anticipation = TRUE, lce = lce_mean()), "mean model")
})

test_that("E18(a): RM holes are reported (finding 12), never the envelope", {
  d <- data.frame(u = "U", g1 = 100, g2 = 100, g3 = 100, yes = 150, no = 150)
  cells <- data.frame(unit = "U", row = "g1", col = "yes", election = c("early", "post"), coef = c(4, 1))
  rs <- list(known = list(list(unit = "U", row = "g1", col = "yes", value = 75, election = "base"),
                          list(unit = "U", row = "g2", col = "yes", value = 25, election = "base"),
                          list(unit = "U", row = "g2", col = "yes", value = 60, election = "post"),
                          list(unit = "U", row = "g3", col = "yes", value = 50, election = "early")),
             linear = list(cells = cells, value = 350),
             total = list(units = "U", row = "g1", col = "yes", election = "post", lower = 10, upper = 90))
  # the count set is one polytope: early (.75 - x/2, .25 + x/2, .50), base (.75, .25, .50), post (.50 + 2x, .60, .40 - 2x)
  f <- eiChange(list(early = d, base = d, post = d), cbind(yes, no) ~ cbind(g1, g2, g3), unit = "u", base = "base",
                row = "g1", col = "yes", integer = "lp", restrictions = rs)
  a <- eiATT(f, treated = "g1", comparison = "g2", population = "fixed", trend = rel_mag(1), spillover = 0, no_anticipation = TRUE)
  pc <- a$results[[1]]$pieces
  expect_equal(unname(pc), rbind(c(-0.9, -0.7), c(-0.3, 0)), tolerance = 1e-9)
  expect_false(a$results[[1]]$rm$collapses_to_dial)
  out <- paste(capture.output(print(a)), collapse = "\n")
  expect_match(out, "[-0.900, -0.700] U [-0.300, 0.000]", fixed = TRUE)
  expect_match(out, "holes are part of the answer", fixed = TRUE)
})

test_that("E18(b): RM increments accumulate over the horizon (ATT_2 = -.30 at h = 2)", {
  h <- function(r, s) data.frame(u = c("T", "C"), tr = c(100, 0), co = c(0, 100), yes = c(r, s), no = 100 - c(r, s))
  f <- eiChange(list(e0 = h(35, 50), b = h(50, 50), p1 = h(50, 50), p2 = h(50, 50)), cbind(yes, no) ~ cbind(tr, co),
                unit = "u", base = "b", row = "tr", col = "yes", integer = "lp")
  a <- eiATT(f, treated = "tr", comparison = "co", population = "fixed", trend = rel_mag(1), spillover = 0, no_anticipation = TRUE)
  e2 <- a$results[[2]]
  expect_equal(c(e2$lower, e2$upper), c(-0.3, 0.3), tolerance = 1e-9)        # contains the true ATT_2 = -.30
  expect_equal(e2$rm$endpointwise, c(-0.3, 0.3), tolerance = 1e-9)            # D_2 +/- 2 A
  expect_equal(e2$rm$horizon, 2L)
  expect_equal(c(a$results[[1]]$lower, a$results[[1]]$upper), c(-0.15, 0.15), tolerance = 1e-9)
  # with a pinned path (b_1 = -.15) the second increment can reach -.30: RM on increments, not levels
  expect_true(e2$lower < -0.15 - 1e-9)
})

test_that("E17(c): the sharp breakdown LP is at least the outer value and reaches 0 exactly", {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10))
  yes <- rbind(c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  f <- eiChange(list(base = mk(1), post = mk(2)), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet",
                col = "yes", integer = "lp", restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
  A <- function(M) eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = trend_dial(M), spillover = 0, no_anticipation = TRUE)
  d <- A(0)$results[[1]]$dial
  expect_true(d$sharp_breakdown >= d$breakdown - 1e-9)
  expect_equal(A(d$sharp_breakdown)$results[[1]]$lower, 0, tolerance = 1e-7)
  expect_output(print(A(0)), "sharp breakdown value")
})

test_that("E33(a) with two post-elections: the local dial chains increments; M = 0 is exact WPT at every post", {
  vet <- c(200, 400, 600); N <- 800
  mk <- function(y) data.frame(u = paste0("P", 1:3), v = vet, o = N - vet, yes = y, no = N - y)
  f <- eiChange(list(b = mk(c(4, 164, 324)), p1 = mk(c(204, 124, 4)), p2 = mk(c(250, 200, 100))), cbind(yes, no) ~ cbind(v, o),
                unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  A <- function(M) att_mean(f, within_pt(M = M), lce_mean(base = TRUE))
  a0 <- A(0); a1 <- A(0.01); a2 <- A(0.02)
  expect_length(a0$results, 2)
  # p1 alone (exact WPT does not couple post-periods)
  f1 <- eiChange(list(b = mk(c(4, 164, 324)), p1 = mk(c(204, 124, 4))), cbind(yes, no) ~ cbind(v, o),
                 unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  expect_equal(iv(a0, 1), iv(att_mean(f1, within_pt(), lce_mean(base = TRUE))), tolerance = 1e-9)
  for (k in 1:2) expect_true(iv(a1, k)[1] <= iv(a0, k)[1] + 1e-12 && iv(a2, k)[1] <= iv(a1, k)[1] + 1e-12)
  # the second post's local gap can drift twice as far: widening at p2 is at least that at p1
  expect_true(iv(a0, 2)[1] - iv(a2, 2)[1] >= iv(a0, 1)[1] - iv(a2, 1)[1] - 1e-9)
})
