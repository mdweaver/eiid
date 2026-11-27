# eiChange / eiATT: reproduce PHASE2e Sec. 9.2 (worked example), N21, and the counterexamples behind
# proof-review issues 1, 2 and 12 (docs/extension/verify_phase2e_output.txt, check W, E17 [F],
# E6 [F], E18). Values printed in 2e Sec. 9.2 to 4 decimals (half away from zero), hence the
# tolerance 5e-5 (plus floating-point slack).

w_data <- function() {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10))
  yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  list(early = mk(1), base = mk(2), post = mk(3))
}
w_fit <- function(rows = NULL, integer = "lp")
  eiChange(w_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet", col = "yes",
           integer = integer, restrictions = rows)
crw <- function(a, b) list(row = "vet", col = "yes", from = a, to = b)
att <- function(f, tr, L = 0, ...) {
  a <- eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = tr, spillover = L,
             no_anticipation = TRUE, ...)
  c(a$results[[1]]$lower, a$results[[1]]$upper)
}
near <- function(x, y) expect_true(all(abs(unname(x) - y) <= 5e-5 + 1e-9), info = paste(format(x, digits = 6), collapse = ", "))

test_that("2e Sec. 9.2 worked example: every row of the ladder", {
  f0 <- w_fit(); fP <- w_fit(list(cr_change = crw("base", "post")))
  fB <- w_fit(list(cr_change = list(crw("base", "post"), crw("early", "base"))))
  # descriptive
  tb <- f0$table
  near(c(tb$dd_lower[2], tb$dd_upper[2]), c(-0.4075, 0.7925))          # D_post, margins
  near(c(tb$dd_lower[1], tb$dd_upper[1]), c(-0.6025, 0.5975))          # D_{-1} = -(D_pre): D_pre in [-0.5975, 0.6025]
  near(c(fP$table$lower[2], fP$table$upper[2]), c(0.1433, 0.2703))     # CR on post changes
  near(c(f0$closed[[2]]$ovr$ER, f0$closed[[2]]$ovr$NM), c(0.2703, 0.1433))
  near(c(f0$closed[[1]]$ovr$ER, f0$closed[[1]]$ovr$NM), -c(0.0080, 0.0043))
  # causal
  near(att(f0, no_trend()), c(-0.5675, 0.7325))
  near(att(f0, parallel_trends()), c(-0.4075, 0.7325))
  near(att(fP, parallel_trends()), c(0.1433, 0.2703))
  near(att(fP, trend_dial(0.05)), c(0.0933, 0.3203))
  near(att(fP, trend_dial(0.05), 0.02), c(0.0733, 0.3403))
  near(att(f0, rel_mag(1)), c(-0.5675, 0.7325))
  near(att(fP, rel_mag(1)), c(-0.4050, 0.7325))
  near(att(fB, rel_mag(1)), c(0.1352, 0.2783))
  near(att(f0, mtr_nd_ps()), c(0, 0.4825))
  near(att(fP, mtr_nd_ps()), c(0, 0.3114))
  # the within-precinct comparison, descriptive: D_W (Prop. E29)
  a <- eiATT(f0, treated = "vet", comparison = "oth", population = "fixed", trend = within_pt(), spillover = 0, no_anticipation = TRUE)
  near(a$results[[1]]$dw$LU, c(-0.9387, 1.4472))
})

test_that("review issue 12: RM with CR on the pre-war change only vs on both changes", {
  fQ <- w_fit(list(cr_change = crw("early", "base")))
  fB <- w_fit(list(cr_change = list(crw("base", "post"), crw("early", "base"))))
  near(att(fQ, rel_mag(1)), c(-0.4155, 0.7325))
  near(att(fB, rel_mag(1)), c(0.1352, 0.2783))
  near(c(fQ$table$lower[1], fQ$table$upper[1]), -c(0.0080, 0.0043))     # pre DiD pinned to [0.0043, 0.0080]
})

n21 <- function(rows = NULL) {
  a21 <- cbind(c(20, 60), c(80, 40)); Y <- rbind(c(30, 45), c(32, 48), c(30, 60))
  mk <- function(t) data.frame(u = c("A", "B"), v = a21[, 1], o = a21[, 2], yes = Y[t, ], no = 100 - Y[t, ])
  eiChange(list(m1 = mk(1), b = mk(2), p = mk(3)), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b",
           row = "v", col = "yes", integer = "lp", restrictions = rows)
}
att21 <- function(f, ...) eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(1),
                                spillover = 0, no_anticipation = TRUE, ...)

test_that("N21 reproduced exactly and corrected (2e Sec. 5.5)", {
  f <- n21()
  expect_equal(c(f$table$dd_lower[1], f$table$dd_upper[1]), -c(1.270833, -1.229167), tolerance = 1e-6)
  expect_equal(c(f$table$dd_lower[2], f$table$dd_upper[2]), c(-1.083333, 1.416667), tolerance = 1e-6)
  a <- att21(f, rate_limits = FALSE)
  expect_equal(c(a$results[[1]]$lower, a$results[[1]]$upper), c(-2.354167, 2.645833), tolerance = 1e-6)  # joint
  expect_equal(a$results[[1]]$rm$endpointwise, c(-2.354167, 2.6875), tolerance = 1e-6)                    # endpoint-wise
  a2 <- att21(f)
  expect_equal(c(a2$results[[1]]$lower, a2$results[[1]]$upper), c(-0.75, 1.00), tolerance = 1e-9)
  expect_true(a2$results[[1]]$rm$no_information)
  expect_equal(c(a2$results[[1]]$none$lower, a2$results[[1]]$none$upper), c(-0.75, 1.00), tolerance = 1e-9)
  expect_output(print(a2), "NO INFORMATION")
  expect_output(print(a2), "implied M")
  expect_output(print(a, digits = 3), "OUTER")
})

test_that("N21 narrow variant: pre-period totals pinned to +/-2 around the Duncan-Davis midpoint", {
  a21 <- cbind(c(20, 60), c(80, 40)); Y <- rbind(c(30, 45), c(32, 48))
  pins <- lapply(1:2, function(t) {
    lo <- sum(pmax(0, a21[, 1] + Y[t, ] - 100)); hi <- sum(pmin(a21[, 1], Y[t, ])); ctr <- lo + (hi - lo) / 2
    list(row = "v", col = "yes", election = c("m1", "b")[t], lower = ctr - 2, upper = ctr + 2, units = c("A", "B"))
  })
  f <- n21(list(total = pins))
  expect_equal(c(f$table$lower[1], f$table$upper[1]), c(-5 / 48, 1 / 16), tolerance = 1e-9)   # pre DiD in [-1/16, 5/48]
  a <- att21(f, rate_limits = FALSE); a2 <- att21(f)
  # 2e N21: [-29/48, 43/48] both without and with the rate limits (they do not bind here)
  expect_equal(c(a$results[[1]]$lower, a$results[[1]]$upper), c(-29 / 48, 43 / 48), tolerance = 1e-9)
  expect_equal(c(a2$results[[1]]$lower, a2$results[[1]]$upper), c(-29 / 48, 43 / 48), tolerance = 1e-9)
})

test_that("review issue 1: E17(b) exact criterion with a coupling row", {
  d <- data.frame(u = "U", tr = 100, co = 100, yes = 100, no = 100)
  cells <- data.frame(unit = "U", row = "tr", col = "yes", election = c("e0", "e1"), coef = c(3, -2))
  f <- eiChange(list(e0 = d, e1 = d), cbind(yes, no) ~ cbind(tr, co), unit = "u", row = "tr", col = "yes", integer = "lp",
                restrictions = list(linear = list(cells = cells, value = 140)))
  a <- eiATT(f, treated = "tr", comparison = "co", population = "fixed", trend = parallel_trends(), spillover = 0,
             no_anticipation = TRUE)
  e <- a$results[[1]]
  expect_equal(e$D[2], -0.4, tolerance = 1e-9)
  expect_equal(c(e$lower, e$upper), c(-14 / 15, -0.8), tolerance = 1e-9)
  expect_false(e$dial$criterion_upper)               # the criterion predicts upper < max D
  expect_true(e$dial$criterion_lower)
  expect_true(e$upper < e$D[2] - 1e-6)
})

test_that("review issue 2: E6 emptiness when Dbar+ < 0 (mix belief refuted by the margins)", {
  aE <- cbind(c(20, 0, 80, 0), c(0, 80, 0, 20)); yE <- rbind(c(10, 40, 40, 4), c(8, 40, 72, 16))
  mk <- function(t) data.frame(u = paste0("U", 1:4), g1 = aE[, 1], g2 = aE[, 2], yes = yE[t, ], no = rowSums(aE) - yE[t, ])
  f <- eiChange(list(a = mk(1), b = mk(2)), cbind(yes, no) ~ cbind(g1, g2), unit = "u",
                targets = list(bin_did("yes", "a", "b", low = c("U1", "U2"), high = c("U3", "U4"))),
                restrictions = list(contrast_change = list(col = "yes", from = "a", to = "b", low = c("U1", "U2"), high = c("U3", "U4"))))
  cc <- f$closed[[1]]$contrast
  expect_equal(f$table$dd_upper, -0.15, tolerance = 1e-12)
  expect_equal(cc$lambda, 0.6, tolerance = 1e-12); expect_equal(cc$dtau, 0.46, tolerance = 1e-12)
  expect_true(cc$empty)
  expect_true(cc$dtau > 0 && cc$cap >= f$table$dd_lower)   # the earlier criterion would not flag it
  expect_true(any(grepl("Dbar\\+ < 0", cc$why)))
  expect_false(f$table$continuous_feasible)
})
