# Prop. E35 (2e Sec. 5): realised within-precinct parallel trends with tolerance tau_i, the default
# reading of within_pt(). Worked example on the Sec. 9.2 data (docs/extension/verify_phase2e.py, E35 checks).

e35_data <- function() {
  aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10))
  yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
  mk <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
  list(early = mk(1), base = mk(2), post = mk(3))
}
e35_fit <- function(rows = NULL) eiChange(e35_data(), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base",
                                          row = "vet", col = "yes", integer = "lp", restrictions = rows)
kc <- function(u, e, v) list(unit = u, row = "oth", col = "yes", value = v, election = e)
e35_aud <- list(known = list(kc("P3", "base", 20), kc("P3", "post", 22), kc("P4", "base", 4), kc("P4", "post", 5)))   # hypothetical audited cells
e35 <- function(f, ...) eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = within_pt(...),
                              spillover = 0, no_anticipation = TRUE)
near4 <- function(x, y) expect_true(all(abs(unname(x) - y) <= 1e-4), info = paste(format(x, digits = 6), collapse = ", "))
ivl <- function(a) c(a$results[[1]]$lower, a$results[[1]]$upper)

test_that("E35 worked example: margins only, audited cells at M = 0 and 0.1, breakdown M*", {
  near4(ivl(e35(e35_fit())), c(-0.5675, 0.7322))
  fa <- e35_fit(e35_aud)
  a0 <- e35(fa, M = 0)
  near4(ivl(a0), c(0.0525, 0.3522))
  near4(ivl(e35(fa, M = 0.1)), c(0.0310, 0.3740))
  near4(a0$results[[1]]$breakdown_M, 0.2438)
  # M = 0 lies inside D_W's count set [-0.0087, 0.3772] (the local rate limits bind)
  near4(a0$results[[1]]$dw, c(-0.0087, 0.3772))
  # at M*, zero is the lower end; below it the sign holds
  Ms <- a0$results[[1]]$breakdown_M
  expect_lt(abs(ivl(e35(fa, M = Ms))[1]), 1e-7)
  expect_gt(ivl(e35(fa, M = 0.9 * Ms))[1], 0)
  # a user tau vector equal to the default tolerance reproduces M = 0.1
  aV <- c(10, 50, 150, 190); aN <- c(190, 150, 50, 10)
  expect_equal(ivl(e35(fa, tau = 0.1 * sqrt(1 / aV + 1 / aN))), ivl(e35(fa, M = 0.1)), tolerance = 1e-9)
})

test_that("E35 nesting in M; default reading is realised; expected reading still available", {
  fa <- e35_fit(e35_aud)
  iv <- t(vapply(c(0, 0.05, 0.1, 0.2, 0.5), function(M) ivl(e35(fa, M = M)), numeric(2)))
  expect_true(all(diff(iv[, 1]) <= 1e-12) && all(diff(iv[, 2]) >= -1e-12))
  expect_identical(within_pt()$target, "realised")
  ae <- eiATT(fa, treated = "vet", comparison = "oth", population = "fixed", trend = within_pt(target = "expected"),
              spillover = 0, no_anticipation = TRUE)
  expect_false(is.null(ae$wpt)); expect_null(ae$rwpt)
  expect_error(e35(fa, target = "expected", tau = 1), "realised")
})

test_that("E35 unsupported precincts: theta_i in [0, 1], pi_U reported", {
  a <- e35(e35_fit(), a_min = 60)                    # P3 (50) and P4 (10) unsupported
  expect_equal(a$rwpt$design$pi_U, 340 / 400)
  out <- capture.output(print(a))
  expect_true(any(grepl("pi_U = 0.85", out, fixed = TRUE)))
  expect_true(any(grepl("P4 10 (unsupported)", out, fixed = TRUE)))
  # bounding the unsupported precincts with no assumption can only widen the set
  i0 <- ivl(e35(e35_fit())); i1 <- ivl(a)
  expect_true(i1[1] <= i0[1] + 1e-12 && i1[2] >= i0[2] - 1e-12)
})

test_that("E35 labels, R4-1 note, LCE only through E13 (confidence set with the coverage line)", {
  out <- capture.output(print(e35(e35_fit(e35_aud))))
  expect_true(any(grepl("[IDENTIFICATION SET, realised ATT, count layer (E35)]", out, fixed = TRUE)))
  expect_false(any(grepl("[CONFIDENCE SET", out, fixed = TRUE)))
  expect_true(any(grepl("not a confidence level", out, fixed = TRUE)))
  expect_true(any(grepl("M* = 0.2438", out, fixed = TRUE)))
  set.seed(3); p <- 40
  x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
  t0 <- .3 + .2 * x + rnorm(p, 0, .03); t1 <- t0 + .02 + .1 * x + rnorm(p, 0, .02); t2 <- t1 + .05 + .15 * x + rnorm(p, 0, .02)
  mkT <- function(l) round(N * pmin(pmax(l, 0), 1))
  mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = mkT(t), no = N - mkT(t))
  f <- eiChange(list(m1 = mk(t0), b = mk(t1), p = mk(t2)), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b",
                row = "v", col = "yes", integer = "lp", model = fd())
  A <- function(...) eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(M = 1),
                           spillover = 0, no_anticipation = TRUE, ...)
  a <- A(); b <- A(lce = TRUE)
  ia <- ivl(a); ib <- ivl(b)
  expect_true(ib[1] >= ia[1] - 1e-9 && ib[2] <= ia[2] + 1e-9)      # the E13 row only intersects
  expect_equal(unname(b$results[[1]]$lce$rho)[b$rwpt$design$aV > 0], 4 * b$rwpt$design$wV[b$rwpt$design$aV > 0], tolerance = 1e-12)
  ob <- capture.output(print(b))
  expect_true(any(grepl("[CONFIDENCE SET, realised ATT (E35 + E13)]", ob, fixed = TRUE)))
  expect_true(any(grepl("under FP, NA, RWPT-tau, |S_i| <= Lambda", ob, fixed = TRUE)))
  expect_false(any(grepl("[IDENTIFICATION SET, realised ATT, count layer (E35)]", ob, fixed = TRUE)))
  expect_error(eiATT(e35_fit(), treated = "vet", comparison = "oth", population = "fixed", trend = within_pt(),
                     spillover = 0, no_anticipation = TRUE, lce = TRUE), "model layer")
  expect_error(A(lce = lce_mean()), "target = \"expected\"")
})
