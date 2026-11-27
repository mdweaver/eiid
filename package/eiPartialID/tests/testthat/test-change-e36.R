# Prop. E36 (2e Sec. 5): baseline-curvature allowance epsilon on the realised route (E35 + LCE by E13).
# Worked example of docs/extension/verify_phase2e.py (E36): three precincts whose observed returns are
# taken as the expected rates, so the model interval is the plug-in with no allowance (internal
# spec fields .plugin and .t_E reproduce that population example).

e36_fit <- function(eps, plugin = TRUE) {
  a <- cbind(c(50, 100, 150), c(150, 100, 50)); y0 <- c(1, 41, 81); y1 <- c(51, 31, 1)
  mk <- function(y) data.frame(p = paste0("P", 1:3), v = a[, 1], o = a[, 2], yes = y, no = 200 - y)
  m <- lce_levels_change(eps)
  if (plugin) { m$.plugin <- TRUE; m$.t_E <- 0 }
  eiChange(list(b = mk(y0), p = mk(y1)), cbind(yes, no) ~ cbind(v, o), unit = "p", base = "b", row = "v", col = "yes",
           integer = "lp", model = m)
}
e36 <- function(f, ...) eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(...),
                              spillover = 0, no_anticipation = TRUE, lce = TRUE)$results[[1]]
near4 <- function(x, y) expect_true(all(abs(unname(as.vector(x)) - y) <= 1e-4), info = paste(format(x, digits = 6), collapse = ", "))

test_that("E36 worked example: causal sets at epsilon = 0, 0.001, 0.01, 1/2 and the breakdown epsilon*", {
  r0 <- e36(e36_fit(0)); r1 <- e36(e36_fit(0.001)); r2 <- e36(e36_fit(0.01)); r5 <- e36(e36_fit(0.5))
  near4(c(r0$lower, r0$upper), c(0.0267, 0.0800))
  near4(c(r1$lower, r1$upper), c(0.0107, 0.0960))
  near4(c(r2$lower, r2$upper), c(-0.1333, 0.2400))
  near4(r5$lce$pieces, c(-0.4933, 0.6867))
  near4(c(r5$lower, r5$upper), c(-0.4933, 0.2733))
  near4(r0$breakdown_eps, 1 / 600)
  expect_equal(r1$breakdown_eps, r0$breakdown_eps, tolerance = 1e-6)    # eps* does not depend on the fitted epsilon
  expect_true(is.na(r2$breakdown_eps))                                   # the set contains 0
  # nesting in epsilon: CI_LCE and the causal set
  es <- c(0, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5)
  iv <- t(vapply(es, function(e) { r <- e36(e36_fit(e)); c(r$lce$pieces, r$lower, r$upper) }, numeric(4)))
  expect_true(all(diff(iv[, 1]) <= 1e-9) && all(diff(iv[, 2]) >= -1e-9) && all(diff(iv[, 3]) <= 1e-9) && all(diff(iv[, 4]) >= -1e-9))
  out <- capture.output(print(eiATT(e36_fit(0), treated = "v", comparison = "o", population = "fixed", trend = within_pt(),
                                    spillover = 0, no_anticipation = TRUE, lce = TRUE)))
  expect_true(any(grepl("epsilon* = 0.001667", out, fixed = TRUE)))
  expect_true(any(grepl("a sensitivity summary, not a test", out, fixed = TRUE)))
  expect_true(any(grepl("[CONFIDENCE SET, realised ATT (E35 + E13)]", out, fixed = TRUE)))
  o5 <- capture.output(print(e36_fit(0.5)))
  expect_true(any(grepl("change LCE plus the base-level link", o5, fixed = TRUE)))
})

test_that("E36 epsilon = 0 reproduces lce_periods(); eiATT pretrends = 'model' refuses it", {
  set.seed(3); p <- 40
  x <- runif(p, .05, .95); N <- sample(200:400, p, TRUE); a1 <- round(N * x)
  t0 <- .3 + .2 * x + rnorm(p, 0, .03); t1 <- t0 + .02 + .1 * x + rnorm(p, 0, .02)
  mkT <- function(l) round(N * pmin(pmax(l, 0), 1))
  mk <- function(t) data.frame(u = seq_len(p), v = a1, o = N - a1, yes = mkT(t), no = N - mkT(t))
  d <- list(b = mk(t0), p = mk(t1))
  G <- function(m) eiChange(d, cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp", model = m)
  cols <- c("plugin_lower", "plugin_upper", "model_lower", "model_upper", "allowance")
  expect_equal(G(lce_levels_change(0))$model$table[, cols], G(lce_periods())$model$table[, cols], tolerance = 1e-8)
  f <- G(lce_levels_change(0.01))
  expect_error(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = rel_mag(1), spillover = 0,
                     no_anticipation = TRUE, pretrends = "model"), "within_pt")
})
