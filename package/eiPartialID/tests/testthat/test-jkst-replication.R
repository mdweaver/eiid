# The R x C engine must reproduce the original 2x2 package (Jiang et al. 2019) exactly.
test_that("2x2 plug-in bound and CI_x reproduce generateBounds / evaluateBounds", {
  skip_if_not_installed("eco")
  census <- NULL
  utils::data("census", package = "eco", envir = environment())
  x <- census$X; t <- census$Y; n <- census$N
  old <- generateBounds(x, t, n, printSummary = FALSE)
  ev <- NULL
  utils::capture.output(ev <- evaluateBounds(old))
  fit <- eiBounds(x = x, t = t, n = n, model = "M0", range_offset = 1e-5, vcov = "HC1",
                  targets = list(B = c(1, 1)), jkst = TRUE, ci = "clr", seed = 1, tests = FALSE)
  expect_equal(fit$table$estimate_lower, unname(old$hbdl0), tolerance = 1e-10)
  expect_equal(fit$table$estimate_upper, unname(old$hbdu0), tolerance = 1e-10)
  expect_equal(fit$table$dd_lower, old$bdl, tolerance = 1e-12)
  expect_equal(fit$table$dd_upper, old$bdu, tolerance = 1e-12)
  jk <- fit$jkst$B
  expect_equal(jk$lower, ev$CI_x_lower, tolerance = 1e-10)
  expect_equal(jk$upper, ev$CI_x_upper, tolerance = 1e-10)
  expect_equal(jk$selected, ev$CI_x_isSelected)
})
