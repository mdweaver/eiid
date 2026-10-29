# Benchmark tables from ei.Datasets (real elections with known truth). The full benchmark is
# benchmarks/run_ei_datasets_benchmark.R in the repository; this test checks the table builder
# and runs the machinery on a few electorates.
test_that("ei.Datasets tables are built with exact truth margins", {
  skip_if_not_installed("ei.Datasets")
  tb <- eiDatasetsTables(elections = c(2002, 2020))
  expect_equal(length(tb), 69 + 72)
  for (e in tb) {
    expect_equal(unname(rowSums(e$truth)), unname(colSums(e$data[, e$rows])))
    expect_equal(unname(colSums(e$truth)), unname(colSums(e$data[, e$cols])))
  }
})

test_that("benchmark on a few electorates: truth within Duncan-Davis, CIs within Duncan-Davis", {
  skip_if_not_installed("ei.Datasets")
  skip_on_cran()
  tb <- eiDatasetsTables(elections = 2002)[1:3]
  b <- suppressWarnings(eiBenchmark(tb, rungs = c("M0", "L6"), seed = 1, ndraw = 30, nz_draw = 1500))
  expect_true(all(b$in_dd))
  ok <- !is.na(b$ci_actual_lower)
  expect_true(all(b$ci_actual_lower[ok] >= b$dd_lower[ok] - 1e-9 & b$ci_actual_upper[ok] <= b$dd_upper[ok] + 1e-9))
  expect_true(all(c("election", "district", "units") %in% names(b)))
  expect_true(mean(b$in_ci_actual[b$rung == "M0"]) >= 0.8)
})
