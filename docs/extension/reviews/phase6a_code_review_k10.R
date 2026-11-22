# K10: every model output of eiBounds() under target = "model" is identical to 0.3.1.
# Usage: Rscript phase6a_code_review_k10.R LIB OUTDIR TESTDIR   (run once with the 0.3.1 library, once with 0.4.0)
a <- commandArgs(TRUE)
if (a[1] == "compare") {
  A <- readRDS(file.path(a[2], "k10_0.3.1.rds")); B <- readRDS(file.path(a[2], "k10_0.4.0.rds"))
  for (nm in names(A)) {
    fa <- A[[nm]]$fit; fb <- B[[nm]]$fit
    tb <- fb$table[, setdiff(names(fb$table), c("dd_exact_lower", "dd_exact_upper"))]
    comps <- setdiff(intersect(names(fa), names(fb)), "table")
    diffc <- comps[!vapply(comps, function(k) identical(fa[[k]], fb[[k]]), TRUE)]
    cat(sprintf("%-7s table identical %s; other components differing: %s; new: %s; print identical %s\n", nm, identical(fa$table, tb),
                if (length(diffc)) paste(diffc, collapse = ",") else "none", paste(setdiff(names(fb), names(fa)), collapse = ","),
                identical(A[[nm]]$print, B[[nm]]$print)))
  }
  quit(save = "no")
}
lib <- a[1]; out <- a[2]; td <- a[3]
library(eiPartialID, lib.loc = if (lib == "default") NULL else lib)
v <- as.character(packageVersion("eiPartialID", lib.loc = if (lib == "default") NULL else lib))
d3 <- read.csv(file.path(td, "rxc_example_3x3.csv")); d2 <- read.csv(file.path(td, "msign_example_3x2.csv"))
X3 <- as.matrix(d3[, c("X1", "X2", "X3")]); T3 <- as.matrix(d3[, c("T1", "T2", "T3")])
X2 <- as.matrix(d2[, c("X1", "X2", "X3")]); T2 <- as.matrix(d2[, c("T1", "T2")])
tg <- list(c11 = c(1, 1), gap = rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 0, 0)), rowsum = rbind(c(1, 1, 0), c(0, 0, 0), c(0, 0, 0)))
tg2 <- list(c11 = c(1, 1), gap = rbind(c(1, 0), c(0, 0), c(-1, 0)))
cfg <- list(
  M0 = function(...) eiBounds(x = X3, t = T3, n = d3$N, model = "M0", targets = tg, seed = 1, ...),
  M0cons = function(...) eiBounds(x = X3, t = T3, n = d3$N, model = "M0", ci = "conservative", targets = tg, seed = 1, ...),
  M1 = function(...) eiBounds(x = X3, t = T3, n = d3$N, group = d3$group, model = "M1", targets = tg, seed = 2, ...),
  M2 = function(...) eiBounds(x = X3, t = T3, n = d3$N, group = d3$group, model = "M2", seed = 3, ...),
  L6 = function(...) eiBounds(x = X3, t = T3, n = d3$N, model = "L6", targets = tg, ...),
  L5 = function(...) eiBounds(x = X3, t = T3, n = d3$N, group = d3$group, model = "L5", ...),
  msign = function(...) eiBounds(x = X2, t = T2, n = d2$N, model = "M0", targets = tg2, seed = 4,
                                 restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1))), ...),
  cr = function(...) eiBounds(x = X2, t = T2, n = d2$N, model = "M0", targets = tg2, seed = 5,
                              restrictions = list(cr = list(list(row = 1, col = 1))), ...),
  popov = function(...) eiBounds(x = X2, t = T2, n = d2$N, model = "M0", fpc = "popoviciu", targets = tg2, seed = 6, ...),
  noseed = function(...) { set.seed(77); f <- eiBounds(x = X2, t = T2, n = d2$N, model = "M0", targets = tg2, ...); f$after <- runif(1); f })
res <- list()
for (nm in names(cfg)) {
  f <- if (v == "0.3.1") cfg[[nm]]() else cfg[[nm]](target = "model")
  pr <- capture.output(print(f))
  f$call <- NULL
  res[[nm]] <- list(fit = f, print = pr, time = NA)
}
saveRDS(res, file.path(out, paste0("k10_", v, ".rds")))
cat("saved", v, "\n")
# compare (run after both): Rscript phase6a_code_review_k10.R compare OUTDIR
