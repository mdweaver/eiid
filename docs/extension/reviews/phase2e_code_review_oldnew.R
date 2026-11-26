# Existing outputs unchanged (0.4.0 vs 0.5.0): Rscript phase2e_code_review_oldnew.R LIB OUT.rds, once per library,
# then compare the captured printouts (identical() per object). Code review of PHASE2e.
args <- commandArgs(TRUE); .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(eiPartialID))
set.seed(7); p <- 80
X <- matrix(runif(p * 3), p); X <- X / rowSums(X); N <- sample(100:300, p, TRUE)
A <- round(X * N); A[, 3] <- N - A[, 1] - A[, 2]
T <- cbind(rbinom(p, N, .4), 0); T[, 2] <- N - T[, 1]
d <- data.frame(a1 = A[, 1], a2 = A[, 2], a3 = A[, 3], y = T[, 1], n = T[, 2])
out <- list()
set.seed(1); out$b0 <- eiBounds(cbind(y, n) ~ cbind(a1, a2, a3), data = d, model = "M0", target = "model", ci = "clr", tests = FALSE)
set.seed(1); out$bc <- eiBounds(cbind(y, n) ~ cbind(a1, a2, a3), data = d, model = "M0", target = "model", ci = "conservative", tests = FALSE,
                   restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1))))
out$c1 <- eiCounts(cbind(y, n) ~ cbind(a1, a2, a3), data = d[1:12, ], targets = list(rate(1, 1), gap(1, 2, 1)),
                   restrictions = list(cr = list(row = 1, col = 1), known = list(row = 2, col = 1, unit = 3, value = 5)))
out$c2 <- eiCounts(cbind(y, n) ~ cbind(a1, a2, a3), data = d, targets = list(rate(1, 1), gap(1, 2, 1)), model = out$b0)
strip <- function(o) { o$call <- NULL; o$P <- NULL; o$fit <- NULL; o }
saveRDS(lapply(out, function(o) capture.output(print(o))), args[2])
