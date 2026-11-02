args <- commandArgs(TRUE); lib <- if (args[1] == "patched") "/tmp/claude-0/-home-user-eiid/420fe422-ca66-5fa7-b642-d29ac1c8c1c2/scratchpad/ecopatch/lib" else .libPaths()[1]
suppressMessages(library(eco, lib.loc = lib)); data(reg, package = "eco")
cat(args[1], "eco from", find.package("eco"), "\n")
pm <- NULL
for (s in 1:2) { set.seed(10 + s); f <- ecoNP(Y ~ X, N = N, data = reg, n.draws = 22000, burnin = 2000, thin = 3, parameter = TRUE)
  pm <- cbind(pm, colMeans(f$W[, 1, ])); cat(sprintf("%s run %d: alpha %.2f clusters %.1f agg_W1 %.4f\n", args[1], s, median(f$alpha), mean(f$nstar),
     mean(f$W[,1,] %*% (reg$N*reg$X))/sum(reg$N*reg$X))) }
saveRDS(rowMeans(pm), paste0("pm_", args[1], ".rds"))
