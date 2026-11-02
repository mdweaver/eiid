suppressMessages(library(eiDPM)); data(reg, package = "eco"); pm <- NULL
for (s in 1:2) { g <- eiDPM(reg$X, reg$Y, reg$N, K = 30, n.draws = 20000, burnin = 2000, thin = 4, seed = 40 + s); pm <- cbind(pm, colMeans(g$W1)) }
saveRDS(rowMeans(pm), "pm_gibbs.rds")
