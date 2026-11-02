suppressMessages({library(eco); library(eiDPM); library(parallel)}); data(reg, package="eco")
d <- dpmData(reg$X, reg$Y)
g <- function(seed, init_k, init_alpha, K=30L, iters=20000L) { set.seed(seed)
  r <- eiDPM:::dpm_gibbs_cpp(2L,K,d$X,d$Y,d$L,d$U,d$Wobs,d$w1obs,d$w2obs,c(0,0),2,4,diag(10,2),1,0.1,iters,2000L,4L,(d$L+d$U)/2,1L,init_k,init_alpha)
  W1 <- matrix(NA, nrow(r$W1), nrow(reg)); W1[, d$lat] <- r$W1; W1[, d$yb] <- d$Wyb
  list(pm = colMeans(W1), alpha = median(r$alpha), occ = median(r$n_clusters), agg = mean(W1 %*% (reg$N*reg$X))/sum(reg$N*reg$X)) }
e <- function(seed, draws) { set.seed(seed); f <- ecoNP(Y~X, N=N, data=reg, n.draws=draws, burnin=2000, thin=3, parameter=TRUE)
  list(pm = colMeans(f$W[,1,]), alpha = median(f$alpha), occ = median(f$nstar), agg = mean(f$W[,1,] %*% (reg$N*reg$X))/sum(reg$N*reg$X)) }
jobs <- list(list("gibbs one-cluster", function() g(1,1,1)), list("gibbs one-cluster s2", function() g(2,1,1)),
             list("gibbs 30-cluster", function() g(3,30,1)), list("gibbs 30-cluster s4", function() g(4,30,1)),
             list("gibbs K=100", function() g(5,1,1,K=100L)),
             list("eco 22000", function() e(11,22000)), list("eco 22000 s2", function() e(12,22000)))
res <- mclapply(jobs, function(j) j[[2]](), mc.cores=4, mc.preschedule=FALSE)
names(res) <- sapply(jobs, `[[`, 1)
ref <- (res[["eco 22000"]]$pm + res[["eco 22000 s2"]]$pm)/2
for (nm in names(res)) cat(sprintf("%-22s agg_W1 %.4f  alpha %.2f  clusters %4.1f  cor(eco) %.3f  rmsd(eco) %.4f\n", nm, res[[nm]]$agg, res[[nm]]$alpha, res[[nm]]$occ, cor(res[[nm]]$pm, ref), sqrt(mean((res[[nm]]$pm-ref)^2))))
