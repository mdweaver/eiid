suppressMessages({library(eco); library(eiDPM); library(parallel)}); data(reg, package="eco")
cl <- function(w) pmin(pmax(w, 1e-4), 1 - 1e-4)
Wobs <- rbind(qlogis(cl(reg$W1)), qlogis(cl(reg$W2)))
# Gibbs: no latent units; all 275 units fully observed
g <- function(seed) { set.seed(seed)
  r <- eiDPM:::dpm_gibbs_cpp(2L,50L,numeric(0),numeric(0),numeric(0),numeric(0),Wobs,numeric(0),numeric(0),c(0,0),2,4,diag(10,2),1,0.1,22000L,2000L,4L,numeric(0),1L,1L,1)
  c(alpha=median(r$alpha), clusters=mean(r$n_clusters)) }
# eco: two latent units (needed by the interface) plus the 275 observed units as survey data
e <- function(seed) { set.seed(seed)
  dd <- data.frame(X=c(0.3,0.6), Y=c(0.5,0.7), N=c(1,1))
  f <- ecoNP(Y~X, data=dd, supplement=cbind(cl(reg$W1), cl(reg$W2)), n.draws=22000, burnin=2000, thin=3, parameter=TRUE)
  c(alpha=median(f$alpha), clusters=mean(f$nstar)) }
res <- mclapply(list(function() g(1), function() g(2), function() e(11), function() e(12)), function(f) f(), mc.cores=4)
print(rbind(gibbs1=res[[1]], gibbs2=res[[2]], eco1=res[[3]], eco2=res[[4]]))
