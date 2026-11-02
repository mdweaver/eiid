# Calibration plan for the `eiDPM` exact-DP sampler (pre-registered)

*This file is committed before any run of `benchmarks/eco_2x2/calibration.R`, and the margins below are fixed in advance. Results go to `benchmarks/eco_2x2/results/calibration_*` and are reported against these criteria whatever they show.*

## Purpose

Check that the collapsed exact-DP Gibbs sampler (`eiDPM(..., K = Inf)`; Neal 2000 Alg. 3 with the Escobar–West $\alpha$ update) draws from the posterior of the ILS model it states. Earlier SBC found two problems: the truncated blocked sampler under-represents the cluster count when $\alpha$ is large, and `ecoNP` is biased in $\alpha$ and the cluster count. This plan uses equivalence criteria rather than failure to reject uniformity.

## Designs

All designs use $n=100$ units, $X\sim U(0.05,0.95)$ for CAR, and `eco`'s default prior ($\mu_0=0$, $\tau_0=2$, $\nu_0=4$, $S_0=10I$, $\alpha\sim\Gamma(1,0.1)$; $\mu\mid\Sigma\sim N(\mu_0,\Sigma/\tau_0)$).

| Design | Data generated from | Fitted with | Replications |
|---|---|---|---|
| **D1 (primary)** | exact DP mixture, CAR (Chinese restaurant process with $\alpha$; a new table draws $(\mu,\Sigma)$ from the base measure) | exact DP (collapsed), $K=\infty$ | 4000 |
| D1b (contrast) | the D1 simulations | truncated blocked Gibbs, $K=50$ (the old sampler) | 2000 |
| D2 | parametric CAR (one component) | parametric Gibbs, $K=1$ | 2000 |
| D3 | exact DP mixture, NCAR ($p=3$: $(\operatorname{logit}W,\operatorname{logit}X)$) | exact DP (collapsed), NCAR | 2000 |

**Data-only rejection.** A replication is discarded and redrawn if any $Y_i$ (or, in D3, any $X_i$) is within $10^{-6}$ of 0 or 1. This conditions on the observed data only, so the posterior is unchanged and SBC remains valid.

**MCMC settings.** 1000 burn-in sweeps, then $L=199$ draws at thinning 25. With $L+1=200$ possible ranks, the central 50/80/90/95% intervals of the draws contain exactly the nominal share of ranks.

## Test quantities

* Unit rates $W_{1,1},W_{1,2},W_{1,3}$ (three fixed units).
* In-sample means $\bar W_1,\bar W_2$.
* In-sample $\operatorname{sd}(W_1)$ and $\operatorname{cor}(W_1,W_2)$ across units. These are joint quantities, sensitive to the dispersion and clustering of the rates (after Modrák et al. 2023).
* $\alpha$ and the number of occupied clusters (DP designs).

The data-given-parameter log-likelihood suggested by Modrák et al. is not used. Here $Y$ is a deterministic function of $W$, so it is degenerate; the joint quantities above stand in for it.

**Rank.** The number of the $L$ draws below the truth, with ties (the cluster count) broken uniformly at random. Write $u=(\text{rank}+0.5)/(L+1)$.

## Metrics and pre-registered margins

All confidence intervals are 90%, i.e. two one-sided tests at 5% each. In the table, $q$ is the nominal coverage level and $\mathrm{Var}_U$ is the variance of the discrete uniform on the rank scale.

| Metric | Estimate | CI | Equivalence region |
|---|---|---|---|
| Central-interval coverage at 50, 80, 90, 95% | share of ranks in the central $q$ block | Clopper–Pearson | $q\pm0.03$ |
| Left and right tail misses at 90% | share of ranks $<10$; share $\ge190$ | Clopper–Pearson | $[0.03,0.07]$ each |
| Mean rank | mean of $u$ | normal | $[0.47,0.53]$ |
| Rank-variance ratio | $\operatorname{var}(u)/\mathrm{Var}_U$ | bootstrap percentile (2000 resamples) | $[0.90,1.10]$ |
| Max ECDF deviation | $D=\sup\lvert F_u-U\rvert$ | bootstrap upper 95% bound | $<0.05$ |

**Classification per quantity × design.**
* *Equivalent:* every CI lies inside its region.
* *Miscalibrated:* some CI lies entirely outside its region.
* *Inconclusive:* otherwise.

**Primary criterion.** The sampler is called **calibrated** on a design if every test quantity is equivalent on 90% coverage, and none is miscalibrated on any metric.

**Specific claim to test: "the cluster-count weakness is fixed."** In D1, $\alpha$ and the cluster count are not miscalibrated on any metric, and both are equivalent on 90% coverage. D1b is expected to fail; it is there to show that the test has the power to detect the old failure.

## Diagnostics reported alongside (no pass/fail)

* **ESS of each test quantity per fit** (`posterior::ess_bulk` over the 199 draws). Report the median and the share of fits below 100. Rank tests assume near-independent draws, and autocorrelation alone produces U-shaped ranks.
* **ECDF-difference plots with simultaneous 95% bands.** The band is the pointwise binomial band at an adjusted level $\gamma$, calibrated by simulation so that a uniform sample stays inside it everywhere with probability 0.95 (Säilynoja, Bürkner & Vehtari 2022).
* **Posterior z-scores** $(\text{posterior mean}-\text{truth})/\text{posterior sd}$: mean and sd (ideal 0 and 1).
* **Contraction** $1-\operatorname{var}_{\text{post}}/\operatorname{var}_{\text{prior}}$, with the prior variance estimated from the simulated truths.
* **For unit rates:** the median ratio of the 90% posterior interval width to the unit's Duncan–Davis width. This is the share of the identified interval the model claims to rule out.
