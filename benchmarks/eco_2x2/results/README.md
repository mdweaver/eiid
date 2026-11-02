# Results

* `summary.csv`, `units.csv`, `report.md`: full run of `run_benchmark.R 4` (eco vs the eiDPM Gibbs sampler; Stan dropped after a dry run showed no speed gain).
* `unit_agreement.csv`: run-to-run versus sampler-to-sampler agreement of unit posterior means for the DP models (`unit_agreement.R`), added when complete.

## Findings (2026-11-02)

* **Aggregates:** all 8 dataset x model combinations pass the pre-registered test (|z| <= 3 against eco; max |z| = 1.8).
* **Unit level, `unit_agreement.csv`** (correlation of unit posterior means of W1; A/B = independent pairs of runs):

  | data | model | eco A vs B | Gibbs A vs B | Gibbs vs eco | verdict |
  |---|---|---|---|---|---|
  | census | CAR DP | 0.979 | 0.989 | 0.995 | agree |
  | census | NCAR DP | 0.777 | 0.995 | 0.929 | limited by eco's own run-to-run noise |
  | reg | CAR DP | 0.994 | 0.998 | 0.909 | **systematic difference** (RMSD 0.065 vs 0.01-0.02 within sampler) |
  | reg | NCAR DP | 0.963 | 0.983 | 0.916 | **systematic difference** (RMSD 0.063 vs 0.03-0.04) |

* **reg CAR DP diagnostics (`../diagnostics/`):** the Gibbs answer is stable across seeds, initial clusterings (1 or 30) and truncation (K = 30, 100): alpha about 1, 5-6 clusters. eco settles at alpha about 0.7, 4 clusters, also at 22000 draws. With the true rates treated as observed, one Gibbs run matches eco (alpha 0.90 vs 0.84-0.88, 5.2 vs 5.0 clusters) and another does not (alpha 1.9, 8.5 clusters), so the blocked sampler's mixing over the number of clusters is seed-dependent. The priors were checked against eco's `gibbsDP.c` and match. **Unresolved:** which sampler is right on reg. The next test is simulation-based calibration of both samplers.

## Calibration of the exact-DP sampler (2026-11-02)

Plan pre-registered in `docs/bayes/calibration_plan.md` (commit 426115a) before any run. Files: `calibration_raw.csv`, `calibration_metrics.csv`, `calibration_report.md`, `calibration_ecdf_<design>.png`. Run: `calibration.R 4 1`, 51 minutes on 4 cores.

| design | data / fit | reps | primary criterion | notes |
|---|---|---|---|---|
| D1 | exact DP CAR / collapsed, K = Inf | 4000 | **calibrated** | all 9 quantities equivalent on 90% coverage; none miscalibrated. alpha and clusters equivalent on every metric. |
| D1b | same data / blocked, K = 50 | 2000 | **not calibrated** (as expected) | clusters miscalibrated on 8 of 9 metrics: 90% coverage 0.74, right-tail misses 0.20, mean rank 0.555 (true count above the posterior). alpha inconclusive. Unit rates and means fine. |
| D2 | parametric CAR / K = 1 | 2000 | **calibrated** | some metrics inconclusive (power), none miscalibrated. |
| D3 | exact DP NCAR / collapsed | 2000 | **not calibrated (inconclusive)** | nothing miscalibrated; W1_u1 and W1_u3 inconclusive on 90% coverage (0.878, 0.876; CIs reach 0.866 and 0.864, below the 0.87 margin). All ECDF curves stay inside the simultaneous 95% bands. Lower ESS (23% of fits < 100 for unit rates, vs 5% in D1). |

**The cluster-count weakness is fixed (pre-registered claim).** In D1, alpha and the cluster count are equivalent on 90% coverage (0.893 and 0.890) and not miscalibrated on any metric. D1b, on the same data, reproduces the old failure clearly, so the test had power to detect it.

**D3 is not a pass.** The data are consistent with mild under-coverage of NCAR unit rates (by about 2 percentage points at 90% and 95%) or with no problem at all; 2000 replications cannot separate the two at a 0.03 margin. Candidate cause (conjecture): the slower mixing in NCAR (lower ESS), since autocorrelated draws give too-narrow central intervals. Not re-run, so as not to add a post-hoc test to the pre-registered verdict.

**Other diagnostics.** Posterior z-scores have mean about 0 and sd 0.98-1.07 for every quantity except D1b clusters (mean -0.36, sd 1.52). The 90% posterior intervals for unit rates are about 96% of the Duncan-Davis width under the DP (CAR and NCAR) and about 83% under the parametric model, so with n = 100 the DP rules out little of the identified interval. ESS for alpha and the cluster count is the weakest point in D1 (median about 85 of 199), but their ranks are uniform.

**Speed (median seconds per fit, n = 100, 5975 sweeps):** collapsed 1.09 (D1), blocked 1.58 (D1b), parametric 0.35 (D2), NCAR collapsed 0.93 (D3).

Note on `ks_D`: the reported 90% CI is a bootstrap percentile interval of a supremum, which is biased upward, so the point estimate can fall below it; only the upper bound is used, as pre-registered.

## Corrected eco with efficient code (2026-11-02)

See `../eco_fixed/README.md`. In short: per sweep, the speed-up over `eco` is entirely implementation (efficient Algorithm 1 is even a little faster than the collapsed sampler); the collapsed sampler's advantage is that it does not get trapped (census). The reg discrepancy is mostly slow mixing of `eco`'s uniform-proposal MH step for the rates, not B1.
