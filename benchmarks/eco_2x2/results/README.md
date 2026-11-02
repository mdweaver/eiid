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
