# Corrected `eco` with efficient code versus the collapsed sampler

`eco_alg1.cpp` re-implements `ecoNP`'s algorithm (Imai, Lu & Strauss 2008, Appendix B; `eco` `src/gibbsDP.c`) with efficient code, and with `eco`'s two departures from its stated model fixed: B1, the Wishart degrees of freedom; and B2, the α mixture weight. The algorithm is Neal's (2000) Algorithm 1 over units, plus the Bush–MacEachern remixing step, an independence Metropolis step with a uniform proposal for each unit's rates (`eco`'s default `rMH`), and the Escobar–West α update.

| Sampler | Configuration update | Rate step |
|---|---|---|
| `eco_original` | `eco` 3.1 as distributed (B1, B2 present) | uniform MH |
| `eco_units` | Algorithm 1, loop over units as `eco` codes it ($n^2$) | uniform MH |
| `eco_grouped` | Algorithm 1, the identical sum taken over clusters ($nk$) | uniform MH |
| `eco_slice` | as `eco_grouped` | slice |
| `collapsed` | `eiDPM`: Algorithm 3 | slice |

**Exactness.** All three corrected variants pass the brute-force test over all 52 partitions of five units (`verify_exact.R`, `verify_exact_*.txt`), as the collapsed sampler does.

## Per-sweep cost (`compare.R`; ms per sweep, one core)

| | eco_original | eco_units | eco_grouped | eco_slice | collapsed |
|---|---|---|---|---|---|
| reg CAR | 13.2 | 1.32 | 0.13 | 0.21 | 0.30 |
| reg NCAR | 19.0 | 1.45 | 0.14 | 0.20 | 0.32 |
| census CAR | 185 | 15.7 | 0.45 | 0.71 | 0.91 |
| census NCAR | 270 | 20.1 | 0.46 | 0.70 | 0.99 |

Decomposition of `eco_original` against `collapsed`:
* `eco`'s implementation (repeated factorisations and allocations): about 10–13×.
* Summing over units instead of clusters: about 10× on reg, 35–45× on census.
* The collapsed sampler is about 1.4× *slower* per sweep than the efficient Algorithm 1 with the same slice step. Its predictive-t evaluation and its re-caching after each move cost more than a normal density.

So, per sweep, the whole speed advantage over `eco` comes from implementation, none from Algorithm 3.

## Mixing (`compare.R`, then `chains.R` with more and longer chains)

* **census: Algorithm 1 gets trapped.**
  * Of 8 chains of 22,000 sweeps, `eco_slice` had 3 (CAR) and 4 (NCAR) spend the entire run in a state with about 2.5 clusters and α ≈ 0.35, never visiting k = 1.
  * The other chains, and all 8 collapsed chains, agree: k ≈ 1.15, α ≈ 0.16, k = 1 about 85% of the time.
  * The trapped chains shift the CAR aggregate from 0.648 to 0.638, and the NCAR aggregate from about 0.72 to 0.67–0.69.
  * Both samplers target the same posterior (exactness tests), so this is metastability of Algorithm 1, not a different answer. Neal (2000, p. 5) describes the mechanism: explicit cluster parameters can only change one unit at a time. In this compare.R run, the same trap gives R-hat 1.13–1.16 and effective samples per second 100–300× below the collapsed sampler for the aggregate and α.
  * The trapped state is very probably low-mass (heuristic: no collapsed chain and no untrapped chain ever visited it), but this is not proved.
* **reg: the uniform MH rate step mixes slowly.**
  * At 22,000 sweeps, `eco_units` and `eco_grouped` differ from the collapsed sampler: unit-mean correlation 0.88–0.91, RMSD 0.06–0.07, α 0.84–0.93 vs 1.11, about 4.5 vs 5.7 clusters.
  * `eco_slice` agrees exactly (correlation 1.000, RMSD 0.004).
  * At 220,000 sweeps, the MH chains approach the answer (pooled correlation 0.993, RMSD 0.018), but one chain still sits at 5.15 clusters against 5.62.
  * On reg, with the slice step, Algorithm 1 is as efficient per second as the collapsed sampler or better: aggregate 86 vs 50 effective samples per second, α 345 vs 265.

## What this changes

* Earlier, the reg discrepancy between `ecoNP` and `eiDPM` (unit correlation about 0.90) was attributed mainly to B1. **Revised: it is mostly slow mixing of `eco`'s uniform MH rate step**, which the corrected sampler reproduces (0.88–0.91 at a similar run length). The B1 patch alone moved `eco` only from 0.901 to 0.932. B1 remains a real bug (proved), but its contribution on reg is the smaller part.
* **Reason to prefer Algorithm 3:** robustness of mixing (no trapping on census), not speed per sweep. An efficient Algorithm 1 with a slice step is a valid and slightly cheaper alternative on data like reg, but can trap on data like census.
