# Phase 4: implementation in eiPartialID 0.2.0

*How the theory of Phases 2a, 2b and 3 maps onto the R code, the defaults chosen, and how the implementation was validated. The original 2×2 functions (`generateBounds`, `evaluateBounds`, `bounds`) are unchanged.*

## Functions

| Function | Purpose | Theory |
|---|---|---|
| `eiBounds()` | bounds and CIs under one rung | 2a Prop. 7 (LP), 2b Thm 2, P3 Thm 5 |
| `ladder()` | all rungs, most to least plausible | 2b Thm 2(c) |
| `sensitivity()` | bounds as a function of $\eta$ | P3 Prop. S1 |
| `breakdown()` | largest $\eta$ at which a conclusion holds | P3 §7 |
| `eiBenchmark()` | coverage and width on data with known cells | — |
| `print`/`summary`/`plot` | reporting; `summary` adds JKST CI_x | — |

Internal files:
* `R/rxc_core.R`: parametrisation of $K$, regression designs, covariance, composition range, constraint builder, dual LP.
* `R/rxc_infer.R`: setup, candidate collection, CLR and conservative CIs, unit-level allowance, tests.
* `R/rxc_api.R`: user interface.

## Models (rungs)

| `model` | Rung | Regression | Unknowns |
|---|---|---|---|
| `"M0"` | pooled LCE | quadratic monomials | $K$ ($d_0$) |
| `"M2"` | group-specific LCE | monomials × group | $K^{(g)}$ for each group |
| `"M1delta"` | shared LCE, deviations ≤ `delta` | as M2 | $\bar K$ and $K^{(g)}$, linked by $\lvert k^{(g)}-\bar k-c_g\rvert\le\delta$ |
| `"M1"` | shared LCE, group baselines | monomials + group × $X_r$ shifts | common $K$ |
| `"L4"` | M1 (or M0) + `row_constant` has no contextual effect | as M1 / M0 | $K$ with equality rows |
| `"L5"` | Goodman within groups | $X$ × group, no intercept | none (point) |
| `"L6"` | Goodman | $X$, no intercept | none (point) |

## Inputs, ranges and options

* **Inputs:** a formula `cbind(cols) ~ cbind(rows)` with counts; `data` with `rows` and `cols`; or shares `x`, `t` and sizes `n`. Unequal row and column totals stop with an error (A1) unless `mismatch = "rescale"`.
* **Composition range:**
  * `"hull"` (default): group-wise sample hull, which is conservative (P3 Prop. P3.2(c)).
  * `"trim"`: by quantile, minimum size or explicit `drop`. Trimmed units enter the actual rate through their Duncan–Davis bounds only, and the model-implied rate refers to the kept units.
  * `"extrapolate"`: pure units or a user matrix, with a warning.
* **Restrictions:** `zeros`, `order` (mean-level group ordering), `sign` (direction of a contextual effect), `slope` (bounded slopes), `linear` (generic), `external` (benchmark interval on a target). Each adds rows to the LP.
* **Variance:** HC3 by default, falling back to HC1, HC0, then classical. The estimator used is reported, with a note when a fallback occurred. `cluster` gives CR1 variances.

## Inference defaults

* **CLR-type intervals**, with candidate endpoints collected as dual vertices at $\hat\beta$ and at `ndraw = 100` Gaussian draws (P3 Thm 5(b)). Critical values come from `nz_draw = 4000` draws. If no draw yields a feasible LP (e.g. an empty plug-in), the conservative inflated set is used.
* **Model-implied rate:** α/2 per side; the interval is intersected with the logical range [0, 1].
* **Actual rate:** α split 0.04 / 0.005 / 0.005 (estimation / unit-level term / variance bound), intersected with Duncan–Davis.
  * **Bhatia–Davis allowance (default):** by concavity, $\sum_i\pi_im_i(1-m_i)\le\bar m(1-\bar m)$, where $\pi_i=w_i^2/\sum w^2$. The weighted mean $\bar m$ is itself bounded by a CLR interval. This is never larger than Popoviciu's 1/4; it is smaller only when the plausible range of $\bar m$ excludes 0.5.
  * `fpc = "popoviciu"` reproduces JKST and returns the variance-bound share of α to estimation.
* **`jkst = TRUE`** adds CI_x with Popoviciu $S_1$ and the revert-to-DD rule.
* **Tests:** the emptiness test (conservative inflated set at level α) and, for M1/L4 with groups, a Wald test of the overidentifying restrictions against the M2 fit.

## Validation

* **Exact replication:** with `range_offset = 1e-5` and `vcov = "HC1"`, the plug-in bound and every CI_x equal `generateBounds()`/`evaluateBounds()` on `eco::census`, to $10^{-10}$.
* **Independent implementation:** plug-in bounds for M0, M2, M1 and M1(δ = 0.05) on a 3×3, 6-group example match the Python code of Phases 2a/2b to $10^{-9}$.
* **Properties:** the three interfaces agree; the closed form of 2a Prop. 5 holds; L6 equals the no-intercept regression; restrictions only narrow the bound; a true structural zero is compatible with the data; the HC3 → HC1 fallback works; the ladder is nested; the M1 test rejects group-specific slopes and not M1; CIs cover in a small simulation; sensitivity bounds widen with η.
* **`R CMD check`:** clean apart from a NOTE inherited from 0.1.2 (MASS and eco are in Imports but used only in examples).
* **Real data with known truth:** NC registration by race (`eiPack::senc`, 212 precincts in 8 counties).
  * With the very small Native American row kept, the plug-in is empty: that row's composition range is tiny and noisy. The emptiness test does not reject.
  * With that row merged, every model-implied and actual-rate CI contains the true white→Dem (0.468) and black→Dem (0.893) rates on every rung tried. The M1 plug-in estimate does not, which illustrates P3 Thm 4.

## Second release (agreed)

* the series-based sharper bound ($\Theta_X$ via series regression of unit-level DD limits), for pooled designs;
* quadratic contextual effects as a weaker rung (L1′);
* poll-level ordering and bounded dispersion;
* difference-in-differences wrappers;
* the optional non-negative-correlation assumption for the unit-level term.
