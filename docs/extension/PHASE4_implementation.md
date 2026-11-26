# Phase 4: implementation in eiPartialID 0.2.0, 0.3.0 and 0.3.1

*How the theory of Phases 2a, 2b and 3 maps onto the R code, the defaults chosen, and how the implementation was validated. The original 2×2 functions (`generateBounds`, `evaluateBounds`, `bounds`) are unchanged. Release 0.3.0 adds M-SIGN (PHASE5a, experimental) and the package side of the response to the batch-B referee reports (`reviews/B_PI_RESPONSE_PLAN.md`); the changes are listed in §0.3.0 below and in `package/eiPartialID/NEWS.md`. Release 0.3.1 is the package side of the round-2 referee reports (`reviews/ROUND2_PACKAGE_PLAN.md`): the projection critical value uses $q=\dim\beta$ (the rank guard is a safeguard, the rank of $\hat V$ a printed diagnostic); `proj_actual` takes $\bar m$ from a projection interval for the rung alone (P3.5(c1)), with a note when that interval lies within $[0,.05]$ or $[.95,1]$; `print()` states which result covers each fit's intervals.*

## Functions

| Function | Purpose | Theory |
|---|---|---|
| `eiBounds()` | bounds and CIs under one rung | 2a Prop. 7 (LP), 2b Thm 2, P3.5 |
| `ladder()` | all rungs, most to least plausible | 2b Thm 2(c) |
| `sensitivity()` | bounds as a function of $\eta$ | P3 Prop. S1 |
| `breakdown()` | largest $\eta$ at which a conclusion holds | P3 §7 |
| `eiBenchmark()` | coverage and width on data with known cells | — |
| `eiRedundancy()` | what each sign / CR restriction buys; associations under the rung alone; CR branches; descriptive EGGHM quantities; weight concentration (0.3.0) | 5a R2, R3 |
| `eiModules()` | registry of restriction types with source, form, sharpness, dial and proof status | PD §1.1 |
| `print`/`summary`/`plot` | reporting; `summary` adds JKST CI_x | — |

Internal files:
* `R/rxc_core.R`: parametrisation of $K$, regression designs, covariance, composition range, constraint builder, dual LP.
* `R/rxc_infer.R`: setup, candidate collection, CLR and conservative CIs, unit-level allowance, tests.
* `R/rxc_api.R`: user interface.
* `R/rxc_modules.R`: registry of restriction types (parse, rows, label, meta).
* `R/rxc_sign.R`: M-SIGN association forms, CR branches and slab, redundancy verdicts, weight concentration.

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
* **Restrictions:** `zeros`, `order` (mean-level group ordering), `sign` (direction of a contextual effect), `slope` (bounded slopes), `linear` (generic), `external` (benchmark interval on a target). Each adds rows to the LP. `assoc` (signs on $\delta_B$, $\delta_W$) and `cr` (contextual reinforcement) are the M-SIGN module (PHASE5a), **experimental in 0.3.0**: a sign is one row; CR is the W1/2×2 slab (one LP) or a union of two branches (W2, pair). W2 signs one pooled sum over groups. Under grouped rungs each sign is one pooled restriction across groups; the rows use each unit's own group's $H$ (the district functional), but R3/R10 are stated for the pooled model only (R16 open), so the fit prints a note and drops the "sharp" line.
* **Variance:** HC3 by default, falling back to HC1, HC0, then classical. The estimator used is reported, with a note when a fallback occurred. *(Review response, WS5: the `cluster` argument and the CR1 branch were removed. Inference is conditional on the design and the group effects, P3 §0.4.)*

## Inference defaults

* **CLR-type intervals (default, algorithmic).** Candidate endpoints are the dual vertices at $\hat\beta$ and at `ndraw = 100` Gaussian draws. Critical values come from `nz_draw = 4000` draws. With the simulated candidate set this is Conjecture P3.5(b3); P3.5(b1) is the theorem with the full vertex set. If no draw yields a feasible LP, candidates come from relaxed rows (M-SIGN).
* **Projection intervals alongside (`proj_model_*`, `proj_actual_*`).** These are built from the same candidates, with critical value $\sqrt{\chi^2_{q,1-\alpha}}$, $q=\dim\beta$ from 0.3.1 ($q=\operatorname{rank}\hat\Sigma$ in 0.3.0; a numerical rank does not give an exact reduction, P3 §5) (P3.5(b2)(ii), proved by weak duality). They are asymptotically valid for any candidate set when $\hat\Sigma$ is nonsingular in the directions the candidates use; on `eco::census` they are about 0.005–0.01 wider per side. Decision D4, option 1.
  * **Rank guard (0.3.0; P-B40).** A candidate whose gradient has a component (relative size above $10^{-6}$) in the numerical null space of $\hat V$ (eigenvalues at most $10^{-10}$ of the largest) is dropped from the projection interval; an end with no candidate left is the logical range. **From 0.3.1** the critical value uses $q=\dim\beta$, not the numerical rank: a numerical cutoff does not establish an exact reduction (P3 §5, remark after (b2)(ii); round-2 PI referee, comment 3), so the guard is a safeguard that can only widen the interval, and the rank of $\hat V$ is printed as a diagnostic. A note reports the number of candidates dropped. Convention: when $\hat V$ is identically zero (all residuals zero) every candidate is exact. The guard never fired in the WS6 Monte Carlo (P3 §8: no numerical rank drop in any design).
  * **`proj_actual` (0.3.0; WS6).** Its Bhatia–Davis variance bound now takes $\bar m$ from the projection interval built from the same candidates (at level $\alpha_v$), not from the CLR interval, so the projection column rests on P3.5(b2) alone. **From 0.3.1** that interval for $\bar m$ is computed on the rung's own constraint system (no optional rows), as P3.5(c1) requires; `ci_actual` (CLR) keeps its bound from the rung plus linear rows and is labelled algorithmic. It is slightly wider than in 0.2.0: by 0.00015 per side on `eco::census` (where it is not already clipped by Duncan–Davis) and 0.0003 on the voter_vax gap.
* **Zero-variance CLR candidates (0.3.0; P-B23).** In `.rxc_clr_min` a candidate with $\hat\sigma=0$ counts as exact only if its gradient is zero; a candidate with a nonzero gradient in a null direction of $\hat V$ is dropped (each candidate is a valid bound at the truth by weak duality, so the minimum over the rest is valid and wider). If none is left the end is $\pm\infty$, i.e. the logical range.
* **Model-implied rate:** α/2 per side; the interval is intersected with the logical range [0, 1].
* **Actual rate:** α split 0.04 / 0.005 / 0.005 (estimation / unit-level term / variance bound), intersected with Duncan–Davis.
  * **Bhatia–Davis allowance (default):** by concavity, $\sum_i\pi_im_i(1-m_i)\le\bar m(1-\bar m)$, where $\pi_i=w_i^2/\sum w^2$. The weighted mean $\bar m$ is itself bounded by an interval of the same version as the reported CI: a CLR interval for `ci_actual` under `ci = "clr"`, a projection interval for `proj_actual`, and a conservative interval under `ci = "conservative"` (0.3.0, decision 3(a), so each column rests only on its own result). This is never larger than Popoviciu's 1/4; it is smaller only when the plausible range of $\bar m$ excludes 0.5.
  * `fpc = "popoviciu"` reproduces JKST and returns the variance-bound share of α to estimation.
* **`jkst = TRUE`** adds CI_x with Popoviciu $S_1$ and the revert-to-DD rule.
* **Tests:**
  * the emptiness test (conservative inflated set at level α);
  * for M1/L4 with groups, a Wald test of the overidentifying restrictions against the M2 fit, df $(G-1)(C-1)R(R-1)/2$. **From 0.3.0** it is calibrated by the same restricted wild bootstrap as the pooling test (null fit M1, refit M2; decision 3(b)); size in design S0 (19 restrictions): 0.054 (bootstrap), 0.057 ($\chi^2$), 0.055 (F);
  * the **pooling test** (P3.6(c)): a Wald test that all group shifts in the M1 fit are zero, df $(G-1)(C-1)R$, with a flag when some group's compositions are not affinely spanning. **From 0.3.0 its p-value is a restricted wild bootstrap** ($B=399$): fit the regression with all shifts zero, draw $Y^*=\hat Y_0+\hat e_0v$ with Rademacher $v_i$ shared across the $C-1$ equations, refit M1 and recompute the same HC Wald statistic; $p=(1+\#\{W^*\ge W\})/(B+1)$. The $\chi^2$ and F ($W/\text{df}$ against $F_{\text{df},n-k}$) p-values are also returned (`p.value_chisq`, `p.value_F`). The random-number stream is restored afterwards, so no other output changes. Size in the WS6 design S0 (38 restrictions): see §0.3.0.
* **Notes added (WS5):**
  * near-boundary flag: strict-feasibility slack below 0.005 at $\hat\beta$ (P3.2's Slater condition);
  * few units per group (median below 30): CLR coverage rests on simulation (P3 §0.5);
  * unbounded estimate: a composition range is not full-dimensional (2a Lemma 3(b));
  * CLR is algorithmic.
* **Rank errors** name the groups that break the rank condition: A5 per group for M2; $R$ affinely independent compositions per group for M1 (2b Lemma D1(d)); for M0, a hint that for $R=3$ six compositions off a common conic are needed.

## Validation

* **Exact replication:** with `range_offset = 1e-5` and `vcov = "HC1"`, the plug-in bound and every CI_x equal `generateBounds()`/`evaluateBounds()` on `eco::census`, to $10^{-10}$.
* **Independent implementation:** plug-in bounds for M0, M2, M1 and M1(δ = 0.05) on a 3×3, 6-group example match the Python code of Phases 2a/2b to $10^{-9}$.
* **Properties:** the three interfaces agree; the closed form of 2a Prop. 5 holds; L6 equals the no-intercept regression; restrictions only narrow the bound; a true structural zero is compatible with the data; the HC3 → HC1 fallback works; the ladder is nested; the M1 test rejects group-specific slopes and not M1; CIs cover in a small simulation; sensitivity bounds widen with η.
* **`R CMD check`:** 0.2.0 was clean apart from a NOTE inherited from 0.1.2 (MASS and eco in Imports but used only in examples). In 0.3.0 `MASS` is dropped (unused) and `eco` moved to Suggests; examples and vignettes that use `eco::census` run only when it is installed.
* **Real data with known truth:** NC registration by race (`eiPack::senc`, 212 precincts in 8 counties).
  * With the very small Native American row kept, the plug-in is empty: that row's composition range is tiny and noisy. The emptiness test does not reject.
  * With that row merged, every model-implied and actual-rate CI contains the true white→Dem (0.468) and black→Dem (0.893) rates on every rung tried. The M1 plug-in estimate does not, which illustrates P3 Thm 4.

## Second release (agreed)

* the series-based sharper bound ($\Theta_X$ via series regression of unit-level DD limits), for pooled designs;
* quadratic contextual effects as a weaker rung (L1′);
* poll-level ordering and bounded dispersion;
* difference-in-differences wrappers;
* the optional non-negative-correlation assumption for the unit-level term.

## 0.3.0: M-SIGN (experimental) and the batch-B response, package side

Decisions U1(a)–U8(a) of `reviews/B_PI_RESPONSE_PLAN.md` §5 and its three WS6 items. No estimate changes; the only numeric change is `proj_actual` (above).

* **Labels (WS-B0).** `assoc` and `cr` have proof status "experimental", so `print()` shows an EXPERIMENTAL line. The refutability wording names the base model: a $\delta_B$ sign and CR are always feasible only for continuous tables without a model; with the linear rung they can be rejected jointly; with integer counts or other restrictions CR can be refuted. An empty plug-in is called a sample diagnostic; the notes say whether the rung alone is rejected (emptiness test) and whether the rung alone has an empty estimate. `eiRedundancy()` says "empty at the point estimate (see the emptiness test)" instead of "refuted", and "cuts the set, not this target" is documented as zero gain for this target, not implication. The projection label names its conditions (nonsingular $\hat V$, fixed dimension). The "sharp for the regression" line is printed for the pooled model only; grouped rungs print "outer set; sharpness not yet proved" and a note that each sign is one pooled restriction across groups (U2(a)). W2 is documented and labelled as one pooled sum (U8(a)).
* **CR intervals (U6(a)).** CLR (branch union) stays the default; the `ci_model`/`ci_actual` lines say "CLR, algorithmic", and the projection branch union is printed beside them.
* **Reporting (U7(a)).** `print()` shows every piece of a union (the `*_gap_*` columns still hold only the widest gap; `fit$pieces` has all). With restrictions, the rung alone is fitted as well (`base_*` columns, printed beside; with `seed` it equals a separate fit of the rung). `eiRedundancy()$concentration` gives, for each signed association, the largest unit weight and the effective number of units $1/\sum_iu_i^2$, with $u_i\propto\sum_s\lvert W_{is}\rvert$ (for 2×2 $\delta_B$, $u_i\propto N_ix_i(1-x_i)$). Later (0.3.x): the drop-one display and the L5/L6 sign check.
* **Tests:** `tests/testthat/test-batchB.R`.
* **Pooling-test size** (`benchmarks/pooling_test_size.R`, design S0 of P3 §8.3, 2,000 replications): rejection at nominal 5% is 0.078 for $\chi^2$, 0.072 for the F form and **0.043 for the wild bootstrap** (Monte Carlo SE 0.005–0.006; Wald mean 39.6 and variance 95.8 against 38 and 76 for $\chi^2_{38}$). Output: `benchmarks/results/pooling_test_size_output.txt`.
* **Known limits (WS6 Monte Carlo, P3 §8; NZ benchmark):** CLR is algorithmic in every regime; with many small groups (R2) the plug-in under-covers and CLR coverage rests on simulation; in the NZ elections pooled by electorate, M0 is rejected by the pooling test and M1 by its test against M2 (the plug-in M1 estimate is empty, but the emptiness test does not reject M1), and non-empty M1 CLR intervals missed the truth in 5 of 63 cells.


## 0.4.0: the count layer (M-COUNT, build step 1)

*Implementation facts only; the theory is `PHASE6a_count_layer.md` (R19, R27, contract §6), not changed here.*

* **Files.** `R/counts_core.R` (margins, rows, exact rounding, TU classes, LP/MILP, witnesses checks), `R/counts_api.R` (`eiCounts()`, targets `rate()`, `gap()`, `ovr_gap()`, `share()`, `combo()`, `witness()`, print, the R27 block); `eiBounds()` in `R/rxc_api.R` embeds the layer (`target`, `counts`, `integer`). Tests: `tests/testthat/test-counts.R`.
* **After the code review** (`reviews/PHASE6a_code_review.md`; all 16 issues fixed): the MILP is GLPK via the suggested package `Rglpk` (time limit per program `time_limit`, objective scaled to integers from the exact weights, exact check, certificate by a second program with $c'n\ge v^*+1$ that must be infeasible; GLPK status 1 counts as infeasible only if the LP relaxation of the same system is infeasible). Rows that cannot be scaled to integers (coefficients above $2^{40}$, e.g. pairwise $\delta_B$ and W2 at realistic $N_i$) are passed normalised with their bounds relaxed outward by $10^{-6}$, so GLPK solves a superset: an exactly feasible optimum of the superset is optimal and an infeasible superset proves infeasibility; otherwise the status is "integer status unresolved". Joint LPs use GLPK's simplex when available (lpSolve otherwise). `HiGHS` (CRAN `highs`) was requested; installing it from CRAN was refused by this session's permission policy, and `Rglpk` was available from the Ubuntu repository; it meets the same requirements (time limit, MIP status). The P3.5(a) maximum now runs over the beta-dependent terms only (proof review round 3, nit a).
* **Solver (first version).** `lpSolve` (already imported): simplex with sparse (`dense.const`) input, so every LP solution is a basic optimal solution, i.e. a vertex (R19(d)); `all.int = TRUE` for the MILP (branch and bound). Problems in which every row involves one unit (margins, boxes, orderings, dials, known cells) are solved unit by unit; margins-only single cells and two-row one-column contrasts use the closed forms of Lemma 2.1(a)-(b) (witnesses then computed on demand). No new solver dependency.
* **Exact arithmetic (C4, RC9).** New import `gmp`. Margins count as integers only if every count in the data is exactly integer (data given as shares, or rescaled margins, are not counts). With integer margins all thresholds ($a b/N$, $Y^{NM}$, $Y^{ER}(\varepsilon)$, slack shifts, pairwise and W2 coefficients) are `bigq` rationals; slack and user bounds enter as the exact rationals of the given doubles. Rounding (§1.4): bounds and totals to $\lceil\ell\rceil,\lfloor u\rfloor$; other rows scaled to coprime integer coefficients first (if a scaled coefficient exceeds $2^{40}$ the row is passed unscaled and every MILP witness is still checked exactly). Every integer witness is verified against the margins and every original row in exact arithmetic; a failed check of a rounded-LP vertex falls back to the MILP; a failed MILP check is an error.
* **Odds-ratio dial.** The continuous LP uses the root of $k n_{22}=\psi n_{12}n_{21}$ (double); with integer margins the integer bound is found by exact search around the root, and the exact checks use it. Defined for 2 x 2 tables only; ignored in units with a zero margin (6a §1.3 convention).
* **"auto" (C4).** Whole row set per CR branch: generic rows (pairwise orderings, pairwise and W2 signs, `linear` rows) → MILP; otherwise totals with at least two units are tested for classes (i) (bipartite crossing graph), (ii) and its transpose, (iii) and its transpose; single-unit totals count as bounds. `"exact"` refuses non-integer margins.
* **Binding rows (C7).** Drop-one test per restriction: a restriction binds an endpoint iff removing it moves the endpoint outward, i.e. iff it removes the whole optimal face (PD R2).
* **Cross-layer block (R27, C8).** Model interval for the model-implied value of the target: the projection interval (P3.5(b2)(ii)) built from fresh candidates, or the conservative interval (P3.5(a)) when the fit has `ci = "conservative"`, at level $\alpha\cdot0.9$, plus Lemma 4.1's allowance at $\alpha\cdot0.1$: Popoviciu $z S_\varphi$ (default) or Hoeffding $t_E(\varphi)$. CLR is not used here, so the inherited status is P (asymptotic) for pooled M0 (§4.3). Representation: for `ovr_gap()` the smaller $S_\varphi$ of the two representations of Remark 4.1a (design only). The random-number stream is saved and restored, so no other output changes. Refused: `ci = "none"` (an estimate), trimmed units, `eta > 0`, a fit on other units; `share()` targets get the count layer only. Wording strings are those of C8 (ASCII fallback "alpha", "-", "<=" in non-UTF-8 locales).
* **`duncan_davis` column (§9 item 3).** `fit$table$dd_lower/dd_upper` and `ci_actual` (intersected with them) are unchanged; new columns `dd_exact_lower/dd_exact_upper` hold the exact margins-only LP, which the printout shows under `target = "actual"`. Reason: `target = "model"` must reproduce earlier outputs, and the cellwise range is still valid (outer).
* **Role of rows (C8).** Default identity for `known` and `total`, belief otherwise; overridable per row with `role =`.

## 0.5.0: changes and DiD (build step 2)

*Implementation facts only; the theory is `PHASE2e_did_general.md` (proof reviewed), not changed here.*

* **Files.** `R/change_core.R` (input and (FP) checks, the stacked problem, targets `did()`, `gap_at()`, `rate_at()`, `bin_did()`, rows on changes, the E27 classifier), `R/change_api.R` (`eiChange()`, closed forms E5/E6/E26, the joint set `eiJoint()`, print, `witness()`), `R/change_model.R` (`fd()`, `lce_periods()`, actual-target intervals, cross-layer, Prop. E25 for `eiATT`), `R/att.R` (`eiATT()`, trend constructors, the joint causal LP, E17/E19/E20/E21 reports, Def. 5.1). Tests: `tests/testthat/test-change-reference.R`, `tests/testthat/test-change-properties.R`.
* **Stacking (Lemma E1).** The T elections are one `eiCounts` problem whose units are the unit-elections; stacked unit $k=(t-1)p_0+i$. Units are matched across elections by label; units with $N_i=0$ (the same in every election under (FP)) are dropped before stacking. One-election problems with the same number type and log are kept for parsing within-election rows and for the exact one-election thresholds.
* **Change rows.** A new row kind `"change"` (coefficients $\pm1$ on one cell at two elections over a unit set; `row$pair`). W1 signs on changes use the difference of the one-election thresholds of `.ec_ovr_thresholds()` (exact `bigq`), with the slack applied once; pair/W2 signs on changes are the difference of the one-election generic rows (kind `"gen"`). Inward rounding is unchanged (coefficients $\pm1$ give scale 1).
* **Integer classes.** `.ec_system()` now calls `P$classify` when a problem supplies one (only the stacked problem does; `eiCounts` is unaffected). `.ed_classify()`: no change rows → 6a's classes on the unit-elections; change rows → E27 test (R = 2; no generic rows; no within-election totals; one column; period pairs of maximum degree 2 without cycles, by union-find; laminar unit sets per pair). The TU vertex is still accepted only if integral and exactly checked (fallback MILP). The E26 lattice is used whenever the target is a one-vs-rest DiD (or the R = 2 pair) and every row is a W1 sign/CR on the same $(r,c,u,t,S)$; its witnesses fix $Y_{rct}(S)$ and $Y_{rcu}(S)$ by two totals (class (i)).
* **Closed forms as checks.** E5 (Z-intervals mapped by $D(Z)$) and E6(b) are computed beside the LP; a disagreement above $10^{-7}$ is a warning.
* **Causal LP.** Rglpk directly, variables $n$ plus $(\theta_t,S_t,b_t)$ per post-period between the base and the last requested post-period (RM couples consecutive post-periods); $S_t$ fixed at 0 when no comparison group is stated (MTR/ND without PS). Continuous tables only (integer causal ends are open, Remark 7.1). E17(b) is evaluated by an LP on the face $\arg\max_F D_t$ (tolerance $10^{-9}$) for each count branch attaining the maximum. E19(b) vacuity: feasibility of the $2L$ LPs, one post-period and $\Lambda=0$ only.
* **Model layer.** Blocks are `.rxc_setup()` fits sharing one design: FD on consecutive pairs ($\Delta T$ as the response, `st$change = TRUE`, which sets the last column's constant to 0 in `.rxc_idlin()`), LCE per election. The stacked $\hat\beta$ gets one joint sandwich covariance from all residuals (`.rxc_vcov`). FD validity rows: $|\Delta m_{rc}(x)|\le1$ at the hull vertices. Intervals: projection (P3.5(b2)(ii), candidates from `.rxc_ci_clr`) or conservative (P3.5(a)) at $0.9\alpha$, plus the Lemma E12 allowance at $0.1\alpha$ computed from the target's representation over stacked cells (smallest $S$ over the available representations, as 6a Remark 4.1a). The random stream is saved and restored.
* **E25.** LP variables $(k,\theta_t,S_t,b_t,e_q)$; quantities $q$: the post DiDs, the consecutive pre DiDs (RM) or the pre DiDs relative to the base (SD), and under LCE also $\theta^{PT}_t$ and $B_{s^*ct}$; $\alpha_{\rm proj}=0.9\alpha$, $\alpha_q=0.1\alpha/|q|$; one projection event (same draws and critical value) for all branches; intervals clipped to $[-1,1]$.
