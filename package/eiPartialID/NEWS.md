# eiPartialID 0.3.1

Package side of the second-round referee reports (Batch B ecological-inference
report; Phases 1-3 report, item 3). Plan: `docs/extension/reviews/ROUND2_PACKAGE_PLAN.md`.

* `proj_actual`: the Bhatia-Davis variance bound now comes from a projection
  interval for the rung alone (no optional rows), so the projection column
  rests only on P3.5(b2)(ii) and a union bound, whatever restrictions are
  added. `ci_actual` (CLR) keeps its CLR-based bound, now labelled as such.
  `fpc = "popoviciu"` gives the bound 1/4 in every column.
* Projection critical value: `sqrt(qchisq(1 - alpha, p))` with
  `p = dim(beta)`, not the numerical rank of the estimated covariance matrix.
  The rank is printed as a diagnostic, and the rank guard stays as a labelled
  safeguard. Unchanged in every full-rank design.
* `print()` states, for each fit, which result covers its intervals: pooled M0
  (P3.5(b2)(ii), or P3.5(a) for the conservative interval; with sign/CR rows
  PHASE5a R10); grouped models only with few, large groups (regime R1); and
  no coverage theorem for sign/CR restrictions under grouped rungs (R16 open).
  The grouped M-SIGN note says the same.
* A note when the Bhatia-Davis allowance's interval for the weighted mean rate
  lies within [0, .05] or [.95, 1]: the allowance relies on that rate being
  bounded away from 0 and 1 (P3.5(c1)), and `fpc = "popoviciu"` is the safe
  choice. An automatic fallback is planned.
* With several optional restrictions, the rung-only block is labelled as
  their joint effect, with a pointer to the added-last gains in
  `eiRedundancy()`.
* The pairwise CR label names both sides: `delta_B(A, B) [pair] and
  Delta_W(A) [W2, pooled over all groups]`.
* Vignette and plain-language guide: the W1 sign is a weighted sum of two
  within associations (either part can be negative; two-precinct example);
  pairwise CR changes only the between side; each coverage figure names its
  procedure (the literal bracket-and-CI intersection covered 304 of 600, the
  reference procedure that widens afterwards 550 of 600); both false-belief
  results (398, 110 and 0 of 400 at s = 0, .08, .25); the projection proof is
  scoped to pooled M0; open items listed.
* No estimate or interval in the package's reference examples changed.

# eiPartialID 0.3.0

## New: signs on associations and contextual reinforcement (M-SIGN), experimental

* `eiBounds(restrictions = list(assoc = ..., cr = ...))` adds signs on the mean
  between association (delta_B) and the within (contextual) association
  (delta_W), and contextual reinforcement (CR: both have the same sign,
  whichever it is). A sign is one linear row; CR is one LP (the W1 / 2 x 2
  slab) or a union of two branches (W2, pair).
* `eiRedundancy()` reports what each restriction buys (face test; added first
  and added last), the associations already settled by the data under the rung
  alone, the CR branches, descriptive EGGHM quantities and, new in this
  release, the weight concentration of each signed association.
* `sensitivity(dial = "slack")` and `breakdown(dial = "slack")` relax the
  beliefs by a slack in probability units. `eiModules()` lists the restriction
  types and their status.
* **The module is experimental.** Two referee reports on it (review batch B)
  are being answered, and `print()` shows an EXPERIMENTAL line for every sign
  or CR restriction. Under grouped rungs (M1, M1delta, M2, L4 with groups) each
  sign is one pooled restriction across groups, sharpness is proved only for
  the pooled model M0, and a note says so at fit time. W2 signs one pooled sum
  over the original groups, not each group's association; a per-group version
  is planned.

## Changes in response to the referee reports

* Wording: a sign on delta_B and CR are always feasible only in the model-free
  problem with continuous tables (no linear rung, no other restriction). With
  the linear rung they can be rejected jointly with it, and with integer counts
  or other restrictions CR can be refuted outright. The old "can never be
  refuted by aggregate data on their own" is gone from the help pages, the
  notes and the module registry.
* An empty estimate is called a sample diagnostic, not a test. Next to a joint
  rejection the notes say whether the rung alone is rejected (emptiness test)
  and whether the rung alone has an empty estimate. `eiRedundancy()` reports
  "empty at the point estimate (see the emptiness test)" instead of "refuted";
  "cuts the set, not this target" means zero gain for this target, not that the
  restriction is implied.
* Reporting: `print()` now shows every piece of a union of intervals (three or
  more included); the `*_gap_*` columns still keep only the widest gap, and
  `fit$pieces` has every piece. With restrictions, the rung alone is fitted as
  well and printed beside the restricted result (`base_*` columns; argument
  `base`). For CR fits the CLR lines say "algorithmic" and the projection
  interval is printed beside them.
* Projection intervals (`proj_*`): a rank guard. Candidate endpoints whose
  gradients use a direction in which the estimated covariance matrix is
  numerically singular are dropped (the interval stays valid and gets wider),
  with a note. The label now names the conditions: a nonsingular covariance
  matrix in a fixed-dimension design (M0, or few large groups). The guard never
  fired in the Phase 3 Monte Carlo.
* CLR intervals: a candidate with zero estimated variance is treated as exact
  only if its gradient is zero; otherwise it is dropped.
* `proj_actual` takes its Bhatia-Davis variance bound from a projection
  interval instead of a CLR interval, so the projection columns rest only on
  the proved result (P3.5(b2)). They are slightly wider (by about 0.0002 per
  side in the examples).
* Pooling test (M0 against M1): the HC3 Wald test against chi-squared
  over-rejected with many groups (9% at a nominal 5% with 38 restrictions in
  the Phase 3 Monte Carlo; the F form, 8%). Its p-value now comes from a
  restricted wild bootstrap (Rademacher weights, B = 399); the chi-squared and
  F p-values are still returned. The bootstrap does not change the random
  number stream used elsewhere. In the same design (2,000 replications) the
  rejection rate at 5% is 0.078 for chi-squared, 0.072 for F and 0.043 for the
  wild bootstrap (Monte Carlo SE about 0.005;
  `benchmarks/pooling_test_size.R`).
* Test of M1's testable part (M1 against M2): calibrated by the same
  restricted wild bootstrap (null fit M1, refit M2); chi-squared and F
  p-values still returned. In design S0 (19 restrictions, 1,000 replications)
  the rejection rate at 5% is 0.057 for chi-squared, 0.055 for F and 0.054 for
  the wild bootstrap (`benchmarks/m1_test_size.R`).
* With `ci = "conservative"`, `ci_actual` now takes its Bhatia-Davis variance
  bound from a conservative interval rather than a CLR interval, so the
  conservative actual-rate interval rests only on the proved result (P3.5(a)).
  Conservative `ci_actual` values change slightly; CLR and projection output is
  unchanged.
* The `cluster` argument was removed in the earlier review response (no
  clustered variance; inference is conditional on the design and the group
  effects).
* `MASS` is no longer imported (it was unused) and `eco` moved to Suggests;
  examples and vignettes that use `eco::census` run only when it is installed.

## Known limits

* The default CLR intervals are algorithmic in every regime: candidate
  endpoints are collected by simulation, and coverage is supported by
  simulation, not proved. The projection and conservative intervals carry the
  proof, for fixed-dimension designs with a nonsingular covariance matrix.
* With many small groups (regime R2 of the Phase 3 Monte Carlo) the plug-in
  estimate under-covers, and the CLR interval's coverage rests on simulation.
* New Zealand elections pooled by electorate (`eiDatasetsTables()`,
  `benchmarks/`): the pooling test rejects M0 and the M1 test rejects M1 (p near
  0; the plug-in estimates are empty, but the emptiness test does not reject M1), and non-empty M1 CLR intervals missed the truth in 5
  of 63 cells. Grouped M-SIGN fits sit on top of these rungs, which is one more
  reason the module is labelled experimental.
* Not yet available: the drop-one display (each belief removed in turn), a
  check of entered signs at the point estimate for the point-identified rungs
  L5 and L6, and per-group W2 directions.

# eiPartialID 0.2.0

* `eiBounds()`: R x C tables, grouped units (M1, M1delta, M2, L4, L5, L6), the
  assumption ladder (`ladder()`), CLR-type and conservative confidence
  intervals, projection intervals alongside CLR, sensitivity and breakdown
  analysis, optional restrictions, `eiBenchmark()` and `eiDatasetsTables()`.
* The original 2 x 2 functions (`generateBounds()`, `evaluateBounds()`,
  `bounds()`) are unchanged; `eiBounds()` reproduces them in the 2 x 2 case.
