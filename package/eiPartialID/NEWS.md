# eiPartialID 0.5.0

Build step 2 (DiD): changes and difference-in-differences for T >= 2 elections on the same
units. Theory: `docs/extension/PHASE2e_did_general.md` (2e; proof reviewed).

* `within_pt()` now defaults to the realised reading (Prop. E35, RWPT-tau): in every supported
  precinct the veterans' untreated change equals the local non-veterans' realised change up to
  `tau_i = M sqrt(1/a_iV + 1/a_iN)` (or a user `tau`), with one spillover `|S_i| <= Lambda` per
  precinct; unsupported precincts get `theta_i` in [0, 1] and `pi_U` is printed. Identification
  set of the count layer (continuous tables), nested in `M`, with the breakdown value `M*` from one
  LP. `lce = TRUE` adds the fit's model layer only through E13 (a row `D_W^S(n)` in the FD-LCE
  interval, E12 allowance `rho_i = 4 w_Vi`) and the result is labelled a CONFIDENCE SET with its
  coverage assumptions. The printout notes that `M` is a tolerance holding in every precinct at
  once, not a confidence level. The expected-rate version is `within_pt(target = "expected")`.
* The E29 confidence interval of `within_pt(target = "expected")` is widened by
  `pi_S (h M + Lambda)` when the local dial `M` or the spillover limit is positive (E29(d)), and
  the printout states the widening. `eiATT()` documents that `$pieces` is the set and
  `$lower`/`$upper` are its envelope ends.
* New `eiChange()`, the descriptive layer. Under a fixed population the T elections stack into
  one count-layer problem over the unit-elections (Lemma E1), so the machinery of `eiCounts()`
  (exact rationals, inward rounding, TU classes, certified GLPK MILP with `time_limit`,
  witnesses, CR unions) applies to rows that couple elections.
  * Targets `did()`, `gap_at()`, `rate_at()`, `bin_did()`; default: the event-study DiDs
    relative to `base`. Margins-only ranges (Prop. E3) and the joint event-study set: the
    printout says whether the corner of the box of separate ranges is attainable (E3(c));
    `eiJoint()` tests any vector of values.
  * Count rows on changes: `sign_change` (W1 one-vs-rest = a half-line in
    Z = Y_rct - Y_rcu with the closed form of Prop. E5(b); pair and W2 as LP rows),
    `cr_change` (a union, E5(d), (f)), `order_change`, `contrast_change` (Prop. E6, with the
    corrected emptiness criterion including Dbar+ < 0, review issue 2), `linear` rows across
    elections, and every `eiCounts()` row at one election (`election =`). Closed forms are
    printed beside the LP and cross-checked.
  * Integer handling: the lattice of Prop. E26 for one-vs-rest DiDs with W1 rows on the same
    changes ("integer-exact (E26 lattice)", any R); the TU class of Thm E27 (R = 2, change rows
    on one column, paths of period pairs, laminar unit sets; class "E27"); 6a's classes on the
    unit-elections; otherwise the certified MILP. Non-integer margins switch the check off.
  * Model layer `model = fd(...)` (FD, Lemma E7; rows R1, R2, R4, delta for R = 2, DiD-S with the
    E9(a) closed form, model signs on changes; plug-in sharp set of E8(d) for 2 x 2) or
    `lce_periods(...)` (per-period pooled LCE with model signs and CR on changes, E11; grouped
    rungs flagged R16). Actual-target intervals: projection interval (P3.5(b2)(ii)) for the
    model DiD plus the two-election allowance of Lemma E12 (rho_i = 2(w_ir + w_is) for a DiD),
    intersected with the count set (Thm E13). Unions are reported whole (Prop. E24).
  * Gap register (2e Sec. 8): different group counts across elections (G1-G3, G5, G9),
    outcome totals that differ from the group totals, and different unit sets (G6) are refused
    with a message; `margins = "estimated"` (G7) warns and switches the integer check off.
* New `eiATT()`, the causal layer, on an `eiChange()` fit.
  * Labelling rule (Def. 5.1): refuses a causal label unless the population, the treated group,
    the comparison group (when the trend assumption uses one), the trend assumption, the
    spillover allowance and no anticipation are stated.
  * Joint causal LP over the realised tables and (theta_t, S_t, b_t) (Thm E16) with the rate
    limits; trend assumptions `no_trend()` (E15), `trend_dial()`/`parallel_trends()` (E17:
    outer form, the exact E17(b) criterion, breakdown value), `rel_mag()` (E18: union of 2L LPs;
    implied M = Mbar sup|D_pre| always printed; degeneration report of E19: collapse to a dial,
    no information, the vacuity criterion, the margins-only artefact), `smoothness()` (Def. 5.3;
    E20 closed form checked), `mtr_nd_ps()` (E21 with the spillover dial; closed form E21(c)).
    `rate_limits = FALSE` gives the outer set (to reproduce NS N21).
  * Spillover dial Lambda; the printout notes that pre-periods carry no information on it (E22).
  * `pretrends = "model"` (Prop. E25): one joint projection region for the stacked parameter,
    box variables with Lemma E12 radii for every actual quantity, the union bound; under FD
    alone the rate limits are dropped and the set is labelled OUTER.
* Within-precinct and pooled parallel trends in expectation (2e Sec. 5.8-5.9; added after the
  independent review): `trend = within_pt(a_min, M)` and `trend = pooled_pt(M)` in `eiATT()`,
  with no default between them or against the realised `parallel_trends()`.
  * Mean-model LP over expected rates (Props E28, E30): adding-up, the averaged Duncan-Davis ends
    of the law (`law =`, default the observed returns), LCE rows on expected rates from
    `lce = lce_mean(base, post, change, epsilon)` (any two imply the third; baseline-curvature
    allowance epsilon), the rate limits m_iV0 + d_iN (+ b_it - S_i) in [0, 1] in every supported
    precinct (aggregate only under PPT, whose witness is one common untreated rate), spillover
    S_i optimised jointly (E28(f)), the local dial |b_it - b_i,t-1| <= M (E33(a)).
  * Support threshold `a_min`; unsupported precincts bounded with theta_i in [0, 1]; pi_U,
    n_eff = 1/sum w_Vi^2, the local comparison counts and the independence unit always printed.
  * One-LP breakdown values of the sign for M (epsilon fixed) and for epsilon (M fixed).
  * E28(e) FD-LCE projection (I_W + pi_S dkappa; I + s dkappa under PPT).
  * E29: the realised target D_W from the count layer (with the fit's rows) and the Hoeffding
    interval [L - pi_U - r_alpha, U + r_alpha] (`independence = "precincts"` or `"persons"`;
    belief rows carry the P(n not in F) penalty).
  * Reporting rule: every printed interval is labelled IDENTIFICATION SET/BOUND, DESCRIPTIVE
    (realised counts) or CONFIDENCE INTERVAL/SET; E29 covers the realised-count route only and
    the printout says that it does not make the identification bound a confidence interval.
  * Deferred: E31, E32, E33(b)-(c), E34.
* Theory changes after the independent review: `rel_mag()` bounds each election's increment of
  the trend gap; the marginal outer set at horizon h is D_h (+) [-h Mbar sup A - Lambda, +];
  RM branch intervals are merged only where they overlap and holes are printed; the system
  always covers every post-period (one shared table, E16(d), E21(e)); the sharp E17(c)
  breakdown LP; the E16(b) whole-person witness check at the LP ends; E8 is computed for one
  pair only (pair-by-pair stacking would be outer).
* Code-review fixes (`docs/extension/reviews/PHASE2e_code_review.md`): GLPK statuses (an
  unbounded LP is infinite, a failure an error, never EMPTY); infeasible count layers with the
  dial and RM; `post =` keeps later post-periods' constraints; E19(a)-(b) only with one
  post-period; the allowance uses the representation of the model-implied form; numeric unit
  selections refer to data rows; a character election variable needs `elections`; smoothness
  with T = 2; model sign rows without `sign`; wording (bin_did under FD, the E26 path's LP
  cross-check, the vignette's E19 sentence, E25 EMPTY and grouped-rung status).
* Def. 5.1 item 3: a spillover allowance is required only when the trend assumption uses a
  comparison group (MTR and ND alone need none).
* Internal: `.ec_system()` accepts a problem-specific integer classifier (unused by
  `eiCounts()`); `.rxc_idlin()` takes a `change` flag (sum over columns of Delta H is 0).
  Existing outputs are unchanged.

# eiPartialID 0.4.0

Build step 1 (M-COUNT): the count layer. Theory: `docs/extension/PHASE6a_count_layer.md`
(R19, R27; proof reviewed); contract §6 (C1-C10).

* New `eiCounts()`: exact bounds on actual (realised) counts and rates from the unit
  margins and rows on the realised table. Targets: `rate()` (merged groups and sums over
  columns allowed), `gap()` (one-column contrasts, differences over unit sets),
  `ovr_gap()` (the one-vs-rest gap D^r_c), `share()` (outcome-margin denominator; count
  layer only) and `combo()`, each with `where =` unit subsets. Rows: known cells,
  structural zeros, cell boxes, published totals, general linear rows, realised unit
  orderings (one-vs-rest = boxes; pairwise for R >= 3), the 2 x 2 odds-ratio dial,
  realised association signs (W1 as totals; pair and W2 with unequal coefficients) and
  realised CR (a union of two systems). Sign and CR rows print EXPERIMENTAL.
* Integer check `integer = c("auto", "lp", "exact")`. With exactly integer margins,
  "auto" rounds every row inward in exact rational arithmetic (new import: `gmp`), tests
  whether the whole row set lies in one total-unimodularity class of R19(b) (two laminar
  families by 2-colouring the crossing graph; one column or row; R = 2 or C = 2) and then
  reports the LP over the rounded rows as "integer-exact (TU class)"; otherwise every
  endpoint (and CR branch) is a MILP ("integer-exact (MILP)"). Integer infeasibility with
  continuous feasibility is reported as "no whole-number table satisfies these rows".
  Non-integer margins (estimates, or data given as shares) switch the check off with the
  note "exact only for the stated margins".
* `witness()`: a table attaining each endpoint, checked against every margin and row
  (exactly for whole-number tables); a certificate of attainability, not an estimate.
* Binding rows per endpoint (drop-one test).
* Cross-layer block (R27): with `model =` an `eiBounds()` fit (or inside `eiBounds()`),
  the count set is intersected with an actual-target model interval: the projection
  interval (or the conservative one) for the model-implied value plus the residual
  allowance of 6a Lemma 4.1 (Popoviciu S_phi, or Hoeffding t_E with
  `allowance = "hoeffding"`), for multi-cell targets too, with the representation of the
  target chosen from the design. Count rows and model rows never share one LP.
  Disjointness is reported as a test, with the wording of contract C8 for identity rows
  only and for belief rows. The status is that of the model interval (6a §4.3).
* `eiBounds()` embeds the count layer (decisions D1/D2). New arguments `target`
  (`"actual"`, default: the actual-rate block first, Duncan-Davis first; `"model"`: the
  old print order), `counts` (count rows) and `integer`. `assoc`/`cr` restrictions with
  `object = "actual"` go to the count layer; `object = c("actual", "model")` is refused.
  `$table` gains `dd_exact_lower`/`dd_exact_upper` (the exact margins-only LP);
  `dd_lower`/`dd_upper` and every model column are unchanged. The printed
  `duncan_davis` column shows the exact range under `target = "actual"`; it differs from
  the old cellwise column only for combination targets, where the old one is outer
  (6a Lemma 2.1(c), §9 item 3).
* Code review (`docs/extension/reviews/PHASE6a_code_review.md`), all 16 issues fixed:
  * MILP: GLPK (package `Rglpk`, in Imports since round 2) instead of lpSolve's branch and bound, with a
    time limit per integer program (`time_limit`, default 10 s), the objective scaled to
    integers from the exact weights, an exact check of each solution, and a certificate
    of optimality (a second program with c'n >= v* + 1 must be infeasible). On timeout or
    without `Rglpk` the endpoint is the LP over the rounded rows with the status
    "integer status unresolved (...): LP outer bound", never "integer-exact". Problems are
    split into connected components of the unit-coupling graph. Pairwise delta_B and W2
    rows at realistic unit sizes now finish (issues 1, 16).
  * The TU path accepts a vertex only if it is integral before rounding and passes the
    exact check; otherwise it falls back to the MILP (issue 6).
  * `known` values follow the user's unit order; duplicated units and wrong value lengths
    are refused; out-of-range row numbers and units dropped for N_i = 0 are errors with
    clear messages (issues 2, 11).
  * The count block can no longer stop `eiBounds()`: a failure gives a warning and the model
    output; units with N_i = 0 on the x/t/n path are aligned (issue 3).
  * Doubles typed by the user are rationalised (denominator <= 1e6 within 1e-12; printed),
    and character fractions are exact (issue 4). An all-zero pairwise or W2 row is vacuous
    and dropped with a note (issue 5).
  * A `where =` expression that uses an outcome variable makes the target count-layer only;
    sets given as labels or numbers print that the set must be fixed from the design (issue 7).
  * A `total` given as a floor or interval defaults to a belief row; `value =` totals stay
    identity rows (issue 8).
  * Cross-layer statuses follow the 6a §4.3 rows added for L4, L5, L6 and M1(delta); for
    point-identified rungs the model interval is the one-functional normal interval
    (issue 9). The cross-layer block reuses the fit's simulation draws (issue 10).
  * Labels: x/t vectors give X/notX, T/notT (issue 12); delta_B sign rows are labelled
    delta_B (issue 13); long `known` labels are shortened, duplicate target labels are made
    distinct, the binding line says jointly binding rows are not shown, the Lindeberg note
    prints its value and is suppressed for the Hoeffding allowance, and `sign = ">=0"` /
    `"<=0"` are accepted (issue 14).
  * Regression tests: `tests/testthat/test-counts-review.R` (issue 15).
  * Performance: `binding = "auto"` (default) is off above 200 units; only restrictions
    active at the optimal LP vertex are refitted; joint LPs use GLPK when available; the
    rounded LP is skipped when rounding changes nothing (500 units with 11 totals: 218 s
    with binding before, 7.5 s by default and 73 s with `binding = TRUE` now).
* Code review, round 2: a connected component with rows but no target cells no longer crashes
  the MILP certificate (R2-1; any exactly checked feasible table is optimal there); the
  integer objective coefficients are capped at 1e7 after scaling (above it: "unresolved");
  the certificate c'n >= v* + 1 is solved in the shifted variables d = n - n* (cut c'd >= 1),
  so GLPK's relative tolerances stay far below its unit gap, and a returned table must
  improve exactly (strict progress); `Rglpk` and
  `slam` moved to Imports, so the integer check works out of the box.
* P3.5(a) conservative interval (proof review round 3, nit a): the studentised maximum is
  taken over the beta-dependent terms only; terms with zero standard error (e.g. the
  reference district's rows under M1(delta)) are held exactly and no longer enter the
  maximum as a constant 0. With no beta-dependent term the critical value is 0. This can
  narrow conservative intervals slightly where such terms exist; it changes no other output.
* Tests: `tests/testthat/test-counts.R` (contract C10: DD reproduction, brute-force
  enumeration per TU class, MILP rows and mixtures, T5, S5, RC11, N10c, RC9, 2c M2/M3
  closed forms, witnesses, the C8 wording, the embedding).

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
