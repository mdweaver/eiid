# PHASE6a code review: the count layer (M-COUNT, build step 1)

*Independent code review of `eiCounts()`, `witness()`, the cross-layer (R27) block and the `eiBounds()` embedding in package **eiPartialID 0.4.0**. Branch `claude/mcount-code` at 792d3f6. The reviewer did not write this code. The proofs (`PHASE6a_count_layer.md`, cited as 6a) are already proof-reviewed (`PHASE6a_proof_review.md`, cited as PR). This review checks the code against them and does not re-review the proofs. Nothing in the package or the theory documents was edited.*

*Scripts, in this folder:*
* *`phase6a_code_review_oracle.py`: an independent oracle in Python. It builds every row from its definition, using 2c Definition 1 with realised rates β = n/a. The affine coefficients come from evaluating each definition exactly with `fractions.Fraction`, not from Lemma C1's closed forms or the package's builders. It enumerates every integer table, and solves the LP and MILP with scipy/HiGHS.*
* *`phase6a_code_review_checks.R`: the R driver for the oracle, the named checks K1–K16 and the issue reproductions I1–I13.*
* *`phase6a_code_review_k10.R`: the 0.3.1 comparison.*

*How to rerun them is in §7.*

---

## 0. Verdict in brief

The core is sound. Five things were checked independently, and all were right except for the solver problem in item 1:
* The rows are built correctly from 6a §1.3 and Lemma C1.
* They are rounded inward in exact rationals.
* The TU-class router follows R19(b) (classes (i)–(iii) and transposes, decided on the whole row set per CR branch).
* The continuous LP is right.
* Every witness is right.

The checks behind this:
* **7,810 integer endpoints** from brute-force enumeration and 8,388 continuous endpoints from an independent LP, on 1,644 random instances. These cover classes (i), (ii), (iii) and their transposes, mixtures, crossing triples, pairwise orderings, pairwise δ_B, W2, generic linear rows, W1 signs with slack, CR with slack, and the odds dial.
* **15,620 witnesses**, checked exactly against the oracle's own rows: all reproduce the margins, satisfy every row, attain the reported value, and are whole tables when the status says integer-exact.
* T5, S5, RC11 and N10c reproduce.
* The thresholds are exact (RC9-type check, 150/150).
* The margins-only line equals the package's Duncan–Davis to 2.2e-16 for cells and one-column contrasts.
* Every model output under `target = "model"` is **identical** to 0.3.1, as are the print and the RNG stream, on 10 configurations.
* The C8 wording is verbatim, and "given the count rows" never prints.
* S_φ and t_E match Lemma 4.1, and the representation is chosen from the design only.
* For a single cell, the cross-layer model interval equals `eiBounds()`'s `proj_actual` to 1e-16.

**R CMD check** (rerun here): Status 1 NOTE (`ei.Datasets` not available). Tests: **1763 pass, 0 fail, 3 skip**, as the developer reported. In this container's POSIX locale there is also a WARNING ("cannot set locale en_US.UTF-8"). It comes from the environment, and a rerun under `C.UTF-8` gives 1 NOTE only. The POSIX run also exercises the ASCII fallback of the C8 wording.

**But five things must be fixed before release** (issues 1–4 and 16; details in §3):
1. **lpSolve's MILP sometimes returns a suboptimal integer endpoint.** The package then reports a range that is too narrow, labelled "integer-exact (MILP)". This happened on 4 of about 7,700 MILP endpoints. In each case the reported range excludes integer tables that satisfy every row. Scaling the objective to integers fixed all 4.
2. **`known` with unit labels not in sorted order gives the values to the wrong units, silently.** This turned the identity-row coverage of my R27 simulation from 1.000 to 0.000.
3. **Regression: `eiBounds(x, t, n)` with a unit of size 0 now stops with an error** under the default `target = "actual"`. It works in 0.3.1 and with `target = "model"`.
4. **Generic `linear` rows take double coefficients as exact rationals of the doubles.** For example, `coef = 1/3, lower = 1` gives the integer bound n ≥ 4 instead of 3, again labelled integer-exact.
5. **Pairwise δ_B and W2 rows hang the default call at realistic unit sizes** (issue 16). With N_i of a few hundred their coefficients cannot be scaled to integers, and lpSolve's MILP then does not return, even with 5 units. There is no time limit.

Two spec deviations in the R27 block should also be fixed (issues 7 and 8):
* a `where=` selection that uses the outcomes still gets a cross-layer interval printed as covering, against C2 and Lemma 4.1;
* a `total` given as a floor or an interval defaults to an identity row, but C8 lists stated floors among belief rows.

The rest are minor bugs, wording and test gaps.

| Contract | Verdict |
|---|---|
| C1 function and placement | **met** (issue 3 is a regression in one input path) |
| C2 targets; DD line | **met**, except issue 7 (design-measurability of `where`) |
| C3 rows | **met** (all row forms verified against definitions); label bug issue 13 |
| C4 integer modes, exact arithmetic | **not met as implemented**: issues 1, 4, 6, 16; the rest of C4 is met |
| C5 infeasibility | **met** |
| C6 witness | **met** |
| C7 print | **met**; wording issues 13, 14 |
| C8 cross-layer block | **met in structure and wording**; issues 7, 8, 9, 10 |
| C9 one vocabulary | **met** |
| C10 tests | **met as listed**; test gaps, issue 15 |

---

## 1. What was run

| Check | What | Result |
|---|---|---|
| R CMD check | `R CMD build` + `R CMD check --no-manual`, R 4.3.3, `_R_CHECK_FORCE_SUGGESTS_=false` | 1 NOTE (ei.Datasets); tests 1763 / 0 / 3 skip; vignette builds |
| Oracle, small | 420 + 1,200 random instances (`gen`, seeds 20261122 and 7): every scenario of §2.1 below, p = 2–4 units, N_i ≤ 7, R, C ∈ {2, 3}, brute force over up to 60,000 joint tables | 7,690 integer and 8,268 continuous endpoints compared; **4 integer mismatches (issue 1)**, 15 crashes (issue 5), 1 refusal (dW with Var 0, as specified); 15,380 witnesses all valid |
| Oracle, large | 24 instances, 60–150 units, 3×3, TU (totals over blocks + W1 sign + orderings), mixture, pairwise; HiGHS LP and MILP over the oracle's rows | 120 + 120 endpoints, 0 mismatches; 240 witnesses valid |
| K1–K9 | T5, S5, RC11 (with every integer of [2,4] attained), N10c with the rows in their natural non-integer form, RC4 routing, RC9 exact thresholds (150 populations with integer Y^ER), linear-row doubles, non-integer margins, DD to 1e-12 on 3 reference datasets | all PASS except K7 (issue 4) and K8b (issue 12) |
| K10 | `eiBounds()` 0.4.0 `target = "model"` vs 0.3.1 (installed from 9f575dd): M0, M0 conservative, M1, M2, L5, L6, M0 + model sign, M0 + CR, fpc = popoviciu, no seed | `$table` (without the two new `dd_exact_*` columns) and every other component `identical()`, printed output identical, RNG stream after the call identical |
| K11–K15 | R27: model interval vs `proj_actual`; S_φ, t_E and the representation; refusals; wording; statuses | PASS, with findings (issues 8, 9, 10) |
| K16 | R27 simulation, 400 replications (§5) | PASS |
| Edge cases, performance | §4 | issues 3, 5, 11, 12, 14 |

### 2.1 Oracle scenarios

Each instance has random common unit rows: a fractional box, and an audited known cell within the Fréchet bounds. On top of these it carries one of the following:

| Scenario | Rows | Expected route |
|---|---|---|
| tu_i | totals (fractional floors and ceilings, or exact values) on one cell over {all, {1,2}, {3,4}} ∪ {{2,3}} (two laminar families), plus a one-vs-rest ordering | (i) |
| tu_i_cross3 | totals on one cell over {1,2}, {2,3}, {1,3} (odd crossing cycle) | MILP |
| tu_ii, tu_ii_t | totals on cells of one column (one row), one laminar family per cell | (ii), (ii, transposed) |
| tu_iii, tu_iii_t | R = 2 (C = 2), totals on any cells, one laminar family per column (row) | (iii), (iii, transposed) |
| mix | totals on (1,1) and (2,2) over overlapping sets, R = C = 3 | MILP |
| pairwise, pair_dB, w2 | pairwise unit orderings, pairwise realised δ_B(r,s), W2 Δ_r (R = 3) | MILP |
| linear | a random generic row with coefficients {1, 2, −1, 0.5, 1.5} | MILP |
| sign_w1 | δ^r_B or δ^r_W (W1), either sign, slack 0 or 0.1 | (i) |
| cr_w1 | W1 CR, slack 0 or 0.05 | (i) per branch |
| odds | odds-ratio dial (2×2), several [ψ_L, ψ_U] | (i) |

The targets are a cell rate, a cell rate over a unit subset, a two-row gap, `ovr_gap`, a merged-row rate (R = 3), a multi-column rate (C = 3), `share`, and a gap across unit sets. The oracle builds them from §1.2 independently of the package.

---

## 2. Contract-by-contract verdicts

**C1. Function and placement: met.**
* `eiCounts()` is separate.
* `print.eiRxC` prints "Actual (realised) rates in these units", then the count block (the margins-only Duncan–Davis LP line first), then the cross-layer block, then the model table.
* `target = "model"` restores the old order. It still stores `$counts` and prints the count block after the model table only when count rows are given.
* K10: every model output is unchanged.
* One input path regressed (issue 3).

**C2. Targets: met, except issue 7.**
* `rate` (merged rows, several columns), `gap` (two-row and across unit sets), `ovr_gap`, `share` and `combo` are all exact against the oracle.
* The margins-only line is the LP. It uses the closed forms of Lemma 2.1(a)–(b) only when each unit's objective has one cell or two opposite-sign cells of one column, and the oracle confirms those closed forms.
* It reproduces the package's `dd_lower`/`dd_upper` to 2.2e-16 on `rxc_example_3x3.csv`, `msign_example_3x2.csv` and `msign_gap_3x2.csv`, and is narrower or equal for combination targets (K9).
* `share()` is count-layer only.
* But `design = TRUE` is set for any `where=`, including selections on the outcomes (issue 7).

**C3. Rows: met.**
* The oracle builds each row from its definition (2c Definition 1 with β = n/a; the one-vs-rest ordering as β_r − β_¬r ≥ −slack), not from Lemma C1. It agrees with the package on every endpoint in the `sign_w1`, `cr_w1`, `pair_dB`, `w2` and `pairwise` scenarios (issue 1 aside).
* So: W1 signs are totals on Y_rc(S) with the right thresholds and slack; pairwise and W2 rows have the right unequal coefficients; CR is the union of the right two systems; ordering and odds rows are the right boxes.
* "EXPERIMENTAL (2c M2; review batch B comment 1 open)" prints on sign and CR rows only.
* Print-label bug: issue 13.

**C4. Integer modes and exact arithmetic: not met as implemented (issues 1, 4, 6, 16). Everything else in C4 is met.**
* Met:
  * margins count as integers only when exactly integer (counts 3 + 1e-12 and share input both switch the check off; K8);
  * thresholds are bigq and rounded inward (ordering boxes, Y^NM, Y^ER with slack, the odds dial by an exact integer search), and RC9 passes 150/150;
  * the rows are rounded, never the endpoints (T5 and S5 pass);
  * the class test is on the whole row set per CR branch, with crossing-graph 2-colouring for class (i), and mixtures and crossing triples go to the MILP (K5);
  * `"exact"` refuses non-integer margins;
  * `"lp"` is labelled "continuous only".
* Not met:
  * "integer-exact (MILP)" is not guaranteed exact with lpSolve (issue 1);
  * pairwise δ_B and W2 rows hang at realistic sizes (issue 16);
  * the generic row takes inexact doubles (issue 4);
  * the integrality check in the TU path is vacuous (issue 6). The R19(b) theorem makes this harmless when the classifier and solver are right, but it removes the safety net the developer's report describes ("every integer witness is checked exactly, with MILP fallback").

**C5. Infeasibility: met.**
* "infeasible: the rows contradict the margins (continuous tables)" and "integer-infeasible: no whole-number table satisfies these rows" are separate. N10c is reported as the second.
* The CR branches print as nonempty or empty, and the union as one interval or two.

**C6. Witness: met.**
* 15,620 witness tables were checked exactly in Python against the oracle's rows (margins, every row, the reported value, integrality). All passed.
* `witness()` also re-checks them itself in exact rationals.
* The print says "a certificate of attainability, not an estimate".

**C7. print(): met.**
* It shows target and denominator, lower/upper/width, the four integer statuses, binding rows and their `object`, and the non-integer note verbatim.
* Wording issues: 13, 14.

**C8. Cross-layer block: met in structure and wording; issues 7–10.**
* Verified:
  * the α split is 0.9α for the model part and 0.1α for the allowance, the same as `eiBounds(fpc = "popoviciu")` (0.045/0.005 at α = .05); K11 matches `proj_actual` to 1.1e-16;
  * S_φ = ½(Σρ_i²)^{1/2} with ρ_i = Σ_r a_ir(max_c v − min_c v), and t_E = (½ log(2/α_E) Σρ_i²)^{1/2}. Both were recomputed independently for 7 target kinds: max error 2.8e-17 for S_φ, and 4.4e-16 for the Hoeffding-minus-Popoviciu width (K12);
  * for `ovr_gap`, the representation with the smaller S_φ is chosen, from the design only;
  * the Bhatia–Davis allowance is never used;
  * count rows and model rows never share an LP (they are separate code paths);
  * plug-in sets (`ci = "none"`), non-`eiRxC` objects and fits on other units are refused;
  * the emptiness of the margins-only intersection is kept;
  * the identity/belief wording is C8's text verbatim (Unicode, and an ASCII fallback in non-UTF-8 locales);
  * "given the count rows" never prints.
* Issues:
  * the outcome-selected `where` (issue 7);
  * the floor or interval `total` defaulting to identity (issue 8);
  * status labels with no §4.3 row (issue 9);
  * fresh Monte Carlo draws (issue 10).

**C9. One vocabulary: met.**
* `object = "model"` is refused in `eiCounts()`, and `object = c("actual","model")` is refused in both.
* In `eiBounds()`, `assoc`/`cr` with `object = "actual"` move to the count layer; the default stays `"model"`, so old calls are unchanged.
* The cross-layer print lists realised and model beliefs separately.

**C10. Tests: met as listed.**
* Every item of C10 and NS §3.1 has a test.
* The gaps are in issue 15: the 4 MILP failures were missed because the package's brute-force tests use designs too small to trigger them.

---

## 2b. Verdicts on the developer's deviations

| # | Deviation | Verdict |
|---|---|---|
| 1 | P3.8 not implemented; t_E is used as the allowance on an asymptotic interval | **Acceptable.** The union bound with an asymptotic model interval gives an asymptotic statement, and the status line correctly stays "P (asymptotic)". The finite-sample row of 6a §4.3 is not claimed. Suggest the print add "(the model part is asymptotic; the allowance holds in finite samples)" when `allowance = "hoeffding"`. |
| 2 | Generic `linear` row, always sent to the MILP | **Acceptable in principle** (R19(c)). Not as implemented: issue 4 (double coefficients) and issue 1 (MILP exactness). |
| 3 | Odds dial 2×2 only | **Consistent with 6a** (§1.3 defines it for 2×2 units). NS's "pairwise 2×C" is outside 6a's proved scope. |
| 4 | Binding rows by drop-one | **Acceptable.** For a single row it is the R2 face test ("removes the whole optimal face" ⇔ removing it moves the endpoint). Duplicated or jointly binding rows all print "none" (issue 14). Cost grows with rows × targets (§4). |
| 5 | Row roles: known/total identity, the rest belief, `role=` overrides | **Change needed for totals given as floors or intervals** (issue 8). The rest is right. `zero` as belief is conservative and acceptable. |
| 6 | Slack units follow eiBounds' convention | **Acceptable and verified.** The oracle implements the same scales independently (Σp X(1−X), Σp X_r X_s, μ Var_{w^r} + (1−μ)Var_{w^¬r}, Σ_g μ_g Var_{w^g}), and all slack instances agree. Note in the docs that the count layer computes the scale over the sign's own unit set S, whereas the model layer uses all units. |

---

## 3. Numbered issues

Severity: **bug** (wrong output or crash); **spec deviation** (output differs from what 6a requires); **wording**; **test gap**. Every reproduction is in `phase6a_code_review_checks.R` (mode `issues`, plus K7/K8b/K11b/K14b/K15 in modes `all`/`r27`).

### Issue 1 (bug, high): lpSolve's MILP returns suboptimal integer endpoints, reported as "integer-exact (MILP)"

* **What.** On 4 of about 7,700 MILP endpoints in the oracle runs (instances 305; 483, 692 and 797 of seed 7), the reported integer range is strictly narrower than the brute-force range. The status says integer-exact. The excluded values are attained by integer tables that satisfy every row exactly. All 4 cases involve fractional objective coefficients: rates over unit subsets and gaps across unit sets.
* **Reproduction (I1).**
  ```r
  d <- data.frame(r1 = c(2,1,1), r2 = c(3,1,2), c1 = c(1,2,1), c2 = c(3,0,1), c3 = c(1,0,1))
  eiCounts(rows = c("r1","r2"), cols = c("c1","c2","c3"), data = d,
           targets = list(gap("r1","r1","c2", where = 1:2, where2 = 3L)),
           restrictions = list(linear = list(cells = data.frame(unit = c(1L,3L), row = "r1",
                                              col = c("c3","c2"), coef = c(3, 0.5)), lower = 1.65)))
  ```
  The package reports the upper end 0 (integer-exact (MILP)). The table with n_{1,r1,·} = (0,1,1) and n_{3,r1,·} = (0,0,1) satisfies the rounded row 6n + n' ≥ 4 and every margin exactly, and has value 1/3. `integer = "exact"` gives the same 0. The other three cases are similar: pair δ_B [−1, 1/3] vs true [−1, 1]; a pairwise ordering {−2/3} vs [−2/3, −1/3]; linear [−1, 0] vs [−1, 1/3].
* **Cause.** Calling `lpSolve::lp(..., all.int = TRUE)` directly on the same system gives 0 with the objective as given, and 1/3 with the objective multiplied by 3 or 300. With the objective scaled by 1/min|c|, all 4 instances agree with brute force (the scaling patch is in the script). This looks like lp_solve's branch-and-bound pruning or its improved-solution threshold misbehaving with fractional objectives. Whatever the exact mechanism, the MILP's optimality is not certified.
* **Fix.**
  1. Before the MILP, scale the objective to coprime integer coefficients, from the exact bigq weights `tm$w` (lcm of denominators).
  2. Certify each MILP endpoint: after the optimum v*, solve one more MILP with the row c'n ≥ v* + 1 (scaled integer objective, maximisation) and require infeasibility. Or cross-check with a second solver (`highs`, which supports time limits; `Rglpk`), and fall back to the LP outer bound with an explicit status when they disagree.
  3. Add these 4 instances as regression tests, plus a randomised brute-force test on designs of the oracle's size (p ≤ 4, N_i ≤ 7) with subset and gap targets.

### Issue 2 (bug, high): `known` with units not in sorted order assigns the values to the wrong units

* **What.** `.ec_units()` returns `sort(unique(match(...)))` for labels and numbers. `.ec_parse_spec()` then pairs `value` with the sorted units by position (`rep_len(value, length(units))`). `known = list(unit = c("U2","U1"), value = c(4, 1))` sets U1 to 4 and U2 to 1. Duplicated labels also shift the recycling. The fit is silently wrong, and it is labelled an identity row.
* **Reproduction (I2):** the output is U1 = 4, U2 = 1. In the R27 simulation (§5), known cells in `sample(300, 30)` units gave identity-row coverage **0.000**. With sorted units it is 1.000 (§5).
* **Fix.** Keep the user's order for per-unit values: match without `sort`/`unique`, refuse duplicated units for `known`, and check `length(value) %in% c(1, length(units))`. Add a test with unsorted labels and unsorted row numbers.

### Issue 3 (bug, medium): regression, `eiBounds(x, t, n)` with a unit of size 0 fails by default

* **What.** `eiBounds(x = c(.3,.5,.4), t = c(.5,.4,.6), n = c(100,0,80), model = "L6")` works in 0.3.1, and in 0.4.0 with `target = "model"`. Under the default `target = "actual"` it stops with "the model fit and the count layer are not on the same units and margins". The cause: `.ec_problem` drops N_i = 0 units, but `st$X_all` keeps them on the x/t/n path, and `.ec_check_model()` is called with `stop()` inside `.rxc_count_block()`. (Formula input drops such units earlier and works.)
* **Reproduction:** I3.
* **Fix.** In `.rxc_count_block()`, align the units: compare only units with N_i > 0. Turn any cross-layer precondition failure into `cnt$cross_absent <- "<reason>"` instead of an error, so the count block can never make `eiBounds()` fail where 0.3.1 succeeded. Add a test.

### Issue 4 (bug, medium): `linear` rows take double coefficients as exact rationals of the doubles

* **What.** `coef = 1/3, lower = 1` (intended: n ≥ 3) is rounded exactly for the double 0.333…3 < 1/3, so the integer row becomes n ≥ 4. It is reported as "integer-exact (MILP)", and the true tables with n = 3 are excluded (K7). The same happens with `coef = 0.1, lower = 3 * 0.1`. Mixed denominators push the scaled coefficients past 2^40. The row is then left unrounded and solved in doubles, and the exact check of the MILP solution can then fail and stop the fit ("internal: ... please report").
* **Why it matters.** C4 requires exact rational thresholds, or a certified, printed tolerance. The N10c test in `test-counts.R` avoids the problem by pre-scaling the coefficients to integers, so the user-facing form is untested.
* **Fix.**
  * Accept exact coefficients: `gmp::bigq`, or a `num`/`den` pair of columns in `cells`.
  * For doubles, rationalise with a bounded denominator (continued fractions; for example, accept p/q with q ≤ 10^6 if |x − p/q| < 1e-12·max(1,|x|)), and print "coefficients rationalised to tolerance …" per C4.
  * Do the same for user-given fractional `lower`/`upper` of boxes and totals, and for `slack` and `psi`.

### Issue 5 (bug, low): W2 (or pairwise δ_B) row with all-zero coefficients crashes

* **What.** If X_r is constant over the sign's unit set, the W2 row's coefficients are all 0. `.ec_round_row` then divides by gcd(0, …) = 0 and gmp stops with "division by zero" (15 of the oracle's w2 instances).
* **Reproduction:** I5.
* **Fix.** Detect an all-zero row: it reads 0 ≥ −slack·scale, which is vacuous. Drop it with a printed note, or refuse it with the same message as δ^r_W with Var_p(X_r) = 0. The oracle treats it as vacuous.

### Issue 6 (bug, low; a safety net that does not work): the integrality check in the TU path is vacuous

* **What.** In `.ec_system`, `x <- round(s$x)` comes before `.ec_check(P, rows, x, TRUE)`, so `any(x != round(x))` is always FALSE. A fractional vertex (from a misclassified system, or solver noise above 0.5) would be rounded silently. If the rounded table happens to be feasible, its value, which can be suboptimal, would be reported as the integer endpoint. This review found no misclassification, so no wrong output came from it. But the report's "every integer witness is checked exactly with MILP fallback" does not hold for fractional vertices.
* **Fix.** Accept the rounded LP solution only if `max(abs(s$x - round(s$x))) <= 1e-6` **and** |c'round(x) − LP optimum| ≤ tolerance. Otherwise fall back to the MILP and set `fallback`.

### Issue 7 (spec deviation, medium): outcome-dependent `where` gets a cross-layer interval

* **What.** C2 and Lemma 4.1 require the target's coefficients to be functions of the design only, and otherwise "the count layer only". `rate("vet","yes", where = yes/(yes+no) > 0.5)` selects units on the outcome, yet gets `design = TRUE`. Its model interval [0.542, 0.854] is printed as "covering the actual value at the model interval's level" (I7). The same holds for `gap(where2=)`, `combo`, and vectors or labels computed from outcomes.
* **Fix.** Set `design = FALSE` when `all.vars(where)` intersects the outcome (column) variables of the formula or `cols`, or any variable other than the declared design variables. For vector, label or numeric selections, print "the unit set must be fixed from the design (6a Lemma 4.1); a set chosen using the outcomes voids the cross-layer interval". Optionally add an argument such as `design_sets = TRUE` that the user must confirm.

### Issue 8 (spec deviation, medium): a `total` given as a floor or interval defaults to an identity row (developer point c)

* **What.** `.ec_default_role` makes every `total` an identity row. 6a §4.1 and C8 define identity rows as rows whose right-hand side is the realised value (margins, a published total, an audited known cell). They list "a stated floor on a count" among belief rows. With `total = list(..., lower = 10)`, the print says "covers … at the model interval's level (identity rows only)" and "rejects (model + its rows) at level α" (K14b). R27(i)–(ii)'s exact statements do not hold for such a row unless it is a fact.
* **Fix.** Default to `identity` only for `value =` (equality) totals, and to `belief` for `lower`/`upper` totals. A user with a published bracket or a rounded published figure can set `role = "identity"`. Mention this in `?eiCounts`.

### Issue 9 (wording; spec gap, low): status labels with no row in 6a §4.3 (developer point b)

* **What.**
  * L5/L6 print "P (asymptotic; point-identified rung, one linear functional)". 6a §4.3 has no row for point-identified rungs.
  * Pooled L4 (and M1(δ) under the M1/M2 wording) print "P (asymptotic)". §4.3's scope for P3.5(a)/(b2)(ii) names "pooled M0; M1/M2 in R1".
  * The normal-plus-Lemma-4.1 construction is very likely valid (a delta-method interval for one linear functional at 0.9α, with the Lemma 4.1 allowance at 0.1α, by a union bound). But the label "P" asserts a reviewed register entry that does not exist.
* **Fix.** Until the theory owner adds §4.3 rows, print for these: "asymptotic (normal interval + Lemma 4.1 by a union bound; not in the 6a §4.3 table)". Raise the §4.3 rows for L4, L5, L6 and M1(δ) with the theory owner (proof item; not code).

### Issue 10 (wording or reproducibility, low): the cross-layer model interval redraws the simulation draws

* **What.** `.ec_model_interval()` calls `.rxc_zdraws()` afresh with the current RNG state, instead of reusing the fit's draws or seed. For `ci = "conservative"`, the cross-layer model interval for a single cell differs from the `ci_actual` printed in the model table of the same `eiBounds()` output, by 5.2e-4 (K11b). It also depends on the global RNG state when `eiCounts(model = fit)` is called. (The CLR projection interval matched to 1e-16 in K11.)
* **Fix.** Store the draws `Zs` (or the seed) in `fit$setup` and reuse them, so that for a single cell with `fpc = "popoviciu"` the two printed intervals are the same numbers.

### Issue 11 (bug, low): unit selections that silently do something else

* **What.**
  * Numeric unit indices that do not exist are dropped silently (`units = c(1, 99)` acts as unit 1; I11).
  * A label of a unit dropped for N_i = 0 gives "unknown unit label(s)", although the code comment says such units are "silently out".
* **Fix.** Error on out-of-range numbers. For dropped units, say "unit Z has N_i = 0 and was dropped", and decide one policy for both labels and numbers.

### Issue 12 (bug, low): `eiCounts(x =, t =)` vectors get labels "x", "" and "t", ""

* **What.** `cbind(x, 1 - x)` names the columns "x" and "". The documented X/notX and T/notT never apply, the complement row cannot be named, and targets print as " -> " (K8b). `eiBounds()` sets the names explicitly and is fine.
* **Fix.** `X <- cbind(x, 1 - x); colnames(X) <- c("X","notX")` (likewise for T), as in `.rxc_parse()`.

### Issue 13 (wording; a print bug, low): every `sign` row prints as "delta_W"

* **What.** `.ec_parse_restrictions()` builds the label from `c(ps, list(type = ps$type_a))`, but `ps` already has `type = "sign"`, so `s$type` in `.ec_label_of()` is "sign". Every δ_B row prints "realised delta_W …" (without "vs rest"). The computation uses `type_a` and is right (I13).
* **Fix.** `ps2 <- ps; ps2$type <- ps$type_a` before labelling, and add a print test for both types.

### Issue 14 (wording, low)

* **Binding rows.** Two identical binding rows, or two rows that bind only jointly, both print "none" while the endpoint is moved. Print "(drop-one: rows that bind only together are not shown)".
* **Labels.** A `known` row over 50 units prints all 50 values in its label. Targets over different unit subsets of the same size get identical labels ("r1 -> c1 [1 units]" twice), which `witness(target = "label")` cannot tell apart.
* **The Lindeberg note** (developer point d) prints even with `allowance = "hoeffding"`, where no Lindeberg condition is needed.
* **The NS interface sketch.** Its form `sign = ">=0"` fails with a bare R error ("missing value where TRUE/FALSE needed"). Accept ">=0"/"<=0", or give a clear message.

### Issue 16 (bug, high): pairwise δ_B and W2 rows hang at realistic unit sizes

* **What.** With unit sizes N_i of 200–1,500, the coefficients a_is/(N·N_i) of a pairwise δ_B row, and the (X_ir − x̄_g)/N of a W2 row, have denominators whose lcm exceeds 2^40 already at 5 units. `.ec_round_row` then leaves the row unrounded (`rounded = FALSE`), with coefficients near 1e-4 in doubles, and `lpSolve`'s all-integer MILP minimisation does not return. It took more than 120 s on 45 variables, and over 600 s at 100 and 500 units; with no time limit, `eiCounts()` hangs. At 3 units it took 1–2 s. The package's tests and my oracle use units with N_i ≤ 7, where the row is rounded to small integers, so neither saw this. Even when the solver does return, the unrounded double row combined with the exact check can stop the fit with "internal: … please report" (issue 4).
* **Reproduction.** `timeout 120 Rscript phase6a_code_review_checks.R hang pair` (or `hang W2`): 5 units, 3×3, N_i from 200 to 1,500, `sign = list(type = "dB", version = "pair", row = "g1", other = "g2", col = "o1")` (or `type = "dW", version = "W2"`), `integer = "auto"`.
* **Fix.**
  * Use a MILP solver with sound tolerances and a time limit (`highs`), on the row scaled to unit norm. Keep the exact rational check of the solution, and certify optimality as in issue 1.
  * On timeout, report the rounded-LP outer bound with the status "integer check not finished: LP outer bound (R19(c) item 1)".
  * Add a test at realistic N_i with `integer = "auto"` and a time budget.

### Issue 15 (test gap)

Missing tests:
1. MILP optimality on randomised designs at the oracle's scale, with subset and gap targets (issue 1);
2. `known` with unsorted or duplicated units (issue 2);
3. `eiBounds(x, t, n)` with N_i = 0 under the default target (issue 3);
4. `linear` rows with non-dyadic coefficients (issue 4);
5. W2 with constant X_r (issue 5);
6. that for a single cell the cross-layer model interval equals the model table's interval (issue 10);
7. that an outcome-based `where` gets no cross-layer interval (issue 7);
8. the role default of floor totals (issue 8);
9. the sign label (issue 13).

---

## 4. Edge cases and performance

| Case | Behaviour | Verdict |
|---|---|---|
| N_i = 0 units | dropped, counted in the header; a `where`/`known` naming one errors as "unknown label" | issue 11 |
| Zero margins (a group with no members) | `rate()` on it errors clearly; default targets (`targets = NULL`) stop with "row r3 has no members" instead of skipping that row | minor: skip such rows or report NA |
| R = 1 or C = 1 | works (rates are 1 or the shares) | ok |
| single unit | works; orderings and dials give boxes | ok |
| ties / duplicated rows | ranges right; binding "none" | issue 14 |
| infeasible rows (box above the cell's range, non-integer known value) | "infeasible (continuous)" or "integer-infeasible", reported separately | ok |
| CR with slack | branches and union as in 2c M2(d); oracle agrees on all 60 cr_w1 instances | ok |
| odds dial ψ = (0, ∞) | no restriction, as the convention says | ok |
| δ^r_W with Var_p(X_r) = 0 | refused with the Lemma C1(d) message | ok (as specified) |
| W2 with constant X_r | gmp crash | issue 5 |
| non-integer margins / shares / mismatch = "rescale" | integer check off, C7 note printed, `"exact"` refused | ok |

**Performance** (500 units, 3×3, R 4.3.3, one core):

| Case | Time |
|---|---|
| margins only, 9 cell targets | 0.2 s |
| 50 known cells, 9 targets, with binding | 6.5 s |
| + 10 county totals and a state total (class (i)), 9 targets, with binding (12 drop-one refits per target) | **218 s** |
| W1 CR, 9 targets, with binding | 24 s |
| mixture of totals (MILP), 2 targets, no binding | 3.7 s |
| pairwise ordering (MILP), 2 targets, no binding | 1.9 s |
| pair δ_B (MILP), 2 targets | **no result in 600 s at 100 and at 500 units; also at 5 units** (issue 16) |
| W2 δ_W (MILP), 2 targets | **no result in 200 s at 5 units** (issue 16) |
| `eiBounds(M0)`, 9 cells: `target = "model"` vs `"actual"` | 3.1 s vs 4.5 s |

The LP calls dominate: one dense `lpSolve` LP over all 4,500 variables per system per endpoint. The drop-one binding test multiplies that by the number of restrictions. Suggestions:
* compute binding rows from the LP duals or slacks first, and refit only the candidates;
* make `binding = FALSE` the default above, say, 200 units, with a note;
* in `.ec_optimise`, split coupled problems into connected components of the unit-coupling graph.

**MILP time limit (developer point e).** lpSolve exposes no time limit. The MILPs for totals and pairwise orderings were fast. Pairwise δ_B and W2 rows at realistic unit sizes hang the default call indefinitely, already at 5 units (issue 16). Recommend a solver with a time limit (`highs`), and on timeout the status "integer check not finished; LP outer bound shown (R19(c) item 1)".

---

## 5. R27 simulation (K16)

Design: 300 units with N_i from 150 to 600, kept fixed across replications. The unit rates are drawn each replication around constant means (0.7 and 0.3, sd 0.12), so Goodman (L6) holds. The model layer is L6, and the target is the realised rate(g1, yes) at α = .05. The scenarios:
* **(I) identity rows only:** the margins, and audited known cells in 30 units chosen from the design alone;
* **(B1) a belief floor:** (I) plus the floor Y_11 ≥ t, which is false in a sizeable share of replications;
* **(B2) a belief sign:** the realised δ_B ≥ 0 (true in essentially every replication).

| Scenario | Coverage of B_φ by C_φ ∩ CI_φ | Rejection rate (disjoint) | Reference |
|---|---|---|---|
| model interval alone | 1.000 | — | — |
| (I) identity rows | **1.000** (= the model interval's) | **0.000** ≤ α | R27(i), (ii) |
| (B1) belief floor Y_11 ≥ 34,982, false in 53.5% of replications | 0.465 | P(reject ∧ row true) = 0.000 ≤ α | (*): coverage ≥ 1 − α − P(false) = 0.415 ✓ |
| (B2) belief δ_B ≥ 0 | 1.000 | 0.000 | — |

400 replications, 300 units, α = .05 (K16 PASS).

What this shows:
* With identity rows only, the intersection covers whenever the model interval does, as R27(i) says: the count set always contains the truth, and the run would expose issue 2 if the units were unsorted.
* The rejection rate is at most α.
* With the false-in-part belief floor, coverage falls below 1 − α but stays above 1 − α − P(row false), and P(reject ∧ row true) ≤ α, as (*) says.

The L6 interval is conservative here, so the sizes are far below α. The theorem claims only "≤ α".

---

## 6. The developer's points

* **a. N17: integer endpoint 37/60 vs continuous 0.633.** Correct.
  * NS N17's [0.500, 0.633] is the continuous range. With whole people the odds ratio at 38 veteran yes votes is 38·28/(22·12) = 4.03 > 4, so the integer upper end is 37/60 = 0.617, reported as integer-exact (TU class), with `lp_upper` 0.633.
  * The package is right. NS N17 should say "continuous" (a note for the NS owner).
* **b. L5/L6 labelled "P (asymptotic; point-identified rung)".** Not supported by the §4.3 table as written. Relabel, and raise the §4.3 rows with the theory owner (issue 9).
* **c. A total given as an interval or floor defaults to identity.** It conflicts with 6a §4.1 and C8, which list stated floors among belief rows. Default to belief for `lower`/`upper` totals and to identity for `value =` totals (issue 8).
* **d. Note when max ρ²/Σρ² > 0.1.** Reasonable as a heuristic flag for Lemma 4.1(a)'s Lindeberg-type condition. Print the value, and suppress the note for `allowance = "hoeffding"` (issue 14). 0.1 is not a theorem threshold, so say "heuristic".
* **e. No MILP time limit.** Confirmed as a real failure, not only a risk: the default call with a pairwise δ_B or W2 row does not return at 5 units of realistic size (issue 16). Use a solver with a time limit, and degrade to the LP outer bound with an explicit status (§4). The same change can fix issue 1.

---

## 7. Reproducing this review

```sh
# R CMD check (from a copy of package/eiPartialID)
R CMD build eiPartialID && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual eiPartialID_0.4.0.tar.gz
# oracle (Python 3, numpy, scipy >= 1.9)
python3 docs/extension/reviews/phase6a_code_review_oracle.py gen inst.json 420 20261122
Rscript docs/extension/reviews/phase6a_code_review_checks.R oracle inst.json res.json
python3 docs/extension/reviews/phase6a_code_review_oracle.py check inst.json res.json
ORACLE_SCEN=mix,pairwise,pair_dB,w2,linear,tu_i_cross3 python3 .../phase6a_code_review_oracle.py gen inst_milp.json 1200 7   # then oracle / check as above
python3 .../phase6a_code_review_oracle.py genlarge inst_large.json 24 3                                                       # then oracle / check
# named checks and issue reproductions (PKG = path to package/eiPartialID)
Rscript docs/extension/reviews/phase6a_code_review_checks.R all
Rscript docs/extension/reviews/phase6a_code_review_checks.R r27
NREP=400 Rscript docs/extension/reviews/phase6a_code_review_checks.R sim
Rscript docs/extension/reviews/phase6a_code_review_checks.R issues
timeout 120 Rscript docs/extension/reviews/phase6a_code_review_checks.R hang pair   # issue 16
# 0.3.1 comparison: install 9f575dd's package into LIB031, then
Rscript docs/extension/reviews/phase6a_code_review_k10.R LIB031 OUT TESTDIR; Rscript ... default OUT TESTDIR; Rscript ... compare OUT
```

---

## Round 2 confirmation (claude/mcount-code at ef47982)

*This pass reviews the developer's fixes on `claude/mcount-code` at ef47982 (fast-forwarded into this branch). That commit merges this review and theory bbf8cf4. Rglpk 0.6-5.1 and slam 0.1-56 are installed from CRAN (the same Rglpk version as Ubuntu's `r-cran-rglpk`).*

### R2.1 What was rerun

| Check | Result |
|---|---|
| R CMD check (`C.UTF-8`) | Status: 1 NOTE (`ei.Datasets`); tests **1866 pass / 0 fail / 3 skip** |
| Oracle, small (420 instances), unpatched | 1,612 integer endpoints agree; **17 pairwise-ordering instances crash** with Rglpk "not conformable" (R2-1); 1 refusal (δ_W with Var 0, as specified) |
| Oracle, MILP suite (1,200 instances), unpatched | **136 crashes** (R2-1: pairwise, linear, crossing-triple rows); every non-crashing endpoint agrees |
| Oracle, large (24 instances), unpatched | 13 crashes (R2-1) |
| All three suites with a review-only workaround for R2-1 (`phase6a_code_review_patch_zero.R`) | small: 1,753 integer + 1,957 continuous endpoints, **0 mismatches**; MILP suite: 6,029 + 6,403, **0 mismatches**, including the 4 instances of issue 1; large: 118 of 120 agree, and the other 2 are reported as "integer status unresolved (certificate not reached): LP outer bound". Both are outer and correctly labelled. All witnesses are valid; the 7 flagged by the checker are the LP tables of unresolved endpoints, as documented. |
| K1–K16, I1–I13 | all PASS (K15 updated to the new §4.3 strings) |
| K10 vs 0.3.1 (10 configurations) | still identical: `$table`, every other component, the print and the RNG stream (new components: `zdraws`, `counts`, `target`) |
| `hang pair`, `hang W2` (5 units, N_i 200–1,500) | 0.1 s, integer-exact (MILP) |
| Realistic scale (3×3, N_i 200–1,500, 2 targets, no binding) | pair δ_B: 0.4 s at 100 units, 9.4 s at 500; W2: 11.1 s at 100 (one endpoint "unresolved (time limit): LP outer bound"), 52.5 s at 500; mixture of totals: 0.3 s / 1.6 s; pairwise orderings: 1.1 s / 4.7 s |
| `round2` mode | R2-1 reproduction; Rglpk-absent behaviour; the 200-unit switch; pruned vs full drop-one binding |

### R2.2 Verdict per issue

| Issue | Verdict | Evidence |
|---|---|---|
| 1 MILP optimality | **resolved**, apart from R2-1 | Integer objective from bigq. Second MILP certificate c'n ≥ v*+1, or ≤ v*−1 for minimisation (checked: for minimisation the row is −c'n ≥ −v*+1, which is right because c'n is an integer). 0 mismatches over 7,900 endpoints once R2-1 is worked around; I1 gives 1/3. |
| 2 known-cell order | **resolved** | I2: U1 = 1, U2 = 4. Duplicates and wrong lengths are refused. |
| 3 N_i = 0 regression | **resolved** | I3 ok under both targets |
| 4 double coefficients | **resolved** | K7: 1/3 means 1/3 (convergent with q ≤ 1e6, tolerance 1e-12, printed). "1/3" as a string is exact. |
| 5 all-zero W2 row | **resolved** | I5: dropped as vacuous, with a note |
| 6 vacuous integrality check | **resolved** | integral to 1e-6 *before* rounding, plus value = LP optimum, otherwise MILP |
| 7 outcome `where` | **resolved** | I7: `design = FALSE`, count layer only |
| 8 floor totals | **resolved** | K14b: lower/upper totals are belief rows, `value =` totals identity |
| 9 statuses | **resolved** | strings follow the new §4.3 rows of bbf8cf4 (K15) |
| 10 Monte Carlo draws | **resolved** | K11b: 2.2e-16 (the fit's `zdraws` are reused) |
| 11 unit selections | **resolved** | I11: out-of-range numbers and dropped units are clear errors |
| 12 x/t labels | **resolved** | K8b: X/notX, T/notT |
| 13 sign label | **resolved** | I13: "delta_B … (W1 vs rest …)" |
| 14 wording | **resolved** | binding note; known label "50 values"; subset labels "[U1]"/"[U2]"; no Lindeberg note under Hoeffding; `sign = ">=0"` accepted |
| 15 tests | **resolved as listed**, with one gap | `test-counts-review.R`, 17 tests. Its randomised MILP test uses a target that touches every unit, so it never reaches R2-1. |
| 16 pair/W2 hang | **resolved** | terminates at realistic sizes; on the time limit the status is "unresolved: LP outer bound", never integer-exact |
| nit (a) (proof review) | **correct; no reference output changes** | §R2.4 |

### R2.3 New issue

**R2-1 (bug, high; blocks release): a connected component with rows but a zero objective crashes the MILP certificate.**

* **What.** `.ec_milp()` splits the problem into connected components and skips a component only when it has *no rows and* a zero objective. A component with rows but a zero objective goes to `.ec_milp_sub()`. There `cc` is all zero, and the certificate row `cbind(..., sense * cc)[cc != 0, ]` has no triplets. So `max(trip[, 1])` is one less than `length(dir)`, and `Rglpk_solve_LP` stops with "Arguments 'mat', 'dir', and/or 'rhs' not conformable."
* **Why it is common.** It happens whenever a non-TU row lives in units that the target does not touch, which is the common case. Examples: pairwise orderings `order_unit(over = )` in every unit with a target over a county; a pair or W2 sign over a different unit set; a linear row elsewhere. The NS interface sketch has exactly this pattern.
* **How often.** 17 of 30 small pairwise instances, 136 of 1,200 MILP-suite instances, and 13 of 24 large ones.
* **Reproduction:** `Rscript phase6a_code_review_checks.R round2`, line `[R2-1]`. It uses 2 units, 3×2, `order_unit = list(row = "r1", over = "r2", col = "c1")` and `targets = list(rate("r1", "c1", where = "U1"))`.
* **Fix.** In `.ec_milp_sub()`, if `all(cc == 0)` return `optimal` with the first exactly checked feasible solution, since any feasible point is optimal (or skip the certificate). Alternatively, in `.ec_milp()`, solve zero-objective components only for feasibility. Add a test with a pairwise ordering in all units and a target on one unit.
* **Effect of a fix.** With a workaround that treats the empty certificate row as infeasible, all three oracle suites agree exactly (R2.1). So this one-line fix is the only thing left between the code and a clean pass.

### R2.4 The points the coordinator asked about

* **Certificate logic.** Sound.
  * It certifies only a solution that is integral to 1e-6, passes the exact rational check against the original rows and margins, and admits no strictly better integer table in the (outward-relaxed) system.
  * A better table found by the certificate MILP becomes the next incumbent, for up to 5 rounds; after that the status is "unresolved (certificate not reached)", which is sound and appeared twice in the large suite.
  * Remaining risk, low: the scaled integer objective may reach 1e9. With n in the thousands, GLPK's relative tolerances could in principle misjudge a cut c'n ≥ v*+1. A false "feasible" only costs a round, but a false "infeasible" (status 4) would certify a non-optimum. Not seen in 7,900 endpoints. Lowering the cap, or checking the certificate cut's slack, would close it.
* **GLPK status 1.** Handled correctly. Status 1 (undefined: presolve found the relaxation infeasible, or the time limit hit with no incumbent) counts as infeasibility only when lpSolve confirms the LP relaxation is infeasible; otherwise "unresolved (time limit)". Status 2 (time limit with an incumbent) is never certified. N10c is reported integer-infeasible (K4).
* **Outward relaxation.**
  * Rows with integer coefficients over 2^40 are normalised to unit max-norm and widened by 1e-6·max(1,|bound|).
  * The MILP therefore solves a superset of the integer set: an infeasible superset proves infeasibility, and a superset optimum is certified only if it passes the exact check on the original rows.
  * Correct. Optimum tables lying in the 1e-6 band come back "unresolved (solver tolerance)", never wrong.
* **Without Rglpk.**
  * Every MILP endpoint becomes "integer status unresolved (no MILP solver (install package Rglpk)): LP outer bound", for example [−1, 0.444] instead of [−1, 1/3]. It is outer and labelled.
  * TU-class systems are unaffected (T5 still integer-exact), and the joint LPs fall back to lpSolve.
  * Acceptable. But since the advertised integer check needs Rglpk, consider moving Rglpk and slam to Imports, or printing the install hint once at load.
* **The 200-unit binding switch.**
  * `binding = "auto"` computes binding rows at 200 units and prints "binding rows not computed (more than 200 units; set binding = TRUE)" at 201.
  * The pruning is correct: by complementary slackness, a restriction that is slack at the optimal LP vertex has a zero dual in every optimal dual solution, so removing it leaves the vertex optimal. The test uses a conservative 1e-7 activity tolerance, and pruning applies only to the LP or TU path, without CR or MILP.
  * Pruned and full drop-one agree on 114 random endpoint pairs.
* **Nit (a)**, the studentised maximum over terms with se > 1e-12·max(se). Correct:
  * A term with se = 0 does not depend on β, so its row holds exactly with any k ≥ 0.
  * k̂ ≥ 0 because the maximum includes |Z_target|, or else the (1−α)-quantile of a maximum of normals, which is ≥ 0 for α ≤ ½.
  * Moreover quantile_{1−α}(max(0, M)) = max(0, quantile_{1−α}(M)), so dropping exact zeros changes nothing at α ≤ ½. Only terms with 0 < se ≤ 1e-12·max(se), which are numerical zeros that the old code studentised as noise, can change k̂, and then only downward.
  * Reference outputs that change: **none found**. All 10 K10 configurations, plus M0, M1, M1(δ) with δ = 0.05, M2 and L4 under both `ci = "conservative"` and `"clr"` on `rxc_example_3x3.csv`, are numerically identical to 0.3.1 (max |diff| = 0). So the "can narrow slightly" in the developer's note applies only to inputs whose zero-SE rows come out as floating-point noise, not to the shipped reference data.

### R2.5 Release readiness for step 1

**No, not yet.** Everything from round 1 is resolved and verified, but R2-1 crashes the default call in a common configuration: a non-TU row in units the target does not touch. With the one-line fix in R2.3, and a test for it, I would expect a clean pass: the oracle suites already agree exactly with a workaround in place. If the fix lands as described, rerun `phase6a_code_review_checks.R round2` and the three oracle suites unpatched. If those pass, step 1 is ready to release.
