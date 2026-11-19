# Response to the referee: Phases 1, 2a, 2b and 3

**Reviewed material:** `PHASE1_foundations.md` and `PHASE2a_rxc.md` (first report); those two plus `PHASE2b_districts.md` and `PHASE3_inference.md` (revised report).
**Revised version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b) on branch `claude/wizardly-volta-lbh4bm` of https://github.com/mdweaver/eiid. All documents are under `docs/extension/`.
**R package:** `eiPartialID` 0.3.0 (`package/eiPartialID`).

**P1**, **2a**, **2b** and **P3** denote the four documents. Each corrected statement carries a note "*(Corrected after review …)*" saying what the earlier text claimed. Each correction of a mathematical claim also has a numerical check, in the `verify_*.py` scripts, that **fails on the earlier statement and passes on the new one**.

---

## 1. Summary

Thank you for two careful reports. We accepted every point.

Beyond the individual corrections, the revision changed four things:
* the inference framework now states where its randomness comes from, and drops clustered variances;
* the regime in which each inference result holds is tabulated;
* the Monte Carlo study was rerun at scale, with Monte Carlo standard errors and new designs;
* a later external review of a module built on these documents (Batch B) found three further corrections to the core documents, which we report here (§4.6).

---

## 2. What changed

### 2.1 Theory

| Document | Main changes |
|---|---|
| [P1](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE1_foundations.md) | <ul><li>New **§1.7**: the three identified sets, three hulls and three targets.</li><li>(1.5) derived from A5.</li><li>A3 without covariates.</li><li>A5 in matrix form, on the support of $P_X$ restricted to $A$.</li><li>A6 rewritten: model-based randomness, no clustering.</li><li>A9: a stacked rank condition.</li><li>A10: regime coverage.</li><li>D2 restriction count corrected.</li><li>D4 conditions on $(X,N)$.</li><li>Turnout marked open.</li></ul> |
| [2a](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2a_rxc.md) | <ul><li>§0.4 citation table.</li><li>S3 condition corrected (Lemma 3(c), Remark 5.1).</li><li>Compactness under $\operatorname{affdim}(A)=R-1$ (Lemma 3(b), Remark 5.2).</li><li>Theorem 1: facts F1–F7 and the "all directions at once" paragraph.</li><li>Remark 4.2: continuous versus integer rates.</li><li>§4.3: what each set needs from the data.</li><li>$P(X\in A)=1$ carried through §6.</li><li>New **Proposition 7′** (units outside $A$).</li><li>Prop. 7: the outer LP and value-function regularity.</li><li>Prop. 5: nonemptiness.</li><li>Cor. 6.1(c): exact condition.</li><li>Prop. 9: its domain, and moved earlier as motivation (Remark 1.0).</li></ul> |
| [2b](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2b_districts.md) | <ul><li>D0: per-cell versus joint sets.</li><li>New **Lemma D1(d)** (stacked M1 rank).</li><li>Theorem 1′: hypotheses H1–H4.</li><li>D3: pure-composition vertices.</li><li>**D4(b)** conditions on $(X,N)$ (Remark D4.1).</li><li>M1(0) labelled a sensitivity benchmark.</li><li>Remark 6.2 promoted to **Proposition D8** (DiD, exact moment condition).</li></ul> |
| [P3](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE3_inference.md) | <ul><li>P3-S1a.</li><li>New **§0.4** (where the randomness comes from) and **§0.5** (regime table).</li><li>Three targets.</li><li>Lemma P3.1: scope.</li><li>P3.2: non-nested hulls, error decomposition, Slater as interiority.</li><li>P3.4: for the sample-hull endpoint, with distinct projected gradients.</li><li>**P3.5** split into (a), (b1), (b2) and (b3).</li><li>New **P3.6(c)** (pooling test) and the test calibration.</li><li>P3.7 split.</li><li>New **Theorem P3.8** (finite-sample interval).</li><li>S1 labelled an outer bound.</li><li>**§8** rerun.</li><li>§10 updated.</li></ul> |

### 2.2 Code (`eiPartialID` 0.3.0)

* **Clustering.** The `cluster` argument and the clustered variance are removed.
* **Projection intervals** (`proj_model_*`, `proj_actual_*`) are reported beside the default CLR interval, which is labelled algorithmic.
  * Candidates whose gradients use a numerical null direction of $\hat V$ are dropped.
  * The variance bound in each column comes from an interval of the same version.
* **Pooling test** (P3.6(c)).
* **Calibration.** The pooling test and the M1-versus-M2 test are calibrated by a restricted wild bootstrap. Their χ² and F p-values are also returned.
* **Diagnostics.**
  * Rank errors name the groups that break the rank condition.
  * New notes for near-boundary fits, few units per group, and unbounded estimates.
* **Docs:** [`NEWS.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/NEWS.md); [`PHASE4_implementation.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE4_implementation.md). **Tests:** `test-review-response.R`, `test-batchB.R`.

### 2.3 Documentation

* **Guide:** the plain-language guide, the vignette and the README describe the new intervals, the tests and the experimental status of the sign module.
* **Benchmarks** in `benchmarks/`:
  * the NZ known-truth benchmark, rerun;
  * the size scripts for the two tests.

### 2.4 Verification scripts

| Script | Before → after | New checks |
|---|---|---|
| `verify_phase2a.py` | 22 → 32 | <ul><li>S3 singleton</li><li>unbounded $\Theta_O$ for a singleton $A$</li><li>bilinear value function</li><li>five versus six compositions</li><li>integer vertices</li><li>Prop. 5 emptiness</li><li>Prop. 7′</li><li>restricted support</li></ul> |
| `verify_phase2b.py` | 14 → 21 | <ul><li>D4(b) counterexample and the $D\perp(X,N)$ case</li><li>M1 restriction counts</li><li>Lemma D1(d)</li><li>D3 with pure-composition vertices</li></ul> |
| `verify_phase3.py` | 8 → 15 | <ul><li>distinct projected gradients</li><li>weak duality</li><li>projection coverage</li><li>Theorem P3.8</li><li>P3.6(c)</li><li>null-direction and residual-CLT examples</li></ul> |
| `phase3_montecarlo.py` | 200 → 1,000 replications per design | MC standard errors; riding-shock, Slater-failure and many-small-district designs; covariance-rank record |

---

## 3. Point-by-point: first report

Tags: **[Theory]**, **[Code]**, **[Docs]**, **[Check]**.

**2.1. The identification decomposition; A5 in matrix form.**
* **[Theory]** P1 (1.5) is now derived from A5, as in 2a Lemma 1(b), rather than from an unstated open-set condition.
* A5 is stated in matrix form: $x^\top Dx=0$ on the support implies $D=0$.
* The R = 3 point is in 2a Remark 1.1: six compositions off a common conic are needed, and three, four or five never suffice.
* **[Code]** Rank errors name the offending groups and give this hint.
* **[Check]** five versus six compositions.

**2.2. Derivative notation in Lemma 2.**
* **[Theory]** The derivative is replaced by "the coefficient of $K_{c,rs}$ in the ambient entries", with the free coordinates of Prop. 7 for the constrained parameterisation (2a Lemma 2; P1 (1.8)).

**2.3. Theorem 1: make the proof self-contained.**
* **[Theory]** 2a §4 lists facts F1–F7: $L^1$, uniform integrability, Mazur, Dunford–Pettis, weak continuity of the conditional expectation, the measurable maximum, and separation.
* A paragraph explains why the support-function inequality must hold for all directions outside one null set.
* The cited theorem numbers are checked in a new table (2a §0.4). The CLR result we use is their **Theorem 4**, not Theorem 2.
* Textbook theorem numbers are marked "standard" and still need checking against the printed editions.

**2.4. S3 for $A\subseteq\operatorname{conv}\operatorname{supp}P_X$ (error).**
* **[Theory]** Lemma 3(c) now requires $A\subseteq\operatorname{conv}(\operatorname{supp}P_X\cap A)$; the default $A=\operatorname{conv}\operatorname{supp}P_X$ satisfies it.
* Your singleton counterexample is Remark 5.1.
* The later review sharpened this to the support of the measure restricted to $A$ (§4.6, item 1).
* **[Check]**

**2.5. Carry $P(X\in A)=1$ through the propositions.**
* **[Theory]** Props 4–7 state it.
* New **Proposition 7′** bounds units outside $A$: the LP part for kept units, plus Duncan–Davis for the rest, with attainment.
* Remark 2.1 links to it.
* **[Code]** `range = "trim"` computes this bound.
* **[Check]**

**2.6. Joint value-function regularity (error).**
* **[Theory]** 2a Prop. 7 now makes three separate statements: concave and piecewise linear in $Q$ for fixed $\mu$; convex in the objective for fixed $Q$; neither jointly. Your example $V=\mu Q$ is included.
* P3 Lemma P3.1 states its fixed-design scope.
* **[Check]** midpoint test.

**2.7. Compactness and dimension.**
* **[Theory]** The condition is $\operatorname{affdim}(A)=R-1$ (2a Lemma 3(b), Theorem 1). Your singleton example is Remark 5.2.
* Remark 5.2 also shows that A5 implies full dimension, so Theorem 1 under its own assumptions was unaffected.
* $\Theta_O$ is called a polyhedron until boundedness is shown.
* **[Code]** A note warns when an estimate is unbounded.
* **[Check]**

**3.1. Covariates in Phase 1.**
* **[Theory]** A3 is now $E[\beta\mid X,N]=E[\beta\mid X]$ (with districts under M1/M2), and $Z$ is removed from the conditioning set.
* Covariate adjustment is labelled a planned extension (Phase 2c) that needs stronger assumptions. It is not claimed.

**3.2. Keep $\Theta_I$, $\Theta_X$ and $\Theta_O$ distinct.**
* **[Theory]** New P1 §1.7 tabulates each set with the information it uses and what it is sharp relative to.
* "Sharp" for $\Theta_O$ always says "for the regression function".
* Prop. 3 is called a sufficient condition.
* Props 4 and 7 are labelled outer-set results.
* New 2a §4.3 says what each set needs from the data.
* **[Docs]** The package prints "estimated outer identified set".

**3.3. Population support versus sample hull.**
* **[Theory]** P1 §1.7 separates the population hull, the hull of the empirical law and the sample hull used as an estimate.
* P3 §3 shows that coverage of the true parameter does not depend on the sample hull: the true parameter satisfies the constraints at the sample-hull vertices.

**3.4. Sharpness for continuous rates.**
* **[Theory]** 2a Remark 4.2 states that sharpness is relative to real-valued rates.
* It notes that unit transportation polytopes with integer margins have integer vertices.
* A randomisation argument for integer counts is given as a sketch only.
* **[Check]** integer vertices.

**4.2. Prop. 5 nonemptiness.**
* **[Theory]** Prop. 5 states that $\Theta_O\ne\emptyset$ iff all $Q_{c,rs}\ge0$.
* **[Check]**

**4.4. Corollary 6.1(c).**
* **[Theory]** The condition is now exact: strict tightening requires the common-sign condition to fail, so some extrapolated pure-poll rate is negative.

**4.5. Proposition 7 is an outer-set LP.**
* **[Theory]** Stated in Prop. 7.

**5. Proposition 9.**
* **[Theory]** Its domain is the support of $P_X$.
* It is also used earlier as motivation (2a Remark 1.0).

---

## 4. Point-by-point: revised report

### 4.1 Phase 2b

**D0: per-cell versus joint.**
* **[Theory]** The joint set is given separately; it is not the product of the cell intervals.

**D1: rank conditions.**
* **[Theory]** New **Lemma D1(d)**: a stacked rank condition for M1. It is sufficient to have A5 in one district plus $R$ affinely independent compositions in each other district; distinct but collinear compositions are not enough.
* Remark D1.1 separates identification from precision.
* **[Code]** Rank errors name the groups.
* **[Check]**

**Theorem 1′: state the conditions.**
* **[Theory]** Hypotheses H1–H4 (measurability, boundedness, affine and continuous in the parameter, unchanged selection set).

**The ladder: M1(0).**
* **[Theory]** Labelled a sensitivity benchmark rather than a substantive rung.
* **[Docs]** `ladder()` help says the same.

**D3: pure-composition districts.**
* **[Theory]** The interior condition is stated, and the limiting constraints for $t\in\{0,1\}$ are derived.
* **[Check]** LP agreement in 195 instances.

**D4(b): conditioning error.**
* **[Theory]** D4(b) now conditions on $(X,N)$: M0 holds iff $E[Z_D\mid X,N]$ is affine in $X$ and free of $N$, so $D\perp(X,N)$ suffices and $D\perp X$ does not (Remark D4.1).
* "Exactly when" is qualified for the outer bounds.
* Conditional on the observed district labels, P3.6(c) shows that only zero shifts make M0 correct, and gives a test.
* **[Check]** counterexample.

**DiD (Remark 6.2).**
* **[Theory]** Promoted to **Proposition D8**. The change is point-identified whatever the identified set for $K$ iff $\mu^{[\text{post}]}_{rs}/\mu^{[\text{post}]}_r=\mu^{[\text{pre}]}_{rs}/\mu^{[\text{pre}]}_r$ for $s\ne r$.
* Equality of the size-and-composition moments across the periods is sufficient; similar share distributions are not.

### 4.2 Phase 3: structure and consistency

**Lemma P3.1 scope.**
* **[Theory]** Stated for a fixed constraint system and objective.

**Prop. P3.2.**
* **[Theory]** Convergence along nested sampling, plus an argument for non-nested Hausdorff-convergent hulls.
* Slater is a strict-interiority condition, which coverage does not need.
* The error decomposition separates root-$n$ estimation error from a one-sided hull error of unknown rate.
* **[Code]** A note flags near-boundary fits.

**Theorem P3.4.**
* **[Theory]** Stated for the sample-hull endpoint.
* The non-Gaussian limit requires active gradients that differ in a direction where the estimator varies.
* Use $C-1$ columns (Remark P3.4.1).
* **[Check]**

### 4.3 Phase 3: confidence intervals

**P3.5(a) and dependence within districts (error).**
* We accepted that the residual allowance assumed independence that the sampling assumption did not grant. Our fix differs from both of your options; P3 §0.4 sets out the argument.
* **[Theory] Framework.**
  * Following Abadie, Athey, Imbens and Wooldridge (2020, 2023), clustering is justified by clustered sampling or clustered assignment. EI has neither: every precinct and riding is observed, and nothing is assigned.
  * The uncertainty is therefore model-based and conditional on the design, since the contextual-effects assumption concerns $E[\beta\mid V]$.
  * The targets $B$ and $\tilde B$ refer to this electorate, so a riding shock is a fixed feature of it. Under M1/M2 it is a district effect. Under M0 it is a misspecification, which the new pooling test (P3.6(c)) detects.
* **[Theory] Assumption.**
  * P3-S1a assumes independence given the design and the district effects.
  * Spillovers within a riding are a stated limitation.
  * A superpopulation target over re-drawn ridings is not pursued.
* **[Code]** `cluster` is removed.

**P3.5(b): CLR.**
* **[Theory]** Split into four parts:
  * (b1): CLR with the full dual-vertex set, a cited theorem (CLR Thm 4, Condition P).
  * (b2): proved here by LP weak duality. Any candidate set works with a fixed-superset or projection critical value.
  * (b3): the simulated candidate collection, a **conjecture**, with a proved reduction to (b1).
* **[Code]** CLR stays the default, labelled algorithmic, with the projection interval beside it.
* **[Check]** weak duality; projection coverage.

**Asymptotic regimes.**
* **[Theory]** New P3 §0.5 tabulates which result holds where:
  * M0 in every regime;
  * M1/M2 with a fixed number of districts;
  * M1/M2 with growing districts: conjecture;
  * many small districts: CLR coverage a simulation finding.
* New **Theorem P3.8**, a finite-sample interval, is valid in every regime. In every design we tried it is as wide as the logical range, and we say so.
* **[Check]**

**Three targets.**
* **[Theory]** P3 §0.2 separates $B^*$, $B$ and $\tilde B$, the last being the finite-design conditional mean.
* Only $\tilde B$ is affine in $K$.
* Inference for $B^*$ is not pursued.

**P3.6 tests.**
* **[Theory]** P1's restriction count is corrected to $(D^\star-1)(C-1)R(R-1)/2$.
* The χ² calibration is stated for fixed dimension.
* New P3.6(c) is the pooling test.
* **[Code]** Both tests are calibrated by a restricted wild bootstrap, which is a simulation-supported calibration, not a proof. Rejection at nominal 5% in design S0:
  * pooling test, 38 restrictions: 0.043, against 0.078 for χ²;
  * M1-versus-M2 test, 19 restrictions: 0.054, against 0.057 for χ².

**P3.7: many districts.**
* **[Theory]** Split into (a), "dropping constraints enlarges the set" (a proposition), and (b), a conditional-validity statement.
* "CLR restores coverage" is now a simulation finding.

**Sensitivity S1.**
* **[Theory]** Labelled an outer bound.

### 4.4 Cross-phase consistency

* **Turnout.** Marked open in P1 A1 and P3 §10, since it changes both the unit sets and the regression identity.
* **P1 A6 and A10.** Aligned with P3-S1a and with the regime table.

### 4.5 Monte Carlo (P3 §8)

* **Scale.** 1,000 replications per design, with Monte Carlo standard errors.
* **Designs:** riding shocks; a design where Slater fails; 30, 93 and 200 districts with 10 polls each; projection intervals throughout; the covariance rank recorded.
* **Findings:**
  * No proved interval under-covered beyond Monte Carlo error.
  * CLR covered at least 0.994, except where Slater fails (0.944, at the nominal level).
  * With many small districts, only CLR is informative.
  * With a non-affine riding shock, M0 fails, M1 covers, and the pooling test rejects every time.
* **NZ benchmark:** M0 and M1 by electorate are both rejected, and the M1 CLR interval missed in 5 of 63 non-empty cells. The package release notes report this.

### 4.6 Further corrections to these documents from the later review (Batch B)

1. **Support of the restricted measure.**
   * Almost-sure equality on $A$ extends only to the support of $P_X$ restricted to $A$, which can be smaller than $\operatorname{supp}P_X\cap A$.
   * **[Theory]** New notation $\operatorname{supp}_AP_X$ is used in 2a §0.1, Lemma 1, Lemma 3(c), Remarks 5.1–5.2, Props 7′ and 9, P1 A5 and 2b Lemma D1.
   * This also removes a tension between Lemma 3(c) and Remark 5.1.
   * **[Check]**
2. **Projection under a singular covariance.**
   * P3.5(b2)(ii) now requires a positive definite $\Sigma$, or a reduced full-rank parameterisation with exactly zero null-direction error, with a new proof.
   * **[Code]** The package drops candidates that use a null direction.
   * In the Monte Carlo, $\hat V$ had full rank in every replication of every design.
   * **[Check]**
3. **Residual CLT.**
   * P3.5(a) step 5 is completed by a Chebyshev/Lindeberg subsequence argument. The allowance is unchanged.
   * **[Check]**

---

## 5. What we now regard as established

| Result | Your assessment | Now |
|---|---|---|
| R×C identification decomposition; $\dim\mathcal K$ | Established | Unchanged; A5 in matrix form on $\operatorname{supp}_AP_X$ |
| Conditional-selection characterisation of $\Theta_I$ | Established, subject to conditions | Unchanged; proof expanded (F1–F7) |
| Compactness | Required a precise condition | $\operatorname{affdim}(A)=R-1$, implied by A5 |
| Joint value-function concavity | Incorrect | Corrected |
| District decomposition; baseline shifts; overidentifying restrictions | Established | Unchanged; restriction count corrected in P1; stacked rank for M1 |
| D4(b) | Incorrect conditioning | Corrected |
| DiD under M1 | Worth elevating | Proposition D8, with the exact moment condition |
| Sample-hull conservatism | Established | Unchanged; coverage shown not to depend on the hull |
| Root-$n$ plug-in under Slater; directional delta method | Established for fixed design | Stated for the sample-hull endpoint; Slater as interiority |
| Conservative CI for $\tilde B$, and for $B$ with independent residuals | Established | Unchanged; residual step completed |
| CI for $B$ with within-district dependence | Not established | Withdrawn; no clustering (P3 §0.4) |
| CLR with the full candidate set | Appropriate, subject to conditions | Cited theorem (CLR Thm 4) |
| Simulated CLR candidates | Outstanding | Conjecture P3.5(b3); the package labels it algorithmic |
| Projection interval | — | Proved (P3.5(b2)), with nonsingularity |
| Inference in R2/R3 | Requires further theory | Regime table; finite-sample interval valid but uninformative; CLR in R2 a simulation finding |
| Bounded nonlinear sensitivity | Established as an outer bound | Labelled so |
| Bounded-turnover sensitivity | Open | Open |

## 6. Not done, by decision

* **Covariate adjustment** (your 3.1, option B). It is planned as Phase 2c at low priority, because it needs stronger assumptions.
* **Bounded turnover.** It remains open.
