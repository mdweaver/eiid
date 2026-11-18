# Response to the referee: Batch B, partial identification and inference

**Reviewed version:** commit [`8b98589`](https://github.com/mdweaver/eiid/tree/8b98589653db71849fc7cdca500312754c26e95d).
**Revised version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b) on branch `claude/wizardly-volta-lbh4bm`.
**All changes since the reviewed version:** [compare `8b98589…e170be9`](https://github.com/mdweaver/eiid/compare/8b98589653db71849fc7cdca500312754c26e95d...e170be9161b9bbc250a17a27ffd39bc523c1c88b).

Abbreviations follow the report: **2c** is `docs/extension/PHASE2c_monotone.md`, **5a** is `PHASE5a_sign_module.md`, **2a** is `PHASE2a_rxc.md`, **P3** is `PHASE3_inference.md`, **1b** is `PHASE1b_monotone_did.md`. Line numbers in quoted comments refer to the reviewed version.

Each corrected statement in the documents carries a note "*(Corrected after review, Batch B, Nov 2026 …)*" saying what the earlier text claimed. Each correction also has a numerical check that **fails on the earlier statement and passes on the new one**.

---

## 1. Summary

Thank you for a careful, constructive and unusually concrete report.

* **Reproduction.** Your appendix script runs unchanged and reproduces every result, including both coverage bounds (0.933525 and 0.940092). It is saved as `reviews/referee_checks_B_partial_identification.py`.
* **Recommendation.** We agree with major revision, above all for 2c.
* **Disputes.** We accepted every point. On two points we argue only about severity (§5, items 1 and 3), and we made the corrections anyway.

**The main changes, in one paragraph each.**

* **Scope of the "exact" bounds.**
  * N1 (Prop. M2) is now stated as sharp for **continuous** unit tables given the margins.
  * A new remark proves that in 2×2 with margins only, the integer gaps form a lattice inside the continuous interval.
  * Integer tables in general are a new open item (R15).
  * The difference-in-differences rung Δ4 inherits the same label.
* **District rungs.**
  * M5 is restated rung by rung. Part (c) holds for the pooled model and within each district, and the pooled district version has an exact between-district term, illustrated by your reversal example.
  * R3 and R10 are now stated for the pooled model. Their grouped versions are a new open item (R16).
  * The package marks the module **experimental** and says so at fit time under grouped rungs.
* **N2.**
  * N2 is now defined on the finite empirical support.
  * The sharp version is withdrawn; your $V$-conditional statement is recorded as **Conjecture N2.2′**, left open for now.
  * N2.3 is kept only as a structural identity; the local-derivative rows are **dropped** from the build plan.
  * Your binning repair is proved as **Lemma N2.4**.
  * The polynomial degree must be prespecified.
* **5a.**
  * R2's feasibility LP now includes the target intercept, and (e) speaks of zero gain for a target.
  * R3's refutation claim names its base system.
  * R10's status is your version table, and R10(e) is now a theorem with your analytic proof.
* **Core documents** (reported also to the original referee of Phases 1–3):
  * the support of the restricted measure (2a);
  * the projection interval under a singular covariance (P3.5(b2)(ii)), with a guard in the package;
  * the residual CLT (P3.5(a) step 5).

---

## 2. What changed

### 2.1 Theory

| Document | Changes |
|---|---|
| [2c](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2c_monotone.md) | <ul><li>**M2** relabelled as sharp for continuous unit tables. New **Remark M2.1**: your and the second referee's integer examples, plus a proved lattice statement for 2×2 with margins only. M2(d) is componentwise, with a breakdown frontier.</li><li>**M3**: hypotheses stated (continuous rates, no unit restrictions, no other sign rows); (d) covers $\gamma<1$, $\gamma=1$ and $ER=0$; the "Reading" names its base system.</li><li>**M5**: restated per rung with a stacked parameterisation; (c) for M0 and within districts; exact decomposition with the between-district term and your reversal example.</li><li>**Remark M5.1**: rate statement replaced.</li><li>**Definition 1** remarks: W2 is a pooled sum; Lemma M1 hypotheses split by part.</li><li>**N2**: support convention; N2.2(a) outer only; **Conjecture N2.2′**; N2.2(b) on the complete support only; N2.2(e) with a chart and weighted covariances.</li><li>**N2.3**: structural identity only, with your smooth example for (d) and the reverse direction.</li><li>New **Lemma N2.4** (person-weighted bins).</li><li>§5 ladder and §6 (fixed degree, separate notation for support size and degree) rewritten.</li></ul> |
| [5a](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md) | <ul><li>**R3**: scoped to pooled M0; "finite LP when $C$ is polyhedral"; the design-weight convention; (C2) noted as stronger than (a) needs. New **Remark R3.1** (full law versus regression, $C=[.2,.4]$) and **Remark R3.2** (LCE excludes CR, $\tilde D=-261/1925$). Common-witness example in §1.4. (e) "rejects the conjunction; does not attribute".</li><li>**R2**: step 4 uses $h^\top k\ge U-\varphi_0$; (b) duplicate-row example; (c) witness objective; (d) qualifications; (e) "zero incremental gain for this target"; §2.5 items 2–4 (joint endpoint inference).</li><li>**R10**: §3.1 points to P3's regime table; §3.2 status by interval version (your table); actual-gap residual envelope; (b) one event, no Slater; (c) point coverage, no branch selection, triangular arrays; (d) pointwise; **(e) Theorem** with your analytic proof.</li><li>§3.4 table: union lengths merged, a "does not test" list and the false-CR table.</li><li>§3.5 items 2 and 4.</li><li>§4: pointer to this review.</li><li>§5: release 0.3.0 table.</li></ul> |
| [2a](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2a_rxc.md) | New notation $\operatorname{supp}_AP_X$ (support of $P_X$ restricted to $A$), defined in §0.1 with your example; [supp] proved; Lemma 1, Remark 1.1, Lemma 3(c), Remarks 5.1–5.2, Props 7′ and 9 updated. |
| [P3](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE3_inference.md) | <ul><li>**P3.5(b2)(ii)**: requires $\Sigma$ positive definite (ii-a) or a reduced full-rank parameterisation (ii-b); new proof; your null-direction counterexample.</li><li>**P3.5(a) step 5**: Chebyshev/Lindeberg subsequence argument; [Lyap] tag corrected.</li><li>P3.6: calibration of both tests by a restricted wild bootstrap.</li><li>§10 updated.</li></ul> |
| 1b, P1, 2b | 1b: theorem-to-table mapping and pure-unit caveat. P1 A5 and 2b Lemma D1 use the restricted support. |
| [`PACKAGE_DESIGN.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PACKAGE_DESIGN.md) | <ul><li>Proof register: R2, R3, R7, R9 and R10 updated; new **R15** (integer tables) and **R16** (grouped rungs).</li><li>P-N scoped to M0; "exact" → "sharp for continuous tables".</li><li>Build step 8 uses the global N2 only.</li></ul> |

### 2.2 Code (`eiPartialID` 0.3.0)

* **Status.** The sign and reinforcement restrictions are **experimental**, and `print()` shows the flag.
* **Grouped rungs.** A fit-time note says that each sign is one pooled restriction across groups and that sharpness is proved for the pooled model only. The printed sharpness claim appears only for pooled fits.
* **Refutation messages.**
  * They name the base model.
  * Next to a joint rejection, they say whether the rung alone is rejected.
  * An empty plug-in is called a sample diagnostic, not a test.
* **`eiRedundancy()`.** "Refuted jointly" became "empty at the point estimate (see the emptiness test)", and "cuts the set, not this target" is documented as zero gain for this target.
* **Projection rank guard (your §5 item 2).** Candidates whose gradients use a numerical null direction of $\hat V$ are dropped, which is case (ii-b). In the CLR interval, a zero-variance candidate is exact only if its gradient is zero.
* **Variance bounds.** The Bhatia–Davis bound in each column now uses an interval of the same version: CLR, projection or conservative.
* **Tests.** The pooling test and the M1-versus-M2 test are calibrated by a restricted wild bootstrap. Size at nominal 5%: 0.043 and 0.054; χ² gave 0.078 and 0.057.
* **Face test (P-B18).** The code already included the target intercept; a regression test now covers it.
* **Test file:** [`tests/testthat/test-batchB.R`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/tests/testthat/test-batchB.R).
* **Changes and known limits:** [`NEWS.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/NEWS.md).

### 2.3 Documentation

* [`PLAIN_LANGUAGE_GUIDE.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PLAIN_LANGUAGE_GUIDE.md) §9 and the [vignette](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/vignettes/eiPartialID_rxc.Rmd) section on signs and reinforcement:
  * refutation stated with its base model;
  * W2 described as a pooled sum;
  * quoted simulation figures name their method;
  * the default reinforcement interval described as algorithmic;
  * a "reading the output" section.
* [`PHASE4_implementation.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE4_implementation.md): updated for 0.3.0.

### 2.4 Verification scripts

| Script | Before → after | New checks (each fails on the earlier statement) |
|---|---|---|
| `verify_phase2c.py` | 12 → 26 checks | <ul><li>C4 rebuilt with valid transition matrices</li><li>C5(c′), C7, C7a, C7b, C8, C8b, C9, C10, C10b, C11, C12, C12b, C13, C14: integer examples, $\gamma=1$, district reversal, grid and size examples, derivative non-identification, binning, smooth (d), W1/W2 reversal</li></ul> |
| `verify_R2.py` | 4 → 6 checks (+ T4) | <ul><li>T1c: target with intercept</li><li>T4b: target-specific redundancy</li><li>T4 records solver status (147/494/1153 of 1794, plus 6 unclassified)</li></ul> |
| `verify_R3.py` | 2 → 4 checks | <ul><li>S6: Remark R3.1</li><li>S7: Remark R3.2</li><li>S2 now separates empty, point and interval branches</li></ul> |
| `verify_R10.py` | MC → MC + 3 checks + V12 | <ul><li>V9–V11: the event inclusions</li><li>union lengths merged (0.357 / 0.242 / 0.480)</li><li>V12: the false-CR table</li></ul> |
| `verify_phase2a.py`, `verify_phase3.py` | 31 → 32, 13 → 15 | the restricted-support counterexample; the null-direction and residual-CLT examples |

---

## 3. Revised verdict table

| Result | Your verdict | Now |
|---|---|---|
| M1 | Established | Unchanged; hypotheses split by part (2c §1) |
| M2 | Established for continuous rates; integer claims incorrect | **Proposition for continuous unit tables**; Remark M2.1 (integer lattice in 2×2 with margins, proved); R15 open |
| M3 | Incorrect as a whole | **Proposition** with stated hypotheses; (d) covers $\gamma<1$, $\gamma=1$, $ER=0$ |
| M4 | Established | Unchanged |
| M5 | Incorrect across the advertised rungs | (a) every rung; (b) district-specific ordering; (c) M0 and within each district, with the exact between-district term |
| N2.1 | Established | Unchanged |
| N2.2 | Established under conditions | (a) outer only; **Conjecture N2.2′** for the sharp version; (b) complete finite support; (e) chart and weighted covariances |
| N2.3 | Not established as an identified restriction | Structural identity only; local rows dropped from the build |
| (new) Lemma N2.4 | — | Proved (your repair): person-weighted bins give an outer relaxation |
| R3(a)–(e) | Established, with qualifications | Pooled model; finite LP for polyhedral $C$; Remarks R3.1–R3.2; (e) reworded |
| R2(a)–(e) | Established under conditions; (e) wording | Step 4 with the intercept; qualifications in (b)–(d); (e) target-specific |
| R10(a)–(d) | Established under the stated assumptions | Status by version: (a) theorem; (b1) cited under the augmented Condition P; (b2) theorem; (b3) algorithm |
| R10(e) | Established; proof completion | **Theorem**, using your analytic proof |

---

## 4. Point-by-point

Tags show where each change lives: **[Theory]**, **[Code]**, **[Docs]**, **[Check]**.

### Major comment 1: a continuous LP and an integer electorate can reach opposite refutation decisions

> *M2's step 106 gives real cell masses, not integer counts; your three-unit example shows that CR can be infeasible for every integer allocation. M3(d) omits the pure-unit boundary case. N1 should be labelled sharp for continuous tables; a scalar breakdown value needs a path.*

**Response.** Accepted in full.
* **[Theory]**
  * M2 is now "sharp for **continuous** unit tables given the margins" (2c §2), with a note recording the earlier claim.
  * New **Remark M2.1** gives your three-unit example and the second referee's two-precinct example. It also proves that in 2×2 with margins only, the integer gap set is the Duncan–Davis interval intersected with a lattice of step $1/(N\mu(1-\mu))$.
  * M2(d) now says nesting is componentwise in the slacks: a vector of slacks gives a breakdown frontier, and a scalar breakdown value needs a stated path.
  * M3 states its hypotheses (continuous rates, no unit-level restrictions, no other sign rows). M3(d) now distinguishes $\gamma<1$ (exactly one branch when $ER\ne0$), $\gamma=1$ (both branches give the same identified gap) and $ER=0$ (gap zero; the singleton is the gap's image).
  * The difference-in-differences rung Δ4 inherits the continuous label (`PACKAGE_DESIGN.md`).
  * Integer tables in general are a new proof item, **R15**.
* **[Code]** N1 is not implemented. The package's wording about feasibility now names the continuous model-free problem.
* **[Check]** C7, C7a, C7b (integer examples; lattice), C8 and C8b ($\gamma=1$, $ER=0$).
* **Severity (P-B1, P-B3).** We accept both corrections, and note two mitigating facts:
  * In the case the package would use first (2×2, margins only), integer and continuous answers differ by less than one lattice step, about $6\times10^{-6}$ at the scale of our Wisconsin application.
  * $\gamma=1$ requires every unit to be pure, in which case the gap is point-identified.

### Major comment 2: district rungs preserve affine rows but can reverse the sign predicted by M5(c)

> *M5(a) holds for every rung via the stacked parameterisation; (b) needs district-specific ordering; (c) fails under M1 because the district intercept is correlated with $X$ (six-unit example, pooled $\delta_W^*=-0.0493$). R3 and R10 extensions to districts must retain the district design.*

**Response.** Accepted.
* **[Theory]**
  * M5 is restated per rung with the stacked district parameterisation (2c §3): (a) holds for every rung; (b) requires ordering in each district; (c) is proved for M0 and within each district.
  * For pooled district rungs there is an exact law-of-total-covariance decomposition, including the between-district term $\sum_g\mu_g\operatorname{Cov}_{w^g}(X,\zeta_{g,D})$. Your example illustrates it: −0.04927 pooled, +0.00064 within each district.
  * R3 and R10 are now stated for the pooled model (5a §1, §3). Their district versions are a new proof item, **R16**.
  * The nesting chain P-N in `PACKAGE_DESIGN.md` is scoped to M0.
* **[Code]** The sign rows were already computed with each unit's own district coefficients. What changed is the labelling:
  * a fit-time note under M1, M1(δ), M2 and L4 says the module is experimental and that each sign is one pooled restriction;
  * `print()` claims regression sharpness only for pooled fits, and otherwise says "outer set; sharpness under grouped rungs not yet proved".
* **[Check]** C9 (your reversal), and a package test for the grouped-rung note.

### Major comment 3: N2 needs a support convention, conditioning on size, and an identified derivative

> *(i) "Support points" must be the finite empirical design; on a grid of a larger support, extendibility can fail. (ii) The proposed sharp version conditions on $X$ only, but the model conditions on $(X,N)$. (iii) N2.3's derivative is not identified on isolated support points.*

**Response.** Accepted.
* **[Theory]**
  * (i) A support convention defines N2 on the finite empirical design law (2c §4). N2.2(b) holds on the complete finite support only. Your grid example is Remark N2.2.2.
  * (ii) N2.2(a) now claims only an outer set. Your $V$-conditional statement is recorded as **Conjecture N2.2′**, with the proof route through 2a Theorem 1. The size example is Remark N2.2.1. We have deliberately left the conjecture open for now, because a separate review of the rungs bears on it.
  * The same caution is recorded for the planned series bound M-SHARP: it is "tighter than the outer set", not sharp, when size is informative (`PACKAGE_DESIGN.md` R9).
  * (iii) N2.3 is now a structural calculus identity for $\tau_m$, with a statement that it is **not** an observable restriction on a finite support, and your $h(x)=\prod_j(x-x_j)$ example.
  * The local-derivative rows are **dropped** from build step 8.
  * N2.3(d) uses your smooth example ($\delta^*_W=-1/47250$) and states the reverse non-implication.
* **[Check]** C10, C10b (grid and size), C11 (derivative), C13 (smooth (d)).

### Major comment 4: binning and cross-validation in 2c §6 require different arguments

> *Ordinary bin averages do not satisfy the accounting identity; person-weighted bin rates do. A CV-selected degree invalidates fixed-degree inference; "J" is used for two things.*

**Response.** Accepted. We adopted your repair.
* **[Theory]**
  * New **Lemma N2.4** proves that person-weighted bin rates reproduce $E_bT$, preserve the row targets and inherit monotonicity, and so give an outer relaxation of the global N2 bounds (2c §4, credited).
  * §6 now requires a prespecified degree, or sample splitting or a selection-aware procedure.
  * Support size is $J_X$ and degree is $q$.
  * Build step 8 in `PACKAGE_DESIGN.md` cites the lemma.
* **[Check]** C12, C12b.

### Major comment 5: R3 is sharp for specified information; refutability changes with that information

> *(a) and (b) are distinct; "LP" needs a polyhedral $C$. R3(e) should say an empty set rejects the conjunction without attributing it. 5a line 309's "unrefutable belief" imports M3 into the wrong base model: LCE alone can fit while both CR branches are empty.*

**Response.** Accepted.
* **[Theory]**
  * R3 now says "a finite LP when $C$ is polyhedral" and records the design-weight convention. It notes that (C2) is stronger than (a) needs.
  * New **Remark R3.1** gives your full-law-versus-regression example ($C=[.2,.4]$: full-law set empty, regression gap set $[1/6,1/3]$).
  * New **Remark R3.2** gives your example where LCE fits but excludes CR ($\tilde D=-261/1925$).
  * R3(e) reads "an empty population set rejects the conjunction; emptiness alone does not attribute the failure". It includes the $ER=-2$ example of a belief refutable model-free, and distinguishes population emptiness, an empty plug-in, and a test rejection.
  * 5a §3.5 item 2 now names its base system.
  * Your common-witness strip example is in §1.4.
* **[Code] [Docs]**
  * Every "never refuted on their own" string was replaced: roxygen, module metadata, the emptiness note, the guide and the vignette. The new text says that a δ_B sign and CR are always feasible only in the model-free problem with continuous tables, and that under the rung the combination can be rejected.
  * Next to a joint rejection, the output now says whether the rung alone is rejected.
* **[Check]** S6, S7 in `verify_R3.py`. S7 fails on the old wording.

### Major comment 6: R2 certifies changes to a specified endpoint

> *Step 4's feasibility LP omits the target intercept ($h^\top k\ge U$ should be $\ge U-\varphi_0$). A positive multiplier at one optimum does not certify tightening (duplicate-row example). (e)'s "redundant" should be "zero incremental gain for this target".*

**Response.** Accepted.
* **[Theory]**
  * Step 4 now uses $h^\top k\ge U-\varphi_0$.
  * (b) includes your duplicate-row example.
  * (c) notes that the separating objective need not be a reported target.
  * (d) states the qualifications (positive variance, $\gamma<1$, non-constant target).
  * (e) reads "zero incremental gain for this target". Redundancy for the feasible set is decided only by the inclusion test of (c).
  * §2.5 item 4: inference on a width gain needs joint treatment of the base and augmented endpoints.
* **[Code]**
  * The package's face test already included the target intercept; a regression test (P-B18) now covers it.
  * The verdict "cuts the set, not this target" is documented as zero gain for this target.
* **[Check]** T1c (intercept) and T4b (target-specific redundancy) in `verify_R2.py`.

### Major comment 7: R10 inherits a coverage theorem only when the augmented procedure has one

> *A version table: conservative (a) supports (a)–(d); full-candidate CLR needs the augmented Condition P; projection controls both ends at once; simulated-candidate CLR remains conjectural. The regime caveat must stay separate from the candidate conjecture. Actual gaps need their own residual envelope. Analytic counterexample for (e).*

**Response.** Accepted. We adopted your table and your proof.
* **[Theory]**
  * 5a §3.2's status note is replaced by your version table: (a) theorem; (b1) cited, under the augmented Condition P; (b2)(i) and (b2)(ii) theorem; (b3) algorithm. Condition P's actual requirements are listed, with the zero-variance-candidate caveat.
  * §3.1 now points to P3 §0.5's regime table instead of "R1/R3".
  * R10(a) adds the actual-gap envelope $\tfrac12\{\sum_i(w_{i1}+w_{i0})^2\}^{1/2}$ and notes that the package's per-cell sum is also valid (Minkowski).
  * (b): one simultaneous event, no Slater.
  * (c): point coverage, no branch selection, and the triangular-array remark.
  * (d): pointwise only.
  * **R10(e) is now a theorem** with your fixed-design proof (0.933525; 0.940092 with 4,000 draws; far-edge 0.927), credited.
  * §3.5 item 4 replaces "a fortiori" with your direct ≤ 1/2 argument and counterexample.
* **[Code]**
  * The default CR interval (CLR branch union) is labelled **algorithmic** on its own line, with the projection branch union printed beside it.
  * A zero-variance CLR candidate is exact only if its gradient is zero.
* **[Check]** V9–V11 in `verify_R10.py`.

### Minor comments and the EGGHM mapping (report §3)

* **W1 versus W2.** The Lemma M1 hypotheses are split by part (invertibility only for $g$; positive variance for $ER^r$). The W2 sentence says it signs a sum. **[Theory]** 2c §1. **[Check]** C14 (W1/W2 sign reversal).
* **Residual order (2c line 160).** Replaced by a root-sum-of-squared-coefficients bound, valid under P3-S1a. **[Theory]** 2c Remark M5.1.
* **The ladder.** "Both signs" now means both weakly positive or both weakly negative. N2-local is removed. A simplex chart is specified for $R\ge3$, and C5(c) uses the weighted covariances. **[Theory]** 2c §5, N2.2(e). **[Check]** C5(c′).
* **Identification versus sample labels.** R3(e) is a population statement. Raw emptiness is a sample diagnostic. Width-gain inference must be joint.
  * **[Theory]** 5a §1, §2.5.
  * **[Code]** `eiRedundancy()` says "empty at the point estimate (see the emptiness test)".
  * **[Docs]** guide §9.
* **EGGHM mapping.**
  * EGGHM's mixed-unit condition ($\gamma<1$) is cited where branch uniqueness is used (2c §0, M3).
  * Integer interpolation is not licensed by EGGHM's sharpness proof (2c M2).
  * Table 1 is attributed to Theorem 1, and Tables 2–3 to Theorems 2–3 (1b §A1).
  * The local results are cited for structural, not identified, derivatives (2c N2.3).
  * W2, pairwise and their CR are marked as extensions (2c Definition 1 remarks).
  * **[Theory]**

### Verification scripts (report §4)

* **`verify_phase2c.py` C4.** Rebuilt with valid transition matrices.
* **`verify_R3.py` S2.** It separates empty, point and interval branches; frequencies are called descriptive.
* **`verify_R2.py` T4.** It records solver status. Our rerun gives 147/494/1153 of 1794, plus 6 unclassified (status 4); your 153 differs, presumably by solver version.
* **`verify_R10.py`.** Union lengths are merged before measuring; 5a §3.4 now reports occupied length, hull span and the unmerged sum, and lists what the Monte Carlo does not test.

**[Check]** All outputs were regenerated.

### Outside scope: core documents (report §5)

These are also reported to the original referee of Phases 1–3 (`REVIEW_RESPONSE.md`).

1. **Restricted support (2a).** New notation $\operatorname{supp}_AP_X$, defined in 2a §0.1 with your example. [supp] is proved for it, and Lemma 1, Lemma 3(c), Remarks 5.1–5.2 and P1 A5 use it. This also resolves the tension between Lemma 3(c) and Remark 5.1. **[Theory] [Check]** `verify_phase2a.py`.
2. **Projection under a singular covariance (P3.5(b2)(ii)).** The proposition now requires a positive definite $\Sigma$ (ii-a), or a reduced full-rank parameterisation with exactly zero null-direction error (ii-b). It has a new proof, and your counterexample.
   * **[Code]** The package drops candidates whose gradients use a numerical null direction of $\hat V$, with a note when this happens.
   * In the Monte Carlo of P3 §8, $\hat V$ had full rank in every replication of every design, so this is a safeguard rather than a fix for an observed failure.
   * **[Check]** `verify_phase3.py`; package tests.
3. **Residual CLT (P3.5(a) step 5).** Completed by your Chebyshev/Lindeberg subsequence argument, credited; the allowance is unchanged. **[Theory] [Check]**
4. **1b attribution and pure units.** Fixed in 1b §§A1–A2. **[Theory]**

### Comparison with the earlier internal review (report §6)

Agreed. 5a §4 now records that the internal review missed:
* the face-LP intercept;
* set-level versus target-level redundancy;
* the empty LCE–CR example;
* the double-counted union length;
* everything in 2c, which it had not been asked to cover.

Our external review plan now gives internal pre-screens the same scope as the external reviewers. We also fixed the "would double the bound" wording in §3.1.

---

## 5. Remaining open items

1. **Conjecture N2.2′** (sharp N2 with conditioning on size). It is stated, with a proof route, and left open by decision for now.
2. **R15:** integer unit tables beyond 2×2 with margins only.
3. **R16:** regression sharpness and inference for the sign restrictions under district rungs. Until then the module is labelled experimental under those rungs.
4. **Width-gain inference.** It must be joint over the base and augmented endpoints (5a §2.5); not yet in the package.
5. **Conjecture P3.5(b3)** (the simulated CLR candidate collection) and uniformity near point identification remain open, as in P3 §10.

We are grateful for the replacement arguments, which we adopted with credit in the documents: the analytic proof for R10(e), the person-weighted binning lemma, the Chebyshev/Lindeberg completion, the null-direction counterexample, and the $V$-conditional formulation.
