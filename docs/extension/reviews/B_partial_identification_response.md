# Response to the referee report on Batch B (partial identification and inference)

*Reviewed version: commit `8b98589`. The revised documents are on branch `claude/wizardly-volta-lbh4bm`. Line numbers in the comments refer to the reviewed commit. Package changes are in release 0.3.0 of `eiPartialID`.*

## Cover note

We thank the referee for a careful and constructive report. We reproduced every counterexample: the referee's appendix script runs unchanged (`reviews/referee_checks_B_partial_identification.py`), and our own checks are in `reviews/checks/`. We agree with the recommendation of major revision for the theory, above all for `PHASE2c_monotone.md`. Of the 43 points we recorded (P-B1 to P-B43), 41 called for a change; all 41 were reproduced, and none is disputed. We argue the **severity** of two points below, while accepting the corrections.

**What changed, in brief.**
* **Scope labels.**
  * N1 (Prop. M2) and the DiD rung Δ4 are now stated as sharp for **continuous** unit tables given the margins. A new remark proves that in 2×2 with margins only, the integer gaps form a lattice of step $1/(N\mu(1-\mu))$ inside the continuous interval. Integer tables in general are a new open item (R15).
  * Theorems R3 and R10 are stated for the **pooled** model. Their grouped versions are a new open item (R16).
  * The package (0.3.0) labels the sign and reinforcement restrictions **experimental**, and under grouped rungs prints that each sign is one pooled restriction across groups and that sharpness is proved for the pooled model only.
* **Results restated or withdrawn.**
  * M3(d) now covers $\gamma<1$, $\gamma=1$ and $ER=0$.
  * M5 is restated rung by rung, with the exact between-district covariance term and the referee's reversal example.
  * N2.2(a) no longer claims a sharp version: the referee's $V$-conditional statement is recorded as **Conjecture N2.2′**, which we leave open for now.
  * N2.3 is kept only as a structural calculus identity; the local-derivative rows are **dropped** from the build plan, because the derivative is not identified on a finite support.
  * The binning step in 2c §6 is replaced by the referee's person-weighted repair, proved as **Lemma N2.4**. The degree of the polynomial is now prespecified.
* **5a errata.**
  * R2: the target intercept in the face-test LP; "zero incremental gain for this target" in place of "redundant".
  * R3(e): "rejects the conjunction; does not attribute"; the LCE-excludes-CR example (Remark R3.2).
  * R10: the referee's version table replaces our status note.
  * R10(e) is now a **theorem**, using the referee's analytic bound (0.933525; 0.940092 with 4,000 critical-value draws), with credit.
* **Core repairs reported to the original referee as well.**
  * the support of the restricted measure (2a);
  * the projection interval under a singular covariance (P3.5(b2)(ii), now requiring a positive definite covariance or a reduced parameterisation; the package drops candidates that use a null direction);
  * the residual CLT, completed by the referee's Chebyshev/Lindeberg argument.
* **Scripts.** Every correction has a check that fails on the old statement and passes on the new one.
  * `verify_phase2c.py`: 12 → 26 checks.
  * `verify_R2.py`: adds T1c and T4b.
  * `verify_R3.py`: adds S6 and S7.
  * `verify_R10.py`: adds V9–V12, with union lengths merged before measuring.
  * `verify_phase2a.py` and `verify_phase3.py`: one and two new counterexamples.

**The two severity points.**
1. **Integer tables (P-B1).** We accept the relabelling. In the case the package would use first (2×2, margins only), the integer and continuous answers differ by less than one lattice step, about $6\times10^-6$ at the scale of the Wisconsin application. The difference matters for small tables and for refutation exactly at the cut point, which the referee's three-unit example shows.
2. **$\gamma=1$ in M3(d) (P-B3).** We accept the correction. $\gamma=1$ requires every unit to be pure, in which case the gap is point-identified, so the omission did not affect any reported interval.

**What remains open.** Conjecture N2.2′ (sharp N2 with size conditioning); R15 (integer tables beyond 2×2 margins); R16 (grouped R3′/R10′); inference on a width gain jointly over the base and augmented endpoints (P-B22).

## Revised verdict table

| Result | Referee's verdict | Status now | Where |
|---|---|---|---|
| M1 | Established | Unchanged; hypotheses split by part | 2c §1 |
| M2 | Established for continuous rates; integer claims incorrect | **Proposition, continuous unit tables.** Remark M2.1: integer lattice in 2×2 with margins (proved); R15 open | 2c §2 |
| M3 | Incorrect as a whole ((d) at a boundary) | **Proposition** with stated hypotheses (continuous rates, no unit restrictions, no other sign rows); (d) covers all three cases | 2c §2 |
| M4 | Established | Unchanged | 2c §2 |
| M5 | Incorrect across advertised rungs | (a) all rungs; (b) district-specific ordering; (c) M0 and within each district, with the exact between-district term for pooled district rungs | 2c §3 |
| N2.1 | Established | Unchanged | 2c §4 |
| N2.2 | Established under conditions | (a) outer only; **Conjecture N2.2′** for the sharp version; (b) on the complete finite support; (e) chart and weighted covariances | 2c §4 |
| N2.3 | Not established as an observable restriction | Structural identity only; local rows dropped | 2c §4 |
| (new) N2.4 | — | **Lemma**, person-weighted bins give an outer relaxation (referee's repair) | 2c §4 |
| R3(a)–(e) | Established (with qualifications) | Pooled model; (b) finite LP for polyhedral $C$; Remarks R3.1, R3.2; (e) reworded | 5a §1 |
| R2(a)–(e) | Established under conditions; (e) wording | Step 4 intercept; (b)–(d) qualifications; (e) target-specific | 5a §2 |
| R10(a)–(d) | Established under stated assumptions | Version table: (a) theorem; (b1) cited under the augmented Condition P; (b2) theorem; (b3) algorithm | 5a §3 |
| R10(e) | Established; proof completion | **Theorem** (referee's analytic proof) | 5a §3 |

## Point-by-point

Each row gives the comment (summarised from the report), our check, the change as planned in `reviews/B_PI_RESPONSE_PLAN.md`, its status, and the check that now fails on the old text. The full ledger, with document and package locations, is in the plan.

### A. 2c: model-free sign bounds (N1), Props M2–M4 (major comment 1)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B1 | M2 "sharp in the finite-population sense" is sharp only for **continuous** unit tables. Integer counts can reject CR when the LP accepts it (3 units of size 4: the continuous CR set is {0}, the integer gaps are ±1/6). | **R\*.** Reproduced, and I1 generalises it. In 2×2 with margins only, every integer sum is attained, so the integer $D$-set is the DD interval on a lattice of step $1/N_1+1/N_0$. Integer CR is empty only when the continuous CR width $(1-\gamma)\lvert ER\rvert$ is below one step (639 of 2,987 tiny populations, all of that kind). Endpoints differ by less than one step. The step is $6\times10^{-6}$ at 370 units of 500–4,000 people. Severity is high only for tiny tables or near-ties. With unit restrictions $U_i$, or joint columns in R×C, the gap is not bounded by this argument (open). | Relabel M2: sharp for continuous tables conditional on the margins. Add a remark with the lattice statement for 2×2 and the referee's example. Integer N1 is decision **U1**. | Done | **C7:** enumerate integer allocations for the referee's example. The integer CR set is empty while the M2 LP returns {0}, contradicting "attained by unit tables". |
| P-B2 | M2(d): a scalar breakdown value needs a specified conclusion, equality treatment and a scalar path. A vector of slacks gives a breakdown frontier. | **R.** | Restate (d) as nesting componentwise in ε, with a breakdown frontier in general and a scalar value along a stated path. | Done | **C7b:** two slacks for which the set of ε pairs preserving a conclusion is not a product of intervals. |
| P-B3 | M3(d) "exactly one of $F_\pm$ is non-empty when $D_{ER}\ne0$" fails at γ = 1 (pure units (0,0), (1,1): both branches contain the solution). | **R\*.** Reproduced. The proof notices the case (l.127) and the statement drops it. Minor, because γ = 1 means every unit is pure and $D$ is point-identified. | Restate as: exactly one branch when $ER\ne0$ and $\gamma<1$; both give the identified gap at γ = 1; gap 0 at $ER=0$. | Done | **C8:** two pure units; both LP branches are feasible, contradicting "exactly one". |
| P-B4 | M3(c)'s witness $\beta_{i1}=\beta_{i2}=T_i$ can be removed by unit restrictions $U_i$ or by previously imposed signs. | **R.** M3 is stated with $F$ from Def. 2, which includes $U_i$, but the proof ignores $U_i$. | Add "with no unit-level restrictions and no other sign rows" to M3. | Done | **C8b:** a unit ordering $\beta_{i1}\ge\beta_{i2}+0.1$ that makes $\{\delta_B\le0\}$ empty. |
| P-B5 | "CR is never refuted by the data" and "Between signs and CR can never be rejected from aggregates" (2c M3(d), Reading) fail for integer tables, and are read into the wrong base model elsewhere (see P-B14). | **R.** The single between sign still has an integer witness. CR does not (P-B1). | Reword everywhere: "in the continuous-rate model-free problem (Duncan–Davis tables, no linear rung), a $\delta_B$ sign and CR are always feasible; with the linear rung they can be jointly rejected; with integer counts CR can be rejected outright." In the emptiness note, say whether the rung alone is rejected (already computed, `.rxc_empty_test`). | Done | **Package test:** the emptiness note no longer contains "never be refuted" unqualified. **C7** (integer CR empty). |
| P-B6 | M4 established. | R | No action. Cite the referee's two-column example in the C3 discussion. | No change needed (result established); example noted | — |

### B. 2c: K-rows under districts, Prop. M5 (major comment 2)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B7 | M5 is incorrect across the advertised rungs. (a) holds per rung only with the full district parameterisation; the single $M_{Q+K}(X_i)$ formula proves only the pooled case. (b) needs district-specific ordering on each district's support. (c) fails under M1: six units in two districts with common positive slopes give pooled $\delta^*_W=-0.0493$. | **R.** Reproduced. I2 adds three things. (i) Within each district $\delta^*_W=+0.00064$. (ii) The 2×2 local rows $xm_1'+(1-x)m_0'\ge0$ hold everywhere. So PD's P-N chain "signed slopes ⊂ weighted-local ⊂ $\delta_W\ge0$" also fails for pooled M1. (iii) The W1 identity $\delta^*_W=\operatorname{Var}_p(X)(ER^*-D^*)$ still holds with the district-fitted τ, so the CR slab rows remain valid. | Restate M5: (a) per rung with the stacked district parameterisation; (b) per district; (c) for M0, or with district-conditional weights; give the between-district covariance term $\sum_g\mu_g\operatorname{Cov}_{w^g}(X,\zeta_{g,D})$ for M1. Scope P-N and the §4.2 cells to M0. New proof item **R16**: R3′/R10′ for grouped rungs via 2b Thm 1′ (decision **U2**). | Done | **C9:** the referee's six-unit design; M5(c)'s hypothesis holds and the conclusion fails. **Package test:** with `group` given, the sign/CR print line does not claim R3(b) until R16 is proved. |
| P-B8 | Remark M5.1: "order $p^{-1/2}$" for $\delta-\delta^*$ needs sampling and weight conditions; (S2) is a mean condition only. | **R.** | Replace with: mean zero given the design; its scale is the root sum of squared design coefficients under P3-S1a, with a weight condition for concentration. Or delete the rate. | Done | — (wording) |

### C. 2c: N2 and §6 (major comments 3 and 4)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B9 | N2.2(a)–(b): "support points (the distinct observed compositions)" is ambiguous. On a subset of a larger support, extendibility subject to the other accounting equations can fail (support {1/4, 1/2, 3/4}, dropping the middle point). | **R.** Reproduced. | Define N2 on the finite empirical design law, where all of (a)–(d) hold. For grids, state only an outer relaxation, using the binning repair (P-B13). | Done | **C10:** the referee's three-point example. The LP on the two extreme points accepts values that no monotone extension can produce at the middle point. |
| P-B10 | N2.2(a): the proposed sharp version $E[S(X,T)\mid X=x]$ is not sharp. The model conditions on $V=(X,N)$. The condition is $m(X)\in E[S(X,T)\mid V]$ a.s., with the shape restrictions. | **R** for the counterexample (two sizes at $X=1/2$). **N** for the replacement: it must be proved via 2a Thm 1 selection for each measurable $m$. | Replace the sentence with the referee's $V$-conditional statement, marked Conjecture until proved. Add the size example as a remark. Carry the caution into R9. | Partly: the claim is withdrawn and the $V$-conditional statement added as Conjecture N2.2′, left open | **C10b:** the size example; the X-conditional LP admits $(0.75,0.25)$, which the full model excludes. |
| P-B11 | N2.3: the derivative is not identified on isolated support points. Adding $h=(x-\tfrac14)(x-\tfrac12)(x-\tfrac34)$ leaves the law unchanged and moves $\tau'(1/2)$ by $-1/16$. An LP using the wrong $\tau'$ need not cover the true gap. | **R.** Reproduced. | Restate N2.3 for the structural $\tau_m$. As an observable restriction, require an identified derivative: an interior point of an interval support with compatible smooth versions, or a maintained identified regression family. An estimated derivative adds sampling error. Decision **U4**. | Done | **C11:** the referee's $h$ example; two compatible models give different local thresholds. |
| P-B12 | N2.3(d)'s last assertion is true even for smooth functions ($X\sim U[1/5,2/5]$, $m_1=x-x^2/2$, $m_0=1-x^2/2$: local zero, aggregate $-1/47250$). The reverse non-implication also holds. | **R.** Reproduced ($-2.1164\times10^{-5}$). | Replace the kinked V6 witness by the smooth example and add the reverse direction. | Done | **C13:** the smooth example (strengthens the check; nothing fails). |
| P-B13 | 2c §6 binning: "averaging Lemma N2.1 within a bin uses [CE-cvx]" is wrong, since $S(x,\tau(x))$ varies with $x$ (example: 0.3125 against 0.25). Repair: person-weighted bin means $\bar m_{1b}=E_b[Xm_1]/\bar x_b$ etc. give an outer relaxation for global monotonicity. Local rows need a separate argument. | **R** for the counterexample (0.3125, 0.25, 0.625 reproduced). **N** for the repair proof. | Replace the binning sentence by the repair, proved as a lemma (with absent groups by zero-weight extension). | Done | **C12:** the two-point bin example. |
| P-B14 | 2c §6: fixed-degree inference after LOO-CV selection of the degree is invalid (conditional coverage about 68%). Degree $J$ with $J$ support points is not identified; $J$ is used twice. | **R.** I3: $1-0.05/P(\lvert Z\rvert>\sqrt2)=0.682$. | Prespecify the degree, or use sample splitting or a selection-aware procedure. Rename the support size (e.g. $J_X$) apart from the degree. | Done | **C12b:** LOO selection simulation in the balanced two-level design; conditional coverage ≈ 0.68. |

### D. 5a: Theorem R3 (major comment 5)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B15 | R3(b): calling the computation "an LP" needs polyhedral $C$ (a ball is admissible under (C1)–(C2)). (C2) is stronger than (a) needs. The dimension condition is sufficient (A5 implies it). $P(X\in A)=1$ is needed for the witness. "Fixed design ⇒ $\tilde B=B^*$" is a convention, not an identity for a superpopulation rate. The full-law versus regression distinction is illustrated by $C=[0.2,0.4]$: full law empty, regression $[1/6,1/3]$. | **R.** Reproduced. | Say "a finite LP when $C$ is polyhedral". Add the full-law versus regression example as a remark under (b). Reword l.22 as a convention. | Done | **S6:** the referee's $C=[0.2,0.4]$ example; regression interval $[1/6,1/3]$, full-law set empty. Shows (a) and (b) differ; nothing fails. |
| P-B16 | R3(e): "never $C$ alone" needs qualifying when a belief has a model-free implication ($m_r=1.5-2x$: ER = −2, positive within association impossible model-free). 5a l.309 "rests on an unrefutable belief (2c M3)" imports M3 into the wrong base: on $X=0,1/4,1$ LCE fits, and both LCE CR branches are empty ($\tilde D=-261/1925$, ER = 0). | **R.** Reproduced. | R3(e): "an empty LP rejects the conjunction; emptiness alone does not attribute the failure". l.309: name the base model, and say that under LCE CR is refutable even when LCE fits. Package: print "LCE alone is / is not rejected" next to the joint rejection. | Done | **S7:** the $X=0,1/4,1$ example; both LCE branches empty while LCE fits, contradicting "unrefutable" as applied in §3.5. |
| P-B17 | R3(c), (d) established. Minor: a common witness is only sufficient. In the square with $g_B=u+0.5$, $g_W=u-0.5$, target $v$ has identical branch images although the common-witness LP is infeasible. | **R.** | Add the referee's example to §1.4. | Done | — |

### E. 5a: Proposition R2 (major comment 6)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B18 | R2 step 4 (l.169): the face LP $h^\top k\ge U$ omits the target intercept. It should be $h^\top k\ge U-\varphi_0$. On $[0,1]^2$ with target $2+x$ the printed row is $x\ge3$. | **R** for the document. **Code is correct:** `.rxc_face_moves` subtracts $c(\hat\beta)+c_0$ (`check_B_pi_code.R`: a cut through the facet keeps $U$; a cut removing it moves $U$; the document's row is infeasible). `verify_R2.py` uses zero intercept, so it cannot see the error. | Fix step 4. | Done | **T1c:** target $2+x$ on the square; the printed face LP is infeasible (so it reports "moves") while the brute-force endpoint is unchanged. **Package test:** `.rxc_face_moves` with $c_0=2$. |
| P-B19 | R2(b): a positive multiplier certifies tightening only with uniqueness or a positive lower bound on all optimal multipliers (duplicate-row example). R2(c): the separating witness need not be a reported cell-rate target. | **R.** | Add the duplicate-row example to (b). Add "the witness objective may not be a reported target" to (c). | Done | — |
| P-B20 | R2(d) holds under γ < 1, a nonconstant gap target, and nonempty intersections. §2.5 item 2 ("in 2×2 the middle verdict cannot occur") needs the nonconstant-target qualification. | **R.** | Add the qualification. | Done | — |
| P-B21 | R2(e): "redundant given the others" is incorrect if it means implied. Square, target $x$, rows $x+y\le0.8$ and $x\le0.5$: zero last-added gain for $x$, yet the row lowers the upper endpoint of $x+y$ from 1.5 to 0.8. | **R** for the document. **Code is correct:** the added-last verdict is "cuts the set, not this target" (checked on the same example, `check_B_pi_code.R`). | Reword (e): "zero incremental gain for this target"; "redundant for the feasible set" only via the implication test (c). | Done | **T4b:** the square example; the old sentence calls the row redundant while it tightens $x+y$. |
| P-B22 | §2.5 item 4: width-gain inference needs joint treatment of the base and augmented endpoints (and any estimated threshold). An interval for one endpoint is insufficient. | **R.** | Restate item 4 and the open item. | Text corrected; the joint-endpoint width-gain inference is still to be built | — (design item) |

### F. 5a: Theorem R10 (major comment 7)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B23 | R10(a) holds per version only under that version's conditions for the **augmented** system. Growing-district M1/M2 (regime R3) is conjectural (P3 §0.5) independently of candidate collection. "Finite vertex set at each $n$" does not verify CLR Condition P (fixed compact index set, gradients bounded away from zero, well-conditioned covariance, rate). A zero-variance **candidate** is different from an exactly known row. Appending row standard errors to the old CLR maximum is insufficient. | **R.** The package already collects candidates from the augmented system (`.rxc_collect(cons, …)` uses the full `cons`). Zero-s.e. candidates are zeroed in the CLR maximum (`W[, se == 0] <- 0`) and kept as exact bounds. That is right only if their gradient is exactly zero. | Replace "R1/R3" by the P3 §0.5 table. List Condition P's actual requirements in the status note. Add the zero-variance-candidate sentence. | Done | — (status) **Package test:** a candidate with nonzero gradient but numerically zero $\hat\sigma$ is not treated as exact. |
| P-B24 | R10(a) for actual targets: carry over P3-S1a, the weight condition and the residual allowance. An actual gap needs its own residual envelope, e.g. $\tfrac12\{\sum_i(w_{i1}+w_{i0})^2\}^{1/2}$. A single-row $S_1$ need not bound the contrast. | **R.** **Code is valid:** `ci_actual` adds the per-cell envelopes $\sum_{rc}\lvert a_{rc}\rvert\,\mathrm{sd}_{rc}$ (`R/rxc_infer.R:388`). By Minkowski this dominates the referee's contrast envelope, so it is valid and more conservative. | State the allowance for linear contrasts (sum of per-cell envelopes, or the joint envelope). | Done | — |
| P-B25 | R10(b), (c), (d) established under conditions. One simultaneous event includes the target constant (no strict feasibility needed). Branch union: point coverage only, no simultaneous coverage of both branch images, no selection on $\widehat{ER}$. Under triangular arrays, branch validity is needed along each subsequence. Emptiness bound $\alpha_U+\alpha_L$. Pointwise only. | **R.** The code's conservative slab includes the target row (`Gall <- rbind(cons$S, DS, tg$cvec)`, two-sided via `abs`). | Adopt the referee's version table (R10 against P3.5(a), (b1), (b2)(i), (b2)(ii), (b3)) as §3.2's status note. Add the triangular-array sentence to (c). | Done | — |
| P-B26 | R10(e): an analytic counterexample replaces the Monte Carlo sketch. A repeated three-point design gives $\limsup P\le\Phi(\tfrac7{12}z_{.995})=0.9335$; with 4,000 critical-value draws, ≤ 0.9401. The supplied far-edge design also satisfies the condition (ratio 0.564, ≤ 0.927). | **R.** Reproduced. The derivation of the 0.564 ratio for the far-edge design is **N** (not rerun here). | Replace the sketch by the fixed-design proof. Keep the MC as illustration. | Done | **V9:** the three-point design; assert the analytic bound is < 0.95. (Nothing fails; it completes a proof.) |
| P-B27 | 5a l.311: "the naive interval contains bracket ∩ JKST CI, so under-coverage transfers a fortiori" is unproved, and in general false (an empty pre-expansion intersection need not stay empty after expansion). Direct argument: the narrower procedure has limiting coverage ≤ 1/2 at the edge. Also distinguish the two sample ERs when weights differ. | **R.** | Replace the a-fortiori sentence with the direct ≤ 1/2 argument and the counterexample $[0,0.1]$, $[0.15,0.2]$, allowance 0.1. | Done | **V10:** the referee's three intervals; the claimed containment fails. |

### G. Minor comments and the EGGHM mapping (report §3)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B28 | W1 versus W2: invertibility of $M$ is needed for $g$ only; ER needs only positive regressor variance; choose a bounded value for the absent rest group at $X_r=1$; W2 signs concern components of a sum and need not sign every summand. | **R** (M1 boundary identities reproduced). | Edit the hypotheses of M1(b)–(c); add the W2 sentence. | Done | — |
| P-B29 | Ladder (2c §5): "both signs" ⇒ CR only when the directions agree; N2 admits steps where a local derivative is undefined; compare realised and model-implied associations only on the same population; for $R\ge3$ specify a simplex chart. The symmetric Dirichlet example gives negative covariance under both W1 weights. | **R.** | Fix the ladder lines. Upgrade C5(c) to the weighted covariances. | Done | **C5(c′):** the weighted covariances under both W1 weights are negative (strengthens the check). |
| P-B30 | Identification versus sample labels: R3(e) is a population statement; raw LP infeasibility at $\hat\beta$ is a sample diagnostic; level-controlled rejection uses R10's rule. R2 §2.5 l.208 and l.210 need their qualifications (P-B20, P-B22). | **R.** | Rename to "empty at the point estimate (see the emptiness test)". | Done | **Package test:** verdict strings contain no "refuted" for point-level emptiness. |
| P-B31 | EGGHM §2, Lemma 1, Props 1, 4, 6, 9–10: correct. EGGHM assumes a mixed neighbourhood (γ < 1); ER also needs positive composition variance; the strict branch distinction at nonzero ER needs a mixed unit. | **R.** | Add the mixed-unit condition where EGGHM is cited for branch uniqueness. | Done | — |
| P-B32 | EGGHM Props 2–3, 5, 7, 8, Cor. 3: sharpness ranges over compatible joint distributions; it does not license integer-table interpolation. sign(ER − NM) = sign(ER) uses γ < 1. | **R.** | Cite this with P-B1. | Done | — |
| P-B33 | EGGHM §4 end, Theorems 1–3: 1b l.69 attributes all three tables to Theorem 1 (the gap is Theorem 1; the group levels are Theorems 2 and 3). W2 and joint columns need their own proofs (M1, M4). | **R.** 1b l.69 checked; EGGHM theorem numbers taken from the report. | Fix the attribution. | Done | — |
| P-B34 | EGGHM §5, Props 11–13 and 14–16: the neighbourhood readings use Assumption 1 and meaningful derivatives. N2 imposes a common conditional mean. Neither identifies a derivative on an atom. At a pure composition EGGHM's local between association is 0; replacing it by the sign of $D(x)$ changes the belief. | **R.** | With P-B11. | Done | — |

### H. The referee's notes on our verification scripts (report §4)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B35 | `verify_phase2c.py` C4's ordering construction changes one column without rebalancing the others, so it does not always build a valid full transition matrix. No check covers integer feasibility, district intercepts or derivative identification. | **N** (C4 not re-read here). | Rebuild C4 with valid transition matrices. Add C7–C13 above. | Done (script) | C7–C13 |
| P-B36 | `verify_R3.py` S2 labels a branch nonempty only if its width exceeds $10^{-9}$ and calls zero positive-width branches "both degenerate" (possibly two empty sets). Frequencies over reused designs are descriptive. | **R.** Code read: `verify_R3.py` l.102. | Separate empty, degenerate and positive-width branches. Call the frequencies descriptive. | Done (script) | — |
| P-B37 | `verify_R2.py` T4 now gives 153/494/1,153 of 1,800, against 147/494/1,153 of 1,794 saved. Unsuccessful base solves are skipped without recording status. The face LP in the tests uses zero intercept. | **N** (not rerun; the earlier session reproduced the committed output under the old runner). | Record solver status for skipped cases; regenerate the output; add T1c (P-B18). | Done (script) | T1c |
| P-B38 | `verify_R10.py` l.189–190 sum the two branch lengths even when they overlap, so "union length" can exceed the Lebesgue length. The runs do not test CLR, residual coverage, the crossing rule or growing dimension. | **R** (code read). The package merges overlaps (`.rxc_union`). | Merge before measuring; regenerate; correct the 5a §3.4 figure. State what the MC does not test. | Done (script) | **V11:** overlapping branch intervals; the script's length exceeds the merged length. |

### I. Outside scope: core documents (report §5)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B39 | 2a tag [supp] (l.57), Lemma 1 (l.144–150), Lemma 3(c) (l.310–314): continuous functions equal a.e. **on** $A$ agree on the support of the measure restricted to $A$, which can be smaller than $\operatorname{supp}P_X\cap A$ ($X\sim U[0,1]$, $A=\{1/2\}$). This also explains the tension between Lemma 3(c) and Remark 5.1. | **R.** Elementary counterexample. Note that Remark 5.2's argument (A5 ⇒ full-dimensional) also uses $\operatorname{supp}P_X\cap A$. | Define $\operatorname{supp}(P_X\vert_A)$ and impose A5 on it. Rewrite [supp]. Check every use. Report it to referee A with the batch-A material (it corrects a reviewed core result). | Done (also reported to referee A) | **2a check:** $X\sim U[0,1]$, $A=\{1/2\}$; the old [supp] claims $0\equiv1$ at $1/2$. |
| P-B40 | P3.5(b2)(ii) (l.179–189): the rank-$q$ pseudoinverse proof bounds gradients in the range of $\Sigma$ only. A singular CLT need not control smaller errors in null directions, e.g. error $(Z/\sqrt n,-1/n)$. | **R** (I4). | Document: require $\Sigma$ positive definite (as (b1) already does), or work in a full-rank reduced parameterisation, or require exact null-coordinate identities. Code: if $q<\dim\beta$, use $q=\dim\beta$ and warn, or set the projection to NA with a message. Add a test. | Done (proof and package guard; also reported to referee A) | **P3 check:** I4's toy; zero-s.e. bound misses always. **Package test:** a rank-deficient $\hat V$ triggers the guard. |
| P-B41 | P3 [Lyap] (l.40) and the residual step (l.201–203): $\max_iw_i^2/\sum_iw_i^2\to0$ does not give a variance-standardised CLT if only one summand has positive variance. The tail bound with the envelope $S_1$ still holds: Chebyshev where $\mathrm{sd}(E)/S_1\to0$, Lindeberg where the ratio is bounded away from zero. | **R** for the gap. **N** for writing the argument. | Replace step 5 by the subsequence argument. Report to referee A with batch A. | Done (also reported to referee A) | — (proof completion) |
| P-B42 | 1b §§A1–A2 (l.59–80, 116–121): fix the theorem-to-table mapping (P-B33). For pure units, unit-level row ordering cannot be equated with the observable implication of a local between sign (the between covariance is identically zero). | **R.** | Fix the mapping. Add the pure-unit caveat. | Done | — |

### J. The referee's comparison with the earlier internal review (report §6)

| Point | Comment (summary) | Our check | Change | Status | Check that fails on the old text |
|---|---|---|---|---|---|
| P-B43 | The internal review (5a §4) covered 5a §§1–3 only and missed the face-LP intercept, set-versus-target redundancy, the empty LCE–CR example and the union-length double count. Its R10(e) sketch left verification to simulation. Combining two tests by a union bound gives an upper bound; the rejection probability need not double. Inference labels should keep P3's growing-district qualification separate from the candidate-simulation conjecture. | **R.** | Reword the §3.1 parenthesis ("would double the **bound**", not the level). Add to 5a §4 a line pointing to this ledger. Process lesson for EXTERNAL_REVIEW_PLAN §3: the fresh-session pre-screen must also cover the scope document (here 2c), not only the new proofs. | Done | — |

## Acknowledgements

We adopted, with credit in the documents, the referee's analytic proof of persistent under-coverage for R10(e), the person-weighted binning repair (Lemma N2.4), the Chebyshev/Lindeberg completion of the residual step, the null-direction counterexample for the projection interval, and the $V$-conditional formulation recorded as Conjecture N2.2′.
