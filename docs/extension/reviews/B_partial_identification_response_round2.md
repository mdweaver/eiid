# Response to the referee, second round: Batch B, partial identification and inference

**Reviewed version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b).
**Revised version:** commit [`1b8b6fb`](https://github.com/mdweaver/eiid/tree/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f) on branch `claude/wizardly-volta-lbh4bm`.
**All changes since the reviewed version:** [compare `e170be9…1b8b6fb`](https://github.com/mdweaver/eiid/compare/e170be9161b9bbc250a17a27ffd39bc523c1c88b...1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f).
**Your report, as received:** [`reviews/B_partial_identification_report_round2.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/reviews/B_partial_identification_report_round2.md).

Abbreviations as before: **2c** is [`PHASE2c_monotone.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE2c_monotone.md), **5a** is [`PHASE5a_sign_module.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE5a_sign_module.md), **P3** is [`PHASE3_inference.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE3_inference.md). Corrections carry the note "*(Corrected after review, round 2, Nov 2026 …)*". Tags: **[Theory]**, **[Code]** (package 0.3.1), **[Docs]**, **[Check]**.

---

## 1. Summary

Thank you for the second report and the recommendation of minor revision. We made all four corrections and every minor change you asked for, with no disputes. We also accept your criticism of how our first response described the new checks. Every check is now tagged by what it actually shows (§4).

## 2. What changed

| Where | Change | Your item |
|---|---|---|
| 2c "Reading" after M3/M4 | Equal-rates witness gives every weak between sign and, for each pairing, the CR **union**; which branch depends on its within association. Your three-row example added | 1 |
| 5a Theorem R10 | Scope is pooled M0 only. Fixed-district P3 intervals are background; district rungs are R16 | 2 |
| P3 §5, remark after (b2)(ii); P3 §8 | (ii-b) only for a documented exact reduction; numerical cutoff → $q=\dim\beta$, dropping null-direction candidates is a safeguard | 3 |
| 2c N2.3(d) | "Nondecreasing **on a neighbourhood of $x$**", with your $.5+h$ example | 4 |
| 5a §1.1 | "Finite LP" narrowed to the regression-only optimisation over $\Theta_O\cap C$ | minor |
| 2c M4 example | Same weak sign on complementary columns forces $D=ER$; opposite signs duplicate one restriction | minor |
| 2c M3(d), proof step 6 | Branch versus its singleton gap image | minor |
| 2c M2(d) | The set of slack values preserving a conclusion versus its boundary, the frontier | minor |
| 2c Remark M2.1 | Lattice step $1/N_1+1/N_0$ needs both group totals; numerical example | disposition table |
| 2c Lemma N2.4, §6 | Model-implied bin mean $\bar\tau_b$, estimated by the observed bin mean; bias allowance or specification-conditional label | disposition table |
| 2c Remark N2.2.1 | Your sufficiency argument added as a credited remark; Conjecture N2.2′ stays open | disposition table |
| 2c, 5a, P3, 2a, 2b, R2/R3/R10 scripts | Every Batch B and round-2 check tagged [F], [V] or [I] | checks paragraph |
| Package 0.3.1 | Projection critical value uses $q=\dim\beta$; `print()` states the CI status per fit (R10 for pooled M0; no theorem for grouped sign fits) | 2, 3 |
| Checks | 2c 26 → 29 (C3b, C15, C16); new V13 in `verify_R10.py`; all pass | |

## 3. Point-by-point

### Correction 1: the equal-rates witness

Agreed; our sentence was wrong. Your example reproduces: $ER=1$, $\gamma=.24$, $D=NM=.24$, $\delta_B=0$, $\delta_W=.0456>0$, and the negative branch is empty.

* **[Theory]** 2c now reads: "The equal-rates table … sets every pairwise between association to zero. It therefore satisfies every weak between sign and, for each specified pairing, the CR **union**, in the continuous problem without unit restrictions. Which branch contains it depends on its within association; both do only if that is also zero." This is followed by your three-row example and a correction note.
* **[Check]** `verify_phase2c.py`, C15 [F].
* **[Docs]** We searched the guide, the vignette and the package strings for the "both branches" claim and found no other copy.

### Correction 2: the scope of R10

Agreed.

* **[Theory]** 5a Theorem R10 now says "**Pooled M0.**" and nothing more. The fixed-district baseline intervals of P3 (M1/M2 in regime R1) are named as background that covers the rung alone. The added-row, slab and union statements are explicitly not claimed for district rungs (R16).
* **[Theory]** 2c M5 gains a sentence on the CI status of each kind of fit, and 5a §5 a status row.
* **[Code]** In 0.3.1, `print()` states which result covers each fit's intervals:
  * pooled M0: P3.5(b2)(ii), and R10 when sign or CR rows are present;
  * grouped fits without signs: only with few, large groups (regime R1);
  * grouped fits with signs or CR: "no coverage theorem for sign/CR restrictions under grouped rungs (R16 open)".

  A grouped M-SIGN fit also prints a note that its confidence intervals, projection included, are not covered.

### Correction 3: an exact rank reduction versus a numerical diagnostic

Agreed. Our first-round letter called the drop rule "case (ii-b)". That was wrong for the reason you give: a numerical cutoff establishes neither $\beta=L\eta+b$ nor $\hat\beta=L\hat\eta+b$, and it does not calibrate a rank-based critical value.

* **[Theory]** P3 §5 has a new remark after (b2)(ii), "exact versus numerical reduction". Case (ii-b) is claimed only for a documented exact reduction, such as the adding-up identity. With a numerical cutoff the theorem applies with $q=\dim\beta$ (case (ii-a), true $\Sigma$ positive definite). Dropping candidates that load on a numerically null direction is a safeguard that can only widen the interval. P3 §8 and the 5a §5 guard row are corrected to match.
* **[Code]** In 0.3.1 the projection critical value uses $q=\dim\beta$ (`R/rxc_infer.R`). The numerical rank of $\hat V$ is printed only as a diagnostic, and the rank guard is labelled a safeguard. In every reference design $\hat V$ has full rank, so no reference output changed. A test with a rank-deficient toy covariance shows the interval widening, from ±3.0 to ±3.4.

### Correction 4: the premise of N2.3(d)

Agreed.

* **[Theory]** 2c N2.3(d): "If $m_{1c}$ and $m_{2c}$ are nondecreasing on a neighbourhood of $x$ … Monotonicity on the finite support (N2) is not enough", followed by your construction (derivative $-1/16$ at $1/2$).
* **[Check]** C16 [F]. The local rows had already been dropped from the package (decision U4).

### Minor points

* **[Theory] 5a §1.1.** Now reads: "the regression-only computation, optimisation over $\Theta_O\cap C$, is a finite LP; the full-law $\Theta_I\cap C$ of R3(a) need not have a finite representation even for polyhedral $C$."
* **[Theory] 2c M4 example.** Thank you for catching the sign description. It now says that for complementary columns $\delta_{W,2}=-\delta_{W,1}$, so the **same** weak sign on both forces both associations to zero and gives $D=ER$, while opposite signs duplicate one restriction. **[Check]** C3b [F].
* **[Theory] 2c M3(d) and proof step 6.** "Both branches nonempty, with the same gap image $\{D_{ER}\}$" (and $\{0\}$ at $ER=0$).
* **[Theory] 2c M2(d).** The set of slack values preserving a conclusion is a sublevel set and need not be a product of intervals. Its boundary is the breakdown frontier. Check C7b relabelled.

### The disposition table

* **[Theory] Lattice step.** 2c Remark M2.1 gives the step as $1/N_1+1/N_0=1/(N\mu(1-\mu))$, with an example ($N=10^6$, $\mu=.3$: $4.8\times10^{-6}$), and says a numerical value for an application needs both group totals.
* **[Theory] Binning.** 2c Lemma N2.4 now defines the model-implied bin mean $\bar\tau_b=E_b[\tau(X)]$. The observed person-weighted bin mean estimates it, and its error needs inference. In §6, the prespecified degree is valid only if $\tau$ is a polynomial of that degree. Otherwise the bound needs an approximation-bias allowance or a label saying it is conditional on the specification.
* **[Theory] N2.2′.** Thank you for the sufficiency argument. It is added to 2c Remark N2.2.1 as a remark credited to you. Under the selection-expectation reading, take for each positive-mass $(x,n)$ a measurable selection with conditional mean equal to the common $m(x)$, then glue across the finitely many atoms. Following an earlier decision we keep the label **Conjecture N2.2′**, open for two reasons:
  * necessity in general;
  * whether any procedure using only $E[T\mid X]$ attains the set. As you note, it does not when size is informative.
* **[Theory] Line 30** is covered under the minor points above.

### The description of the checks

> *The response describes every new check as failing on an earlier assertion. That overstates the evidence.*

Agreed, and we withdraw that sentence. The 2c and 5a headers now define three tags, and every Batch B and round-2 check carries one:

* **[F] fails on an earlier statement:**
  * 2c: C3b, C7, C8, C8b, C9, C10, C10b, C11, C12, C12b, C14, C15, C16;
  * 2b: D8b;
  * 2a: [supp];
  * P3: the P3.5(b2)(ii) null direction, the P3.5(a) residual CLT, and P3.4′;
  * R2: T1c and T4b;
  * R3: S7;
  * R10: V10.
* **[V] verifies a corrected or new statement:**
  * 2c: the rebuilt C4, C5(c′), C7a and C7b;
  * P3: P3.5(c);
  * R3: S6;
  * R10: V9, V11 and V13.
* **[I] illustration:** C13 and V12.

As you say, the proofs carry the generality; the checks support individual instances.

## 4. Checks run on the revised version

| Script | Result |
|---|---|
| `verify_phase2c.py` | 29/29 PASS (was 26) |
| `verify_phase2a.py` | 32/32 PASS |
| `verify_phase2b.py` | 22/22 PASS |
| `verify_phase3.py` | 17/17 PASS |
| `verify_R2.py`, `verify_R3.py` | all asserted checks PASS; descriptive tables unchanged |
| `verify_R10.py` | V9, V10, V11, V13 PASS; limiting bound 0.933525 and 0.940092 as you found; main table unchanged |
| Package tests (`R CMD check`) | 0 failures |

## 5. Still open

These are unchanged in status, and each is visible in the documents and in the package's printed notes:

* general integer tables (R15);
* grouped-rung proofs for M-SIGN (R16);
* inference on width gains;
* Conjecture P3.5(b3);
* Conjecture N2.2′.

The package keeps M-SIGN labelled experimental.
