# Plan: responding to the batch-B referee reports (partial identification; ecological inference)

*Drafted 18 Nov 2026 from branch `claude/wizardly-volta-lbh4bm` at `27557df`. The report is `B_partial_identification_report.md` in this folder. It reviewed commit `8b98589`; its line numbers refer to that commit. The files it reviewed (2c, 5a, and the background parts of 2a, 2b and P3) are unchanged at `27557df`. This is an assessment and plan only: no theory document, R file or test was changed. The ledger format follows `REVIEW_RESPONSE_PLAN.md`.*

*Updated the same day with the second batch-B report, from an ecological-inference specialist: `B_ecological_inference_report.md` (also major revision; same commit `8b98589`). Its appendix is saved as `referee_checks_B_ecological_inference.py`. §1A summarises it, §2 part II (rows P-E1 to P-E33) is its ledger, and §§3–5 now cover both reports. The file name is kept so earlier links work.*

**Abbreviations.** **2c** = `PHASE2c_monotone.md`, **5a** = `PHASE5a_sign_module.md`, **P3** = `PHASE3_inference.md`, **1b** = `PHASE1b_monotone_did.md`, **PD** = `PACKAGE_DESIGN.md`, **PLG** = `PLAIN_LANGUAGE_GUIDE.md`. Package files are under `package/eiPartialID/`.

**Checks run for this plan.**

| Script | What it does | Result |
|---|---|---|
| `referee_checks_B_partial_identification.py` (the referee's appendix) | all counterexamples and bounds in the report | **reproduces exactly** under `python3 -I` (NumPy 2.x, SciPy 1.x): M1–M5, N2, R3, R2, R10 lines and the bounds 0.933525 and 0.940092 |
| `checks/check_B_pi_claims.py` (new) | I1: integer versus continuous N1 in 2×2 by brute force; I2: district reversal against the local rows and the W1 slab identity; I3: conditional coverage after CV selection; I4: singular-covariance projection; I5, I6 (added for the EI report): the W1/W2 sign reversal and the two-precinct integer example, recomputed independently | all pass (`python3 -I`) |
| `referee_checks_B_ecological_inference.py` (the EI referee's appendix; run from `docs/extension/`) | union lengths in `verify_R10.py`; the false-CR experiment; the W1/W2 table; integer gaps | **reproduces exactly** (28 s): union medians 0.357, 0.242, 0.480 with 134, 591, 0 overlapping runs; false CR at 600 units covers 398, 110, 0 of 400; W1 −0.02633, W2 +0.00110; integer gaps −2/5, 2/15, 2/3 |
| `checks/check_B_pi_code.R` (new; `Rscript` from the repository root, needs only `lpSolve`) | calls the package's internal face test, implication test and union helpers on the referee's toy problems; reads the M-SIGN proof-status labels | all pass; confirms the code findings below |

---

## 1. Summary

**Do I agree with "major revision"?** Yes, for the theory documents, and above all for 2c. Every counterexample in the report reproduces. I found no mathematical claim in the report that is wrong. I dispute the severity of two points:
* The integer-table point (P-B1). I accept it as a theorem-level correction. In 2×2 with margins only, though, the integer and continuous N1 sets differ by less than one lattice step in $D$. That step is $1/(N\mu(1-\mu))$, where $N$ is the total population. It is about $6\times10^{-6}$ at Wisconsin scale (check I1). The referee's "major for actual count tables" is right in principle. In practice it matters only for tiny tables, at near-ties, or (not yet quantified) with unit-level restrictions.
* The γ = 1 branch count (P-B3). It is a real omission, but it arises only when every unit is pure, so it is minor.

**Where the errors are.**
* **2c has three results that are wrong as stated:**
  * M3(d), which omits γ = 1;
  * M5(c), which fails for district rungs;
  * N2.3, which uses a derivative that is not identified.
* **2c needs a new scope for three more:**
  * M2, which is sharp only for continuous tables;
  * N2.2, which needs a support convention and conditioning on size;
  * §6, where binning and CV selection need new arguments.
* **5a needs local fixes:**
  * the face LP omits the target intercept (R2 step 4);
  * "redundant" must mean "for this target" (R2(e));
  * the "unrefutable" sentence imports M3 into the wrong base model (l.309);
  * the "a fortiori" containment (l.311) is unproved;
  * the inference statements must carry the regime and CLR conditions (R10).
* **Core documents.** The outside-scope items hit 2a (the support of the restricted measure) and P3 (projection with a singular covariance; the residual CLT). They are proof repairs. They do not change any default.

**Real exposure for users today: limited, but real.** The code turns out to be more careful than the documents:
* **Not exposed.**
  * N1 (integer or continuous) and N2 are not implemented. The M-SIGN build left them out of scope (5a §5).
  * The face test includes the target intercept (`.rxc_face_moves`, `R/rxc_sign.R:349`; check R2 in `check_B_pi_code.R`).
  * The added-last verdict is target-specific. It says "cuts the set, not this target", not "redundant".
  * Package union widths merge overlapping intervals. Only `verify_R10.py` double-counts.
  * The CLR candidates come from the augmented system, not from appended rows.
  * The conservative slab includes the target constant in the critical value.
  * The actual-gap allowance adds the per-cell envelopes. By Minkowski's inequality this is valid, and more conservative than the referee's single-contrast envelope.
* **Exposed.**
  1. **M-SIGN is labelled "proved".** `assoc` and `cr` carry `proof_status = "proved (...)"` (`R/rxc_modules.R:176, 190`). So `print()` never shows the EXPERIMENTAL flag that NEXT_BUILD and EXTERNAL_REVIEW_PLAN §4 promise until batch B is resolved.
  2. **The "never refuted" wording.** The roxygen (`R/rxc_api.R:152–154`), the emptiness note (`:250`), the module meta (`R/rxc_modules.R:174, 187`), the vignette and PLG §9 all say that with two rows a $\delta_B$ sign and CR "can never be refuted by aggregate data on their own". That holds only in the continuous, model-free table problem (P-B5, P-B14). The message is printed right after an LCE-based rejection, where it invites the reading that the belief survived.
  3. **M-SIGN under grouped rungs.** It runs under M1, M1(δ), M2 and L4.
     * The rows are the correct full-district functionals: `.rxc_mlin` uses each unit's own group's $H$, and the W1 slab identity holds with the district-fitted τ (I2).
     * But R3 and R10 are stated only for the pooled model.
     * The output still prints "sharp for the regression (PHASE5a R3(b))".
     * No test covers M-SIGN with groups.
     * The intuition the documents teach, "positive slopes imply $\delta_W\ge0$" (2c M5(c), 1b l.120, PD P-N), is false under M1 (I2).
  4. **The projection label.** It says "asymptotic guarantee" without qualification. `.rxc_projection` (`R/rxc_infer.R:232`) uses the numerical rank of $\hat V$ for the χ² degrees of freedom and gives zero standard error in null directions. That is the case P3 §5 item 2 shows can fail. With $C-1$ columns and a full-rank design, $\hat V$ is nonsingular, so the risk is limited to degenerate designs.
  5. **Point-level verdicts called "refuted".** `eiRedundancy()` labels an empty system at $\hat\beta$ "empty: refuted jointly with the base" (`R/rxc_sign.R:373`). That is a sample diagnostic, not a level-controlled rejection (minor comment 4).

**Consequences.** The 0.3.0 release does not need to wait for new theory. It needs M-SIGN flagged experimental, the wording fixed, and a guard in `.rxc_projection` (§3(b)). The build order needs three changes:
* the nesting claims that rest on M5(c) or P-N hold for pooled M0 only;
* N2's local rows (N2.3) should be dropped or redesigned before step 8;
* "exact" for N1 and Δ4 must read "exact for continuous tables".

---

## 1A. The second report (ecological inference)

**Agreement.** I agree with "major revision" here too, for the user-facing material (PLG §9, the vignette section, 5a §3.4–§3.5 and §5) rather than the theory. Every numerical claim in the report reproduces. The two reports agree on four points, so these are now confirmed twice:
* N1 is sharp only for continuous rates (P-B1 = P-E11). The EI example also shows that an **endpoint** moves, not just CR emptiness: the positive-between lower end is 0 continuously and 2/15 in integers. That is consistent with my lattice bound (step 8/15 here).
* M3(d) at γ = 1 (P-B3 = P-E13).
* The `verify_R10.py` union length (P-B38 = P-E21). The EI referee quantified it: in the near-zero design the branch union is 0.242, not 0.430, the same as the slab. The headline slab comparison (0.359 against 1.106) survives.
* N2 support, size conditioning, derivative, binning and CV selection (P-B9 to P-B14 = P-E23).

**What is new.**
1. **W2 is described wrongly to users.** PLG §9 calls it "the contextual effect on each group separately". It is one pooled sum, and W1 and W2 can have opposite signs while every group's rate rises (W1 −0.026, W2 +0.001; checked independently, I5). 1b §D.1 has the same error, and 1b l.309 ("exactly what biases Goodman's multiple regression") overstates: signing $\Delta_A$ does not sign the bias of A's coefficient (P-E1, P-E5).
2. **Reporting.** False CR gives narrow, nonempty, reproducibly wrong intervals: at 600 precincts, 0 of 400 cover with a moderate-to-large departure, while the base LCE covers 400 of 400 (P-E22). This is the expected behaviour of a false restriction, not a bug. The point is that the output does not show the price. The referee asks for the base result beside the restricted one, belief-by-belief removal, weight diagnostics, and a misspecification example next to the success example.
3. **Print loses pieces.** With several non-W1 CRs the union can have more than two pieces. The print and the `*_gap_*` columns keep only the widest gap (`.rxc_union_summary`), so the other gaps disappear silently from the printed line. `fit$pieces` has the full list (P-E15, a code finding confirmed here).
4. **Default CI for CR.** The default is the CLR branch union (algorithmic). The quoted 0.359/0.992 are the conservative slab; the package's gated CLR run gave 0.230/0.990. PLG §9 and the vignette mix the two (P-E17).
5. **Meaning of each belief.** Who is averaged, with which weights (pure precincts get zero weight in $\delta_B$); "own share" is the wrong predictor for the complement in W1; pairwise CR mixes a pair $\delta_B$ with an all-rows $\Delta_A$; turnout and denominators (P-E2 to P-E7, P-E24).

**One disagreement.** The referee writes "if a slack parameter is later implemented". It **is** implemented: `sensitivity(dial = "slack")`, per-restriction `slack` in probability units, and `test-msign-dial.R`. The referee was misled by 5a §5, which still lists "slack η on δ" as out of scope. So this is a stale-document error on our side (P-E22, P-E28), and the slack dial already meets part of the request.

**Effect on my earlier conclusions.**
* The real exposure is larger than I put it in §1, on the reporting side: the W2 wording, the hidden extra gaps, and the absence of the base result next to the restricted one.
* The case for flagging M-SIGN experimental in 0.3.0 (U3(a)) is stronger. The EI referee would not recommend the module for election analysis until the claims, examples and default reporting are revised.
* Nothing in the EI report changes the theory corrections or the build-order consequences in §3(c).

---

## 2. Ledger

**Verification codes.**
* **R:** reproduced; I agree with the claim.
* **R\*:** reproduced, but I dispute the severity (reason given).
* **D:** disputed.
* **N:** needs work. Either the referee's replacement statement or proof must still be written and checked, or the claim concerns material I could not verify here.

**Types** are as in `REVIEW_RESPONSE_PLAN.md`: **E** error, **S** scope, **P** proof, **X** exposition, **C** code.

The last column names a check that must **fail on the current statement**. Checks named `Cn`, `Tn`, `Sn`, `Vn` are proposed additions to the existing `verify_*.py` scripts.

### A. 2c: model-free sign bounds (N1), Props M2–M4 (major comment 1)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B1 | M2 "sharp in the finite-population sense" is sharp only for **continuous** unit tables. Integer counts can reject CR when the LP accepts it (3 units of size 4: the continuous CR set is {0}, the integer gaps are ±1/6). | **R\*.** Reproduced, and I1 generalises it. In 2×2 with margins only, every integer sum is attained, so the integer $D$-set is the DD interval on a lattice of step $1/N_1+1/N_0$. Integer CR is empty only when the continuous CR width $(1-\gamma)\lvert ER\rvert$ is below one step (639 of 2,987 tiny populations, all of that kind). Endpoints differ by less than one step. The step is $6\times10^{-6}$ at 370 units of 500–4,000 people. Severity is high only for tiny tables or near-ties. With unit restrictions $U_i$, or joint columns in R×C, the gap is not bounded by this argument (open). | S (E for "exact") | 2c §0.1, M2(a)(b), Reading; 1b §A2 table l.145, §D decision 1; PD §3 T1.6 "exact (M2)", §4.1 M-SIGN "sharp (M2)", M-PER "Δ0, Δ4 exact"; PLG §8 | None today: N1 is not implemented. Future N1 output must say "sharp for continuous unit tables given the margins". | Relabel M2: sharp for continuous tables conditional on the margins. Add a remark with the lattice statement for 2×2 and the referee's example. Integer N1 is decision **U1**. | **C7:** enumerate integer allocations for the referee's example. The integer CR set is empty while the M2 LP returns {0}, contradicting "attained by unit tables". |
| P-B2 | M2(d): a scalar breakdown value needs a specified conclusion, equality treatment and a scalar path. A vector of slacks gives a breakdown frontier. | **R.** | S | 2c M2(d), §6 bullet 1 | `breakdown()` uses one scalar slack per module. Document it as a scalar path. | Restate (d) as nesting componentwise in ε, with a breakdown frontier in general and a scalar value along a stated path. | **C7b:** two slacks for which the set of ε pairs preserving a conclusion is not a product of intervals. |
| P-B3 | M3(d) "exactly one of $F_\pm$ is non-empty when $D_{ER}\ne0$" fails at γ = 1 (pure units (0,0), (1,1): both branches contain the solution). | **R\*.** Reproduced. The proof notices the case (l.127) and the statement drops it. Minor, because γ = 1 means every unit is pure and $D$ is point-identified. | E (minor) | 2c M3(d); 5a R3(d) inherits ("this is the only CR") | None: `.rxc_cr_parts` builds both branches as a slab; at γ = 1 the slab is a point. | Restate as: exactly one branch when $ER\ne0$ and $\gamma<1$; both give the identified gap at γ = 1; gap 0 at $ER=0$. | **C8:** two pure units; both LP branches are feasible, contradicting "exactly one". |
| P-B4 | M3(c)'s witness $\beta_{i1}=\beta_{i2}=T_i$ can be removed by unit restrictions $U_i$ or by previously imposed signs. | **R.** M3 is stated with $F$ from Def. 2, which includes $U_i$, but the proof ignores $U_i$. | S | 2c M3 hypotheses, M3(c) | None (N1 not implemented). | Add "with no unit-level restrictions and no other sign rows" to M3. | **C8b:** a unit ordering $\beta_{i1}\ge\beta_{i2}+0.1$ that makes $\{\delta_B\le0\}$ empty. |
| P-B5 | "CR is never refuted by the data" and "Between signs and CR can never be rejected from aggregates" (2c M3(d), Reading) fail for integer tables, and are read into the wrong base model elsewhere (see P-B14). | **R.** The single between sign still has an integer witness. CR does not (P-B1). | E (wording, user-facing) | 2c M3, Reading; 5a §3.5 item 2; PLG §9; vignette "What the data can and cannot check" | **Yes.** Roxygen `R/rxc_api.R:152–154`, note `:250`, meta `R/rxc_modules.R:174, 187`, vignette, PLG §9. | Reword everywhere: "in the continuous-rate model-free problem (Duncan–Davis tables, no linear rung), a $\delta_B$ sign and CR are always feasible; with the linear rung they can be jointly rejected; with integer counts CR can be rejected outright." In the emptiness note, say whether the rung alone is rejected (already computed, `.rxc_empty_test`). | **Package test:** the emptiness note no longer contains "never be refuted" unqualified. **C7** (integer CR empty). |
| P-B6 | M4 established. | R | — | — | — | No action. Cite the referee's two-column example in the C3 discussion. | — |

### B. 2c: K-rows under districts, Prop. M5 (major comment 2)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B7 | M5 is incorrect across the advertised rungs. (a) holds per rung only with the full district parameterisation; the single $M_{Q+K}(X_i)$ formula proves only the pooled case. (b) needs district-specific ordering on each district's support. (c) fails under M1: six units in two districts with common positive slopes give pooled $\delta^*_W=-0.0493$. | **R.** Reproduced. I2 adds three things. (i) Within each district $\delta^*_W=+0.00064$. (ii) The 2×2 local rows $xm_1'+(1-x)m_0'\ge0$ hold everywhere. So PD's P-N chain "signed slopes ⊂ weighted-local ⊂ $\delta_W\ge0$" also fails for pooled M1. (iii) The W1 identity $\delta^*_W=\operatorname{Var}_p(X)(ER^*-D^*)$ still holds with the district-fitted τ, so the CR slab rows remain valid. | E | 2c M5, §5 ladder (l.221); 5a §1.1 (sign rows "2c M5(a)"), R3 setting (pooled only); 1b l.97, l.120 ("proved under LCE"); PD §0 P-N (l.39), §4.2 cells ORD⇒SIGN-B, LOC⇒SIGN-W (footnote d), N2⇒SIGN-W | **Code correct, labels not.** `.rxc_mlin` uses each unit's own group's $H$ (full parameterisation, M5(a) per rung). The `order` rows are imposed per group at each group's vertices (M5(b) per district). But M-SIGN under M1/M1(δ)/M2/L4 prints "sharp for the regression (PHASE5a R3(b))", R3/R10 are stated only for the pooled model, and no test covers M-SIGN with groups. | Restate M5: (a) per rung with the stacked district parameterisation; (b) per district; (c) for M0, or with district-conditional weights; give the between-district covariance term $\sum_g\mu_g\operatorname{Cov}_{w^g}(X,\zeta_{g,D})$ for M1. Scope P-N and the §4.2 cells to M0. New proof item **R16**: R3′/R10′ for grouped rungs via 2b Thm 1′ (decision **U2**). | **C9:** the referee's six-unit design; M5(c)'s hypothesis holds and the conclusion fails. **Package test:** with `group` given, the sign/CR print line does not claim R3(b) until R16 is proved. |
| P-B8 | Remark M5.1: "order $p^{-1/2}$" for $\delta-\delta^*$ needs sampling and weight conditions; (S2) is a mean condition only. | **R.** | E (minor) | 2c Remark M5.1 (l.160) | None (N1 and K-rows are reported separately in any case). | Replace with: mean zero given the design; its scale is the root sum of squared design coefficients under P3-S1a, with a weight condition for concentration. Or delete the rate. | — (wording) |

### C. 2c: N2 and §6 (major comments 3 and 4)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B9 | N2.2(a)–(b): "support points (the distinct observed compositions)" is ambiguous. On a subset of a larger support, extendibility subject to the other accounting equations can fail (support {1/4, 1/2, 3/4}, dropping the middle point). | **R.** Reproduced. | S | 2c N2 outer set (l.177), N2.2(b), proof step 2; PD §4.1 M-N2 "grid LP" | None today (N2 not built). Shapes step 8 (grid LP), where a grid coarser than the support is a relaxation, not "nothing lost". | Define N2 on the finite empirical design law, where all of (a)–(d) hold. For grids, state only an outer relaxation, using the binning repair (P-B13). | **C10:** the referee's three-point example. The LP on the two extreme points accepts values that no monotone extension can produce at the middle point. |
| P-B10 | N2.2(a): the proposed sharp version $E[S(X,T)\mid X=x]$ is not sharp. The model conditions on $V=(X,N)$. The condition is $m(X)\in E[S(X,T)\mid V]$ a.s., with the shape restrictions. | **R** for the counterexample (two sizes at $X=1/2$). **N** for the replacement: it must be proved via 2a Thm 1 selection for each measurable $m$. | E (proposal) | 2c N2.2(a); also bears on M-SHARP (PD step 1b, R9): the series bound conditions on $X$ only, so call it "tighter than $\Theta_O$", never sharp | M-SHARP's documentation must not claim sharpness when size is informative. | Replace the sentence with the referee's $V$-conditional statement, marked Conjecture until proved. Add the size example as a remark. Carry the caution into R9. | **C10b:** the size example; the X-conditional LP admits $(0.75,0.25)$, which the full model excludes. |
| P-B11 | N2.3: the derivative is not identified on isolated support points. Adding $h=(x-\tfrac14)(x-\tfrac12)(x-\tfrac34)$ leaves the law unchanged and moves $\tau'(1/2)$ by $-1/16$. An LP using the wrong $\tau'$ need not cover the true gap. | **R.** Reproduced. | E | 2c N2.3, C6, §5 ladder line 2, §6 ($\tau'$); 1b A2 (N2-local); PD §4.1 (N2-local), §4.2 cell LOC⇐N2; EGGHM mapping (P-B33) | None today. **Blocks the N2-local rows** of PD step 8. | Restate N2.3 for the structural $\tau_m$. As an observable restriction, require an identified derivative: an interior point of an interval support with compatible smooth versions, or a maintained identified regression family. An estimated derivative adds sampling error. Decision **U4**. | **C11:** the referee's $h$ example; two compatible models give different local thresholds. |
| P-B12 | N2.3(d)'s last assertion is true even for smooth functions ($X\sim U[1/5,2/5]$, $m_1=x-x^2/2$, $m_0=1-x^2/2$: local zero, aggregate $-1/47250$). The reverse non-implication also holds. | **R.** Reproduced ($-2.1164\times10^{-5}$). | X | 2c N2.3(d), C6; 1b V6 | — | Replace the kinked V6 witness by the smooth example and add the reverse direction. | **C13:** the smooth example (strengthens the check; nothing fails). |
| P-B13 | 2c §6 binning: "averaging Lemma N2.1 within a bin uses [CE-cvx]" is wrong, since $S(x,\tau(x))$ varies with $x$ (example: 0.3125 against 0.25). Repair: person-weighted bin means $\bar m_{1b}=E_b[Xm_1]/\bar x_b$ etc. give an outer relaxation for global monotonicity. Local rows need a separate argument. | **R** for the counterexample (0.3125, 0.25, 0.625 reproduced). **N** for the repair proof. | E | 2c §6 bullet 3; PD step 8 | Step 8 must implement the person-weighted bin rates. | Replace the binning sentence by the repair, proved as a lemma (with absent groups by zero-weight extension). | **C12:** the two-point bin example. |
| P-B14 | 2c §6: fixed-degree inference after LOO-CV selection of the degree is invalid (conditional coverage about 68%). Degree $J$ with $J$ support points is not identified; $J$ is used twice. | **R.** I3: $1-0.05/P(\lvert Z\rvert>\sqrt2)=0.682$. | E | 2c §6 bullets 1–2; 5a none | Step 8: no CV-selected degree in the default. | Prespecify the degree, or use sample splitting or a selection-aware procedure. Rename the support size (e.g. $J_X$) apart from the degree. | **C12b:** LOO selection simulation in the balanced two-level design; conditional coverage ≈ 0.68. |

### D. 5a: Theorem R3 (major comment 5)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B15 | R3(b): calling the computation "an LP" needs polyhedral $C$ (a ball is admissible under (C1)–(C2)). (C2) is stronger than (a) needs. The dimension condition is sufficient (A5 implies it). $P(X\in A)=1$ is needed for the witness. "Fixed design ⇒ $\tilde B=B^*$" is a convention, not an identity for a superpopulation rate. The full-law versus regression distinction is illustrated by $C=[0.2,0.4]$: full law empty, regression $[1/6,1/3]$. | **R.** Reproduced. | X/S | 5a §1.1 (l.22, l.24–32), R3(b) | None. All package restrictions are finite row sets. | Say "a finite LP when $C$ is polyhedral". Add the full-law versus regression example as a remark under (b). Reword l.22 as a convention. | **S6:** the referee's $C=[0.2,0.4]$ example; regression interval $[1/6,1/3]$, full-law set empty. Shows (a) and (b) differ; nothing fails. |
| P-B16 | R3(e): "never $C$ alone" needs qualifying when a belief has a model-free implication ($m_r=1.5-2x$: ER = −2, positive within association impossible model-free). 5a l.309 "rests on an unrefutable belief (2c M3)" imports M3 into the wrong base: on $X=0,1/4,1$ LCE fits, and both LCE CR branches are empty ($\tilde D=-261/1925$, ER = 0). | **R.** Reproduced. | E (wording) | 5a R3(e), §1.6 item 4, §3.5 item 2 (l.309); PLG §9; vignette | **Yes:** same strings as P-B5. The emptiness note should state whether LCE alone fits. The code computes this (`.rxc_empty_test` reports whether the rung alone is rejected) but the note does not use it. | R3(e): "an empty LP rejects the conjunction; emptiness alone does not attribute the failure". l.309: name the base model, and say that under LCE CR is refutable even when LCE fits. Package: print "LCE alone is / is not rejected" next to the joint rejection. | **S7:** the $X=0,1/4,1$ example; both LCE branches empty while LCE fits, contradicting "unrefutable" as applied in §3.5. |
| P-B17 | R3(c), (d) established. Minor: a common witness is only sufficient. In the square with $g_B=u+0.5$, $g_W=u-0.5$, target $v$ has identical branch images although the common-witness LP is infeasible. | **R.** | X | 5a §1.6 item 3 (already says "sufficient") | None (the certificate is reported as sufficient). | Add the referee's example to §1.4. | — |

### E. 5a: Proposition R2 (major comment 6)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B18 | R2 step 4 (l.169): the face LP $h^\top k\ge U$ omits the target intercept. It should be $h^\top k\ge U-\varphi_0$. On $[0,1]^2$ with target $2+x$ the printed row is $x\ge3$. | **R** for the document. **Code is correct:** `.rxc_face_moves` subtracts $c(\hat\beta)+c_0$ (`check_B_pi_code.R`: a cut through the facet keeps $U$; a cut removing it moves $U$; the document's row is infeasible). `verify_R2.py` uses zero intercept, so it cannot see the error. | E (document only) | 5a §2.3 step 4; `verify_R2.py` | None. Add a regression test with a nonzero intercept. | Fix step 4. | **T1c:** target $2+x$ on the square; the printed face LP is infeasible (so it reports "moves") while the brute-force endpoint is unchanged. **Package test:** `.rxc_face_moves` with $c_0=2$. |
| P-B19 | R2(b): a positive multiplier certifies tightening only with uniqueness or a positive lower bound on all optimal multipliers (duplicate-row example). R2(c): the separating witness need not be a reported cell-rate target. | **R.** | X | 5a §2.2(b), (c), §2.5 | None (the package never uses duals for verdicts). | Add the duplicate-row example to (b). Add "the witness objective may not be a reported target" to (c). | — |
| P-B20 | R2(d) holds under γ < 1, a nonconstant gap target, and nonempty intersections. §2.5 item 2 ("in 2×2 the middle verdict cannot occur") needs the nonconstant-target qualification. | **R.** | S | 5a §2.2(d), §2.5 l.208 | Package test "middle verdict never occurs in 2×2" uses gap targets. Check it also holds for $\tilde B_{rc}$ targets (nonconstant in $w_1$ unless a degenerate design). | Add the qualification. | — |
| P-B21 | R2(e): "redundant given the others" is incorrect if it means implied. Square, target $x$, rows $x+y\le0.8$ and $x\le0.5$: zero last-added gain for $x$, yet the row lowers the upper endpoint of $x+y$ from 1.5 to 0.8. | **R** for the document. **Code is correct:** the added-last verdict is "cuts the set, not this target" (checked on the same example, `check_B_pi_code.R`). | E (wording) | 5a §2.2(e) l.161, §2.5 item 3; PLG §9 | None in verdict strings. Column names `gain_last` are fine. | Reword (e): "zero incremental gain for this target"; "redundant for the feasible set" only via the implication test (c). | **T4b:** the square example; the old sentence calls the row redundant while it tightens $x+y$. |
| P-B22 | §2.5 item 4: width-gain inference needs joint treatment of the base and augmented endpoints (and any estimated threshold). An interval for one endpoint is insufficient. | **R.** | S | 5a §2.5 item 4 | The open item "width change with its CI" (5a §5) must be built jointly. The current `near_tie` uses a one-sided CI on the association margin: a different, valid diagnostic. | Restate item 4 and the open item. | — (design item) |

### F. 5a: Theorem R10 (major comment 7)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B23 | R10(a) holds per version only under that version's conditions for the **augmented** system. Growing-district M1/M2 (regime R3) is conjectural (P3 §0.5) independently of candidate collection. "Finite vertex set at each $n$" does not verify CLR Condition P (fixed compact index set, gradients bounded away from zero, well-conditioned covariance, rate). A zero-variance **candidate** is different from an exactly known row. Appending row standard errors to the old CLR maximum is insufficient. | **R.** The package already collects candidates from the augmented system (`.rxc_collect(cons, …)` uses the full `cons`). Zero-s.e. candidates are zeroed in the CLR maximum (`W[, se == 0] <- 0`) and kept as exact bounds. That is right only if their gradient is exactly zero. | S | 5a §3.1 (l.218 "R1/R3"), §3.2 status note, step 4; P3 §0.5; PD §8 R10 | Labels: the CLR branch union (default for CR) is algorithmic in every regime. The conservative slab and projection are proved for M0, and for M1/M2 with few large groups. Add a candidate-gradient check: treat a candidate with $\hat\sigma=0$ as exact only when its gradient is exactly 0. | Replace "R1/R3" by the P3 §0.5 table. List Condition P's actual requirements in the status note. Add the zero-variance-candidate sentence. | — (status) **Package test:** a candidate with nonzero gradient but numerically zero $\hat\sigma$ is not treated as exact. |
| P-B24 | R10(a) for actual targets: carry over P3-S1a, the weight condition and the residual allowance. An actual gap needs its own residual envelope, e.g. $\tfrac12\{\sum_i(w_{i1}+w_{i0})^2\}^{1/2}$. A single-row $S_1$ need not bound the contrast. | **R.** **Code is valid:** `ci_actual` adds the per-cell envelopes $\sum_{rc}\lvert a_{rc}\rvert\,\mathrm{sd}_{rc}$ (`R/rxc_infer.R:388`). By Minkowski this dominates the referee's contrast envelope, so it is valid and more conservative. | P | 5a R10(a) statement and proof step 3 | None for validity. Optional: use the contrast envelope for gap targets (narrower). | State the allowance for linear contrasts (sum of per-cell envelopes, or the joint envelope). | — |
| P-B25 | R10(b), (c), (d) established under conditions. One simultaneous event includes the target constant (no strict feasibility needed). Branch union: point coverage only, no simultaneous coverage of both branch images, no selection on $\widehat{ER}$. Under triangular arrays, branch validity is needed along each subsequence. Emptiness bound $\alpha_U+\alpha_L$. Pointwise only. | **R.** The code's conservative slab includes the target row (`Gall <- rbind(cons$S, DS, tg$cvec)`, two-sided via `abs`). | S/P | 5a §3.2(b)–(d), §3.3 steps 5–17 | None. | Adopt the referee's version table (R10 against P3.5(a), (b1), (b2)(i), (b2)(ii), (b3)) as §3.2's status note. Add the triangular-array sentence to (c). | — |
| P-B26 | R10(e): an analytic counterexample replaces the Monte Carlo sketch. A repeated three-point design gives $\limsup P\le\Phi(\tfrac7{12}z_{.995})=0.9335$; with 4,000 critical-value draws, ≤ 0.9401. The supplied far-edge design also satisfies the condition (ratio 0.564, ≤ 0.927). | **R.** Reproduced. The derivation of the 0.564 ratio for the far-edge design is **N** (not rerun here). | P | 5a R10(e), §3.4; vignette ("about 0.90 … does not improve with $n$") | None in code (the naive bracket is never imposed). The vignette and PLG §9 numbers can cite the analytic bound. | Replace the sketch by the fixed-design proof. Keep the MC as illustration. | **V9:** the three-point design; assert the analytic bound is < 0.95. (Nothing fails; it completes a proof.) |
| P-B27 | 5a l.311: "the naive interval contains bracket ∩ JKST CI, so under-coverage transfers a fortiori" is unproved, and in general false (an empty pre-expansion intersection need not stay empty after expansion). Direct argument: the narrower procedure has limiting coverage ≤ 1/2 at the edge. Also distinguish the two sample ERs when weights differ. | **R.** | E (minor) | 5a §3.5 item 4 | None. The descriptive print already warns. | Replace the a-fortiori sentence with the direct ≤ 1/2 argument and the counterexample $[0,0.1]$, $[0.15,0.2]$, allowance 0.1. | **V10:** the referee's three intervals; the claimed containment fails. |

### G. Minor comments and the EGGHM mapping (report §3)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B28 | W1 versus W2: invertibility of $M$ is needed for $g$ only; ER needs only positive regressor variance; choose a bounded value for the absent rest group at $X_r=1$; W2 signs concern components of a sum and need not sign every summand. | **R** (M1 boundary identities reproduced). | X | 2c Lemma M1, Def. 1 | None. | Edit the hypotheses of M1(b)–(c); add the W2 sentence. | — |
| P-B29 | Ladder (2c §5): "both signs" ⇒ CR only when the directions agree; N2 admits steps where a local derivative is undefined; compare realised and model-implied associations only on the same population; for $R\ge3$ specify a simplex chart. The symmetric Dirichlet example gives negative covariance under both W1 weights. | **R.** | X/S | 2c §5, N2.2(e), C5(c) | None. | Fix the ladder lines. Upgrade C5(c) to the weighted covariances. | **C5(c′):** the weighted covariances under both W1 weights are negative (strengthens the check). |
| P-B30 | Identification versus sample labels: R3(e) is a population statement; raw LP infeasibility at $\hat\beta$ is a sample diagnostic; level-controlled rejection uses R10's rule. R2 §2.5 l.208 and l.210 need their qualifications (P-B20, P-B22). | **R.** | X/C | 5a §1.6 item 4, §2.5 | **Yes, small:** `eiRedundancy()` labels an empty system at $\hat\beta$ "empty: refuted jointly with the base" (`R/rxc_sign.R:373`). | Rename to "empty at the point estimate (see the emptiness test)". | **Package test:** verdict strings contain no "refuted" for point-level emptiness. |
| P-B31 | EGGHM §2, Lemma 1, Props 1, 4, 6, 9–10: correct. EGGHM assumes a mixed neighbourhood (γ < 1); ER also needs positive composition variance; the strict branch distinction at nonzero ER needs a mixed unit. | **R.** | S | 2c §0, M3; 1b §0 | None. | Add the mixed-unit condition where EGGHM is cited for branch uniqueness. | — |
| P-B32 | EGGHM Props 2–3, 5, 7, 8, Cor. 3: sharpness ranges over compatible joint distributions; it does not license integer-table interpolation. sign(ER − NM) = sign(ER) uses γ < 1. | **R.** | S | 2c M2, M3; 1b §A2 | None. | Cite this with P-B1. | — |
| P-B33 | EGGHM §4 end, Theorems 1–3: 1b l.69 attributes all three tables to Theorem 1 (the gap is Theorem 1; the group levels are Theorems 2 and 3). W2 and joint columns need their own proofs (M1, M4). | **R.** 1b l.69 checked; EGGHM theorem numbers taken from the report. | X | 1b l.69 | None. | Fix the attribution. | — |
| P-B34 | EGGHM §5, Props 11–13 and 14–16: the neighbourhood readings use Assumption 1 and meaningful derivatives. N2 imposes a common conditional mean. Neither identifies a derivative on an atom. At a pure composition EGGHM's local between association is 0; replacing it by the sign of $D(x)$ changes the belief. | **R.** | S | 2c N2.3; 1b §A2 "EGGHM local results" (l.162) | None. | With P-B11. | — |

### H. The referee's notes on our verification scripts (report §4)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B35 | `verify_phase2c.py` C4's ordering construction changes one column without rebalancing the others, so it does not always build a valid full transition matrix. No check covers integer feasibility, district intercepts or derivative identification. | **N** (C4 not re-read here). | C (script) | `verify_phase2c.py` C4 | None. | Rebuild C4 with valid transition matrices. Add C7–C13 above. | C7–C13 |
| P-B36 | `verify_R3.py` S2 labels a branch nonempty only if its width exceeds $10^{-9}$ and calls zero positive-width branches "both degenerate" (possibly two empty sets). Frequencies over reused designs are descriptive. | **R.** Code read: `verify_R3.py` l.102. | C (script)/X | `verify_R3.py` S2; 5a §1.5 | None. | Separate empty, degenerate and positive-width branches. Call the frequencies descriptive. | — |
| P-B37 | `verify_R2.py` T4 now gives 153/494/1,153 of 1,800, against 147/494/1,153 of 1,794 saved. Unsuccessful base solves are skipped without recording status. The face LP in the tests uses zero intercept. | **N** (not rerun; the earlier session reproduced the committed output under the old runner). | C (script) | `verify_R2.py` T4; 5a §2.4 | None. | Record solver status for skipped cases; regenerate the output; add T1c (P-B18). | T1c |
| P-B38 | `verify_R10.py` l.189–190 sum the two branch lengths even when they overlap, so "union length" can exceed the Lebesgue length. The runs do not test CLR, residual coverage, the crossing rule or growing dimension. | **R** (code read). The package merges overlaps (`.rxc_union`). | C (script)/X | `verify_R10.py`; 5a §3.4 table ("CR union" 0.386) | None. | Merge before measuring; regenerate; correct the 5a §3.4 figure. State what the MC does not test. | **V11:** overlapping branch intervals; the script's length exceeds the merged length. |

### I. Outside scope: core documents (report §5)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B39 | 2a tag [supp] (l.57), Lemma 1 (l.144–150), Lemma 3(c) (l.310–314): continuous functions equal a.e. **on** $A$ agree on the support of the measure restricted to $A$, which can be smaller than $\operatorname{supp}P_X\cap A$ ($X\sim U[0,1]$, $A=\{1/2\}$). This also explains the tension between Lemma 3(c) and Remark 5.1. | **R.** Elementary counterexample. Note that Remark 5.2's argument (A5 ⇒ full-dimensional) also uses $\operatorname{supp}P_X\cap A$. | E (core) | 2a §0.3 [supp], §0.2 (A5 wording), Lemma 1(a) step 7, (b) step 3, Lemma 3(c) step 8, Remarks 5.1–5.2; P1 A5; 2b Lemma D1 if it uses the same tag | None. The package works with the finite empirical design, where the restricted support is the set of observed units in $A$; `range = "trim"` already uses those. | Define $\operatorname{supp}(P_X\vert_A)$ and impose A5 on it. Rewrite [supp]. Check every use. Report it to referee A with the batch-A material (it corrects a reviewed core result). | **2a check:** $X\sim U[0,1]$, $A=\{1/2\}$; the old [supp] claims $0\equiv1$ at $1/2$. |
| P-B40 | P3.5(b2)(ii) (l.179–189): the rank-$q$ pseudoinverse proof bounds gradients in the range of $\Sigma$ only. A singular CLT need not control smaller errors in null directions, e.g. error $(Z/\sqrt n,-1/n)$. | **R** (I4). | E (core) | P3.5(b2)(ii), step 4; 5a §3.2 status note; PD §4.1 M-LCE | **Yes, limited:** `.rxc_projection` (`R/rxc_infer.R:232`) sets $q$ to the numerical rank of $\hat V$ and gives zero s.e. to null-direction candidates. With $C-1$ columns and a full-rank design $\hat V$ is nonsingular, so only degenerate designs are exposed. The print line "projection intervals (asymptotic guarantee)" is unqualified. | Document: require $\Sigma$ positive definite (as (b1) already does), or work in a full-rank reduced parameterisation, or require exact null-coordinate identities. Code: if $q<\dim\beta$, use $q=\dim\beta$ and warn, or set the projection to NA with a message. Add a test. | **P3 check:** I4's toy; zero-s.e. bound misses always. **Package test:** a rank-deficient $\hat V$ triggers the guard. |
| P-B41 | P3 [Lyap] (l.40) and the residual step (l.201–203): $\max_iw_i^2/\sum_iw_i^2\to0$ does not give a variance-standardised CLT if only one summand has positive variance. The tail bound with the envelope $S_1$ still holds: Chebyshev where $\mathrm{sd}(E)/S_1\to0$, Lindeberg where the ratio is bounded away from zero. | **R** for the gap. **N** for writing the argument. | P (core) | P3 §0.3 [Lyap], P3.5(a) step 5, P3.6 | None. The allowance is unchanged. | Replace step 5 by the subsequence argument. Report to referee A with batch A. | — (proof completion) |
| P-B42 | 1b §§A1–A2 (l.59–80, 116–121): fix the theorem-to-table mapping (P-B33). For pure units, unit-level row ordering cannot be equated with the observable implication of a local between sign (the between covariance is identically zero). | **R.** | X | 1b §§A1–A2 | None. | Fix the mapping. Add the pure-unit caveat. | — |

### J. The referee's comparison with the earlier internal review (report §6)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-B43 | The internal review (5a §4) covered 5a §§1–3 only and missed the face-LP intercept, set-versus-target redundancy, the empty LCE–CR example and the union-length double count. Its R10(e) sketch left verification to simulation. Combining two tests by a union bound gives an upper bound; the rejection probability need not double. Inference labels should keep P3's growing-district qualification separate from the candidate-simulation conjecture. | **R.** | X (process) | 5a §4, §3.1 (l.224, "would double the bound"); MSIGN_REVIEW_LOG; EXTERNAL_REVIEW_PLAN §3 | None. | Reword the §3.1 parenthesis ("would double the **bound**", not the level). Add to 5a §4 a line pointing to this ledger. Process lesson for EXTERNAL_REVIEW_PLAN §3: the fresh-session pre-screen must also cover the scope document (here 2c), not only the new proofs. | — |

### Part II. Ecological-inference report (`B_ecological_inference_report.md`)

Same columns and codes. "= P-Bn" marks a point the first referee also made; the action is shared and not repeated.

#### K. Meaning of the beliefs and weights (major comments 1–2)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-E1 | W2 is described as "the contextual effect on each group separately" (PLG §9 l.214; 1b §D.1). It signs one pooled sum $\Delta_A=\sum_s\mu_s\operatorname{Cov}_{w^s}(X_A,m_s)$; summands can have either sign. W1 and W2 can disagree: two units in which every group's rate rises 2 points give W1 −0.02633, W2 +0.00110. | **R** (I5, recomputed independently). Same as minor point P-B28 (W2 signs a sum), with a user-facing consequence. | E (wording) | PLG §9; vignette selector text; 1b §D.1 (l.305–311), §A3; 2c Def. 1 | **Yes, wording.** Label W2 "average within the original groups (pooled over groups)" in roxygen and the `assoc`/`cr` labels. A per-group direction needs separate restrictions, which the package cannot express today (decision **U8**). | Fix the wording. Put the two-unit example next to the W1/W2 selector in the vignette. | **C14:** the two-unit table; old wording "each group" is contradicted by W2 > 0 with W1 < 0. |
| P-E2 | W1: the predictor in **both** covariances is A's share. The vignette's "rises with the group's own share" is the wrong predictor for the complement. Either covariance can be negative when the sum is positive. | **R.** | E (wording) | vignette l.155–156; PLG §9 | Vignette and roxygen text only. | Reword: "rates of A, and of everyone else, rise on average with A's share". | — |
| P-E3 | $\delta_B$ (one-vs-rest) is a weighted average with weights $p_ix_i(1-x_i)$: pure precincts get zero weight; it is not an ordering in every precinct or against each rival; a large mixed precinct can dominate. | **R** (Def. 1). | X | PLG §9; vignette; 2c §1 | Diagnostics: report the largest weight and the effective number of units for each signed association (P-E25). | Add the weighting sentence wherever $\delta_B$ is explained. | — |
| P-E4 | Within covariances also reflect sorting, campaign targeting and other between-precinct differences; a causal reading needs a separate design. | **R.** | X | PLG §9; vignette; 1b §A1 | — | Add one sentence. | — |
| P-E5 | $g-B=M^{-1}\Delta$: signing $\Delta_A$ generally does not sign the bias of A's Goodman coefficient. | **R** (2c M1(a)). | E (wording, background) | 1b l.309 ("exactly what biases Goodman's multiple regression") | — | Reword 1b l.309: $\Delta$ enters the bias through $M^{-1}$. | — |
| P-E6 | Pairwise CR (5a l.94) pairs a pair $\delta_B(A,B)$ with $\Delta_A$, which involves all rows. An A-versus-B argument supports the first, not the second. Print both terms separately with every CR result. | **R.** | S/C | 5a §1.4; PLG §9; vignette | **Yes, small:** the CR print should show the two terms (point interval and CI under the rung) next to the CR result. The "already implied" table has them but is printed separately. | Add the terms to the CR table and the print. | **Package test:** a pair-CR fit prints both terms. |
| P-E7 | Applied credibility: an argument for a between sign rarely supplies one for CR (racial bloc voting, partisan transitions with abstention, veterans). Pre-specify each restriction (row, column, population, geography, pairing, reason); show results without the belief and under the alternative. Under district regressions the sign is **one pooled restriction**: not per district, and cross-district composition differences remain. | **R.** The pooled-restriction point = P-B7. | X/C | PLG §9; vignette; 5a §5 l.353 | **Yes:** print the base (rung-only) interval next to the restricted one by default (decision **U7**); say "one pooled restriction across groups" in the grouped-rung note (P-B7 label). | Add a reporting checklist to the vignette and PLG §9. | — |

#### L. Targets, sharpness and refutation (major comments 3–4)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-E8 | EGGHM mapping faithful for a bounded column indicator; W2, pair weighting and CR with those pairings are extensions; EGGHM's mechanism is a reason to believe a sign, not a causal claim; local Props 12–13 versus 14–16. | **R.** = P-B31 to P-B34. | X | 2c §0, 1b §0 | — | Mark W2/pair/their CR as extensions of EGGHM in 2c and PLG. | — |
| P-E9 | JKST plus a sign is a coherent intersection; a sign already implied by the base adds nothing; inferring the realised sign from the model sign, demanding it per district, or reading it as an order at every composition adds assumptions. The base also assumes no size dependence (2a A3); weighting does not repair a violation. | **R.** | X | PLG §9; vignette | The "already implied" table already shows the first point. | Add the list of "what a sign does not mean" to PLG §9. | — |
| P-E10 | R3 delivers an *estimated bound for the model-implied gap using the regression and design*; the full-law set can be smaller ($T\in\{0,1\}$ example); report whether the CI targets $\tilde D$ or $D$. | **R.** = P-B15 (same example). | X | PLG §9 l.215 ("exact for what the regression knows"); vignette | **Yes, wording:** PLG and vignette drop "exact"; print already labels `ci_model`/`ci_actual`. | Reword as the referee suggests. | — |
| P-E11 | N1 integer: two four-voter precincts (1 and 2 A voters, one P vote each). Continuous CR set {0}; integer gaps −2/5, 2/15, 2/3, so CR is empty, and under the positive between sign the integer lower end is 2/15, not 0. | **R** (I6). = P-B1. Adds endpoint failure, consistent with the lattice bound (step 8/15). | S | as P-B1 | as P-B1 | as P-B1 (U1). Cite both examples. | **C7** with this example as a second case. |
| P-E12 | Model-free non-refutability: the equal-rates witness ($\beta_{ir}=T_i$ for all $r$) works for **all rows and columns**, so a weak $\delta_B$ sign and CR are model-free feasible in R×C too. But known cells, zeros, integer counts or other restrictions can remove the witness, so "only with straight lines" is too narrow. An empty plug-in, an empty population set and a test rejection are different events. Show the base-only result beside the joint result. | **R.** The R×C extension is new and correct (all pairwise $\delta_B=0$ at the witness, so both weak CR branches hold). = P-B4, P-B5, P-B16, P-B30. | E (wording) / S | 2c M3(c), Reading; 5a R3(e); PLG §9; vignette; package strings | as P-B5, P-B16, P-B30 | Extend M3(c) to R×C as a remark (with "no unit restrictions"). Use the referee's three-way distinction in the emptiness note. | **C8b** (unit restriction removes the witness). |
| P-E13 | M3(d) at γ = 1; with $ER=0$ the singleton is the gap's image, not a single latent table. | **R.** = P-B3. | E (minor) | 2c M3(d) | — | as P-B3, plus the image sentence. | **C8** |

#### M. Reporting, uncertainty and the simulation (major comments 5–6)

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-E14 | Vignette l.163 writes the W1 slab as $[\gamma_rER^r, ER^r]$, reversed for negative ER (voter_vax, 5a l.373). The single-LP result does not license choosing a branch from the sign of $\widehat{ER}$. | **R** for the text. **Code correct:** `.rxc_descriptive` uses min/max, and the slab is a hull. | X | vignette l.163; 5a §1.2(d) already uses "hull" | Vignette text only. | Write $\operatorname{hull}\{\gamma_rER^r,ER^r\}$. | — |
| P-E15 | One W2/pair CR can give two pieces, one or none; a failed certificate leaves overlap undecided; several CRs can give more than two pieces. Say what the summary columns keep when the piece list is longer. | **R**, and confirmed in code: `.rxc_union_summary` keeps the hull, the piece count and only the **widest** gap; `print.eiRxC` prints `[a, b] U [c, d]` from those, so with three or more pieces the other gaps are not shown. `fit$pieces` keeps all pieces. | C | 5a §5 table; vignette; PLG §9 | **Yes:** print every piece when `n_pieces > 2` (or print "k pieces; see fit$pieces"); document the columns. | Fix print; add a test with two W2 CRs producing three pieces (or a mock union). | **Package test:** a three-piece union prints all three pieces. |
| P-E16 | Print the estimated set and the confidence set on separate lines; keep every component, including point components; state units; values in an estimated gap are excluded by the assumptions, but the CI may still cover them; branches have no probability weights; a hull includes the gaps; marginal endpoints of different cells need not be jointly attainable. | **R.** | C/X | print, vignette, PLG §9 | **Yes, print layout** (with P-E15). | Add the interpretive sentences to the print footer and vignette. | — |
| P-E17 | The default CR CI (CLR branch union) is conjectural even under a correct model; say so on the result line; show the conservative result by default or alongside. The 0.359 width and 0.992 coverage (5a l.300) are the conservative method; the default-CLR experiment (l.377) gave different numbers. Uniformity open. | **R.** PLG §9 quotes "0.36 … 99%" next to "the default interval is algorithmic", which mixes the two. | X/C | PLG §9; vignette ("coverage 0.99"); 5a §3.4, §5 | **Yes:** print labels CLR as algorithmic in the footer only. For CR fits, show the conservative slab (or projection) beside the CLR union (decision **U6**). | Correct the quoted numbers (method named with each). | — |
| P-E18 | Descriptive values: name them "weighted aggregate regression slope" and "neighbourhood gap", name the focal group, give γ, say fitted τ or observed $T$; equality with EGGHM holds only under size-weighted LS; voter_vax is an arithmetic check; the 8–10% miss rates are design-specific. | **R.** | X/C | `print.eiRedundancy`; vignette; PLG §9 ("8–10% … however large the sample") | **Yes, labels** in the descriptive table header. | Rename columns or add a header line; reword PLG's "8–10%" as design-specific, citing the analytic bound of P-B26 for persistence. | — |
| P-E19 | R2 verdicts need a nonempty convex base and a feasible added restriction; label each with its base, target and order; zero last-added gain ≠ logical redundancy; with unions compare membership, not hull width (removing the middle of [0,1] keeps the endpoints); report occupied length and hull span; label the width-only fallback (5a l.357) inconclusive about implication; near ties, mark verdicts as estimated. | **R.** = P-B21, P-B22, plus the union case. The code's union fallback compares occupied length (`.rxc_union_len`), so it sees the [0,1] example; but its `method` label does not say "inconclusive about implication". | C/X | 5a §2.5, §5 | **Yes, small:** add hull span next to the occupied-length gain; relabel the fallback method; add base and order to the verdict print header. | as stated | **Package test:** verdict print names base and order. |
| P-E20 | Report table: DD for the realised target, unrestricted LCE, LCE plus the belief; target per column; sign direction, comparator, weighting, rejection flags; number of units, composition range, size concentration. For point-identified rungs that ignore an entered sign (5a l.363), check whether the fitted point violates the sign or mark the combination unevaluated. | **R.** Code: `print.eiRxC` shows DD, estimate and CIs, not the unrestricted LCE interval. L5/L6 ignore signs with a note (`R/rxc_api.R:254`) and do not check the fitted point. | C | print; `ladder()` | **Yes:** base LCE column (U7); evaluate the sign at the point for L5/L6 and report "violates"/"satisfies". | Implement in WS-B9. | **Package test:** an L6 fit with a sign its point violates reports the violation. |
| P-E21 | `verify_R10.py` sums branch lengths before merging: near-zero design 0.430 → 0.242 (591/600 runs overlap); edge 0.386 → 0.357; the slab comparison survives. | **R** (reproduced exactly). = P-B38. | C (script)/X | `verify_R10.py`; 5a §3.4 table | None (package merges). | as P-B38; replace the 5a table column with the referee's three statistics (occupied length, hull span, unmerged sum). | **V11** |
| P-E22 | The design is narrow (fixed composition and sizes, Beta rates, cancelling slopes). In a matched-regression experiment where only the truth of CR changes, the CR slab covers 398, 110 and 0 of 400 at 600 precincts ($s=0, .08, .25$), always nonempty, while base LCE covers 400. Show a misspecification comparison beside the success example; expose one-belief-at-a-time removal now; calibrate any slack dial in interpretable units. | **R** (reproduced exactly). **D** on one point: a slack dial already exists in probability units (`sensitivity(dial = "slack")`, `test-msign-dial.R`); 5a §5 wrongly lists it as out of scope, which misled the referee. | X/C | 5a §3.4, §5 (stale "slack" line); vignette; PLG §9 | **Yes:** a drop-one display (each restriction removed in turn; `eiRedundancy` added-last gains are close, but not a display of intervals) (U7). Vignette: add the false-CR example with its numbers. | Fix 5a §5. Add `verify_R10.py` V12: the false-CR design. | **V12:** false CR at 600 units; coverage far below 0.95 while every interval is nonempty (documents the price; nothing in the theory fails). |
| P-E23 | N2: support convention, size conditioning, derivative non-identification (adding $a\prod_j(x-x_j)$), binning with row-specific weights, post-selection after CV. | **R.** = P-B9 to P-B14. | E/S | 2c §§4, 6 | as P-B9 to P-B14 | as P-B9 to P-B14 | C10–C12 |

#### N. Minor comments, outside scope and process

| ID | Referee claim | Verification | Type | Documents | Package impact | Proposed action | Check that fails now |
|---|---|---|---|---|---|---|---|
| P-E24 | A hand-worked turnout example: P ballots = turnout × support; loyalty among participants can coexist with a negative eligible-population gap; abstention column; who enters the denominator; new entrants and deaths. | **R.** | X | vignette; PLG §9 | — | Add the example. | — |
| P-E25 | Report the largest precinct weight and the effective number $1/\sum_i(w^r_i)^2$ per row, total and largest overlap weight for between comparisons, and a leave-largest-unit-out display. The $p^{-1/2}$ residual claim needs conditions. | **R.** Residual claim = P-B8. | C/X | 2c Remark M5.1; print | **Yes (new diagnostic):** weight concentration in `eiRedundancy()` (cheap: weights are design-only) (U7). Leave-largest-out later. | Implement the weight table in WS-B9. | — |
| P-E26 | Position M-SIGN beside King's EI and hierarchical R×C; a posterior mean inside a sign bound is no independent support. Benchmark on aggregated individual panels or linked surveys, with sign choices locked on other units. | **R.** | X | vignette; PD step S11 (real data) | Benchmark design for S11. | Add a paragraph; add the held-out benchmark to S11. | — |
| P-E27 | 2×3 versus 3×2; restrictions on all columns must respect adding-up; signs on complementary binary outcomes reverse; "both signs" in the ladder means both weakly positive or both weakly negative. | **R.** Ladder = P-B29. | X | 2c §5; PLG §9 | Possible input check: contradictory signs on complementary columns (C = 2) → error. | Add the sentences; optional input check. | — |
| P-E28 | One status key: N1, N2, local monotonicity and sign slack are outside the implementation (5a l.385); the planned API (2c l.7) differs from the shipped `restrictions` syntax. | **R** for the key and API; **D** for slack, which is implemented (see P-E22). | X | 2c header; 5a §5 | — | One status key in 5a §5 and PLG §9; fix the slack line; update 2c's API sentence. | — |
| P-E29 | 1b §A2 l.91–98: "CR is nonconvex and needs two LPs" is superseded by R3(d) (W1 slab); "different moments" does not prove non-nesting. | **R.** | X | 1b §A2 | — | Mark superseded. | — |
| P-E30 | 1b §D.1: W2 per-group description (see P-E1); §A1 l.66–69 attribution (= P-B33); l.84–87 vaccination evidence caveats. | **R.** | X | 1b | — | Fix. | — |
| P-E31 | PLG §2: A3 and A5 are not shared requirements; DD needs neither; N1 needs no quadratic regression. State denominators and rank requirements per rung. | **R.** | X | PLG §2 l.23–29 | — | Fix. | — |
| P-E32 | 1b §D.2: a reversal in the aggregate yes-share is not evidence that a group's conditional mean reverses. | **R.** | X | 1b §D.2 | — | Fix. | — |
| P-E33 | Agrees with the internal review on slab versus CLR and conditional coverage; conditioning protects nothing against choosing a belief because it gives an attractive result. Differs from "accept with fixes". | **R.** = P-B43. | X (process) | 5a §4 | — | Add the belief-selection sentence to R10(d). | — |


**Totals, part I.** There are 43 rows, P-B1 to P-B43. Two are "established, no action" (P-B6, P-B17), which leaves 41 substantive points:
* **reproduced:** 39. Of these:
  * 37 are agreed in full;
  * 2 have their severity disputed (P-B1, P-B3);
  * 4 also contain a replacement statement or proof that is not yet written or checked (P-B10, P-B13, P-B26, P-B41).
* **disputed outright:** 0;
* **not verified here (open):** 2 (P-B35, P-B37: script findings not rerun).

**Totals, part II.** 33 rows, P-E1 to P-E33. All 33 reproduce or are agreed. 13 duplicate a part-I point (P-E7 partly, P-E8, P-E10 to P-E13, P-E19, P-E21, P-E23, P-E25 partly, P-E27 partly, P-E30 partly, P-E33). One point is disputed in part: the referee's "if a slack parameter is later implemented" (P-E22, P-E28), because the slack dial exists; the fault is our stale 5a §5. New package findings confirmed in code: P-E1 (W2 wording), P-E6, P-E15 (print hides extra gaps), P-E17, P-E19, P-E20 (no base column; L5/L6 do not check the sign), P-E25.

---

## 3. Implications

### (a) The M-SIGN module

**Keep, as is (the code was checked against the referee's toy problems):**
* sign rows and the W1/2×2 slab for M0;
* the conservative slab CI;
* the CLR branch union (already labelled algorithmic);
* the face-test redundancy diagnostic and its three verdicts, added first and added last;
* the CR table and its certificate;
* the descriptive ER/NM with the "never intersect" warning;
* the union printing.

**Relabel now (no change to numbers):**
* `proof_status` of `assoc` and `cr` becomes `"experimental (referee report B under revision; see docs/extension/reviews/B_PI_RESPONSE_PLAN.md)"`. The existing print code then shows the EXPERIMENTAL line, as NEXT_BUILD and EXTERNAL_REVIEW_PLAN §4 intended.
* The "never refuted" strings (P-B5, P-B16):
  * name the base model;
  * add "with integer counts, CR can be refuted";
  * print whether LCE alone is rejected next to a joint rejection.
* "sharp for the regression (PHASE5a R3(b))" in `print.eiRxC`: keep it for M0. For grouped rungs say "outer for the regression; sharpness under district rungs not yet proved" until R16.
* `eiRedundancy()`: "empty: refuted jointly with the base" becomes "empty at the point estimate (see the emptiness test)".
* Vignette and PLG §9: the same wording, plus "redundant" meaning "for this target".
* From the EI report:
  * W2 is "the average within the original groups, pooled", never "each group separately" (P-E1);
  * W1 uses A's share as the predictor for both A and the rest (P-E2);
  * the slab is $\operatorname{hull}\{\gamma_rER^r,ER^r\}$ (P-E14);
  * PLG §9 drops "exact" (P-E10);
  * quoted coverage and width numbers name their method (P-E17);
  * the 8–10% miss rate is design-specific (P-E18);
  * under groups, the note says "one pooled restriction across groups" (P-E7).

**Disable?** I do not recommend disabling anything. The one candidate is M-SIGN under grouped rungs: its rows are right, but its sharpness and inference statements are unproved and untested. Decision **U2**. My recommendation is to keep it, flagged experimental, with a note at fit time.

**Change (code), small:**
* a regression test for the face test with a target intercept;
* tests of M-SIGN under M1 (row forms equal the direct district computation; the W1 identity with district τ);
* the zero-variance-candidate check (P-B23);
* the projection rank guard (P-B40, core code that M-SIGN also calls).

**Change (code), reporting, from the EI report.** These are output changes; they add information and change no estimate.
* Print every piece of a union with more than two pieces (P-E15). This is the only point where the current print can mislead silently.
* Print the base (rung-only) estimate and CI next to a restricted result (P-E7, P-E20; U7).
* For CR fits, show the conservative slab or the projection interval next to the CLR union (P-E17; U6).
* Print the two terms of a pair or W2 CR with the CR result (P-E6).
* Weight concentration for each signed association: largest weight and effective number of units (P-E25).
* For L5/L6, evaluate an entered sign at the fitted point and report whether it holds (P-E20).
* `eiRedundancy()`: name base, target and order; give hull span with the occupied-length gain; label the union fallback "inconclusive about implication" (P-E19).
* Later: a drop-one display (each restriction removed in turn, intervals shown) (P-E22).

### (b) The package overall

* **Core 2a/P3 issues (P-B39 to P-B41)** are proof repairs. None changes an estimator or a default.
  * [supp] does not touch the code, which works on the empirical design.
  * The residual CLT leaves the allowance unchanged.
  * The singular-covariance point needs a guard in `.rxc_projection` and a qualified label ("asymptotic guarantee for fixed-dimension designs with nonsingular $\hat V$: M0, or few large groups").
* **These corrections touch documents that referee A reviewed.** They should go into the batch-A response (WS7) as new corrections, so that referee A sees them. That keeps the two review tracks consistent.
* **Inference status labels.** Keep D4 (CLR default, algorithmic; projection alongside). Qualify "projection: proved" with the regime (P3 §0.5) and nonsingularity. For CR in particular:
  * the default (CLR branch union) is algorithmic in every regime;
  * the conservative slab is proved for M0, and for grouped rungs with few large groups once R16 is written.
* **The 0.3.0 plan in NEXT_BUILD.** The order (WS6 → WS7 → release hygiene) stands. Release hygiene gains gated items: **WS-B0** (relabel and rewording, now including the EI wording), the projection guard from WS-B3, and from WS-B9 at least the print of all union pieces (P-E15) and the base column (P-E20). The rest of WS-B9 can follow in 0.3.x. With both reports in, M-SIGN should ship flagged experimental (decision **U3**). NEWS.md should list the M-SIGN status, the wording changes and the new output columns.

### (c) Proof register and build order (PD §6, §8)

| Register item / step | Effect of this report | New status or action |
|---|---|---|
| R2 (redundancy) | Two document errors (step 4 intercept; (e) wording); the code is right | stays **P** after WS-B1, with T1c added |
| R3 (sharpness) | Correct for pooled LCE; "LP" needs polyhedral $C$; full-law versus regression example | **P** (pooled), after WS-B1 |
| R10 (inference) | Version table; regime and Condition P qualifications; (e) now analytic | **P** for version (a) and (b2) in fixed dimension; (b1) under the augmented Condition P; (b3) **J**; R3 regime **J** |
| **R15 (new)** | Integer N1: sharpness, refutation and target sets for integer tables | **O** (decision U1); blocks any "exact" label for N1 and Δ4 |
| **R16 (new)** | R3′/R10′ for district rungs (M1, M1(δ), M2, L4), via 2b Thm 1′ | **O**; blocks removing "experimental" from M-SIGN under groups, and the M5/P-N cells for grouped rungs |
| R1 / C1 and P-N (step 2, M-LOC) | P-N's chain fails for pooled M1 (I2); M-LOC rows are linear rows and stay valid | **Not blocked.** Scope P-N and the LOC⇒SIGN-W cell to M0 (pooled) or district-conditional weights; build M-LOC as planned |
| R7 (step 8, N2 grid LP) | Support convention; V-conditioning; derivative non-identification; binning repair; selection | **Blocked** until 2c §4 and §6 are corrected. N2-local rows dropped or redesigned (U4). Global N2 on the empirical support is feasible after WS-B2 |
| R11 / §7 DiD (Δ-rungs, M-SIGN-Δ = Δ4) | "Δ4 exact" inherits P-B1 (continuous only); M-SIGN-Δ under district FD inherits P-B7 | **Not blocked**, but relabel "exact for continuous tables" in PD §4.1, §7 and 2d. Batch C must check Δ4 against P-B1 and P-B7 |
| R9 (step 1b, M-SHARP) | Size-conditioning example (P-B10): conditioning on $X$ alone is not sharp when size is informative | Not blocked. M-SHARP must not be labelled sharp; it is "tighter than $\Theta_O$" |
| Step 3 (2×C joint signs, R6) | M4 confirmed; the integer remark applies to N1 versions only | not affected |
| §4.1 rows M-SIGN, M-ORD, M-N2, M-PER | "sharp (M2)", "proved (M2–M5)", "never (model-free)", "Δ0, Δ4 exact" | rewrite per P-B1, P-B5, P-B7 |
| §4.2 matrix | ORD⇒SIGN-B (per district); LOC⇒SIGN-W (M0 only); SIGN-W⇐N2 (R = 2, finite support, superpopulation reading) | annotate with conditions and proof references |

---

## 4. Workplan

Effort: **S** = under a day of session work; **M** = one to three sessions; **L** = more. Workstreams are listed in the order I recommend.

| WS | Content | Ledger items | Depends on | Effort | Gate (done when) |
|---|---|---|---|---|---|
| **WS-B0** Package labels and wording (now) | `proof_status` experimental for `assoc` and `cr`; the "never refuted" strings (roxygen, note, meta, vignette, PLG §9); LCE-alone status next to a joint rejection; grouped-rung print line; `eiRedundancy` point-level emptiness wording; projection label qualified. **Tests:** EXPERIMENTAL appears in `print()` for sign/CR fits; no unqualified "never be refuted"; no "refuted" in point-level verdicts; face test with intercept (P-B18). | P-B5, P-B7 (labels), P-B16, P-B18 (test), P-B21, P-B30, P-B40 (label) | user approval of this plan | S | full test suite passes; no numeric output changes (compare `test-msign.R` reference values) |
| **WS-B1** 5a errata | R2 step 4; R2(b), (c), (d), (e); §2.5 items 2–4; R3(b) "LP" and examples; R3(e) and l.309; l.311; §3.1 "R1/R3" → regime table; the version table as the R10 status note; R10(e) analytic proof; the §3.1 "double" wording; 5a §5 status rows updated. **Scripts:** T1c, T4b in `verify_R2.py` (plus T4 status recording); S6, S7 in `verify_R3.py` (plus S2 labels); V9–V11 in `verify_R10.py` (plus union merge); regenerate outputs. | P-B15 to P-B27, P-B36 to P-B38, P-B43 | — (can run alongside WS-B0) | M | each new check fails on the old text and passes on the new; outputs regenerated |
| **WS-B2** 2c corrections | M2 relabel and lattice remark; M2(d) frontier; M3 hypotheses, (c), (d); Reading; M5 restated per rung, with the between-district term; Remark M5.1; N2 support convention; N2.2(a) $V$-conditional statement (Conjecture until proved); N2.3 restated with an identified derivative; smooth (d) example; §5 ladder; §6 binning lemma, fixed degree, notation. Plain version (`PLAIN_PHASE2c`) to match. **Scripts:** C7–C13, C5(c′), rebuilt C4 in `verify_phase2c.py`. | P-B1 to P-B5, P-B7 to P-B14, P-B28, P-B29, P-B31, P-B32, P-B34, P-B35 | U1, U4 (wording of the N1 and N2-local parts) | M–L | every **E** row in sections A–C has a corrected statement and a failing-then-passing check |
| **WS-B3** Core 2a/P3 repairs | [supp] with the restricted measure (2a, P1 A5, 2b where used); P3.5(b2)(ii) nonsingularity or reduced parameterisation; residual CLT via Chebyshev/Lindeberg; 1b mapping. **Code:** `.rxc_projection` rank guard and test; zero-variance-candidate check in `.rxc_clr_min`. | P-B23 (code), P-B39 to P-B42, P-B33 | — | M | `verify_phase2a.py`, `verify_phase3.py` gain the two counterexamples; package tests pass; items added to the batch-A response (WS7) |
| **WS-B4** District rungs (R16) | Write R3′ (regression sharpness under M1/M1(δ)/M2 via 2b Thm 1′, with $A_d$ and all free coordinates) and R10′ (regime caveat). **Tests:** M-SIGN under M1 (forms equal direct district computation on a noise-free design; W1 identity with district τ; the referee's design as a false-belief example). | P-B7 | WS-B2 (M5 restated); U2 | M | R16 proved and independently pre-screened; then the grouped-rung print line can return to "sharp for the regression" |
| **WS-B5** Register and design documents | PD §0 (P-N scope), §3 T1.6, §4.1 rows, §4.2 cells, §6 step 8 prerequisites, §8 (R2, R3, R10 statuses; new R15, R16; R7 and R9 notes); NEXT_BUILD (release item); EXTERNAL_REVIEW_PLAN (batch B status, pre-screen lesson); MSIGN_REVIEW_LOG (close M-9 to M-12 with pointers). | §3(c) table, P-B43 | WS-B1, WS-B2 | S | every cell that cites M5(c), P-N, N2.3 or "exact" carries its condition |
| **WS-B6** Decision-dependent work | If U1 = implement: integer N1 (ILP with discrete target sets, separate integer branch feasibility, R15). If U4 = keep N2-local: an identified-derivative design (maintained family plus derivative inference). | P-B1, P-B11 | U1, U4 | L each | — |
| **WS-B7** Vignette and PLG §9 rewrite (EI report) | The substantive rewrite: meaning and weights of each belief (P-E2 to P-E5, P-E9); W2 example at the selector (P-E1); reporting checklist and pre-specification (P-E7); a hand-worked turnout example (P-E24); the false-CR example beside the success example (P-E22); the King/R×C positioning (P-E26); one status key (P-E28); method-named numbers (P-E17). Also 1b fixes P-E29 to P-E32 and PLG §2 (P-E31). | P-E1 to P-E7, P-E9, P-E10, P-E17, P-E22, P-E24, P-E26, P-E28 to P-E32 | WS-B0; WS-B9 (to show the new output) | M | the vignette builds; every user-facing sentence in PLG §9 matches a ledger action |
| **WS-B9** Reporting and output (EI report) | Print all union pieces; base column; conservative or projection beside CLR for CR (U6); pair/W2 CR terms; weight concentration; L5/L6 sign check; verdict labels and hull span; later the drop-one display. **Tests:** three-piece print; base column present; L6 violation reported; verdict header. **Script:** V12 (false CR) in `verify_R10.py`. | P-E6, P-E15, P-E16, P-E17, P-E19, P-E20, P-E22, P-E25 | U6, U7 | M | tests pass; existing numeric outputs unchanged |
| **WS-B8** Responses to both referees | See below | all | WS-B1 to WS-B5, WS-B7 | M | the two response files cover P-B1 to P-B43 and P-E1 to P-E33 |

**Relation to the Monte Carlo session (`claude/ws6-montecarlo`).** No overlap: WS6 reruns the P3 Monte Carlo and the NZ benchmark, not M-SIGN. Three of its results bear on this plan:
* **projection coverage everywhere:** ask that the runs record the rank of $\hat V$ and whether the guard of P-B40 would have fired. A rank drop in any design would raise P-B40's priority;
* **many small ridings (R2):** feeds the regime row added to the R10 status (P-B23);
* **truth on the boundary:** the CR edge is a boundary case of the same kind; the WS6 design informs, but does not replace, the analytic R10(e).

The gated M-SIGN coverage MC (MSIGN_REVIEW_LOG M-17) is not part of WS6. After WS-B1 it matters less, because R10(e) is now analytic. Rerun it with Monte Carlo standard errors only if U3 keeps M-SIGN unflagged.

**How the EI report changed the plan.** The items that were waiting for it are now planned:
* the substantive vignette and PLG §9 rewrite is WS-B7;
* the final wording of the refutation messages uses the referee's three-way distinction (empty plug-in, empty population set, test rejection) and the R×C equal-rates witness (P-E12);
* the 1b EGGHM text is fixed in WS-B7 (P-E29 to P-E32), together with P-B33 and P-B42;
* U2 stands: the referee does not object to grouped rungs, but asks that the output say the sign is one pooled restriction (P-E7);
* reporting changes are a new workstream, WS-B9.

No theory item in part II is new: every theory point the EI referee makes is already in part I.

**How the responses will be written (WS-B8).**
* Two files, one per referee: `B_partial_identification_response.md` (P-B rows) and `B_ecological_inference_response.md` (P-E rows), written after WS-B1 to WS-B5 and WS-B7. Code changes (WS-B0, WS-B9) are described in them but need not be finished first. Where a point appears in both reports, each response answers it in full and notes that the other referee raised it too.
* For the EI referee: acknowledge the false-CR experiment and the union-length correction by name; explain, politely, that the slack dial exists and that our stale 5a §5 misled the referee; report the reporting changes with an example of the new print.
* **Structure:**
  * a short cover note (what changed; what is now labelled continuous-only, M0-only, conjecture or experimental);
  * a point-by-point section in the report's order (major comments 1–7, minor, EGGHM table, scripts, outside scope). Each point quotes the comment, gives the response, and names the new location (file, section, result number) and the check that now fails on the old statement;
  * a revised verdict table using the status labels of `REVIEW_RESPONSE_PLAN.md`. Every row the referee marked *Incorrect* or *Not established* shows either the corrected result or a downgrade to Conjecture or Open.
* **The two severity disputes** are argued, not asserted:
  * P-B1 with the lattice statement and check I1, while accepting the relabelling;
  * P-B3 by noting that γ = 1 forces point identification.
* **The core items (P-B39 to P-B41)** go to referee A as well, in the batch-A response.
* **Credit.** Acknowledge the referee's replacement arguments (R10(e) analytic bound, the binning repair, the residual-CLT argument) by name where they are adopted.
* No model names, and the repository's attribution conventions for commits.

---

## 5. Decisions needed from the user

1. **U1, integer tables (N1, and Δ4 by inheritance).**
   * (a) Recommended: relabel N1 and Δ4 as sharp for **continuous** unit tables. Add the 2×2 lattice remark and an "integer check" item to the backlog. N1 is not built, and in 2×2 with margins only the difference is below one lattice step (about $10^{-5}$ at realistic sizes).
   * (b) Implement integer N1 as an integer program with discrete target sets (new R15; L effort; integer programs are much slower at Wisconsin scale).
   * (c) Both: continuous by default, with an integer feasibility check for small tables.
2. **U2, M-SIGN under grouped rungs (M1, M1(δ), M2, L4).**
   * (a) Recommended: keep it enabled, flagged experimental with a fit-time note saying the sign is one pooled restriction across groups (P-E7), until R16 (R3′/R10′) is proved and tested.
   * (b) Restrict M-SIGN to M0 (error under groups) until then.
   * (c) Leave it as is.
3. **U3, M-SIGN status in release 0.3.0.**
   * (a) Recommended: ship 0.3.0 with M-SIGN flagged experimental after WS-B0, the projection guard, the full union print and the base column (P-E15, P-E20). The EI referee's recommendation (revise claims, examples and default reporting before use in election analysis) supports the flag.
   * (b) Hold 0.3.0 until both batch-B reports are resolved (WS-B1 to WS-B9).
4. **U4, N2's local rows (N2.3) in build step 8.**
   * (a) Recommended: drop N2-local; build global N2 on the finite empirical support with a prespecified degree and person-weighted bins.
   * (b) Keep N2-local with a maintained smooth family and derivative inference (WS-B6, L).
   * (c) Defer N2 entirely.
5. **U5, approval of the immediate work.** Approve WS-B0 (labels and wording only, no numeric change) and WS-B1 to WS-B3 to start now. WS-B4 waits for U2; WS-B7 and WS-B9 wait for U6 and U7.
6. **U6, the interval shown for CR restrictions.**
   * (a) Recommended: keep CLR (branch union) as the default, and always print the conservative slab (W1/2×2) or the projection interval beside it for CR fits, with "algorithmic" on the CLR line itself. This keeps decision D4 and answers P-E17.
   * (b) Make the conservative slab the default for CR fits, with CLR as an option.
   * (c) No change (CLR plus projection, as now).
7. **U7, reporting additions (EI report).**
   * (a) Recommended now (0.3.0): print all union pieces; the base (rung-only) column; weight concentration per signed association. Later (0.3.x): the drop-one display and the L5/L6 sign check.
   * (b) All of WS-B9 before 0.3.0.
   * (c) Only the union-print fix; the rest as vignette advice.
8. **U8, per-group directions (W2 summands).** The EI referee notes that a user who wants each named group to satisfy a direction needs separate restrictions. Each summand $\mu_s\operatorname{Cov}_{w^s}(X_A,m_s)$ is linear in $m$, so a per-group version is one more LP row type.
   * (a) Recommended: document now that W2 is pooled; add a per-group version (`version = "group"`, with `group = s`) in a later build, with its own R3/R10 check.
   * (b) Add it now.
   * (c) Do not add it.
