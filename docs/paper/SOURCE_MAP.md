# Source map (internal; not part of the paper)

Links every numbered result of `eiPartialID_paper.tex` to its label in the repository at
tag `v0.3.1` (commit `c02954e`). The supplement (`eiPartialID_supplement.tex`) proves the results
under the same numbers. Abbreviations: **2a** = `docs/extension/PHASE2a_rxc.md`, **2b** =
`PHASE2b_districts.md`, **2c** = `PHASE2c_monotone.md`, **P3** = `PHASE3_inference.md`, **5a** =
`PHASE5a_sign_module.md`, **PD** = `PACKAGE_DESIGN.md`, **P1** = `PHASE1_foundations.md`. Paper
numbers come from the compiled `.aux`; theorems, propositions, lemmas, corollaries, conjectures and
remarks have separate counters.

## Notation choices

* Number of units: `n` (as in P3), not `p` (as in `notation.tex`), to avoid a clash with the unit
  weights `p_i = N_i / sum N_j` of 2c.
* Number of districts: `D^*` (2b, P3); `D_i` the district of unit i. In Section 5, a bare `D` is a
  gap between two groups' rates (2c, EGGHM).
* Target coefficient vector: `omega` (2a and P3 use `gamma`, 5a uses `h`); `gamma_r` is kept for the
  variance ratio of 2c/5a. Target constant: `varphi_0(beta)` (P3 writes `c(beta)`). Critical
  values: `c-hat`, `c-hat_p` (P3 writes `k-hat`, `k_p`, which clash with the free coordinates `k`).
* Undetectable nonlinear deviation: `nu(x)` (2a Prop. 9 and P3 Prop. S1 write `Delta(x)`, which
  clashes with the W2 association `Delta` and the simplex).
* Symmetric matrix in A5: `Psi` (sources: `D`).
* Sampling conditions P3-S1a and P3-S2 are called (I1) and (I2).

## Main text results

| Paper | Repository label | Source location |
|---|---|---|
| Assumptions A1--A10 | A1--A10 | P1 §2 (A3, A5, A6, A9 as corrected there); readings from PLAIN_LANGUAGE_GUIDE §2 and report §2.5 |
| Lemma 1 (unit polytope) | Lemma 0 | 2a §1 |
| Lemma 2 (identification decomposition) | Lemma 1 (a)--(c) | 2a §2 |
| Remark 1 (A5 in practice) | Remarks 1.1, 1.2 | 2a §2 |
| Lemma 3 (estimand) | Lemma 2 (and Remark 2.1, in the supplement) | 2a §3 |
| Theorem 1 (sharp set) | Theorem 1, Remark 4.1 | 2a §4 |
| Corollary 1 | Corollary 1 | 2a §4 |
| "Continuous versus integer rates" paragraph | Remark 4.2 (1) proved, (2) sketch | 2a §4 |
| Lemma 4 (vertex reduction, boundedness, validity) | Lemma 3 (a)--(c); Remarks 5.1, 5.2 | 2a §5 |
| Proposition 1 (nesting) | Proposition 1 | 2a §5 |
| Proposition 2 (sharp for the regression) | Proposition 2 | 2a §5 |
| Proposition 3 (when the outer set is sharp) | Corollary 2 (a), Proposition 3 (b), Remark 3.1 | 2a §5 |
| Proposition 4 (JKST equivalence) | Proposition 4 | 2a §6.1 |
| Proposition 5 (whole simplex) | Proposition 5 | 2a §6.2 |
| Lemma 5 (box ∩ hyperplane) | Lemma 4 | 2a §6.3 |
| Proposition 6 (two rows) | Proposition 6 | 2a §6.4 |
| Corollary 2 (joint vs collapsed) | Corollary 6.1 | 2a §6.4 |
| Proposition 7 (general LP) + value-function paragraph | Proposition 7 and "Regularity of the value function" | 2a §6.5 |
| Proposition 8 (units outside A) | Proposition 7′ | 2a §6.6 |
| Proposition 9 (non-emptiness) | Proposition 8 | 2a §7 |
| Proposition 10 (undetectable violations) | Proposition 9 (with Remark 1.0) | 2a §7 |
| §3.6 illustration; Cross–Manski paragraph | BC 3×3 case; relation to Cross & Manski | 2a §8, §12 |
| Proposition 11 (grouping adds nothing) | Proposition D0 and the joint DD set | 2b §1 |
| Lemma 6 (district identification) | Lemma D1 (a)--(d), Remark D1.1 | 2b §2 |
| Lemma 7 (district estimand) | Lemma D2 | 2b §2 |
| Theorem 2 (affine mean models) | Theorem 1′ | 2b §3 |
| Theorem 3 (rungs and ladder) | Theorem 2 (a)--(c) | 2b §3 |
| Proposition 12 (M1, two rows) | Proposition D3 (incl. pure-composition vertices) | 2b §4 |
| Proposition 13 (when baselines matter) | Proposition D4 (a), (b), Remark D4.1 | 2b §4 |
| Proposition 14 (M2) | Proposition D5 | 2b §4 |
| Proposition 15 (M1(δ)) | Proposition D6 | 2b §4 |
| Proposition 16 (L5, L4) | Proposition D7, Remark D7.1 | 2b §4 |
| §4.5 worked examples | Worked examples | 2b §6 |
| Definition 1 (associations) + reading | Definition 1 and its remarks | 2c §1 |
| Lemma 8 (accounting identities) | Lemma M1 (a)--(d) | 2c §1 |
| Definition 2 (N1 set) | Definition 2 | 2c §2 |
| Proposition 17 (N1 sharp, continuous tables) | Proposition M2 (a)--(d) | 2c §2 |
| Proposition 18 (integer counts, 2×2) | Remark M2.1 (proposition inside it) | 2c §2 |
| Proposition 19 (two rows) | Proposition M3 (a)--(d) and "Reading" | 2c §2 |
| Proposition 20 (joint restrictions) | Proposition M4 | 2c §2 |
| Proposition 21 (sign rows under linear rungs) | Proposition M5 (a)--(c) | 2c §3 |
| Remark 2 (realised vs model-implied) | Remark M5.1 | 2c §3 |
| Theorem 4 (sharpness with added restrictions) | Theorem R3 (a)--(e) | 5a §1.2--1.3 |
| Remark 3 (full law vs regression) | Remark R3.1 | 5a §1.2 |
| Remark 4 (LCE can exclude CR) | Remark R3.2 | 5a §1.4 |
| Paragraph after Remark 4 (other CR pairings) | §1.4 "What fails" | 5a §1.4 |
| Proposition 22 (tightening) | Proposition R2 (a)--(e) | 5a §2 |
| Proposition 23 (JKST plus signs) | G2 (proof §2.2); audit rows G1, G3 in the reading | PD §2.1--2.2 |
| Proposition 24 (sandwich) | Proposition P-S | PD §5.1 |
| Proposition 25 (pairwise directions) | Proposition P-R (a), (b), (d); (c) stated as checked only | PD §5.2 |
| Proposition 26 (nesting in 2×2) | Proposition P-N (a), (b) | PD §5.3 |
| Conjecture 1 | Conjecture C1 (register R1) | PD §5.3 |
| Lemma 9 | Lemma N2.1 | 2c §4 |
| Proposition 27 (N2 outer set) | Proposition N2.2 (a)--(e), Remarks N2.2.2 | 2c §4 |
| Conjecture 2 | Conjecture N2.2′ (Remark N2.2.1 and the sufficiency remark) | 2c §4 |
| Proposition 28 (local identity) | Proposition N2.3 (a)--(d) and "Not an observable restriction" | 2c §4 |
| Lemma 10 (person-weighted bins) | Lemma N2.4 | 2c §4 |
| "What inference for N2 needs" | §6 (degree selection) | 2c §6 |
| (I1), (I2); "Where the randomness comes from" | P3-S1a, P3-S2; §0.4 | P3 §0.1, §0.4 |
| Table 2 (regimes) | §0.5 table | P3 §0.5 |
| Lemma 11 | Lemma P3.1 | P3 §2 |
| Proposition 29 | Proposition P3.2 (and remarks) | P3 §3 |
| Theorem 5 | Theorem P3.4 (i)--(iii) (fixed constraint system) | P3 §4 |
| Remark 5 | Remark P3.4.1 | P3 §4 |
| Proposition 30 | Proposition P3.4′ (register R17) | P3 §4 |
| Theorem 6 (a) | Theorem P3.5(a) | P3 §5 |
| Theorem 6 (b1) | P3.5(b1) (cited, CLR Thm 4) | P3 §5 |
| Theorem 6 (b2) | P3.5(b2)(i), (ii) with the exact-reduction remark | P3 §5 |
| Conjecture 3 | Conjecture P3.5(b3) (algorithm) | P3 §5 |
| Lemma 12 | proved reduction in Conjecture P3.5(b3) | P3 §5 |
| Proposition 31 | Proposition P3.5(c): (c1) proved; (c2) in the reading (register R18) | P3 §5 |
| Boundary-cases paragraph | "Boundary cases" | P3 §4 |
| Proposition 32 (tests) | Proposition P3.6 (a)--(c) | P3 §5 |
| Wild-bootstrap paragraph | "Calibration in the package" | P3 §5 |
| Theorem 7 (finite-sample interval) | Theorem P3.8 (a)--(d) and remarks | P3 §5.1 |
| Proposition 33 | Proposition P3.7 (a), (b) | P3 §6 |
| Proposition 34 | Proposition S1 and breakdown value | P3 §7 |
| Theorem 8 | Theorem R10 (a)--(e); Table 3 = 5a §3.2 status table | 5a §3 |
| §6.9 recipe | Practical recipe | P3 §9 |
| §7 package | NEWS.md 0.3.0/0.3.1, `R/rxc_api.R` (`eiBounds`, `.rxc_ci_status`, `print.eiRxC`), `R/rxc_sign.R` (`eiRedundancy`), PHASE4, vignettes | package/eiPartialID |
| §8.2--8.3 Monte Carlo (Tables 4--6) | §8.1--8.4 (projection-for-B column from the rerun with P3.5(c1)) | P3 §8; `phase3_montecarlo_output.txt` |
| §8.3 test calibration table | calibration table | P3 §5; `benchmarks/results/*_size_output.txt` |
| §8.4 NZ benchmark | benchmark summary | `benchmarks/results/ei_datasets_nz_summary.md` (run with 0.2.0) |
| §8.5 sign-module simulations (Table 7, false CR) | §3.4 (V11--V13), PD §2.2 [V2], 5a §2.4 T4, 5a §5 gated run | 5a; `verify_R10_output.txt` |
| §8.6 worked example | computed with eiPartialID 0.3.1 on `tests/testthat/voter_vax_MonotoneEI.csv` (seed 1); matches `test-msign.R` ("about [-0.470, -0.050]") | this paper |

## Coverage of the technical report's status table (§12), DiD items excluded

Every row of `docs/report/sections/12_status.tex` is either in the paper or omitted deliberately.

| Status-table row | In the paper | Note |
|---|---|---|
| P1 Lemma 1; Prop. 2-RC; closed form R=2 and joint-vs-collapsed corollary; Theorem 1 | Lemma 2; Lemma 4, Props 5, 7; Prop. 6, Cor. 2; Thm 1 | stated in their 2a form (where they are proved) |
| P1 nesting (1.13); equality conjecture | Prop. 1; Prop. 3 | only the proved sufficient condition is stated |
| P1 (D0)--(D5) | Prop. 11, Lemmas 6--7, Props 13--14 | stated in their 2b form |
| P1 assumptions A1--A10 | §2.5 | |
| 2a Lemma 0 -- Prop. 9, Remarks 1.0--1.2, 2.1, 4.1, 4.2, 5.1, 5.2, 3.1, Cor. 6.1, value-function regularity | §3 (see table above) | Remark 4.2(2) stated as a sketch, not used |
| 2a Cor. 6.1 R≥3 characterisation | §3 reading of Cor. 2; §9 | open |
| 2b D0, joint DD set, D1, D1.1, D2, Thm 1′, Thm 2, D3, D4, D4.1, D5, D6, D7, D7.1 | §4 | |
| 2b Prop. D8 | **omitted** | difference-in-differences material, excluded by scope; one future-work sentence in §9.3 |
| 2b Remarks 6.1--6.3 | 6.1: readings of A9/A10 and Prop. 33; 6.3: §4 opening sentence | 6.2 (periods as districts) only as the future-work sentence |
| 2b strict-gain characterisation R≥3; M1′ | §9 (strict gains) | M1′ (covariate-dependent effects) omitted: covariates are not developed in the sources |
| 2c Lemma M1, Props M2, M3, M4, M5, Remarks M2.1, M5.1, Lemma N2.1, Prop. N2.2, Conj. N2.2′, Prop. N2.3, Lemma N2.4, degree selection | §5 | |
| P3 Lemma P3.1, P3.2, P3.4, P3.4′, P3.5(a), (b1), (b2), (b3), (c1), (c2), P3.6, wild bootstrap, P3.8, P3.7, S1, breakdown | §6 | (b3) as Conjecture 3, (c2) in the reading of Prop. 31 |
| P3 M1/M2 in R3; CLR in R2 | Table 2; §9 | conjecture; simulation finding |
| 5a Theorem R3, Remarks R3.1, R3.2, Prop. R2, Theorem R10 (a)--(e) | §5.4, §6.8 | pooled M0 |
| PD G2, P-S, P-R, P-N, Conjecture C1 | Props 23--26, Conj. 1 | |
| PD DiD-S | **omitted** | difference-in-differences, excluded by scope |
| PHASE1b, PHASE2d items | **omitted** | difference-in-differences proposals, excluded by scope |

## Proof register (PD §8)

| ID | In the paper |
|---|---|
| R1 | Conjecture 1 |
| R2 | Proposition 22 |
| R3 | Theorem 4 |
| R4 | §9 (refutability for R≥3) |
| R5 | Proposition 25 (star proved; general trees in §9) |
| R6 | §9 (joint sign restrictions across columns) |
| R7 | §9 (outer N2 LP with an estimated regression) |
| R8 | **omitted**: modules B1(M), κ, CAR-δ belong to NEW_AVENUES and are outside the scope list |
| R9 | §9 (series bound towards Θ_X) |
| R10 | Theorem 8 |
| R11 | **omitted**: difference-in-differences |
| R12 | Proposition 24 |
| R13 | Proposition 23 |
| R14 | §9 (single-index dimension) |
| R15 | Proposition 18; §9 |
| R16 | §5.4 opening; §9 |
| R17 | Proposition 30 |
| R18 | Proposition 31 |

## Other material deliberately left out

* The audit rows G4--G9 of PD §2.1 (verdicts on an outside proposal): their content appears where
  it is proved (Props 21, 23--26, Thm 8); the audit itself is not reproduced.
* NEW_AVENUES modules, M-LOC/M-ORD build details, the build order of NEXT_BUILD: not in scope.
* Verification-check labels ([F], [V], [I], C1--C16, S1--S7, T1--T4, V1--V13): the paper reports the
  checks' content (in readings and §8.1) without the labels.
