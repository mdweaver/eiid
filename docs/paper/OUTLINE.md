# Paper outline and page budget

Working title: *Partial Identification for R×C Ecological Inference: Identified Sets, Inference and Software*.

Files: `eiPartialID_paper.tex` (main text), `eiPartialID_supplement.tex` (supplementary appendix with the proofs),
`references.bib`. Notation follows `docs/extension/notation.tex` (V = (X, N, D)). Results are numbered
with a separate counter for each kind (Theorems 1–8, Propositions 1–34, Lemmas 1–12, Corollaries 1–2,
Conjectures 1–3, Remarks 1–5); the supplement proves them under the same numbers. `SOURCE_MAP.md` links each number to its repo label.

Target: main text about 42 pages (11pt, 1in margins, 1.15 spacing), supplement about 25 pages.

| § | Section | Content | Pages |
|---|---|---|---|
| 1 | Introduction | EI and R×C; why point identification needs strong assumptions (Goodman, King, hierarchical models); aim of the bounding approach (what each assumption buys, honest uncertainty); contributions; status conventions (proved / cited / conjecture / algorithm); plan | 4 |
| 2 | Setting and assumptions | data, accounting identity, DD polytope; three targets B, B̃, B*; LCE in homogeneous form, Q and K; three sets and three hulls; assumptions A1–A10 each with a plain reading; what is not assumed | 6 |
| 3 | Identification in R×C under LCE | decomposition; estimand; sharp set (theorem) and linear targets; outer set: reduction, nesting, sharp for the regression; when sharp = outer; closed forms: JKST equivalence, A = simplex, box∩hyperplane, two rows, joint vs collapsed; the LP; units outside A; refutability and undetectable violations; illustration | 7 |
| 4 | Districts and the rung ladder | grouping alone adds nothing; what each rung identifies; sharpness for affine mean models; the ladder theorem; M1 with two rows; when dummies matter; M2, M1(δ), L4/L5/L6; worked examples | 5 |
| 5 | Signs, reinforcement and monotonicity | associations; identities; N1 (sign rows on realised tables), two-row case, joint restrictions; integer counts; K-rows; **M-SIGN module (experimental)**: sharpness R3, redundancy R2, JKST+signs closed form, sandwich, rank and nesting results (C1 conjecture); N2 outer set, Conjecture N2.2′, local identity, binned LP | 7 |
| 6 | Inference | sampling and the source of randomness; estimators; regime table; endpoint structure; consistency; limit theory (fixed system, triangular arrays); confidence intervals: conservative, CLR (cited / algorithm), projection, fixed superset, variance-adaptive allowance; tests (emptiness, M1 vs M2, pooling; wild bootstrap); finite-sample interval; many small districts; sensitivity; inference with sign rows (R10); recipe | 8 |
| 7 | The R package | design; `eiBounds()` inputs/outputs; defaults; rungs, ranges, restrictions; intervals reported; tests; `eiRedundancy()`; the `status:` line of `print()`; what is not implemented | 3.5 |
| 8 | Numerical evidence and an application | verification scripts; Monte Carlo (model holds, many small districts, riding shocks, projection for B under the variance-adaptive allowance); NZ benchmark (run with 0.2.0); M-SIGN simulations (CR at the edge, false CR, redundancy); worked example: voter_vax 2×2 | 5 |
| 9 | Open questions | conjectures and algorithms (CLR default; N2.2′; C1; single-index dimension); open problems (grouped-rung M-SIGN; integer tables; width-gain inference; many small groups; uniform limit theory under changing active sets; Slater failures; hull rate; B*; dependence; turnover); a pointer to difference-in-differences as future work | 2 |
| | References | | 1.5 |

Supplement (proofs, one consistent notation):

| § | Content |
|---|---|
| S1 | Conventions; cited theorems and where they are used |
| S2 | Proofs for Section 3 (identification) |
| S3 | Proofs for Section 4 (districts) |
| S4 | Proofs for Section 5 (signs, M-SIGN, N2) |
| S5 | Proofs for Section 6 (inference, incl. R10) |
| S6 | Status of each result (proved / cited / sketch / conjecture / algorithm) |

Excluded: difference-in-differences (PHASE1b, PHASE2d, Proposition D8, DiD-S): one sentence in §9.
