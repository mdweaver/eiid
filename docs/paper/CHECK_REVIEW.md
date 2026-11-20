# Check of the paper against the repository

**Paper checked:** branch `claude/paper` at `bc9bc6c`, folder `docs/paper/`.
**Sources of truth:** branch `claude/wizardly-volta-lbh4bm` at `44db8f1`. That commit fixes items 1–4 of the writer's `ISSUES.md`. The `docs/report/` folder on the same branch was also consulted (its §12 status table and register).
**What was done:**
- read every section of the main text and the supplement, and compared each numbered result with its source statement;
- rebuilt both PDFs with `build.sh` (TeX Live, pdflatex + bibtex);
- installed `package/eiPartialID` 0.3.1 from the source branch and reran the `voter_vax` worked example with seed 1;
- recomputed the Monte Carlo, test-size, NZ and sign-module numbers from their source tables;
- searched the `.tex` sources and the compiled PDF text for process language.

**Scope of this check:** it compares the paper with the sources. It does not re-derive the sources' proofs and does not judge whether the paper is publishable.

---

## 1. Verdict on each point of the owner's specification

| # | Specification point | Verdict | Notes |
|---|---|---|---|
| 1 | Journal-style methods paper, about 35–50 pp. of main text, with a supplement | **Met** | The main-text PDF has 52 pages. Of these, the text runs from p. 4 (after the title page and contents) to p. 51, where the references start: about 47–48 pages, near the top of the range. The supplement has 25 pages. |
| 2 | Process-free, with every correction present as the current result | **Met** | There are no references to referees, reviews, rounds, responses, ledgers, workstreams, sessions or revision history (§D). Every corrected statement checked appears in its corrected form (§A). One sentence gives software-version history, which is not review process (§D, optional fix). |
| 3 | Opens with the aims, then the methods in a logical order | **Met** | §1 runs problem → aim of the bounding approach → contributions → related work → plan. The methods then follow in order: setting → R×C → districts → signs → inference → package → evidence → open questions. |
| 4 | Scope | **Met** | It covers the core bounds, rungs, districts and inference; the monotone results and M-SIGN, labelled experimental; the package; and the evidence. Difference-in-differences appears only as the future-work paragraph in §9.3. Prop. D8, DiD-S, PHASE1b and PHASE2d are omitted, as the specification requires. |
| 5 | Formal statements, plain readings and sketches in the main text; detailed proofs in one notation in the supplement | **Met** | Every key result has a "Reading", and the main results have a "Sketch". The supplement proves every numbered result in the paper's notation (`n` for units, `ω` for target coefficients, `ν` for the undetectable deviation, `ĉ` for critical values). It is rewritten, not a copy of the design documents. Gaps are listed under §C. |
| 6 | Conjectures and algorithms labelled where they appear and collected at the end; nothing unproved presented as proved | **Met** | The labelled items are: Conjectures 1–3; the CLR procedure ("algorithm"); the R3 regime ("conjecture" in Table 2); the CLR result in R2 ("simulation finding"); the wild bootstrap ("supported by simulation"); P-R(c), the spanning-tree claim ("checked … not proved"); and Remark 4.2(2) ("only a sketch"). All are collected in §9.1–9.2. The status tags agree with the report's §12 table, and no unproved result carries a "proved" tag. |
| 7 | Compiles with pdflatex + bibtex (`build.sh`) | **Met** | The build ends with "build ok". It produces 52 + 25 pages, with no undefined references or citations, no overfull boxes, and only hyperref "Token not allowed in a PDF string" warnings (§H). |

**Overall:** the paper is faithful to the corrected sources. **Nothing must be fixed for correctness.** Five fixes are recommended (§I): two wording corrections and three source-tracking updates. The rest are optional.

---

## A. Faithfulness of every numbered result

The counts are 8 theorems, 34 propositions, 12 lemmas, 2 corollaries and 3 conjectures: **59 numbered results**, plus Assumptions A1–A10, Definitions 1–2 and Remarks 1–5. Sources are abbreviated as in `SOURCE_MAP.md`; all were read at `44db8f1`.

| Paper | Source | Verdict | Note or fix |
|---|---|---|---|
| A1–A10 | P1 §2; report §2.5 | OK | A5 uses the support of the restricted measure, `supp_A P_X`, with the uniform-on-[0,1], `A={1/2}` example. A6 makes independence a substantive assumption (P3 §0.4). |
| Lemma 1 | 2a Lemma 0 | OK | |
| Lemma 2 | 2a Lemma 1 | OK | |
| Remark 1 | 2a Remarks 1.1, 1.2 | OK | |
| Lemma 3 | 2a Lemma 2 | OK | Remark 2.1 appears in the supplement after the proof. |
| Thm 1 | 2a Thm 1 | OK | It states one `P_V`-null set for all `u`, and compactness when `affdim(A)=R-1`, which A5 implies. |
| Cor 1 | 2a Cor 1 | OK | |
| "Continuous vs integer" paragraph | 2a Remark 4.2 | OK | (2) is labelled as a sketch that is not used. |
| Lemma 4 | 2a Lemma 3, Remarks 5.1, 5.2 | OK | (c) has the hypothesis `A ⊆ conv(supp_A P_X)`. The reading covers the singleton-`A` cases. |
| Prop 1 | 2a Prop 1 | OK | |
| Prop 2 | 2a Prop 2 | OK | |
| Prop 3 | 2a Cor 2, Prop 3, Remark 3.1 | OK | It says the condition is "sufficient, not necessary", and gives 46/400. |
| Prop 4 | 2a Prop 4 | OK | |
| Prop 5 | 2a Prop 5 | **MISMATCH (minor)** | The displayed set omits the diagonal condition. Fix: write `Θ_O = {K: |K_{c,rs}| ≤ Q_{c,rs} ∀c, r<s}`, "which also requires `Q_{c,rr} ≥ 0`". The non-emptiness sentence and the supplement proof already include it. |
| Lemma 5 | 2a Lemma 4 | OK | |
| Prop 6 | 2a Prop 6 | OK | |
| Cor 2 | 2a Cor 6.1 | OK | It includes the open characterisation for R≥3 and the 20 strict 3×3 cases. |
| Prop 7, value-function paragraph | 2a Prop 7 | OK | `γ` is renamed `ω`. |
| Prop 8 | 2a Prop 7′ | OK | It has `P(X∈A)>0`, (S2)–(S3) on `A`, A5 on `supp_A P_X`, `μ^A`, and the outside terms; each end is attained. |
| Prop 9 | 2a Prop 8 | OK | |
| Prop 10 | 2a Prop 9 | **MISMATCH (nit)** | The last sentence says "Under A5". The source says "Under (A5) on `supp_A P_X` (with `A ⊇ supp P_X`)". Fix: add that qualifier. |
| Prop 11 | 2b Prop D0, joint DD set | OK | |
| Lemma 6 | 2b Lemma D1 (a)–(d), Remark D1.1 | OK | (d) has the stacked rank condition, sufficiency, and the requirement of affinely independent rather than merely distinct compositions. This is item 2b in the list to check. |
| Lemma 7 | 2b Lemma D2 | OK | |
| Thm 2 | 2b Thm 1′ | OK | |
| Thm 3 | 2b Thm 2 | OK | The reading frames M1(0) as a sensitivity benchmark. |
| Prop 12 | 2b Prop D3, pure vertices | OK | |
| Prop 13 | 2b Prop D4, Remark D4.1 | OK, **incomplete (nit)** | (b) is correct and has the `(X,N)` conditioning; the counterexample to D⊥X alone is in the reading. It omits two parts of the source: the second sufficient condition ("`D⊥N | X` with `E[Z_D | X]` affine"), and the sentence that M0's (S3) is automatic only on the population hull. Optional: add both. |
| Prop 14 | 2b Prop D5 | OK | |
| Prop 15 | 2b Prop D6 | OK | |
| Prop 16 | 2b Prop D7, Remark D7.1 | OK | |
| Def 1 and its reading | 2c Def 1 and remarks | OK | It includes the W1/W2 example (−0.02633 / +0.00110). |
| Lemma 8 | 2c Lemma M1 | OK | The hypotheses are now part by part, as corrected. |
| Def 2 | 2c Def 2 | OK | |
| Prop 17 | 2c Prop M2 | OK | It is sharp for continuous tables only. (d) has componentwise nesting, the breakdown frontier and the scalar path. This is item 2c. |
| Prop 18 | 2c Remark M2.1 | OK | It covers the integer lattice: the step `1/N₁+1/N₀`, 4.8×10⁻⁶ for `N=10⁶, μ=.3`, the three-unit example, and the open cases. |
| Prop 19 | 2c Prop M3 | OK | It assumes no unit restrictions; (d) has the cases γ<1, γ=1 and `D_ER=0`. The reading says the equal-rates witness satisfies the CR **union** and only the branch that matches its within association. This is item 2c (M3). |
| Prop 20 | 2c Prop M4 | OK | The reading has the "same weak within sign" version. |
| Prop 21 | 2c Prop M5 | OK | It has the stacked `m_id + L_i k` form and "(c) for M0 only". The district decomposition and the six-unit reversal example (+0.00064 / −0.04927) are in the reading and the proof. This is item 2c (M5). |
| Remark 2 | 2c Remark M5.1 | OK | The standard-deviation bound is under A6, and "a mean condition alone gives no rate". |
| Thm 4 | 5a Thm R3 | OK | It is for pooled M0 only, with `affdim(A)=R-1`, the (C1)–(C2) admissibility conditions, "finite LP when C is polyhedral", and the certificate as sufficient only. |
| Remark 3 | 5a Remark R3.1 | OK | |
| Remark 4 | 5a Remark R3.2 | OK | |
| "Other CR pairings" paragraph | 5a §1.4 | OK | It has the gaps of up to about 31% and the pairwise pooling caveat. |
| Prop 22 | 5a Prop R2 | OK | The face LP carries the intercept: R2 face test with intercept, as listed in item 5a. It includes the dual-multiplier caveat and (e) on zero incremental gain. |
| Prop 23 | PD G2; 5a R2(d) | OK | |
| Prop 24 | PD P-S | OK | |
| Prop 25 | PD P-R (a), (b), (d) | OK | (c) is stated in the reading as checked numerically and not proved. |
| Prop 26 | PD P-N | OK | |
| Conj 1 | PD C1 (R1) | OK | |
| Lemma 9 | 2c Lemma N2.1 | OK | |
| Prop 27 | 2c Prop N2.2 | OK | It uses the finite-support convention, and (b) fails on a grid. |
| Conj 2 | 2c Conj N2.2′ and the sufficiency remark | OK | |
| Prop 28 | 2c Prop N2.3 | OK | (d) says "on a neighbourhood of x". Finite-support monotonicity is not enough, and the local rows are not observable. This is item 2c (N2.3(d)). |
| Lemma 10 | 2c Lemma N2.4 | OK | It uses person weights and defines `τ̄_b` as the model-implied bin mean, with the estimation caveat in the reading. This is item 2c (N2.4). |
| (I1), (I2), "Where the randomness comes from" | P3 §0.1, §0.4 | OK | Independence is a substantive assumption, the pooling test is not a test of dependence, and nothing is covered under dependence. This is item P3 (§0.4). |
| Table 2 | P3 §0.5 | OK | |
| Lemma 11 | P3 Lemma P3.1 | OK | It is tagged "fixed constraint system". |
| Prop 29 | P3 Prop P3.2 | OK | |
| Thm 5 | P3 Thm P3.4 | OK | It is "cited structure" for a fixed constraint system, and (ii) has the variance condition. This is item P3 (P3.4). |
| Remark 5 | P3 Remark P3.4.1 | OK | |
| Prop 30 | P3 Prop P3.4′ | OK | It has the localisation condition and the `min{Z,1−Z}` example. This is item P3.4′. |
| Thm 6 (a) | P3.5(a) | OK | It has the weight condition. The proof uses the Chebyshev/Lindeberg subsequence argument. |
| Thm 6 (b1) | P3.5(b1) | OK | It is cited as CLR Thm 4 with Condition P. |
| Thm 6 (b2) | P3.5(b2)(i)/(ii) | OK | It has `q = dim β` or an exact reduction. The reading says a numerical rank cutoff is not an exact reduction and gives the null-direction counterexample. This is item P3.5(b2). |
| Conj 3 | P3.5(b3) | OK | It is labelled an algorithm. |
| Lemma 12 | P3.5(b3), the proved reduction | OK | |
| Prop 31 | P3.5(c) | OK | (c1) requires the input interval from the rung alone and `max w_i²/s̄² → 0`, the condition on `s̄`. (c2) is in the reading and inherits Conjecture 3. This is item P3.5(c1)/(c2). |
| Boundary-cases paragraph | P3 §4 | OK | |
| Prop 32 | P3 Prop P3.6 | OK | (c) is conditional on the design, needs R affinely independent points, and has `(D*−1)(C−1)R` restrictions. The reading says it has no power against dependence. This is item P3.6(c). |
| Wild-bootstrap paragraph | P3 "Calibration in the package" | OK | |
| Thm 7 | P3 Thm P3.8 | OK | It covers all four parts, (a)–(d). This is item P3.8. |
| Prop 33 | P3 Prop P3.7 | OK | (b) transfers only under the regime conditions. |
| Prop 34 | P3 Prop S1 and breakdown | OK | It is labelled an outer bound, with η* = 0.060. |
| Thm 8 and Table 3 | 5a Thm R10 and the §3.2 table | OK | It is for pooled M0 only, with status by interval version. (a) has the contrast envelope, (b) is conservative only, (c) is point coverage, (d) is pointwise, and (e) gives 0.933525. This is item 5a (R10 pooled). |
| §8 "plug-in bracket" text | 5a §3.5 item 4 | OK | It gives 0.917 = 550/600 for `CR-naive` and 304/600 (51%) for the literal intersection, with "at most one half". This is item 5a §3.5. The paper gives 0.917 rather than the count 550/600; optionally add the count. |

**Summary:** 59 numbered results checked. There are **0 substantive mismatches** and **3 minor ones**: Prop 5 (diagonal condition), Prop 10 (scope of A5) and Prop 13 (an omitted sufficient condition). Every corrected statement in the list to check matches its current source.

**The four source fixes in `44db8f1`:**
1. **P3 §8.1 projection-B values:** the paper's Table 4 has 0.251 (design E) and 0.378 (demobilisation). Both match.
2. **NZ results stated as tests:** see §E. The paper states the test results, but it puts the empty plug-in under "the data reject the model". **This should be reworded.**
3. **The [Lyap] citation:** the supplement cites "Lindeberg CLT … Billingsley (1995), Section 27". This agrees with the source. Optionally add "Theorem 27.2", as the source now does.
4. **The M-SIGN status string:** §7 paraphrases it ("a line marking the module as experimental, with its proof status"), which is consistent with the new string ("experimental (sharpness and inference proved for the pooled model M0 only; …)").

---

## B. Coverage (own list, from report §12 and PD §8)

### Report §12 status table

| Item | Location in the paper | Verdict |
|---|---|---|
| P1 Lemma 1; Prop 2-RC; R=2 closed form and corollary; Thm 1; nesting and the equality conjecture; (D0)–(D5); A1–A10 | Lemma 2; Lemma 4 and Props 5, 7; Prop 6 and Cor 2; Thm 1; Props 1 and 3; Prop 11, Lemmas 6–7, Props 13–14; §2.5 | Covered, in their proved 2a/2b forms |
| 2a Lemma 0 – Prop 9, Remarks 1.0–5.2, Cor 6.1, value-function regularity | §3 | Covered |
| 2a, characterisation for R≥3 | Reading of Cor 2; §9.2 | Covered (open) |
| 2b D0 – D7, D1.1, D4.1, D7.1 | §4 | Covered |
| 2b Prop D8 | — | **Omitted, with sound reason** (DiD scope) |
| 2b Remarks 6.1–6.3 | 6.1 in the readings of A9/A10 and Prop 33; 6.3 in the §4 opening; 6.2 only as future work | Covered as appropriate |
| 2b strict-gain characterisation R≥3; M1′ | §9.2 (strict gains) | Covered. M1′ is omitted: covariates are not developed in the sources. Sound. |
| 2c M1–M5, M2.1, M5.1, N2.1–N2.4, Conj N2.2′, degree selection | §5; "What inference for N2 needs" | Covered |
| P3 P3.1, P3.2, P3.4, P3.4′, P3.5 (a), (b1), (b2), (b3), (c1), (c2), P3.6, wild bootstrap, P3.8, P3.7, S1, breakdown | §6 | Covered |
| P3, M1/M2 in R3 and CLR in R2 | Table 2; §9.1 | Covered |
| 5a R3, R3.1, R3.2, R2, R10 (a)–(e) | §5.4, §6.8 | Covered |
| PD G2, P-S, P-R, P-N, C1 | Props 23–26, Conj 1 | Covered |
| PD DiD-S; PHASE1b; PHASE2d | — | **Omitted, with sound reason** (DiD scope) |

### PD §8 register (R1–R18)

| ID | Location in the paper | Verdict |
|---|---|---|
| R1 | Conj 1 | Covered |
| R2 | Prop 22, and Prop 23 for the 2×2 closed forms | Covered |
| R3 | Thm 4 | Covered |
| R4 | §9.2 ("refutability of sign, ordering and local beliefs when R≥3") | Covered (open) |
| R5 | Prop 25; the chain is in the reading; §9.1 | Covered. The source's sub-item "behaviour with lower-dimensional support" is not mentioned; optional to add. |
| R6 | §9.2 | Covered |
| R7 | Lemma 10; Conj 2; §9.2 ("an outer N2 LP with an estimated regression") | Covered |
| R8 | — | **Omitted, with sound reason**: B1(M), κ and CAR-δ are NEW_AVENUES modules, outside the scope list |
| R9 | §9.2 (series bound) | Covered |
| R10 | Thm 8, Table 3 | Covered |
| R11 | — (future-work pointer only) | **Omitted, with sound reason** (DiD) |
| R12 | Prop 24 | Covered |
| R13 | Prop 23 | Covered |
| R14 | §9.1 ("a related conjecture that single-index contextual effects have dimension C−1") | Covered |
| R15 | Prop 18; §9.2 | Covered |
| R16 | §5.4 opening; §7; §9.2 | Covered |
| R17 | Prop 30 | Covered |
| R18 | Prop 31 | Covered |

**Missing proofs and inference results:** no main-text result is left without a proof in the supplement. Every method has its inference result:
- linear rungs: Thms 6–7;
- districts: Table 2 with the regime caveat;
- sign rows: Thm 8, for pooled M0. For district rungs with sign rows the paper says that no coverage theorem exists, which is correct.
- N1: needs no sampling inference.
- N2: inference is stated as conditional on a prespecified degree, which is appropriate because N2 is not implemented.

---

## C. Proofs in the supplement

A script cross-referenced every main-text label against the supplement. **Every numbered lemma, proposition, theorem and corollary has a proof environment.** Theorem 1 is proved as "Proof of necessity" and "Proof of sufficiency", Theorem 6 in two parts. The proofs read in full were compared with the source proofs step by step: S2 (all of §3), S3 (all of §4), and in S4/S5 the proofs of Prop 21, Prop 28, Lemma 10, Lemma 11, Prop 29, Thm 5, Prop 30, Thm 6(a)/(b2), Prop 31, Prop 32 and Thm 8. The other S4 proofs were checked for presence and length only.

**Complete and faithful.** These proofs are in one notation and drop no step:
- Thm 1: necessity has the countable dense set; sufficiency has steps F1–F7 and the contradiction at a random direction.
- Lemma 4(b), together with the proof that A5 implies `affdim(A)=R−1`.
- Prop 5, Cor 2, Prop 8.
- Lemma 6(d), with the stacked-rank sufficiency argument.
- Prop 13(b), with the counterexample.
- Thm 6(a), with the Chebyshev/Lindeberg subsequence argument.
- Thm 6(b2), with weak duality and the Cauchy–Schwarz step for (ii).
- Prop 31: the (c1) argument with `s̄`, and the (c2) remark.
- Thm 8(e): the 7/12 ratio, the union bound over eight rows and two tails, and the simulated-critical-value bound 0.940092.
- The proof that the literal intersection covers at most one half.

**Sketch-only sources are flagged:**
- Remark 4.2(2) (randomising over integer vertices): the supplement says "only a sketch … not used".
- Conjecture 2: the supplement gives the sufficiency sketch and says that necessity "is not written".
- Thm 2 is proved by step-by-step reduction to Thm 1, as in the source, and says so.
- Thm 5 is labelled "cited structure".

**Small gaps (optional):**
1. **Prop 28(d)** (`supplement/d_signs.tex`, proof of N2.3): "non-negative derivatives give `δ_W(x) ≥ 0`" skips the step "differentiable at x and non-decreasing on a neighbourhood of x ⇒ `m'_{rc}(x) ≥ 0`". It is one line, and the hypothesis on the neighbourhood is what the corrected statement needs, so state it.
2. **Thm 8 proof** (`supplement/e_inference.tex`, from l. 252): the step labels jump from (1)–(7) to (8)–(23), interleaved with references to steps of Thm 6. That is the source's continuous numbering. It reads oddly but is not wrong. Optional: renumber within each part.
3. **Prop 13 (D4)**: as in §A, the proof does not cover the omitted second sufficient condition, because the statement omits it.

**No proof is missing or truncated** relative to the source.

---

## D. Process-free check

**Searches run:**
- in the `.tex` files of `sections/`, `supplement/` and the two root files: referee, review, round, respons*, decision, workstream, ledger, session, revision, revised, earlier version, earlier, previously, originally, "we correct", corrected, "now ", "after ", batch, WS[0-9], P-B, P-E, U[0-9](, R1–R18 IDs, footnote, superseded, stale, "had been", "used to";
- in the compiled PDF text (pdftotext): referee, review, round, response, ledger, session, revision, earlier version, previously, originally, corrected, superseded, "later change", "before the".

**Hits:** every hit is benign:
- "earlier election" (s1 l.9, s2 l.10);
- "after relabelling", "after clipping", "after dividing";
- `thm:R10` and `tab:R10`, which are labels, not register IDs;
- "Review" in journal titles in the references;
- "now includes" in the supplement (s5 l.273), a step in a proof.

There are no footnotes.

**One borderline hit**, about software rather than review: §8.4 (s8, "A benchmark with known truth") says *"The benchmark was run with version 0.2.0 of the package, before the later changes to the tests and to the variance bound of the projection interval for B …"*. This is a provenance statement that a reader needs, not review process, but "before the later changes" reads like revision history. Optional rewording: *"The benchmark used version 0.2.0 of the package. Version 0.3.1 calibrates the tests by a wild bootstrap and takes the variance bound of the projection interval for B from Proposition 31; these differences can move the projection and test columns slightly."*

**Superseded statements:** none found. Every place where a source carries a "Corrected after review" note shows only the corrected text in the paper. The ones checked are:
- 2a Remark 5.1/5.2, Prop 7 value function;
- 2b D1, D1(d), D4(b), D6(d), Remark 6.1;
- 2c M2, M3, M4, M5, M5.1, N2.2, N2.3, N2.4, §6;
- P3 §0.1, §0.4, P3.2 remarks, P3.4, P3.5(a) step 5, P3.5(b2)(ii), P3.6(c), P3.7;
- 5a R3 setting, R3(b), R3(e), §1.4, R2(a)/(b)/(e), R10 scope, R10(e), §3.4 union length, §3.5 item 4.

**Counts:** 0 hits for review process; 1 borderline hit for software-version history.

---

## E. Numbers

| Block | Source | Verdict |
|---|---|---|
| Table 4 (Monte Carlo, model holds) | P3 §8.1 at `44db8f1`, with the rerun table | **OK**: all 48 coverage/width cells, `q` and DD agree, including the corrected 0.251 and 0.378. One caption issue: "`q = dim β`" is wrong for design B, where the Monte Carlo fitted all three columns (dim β = 18) and `q = rank Σ̂ = 12`. Fix: "`q` is the number of free coefficients (dim β on the C−1 free columns)". |
| Table 5 (R2) | P3 §8.2 | OK: every cell agrees. |
| Items 1–8 of §8.3 | P3 §8.4 | OK. Recomputed: CLR ≥ 0.994 where Slater holds (min 0.994, design C); critical values √χ²: 2.80 (q=3), 4.59 (12), 14.8 (187), 21.2 (401); the B̃-versus-B narrowing ranges 0.06–0.10 (3×3) and 0.06–0.19 (2×2 M1); the Bhatia–Davis allowance gives 0.086 → 0.057. |
| Table 6 (district shocks) | P3 §8.3 | OK: every cell agrees. |
| Test sizes | P3 §5 table; `benchmarks/results/*_size_output.txt` | OK: 0.078/0.072/**0.043** (pooling, 2000 replications); 0.057/0.055/**0.054** (curvature, 1000); 0.090 and mean 39.7, variance 100 from the 1000-replication study. |
| NZ, per electorate | `ei_datasets_nz_summary.md` | OK: 492 tables, 4302 cells, median 64 units; CLR width 0.833 against DD 0.843; M0 and M0+order cover every cell; L6 0.939, 2.2%, 0.415; plug-in empty 48% (0.476). |
| NZ, electorates as districts | same, Part 2 | **Numbers OK, wording to fix.** The figures are right: 3398–6550 units, 366–408 restrictions, emptiness test rejecting M0 in 2002 and 2020 and M1 never, 21% rejected, 8% (5/63) missed, projection never missing and as wide as DD. But the sentence *"the data reject the model: the M1 plug-in estimate is empty in all seven elections, and the curvature and pooling tests have p≈0"* presents the empty plug-in as part of the rejection. Corrected source (P3 §8.4 item 10; report §10): the M1 test rejects M1 and the pooling test rejects M0 in all seven elections (p≈0); the empty plug-in is a sample diagnostic only. Suggested text: *"Pooling each election's 62–69 electorates (3398–6550 polling places) and fitting M1 with electorates as districts, the curvature test rejects M1 and the pooling test rejects M0 in all seven elections (p≈0; the pooling test has 366–408 restrictions). The M1 plug-in estimate is empty in all seven elections, but that is a sample diagnostic: the emptiness test rejects the pooled model in two elections and M1 in none."* |
| CR table (Table 7) | 5a §3.4 | OK |
| 304/600, 0.917, 0.990, 0.230 | 5a §3.5 item 4; §5 gated run | OK. Optionally write "550/600" next to 0.917. |
| False-belief table (398/110/0 of 400 at 600 units) | 5a §3.4 (V12) | OK: every cell agrees. |
| Sign-module frequencies (83%/0.70, 89%/0.67, 80% of 1060/0.69; 8/28/64%) | PD §2.2; 5a §2.4 T4 | OK |
| Readings: 46/400; 20 strict 3×3 cases; 75/139; 4.8×10⁻⁶; 0.68; −1/√π; η* = 0.060 | 2a, 2c, P3 | OK |
| §3.6 illustration; §4.5 ladder tables | 2a §8; 2b §6 | OK (the paper's 3×3 ladder table drops the M0, M1(.10) and M1(.02) columns, which is fine) |
| **voter_vax (Table 8 and the text)** | **rerun**: eiPartialID 0.3.1 from `44db8f1`, `voter_vax_MonotoneEI.csv`, `weights = "n"`, seed 1 | **OK: every number reproduces exactly.** M0: estimates [−0.486, −0.001] and [0.345, 0.597]; CLR B̃ [−0.571, 0.056] and [0.301, 0.623]; projection B̃ [−0.607, 0.080] and [0.282, 0.634]; conservative B̃ [−0.599, 0.112] and [0.291, 0.651]; CLR B [−0.750, 0.233] and [0.229, 0.694]; projection B [−0.763, 0.257] and [0.210, 0.705]; DD. M0+CR: every cell as printed. Descriptive: μ = 0.482, variance 0.0267, γ = 0.107, ER −0.470, NM −0.050. The four single-sign sets and the CR set reproduce; δ_B ∈ [−0.488, 0.055] and δ_W ∈ [−0.525, 0.018]; the CR gain is 0.066 added first and last; n_eff is 271 and 188; the joint δ_W≥0, δ_B≥0 fit is empty, the emptiness test rejects, and the rung alone is not rejected. |

**Worked-example code (should-fix):** the displayed call does not reproduce Table 8 as written.
- It has one target (`D`), while the table also reports the Republican rate.
- It omits `seed = 1`, which the text mentions.
- The M0+CR columns need `restrictions = list(cr = list(list(row = 1, col = 1)))`, which is not shown.

Fix: show `targets = list(D = matrix(c(1,-1,0,0),2,2), R = matrix(c(1,0,0,0),2,2))` and `seed = 1`, and give the restrictions argument for the CR columns.

---

## F. Package section (§7)

Checked against `R/rxc_api.R` (`eiBounds` signature, `.rxc_ci_status`, `print.eiRxC`, the notes), `R/rxc_infer.R`, `R/rxc_modules.R`, `NEWS.md` 0.3.1 and `PHASE4_implementation.md`.

- **Defaults: OK.** These match the signature: `model = "M0"`; `range = "hull"`; `weights` "none"/"n"; `vcov = "HC3"`; `ci = "clr"` with `ndraw = 100`; `alpha = 0.05`; `fpc = "bhatia-davis"`; `union = "pieces"`; `jkst = FALSE`; `tests = TRUE`; bootstrap `B = 399`.
- **Projection: OK.** `q = dim β` (`sqrt(qchisq(1-alpha, q))`, with `q = ncol(V)`). The rank is printed as a diagnostic; the rank guard is labelled a safeguard. For `proj_actual`, `m̄` comes from a projection interval on the rung alone.
- **Status line: OK.** The paper's three cases match `.rxc_ci_status` word for word in substance: pooled (P3.5(b2)(ii), or P3.5(a) when conservative; R10 with sign rows); grouped without sign rows (regime R1 only; otherwise a simulation finding); grouped with sign rows (no coverage theorem).
- **Notes: OK.** The CLR labels, the Bhatia–Davis-from-CLR label, the CR branch-union label, the EXPERIMENTAL module line, the grouped-sign "outer set" line, the pooled "sharp for the regression" line, Slater slack < 0.005, median group size < 30, unbounded estimates, the [0, .05]/[.95, 1] Bhatia–Davis note, and the joint-rejection notes all match the code.
- **One inaccuracy (should-fix).** "What a fit reports" says that `proj_model` and `proj_actual` are always returned. The rerun shows they are `NA` when `ci = "conservative"`: the projection intervals use the CLR candidates and are computed only with `ci = "clr"`. Fix: add "(with `ci = "clr"`)".
- **Not implemented: OK.** The list agrees with NEWS and PHASE5a §5. The drop-one display and per-group W2 are "not yet" in NEWS; the paper does not claim them.

---

## G. Readability

- **Aims-first structure:** it holds throughout.
- **Undefined or clashing symbols:**
  - `D` is used for districts in §§2, 4 and 6 and for a gap in §5. This is declared in §5 ("In this section `D` without a unit index denotes a gap"), but §8.6 (worked example) uses `D` for the gap again outside §5. Optional: write "the gap `D = B₁ − B₂`" there, or use `Δ_B`.
  - `ĉ` in Thm 6(a) and `ĉ_p` in (b1) are defined. In Prop 31, `α_est` is defined only in words.
  - `S_1` is defined before Thm 6. `t_E`, `h_j` and `h_0` are defined before Thm 7.
  - The reading of Thm 7 has an awkward clause, "(about 0.01 was the truth, against ¼)". Suggest "(the true variance was about 0.01, against the worst case ¼)".
- **Paper versus supplement notation:** consistent (`n`, `ω`, `ν`, `ĉ`, `Ψ`, (I1)/(I2), `supp_A`). The supplement's Lemma S1 ("agreement on the restricted support") supports the A5 wording.
- **Internal files are out of date:**
  - `OUTLINE.md` says results "share one counter". The preamble uses separate counters, and `SOURCE_MAP.md` says so.
  - `SOURCE_MAP.md` and `ISSUES.md` cite tag `v0.3.1` (`c02954e`) as the source.
  - `ISSUES.md` items 1–4 are now resolved in `44db8f1`.

  These files are not part of the paper.

---

## H. Compilation

Command: `build.sh` on a fresh copy, with TeX Live from Ubuntu apt.
- **Result:** "build ok". The paper is 52 pages and the supplement 25, the same as the committed PDFs.
- **Undefined references or citations:** none (0 "undefined" in either log). The `xr-hyper` cross-references from the supplement to the paper (`P-…`) resolve.
- **Overfull boxes:** none in either log. The text extracted from both PDFs has no "??".
- **Warnings:** only hyperref "Token not allowed in a PDF string (Unicode)" (3 in each log), from math in section titles going into PDF bookmarks. They are harmless; optionally use `\texorpdfstring` for R$\times$C.

---

## I. Prioritised fix list

### Must-fix
None. No result is misstated, no proof is missing, there is no process language, and the paper compiles.

### Should-fix
1. **§8.4, NZ electorates as districts:** reword so that the rejection comes from the curvature test (M1) and the pooling test (M0), and the empty M1 plug-in is called a sample diagnostic only. Suggested text is under §E. Corrected source: P3 §8.4 item 10.
2. **§7, "What a fit reports":** `proj_model` and `proj_actual` are computed only with `ci = "clr"`; they are `NA` with `ci = "conservative"`.
3. **§8.6 worked example:** make the displayed call reproduce Table 8: add the Republican-rate target, `seed = 1`, and the `restrictions = list(cr = …)` call for the M0+CR columns.
4. **Table 4 caption:** "`q = dim β`" is wrong for design B (rank 12 of 18). Use "the number of free coefficients (dim β on the C−1 free columns)".
5. **Internal files:**
   - update `SOURCE_MAP.md` and `ISSUES.md` to cite `44db8f1`;
   - mark `ISSUES.md` items 1–4 as resolved in the sources;
   - correct the counter sentence in `OUTLINE.md`.

### Optional
6. **Prop 5:** add "which also requires `Q_{c,rr} ≥ 0`" to the displayed set.
7. **Prop 10:** "Under A5 on `supp_A P_X` (with `A ⊇ supp P_X`)".
8. **Prop 13(b):** add the second sufficient condition (`D⊥N | X` with `E[Z_D | X]` affine), and add that (S3) is automatic only on the population hull.
9. **Supplement, proof of Prop 28(d):** add the step "differentiable and non-decreasing on a neighbourhood ⇒ `m'(x) ≥ 0`".
10. **Supplement, proof of Thm 8:** renumber the steps within each part.
11. **§8.4, provenance sentence:** reword to drop "before the later changes" (text under §D).
12. **Supplement S1 citation table:** write "Billingsley (1995), Theorem 27.2" to match the corrected source.
13. **Prop 25 reading:** mention R5's "lower-dimensional support" sub-item.
14. **§8.5:** add the count 550/600 next to 0.917.
15. **Theorem 7 reading:** fix the clause about the true variance (text under §G).
16. **Worked example (§8.6):** avoid the bare `D` for the gap outside §5 (§G).
17. **Section titles:** use `\texorpdfstring` for R$\times$C to silence the hyperref warnings.

### Counts
- Numbered results checked: **59**, plus A1–A10, Definitions 1–2, Remarks 1–5 and the unnumbered paragraphs listed in §A.
- Mismatches: **0 substantive, 3 minor** (Props 5, 10, 13), plus 1 wording issue in the evidence section (NZ) and 1 table caption.
- Process-language hits: **0** (1 borderline software-provenance sentence).
- Proofs missing or truncated: **0** (2 small gaps for clarity).
- Build: OK (52 + 25 pages; 0 undefined references; 0 overfull boxes).
