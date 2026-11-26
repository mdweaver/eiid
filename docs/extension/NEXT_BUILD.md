# What to build next in eiPartialID

*Drafted 18 Nov 2026, after the review-response package changes (WS5) and the merge of the M-SIGN branch. This combines three sources: the build ladder in `PACKAGE_DESIGN.md` §6, the "second release" list in `PHASE4_implementation.md`, and the open items of the review response.*

## Where the package stands

* **Rungs:** DD, M0, M2, M1(δ), M1, L4, L5, L6 (Phases 2a–2b), plus restriction rows and M-SIGN sign/CR restrictions with the redundancy diagnostic (PHASE5a).
* **Inference:** CLR-type intervals by default (algorithmic), with projection intervals alongside (proved), and the conservative interval as an option. Tests: emptiness, M1's testable part, pooling (M0 vs M1). Sensitivity (η, slack) and breakdown.
* **Verification:** full test suite passes (0 failures); `R CMD check` has no notes, and only a container locale warning.
* **Version:** **0.4.0** on the working branch (M-COUNT, step 1); **0.3.1** merged to `master` (Nov 2026; mdweaver/eiid#3). 0.3.0 was the Batch B response; 0.3.1 the round-2 response (per-fit inference status in `print()`, $q=\dim\beta$, `proj_actual` per P3.5(c1)). M-SIGN is flagged experimental.
* **Review:** both rounds answered for all three referees (`REVIEW_RESPONSE.md`, `REVIEW_RESPONSE_ROUND2.md`, `reviews/`). Round-2 verdicts: favourable (Phases 1–3), minor revision (both Batch B referees).
* **Technical report:** `docs/report/eiPartialID_report.pdf` (150 pp.), with a completeness check against the repo.

## Recommended order

### Now: finish the review response (no new features)
1. **WS6, Monte Carlo rerun — done (Nov 2026).** 1000 replications per design with MC standard errors; riding shocks, Slater failure, R2 up to 200 districts, projection everywhere; NZ benchmark rerun with M1 by electorate (P3 §8, `REVIEW_RESPONSE_PLAN.md` WS6 status). Items it raises for later:
   * the pooling test over-rejects (9% at 5% with 38 restrictions). **Done in 0.3.0:** both tests use a restricted wild bootstrap (0.043 and 0.054 at 5%);
   * `proj_actual` took its Bhatia–Davis variance bound from a CLR interval. **Done in 0.3.0:** projection-based for `proj_actual`, conservative-based for conservative `ci_actual`;
   * on the NZ data, M0 (pooling test) and M1 (test against M2) by electorate are both rejected, though the emptiness test does not reject M1; a non-empty CLR interval missed the truth in 8% of cells of the pooled M1 fit.

2. **WS7, response letter — done** (`REVIEW_RESPONSE.md`), to go with batch A of `EXTERNAL_REVIEW_PLAN.md`. It now includes three core corrections found by the Batch B referee (2a support of the restricted measure; projection under a singular covariance; residual CLT).
2c. **Round 2 (Nov 2026), theory side — done** (`reviews/ROUND2_THEORY_PLAN.md`): P3.4 restated for a fixed constraint system and new Proposition P3.4′ (triangular arrays); new Proposition P3.5(c) for the Bhatia–Davis allowance ((c1) valid with a projection or conservative interval for $\bar m$ on the rung alone; (c2) with CLR, algorithmic); the projection theorem is claimed with $q=\dim\beta$ (numerical rank cutoff = safeguard); 2b Prop. D8 split into cancellation and given-the-set identification; P3 §0.4 independence as a substantive assumption; 2c/5a round-2 wording. P3 §8's projection column for $B$ rerun with the P3.5(c1) bound. **Package (0.3.x), for the package session:** $q=\dim\beta$; $\bar m$ projection interval on the rung alone for `proj_actual`, with a note near 0 or 1 (Popoviciu fallback later); per-fit CI status in `print()`.
2b. **WS-B8, responses to the two Batch B referees — done** (`reviews/B_partial_identification_response.md`, `reviews/B_ecological_inference_response.md`). The corrections themselves are done (WS-B0 to WS-B3, Nov 2026).
3. **Release hygiene — done (0.3.0):** version and `NEWS.md`; M-SIGN flagged `experimental`; `MASS` removed and `eco` moved to Suggests (the NOTE is gone). Pull request to `master` merged (0.3.1); version tag `v0.3.1` to follow.

### Done since: technical report and paper
3b. **Technical report** (`docs/report/`) and **methods paper** (`docs/paper/`): done and merged (Nov 2026).

### Build order from Nov 2026 (supersedes the earlier items 4–9)
*Decided by the owner on 21 Nov 2026 from `next_steps/NEXT_STEPS_RECOMMENDATION.md`, which has the plain-language explanation, worked example, ordering reason and technical notes for every step, the combination matrix (Part 3A) and the avenue-by-avenue review of `NEW_AVENUES.md` (Part 4). **Scope rule:** this repository develops general theory and methods; application-specific work (data assembly, calibration, covariate choice) is outside it (recommendation Part 9). **Review rule:** every step gets a proof review and a separate code review.*

0. **Housekeeping:** printout fixes; the planned batches D and E become the proof and code reviews of existing modules.
1. **M-COUNT / `eiCounts()`**: exact bounds on actual counts and rates (Duncan–Davis first), realised ordering, sign and CR rows, known cells and published totals, integer check ("auto"; off for non-integer margins), witness tables. `eiBounds()` prints it first (`target = "model"` restores the old order). Combination rule: count rows and model rows never share one LP; the actual-rate answer intersects the two, and disjointness is a test of the model. Proof items R19, R27: **theory done and proof reviewed** (`PHASE6a_count_layer.md`, `reviews/PHASE6a_proof_review.md`; `verify_phase6a.py` 36/36). Note for the code: the package's `duncan_davis` column is outer for combination targets (PHASE6a §9 item 3). **Code done and code reviewed (Nov 2026; 0.4.0):** `eiCounts()`, `witness()`, the count block in `eiBounds()` (actual-rate block first; `target = "model"` unchanged from 0.3.1), certified GLPK integer programs with a time limit and an "unresolved" status, exact rational rounding (`reviews/PHASE6a_code_review.md`, three rounds; 1897 tests pass). Suggested follow-up: keep a sample of the review oracle instances as a periodic regression check of the GLPK integer path.
1b. *Moved to the end of the plan (owner, 22 Nov 2026): see "Uncertain totals and compositions" after step 11.*
2. **DiD theory track** (general; T ≥ 2 with pre-periods): descriptive and causal (trend dial, bounds relative to set-valued pre-trends, MTR/ND/PS), with a register of where changing populations need assumptions or data-based bounds. Implementation (`eiChange()`, `eiATT()`) after its proof review. R11a–d, R11b′. **Theory done and proof reviewed (Nov 2026):** `PHASE2e_did_general.md`, `reviews/PHASE2e_proof_review.md` (two rounds; 19 issues and four residual items fixed), `verify_phase2e.py` 76/76; register rows R11a–d, R11b′, R31 (gap register), R32 (margin sets over elections); errata to 2d, 1b, PD §7 and NS from 2e §11 applied. **Code: in progress.**
   *First DiD release scope (owner, 26 Nov 2026).* Include: the 9 code-review fixes; within-precinct parallel trends in expectation (no default; report $n_{\rm eff}$, $\pi_U$, local comparison counts); RM on increments; joint baseline/change LCE with the baseline-curvature setting (E30); the local parallel-trends dial (E33(a)); the Hoeffding allowance (E29). Either setting can be varied with the other held fixed, each with a one-LP breakdown value. **Reporting rule:** every printed interval says whether it is an identification bound under the mean model or a confidence interval from one observed election pair; E29 covers the realised-count route only and does not make the sharper LCE interval a confidence interval. Deferred, with reasons: the variance-based allowance (E31(b) is cheap and can help, e.g. 45.9 → 30.2 points with 500 precincts and a rare outcome, but its gain on the owner's data is unestablished); the moment-projection confidence set (E32; keep the prototype; very wide in the one synthetic design tested, still including zero at 5,000 precincts); stratum relaxations and the fast per-precinct RM method (E33(b)–(c), E34; not needed for the local dial).
3. **M-CONTRAST**: finite comparisons between prespecified kinds of places, in levels and changes; M-LOC becomes its `at = "hull"` option. R22.
4. **B1(M) + κ**: context at a second scale, per column, with κ-annotated endpoints. R8.
5. **M1(δ) with named neighbours**: adjacency or distance matrix (Lipschitz default δ_de = λ·d_de). R30.
6. **CAR(Z)-δ**, then **covariate-benchmarked κ/δ** with three diagnostics (benchmark size, composition R², overlap factor). R28. **Widened (owner, 26 Nov 2026) to covariates in both ordinary EI and DiD:** conditional parallel trends (PHASE2e gap item G8) and the age adjustments proposed in `reviews/PHASE2e_independent_review.md` (option 1, an age-gap sensitivity slope; option 2, parallel trends within precinct × age band), designed once for both layers.
7. **M-EXCEPT**: exception-mass frontiers. R23.
8. **M-SHARP diagnostic** (kink mass, (X,N) conditioning); the estimator only if polarised precincts are found. R24, R9′.
9. **M-BERN** polynomial robustness ladder, with L1′.
10. **C2, contests as districts**, with its overidentification test. R29.
11. *E4 table transfer: moved to the very end of the plan (owner, 22 Nov 2026); see below.*

General-user track: M-AUDIT, M-CVR (exact published-table rows ship inside step 1).

**After step 10: uncertain totals and compositions** (owner's decision, 22 Nov 2026; replaces step 1b). One general module for margins and compositions that are not known exactly but lie in a set with stated confidence: match probabilities (the former step 1b, R25/R26: Poisson-binomial margin sets, calibration slack, the Berkson correction in the LCE layer), survey-based margins, and changes in composition between elections (entry, exit, aging-in, from the DiD gap register). Placed late so that the needs found along the way (above all in the DiD theory's gap register) feed its specification. Until then, non-integer margins switch the integer check off and `print()` says the bounds are exact only for the stated margins.

**Then: open items from M-COUNT** (owner's decision, 22 Nov 2026; `PHASE6a_count_layer.md`). Left open until the end of this plan, in this order of priority:
1. **A smaller valid set than the intersection** (R27(iii-b); computed example X5 only). Use the count rows when computing the model-to-actual allowance, so the combined answer can be tighter than count set ∩ model interval. Highest priority of the three: the owner considers it the most useful.
2. **Integer theory outside the TU classes** (R15 remainder; PHASE6a §3.4). A theorem bounding the integer–continuous gap for pairwise orderings with R ≥ 3, W2 rows, totals across rows and columns, and mixtures. Lower priority: the certified MILP already gives exact answers.
3. **Bhatia–Davis allowance for combined targets** (PHASE6a §9 item 7; P3.5(c1) covers single cells only). Lower priority: Popoviciu is valid, and the difference matters mainly with rates near 0 or 1 and small samples.

**Last: E4 table transfer** (NZ leave-one-out; owner's decision, 22 Nov 2026). Moved after everything else, because the owner expects tables of this kind (a known cross-tabulation from a comparable setting) to exist in too few applications for the tool to be widely useful. Build only if such tables turn up.

### Later or research-dependent
* **From the Batch B review:** grouped-rung proofs for M-SIGN (R16, needed to drop "experimental" under M1/M2); integer unit tables (R15); per-group W2 restrictions (U8(a)); the drop-one display and the L5/L6 sign check; the N2 global LP only (local rows dropped, U4(a)).
* Research only (recommendation Part 4): Avenue A (unit size, with the dispersion ceiling and the covariance route), single index (R14), KL-ball reference dial.
* **Covariate adjustment** (Phase 2c; review-response WS8, low priority by decision).
* **Informative finite-sample inference for many small groups.** Theorem P3.8 is valid but uninformative; this is open (P3 §10).
* **Closing Conjecture P3.5(b3)**, which would make the CLR default a theorem.

## Decisions

All settled (recommendation Part 6, D1–D12); D9 (match calibration) and D10 (covariate set) are user inputs, outside this repository.
