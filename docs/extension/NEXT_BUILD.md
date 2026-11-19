# What to build next in eiPartialID

*Drafted 18 Nov 2026, after the review-response package changes (WS5) and the merge of the M-SIGN branch. This combines three sources: the build ladder in `PACKAGE_DESIGN.md` §6, the "second release" list in `PHASE4_implementation.md`, and the open items of the review response.*

## Where the package stands

* **Rungs:** DD, M0, M2, M1(δ), M1, L4, L5, L6 (Phases 2a–2b), plus restriction rows and M-SIGN sign/CR restrictions with the redundancy diagnostic (PHASE5a).
* **Inference:** CLR-type intervals by default (algorithmic), with projection intervals alongside (proved), and the conservative interval as an option. Tests: emptiness, M1's testable part, pooling (M0 vs M1). Sensitivity (η, slack) and breakdown.
* **Verification:** full test suite passes (0 failures); `R CMD check` has no notes, and only a container locale warning.
* **Version:** **0.3.1**, merged to `master` (Nov 2026; mdweaver/eiid#3). 0.3.0 was the Batch B response; 0.3.1 the round-2 response (per-fit inference status in `print()`, $q=\dim\beta$, `proj_actual` per P3.5(c1)). M-SIGN is flagged experimental.
* **Review:** both rounds answered for all three referees (`REVIEW_RESPONSE.md`, `REVIEW_RESPONSE_ROUND2.md`, `reviews/`). Round-2 verdicts: favourable (Phases 1–3), minor revision (both Batch B referees).
* **Technical report:** `docs/report/eiPartialID_report.pdf` (150 pp.), with a completeness check against the repo.

## Recommended order

### Now: finish the review response (no new features)
1. **WS6, Monte Carlo rerun — done (Nov 2026).** 1000 replications per design with MC standard errors; riding shocks, Slater failure, R2 up to 200 districts, projection everywhere; NZ benchmark rerun with M1 by electorate (P3 §8, `REVIEW_RESPONSE_PLAN.md` WS6 status). Items it raises for later:
   * the pooling test over-rejects (9% at 5% with 38 restrictions). **Done in 0.3.0:** both tests use a restricted wild bootstrap (0.043 and 0.054 at 5%);
   * `proj_actual` took its Bhatia–Davis variance bound from a CLR interval. **Done in 0.3.0:** projection-based for `proj_actual`, conservative-based for conservative `ci_actual`;
   * on the NZ data, M0 and M1 by electorate are both rejected; a non-empty CLR interval missed the truth in 8% of cells of the pooled M1 fit.

2. **WS7, response letter — done** (`REVIEW_RESPONSE.md`), to go with batch A of `EXTERNAL_REVIEW_PLAN.md`. It now includes three core corrections found by the Batch B referee (2a support of the restricted measure; projection under a singular covariance; residual CLT).
2c. **Round 2 (Nov 2026), theory side — done** (`reviews/ROUND2_THEORY_PLAN.md`): P3.4 restated for a fixed constraint system and new Proposition P3.4′ (triangular arrays); new Proposition P3.5(c) for the Bhatia–Davis allowance ((c1) valid with a projection or conservative interval for $\bar m$ on the rung alone; (c2) with CLR, algorithmic); the projection theorem is claimed with $q=\dim\beta$ (numerical rank cutoff = safeguard); 2b Prop. D8 split into cancellation and given-the-set identification; P3 §0.4 independence as a substantive assumption; 2c/5a round-2 wording. P3 §8's projection column for $B$ rerun with the P3.5(c1) bound. **Package (0.3.x), for the package session:** $q=\dim\beta$; $\bar m$ projection interval on the rung alone for `proj_actual`, with a note near 0 or 1 (Popoviciu fallback later); per-fit CI status in `print()`.
2b. **WS-B8, responses to the two Batch B referees — done** (`reviews/B_partial_identification_response.md`, `reviews/B_ecological_inference_response.md`). The corrections themselves are done (WS-B0 to WS-B3, Nov 2026).
3. **Release hygiene — done (0.3.0):** version and `NEWS.md`; M-SIGN flagged `experimental`; `MASS` removed and `eco` moved to Suggests (the NOTE is gone). Pull request to `master` merged (0.3.1); version tag `v0.3.1` to follow.

### Next: a process-free paper version
3b. **Paper version of the technical report** (on hold until the owner starts it): a journal-style methods paper, 35–50 pp., with a rewritten proof appendix; covers the core bounds and inference, the monotone results and M-SIGN (experimental), the package and the evidence; no DiD; no reference to the review process; conjectures and open items labelled and collected in a closing section.

### Then: the first new module
4. **M-SHARP, the series-based sharper bound** (PACKAGE_DESIGN step 1b; also first on PHASE4's agreed second-release list).
   * It estimates the tighter set $\Theta_X$ by series regression of the unit-level Duncan–Davis limits, for pooled designs. It conditions on composition only, so when precinct size is informative it is "tighter than the outer set", not sharp (2c Remark N2.2.1; Conjecture N2.2′ is left open by decision, pending another reviewer's work on the rungs).
   * It tightens every rung with **no new assumption**. It also answers the reviewer's point (2a §4.3) that the sharp set needs smoothing or pooling to be estimable.
   * Proof item: R9 (consistency; never wider than the current outer set).
   * The M-SIGN session's queued work can run alongside, on this branch: real-data runs (scripts only) and the slack documentation.

### Then, in order of proof readiness (PACKAGE_DESIGN §6)
5. **M-LOC**, weighted-local contextual monotonicity, in 2×2 and 2×C (step 2; needs conjecture C1, R1). Not blocked by Batch B, but P-N is now scoped to the pooled model (M0) or within a district.
6. **2×C joint column signs** and ordinal cumulative columns (step 3; R6).
7. **R×2 pairwise ordering plus local restrictions on spanning trees**, with hull diagnostics (step 4; R4, R5).
8. **DiD wrappers**: the Δ-rungs and FD (PHASE4 second-release list). Gated on batch C of the review plan and the Wisconsin/Iowa data decisions. Δ0 and Δ4 are sharp for continuous tables (integer tables: R15).
9. **Quadratic contextual effects (L1′)** as a weaker rung, and the optional non-negative-correlation assumption for the unit-level term (both on PHASE4's list).

### Later or research-dependent
* **From the Batch B review:** grouped-rung proofs for M-SIGN (R16, needed to drop "experimental" under M1/M2); integer unit tables (R15); per-group W2 restrictions (U8(a)); the drop-one display and the L5/L6 sign check; the N2 global LP only (local rows dropped, U4(a)).
* PACKAGE_DESIGN steps 5–11: shared $K$ across validation tables (E4/E5), B1(M), κ ceilings, N2 grid LP, multiple contests, CAR-δ, dispersion. These are mostly data-dependent.
* **Covariate adjustment** (Phase 2c; review-response WS8, low priority by decision).
* **Informative finite-sample inference for many small groups.** Theorem P3.8 is valid but uninformative; this is open (P3 §10).
* **Closing Conjecture P3.5(b3)**, which would make the CLR default a theorem.

## Decision for you

The main choice is whether **M-SHARP** (my recommendation: no new assumption, fixes a reviewer concern, and was already agreed for the second release) comes before **M-LOC** (PACKAGE_DESIGN's step 2). Either way, items 1–3 come first.
