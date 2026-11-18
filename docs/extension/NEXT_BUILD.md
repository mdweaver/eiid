# What to build next in eiPartialID

*Drafted 18 Nov 2026, after the review-response package changes (WS5) and the merge of the M-SIGN branch. This combines three sources: the build ladder in `PACKAGE_DESIGN.md` §6, the "second release" list in `PHASE4_implementation.md`, and the open items of the review response.*

## Where the package stands

* **Rungs:** DD, M0, M2, M1(δ), M1, L4, L5, L6 (Phases 2a–2b), plus restriction rows and M-SIGN sign/CR restrictions with the redundancy diagnostic (PHASE5a).
* **Inference:** CLR-type intervals by default (algorithmic), with projection intervals alongside (proved), and the conservative interval as an option. Tests: emptiness, M1's testable part, pooling (M0 vs M1). Sensitivity (η, slack) and breakdown.
* **Verification:** 439 tests pass; `R CMD check` is clean apart from pre-existing notes.
* **Version:** still labelled 0.2.0, although M-SIGN and the WS5 changes are substantial. Before any release it should become **0.3.0**, with a `NEWS.md`.

## Recommended order

### Now: finish the review response (no new features)
1. **WS6, Monte Carlo rerun — done (Nov 2026).** 1000 replications per design with MC standard errors; riding shocks, Slater failure, R2 up to 200 districts, projection everywhere; NZ benchmark rerun with M1 by electorate (P3 §8, `REVIEW_RESPONSE_PLAN.md` WS6 status). Items it raises for later:
   * the pooling test over-rejects (9% at 5% with 38 restrictions): consider the F form or a wild bootstrap;
   * `proj_actual` takes its Bhatia–Davis variance bound from a CLR interval: use the projection interval (or Popoviciu) there so the projection column is fully covered by P3.5(b2);
   * on the NZ data, M0 and M1 by electorate are both rejected; a non-empty CLR interval missed the truth in 8% of cells of the pooled M1 fit.

2. **WS7, response letter** (`REVIEW_RESPONSE.md`), with batch A of `EXTERNAL_REVIEW_PLAN.md`.
3. **Release hygiene:**
   * version 0.3.0 and `NEWS.md`;
   * keep M-SIGN modules flagged `experimental` until batch B is reviewed;
   * resolve the long-standing NOTE about unused imports (`MASS`, `eco` are used only in examples).

### Next: the first new module
4. **M-SHARP, the series-based sharper bound** (PACKAGE_DESIGN step 1b; also first on PHASE4's agreed second-release list).
   * It estimates the tighter set $\Theta_X$ by series regression of the unit-level Duncan–Davis limits, for pooled designs.
   * It tightens every rung with **no new assumption**. It also answers the reviewer's point (2a §4.3) that the sharp set needs smoothing or pooling to be estimable.
   * Proof item: R9 (consistency; never wider than the current outer set).
   * The M-SIGN session's queued work can run alongside, on this branch: real-data runs (scripts only) and the slack documentation.

### Then, in order of proof readiness (PACKAGE_DESIGN §6)
5. **M-LOC**, weighted-local contextual monotonicity, in 2×2 and 2×C (step 2; needs conjecture C1, R1).
6. **2×C joint column signs** and ordinal cumulative columns (step 3; R6).
7. **R×2 pairwise ordering plus local restrictions on spanning trees**, with hull diagnostics (step 4; R4, R5).
8. **DiD wrappers**: the Δ-rungs and FD (PHASE4 second-release list). Gated on batch C of the review plan and the Wisconsin/Iowa data decisions.
9. **Quadratic contextual effects (L1′)** as a weaker rung, and the optional non-negative-correlation assumption for the unit-level term (both on PHASE4's list).

### Later or research-dependent
* PACKAGE_DESIGN steps 5–11: shared $K$ across validation tables (E4/E5), B1(M), κ ceilings, N2 grid LP, multiple contests, CAR-δ, dispersion. These are mostly data-dependent.
* **Covariate adjustment** (Phase 2c; review-response WS8, low priority by decision).
* **Informative finite-sample inference for many small groups.** Theorem P3.8 is valid but uninformative; this is open (P3 §10).
* **Closing Conjecture P3.5(b3)**, which would make the CLR default a theorem.

## Decision for you

The main choice is whether **M-SHARP** (my recommendation: no new assumption, fixes a reviewer concern, and was already agreed for the second release) comes before **M-LOC** (PACKAGE_DESIGN's step 2). Either way, items 1–3 come first.
