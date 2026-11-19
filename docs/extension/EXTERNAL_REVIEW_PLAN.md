# Plan: external review of the unreviewed material

*Drafted 18 Nov 2026. The external review so far covers Phases 1, 2a, 2b and 3 (two reports; response in `REVIEW_RESPONSE_PLAN.md`). This plan covers everything else in `docs/extension/` and the package code that has had no external review. Open consistency items are tracked in `MSIGN_REVIEW_LOG.md`; this plan is about substantive review.*

## 1. What is and is not reviewed

| Material | Review so far | Used by the package? |
|---|---|---|
| P1, 2a, 2b, P3 | external (two reports); corrections applied (WS1–WS4) | yes (core) |
| **New results added in the response** (listed in batch A) | none | yes |
| `PHASE5a_sign_module.md` §§1–3 (R3, R2, R10) | fresh-session adversarial review only (5a §4) | **yes** (M-SIGN, exported) |
| `PHASE2c_monotone.md` (M1–M5, N2.1–N2.3) | none | **yes** (the sign rows rest on M5) |
| `PHASE1b_monotone_did.md` (proposal), `PHASE2d_did.md` (DiD derivations) | none | not yet (DiD wrappers planned) |
| `PACKAGE_DESIGN.md` (architecture, module rules, proof register R1–R14) | none | governs what is built next |
| `NEW_AVENUES.md` | none (research agenda) | only through results quoted elsewhere |
| Package code: `rxc_sign.R`, `rxc_modules.R`, M-SIGN and WS5 changes to `rxc_core.R`/`rxc_infer.R`/`rxc_api.R` | none (tests only: 439 pass) | — |

## 2. Review batches

Each batch is sent as a frozen commit (a git tag, e.g. `review-B-2026-12`). It includes the documents, their verify scripts and outputs, and a one-page cover note listing the specific questions below. Reviewers are asked for the same format as the first two reports: a verdict per result (established / established under conditions / incorrect / not established), with major and minor comments.

### Batch A: new results introduced by the review response (return to the original reviewer)

**Status (Nov 2026): ready to send.** The response letter is `REVIEW_RESPONSE.md`. It includes the three core corrections found in Batch B. **Round 2 (Nov 2026):** the referee's second report (`reviews/A_phases1to3_report_round2.md`, favourable, targeted revisions) is answered in `REVIEW_RESPONSE_ROUND2.md` (P3.4 fixed-constraint + P3.4′; P3.5(c); D8 split; independence as an assumption).

*Send with the response letter (WS7). The reviewer knows the context, so this is the cheapest batch and closes the current review.*

* **2a:** Remarks 5.1–5.2 (S3 and dimension counterexamples); Remark 4.2 (continuous vs integer rates; part 2 is a sketch); §4.3; **Prop. 7′** (units outside $A$); Theorem 1 facts F1–F7; the §0.4 citation table.
* **2b:** **Lemma D1(d)** (stacked rank under M1); Prop. D4(b) as corrected, and Remark D4.1; D3 with pure-composition vertices.
* **P3:**
  * §0.4 (model-based randomness, no clustering) and §0.5 (regime table);
  * **P3.5(b2)** (weak-duality proof) and the reduction in **Conjecture P3.5(b3)**;
  * **P3.6(c)** (pooling test, conditional on the design);
  * **Theorem P3.8** (finite-sample CI: valid but uninformative);
  * P3.2's non-nested hull argument.
* **Questions:**
  1. Is the D3/clustering argument (P3 §0.4) acceptable as the answer to the reviewer's §4.5?
  2. Is P3.6(c) right that D4(b)'s weaker condition is a superpopulation statement?
  3. Can the (b3) conjecture be closed, or is (b2)(ii) the right formal default?
  4. Is the Remark 4.2(2) sketch worth completing?

### Batch B: M-SIGN theory (highest priority among the unreviewed)

**Status (Nov 2026): reviewed.** Two external reports (partial identification; ecological inference), both major revision. Corrections were applied in theory and code (0.3.0); see `reviews/B_PI_RESPONSE_PLAN.md`. Responses sent (WS-B8). **Round 2 (Nov 2026):** both referees recommend minor revision; reports `reviews/B_*_report_round2.md`, answered in `reviews/B_*_response_round2.md` (package 0.3.1; theory plan `reviews/ROUND2_THEORY_PLAN.md`, package plan `reviews/ROUND2_PACKAGE_PLAN.md`). **Lesson for later batches:** the internal fresh-session review of PHASE5a §4 had accepted it with fixes, and missed the integer-table, district-slope, derivative-identification and binning problems, all in 2c, which it had not been asked to cover. Every batch should give the pre-screen the full scope that the external reviewers get.

*Already shipped and exported, so it is the largest exposure.*

* **Documents:** `PHASE2c_monotone.md` (all); `PHASE5a_sign_module.md` §§1–3.
* **Results:**
  * 2c Lemma M1; Props M2–M5; Lemma N2.1; Props N2.2–N2.3;
  * 5a Theorem R3 (sharpness of LCE ∩ restrictions; CR as union or slab);
  * Prop. R2 (when an added row tightens: the face test);
  * Theorem R10 (inference with added rows, slab, union of branches, emptiness), with the status note added after the review.
* **Questions** (from `MSIGN_REVIEW_LOG.md`):
  * M-9: does R3(c) handle gaps correctly? The certificate is now stated as sufficient, not necessary.
  * M-10: the face-test proof, and the T1b addendum.
  * M-11: R10(c) step 9. The M-SIGN session says all branch rows are exact; please confirm.
  * R10(e): the persistent under-coverage of the naive bracket.
  * Does M5's mapping from sign beliefs to rows in $k$ hold for every LCE rung?
* **Verify scripts:** `verify_phase2c.py`, `verify_R2.py`, `verify_R3.py`, `verify_R10.py`.

### Batch C: difference-in-differences

*Needed before the DiD wrappers are built.*

* **Documents:** `PHASE1b_monotone_did.md` (proposal: assumptions and structure only); `PHASE2d_did.md` (derivations).
* **Results:** 2d Lemma G; Props Δ0, FD1–FD4, Δ4, Δ6, C1, C2, CL, P2, P3; the nesting display (§10).
* **Questions:**
  * Whether the Δ4 rung (monotone EI on changes) is sharp.
  * The "repeated unit across periods" treatment in place of clustering (M-5).
  * The aging-in decomposition (C2).
  * The Athey–Imbens citation flagged "to verify" in 1b.
  * The Phase-3 inference items in 2d §11.
* **Verify scripts:** `verify_phase1_mei_did.py`, `verify_phase2d.py`.

### Batch D: design and agenda

* **Documents:** `PACKAGE_DESIGN.md` (all); `NEW_AVENUES.md` (skim for results quoted elsewhere).
* **Questions:**
  * Module combination rules (§1.2, §4.2 layering matrix).
  * The lemmas P-S (sandwich), P-R (rank) and P-N (nesting), with conjecture C1.
  * G2 (§2.2) and DiD-S.
  * Whether each proof-register status (R1–R14) is justified.
  * Re-derive the two cluster-based figures in `NEW_AVENUES.md` (M-15).
* **Verify scripts:** `verify_package_design.py`, `new_avenues/*.py`.
* **Order:** after batch B, since several statuses depend on it.

### Batch E: package code (in parallel with B)

* **Scope:**
  * `rxc_sign.R`, `rxc_modules.R`;
  * the M-SIGN changes and the WS5 changes in `rxc_core.R`, `rxc_infer.R`, `rxc_api.R`;
  * `test-msign.R`, `test-msign-dial.R`, `test-review-response.R`.
* **Questions:**
  * Does the code implement PHASE5a §5 and PHASE4 as documented?
  * The $2^m$ growth of branch systems (M-16).
  * That `.rxc_ci_clr` refuses disjunctive rows.
  * The projection interval and the pooling test against P3.5(b2) and P3.6(c).
  * Are the tests adequate, including the gated coverage MC (M-17)?
* **Reviewer:** an R programmer with numerical-optimisation experience. This need not be an econometrician.

## 3. Reviewer profiles and pre-screening

* **A–D:** an econometrician working on partial identification and moment-inequality inference (A is the original reviewer). B and C benefit from familiarity with ecological inference (JKST; Elzayn et al.).
* **Pre-screen (PACKAGE_DESIGN S12).** Before a batch goes out, a fresh Claude session re-derives the proofs adversarially and reports per result, as was done for PHASE5a §4. This catches the cheap errors. It does not replace the external review: none of the review response's corrections was caught by the earlier internal checks.

## 4. Handling the reports

* Store each report as `docs/extension/reviews/<batch>_<date>.md`.
* Answer it with a ledger in the same format as `REVIEW_RESPONSE_PLAN.md` (IDs, type E/S/P/X/C, location, action, package impact, a check that fails on the old statement).
* Until a result's batch is reviewed and resolved, the package keeps `experimental` in `meta$proof_status` for modules that depend on it (PACKAGE_DESIGN §8).

## 5. Order and gates

| Order | Batch | Gate before sending | Blocks |
|---|---|---|---|
| 1 | A | WS6 (Monte Carlo rerun) and WS7 (response letter) done | closing the current review |
| 2 | B and E (in parallel) | fresh-session pre-screen of B | removing "experimental" from M-SIGN; the 0.3.0 release notes |
| 3 | C | pre-screen; DiD data decisions (1b §D) | DiD wrappers |
| 4 | D | B resolved | the build order beyond step 1 (PACKAGE_DESIGN §6) |
