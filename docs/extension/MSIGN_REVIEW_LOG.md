# Review log: material merged from `claude/msign-build`

*Started 18 Nov 2026. On that date `claude/msign-build` (head `b3f709b`, 33 commits ahead of `master`, including all of `claude/package-design`) was merged into `claude/wizardly-volta-lbh4bm`, so that the review-response package changes (WS5) could be built on the newest code. The merge had no conflicts. This log records what in the merged material still needs review, and why. It is a list of open items, not a review.*

Status labels:
* **open**: needs review;
* **fix needed**: the text contradicts a corrected result and must change;
* **checked**: reproduced in this session, which is not a full review;
* **done**: fixed in this branch, with a pointer.

## 1. What was merged

* **Docs:**
  * `PACKAGE_DESIGN.md`, `PLAIN_PACKAGE_DESIGN.md`
  * `PHASE1b_monotone_did.md`, `PHASE2c_monotone.md` (+ plain), `PHASE2d_did.md` (+ plain)
  * `PHASE5a_sign_module.md`, `NEW_AVENUES.md` and `new_avenues/`
  * additions to `ASSUMPTION_MENU.md` and `PLAIN_LANGUAGE_GUIDE.md` (§9)
* **Verify scripts:** `verify_R2.py`, `verify_R3.py`, `verify_R10.py`, `verify_package_design.py`, `verify_phase1_mei_did.py`, `verify_phase2c.py`, `verify_phase2d.py`.
* **Package:**
  * new `R/rxc_sign.R` and `R/rxc_modules.R`;
  * changes to `rxc_api.R`, `rxc_core.R`, `rxc_infer.R`;
  * new exports `eiModules`, `eiRedundancy`;
  * tests `test-msign.R` and `test-msign-dial.R` with data files;
  * a vignette section.

## 2. What was reproduced here (not a review)

| Item | Result |
|---|---|
| Package test suite after the merge (R 4.3.3) | 425 passed, 0 failed, 3 skipped (two need `ei.Datasets`; one slow MC gated by `EIPARTIALID_SLOW`) |
| `verify_phase1_mei_did.py`, `verify_phase2c.py`, `verify_phase2d.py` | 11, 13, 16 PASS; output identical to the committed files |
| `verify_R10.py`, `verify_package_design.py`, `verify_R2.py` | output identical to the committed files |
| `verify_R3.py` | see §5 |
| Running convention | `verify_R2.py` and `verify_R3.py` import `verify_package_design` from their own folder, so they **fail under `python3 -I`**, which the phase documents use. Run them without `-I`, or fix the import. **open** |

The phase documents 2a, 2b and P3 have an external review (see `REVIEW_RESPONSE_PLAN.md`). Of the merged material, only PHASE5a §§1–3 records a review: its §4 is a fresh-session adversarial review. **No external review exists** for `PHASE1b`, `PHASE2c`, `PHASE2d`, `PACKAGE_DESIGN.md`, `NEW_AVENUES.md`, or the M-SIGN package code.

## 3. Items that conflict with corrections made in response to the review

| # | Where | Issue | Status |
|---|---|---|---|
| M-1 | `PHASE5a` Theorem R3 | Assumes "$A$ has nonempty relative interior". The reviewer showed that this is too weak for compactness (2a Remark 5.2): the condition is $\operatorname{affdim}(A)=R-1$, which A5 implies. Restate the hypothesis, or note that A5 already gives it. | **done** (WS5, see §6) |
| M-2 | `PHASE5a` Theorem R10(a), (c), (d), steps 13 and 15 | They say "P3.5(b) holds unchanged" and use "CLR's one-sided validity". After the review, P3.5(b) is split. (b1) is a theorem only with the **full** dual-vertex set. The package's simulated candidate collection is **Conjecture P3.5(b3)**. The projection interval (b2)(ii) and the conservative interval (a) are proved. So R10's statements for version (b) inherit conjecture status. R10 for version (a) stands. | **done** (WS5, see §6) |
| M-3 | `PHASE5a` §5, "CLR no longer falls back…" | Candidates from rows relaxed by the smallest uniform slack are dual-feasible for the original problem (the dual feasible set does not depend on $s$), so they are valid upper bounds by weak duality (P3.5(b2) step 1). That supports the device. Coverage of the CLR critical value computed over the collected set is still Conjecture (b3). | **done** (wording, §6) |
| M-4 | `PHASE1b` l.166 | "That is exactly P3 Thm 5(b)'s finite-candidate CLR setting, so it is valid uniformly over near-ties." Two problems: uniformity near point identification is an **open** item in P3 §10, and the theorem number is outdated (now P3.5(b1), with CLR Thm 4). | **fix needed** |
| M-5 | `PHASE1b` l.255–256, l.396; `PHASE2d` l.327 | "Clustered by county". **Case 1:** clustering the *same* county across two stacked periods is within-unit dependence, which P3-S1a allows and which is consistent with P3 §0.4. **Case 2:** clustering *townships within counties* ("52 county clusters", groups "as clusters for inference") is the cluster-of-units case. The response to the reviewer rejects that case and replaces it with group effects in the mean (M1/M2) and P3.6(c). Decide which case is meant and align. | **open** (needs your decision) |
| M-6 | `PACKAGE_DESIGN.md` l.87, 202, 404; `PHASE2d` l.329; `PHASE1b` l.256 | References to "P3 Thm 5" and "CLR stack" as proved. Renumber to P3.5(a)/(b1)/(b2)/(b3), and mark CLR-as-implemented as algorithmic. The G9 row ("already handled if all rows… enter one CLR stack") should say: proved for (a) and (b2), conjecture for the implemented CLR. | **fix needed** |
| M-7 | `PACKAGE_DESIGN.md` l.202, "M-LCE … proved (2a, 2b, P3)" | Qualify, per P3 §0.5: formal for M0 in any regime and M1/M2 with few large districts; CLR as implemented is a conjecture; R2 coverage is a simulation finding. | **fix needed** |
| M-8 | Package code: `cluster` argument | M-SIGN's API still exposes `cluster`. It is removed in WS5 (decision D3). Check that no M-SIGN path relied on it. | **done** (WS5, see §6) |

## 4. Items in the merged material that were never externally reviewed

| # | Item | What a reviewer should check |
|---|---|---|
| M-9 | `PHASE5a` Theorem R3 (sharpness of LCE ∩ restrictions) | the full proof; adversarial search S5 found 3 gaps where the certificate failed, so check the stated conditions exclude them |
| M-10 | `PHASE5a` Proposition R2 (when a restriction tightens) | the face-test argument; the T1b addendum |
| M-11 | `PHASE5a` Theorem R10(b) (parallel disjunction, version (a) only) and (c) (union of branch CIs at full level) | the step "$j^*$ is non-random" (9) when the restriction set itself is estimated; R10(e)'s counterexample for the naive bracket (coverage ≈ 0.91, persisting in $n$) |
| M-12 | `PHASE2c_monotone.md` (monotone EI, sign restrictions on $\delta$) | all propositions; 13 numeric checks pass |
| M-13 | `PHASE2d_did.md` and `PHASE1b_monotone_did.md` (DiD: FD, D0, D4, R1–R4) | all propositions; that "D4" here is a DiD result, not 2b Prop. D4 (name clash, see M-14) |
| M-14 | Naming | "D4" means 2b Prop. D4 (district dummies) in the review response, and the monotone-EI-on-changes rung in `PHASE1b`/`PHASE2d`/`PACKAGE_DESIGN`. Rename one, e.g. DiD rungs → "Δ4". |
| M-15 | `NEW_AVENUES.md` | it is a research agenda; any result quoted from it (e.g. cluster SEs l.173, cluster-bootstrap l.720) should be re-derived before use, and is subject to M-5 |
| M-16 | Package: `rxc_sign.R`, `rxc_modules.R`, M-SIGN changes in `rxc_core.R`/`rxc_infer.R` | code review against PHASE5a §5; the $2^m$ growth of branch systems; that `.rxc_ci_clr` refusing disjunctive rows is enforced |
| M-17 | M-SIGN coverage MC (gated) | was run once in the build session (conservative slab 1.000, CLR branch union 0.990); rerun with MC standard errors as in R-37 |

## 5. `verify_R3.py`

Rerun without `-I` (about 12 minutes): output identical to the committed `verify_R3_output.txt`, with T1, T2, T3 and T1b PASS and T4 reporting its counts. **checked**

## 6. Changes made in this branch to merged material

| Item | Change | Where |
|---|---|---|
| M-1 | Theorem R3: hypothesis is now $\operatorname{affdim}(A)=R-1$, with a note that A5 implies it | `PHASE5a_sign_module.md` §1 |
| M-2 | Theorem R10: status note. Version (b) statements are theorems only for P3.5(b1); with simulated candidates they inherit Conjecture P3.5(b3). Version (a) and projection stand. **The proofs of R10 were not re-reviewed.** | `PHASE5a_sign_module.md` §3.2 |
| M-3 | The relaxed-rows device is justified by weak duality; coverage over the collected set is the conjecture | `PHASE5a_sign_module.md` §5 table |
| M-8 | `cluster` removed from `eiBounds()`, `.rxc_setup()` and `.rxc_vcov()`. No M-SIGN code (`rxc_sign.R`, `rxc_modules.R`) referred to it. The M-SIGN tests pass unchanged (`test-msign.R`, `test-msign-dial.R`). | `R/rxc_api.R`, `R/rxc_infer.R`, `R/rxc_core.R` |
| (code) | `.rxc_ci_clr()`, which M-SIGN's `rxc_sign.R` also calls, gained an optional `proj_alpha` argument. With the default `NULL` its return value is unchanged, so M-SIGN's two calls are unaffected. | `R/rxc_infer.R` |


## 7. Consultation with the M-SIGN session (18 Nov 2026) and resolutions

At the user's request, the questions were put to the session that built M-SIGN (`session_01S2ZdrKuNjkb9Qi2eioNYtH`). It made no changes. Its answers were checked against the files, and the user approved the resulting list, decided M-5 as (b), and moved future M-SIGN work to this branch.

| # | Answer (summary) | Resolution | Status |
|---|---|---|---|
| M-1 | The affdim fix also had to go in §1.4 (degenerate hull) and in Prop. R2's setting, whose step 3 uses compactness | both fixed | **done** |
| M-2 | Agrees with the status note, and itemised the steps: (a) step 4; steps 13 and 15; (d); union steps 9–12 version-free; slab (b) proved. The CLR uses in `rxc_sign.R` ("already implied", near-tie margin) are also (b3) | itemised note added to R10 | **done** |
| M-3 | Wording correct | — | **done** |
| M-4 | Drop "uniformly". The sentence is about smooth nonlinear candidates without LCE, so cite the finite-candidate case of CLR Thm 4 | fixed in `PHASE1b` | **done** |
| M-5 | `PHASE1b` l.257 = same county in two stacked periods (within-unit). l.398 and `PHASE2d` l.329 = townships within counties. No reason specific to DiD to cluster | **user decision (b)**: l.257 reworded as "one unit observed twice"; l.398 and 2d l.329 now use county effects in the mean plus the pooling test, with no clustered variance | **done** |
| M-6, M-7 | No objection; extend the "algorithmic" label to wherever CLR is called the CR default | `PACKAGE_DESIGN` G9, M-LCE, R10 and item 5; `PHASE5a` §1.6 and §3.5; `PHASE2c` l.235; vignette; `PLAIN_LANGUAGE_GUIDE` §6 and §9 | **done** |
| M-9 | The certificate is sufficient, not necessary; R3(c) allows gaps | sentence added to §1.6 | **done** |
| M-11 | All branch rows are exact linear-in-β rows with design-only coefficients, so $j^*$ is fixed | no change needed | **done** |
| M-14 | Rename the DiD family D0–D8 → Δ0–Δ8; keep 2b's names | renamed in `PHASE2d`, `PHASE1b` (except l.124, which is P1/2b's D0), `PACKAGE_DESIGN`, `PLAIN_PHASE2d`, `ASSUMPTION_MENU` l.72, `PLAIN_LANGUAGE_GUIDE` l.198–199, `verify_phase2d.py` (output regenerated: 16 PASS) | **done** |
| M-15 | — | `NEW_AVENUES.md` cluster SE (l.173) and cluster bootstrap (l.720) marked for re-derivation | **done** (marking only) |
| scripts | Not intended; use a path insert | path insert added to `verify_R2.py`, `verify_R3.py` and `make_msign_reference.py`. All three now run under `python3 -I`: `verify_R2.py` and `verify_R3.py` reproduce their committed outputs exactly, and `make_msign_reference.py` regenerates the three test CSVs byte for byte | **done** |
| future work | Nothing unpushed. Planned: S11 real-data scripts; PHASE5a §5 proofs on slack monotonicity and the slack slab; slack docs. None touches `rxc_infer.R`/`rxc_core.R` | **user decision**: future M-SIGN work moves to `claude/wizardly-volta-lbh4bm` | **done** |

**Still open (substantive review, not consistency):** M-9 to M-13 and M-16 to M-17 need an external reviewer. M-15 needs re-derivation. See `EXTERNAL_REVIEW_PLAN.md`.
