# PHASE2e code review: changes and DiD (build step 2)

*Independent code review of `eiChange()`, `eiJoint()`, `eiATT()`, the model layers `fd()` / `lce_periods()` and the two one-line edits in package **eiPartialID 0.5.0**. Branch `claude/did-code` at 43c3acc. The reviewer did not write this code. The theory (`PHASE2e_did_general.md`, cited as 2e) is proof-reviewed (`PHASE2e_proof_review.md`). This review checks the code against it and does not re-review the proofs. Nothing in the package or the theory documents was edited.*

*Scripts, in this folder:*
* *`phase2e_code_review_oracle.py`: an independent oracle in Python (run with `python3.12 -I`; scipy 1.11 / HiGHS). It builds every row from its definition (2c Definition 1 with realised rates β = n/a; a row on changes is the association at t minus the one at u; orderings, contrasts, CR unions with slack). It does not use the package's thresholds or Lemma E4's closed forms. It enumerates every integer table of every unit-election for the brute force, solves the continuous LP per union branch, and solves the causal LP of Thm E16 in (n, θ_t, S_t, b_t), with Lemma E14's identity, the rate limits, the dial, RM as a union of 2L polyhedra on consecutive pre-period differences, SD, and MTR/ND/PS with spillover. Output: `phase2e_code_review_oracle_output.txt`.*
* *`phase2e_code_review_run.R`: the R driver for the oracle (runs `eiChange` with `integer = "auto"` and `"lp"`, and `eiATT`).*
* *`phase2e_code_review_probes.R`: the probes for edge cases, alignment, the model layer, E25 and performance (output `phase2e_code_review_probes_output.txt`).*
* *`phase2e_code_review_oldnew.R`: the 0.4.0 vs 0.5.0 comparison of existing outputs.*

*How to rerun them is in §6.*

---

## 0. Verdict in brief

The core is sound. Every set computed with default arguments that the oracle could check came out right:
* **1,831 comparisons, 0 failures** (`phase2e_code_review_oracle_output.txt`), on 261 instances:
  * 401 integer endpoints by brute force (2×2 with T = 2 and T = 3; 3×2 and 2×3 with T = 2). They cover margins, W1 δ_B/δ_W signs on changes of either sign, CR with and without slack, CR W2 and pair, pairwise δ_B, W2 δ_W, unit orderings vs the rest and vs a named row (with slack), the E6 contrast (both beliefs, rate alone), known cells, within-election totals, and rows on two period pairs (consecutive, and sharing a base).
  * 485 continuous endpoints, including 3×3 with T = 3 and 2 × 2 CR branches.
  * 94 causal sets (Thm E16) for no trend, the dial with M ∈ {0, .05, .3}, RM with M̄ ∈ {.3, 1}, SD, and MTR/ND/PS (all three; without PS; without ND). They use Λ ∈ {0, .03}, with and without the rate limits, a comparison row or the rest (R = 3), T = 3 and 4, CR on the pre or post change, pinned pre-period totals, and coupling rows across elections.
  * The E15 set and its closed form; E17(a), and the E17(b) criterion against the LP's ends (27/27, including coupling rows); E19's sup A, the implied dial and the inclusion E19(a); and E19(b) vacuity when there is one post-period (15/15).
* **§9.2, N21 (with and without rate limits, and the narrow variant), and review issues 1, 2 and 12 reproduce exactly** from the oracle, independently of the package (§1.2).
* The joint event-study set: `eiJoint()` agrees with the interval criterion of E3(c) on 300 random vectors.
* Non-integer margins: sets are scale-invariant (15/15).
* Units are aligned by label across elections, and shuffled rows and long format give identical results. N_i = 0 units are dropped without changing any set (but see issue 6).
* **Existing outputs are unchanged.** All 0.4.0 test files are unchanged and pass. `eiBounds(target = "model")` (CLR; conservative with an M-SIGN row), `eiCounts()` with CR and known cells, and the eiCounts cross-layer block print byte-identically under 0.4.0 and 0.5.0. The two one-line edits are inert for existing callers:
  * `P$classify` is set only by `eiChange`;
  * `st$change` is set only by `fd()`, and `isTRUE(NULL)` is `FALSE`.
* R CMD check reproduces the developer's result: 1 NOTE (eco, eiPack, ei.Datasets). The one WARNING is this container's missing en_US locale, not the package. Tests: 2,328 pass, 0 fail, 4 skip (CRAN, eco, ei.Datasets ×2).

Problems found: nine bugs or spec deviations that should be fixed before merging, and six minor items. None changes a default-argument result that the oracle checked, except issue 4 (a printed flag). The serious ones are:
* **issue 1**: an unbounded or failed causal LP is printed as an EMPTY set that "refutes" the assumptions;
* **issue 3**: `post =` silently drops the later post-periods' constraints, so the set is valid but not sharp for the stated RM/SD assumption;
* **issue 5**: for one-vs-rest targets with R ≥ 3, the allowance is computed for one representation of the target and added to an interval for another, which Lemma E12 does not cover;
* **issue 6**: numeric unit selections are shifted by dropped N_i = 0 units, so the wrong unit is selected silently.

| Component | Verdict |
|---|---|
| Input, (FP) and gap-register refusals (G1–G3, G5, G6, G7, G9, common base) | **met**, except issue 6 (numeric unit selections after N_i = 0 units) and issue 7 (election order in long data) |
| Stacking (Lemma E1), targets `did`, `gap_at`, `rate_at`, `bin_did` | **met** |
| Rows on changes (E4, E5(e)–(f), E6, orderings, multi-election `linear`) | **met** (0 mismatches against the definitions) |
| Closed forms E5(a)–(d), E6(b), E26 lattice | **met**; issue 11 (in integer mode the "closed form = LP" check compares the closed form with itself) |
| Integer classes (E26, E27, 6a classes on unit-elections, MILP) | **met**: the E27 router checks R = 2, one column, degree ≤ 2 with no cycle, laminar sets per pair, and boxes only (a total or a generic row routes to the MILP) |
| Joint set E3(c), `eiJoint()` | **met** |
| Causal LP (Lemma E14, Thm E16, E15, E17, E18, SD, E21, E22) | **met numerically**; issues 1, 2, 3, 8 |
| Reports E17(a)–(c), E19(a), (c2), E20, E21(c) | **met**; E19(b) issue 4; E20's display condition, issue 3 |
| Def. 5.1 labelling | **met** (ambiguity a needs the theory owner) |
| Model layer (E7, E8 plug-in, E9 closed form, E11, E12, E13/R27 wording, statuses) | **met**, except issue 5 (R ≥ 3 rest representation), issue 9 (sign rows without `sign`), issue 10 (wording) |
| E25 (H1) stacking, (H2) OUTER, (H3) boxes and α accounting | **met**; issues 3 and 5 carry over; issue 13 (wording) |
| Unions without multiplicity correction; no plug-in branch selected | **met** (count, model and E25 unions are reported whole; plug-ins are labelled and never intersected) |
| Existing outputs after the two one-line edits | **unchanged** |
| Tests | good coverage of the reference numbers and closed forms; gaps listed in issue 15 |

---

## 1. What was run

| Check | What | Result |
|---|---|---|
| R CMD check | `R CMD build` + `R CMD check --no-manual`, R 4.3.3, Rglpk/gmp/slam from Ubuntu `r-cran-*`, `_R_CHECK_FORCE_SUGGESTS_=false` | 1 NOTE (eco, eiPack, ei.Datasets unavailable); 1 WARNING from the container's locale (`en_US.UTF-8` cannot be set; not the package); tests **2,328 / 0 / 4 skip**; vignette builds |
| Oracle | 261 instances (seed 20261126), sets A–F plus reference sets W, N21, I1, I2 (§1.1) | **1,831 pass, 0 fail** |
| Old vs new | `phase2e_code_review_oldnew.R` under 0.4.0 (installed from bf011ea) and 0.5.0 | 4/4 printouts identical |
| Probes | `phase2e_code_review_probes.R` | issues 1–10, 13; checks of H1–H3, alignment, scale invariance, performance |

### 1.1 Oracle scenarios

| Set | Design | Rows | Compared |
|---|---|---|---|
| A (80) | 2×2, p = 3, N_i ≤ 6, T = 2 | none; δ_W±, δ_B± (W1); CR (slack 0, .05); CR W2; W2 δ_W; pair δ_B; order vs rest (±, slack .1); contrast (both; rate only); CR + known cell; δ_W + an interval total | DiD, a rate, a pairwise gap, bin DiD: integer (brute force), continuous, margins only |
| B (24) | 2×2, p = 2, T = 3, base 1 | CR, δ_W sign, or orderings on (0,1) and (1,2), or on (0,1) and (0,2) | event-study DiDs and a pairwise DiD (0 → 2) |
| C (40) | 3×2 and 2×3, p = 2, T = 2 | W1, pair, W2 signs; CR W1/W2/pair; orderings vs a row and vs the rest | one-vs-rest and pairwise DiD |
| D (28) | 3×3, p = 5, N 20–60, T = 3 | as C, plus CR W2 with slack, pairwise orderings with slack, CR in two pairs | continuous only |
| E (60) | 2×2 and 3×2, p = 4, T = 3/4 | none, CR on the pre or post change | causal sets for 10 trend forms, Λ ∈ {0, .03}, rate limits on/off, comparison row or rest |
| F (12) | 2×2, p = 3, T = 3 | a pre-period total pinned ±2, plus a coupling row 3 y₁ − 2 y₂ ∈ [y − 3, y + 3] | dial, RM, SD causal sets; E17(b) with coupling rows |
| W, N21, I1, I2 | §9.2 data; N21; review issues 1 and 2 | as in 2e | every row of the §9.2 ladder; N21; issues 1, 2, 12 |

### 1.2 Reference values (oracle, independent of the package)

| Item | Oracle (exact) | 2e (3 d.p.) |
|---|---|---|
| §9.2 DiD, margins only | [−0.4075, 0.7925] | [−0.407, 0.792] |
| §9.2 effect, no assumption | [−0.5675, 0.7325] | [−0.568, 0.732] |
| CR on the post change | [0.14325, 0.270283] | [0.143, 0.270] |
| + PT; dial .05; + Λ = .02 | [0.14325, 0.270283]; [0.09325, 0.320283]; [0.07325, 0.340283] | [0.143, 0.270]; [0.093, 0.320]; [0.073, 0.340] |
| RM, margins only | [−0.5675, 0.7325] | no information |
| RM, CR on the pre change only (issue 12); on both | [−0.415519, 0.7325]; [0.135231, 0.278302] | [−0.416, 0.732]; [0.135, 0.278] |
| pre-period DiD with CR | [−0.008019, −0.00425] (i.e. D_pre ∈ [0.00425, 0.008019]) | [0.004, 0.008] |
| MTR+ND+PS; with CR | [0, 0.4825]; [0, 0.311392] | [0, 0.482]; [0, 0.311] |
| N21 DiDs; RM without / with rate limits | [−1.229167, 1.270833], [−1.083333, 1.416667]; [−2.354167, 2.645833] / [−0.75, 1] | as 2e |
| N21 narrow: pre DiD; RM without / with rate limits | [−0.0625, 0.104167]; [−0.604167, 0.895833] / **the same** | [−0.06, 0.10]; [−0.60, 0.90] / [−0.604, 0.896] |
| issue 1 (coupling row) | max D = −0.4; PT set [−14/15, −0.8] | as 2e |
| issue 2 | margins-only D̄ = −0.15; contrast set empty | as 2e |

The exact §9.2 ends are halves of 10⁻³ (−0.4075, 0.7925, −0.5675, 0.7325, 0.4825). 2e rounds them in both directions, which is floating-point representation, not an error (ambiguity b). In N21's narrow variant the two sets are identical. 2e's sentence "[−0.60, 0.90] without rate limits, [−0.604, 0.896] with them" reads as if the rate limits changed the set, and they do not. This is a wording point for the theory owner, not a code issue.

---

## 2. Checks by topic

**Descriptive rows.** Every row type agrees with its definition on integer and continuous tables:
* W1 one-vs-rest rows reduce to half-lines in Z (E4(b)). The package computes them as differences of the one-election thresholds, with the slack applied once, in scale units of election t. The oracle builds the same rows from 2c Def. 1 with slack × scale on the change, and they agree, including CR with slack.
* Pair and W2 rows on changes have no constant term (they are linear in n). This is correct, and `gen_change` relies on it.
* `order_change` rows: Z ≥ (ΔB/(N−a) − slack)/κ_i vs the rest, and a_q Δn_r − a Δn_q ≥ −slack·a·a_q vs a row. Both are correct, and units with a_ir ∈ {0, N_i} are skipped.

**Integer classes.** The E26 lattice is used only when every row is a W1 sign/CR on the same (r, c, u, t, S) tuple, which is E26's hypothesis. `.ed_classify` implements E27's hypotheses. A within-election total, any generic row, R > 2, several columns, degree 3 or a cycle all route to the MILP. That is E27's remark ("LP outer, MILP exact"), so ambiguity f is resolved by the spec. Rows of kind "bound" (known cells, boxes, one-vs-rest unit orderings within an election, the odds dial) are single-cell boxes, which E27 allows. In set B (two pairs) and in the performance run (order_change in all 500 units), the TU path was taken (class "E27") and agreed with the brute force (set B).

**Causal layer.** The code's LP matches Thm E16 term by term:
* Lemma E14: b = θ − θ^PT + S.
* Rate limits θ ∈ [0, 1] and B_s − S ∈ [0, 1], with |S| ≤ Λ, and S = 0 when no comparison group is stated.
* RM: one branch per (pre pair, sign), with consecutive differences δ_t − δ_{t−1} for every post-period in the system. A(n) is taken over the L consecutive pre pairs, ending with (−1, 0).
* SD: second differences over all elections, with pre-period δ_s = D_s.
* MTR, ND and PS as in E21.

The oracle's independent LP agrees on all 94 sets.

**Def. 5.1.** Every refusal fires with the item named: population, treated group, comparison group when the trend uses one (dial, RM, SD, PS), trend, spillover, and no anticipation (missing or `FALSE`). A population that differs from the fit's is refused. MTR+ND without a comparison group is accepted; with a spillover still required (ambiguity a). A negative Λ, or a comparison equal to the treated group, is refused.

**Model layer.**
* FD blocks are `.rxc_setup` fits on ΔT with `st$change = TRUE`. The last column's identified constant is then 0, so Σ_c ΔH_c = 0. The joint sandwich's diagonal blocks equal each block's own V, the β ordering matches, and the stacked V is positive definite (H1 (ii-a)). Overlapping or non-consecutive pairs are refused.
* ρ_i = 2(w_ir + w_is) for a DiD with a named row, or with R = 2 (dev test, and re-derived here).
* The Popoviciu allowance is z_{1−α_E/2}·½(Σρ²)^{1/2}, and the Hoeffding radius is (½ log(2/α_E) Σρ²)^{1/2}, as in E12.
* The R27 wording is carried over. The identity-only line prints for known cells and value totals. With a `cr_change` (or any sign, order, contrast, or interval total) the belief line "coverage at least 1 − α − P(belief rows false)" prints.
* Statuses: pooled FD and pooled LCE get "asymptotic guarantee (P3.5(b2)(ii) …)"; P3.5(a) when `ci = "conservative"`; a grouped rung with rows gets "R16 open", and without rows "regime R1". These are correct.
* The E9(a) closed form equals the LP (dev test). E8 is a plug-in at the empirical law. With unique compositions it equals FD4's outer set, as E8(c) predicts; with 4 support points it is strictly inside ([−0.917, 0.889] vs [−0.938, 0.954]).

**E25.**
* (H1) Consecutive FD blocks or per-election LCE blocks, one joint V, one projection event (the same draws and critical value) for every branch and post-period.
* (H2) FD gives `outer = TRUE` with no rate-limit rows, and prints "OUTER"; LCE adds θ^PT and B_s quantities with rate limits.
* (H3) One box e_q ∈ [−t_q, t_q] per actual quantity: post DiDs, consecutive pre DiDs (RM) or pre DiDs relative to the base (SD), and under LCE also θ^PT_t and B_s,t. α_proj = 0.9α and α_q = 0.1α/|q|, so α_proj + |q|·α_q = α in every probe.
* MTR and no-trend are refused. Clipping to [−1, 1] is valid (an ATT of a rate).

---

## 2b. Verdicts on the developer's deviations and ambiguities

| # | Deviation | Verdict |
|---|---|---|
| 1 | Causal LP for continuous tables only | **Consistent with 2e** (Remark 7.1: integer causal ends are O) and printed. |
| 2 | E8 sharp set: 2×2, one pair, plug-in at the empirical law, no inference | **Consistent** (E23(c): sharp-set inference is O), and labelled a plug-in. "Nearly empty with continuous designs" did not reproduce here: with unique compositions it equals FD4's outer plug-in (E8(c), degenerate case). Suggest printing that equality when the support is the units themselves. |
| 3 | FD R1/R2/R4/δ for R = 2 only | **Consistent** (2d FD2/FD3 are for R = 2; E7 applies them per column in 2×C, which the code allows). |
| 4 | FD validity rows \|Δm\| ≤ 1 at the hull vertices | **Acceptable.** They are valid (Δβ ∈ [−1, 1], and the hull vertices are data compositions), and they are linear rows of E23(b)'s kind, so they are covered. They are not the tightest valid rows. |
| 5 | Projection or conservative only; α 0.9/0.1; α_q = 0.1α/\|q\| | **Acceptable** (E23(b): CLR algorithmic; the union bound is E25(a) as stated). Verified. |
| 6 | E25 for the dial, RM and SD; MTR refused; clip to [−1, 1] | **Consistent** (E25 is stated for "each branch of E18(a), or for SD or the dial"). Clipping is valid. Issue 13 covers an interval wholly outside [−1, 1]. |
| 7 | E6 for R = 2 | **Consistent** with E6 ("with R ≥ 3 only named pairs"). |
| 8 | E19(b) for one post-period with Λ = 0 | **Consistent in intent, not as implemented**: the flag is computed whenever k = 1, even with later post-periods in the system (issue 4). |
| 9 | E3(c) automatic test on the first two DiDs | **Acceptable** (illustrative; `eiJoint()` is general and verified). |
| 10 | Interface; `rate_limits = FALSE` labelled OUTER | **Acceptable.** But `rate_limits = FALSE` with forms that leave θ unbounded triggers issue 1. |
| 11 | Model signs allowed under FD; grouped rungs flagged R16 | **Consistent** (E23(b) lists model signs on changes; E11 and the matrix give "O grouped (R16)"). Issue 9 is a parsing bug in these rows. |

| | Ambiguity | View |
|---|---|---|
| a | Def. 5.1 item 3 requires a spillover assumption even for MTR+ND without PS (no comparison group) | **Needs a theory decision.** As written, Def. 5.1 lists item 3 unconditionally, while item 1 is conditional ("the comparison group when the trend assumption uses one"). S_t is defined only relative to a comparison group, so without one the spillover statement is vacuous: the code requires Λ and then ignores it (S fixed at 0). The code follows the letter and is safe. Recommendation to the theory owner: amend item 3 to "when a comparison group is used", mirroring item 1 (E15's remark already says MTR and ND need none). If item 3 is kept, the print should say "spillover: not used (no comparison group)" rather than "|S_t| ≤ Λ". |
| b | §9.2 rounding | No code issue. The exact values are halves of 10⁻³ (§1.2), and the tests' 6e−4 tolerance is right. Suggest the theory owner print 4 decimals, or the exact values, in §9.2 and fix the N21 narrow-variant sentence (§1.2). |
| c | E17(c) breakdown only in the outer form | **Consistent**: E17(c) is stated "on the Minkowski form". A sharp breakdown (smallest M with the sharp lower end ≤ 0, by bisection over the LP) would be an optional addition, not a fix. |
| d | E19(c2) uses the width at base − 1 when L ≥ 2 | **Correct**: 2e defines w_{−1} as election −1's gap width, and the bound Ā(g₀) ≥ ½w_{−1} holds for every L. Nothing to change. |
| e | Level targets get no model interval under FD | **Correct** (E12(i)). But `bin_did`, a change target, gets the same message (issue 10). |
| f | Totals mixing elections with change rows go to the MILP | **Correct** (E27's remark: "change rows mixed with within-election totals of 6a's classes are outside the theorem: LP outer, MILP exact"). |

**Theory-owner decisions needed:** ambiguity a (Def. 5.1 item 3), and the wording in ambiguity b. Issues 3 and 4 follow from 2e as written (Defs 5.2–5.3 hold "for every post-period" or "every t", and E19 is stated for one post-period), so they are code fixes, not theory decisions. Optionally, the theory owner can confirm that a user who requests one post-period still means the assumption for all elections in the data.

---

## 3. Numbered issues

### Issue 1 (bug, medium): an unbounded or failed causal LP is printed as EMPTY, "refutes …"

`.att_lp` (R/att.R:114–116) calls `Rglpk_solve_LP` with canonical statuses and returns `NULL` for any status ≠ 0. Rglpk's canonical status is 1 for "infeasible", for "unbounded", and for any failure or time-out (checked: both an infeasible and an unbounded LP return 1). `.att_range` treats `NULL` as an empty branch. If every branch fails, the print says "EMPTY -- refutes (FP) + no anticipation + count rows + the trend and spillover assumptions (Thm E16(d))".

*Reproduction* (probes, section probe1):
* `eiATT(fit, …, trend = no_trend(), rate_limits = FALSE)` gives `lower = upper = NA`, printed EMPTY. The true set is (−∞, ∞).
* `trend = mtr_nd_ps(nd = FALSE, ps = FALSE), rate_limits = FALSE` gives the same. θ is unbounded below, so the ATT is unbounded above.

The same path would turn a numerical failure of any causal LP into a false refutation. 6a's `.ec_optimise` already handles this correctly (`canonicalize_status = FALSE`; 5 optimal; 3/4 infeasible; anything else falls back to lpSolve).

*Fix.* Mirror `.ec_optimise`: `control = list(canonicalize_status = FALSE, tm_limit = …)`. Treat status 5 as optimal, 3/4 as infeasible (empty branch), 6 (unbounded) as ±∞, and anything else as an error or a fallback, never "empty". Apply the same to `.att_nrange` and the E17/E19 feasibility LPs.

### Issue 2 (bug, low–medium): eiATT crashes on an infeasible count layer for the dial and RM

When the count rows contradict the margins, the fit's status is "infeasible …". `eiATT` then reports `NA` for no-trend and MTR, which is correct (E16(d): the conjunction is refuted). For the dial it stops with "missing value where TRUE/FALSE needed" in `.att_dial_report` (`if (D[1] - Lambda > 0)` with D = NA). For RM it stops with "trend_dial(M) must be one nonnegative number" in `.att_rm_report` (sup A = −Inf).

*Reproduction:* probe1 P6 (a known cell of 11 for a unit with 10 veterans).

*Fix.* In both report functions, return `NULL` when D or the pre-period ranges are NA. The main set is then EMPTY with the E16(d) message.

### Issue 3 (spec deviation, sharpness): `post =` drops the constraints of later post-periods

`eiATT` builds the system over post-periods `seq(b0 + 1, max(tgt))` (R/att.R:356; also `cs$posts` in `.att_model`). Requesting `post = "p1"` in data with a later post-period p2 therefore drops:
* p2's RM constraint |b₂ − b₁| ≤ M̄A(n), or SD's second differences involving b₂;
* and p2's rate limits.

The printed set for ATT_{p1} is then the set under a weaker assumption than the one stated (Def. 5.2: "for every post-period t ≥ 1"; Def. 5.3: "for every t with t ± 1 ∈ 𝒯"; Thm E16 in vector form). It is still valid, but not sharp, and it differs from the default output for the same target.

*Reproduction:* probe2. In 31 of 60 random designs (T = 4, base 2) the set for ATT_{p1} changes with `post = "p1"`. Design 4, SD(.02): [−0.060, 0.800] with all posts, [−0.530, 0.800] with `post = "p1"`. The dial is unaffected (its posts are not coupled).

*Fix.* Always build the system over every post-period in the data and report only the requested ones. Do the same in `.att_model`. Then show E20's closed form only when the system has one post-period (T = b0 + 1). The current condition `length(tgt) == 1` holds only because of this truncation. If truncation is ever wanted, it should be an explicit argument with the label "assumption imposed through election …".

### Issue 4 (spec deviation): the E19(b) vacuity flag is computed outside E19's hypothesis (one post-period)

`.att_rm_report` evaluates the E19(b) criterion whenever `k == 1 && Lambda == 0` (R/att.R:435), even when the system contains later post-periods. With later post-periods the RM constraints chain |b₂ − b₁| ≤ M̄A(n), and p2's rate limits can exclude θ₁ = 0 or θ₁ = 1. The criterion, which looks only at t = 1, then mis-reports.

*Reproduction:* probe9 (T = 4, base 2, default posts). In 22 of 40 designs the flag disagrees with the computed set. For example, k = 1: the set is [−0.684, 0.297], no-assumption [−1, 0.297], yet "lower end uninformative (= B^L − 1)".

*Fix.* Compute E19(b) only when the system has one post-period (after issue 3's fix: T = b0 + 1); otherwise print nothing. The alternative is to add the other post-periods' constraints to the 2L feasibility LPs, which is outside E19 as stated.

### Issue 5 (bug, validity gap): the allowance and the model-implied target use different representations (one-vs-rest, R ≥ 3)

For a one-vs-rest DiD (or gap) with R ≥ 3, `tg$reps` has two representations. Rep 1 is (1/A_r)Y_rc − (1/(N−A_r))Σ_{q≠r}Y_qc. Rep 2 is κ_r·Y_rc − Y_c/(N−A_r); its constant is built from the observed outcomes.
* `.ed_allowance` (R/change_model.R:331–337) chooses the representation with the smaller S_φ.
* The model-implied target always comes from `.ed_model_target` → `.ed_gapW`, which is rep 1.

Lemma E12 (6a Lemma 4.1) bounds φ(n°) − φ̃ for the representation whose φ̃ is targeted. 6a's `.ec_model_interval` does this correctly: it builds `tgm` from the chosen representation, with `c0 = rep$const`. Here, when rep 2 is chosen, the interval covers φ̃₁ while the allowance bounds φ − φ̃₂. The difference

  φ̃₂ − φ̃₁ = (1/(N − A_r)) Σ_i N_i (ΔT_ic − E[ΔT_ic | X_i])

is a residual sum that no allowance accounts for. The same mismatch is in `.att_model` (E25), whose `gap_reps` also carry `terms2`.

*Reproduction:* probe3 (R = 3, 300 units, FD, versus = "rest"). S_φ(rep 1) = 0.1187 and S_φ(rep 2) = 0.0976; the package uses 0.0976, but the model form's constant is 0 (rep 1). At the true rates φ̃₁ = 0.13997 and φ̃₂ = 0.14151.

The practical effect was small here: the difference is 0.0015 against an allowance of 0.274, and Monte Carlo coverage of the actual DiD was 150/150, because FD intervals are wide when ΔK is free. But the coverage claim printed with the interval is not covered by E12 as applied. R = 2 is unaffected, because both representations coincide.

*Fix.* Build the model-implied form from the chosen representation, as 6a does: terms → weights through `.rxc_mlin`, plus `c0 = rep$const`. Under FD the rep-2 form is a change target, since its constant uses the observed ΔY_c. Alternatively, restrict the choice to the representation the model form uses.

### Issue 6 (bug, medium): numeric unit selections are shifted by dropped N_i = 0 units

`.ed_problem` builds the per-election problems `P$Pt` from `P$A0`, which already excludes N_i = 0 units, so their `orig` is 1..p₀. Every numeric unit selection is then matched against positions among the kept units, not the user's data rows:
* `where = <row numbers>` in targets;
* `unit = / units =` numbers in `known`, `total`, `sign_change`, `order_change` and other rows;
* numeric `unit` in multi-election `linear` cells (`match(as.integer(cl$unit), P$Pt[[1]]$orig)`).

With an N = 0 unit in data row 1, `where = 2` (unit P1) selects P2. `known` with `unit = 2` is imposed on P2. Both are silent. A logical `where` vector of full length fails ("must have one value per row of the data").

*Reproduction:* probe1 P3, P3d. Labels and expressions work (P3 "by label", P3c).

*Fix.* Build `P$Pt[[t]]` with the full aligned data and let `.ec_problem` drop the N = 0 units itself (so `orig` and `n_orig` refer to data rows). Alternatively, map numeric selections through `which(keep)`. Document that numeric positions refer to the first election's row order, since later elections are realigned to it by label.

### Issue 7 (bug, medium): long-format data with a character election variable is ordered by first appearance

In `.ed_input`, a character election variable without `elections` gives `lev <- unique(as.character(ev))`, which is the order of first appearance in the data. Time order matters for every pre-period quantity, RM, SD, the FD pairs and the event-study sign. Shuffled long data silently reorders the elections.

*Reproduction:* probe1 P4c, where shuffled rows give "post, early, base". Numeric and factor election variables are handled correctly.

*Fix.* Require `elections` (or a factor with ordered levels) when the election variable is character, and stop otherwise. At minimum, warn.

### Issue 8 (bug, low): `smoothness()` with T = 2 crashes; the same loop is in E25

`.att_system` and `.att_model` loop `for (j in 2:(length(els) - 1))`. With T = 2 that is `2:1`, which iterates j = 2, 1, and `dl[[3]]` fails with "subscript out of bounds" (probe1 P2).

*Fix.* Use `for (j in seq_len(max(0, T - 2)) + 1)`. With no triple, SD imposes nothing, so the set is E15's: print that, or refuse SD with fewer than three elections.

### Issue 9 (bug, low): model sign rows without an explicit `sign` fail (partial matching)

In `.ed_fd_rows` and `.ed_lce_rows`, `s$sign` partially matches `s$sign_type` when `sign` is omitted. `.ed_sign_arg("dW")` then stops with "sign must be +1 or -1". The documented default (sign ≥ 0) is therefore unreachable whenever `sign_type` is given (probe7).

*Fix.* Use exact extraction, `s[["sign"]]`. Audit the other `$` reads on user lists for the same hazard; for example `s$unit` after `s$units` is safe only by order.

### Issue 10 (wording): `bin_did` under FD gets the level-target message

`.ed_model_target` returns "FD leaves levels free: a level target has no model-implied value (2e Lemma E12(i))" for `bin_did`, which is a change target (probe1 P10). The true reason is that no model-implied version is implemented. E10 covers compositions; bin aggregates are evaluated at row-specific weighted compositions, which E10's own remark excludes.

*Fix.* Give `bin_did` its own note.

### Issue 11 (test gap / wording): in integer mode the "closed form = LP" cross-check is vacuous

On the E26 path, `.ed_e26_range` sets `out$lp` to the closed form's continuous hull. `.ed_fit` then compares `closed_*` with `lp_*`, which is the closed form with itself, and the `lp_lower/lp_upper` columns hold no LP. The oracle confirms the values are right, so no number is wrong.

*Fix.* Run the continuous LP on that path too (it is cheap), or rename the columns and skip the check there.

### Issue 12 (wording): vignette overstates E19

"Relative magnitudes … On the margins alone they collapse to no information." E19 shows that they reduce to a trend dial whose M is at least M̄w_{−1}/2, and *usually* to no information. The sentence is true for the example.

*Fix.* "Here, on the margins alone, they give no information" (or "usually").

### Issue 13 (wording, minor): E25 print details

(a) `.att_model` clips each branch interval to [−1, 1], but leaves an interval that lies wholly outside [−1, 1] unclipped. That case is a rejection at the stated level and should print as EMPTY. (b) For a grouped LCE rung without model rows the E25 status says "R16 open". The eiChange block correctly distinguishes "regime R1" (no rows) from "R16" (rows); E25 should do the same.

### Issue 14 (wording, theory documents; not code): §9.2 rounding and N21 narrow variant

See §1.2 and ambiguity b. For the theory owner; the code is right.

### Issue 15 (test gap)

The tests cover the reference numbers, closed forms, labelling and unions well. Add tests for:
* unbounded forms with `rate_limits = FALSE` (issue 1);
* eiATT on an infeasible fit (issue 2);
* `post =` invariance (issue 3);
* E19(b) with two post-periods (issue 4);
* R ≥ 3 one-vs-rest model intervals (issue 5);
* numeric unit selections with N = 0 units (issue 6);
* long data with character elections (issue 7);
* SD with T = 2 (issue 8);
* model sign rows without `sign` (issue 9);
* a brute-force check of R = 3 rows on changes (pair, W2, CR W2), as in oracle set C.

---

## 4. Edge cases and performance

| Case | Result |
|---|---|
| Unit order differs across elections | aligned by label; identical sets (probe1 P4) |
| Long format, numeric or factor election variable | correct; character → issue 7 |
| N_i = 0 units | dropped, same sets (P5); numeric selections → issue 6 |
| Group counts change between elections | refused as G1/G2/G3/G5/G9 (dev test, verified) |
| Outcome totals ≠ group totals (voters-only base) | refused ("common base") |
| Different unit sets | refused as G6, with Lemma G6's advice |
| Non-integer margins | integer check off, sets scale-invariant (probe8: 15/15) |
| `margins = "estimated"` | G7 warning, integer check off; the model layer warns that R26 is deferred |
| Infeasible rows | count layer: "infeasible: the rows contradict the margins"; eiATT: issue 2 |
| CR unions with slack | oracle sets A and D (slack .05, .02): exact |
| T ≥ 3 | oracle sets B, D, E (T = 3, 4); SD with T = 2 → issue 8 |
| 500 units, 3 elections | margins only, integer auto: 4.5 s; CR on both changes: 12 s (LP) / 16 s (auto, TU); order_change in every unit (E27): 12 s; eiATT RM on the CR fit: 49 s; SD: 18 s; FD + DiD-S: 2 s; E25 RM: 0.3 s. Acceptable. Most of RM's time is the report LPs (E19 dial, vacuity), which could be skipped on request. |
| Silent fallbacks | issue 1 (causal LP statuses). Elsewhere the count layer inherits 6a's careful status handling, and closed-form/LP disagreements warn. |

---

## 5. Suggested order of fixes

1. Issues 1, 6, 7 and 5: wrong or unsupported output without a warning.
2. Issues 3 and 4: sharpness, and a wrong printed flag.
3. Issues 2, 8 and 9: crashes and parse errors.
4. Issues 10–13 and the tests in issue 15.

Ambiguity a goes to the theory owner. None of the fixes needs a change to 2e.

---

## 6. How to rerun

```
# R 4.3 with Rglpk, gmp, slam (Ubuntu r-cran-*); scipy for python3.12 (python3-scipy)
R CMD INSTALL -l LIB package/eiPartialID
cd docs/extension/reviews
python3.12 -I phase2e_code_review_oracle.py gen /tmp/inst.json            # instances + oracle answers
Rscript phase2e_code_review_run.R /tmp/inst.json /tmp/res.json LIB          # the package on each instance
python3.12 -I phase2e_code_review_oracle.py check /tmp/inst.json /tmp/res.json
Rscript phase2e_code_review_probes.R LIB                                    # SKIP_PERF=1 to skip the 500-unit run
Rscript phase2e_code_review_oldnew.R LIB_0.4.0 old.rds; Rscript phase2e_code_review_oldnew.R LIB new.rds   # then identical() per printout
```
