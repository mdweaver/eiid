# What to build next in eiPartialID: a recommendation

*Planning note, 20 November 2026; revised 21 November after your decisions D1–D12. Written against `claude/wizardly-volta-lbh4bm` at 7077dfd (v0.3.1 on `master`, the round-2 corrections, the technical report and the methods paper). Nothing in the package, the theory documents, the report or the paper is changed. Checks are in this folder: `verify_next_steps.py` (labels N1–N24, 34/34 PASS, output in `verify_next_steps_output.txt`) and two population studies, `sharp_gain_study.py` and `sharp_gain_study_2xC.py` (outputs alongside).*

*Scope (your rule, 21 Nov): this repository develops **general theory and methods**. Wisconsin/Iowa and NZ appear only as **illustrations** of the data shapes the methods are built for:
* aggregates with a second-scale composition and repeated elections;
* aggregates with published tables;
* matched (estimated) group counts.

*Project-specific choices are not deliverables of this plan: data assembly, match calibration, covariate sets, township selection and Wisconsin-specific dial settings. They are collected in a short note at the end (Part 9).*

*How to read it. Part 1 is the one-page summary. Part 3 has one section per recommended step. Each section starts with three plain-language parts (what it is, a worked example, why it comes at this point), which you can read without the technical notes that follow them. Part 4 works through NEW_AVENUES.md avenue by avenue. Parts 5–8 cover what to leave for later, the external review, the decisions needed from you, and what was verified.*

---

## Part 1. One-page summary

**The recommendation in one sentence.** Build next a **count-table (“realised electorate”) layer** (`eiCounts()`, printed first by `eiBounds()`), with honest handling of **estimated group margins**. Then a **general DiD theory track** (descriptive and causal, proofs first). Then, on top: finite comparisons, the second-scale context dial B1(M) with κ annotation, covariate adjustment CAR(Z)-δ (sharing B1's model), exception frontiers, neighbour-bounded district deviations (an M1(δ) option), the M-SHARP diagnostic, the Bernstein robustness ladder, and contests as districts (C2). NZ table transfer (E4) comes last. Audits and ballot patterns form a general-user track. Each step gets a separate proof review and code review. Part 3A answers your layering question: everything combines, either in one LP within a layer or by intersecting the two layers' answers, apart from a few named incoherent cases.

**Why.** Three findings drive the order.

1. **Most applied questions are about the people who actually voted**, for example “what share of these veterans voted yes?”. Today the package answers that question only through the Duncan–Davis column and `ci_actual`, which come last in `print()`. Every optional restriction acts on the *model-implied* rate (`eiModules()`: zeros, order, sign, slope, linear, external, assoc and cr all restrict the mean). The paper (§7) lists sign rows on realised tables as not implemented. The realised-table results are already proved (2c Prop. M2; DD; PHASE2d Δ0/Δ4), so the gap is in the code, not the theory.
2. **M-SHARP's gain depends on one feature of the data: polarised precincts.** In an oracle study that knows the law of returns exactly (`sharp_gain_study*.py`), using the whole law instead of the regression narrowed the model-implied group rates by a median of **0–1%** (at most 5%) when precinct rates are spread unimodally around their means (2×2, 2×3 and 2×4; SD up to about .25). There, the review's subset rows added nothing beyond averaged cell limits. When precinct rates are **bimodal** (some precincts near 0 and others near 1 within a group), the gain is **8–37%** (median 17%; 2×2, small or large precincts). Small precincts alone add 1–8%. So M-SHARP is neither worthless nor generally useful. A cheap "kink mass" diagnostic on your NZ and Wisconsin data should decide whether to build the estimator. The review's subset identity is correct (verified, N1–N3), and it is the right tool if the estimator is built.
3. **The applications in view have aggregate margins only** (illustrated by the Wisconsin/Iowa suffrage data and the NZ tables). The count layer serves such data directly: the exact DiD rungs are realised-table results, and comparisons and exception frontiers need no extra data. The modules with the largest gains for *other* users bring in information (an audited precinct, a published total, a cast-vote record); in the count layer each is an exact equality or extra columns, so they are cheap to add later (verified, N6–N11).

**What changed after your decisions of 21 Nov.**
* D3 accepted, with E4 moved last (Step 11).
* D1: `print()` leads with the actual-rate block, DD first.
* D2: a separate `eiCounts()`, embedded and printed by `eiBounds()`.
* D4: DiD becomes a **general theory track** (descriptive and causal), with changing populations deferred and named in a gap register. Implementation follows the proof review.
* D5: covariates are potentially available, so **CAR(Z)-δ is Step 6**, sharing B1(M)'s model.
* D6: integer check "auto", and the matched (non-integer, estimated) veteran counts lead to a new **Step 1b**. M-MISCLASS moves into your track, with margin sets from the match probabilities and a Berkson correction for the LCE layer (R25, R26).
* D7: two separate reviews (proof, code) per step.
* D8: the M-15 flags go to the NEW_AVENUES author.

**What changed after working through NEW_AVENUES (Part 4).** Avenue B's B1(M) enters as Step 4: it serves your Wisconsin/Iowa data directly (the county share is the second scale), and it is the DiD's yardstick. κ annotation joins it. The odds-ratio dial joins Step 1, outcome nesting joins Step 2, and the Bernstein robustness ladder becomes Step 9. Avenue A, the single index and the covariance route stay research.

**What changed after the owner's data answers (20 Nov).** No poll books, linked records, CVRs or sub-tables exist for your projects; your applications use aggregate margins only (NZ electorate tables; the Wisconsin/Iowa DiD). The order below therefore ranks modules by their value for that work first, and moves the modules that need external, audit or ballot-level data into a **general-user track** (§3.8). They stay in the plan, because other users have such data, but they are labelled as serving other users.

**Recommended order**

| Step | Module (working name) | What the user gets | Typical data (illustration) | New proofs | Data | Effort |
|---|---|---|---|---|---|---|
| 0 | Finish in-flight work | tag v0.3.1; print fixes (§3.0); the planned batches D/E become the proof and code reviews of existing modules | all | none | none | S |
| 1 | **M-COUNT**: count-table layer, `eiCounts()`, printed first by `eiBounds()` | exact bounds for the actual rates, gaps and sub-area totals; realised signs, orderings, odds-ratio dial and CR (one vocabulary with `object =`); known cells and totals; integer check "auto"; the cross-layer combination rule and free LCE test | any aggregates | R19, R27 | margins | M |
| 1b | **Estimated margins** (M-MISCLASS, moved in): margin sets from match probabilities; Berkson correction in the LCE layer | honest ranges when group counts are expected counts from matching | matched group counts (e.g. veterans matched to a census) | R25, R26 | match probabilities | M |
| 2 | **DiD theory track**: general descriptive and causal theory with a gap register for changing populations (PHASE2e), proof review; then `eiChange()`/`eiATT()` and code review | a general, proved DiD layer; ATT only with explicit trend dials | repeated elections, same units (e.g. Wisconsin/Iowa referenda) | R11a–d; R25/R26 across periods | two elections' margins | L (theory), then M |
| 3 | **M-CONTRAST**: finite-composition comparisons, in levels and in changes (M-LOC as its "hull" option under LCE; the finite form of NA's no-offsetting ceiling) | signs on *prespecified* comparisons, with an upper cap on the gap | any aggregates | R22 | margins | S–M |
| 4 | **M-CTX2**: B1(M) with κ-annotated endpoints, per column, range form | a dial benchmarked by the identified county-level contextual effect; the DiD's yardstick | a composition at a second scale (e.g. county share) | R8 (B1) | second-scale composition | M |
| 5 | **M1(δ) with named neighbours**: bounded differences between neighbouring districts, defined by an adjacency or a distance matrix | district deviations tied to well-identified neighbours instead of one common centre | grouped data with a geography (districts, counties) | R30 (R16-type status) | adjacency or distance matrix | S–M |
| 6 | **CAR(Z)-δ**: covariate adjustment with a dial, sharing Step 4's model | covariate-adjusted bounds with δ and κ curves; covariate-benchmarked κ_k/δ_k with the bridge "remainder ≤ m × covariate k"; diagnostics R²(X~Z) and overlap; shown beside B1(M) | unit covariates with overlap | R8 (CAR), R28, overlap factor | covariates with overlap | M |
| 7 | **M-EXCEPT**: exception-mass frontiers | "how many exceptional precincts before my conclusion fails?" | any aggregates | R23 | margins | M |
| 8 | **M-SHARP diagnostic** (estimator only if polarised precincts are found) | when the full law could tighten an LCE bound | any aggregates (diagnostic) | R24, R9′ | margins | M |
| 9 | **M-BERN**: polynomial outer/inner ladder (with L1′) | how much of an LCE answer rests on straight lines | any aggregates | designed (NA §4) | margins | M |
| 10 | **C2: contests as districts**: several contests on one ballot share contextual structure up to level and scale | lopsided contests pin the hidden pattern for close ones; testable restrictions | several contests, same units | R29 | contest returns per unit | M |
| 11 | **M-DATA E4**: NZ leave-one-electorate-out transfer (last, by your decision) | validation of the model layer; a reference distribution for Steps 4–6's dials | published area tables (e.g. NZ) | R21 + allowance | published tables | M |
| G1–G2 | **General-user track**: M-AUDIT, M-CVR; E5 rows ship inside Step 1 | features for users with audits or CVRs | audits, CVRs | R20; none | external data | S, S |
| later | KL-ball reference dial; R16/R15 for M-SIGN (theory track); 2×C joint signs (R6), R×2 trees (R5). Research: Avenue A on pooled NZ (with the dispersion ceiling), single index (R14) | — | — | — | — | — |

Every step gets a **proof review** and a separate **code review** (D7; table at the start of Part 3).

**Not now** (Part 5): a series-based M-SHARP estimator (unless the diagnostic finds polarised precincts); M-LOC as a separate module; Avenue A and the dispersion ceiling (research); the covariance diagnostic (dropped by your decision); DiD code before its theory has passed proof review; switching the default interval (this would reverse a decision of yours); general R×C full-law estimation.

**Decisions** (Part 7). D1–D12 are settled or deferred:
* D9 and D10 are deferred to your project; in the plan they are user inputs.
* D11: T ≥ 2 and all three causal layers.
* D12: CAR(Z) before M-EXCEPT, as ordered.

Added by your decision: C2 (contests as districts) as Step 10, just before E4; bounded differences between named neighbours as an option of M1(δ), Step 5 (adjacency or distance matrix). Dropped by your decision: the covariance diagnostic (one line kept inside Avenue A). No decisions remain open.

---

## Part 2. Where the project stands, and what the external review got right

### 2.1 Status in brief

* **Package 0.3.1**:
  * the LCE ladder (DD, M0, M1, M1δ, M2, L4–L6);
  * CLR (algorithmic, the default), projection (proved for pooled M0 and for few large groups) and conservative intervals;
  * emptiness, pooling and M1 tests with a wild bootstrap;
  * sensitivity and breakdown;
  * M-SIGN (experimental; proved for pooled M0 only).
* **Proof register**:
  * open: R1 (C1), R4, R7, R9, R11 (beyond DiD-S), R14, R15, R16;
  * checked only: R5(c), R6;
  * partly done: R8.
* **Not reviewed**:
  * PHASE1b and PHASE2d (batch C);
  * PACKAGE_DESIGN and NEW_AVENUES (batch D);
  * the M-SIGN and WS5 code (batch E).

  So the build order in PACKAGE_DESIGN §6 has not had an outside check. None of batches C, D and E has been sent: C is gated on the DiD data decisions (now partly answered: no poll books exist); D is unblocked but unscheduled; E has no reviewer.
* **The repository's own next step** (`NEXT_BUILD.md` item 4) is M-SHARP as a series bound, followed by M-LOC (C1), 2×C joint signs (R6), R×2 trees (R5) and then the DiD wrappers. PACKAGE_DESIGN §6 justifies this order by proof readiness.

### 2.2 Which of the review's points (written at 27557df) later work already addresses

| Review point | Status after v0.3.1, round 2 and the paper | Still to do |
|---|---|---|
| Projection vs CLR labelling | **Done.** CLR is printed as "algorithmic", with projection beside it; `print()` prints a status line for each fit (0.3.1) | — |
| Coverage claim scoped by regime | **Mostly done.** The status line names P3.5(b2)(ii), R1 vs R2, and "no coverage theorem" for grouped M-SIGN | an R2-regime line (effective *n* for the focal group, the finite-sample width beside DD), Part 3.0 |
| "Rejection" wording | **Done** (0.3.0). An empty estimate is a sample diagnostic; "rejected" is kept for the calibrated emptiness test | — |
| Integer tables | **Partly.** 2×2 margins-only proved (2c Remark M2.1); R15 open | an integer check in the count layer (Step 1) |
| M-SIGN scope | **Done.** EXPERIMENTAL lines, pooled-only proofs, W2 described as a pooled sum, pairwise-CR label names both sides | R16 |
| Lead with the realised-electorate target | **Partly.** `ci_actual` is intersected with DD, and `proj_actual` rests on P3.5(c1). But `print()` puts `estimate`/`ci_model` first and `duncan_davis` last, and no restriction acts on the actual table | Step 1 and decision D1 |
| Count-table interface | **Not done** (paper §7, "Not implemented") | Step 1 |
| M-SHARP specification | **Partly.** NEXT_BUILD now says "tighter than Θ_O, not sharp" when size is informative | Step 8: a narrower scope and a different place in the order |
| Short proved chains instead of one master ladder | **Not done** | `ladder()` prints proved implications only (Part 3.0) |
| ORD-Δ labelled "MTR" | **Not done** (PACKAGE_DESIGN §7 still says "M-ORD-Δ (MTR)") | rename to "ordering of changes" (Part 3.2) |

---

## Part 3. The recommended steps

**Review process for every step (your rule, D7).** Each module gets **two separate reviews**:
1. a **proof review** of its theory and proof items, before its status changes from "experimental" to "proved";
2. a **code review** of its R implementation against the theory document, before release.

This replaces the earlier batch C/D/E framing for future work. The batches already planned for existing material (D: PACKAGE_DESIGN and NEW_AVENUES; E: the 0.3.x code) are folded in as the proof and code reviews of the existing modules.

| Step | Proof review covers | Code review covers |
|---|---|---|
| 1 M-COUNT | R19; realised rows from 2c M2 (Batch B comment 1) | `eiCounts()`, integer/MILP switch, witnesses, the `object =` vocabulary, the combination rule |
| 1b estimated margins | R25, R26 (nondifferential matching) | margin sets, Charnes–Cooper/profiling, the Berkson correction |
| 2 DiD theory | PHASE2e: R11a–d, the gap register | `eiChange()`, `eiATT()` (after the proof review) |
| 3 M-CONTRAST | R22 | contrast rows in both layers; bin prespecification |
| 4 B1(M)+κ | R8 (B1 rank, monotonicity) | z handling, range-form M, κ annotation |
| 5 M1(δ) neighbours | R30 (neighbour form: nesting, components; R16-type sharpness) | adjacency/distance input, graph checks |
| 6 CAR(Z)-δ | R8 (CAR), overlap factor | covariate cells or η(Z), benchmarks, overlap printout |
| 7 M-EXCEPT | R23 | MILP frontier |
| 8 M-SHARP diagnostic | R24, R9′ | `eiSharpCheck()` |
| 9 M-BERN | outer/inner validity, dimension formula | Bernstein certificates, grids |
| 10 C2 contests as districts | R29 (identification of a_j, b_j; sharpness; the test) | contest stacking, Wald/wild-bootstrap test |
| 11 E4 | R21, the allowance | table rows, leave-one-out driver |

### 3.0 Step 0: finish what is in flight (S; no new theory)

* **Plain language.** Release v0.3.1 as a tag and get the code read by someone outside the project. Then make three small changes to what `print()` says, so that what the package claims matches what is proved.
* **Contents.**
  1. Tag `v0.3.1` (NEXT_BUILD item 3 says the tag is still to follow).
  2. Under your two-review rule (D7), the planned batches become reviews of existing modules: batch E becomes the **code review** of the 0.3.x code (M-SIGN, WS5), and batch D the **proof review** of PACKAGE_DESIGN's results (P-S, P-R, P-N, G2). Choose reviewers now; the count layer reuses the LP code that E covers.
  3. In regime R2 (many small groups), print the effective number of units for the focal group, $1/\sum_iw_{ir}^2$, and the width of the finite-sample interval next to the DD width. This adds information without changing the default.
  4. Have `ladder()` print only proved implications (for example M1 ⊂ M1δ ⊂ M2 and "signed slopes ⇒ weighted-local") and list the other rungs as side by side.
  5. Rename "M-ORD-Δ (MTR)" in PACKAGE_DESIGN §7 to "ordering of changes". MTR is a causal assumption.
* **Worked example.** A fit with 93 ridings of 10 polls (regime R2). *Today* `print()` shows the CLR interval with the line "algorithmic", and nothing tells the reader that the proved alternative is uninformative here. *With Step 0* it adds, for instance: "effective number of units for veterans: 41; finite-sample interval [0.08, 0.97] vs Duncan–Davis [0.10, 0.95]". The reader then sees directly that, in this regime, only the algorithmic interval is informative. (Illustrative numbers.)
* **Why first.** It is cheap, it uses no new theory, and the review batches it schedules take weeks to come back, so they should start now. Batch D's gate is met but it is unscheduled; batch E has no reviewer; batch C waits for the DiD data decisions (D4).
* **Proofs**: none. **Tests**: snapshot tests of `print()`.

---

### 3.1 Step 1: M-COUNT, the count-table (realised-electorate) layer

**What it does (plain language).** When a researcher asks "how many of these veterans voted yes?", the honest answer starts from arithmetic. In each precinct the veterans' yes votes must fit the precinct's totals, and nothing else is assumed. M-COUNT computes that range for any question that adds up votes:
* a group's rate;
* a gap between two groups;
* a rate for one county, or for precincts above some size.

Beliefs about the people who actually voted then go in directly, for example:
* "in every precinct, veterans were at least as supportive as other men";
* "no Democrat voted yes";
* "this township's poll book shows 12 veteran yes votes".

Each one narrows the range, and the package checks that whole-number vote tables can still produce the answer. No straight-line model, no regression and no sampling theory is involved. If the counts are right, the range is exact.

**Information or assumption used.** Only the precinct margins, plus whatever counts or beliefs about realised votes the user supplies.

**The question it answers.** "Given these returns, and given what I am willing to believe about these particular voters, what could the true count or rate be?"

**Worked example.** Two precincts:

| | Men | Veterans | Yes votes |
|---|---|---|---|
| P1 | 100 | 20 | 30 |
| P2 | 100 | 60 | 50 |

Target: the veterans' yes rate across both precincts (80 veterans).

* **Without the module.** The package prints a DD range for the actual rate, [0.125, 0.875]. In P1 the veteran yes count can be anything from 0 to 20. In P2 it is at least 10 (only 40 non-veterans can account for the 50 yes votes) and at most 50. The model rungs say nothing useful with two precincts. A belief such as "veterans were more supportive *in every precinct*" can enter today only as a statement about the model's averages.
* **With the module.** The same DD range comes first. The belief "in each precinct veterans voted yes at least as often as others" then raises P1's floor to 6 and P2's floor to 30, so the rate is in **[0.45, 0.875]** (N6, N7). The output also shows a witness table: an actual pair of precinct tables with whole numbers that reaches each end.

**Why now.**
* Most later modules (descriptive DiD, comparisons, exceptions, and the general-user features: audits, known tables, ballot patterns) are a few rows or columns added to this layer.
* It is the base of the DiD layer: the exact DiD rungs (Δ0, Δ4) are realised-table results.
* Its theory is done: DD, 2c Prop. M2 for realised signs and CR, and PHASE2d Δ0/Δ4.
* It answers the question applied users ask first.
* It needs only aggregate margins, which is all many applications have.
* **Cost:** moderate code (a new LP builder, targets, witnesses and an optional MILP), one small proof item, no data.

**Technical notes.**
* **Object.** Unknown counts $n_{irc}\ge0$ with the observed row and column margins in each unit. The target is any $\sum_i c_i^\top n_i$, divided by a known denominator (or by a sum of known denominators). Restrictions are linear rows on the $n$:
  * structural zeros and known cells;
  * published totals over unit sets (§3.7c; a general-user feature);
  * the realised unit ordering $n_{irc}N_{is}\ge n_{isc}N_{ir}$ (N1/B5);
  * the realised association signs $\delta_B,\delta_W$ of 2c M1 (rows in $n$);
  * CR as a union of two LPs;
  * slack dials, as in 2c M2(d);
  * the odds-ratio dial of NEW_AVENUES §12.3(c) for 2×2 and pairwise 2×C: each unit's group-by-outcome odds ratio in $[\psi_L,\psi_U]$, a unit-level interval because the cell is monotone in the odds ratio. Ordering is $\psi_L=1$. In P2 of the example, an odds ratio between 1 and 4 gives the veteran rate [0.500, 0.633] instead of [0.500, 0.833] (N17). Label it untestable when the dial contains 1.
* **Proof items.**
  * Proposed **R19**:
    * (a) the bound is the LP over the product of unit transportation polytopes, intersected with the rows. This is standard.
    * (b) For margins, known cells and per-unit box rows the constraint matrix is totally unimodular. The LP endpoints of integer-coefficient count targets are integers once the box bounds (and any totals) are rounded to integers first; rounding the LP endpoint afterwards is not enough (PHASE6a R19(b), check T5). Checked in N10a (300 tables, 0 fractional optima). *Corrected after PHASE6a (proof reviewed).*
    * (c) Rows that couple cells across units or within a unit (pooled realised signs, ordering with fractional coefficients) can break this. Then the LP endpoint is an outer bound, and an integer program decides. N10b gives a one-unit fractional endpoint. N10c gives a three-unit system that is LP-feasible but has **no** integer table, which matches Batch B comment 1.
  * Realised signs and CR for continuous tables: 2c M2 (proved). The integer version: the default one-vs-rest (W1) signs and CR are rows on the single aggregate $Y_{rc}(\mathcal S)$, i.e. totals, so they are integer-exact by the rounded LP (PHASE6a Lemma C1, Prop. C2(b)). MILPs are needed for pairwise $\delta_B(r,s)$, W2, rows on cells in different rows and columns ($R\ge3$ and $C\ge3$) and mixed row sets (PHASE6a §3.4); for those, R15 remains a computed check rather than a theorem. *Corrected after PHASE6a: an earlier draft said all pooled realised signs and CR break total unimodularity.* Any "never refutable" label must say "continuous tables" (Batch B comment 1).
  * **Review gates** (from the PACKAGE_DESIGN author, §4.13). The realised-sign and CR rows rest on 2c M2, where Batch B comment 1 is still being worked through. The DiD rungs (Step 2) rest on Δ0/Δ4, which are unreviewed (Batch C). Those rows carry "experimental" until the batches resolve. Margins, known cells and totals do not.
* **Combination rule with the LCE layer.** Rows on counts and rows on the model's coefficients never share one LP. Under LCE, a finite realised table puts no per-unit constraint on the counts, because LCE restricts only a conditional mean. The actual-rate answer is the **intersection** of the count-layer set (exact) with the LCE interval for $\tilde B$ plus P3's allowance (level $1-\alpha$). That intersection covers $B$ with probability at least $1-\alpha$. If the two are disjoint, this **rejects LCE (with its rows) together with the count rows**: at level α exactly when the count rows are identity rows (margins, published totals, audited cells); with belief rows (ordering, signs, floors) the probability of rejecting while those rows are true is at most α (PHASE6a R27). A free test for PACKAGE_DESIGN §4.3. Footnote f of §4.2 becomes this rule instead of a warning.
* **One restriction vocabulary.** "Veterans more supportive" can mean a realised sign (2c M2, rows on counts) or a model sign (M5, rows on $k$). The API uses one name with an `object = c("actual", "model")` tag. It refuses to impose both under one label, and `meta$object` drives the combination rule. The R2 face test carries over unchanged, so `eiRedundancy()` extends to the count layer.
* **Inference status.** None is needed when the margins are exact counts. If a margin is an estimate (matched counts, survey margins), the exact claim fails. Step 1b handles this with margin sets and a stated coverage. `integer = "auto"` switches the integer check off for non-integer margins, and `print()` says which case applies.
* **Interface sketch.**
  ```r
  fit <- eiCounts(formula = cbind(yes, no) ~ cbind(vet, other), data = d, unit = "precinct",
                  targets = list(rate("vet", "yes"), gap("vet", "other", "yes"),
                                 rate("vet", "yes", where = county == "Dane")),
                  restrictions = list(order_unit = list(row = "vet", over = "other", col = "yes"),
                                      known = list(unit = "P7", row = "vet", col = "yes", value = 12),
                                      total = list(units = c("P1","P2"), row = "vet", col = "yes", value = 45),
                                      sign = list(type = "dB", sign = ">=0", object = "actual"),
                                      cr = list(object = "actual")),
                  integer = c("auto", "lp", "exact"))   # auto: MILP only when rows break TU
  print(fit)        # target, lower, upper, width, integer status, binding rows
  witness(fit)      # unit tables attaining each endpoint
  ```
  Also (decided, D1/D2): `eiBounds()` embeds `eiCounts()` and prints its actual-rate block first, Duncan–Davis first. `target = "model"` restores the old order.
* **Tests.**
  * Reproduce `eiBounds()`'s `dd_lower`/`dd_upper` to 1e-12 for single cells and one-column contrasts. For other combination targets, including the one-vs-rest gap $D^r_c$ with $R\ge3$, the package's cellwise interval arithmetic is outer and `eiCounts()` is correctly narrower (PHASE6a Lemma 2.1(c)).
  * Brute-force enumeration of all integer tables in tiny examples.
  * Closed forms of 2c M2 and M3 in 2×2 (reuse the C7b/C8 checks in `verify_phase2c.py`).
  * Every witness satisfies every row.
  * The N10c system is reported as integer-infeasible.
  * Restrictions only narrow the bound.

---

### 3.1b Step 1b: estimated margins (probabilistically matched group counts): M-MISCLASS moves into your track

**The problem (plain language).** Group counts are sometimes not counts. In the Wisconsin illustration, the number of veterans in a township It is an *expected* count from matching soldier lists to the census, so each man carries a probability of being a veteran. Two consequences follow:
* The package must accept non-integer margins, and the integer check must switch off for them.
* More fundamentally, the arithmetic bounds are exact only if the margins are right. An estimated margin is not right; it is uncertain.

Treating the expected count as exact gives a range that is too narrow and can miss the truth.

**What the module does.** It replaces "the township had 60 veterans" by "the township had between 56 and 64 veterans, with 95% probability, given the match probabilities". The arithmetic bounds are then computed over every table consistent with any count in that range. The answer becomes a confidence statement: "with probability at least 95%, the true veteran yes rate lies in this range". It holds *if* the match probabilities are calibrated (a man given 0.5 is a veteran half the time) and independent across men. If calibration is doubtful, a slack dial widens the margin sets. An audit of some matches, if one ever exists, would narrow them. This is the misclassification problem (formerly G3) with the error law supplied by the matcher. It moves from the general-user track into the core plan, because estimated margins are common. The calibration slack and the match probabilities are user inputs; choosing them for a given project is outside this plan.

**Worked example** (N20). Precinct P2: 100 men, 50 yes votes. Matching makes 50 men certain veterans and gives 20 men probability 0.5 each, so the expected count is 60.
* *Treating 60 as exact* (what M-COUNT would do without this step): the veteran yes rate is in **[0.167, 0.833]**.
* *With the module*: the true count lies in [56, 64] with probability 0.95, so the rate is in **[0.107, 0.893]**. The range is wider, and it is honest about the matching.

**Why here (Step 1b, immediately after Step 1).**
* It decides whether Step 1's "exact" claim holds whenever a group margin is estimated (matched records, survey-based margins), and it changes the DiD theory in Step 2, where a matched composition enters both periods.
* It is mostly linear programming over margin sets plus a probability calculation. The theory is short but new (below).
* **Cost:** M. New proof items R25 and R26. Data: the match probabilities (or at least, per township, the number of certain matches and the probabilities of the uncertain ones).

**Technical notes.**
* **Two inputs, two treatments.**
  * (i) *Non-integer margins as given* (point estimates only). `eiCounts()` accepts them. The integer check switches off automatically (`integer = "auto"` detects non-integer margins). `print()` says the range is exact only for the stated margins, which are estimates.
  * (ii) *Match probabilities available.* For each unit, a set $\mathcal V_i$ for the true group count. One option is a central Poisson-binomial interval at level $1-\alpha_i$. The other is a joint set across units: $\sum_i$-type constraints with a Bernstein or Hoeffding radius, which is tighter for aggregate targets than per-unit intervals plus a union bound.
* **The bound.** Optimise over counts and margins jointly: $n_{iVc}\le V_i$, $n_{iVc}\ge Y_{ic}-(N_i-V_i)$, $V\in\mathcal V$. The target "rate among true veterans" $\sum n/\sum V$ is linear-fractional (Charnes–Cooper, or profile over integer $\sum V$). Count targets stay linear. The true $V_i$ are integers, so the integer check applies to the *latent* counts, not to the expected ones.
* **Proof items.**
  * **R25** (new): (a) if $P(V\in\mathcal V)\ge1-\alpha$, the bound covers the actual target with probability at least $1-\alpha$. The argument is one line: the count-layer bound is exact given $V$. (b) The joint-set construction and its level. (c) Monotonicity in a calibration slack $\varepsilon$, with probabilities moved by up to $\varepsilon$.
  * **R26** (new): the LCE layer with matched composition. If matches are calibrated and independent, the true share is $X=\bar q+u$ with $E[u\mid q]=0$ (Berkson error) and $\operatorname{Var}(u\mid q)=\sum_jq_j(1-q_j)/N^2$. Then $E[T\mid q]=\bar q^\top Q\bar q+\operatorname{tr}(Q\Sigma_u)$: a **known** correction that involves only the identified symmetric part $Q$, so $K$ remains the only unknown (checked, N19). The weights of the model-implied target also use the true $X$, which is uncertain, and that needs its own allowance.
  * Both rest on a **nondifferential matching** assumption: given the true composition, the match probabilities carry no information about how anyone voted. If the matcher used variables that predict votes, this fails, and the theory must say so explicitly.
* **Relation to G3/M-MISCLASS.** It is the same object, counts by true group × recorded label × outcome, with the error law given by the match probabilities rather than by an audit. The general-user version, audit bounds without probabilities, stays in the general track as a variant of the same code.
* **Inference.** Coverage comes from the margin set (R25). With the LCE layer, intersect as in Step 1's combination rule, splitting α between the margin set and the LCE interval.
* **Interface.**
  ```r
  eiCounts(..., margins = list(row = "vet", certain = d$vet_sure, probs = match_probs_by_unit,
                               level = 0.95, joint = TRUE, calib_slack = 0))
  eiBounds(..., composition_error = list(probs = match_probs_by_unit))   # LCE layer, R26 correction
  ```
* **Tests.**
  * N19 and N20.
  * Exact Poisson-binomial coverage by simulation.
  * Exactness when all probabilities are 0 or 1.
  * Monotone widening in the slack.
  * Brute force over small tables.
* **Reviews**: proof review (R25, R26), then code review.

### 3.2 Step 2: a general DiD theory track (descriptive and causal), then implementation

**What it is (plain language).** Before writing code, write down and prove a general theory of difference-in-differences with ecological data. It is meant to work for any units, groups, outcomes and pair (or sequence) of elections, not only Wisconsin.
* **Descriptive part.** How much more did one group's support change than another's, among the same people? Bounds come from arithmetic, from sign beliefs on changes, from comparisons between kinds of places, and from straight-line models of changes.
* **Causal part.** What would the treated group have done without treatment? Each causal assumption is stated as an explicit, adjustable belief: a bound on how different the untreated trends could have been, a relative-magnitude rule from pre-periods, or the monotone-response assumptions.
* **The gap register.** People enter, leave and age into the electorate between elections. You asked to defer this. So the theory marks **every place** where an assumption about the changing population, or a data-based bound on it, has to be added, and leaves those places explicitly open.

You will supply data later. Implementation follows once the theory has passed its proof review.

**Worked example** (arithmetic only, unchanged).
* The descriptive layer finds veterans' support rose 5 to 25 points more than others'.
* If untreated trends could differ by up to 10 points, the effect of service is −5 to +35 points. If by at most 2, it is +3 to +27.
* *Without the theory*, users combine these by hand and can mix descriptive and causal claims. PACKAGE_DESIGN §7's "MTR" label for orderings of changes is one example of that mix.

**Why here.** DiD is central to the applications in view, and you asked for general theory and proofs first. The theory reuses PHASE2d, which is written but unreviewed, and generalises it. Doing it second lets it build on Step 1's count layer and Step 1b's estimated margins (a matched composition enters both periods). **Cost:** theory L; implementation M.

**Deliverables.**
1. **Document** `docs/extension/PHASE2e_did_general.md`, built from PHASE2d and generalising it. Contents:
   * **Setting**: units, R groups, C outcomes, **T ≥ 2 elections with pre-periods** (your decision D11); a fixed population as the base case.
   * **Descriptive rungs**:
     * Δ0 (product of periods);
     * Δ4 (signs and CR on changes) beyond 2×2;
     * finite contrasts on changes (Step 3's changes version);
     * FD (first-difference LCE) and DiD-S, with their sharpness;
     * estimated margins in both periods (Step 1b).
   * **Causal layer**, with **all three** forms developed (D11):
     * a trend dial: $ATT=D-b$ with $\lvert b\rvert\le M$;
     * bounds relative to pre-period trends (Rambachan–Roth relative magnitudes, and smoothness), over **one feasible sequence** of latent pre-trends taken jointly with the post-period violation. The precise treatment is in "Set-valued pre-trends" below;
     * MTR with ND and PS (2d Prop. CL);
     * spillovers onto the comparison group as an explicit sensitivity term;
     * the rule that a causal label requires a stated comparison group and trend assumption.
   * **Gap register** (deferred, named): entry; exit; aging-in; stayers vs leavers; interval bases from censuses; changing unit boundaries. For each, the theorem that holds without it, the assumption or data bound that would close it, and the rung it would modify.
2. **Proof items.** R11 split into:
   * R11a: sharpness of the descriptive rungs, Δ4 beyond 2×2, DiD-S within FD;
   * R11b: validity of the causal layer for each trend dial;
   * R11c: inference for unions across periods and for FD;
   * R11d: integer versions (with R15).
   * Plus R25/R26 applied to two periods.
3. **Checks.** A `verify_phase2e.py` that:
   * ports 2d's 15 checks;
   * brute-forces small two-period integer tables for every descriptive rung;
   * checks the causal layer's closed forms against LPs;
   * confirms each gap-register item fails as stated when its assumption is dropped, each with a check that fails on the naive statement.
4. **Reviews**: an external proof review of PHASE2e (replacing batch C), then implementation (`eiChange()`, `eiATT()`) and a separate code review.

**Set-valued pre-trends: why relative-magnitude bounds may be uninformative, and how the theory must handle it.** *Your caution, made precise.*

*Plain language.* In ordinary DiD, the pre-period trends are estimated as numbers. Rambachan–Roth then say: "the post-period violation of parallel trends is at most M̄ times the largest pre-period violation". Here the pre-period trends are themselves EI bounds, ranges rather than numbers. The rule must then allow for the *largest* pre-period trend the data cannot rule out. When the pre-period ranges are wide, that allowance swamps everything.

*Theory to write (PHASE2e).*
* **Joint object.** The latent tables of all periods determine every DiD at once. Let $\mathcal D$ be the joint set of $(D_{\text{pre}},D_1)$ they allow. It is not a box: adjacent DiDs share a period's latent table.
* **Treatment of the pre-periods.** With no treatment before $t=1$, $D_t=b_t$ for $t\le0$. The relative-magnitude set is
  $$\{D_1-b_1:\ (D_{\text{pre}},D_1)\in\mathcal D,\ \lvert b_1\rvert\le\bar M\max_{t\le0}\lvert D_t\rvert\},$$
  optimised jointly. Taking each pre-period's extreme separately can pair incompatible histories, so it is weakly wider; the review insists on the joint form.
* **Degeneration** (a proposition to prove).
  * (a) The relative-magnitude bound is contained in the trend-dial bound with $M=\bar M\cdot\sup_{\mathcal D}\lvert D_{\text{pre}}\rvert$. It equals it when the maximisers are compatible, which is the endpoint-wise case.
  * (b) It is uninformative, covering the whole logical range, once $\bar M\sup\lvert D_{\text{pre}}\rvert$ exceeds the gap between the post-period DiD range and the logical limits.
  * (c) It is informative only when a rung narrows the pre-period DiDs: a model rung (FD/R1), Δ4 signs, or contrasts on changes. With point-identified pre-trends it reduces to standard Rambachan–Roth.
  * (d) Sampling uncertainty in model-based pre-trends enters through one joint confidence region, not period-by-period intervals.
* **Proof items.** R11b′: the joint-set formula and (a)–(d). R11c: inference with a joint region.

*Worked example* (N21; two precincts, three elections, arithmetic only).
* The pre-period DiD is only known to lie in **[−1.23, 1.27]**, and the post-period DiD in [−1.08, 1.42]. Both are wide, because Duncan–Davis alone is wide here.
* With M̄ = 1, the ATT bound is:
  * **[−2.35, 2.65]** jointly;
  * [−2.35, 2.69] endpoint-wise.

  Both extend beyond the logical range of a DiD of rates ([−2, 2]), so the result is **no information**.
* A plain trend dial with M = 0.10 gives [−1.18, 1.52]. This is wide too, but only because the post-period DiD itself is wide.
* If a model rung pinned the pre-period totals to within ±2 votes, the pre-period DiD would shrink to [−0.06, 0.10]. The relative-magnitude bound becomes [−0.60, 0.90]: the allowance falls from ±1.27 to ±0.10.

The lesson for the plan: relative-magnitude bounds are worth computing only on top of rungs that narrow the pre-period DiDs, and `eiATT()` should print the implied M = M̄·sup|D_pre| beside them so the reader sees when they collapse to a trend dial or to nothing.

**Technical notes.**
* **Descriptive.** The rows on changes are count-layer rows on the two-period product, with exact coverage. FD and DiD-S are LCE-layer rows on changes. They combine by intersection (Part 3A).
* **Causal.** The ATT set is the descriptive set shifted by the trend dial's set. With relative magnitudes it is the projection of a joint LP over latent pre-trends and the post-period violation. With MTR/ND/PS it is 2d Prop. CL.
* **Interface** (after the reviews): `eiChange(d0, d1, unit = "township", restrictions = list(cr = list(object = "actual", on = "changes")))`; `eiATT(change_fit, trend_gap = 0.1, pre = NULL, mtr = FALSE)`.

### 3.3 Step 3: M-CONTRAST, finite comparisons between observed kinds of places (M-LOC folded in)

**What it does (plain language).** Researchers often believe something about two *real* kinds of places, not about a curve. Take places with few veterans and places with many. Many researchers would say two things:
* veterans were at least as supportive as other men on average;
* both groups were at least as supportive where veterans were common.

Together these two directional beliefs **cap** the gap between groups at the observed difference in overall support between the two kinds of place. A plain "veterans were more supportive" belief gives no cap at all. The module needs no straight line between the two kinds of place and no slope at any composition the data never show.

**Information or assumption used.** Two directional beliefs about a prespecified comparison (which places count as "low" and "high", and which groups are compared).

**The question it answers.** "If veterans were more supportive and support rose where veterans were common, how big can the veteran gap be?"

**Worked example.** Two bins of 100 men each:

| Bin | Veteran share | Yes rate |
|---|---|---|
| Low | 25% | 45% |
| High | 75% | 55% |

The target is the average veteran–other gap across the two bins.

* Arithmetic alone: the gap is in [−0.60, 0.73].
* "Veterans more supportive" alone raises the floor to 0 but leaves the ceiling high.
* With both beliefs, the 10-point difference between bins is split into a "who lives there" part and a "how everyone votes there" part. Each part must lie between 0 and 10 points, so the gap is in **[0, 0.20]** (N12b).

**Why here.**
* It is the most useful assumption-based module and costs little.
* Its identity is algebra (N12a). Under LCE the same rows are exactly the existing sandwich objects at the midpoint of the comparison (N12c). M-LOC is therefore the special case "every comparison along the hull", and its proof debt (C1) concerns only the implication to $\delta_W$, not its computation.
* It works both in the count layer (bins of actual voters) and under LCE.
* It needs no data beyond the margins, so it applies to any aggregate data. Its version on *changes* (comparing low- and high-veteran townships' changes) is the finite form of DiD-S (PACKAGE_DESIGN §7) and needs no linearity in either period.
* It comes after Steps 1–2 because those are exact; this one adds beliefs.
* **Cost:** S–M, one small proof item.

**Technical notes.**
* **Identity** (Kitagawa 1955):
  $$\tau_z-\tau_x=(z-x)^\top\bar m+\bar x^\top(m_z-m_x).$$
  For $z-x=t(e_r-e_s)$, the two signs confine each term to $[0,\tau_z-\tau_x]$.
* **Count layer.** Bins are user-chosen sets of units. They must be **prespecified**, or reported as a family of beliefs. Several contrasts intersect, so trying bins until the bound is narrowest is belief-shopping, even though no inference is involved (PACKAGE_DESIGN author, §4.13). At `at = "hull"` the module inherits the proved sandwich (P-S, R12). Bin-level group rates are linear-fractional in counts with known denominators, so they are linear. The rows are exact for those bins' actual voters. In R≥3, use only comparisons that the user names. A general displacement has a well-defined mix term, but attributing it to one pair is an assumption.
* **LCE layer.** The mix term equals $t\,O_{rs}(\bar x)$ and the rate term equals $t\,L_{rs}(\bar x)$ (N12c). Options:
  * `at = "pairs"`: chosen midpoints;
  * `at = "hull"`: every point, which is M-LOC (P-S, R12).
* **Proof item.** Proposed **R22**: (a) the identity (checked); (b) the LCE midpoint form (checked); (c) the realised bin version as LP rows (trivial); (d) relation to M-SIGN: M-LOC on the hull implies the midpoint rows, but not conversely.
* **Inference.**
  * Count layer: exact.
  * LCE pooled M0: linear rows in $k$, so R3 and R10 apply (sharp; projection proved, CLR algorithmic).
  * Grouped rungs: R16, as for M-SIGN, and experimental until then.
* **Interface.**
  ```r
  restrictions = list(contrast = list(rows = c("vet","other"), col = "yes", from = d$vetshare < .3,
                                      to = d$vetshare > .6, mix = ">=0", rate = ">=0", slack = 0))
  eiBounds(..., restrictions = list(contrast = list(..., at = c("pairs", "hull"))))   # "hull" = M-LOC
  ```
* **Tests.**
  * N12a–c.
  * Brute force on bins.
  * Rows only narrow the bound.
  * The truth is covered when the signs hold, and is excluded at the expected rate when they are false (the 0.3.1 false-belief design).
  * The test for the sign of the observed contrast ($\tau_z-\tau_x<0$ contradicts the two signs jointly).

---

### 3.4 Step 4: M-CTX2, context at a second scale (B1(M)), with κ-annotated endpoints

**What it does (plain language).** A township's own veteran share has effects the data cannot assign to veterans or to non-veterans. That is the core ambiguity of EI. The *county's* veteran share is different: its effect on each group's rate **can** be read off the data. B1(M) assumes the township's own make-up moves each group's rate at most *M* times as much as the county's make-up does. The researcher then reports results across a range of *M*. Beside it, each end of every bound is labelled with κ, the amount of contextual variation it requires ("this end needs group rates to swing by at least so much with composition"). Readers see what an extreme answer demands without anyone committing to a ceiling.

**Information or assumption used.** A composition measured at a second scale (the county), and a judgement about how large own-township effects can be relative to the identified county effects. The κ annotation needs no new assumption; the κ *cap* is an assumption with a dial.

**The question it answers.** "If context inside a township matters no more than M times what the county-level pattern shows, what is the veteran gap?" And: "how much contextual variation does each end of my bound require?"

**Worked example** (illustration on the Wisconsin data; [NA, `apply_B_suffrage`, `apply_kappa_suffrage`]; trusted; plug-in, no confidence intervals). Wisconsin 1865, trimmed units, veteran minus non-veteran yes rate.
* *Without the module*: pooled LCE gives **[−0.066, 0.452]**, which does not settle the sign.
* *With B1(M = 1)* on the yes column: **[0.049, 0.297]**, a positive gap. With M = 0.5: [0.111, 0.235].
* *κ annotation* (same data, pooled LCE; see the κ table in §4.8(a) for the cap's effect): the lower end of −0.065 needs κ above 0.25. Capping κ at 0.25 lifts the lower end to 0.019. The data's own floor is κ = 0.015.
* *A warning that belongs in the output*: imposed on all three columns jointly, B1 is empty at the plug-in unless M ≥ 3.85 (1857) or 6.65 (1865). Iowa 1868 is empty even under plain M0. This is why the step includes the emptiness test and proper intervals, not plug-in sets.

**Why here.**
* It is NEW_AVENUES' first choice for near-term use (§9), and it applies to any data with a composition measured at a second scale (in the Wisconsin/Iowa illustration, the county share).
* It gives the DiD its yardstick. M-CTX2-Δ replaces FD3's free δ with a measured contextual effect (PACKAGE_DESIGN §7). It is also the benchmark of the right kind for any CAR-style sensitivity (NEW_AVENUES §14.5).
* It comes after Steps 1–3 because those are exact, or use only directional beliefs. This step adds a magnitude belief on the LCE layer, and it needs a proof first (R8; until then the module carries "experimental").
* Calibration. NZ (Step 11) is the only data with known truths, and the NEW_AVENUES author argued that it should calibrate these dials first. You moved E4 to last, because most applications will lack the tables it needs. So Step 4 ships with its dials uncalibrated: M_min and κ annotation, reported with intervals, plus breakdown values. When Step 11 is built, its NZ output (the κ, δ and M values the truths needed) becomes a reference distribution that is printed beside these dials.
* **Cost:** M. The B1 rank condition (R8), the parametrisation-invariant M, κ's ellipsoid in R×C, and inference through the existing stack.

**Technical notes.**
* **Model.** $E[\beta_{rc}\mid x,z]=e_r^\top H_cx+e_r^\top\Gamma_cz$; $\Gamma$ is identified up to the intercept normalisation (NA Prop. B1, a sketch). B1(M): $\text{Own}_{rc}\le M\,\text{Area}_{rc}$, with ranges over the *observed* supports (§2.3, the invariant form). These are LP rows. B2 (proportionality) and B3 (no own effect) are points; report them only as diagnostics (heavy-tailed ratio, §2.4).
* **κ.** The annotation computes, at each LP endpoint, the implied $\operatorname{Var}_X(E[\beta_r\mid X])/\operatorname{Var}(Y)$; this is cheap. The cap is LP rows in 2×2 and an ellipsoid (SOCP) in R×C. Data-implied κ never tightens (§12.1).
* **Proof items.** R8, split into: (a) the B1 rank condition, with $z$ not affine in $x$ on the support; (b) monotonicity in $M$ and κ (elementary); (c) the invariant M's LP form.
* **Inference.** The rows enter the CLR and projection stacks. Pooled M0 is covered by R3/R10 (linear rows in $k$ and $\Gamma$); grouped rungs follow R16. $M_{\min}$ (the smallest M with a non-empty set) is reported with the emptiness test, not as a plug-in.
* **Interface.**
  ```r
  eiBounds(..., z = ~ county_vetshare, restrictions = list(ctx2 = list(M = c(0.5, 1, 2), form = "range", cols = "yes")))
  sensitivity(fit, dial = "M", M = seq(0, 4, .25)); breakdown(fit, dial = "M", target = "gap", threshold = 0)
  print(fit, annotate = "kappa")   # each endpoint with the kappa it requires
  ```
* **Tests.**
  * Reproduce `sim_two_scale_context` (identification column) and the Wisconsin tables from `apply_B_suffrage`. The data are not in the repository, so this must run with the owner's copy.
  * Truth excluded at λ = 1.5, M ≤ 1.
  * Nesting in M.
  * κ annotation agrees with `apply_kappa_suffrage`.

### 3.4b Step 5: M1(δ) with named neighbours (bounded differences between neighbouring districts)

**What it does (plain language).** With grouped data (districts, counties), the existing M1(δ) rung lets each district's hidden contextual pattern differ from one common pattern by at most δ. The neighbour option says something more local and often more believable: **neighbouring districts' patterns differ from each other by at most δ**. The user names the neighbours, either with an adjacency matrix ("these counties share a border") or with a distance matrix ("patterns may differ more the further apart two districts are"). The payoff is that a district whose own data say little borrows strength from a well-identified neighbour, not from an average of everyone.

**Information or assumption used.** A geography of the districts and a judgement δ about how different neighbours can be.

**The question it answers.** "If neighbouring counties' contextual patterns differ by no more than δ, what is this county's group rate?"

**Worked example** (N23). Three districts on a line, A – B – C. On its own (M2), A's hidden deviation could be anywhere in [−0.10, 0.10], B's is well identified at [−0.02, 0.02], and C's is in [−0.12, 0.12].
* **The existing centre form**, M1(δ = 0.05), leaves A at **[−0.10, 0.10]**: no gain, because the common centre can sit anywhere B allows.
* **The neighbour form** with adjacency A~B, B~C and δ = 0.03 gives A **[−0.05, 0.05]**: A is tied to its well-identified neighbour.
* A distance matrix with the Lipschitz default (λ = 0.03 per unit distance) gives the same here.
* If A has no neighbour (a disconnected graph), A stays at [−0.10, 0.10].

**Why here (Step 5).** It is a variant of a rung that already exists, so its theory is a small extension of 2b. It is placed after Step 4 because both concern how context varies across space, and the two must be specified together when both are used (see below). **Cost:** S–M; proof item R30.

**Technical notes.**
* **Object.** District deviations $v_d=k^{(d)}-c_d$ (2b §0.2). The neighbour rows are $\lVert v_d-v_e\rVert_\infty\le\delta_{de}$ for each named pair. They are linear rows in the stacked $k$ (the same LP as M1(δ)).
* **Adjacency input**: $\delta_{de}=\delta$ for adjacent pairs.
* **Distance input.**
  * **Default: Lipschitz**, $\delta_{de}=\lambda\,d_{de}$ for **all** pairs. Reasons:
    * it has one dial (λ, in probability units per unit distance);
    * no threshold has to be chosen;
    * it never disconnects the graph;
    * λ = 0 gives M1(0) and λ → ∞ gives M2;
    * it is consistent with the metric, so for nearby pairs only the tightest constraints bind.
  * **Alternatives**: a radius ($\delta_{de}=\delta$ if $d_{de}\le r$, otherwise none), which reduces to an adjacency graph; or a user function $\delta_{de}=\delta\,f(d_{de})$ with $f$ non-decreasing and $f(0)\ge0$.
* **Relation to the centre form.** If the neighbour graph is connected with diameter $D$, the neighbour rows imply the centre form with $\delta_c=\delta D/2$ (coordinate-wise midpoints). Conversely, the centre form with $\delta_c$ implies the complete-graph neighbour form with $2\delta_c$. So the two forms are nested within a factor of the diameter. They can be imposed together; this is an intersection in one LP. Report which one binds.
* **Incoherences and warnings.**
  * **A disconnected graph**: each component drifts independently. Across components there is no restriction (M2 across components) unless the centre form is also imposed. Isolated districts are unrestricted. `print()` reports the components.
  * **Drift along chains**: with δ per edge, a long chain allows a total drift of up to δ × (path length). This is the district-level analogue of NEW_AVENUES §12.1 rule 2. The Lipschitz default bounds it through the metric.
  * **With B1(M)**: B1 needs its second-scale composition to vary *within* districts (Step 4, footnote 3 of Part 3A). If districts are the same units as B1's areas (counties with county share), B1 and district rungs cannot be combined at all, neighbours or not. With districts coarser than B1's areas, both are rows in one LP. B1's county effect Γ is shared, and the neighbour rows act on the district deviations; that combination's status is open.
* **Scope**: grouped rungs only (M1(δ) and its variants). It is meaningless for pooled M0.
* **Status.**
  * Validity as an outer bound and monotonicity in δ or λ: elementary (Prop. D6 extended).
  * Sharpness for the regression: R16-type, open (2b Thm 1′ with these rows).
  * Inference: projection or CLR in the few-large-districts regime (R1). There is no theorem for many small districts.
* **Proof item R30**: (a) nesting and monotonicity; (b) the diameter relation with the centre form; (c) component structure; (d) sharpness via 2b Thm 1′ (open).
* **Interface**:
  * `eiBounds(..., model = "M1delta", neighbours = list(adjacency = A, delta = 0.03))`;
  * `neighbours = list(distance = D, lambda = 0.03)` (default form);
  * `neighbours = list(distance = D, radius = r, delta = 0.03)`;
  * `neighbours = list(distance = D, decay = f, delta = 0.03)`.

  Combine with the centre form via `delta = ` as now.
* **Tests**:
  * N23;
  * δ = 0 on a connected graph equals M1(0);
  * δ → ∞ equals M2;
  * the diameter inequality;
  * a disconnected graph gives per-component results;
  * a symmetric A is checked, and an asymmetric one is refused.
* **Reviews**: proof review (R30), then code review.

### 3.5 Step 6: CAR(Z)-δ, covariate adjustment with a dial (shares Step 4's model)

**What it does (plain language).** Among places alike on covariates (for example wealth, nativity, denomination), assume the veteran share tells you nothing more about how each group voted. That is "no contextual effect once the covariates are held fixed". The dial δ lets that fail: group rates may still depend on composition by up to δ. At δ = 0 the answer is a point (Goodman within covariate cells). As δ grows it widens back toward the no-model bounds. The package reports the whole curve, a benchmark for δ from each covariate ("dropping wealth would create a violation of size δ_wealth"), and an overlap warning: where composition varies little within covariate cells, every δ is magnified.

**Information or assumption used.** Township covariates, with veteran share varying within covariate cells (positivity), and a judgement about δ.

**How it relates to B1(M)+κ (Step 4).** Both relax the same strong claim, that a place's own composition has no effect once something else is held fixed. They differ in what is held fixed and in how the dial is benchmarked:
* **CAR-δ** conditions on covariates and benchmarks δ against what dropping a covariate would do. NEW_AVENUES §14.4 shows such covariate benchmarks are **optimistic about contextual effects**, by up to 12× in its check.
* **B1(M)** conditions on composition at a second scale and benchmarks M against the *identified* county-level contextual effect: a benchmark of the same kind as the threat.

They combine in **one LP** on the shared model "each group's mean rate = effect of covariates + effect of own composition + effect of county composition", with both caps imposed. If the county share is itself one of the covariates, CAR(Z) is B1(M = 0), so report them as one model with two dials, not as two independent checks. NEW_AVENUES §14.7 recommends reporting CAR-δ curves, the overlap factor, covariate benchmarks and B1(M) side by side, and that is what Steps 4–6 produce together.

**Worked example** [NA, `check_car_sensitivity`, reproduced].
* A contextual violation that leaves the regression exactly unchanged biases a group's rate by 0.080.
* A covariate-benchmarked sensitivity bound would report 0.006–0.032 for the same bias.
* With compositions confined to [0.3, 0.7] within covariate cells, a δ of 0.05 becomes a half-width of 0.125 at the vertex rate (overlap factor 2.5).
* *Without the module*: no covariate adjustment at all, or a point estimate (Goodman within cells) that silently assumes δ = 0. *With it*: an interval that widens honestly with δ and with poor overlap, beside B1(M)'s like-for-like benchmark.

**Why here.** Covariates are commonly available in applications of these models, so it is a real candidate, not a rare case. The covariate set is a user input. It shares Step 4's model, so it comes right after it. It is placed after B1(M) because its benchmark is the weaker one for EI's signature threat. Covariate adjustment was earlier low priority by your decision. Covariate availability makes it worth a numbered step, but it can trade places with Step 7 if you prefer. **Cost:** M (LP rows per covariate cell or a linear η(Z); the overlap factor; the benchmarks).

**Technical notes.**
* **Model.** $E[\beta\mid X,Z,z]=\eta(Z)+\Delta$, with $\lVert\Delta\rVert_\infty\le\delta$ (LP) or $E[\Delta_j^2]\le\kappa_j\operatorname{Var}(Y)$ (SOCP). With linear η, CAR(Z) is the LCE+Z model with the own-composition effect set to zero (NEW_AVENUES §13.4).
* **Proof items.**
  * R8 (CAR part): validity and monotonicity in δ.
  * The invisible-violation class in R×C (§14.3; proved in its simple form).
  * The overlap factor beyond 2×2 (a conjecture).
  * The benchmark $\delta_k$ and its caveats (§14.5).
* **Inference.** η(Z) is estimated, by series or within cells. Rows enter the CLR and projection stacks for linear η. Semiparametric nuisance estimation is a Phase-3 question (§13.4). Positivity diagnostics are printed.
* **Do not build** the Cinelli–Hazlett/Chernozhukov contours. Point to `seine` with §14.4's caution.
* **Interface.** `eiBounds(..., covariates = ~ Z1 + Z2, restrictions = list(car = list(delta = seq(0, .1, .01), norm = "sup", bench = list(kappa = TRUE, m = c(0.5, 1, 2)))), z = ~ area_share)`; `print(fit, overlap = TRUE, benchmarks = TRUE)` prints κ_k, δ_k, $R^2(X\sim Z)$, $R^2(Z_k\sim X\mid Z_{-k})$, the overlap factor and, when z is given, the B1 benchmark beside them.
* **Tests.**
  * Reproduce `check_car_sensitivity` Parts 1–3.
  * δ = 0 equals Goodman within cells.
  * δ → ∞ equals DD.
  * Monotone in δ.
  * B1(M = 0) coincides with CAR(Z) when z ∈ Z.
* **Reviews**: proof review, then code review.

**Covariate-benchmarked κ and δ (a general part of this step, built on the CAR(Z)-δ construction above; kept as written by your decision).** Following NEW_AVENUES §14.8, at your suggestion.
* **Construction.** Under CAR(Z), dropping covariate $k$ creates an *identified* violation of CAR($Z_{-k}$): $\Delta^{(-k)}(x,z_{-k})=E[\eta(Z)\mid X=x,Z_{-k}]-E[\eta(Z)\mid Z_{-k}]$. Its size gives two benchmarks: $\kappa_k=E[(\Delta^{(-k)}_j)^2]/\operatorname{Var}(Y)$ (a variance share) and $\delta_k=\sup\lvert\Delta^{(-k)}_j\rvert$ (sup norm). A joint benchmark $\kappa_Z$ drops all covariates at once.
* **The bridge assumption, in plain language.** "Whatever contextual violation remains after adjusting for the covariates is at most *m* times what covariate *k* accounts for": $\kappa\le m\,\kappa_k$, or $\delta\le m\,\delta_k$. **m is a judgement.** The benchmark only sets the scale on which that judgement is expressed. It does not identify m: the multiple the truth actually needs, $m^*$, is never observed.
* **When it is informative: your point, made precise.** Three printed quantities decide it.
  1. **The size of the benchmark, κ_k (and δ_k).** A small κ_k at a believable m gives narrow bounds. A large κ_k gives a large allowance even when it is conservative.
  2. **How much of the composition the covariates absorb**, $R^2(X\sim Z)$, and for each covariate $R^2(Z_k\sim X\mid Z_{-k})$. When covariates absorb most composition variation, little composition variation is left for a contextual remainder. Its variance share is then small. This is your "remainder from contextual bias may be highly limited", and it holds **in variance terms**.
  3. **The overlap factor.** The same high $R^2(X\sim Z)$ means composition varies little *within* covariate cells. Each group's rate is then an extrapolation to an all-one-group place, and any remainder is magnified on the way. So a small variance share does **not** imply a small error in the group rate.
* **Small examples** (N22; one covariate, an invisible contextual violation of known size).
  * *Covariate shares composition, small contextual effect* ($R^2(X\sim Z)=0.26$). The truth needs only a tiny fraction of the benchmark ($m^*\approx0$), and the CAR(Z) answer for group 1 is off by only −0.026. The benchmark is safe and informative if κ_Z is modest.
  * *Covariate absorbs most composition variation* ($R^2(X\sim Z)=0.92$). The remainder's variance share is small, κ_true = 0.006. But the overlap factor doubles (2.4 → 5.0), and the CAR(Z) answer for group 1 is still off by **−0.123**. Moreover κ_Z = 0.99, so a bridge at m = 1 allows almost any violation: safe but uninformative.
  * *Covariate explains the outcome but not composition* ($R^2(Y\sim X,Z)=0.88$, $R^2(X\sim Z)=0.00$). κ_Z ≈ 0, so any m calibrated to it allows nothing, and the truth needs **m* ≈ 2,800**. This is the §14.4 understatement case in its extreme form.
* **Reconciling with §14.4** (benchmarks can understate contextual bias by up to 12×).
  * Understatement happens when the true violation acts **through composition itself** (an own-composition contextual effect) and the covariates are only weakly tied to composition. κ_k is large only for covariates that both move rates *and* are associated with composition given the other covariates. No observed covariate is as collinear with composition as composition itself.
  * The three cases separate two notions that are easy to conflate:
    * covariates that explain a lot of the **outcome**: $R^2(Y\sim X,Z)$ was 0.88–0.92 in all three cases, and it predicted nothing about $m^*$;
    * covariates that explain a lot of the **composition**: $R^2(X\sim Z)$. This one matters, but in two opposite directions (point 2 against point 3).
  * Outcome R² cannot bound the contextual remainder at all, because that remainder is invisible in $E[Y\mid X,Z]$ (§14.3).
* **Relation to Step 4's B1 benchmark.**
  * B1(M) benchmarks the own-composition effect against the *identified* effect of composition at another scale, a benchmark of the **same nature** as the threat.
  * Covariate benchmarks calibrate omitted-covariate-type violations, and they are optimistic for contextual ones.
  * **Show both whenever a second-scale composition exists.** When they disagree, the B1 allowance should govern statements about contextual bias.
  * Without a second-scale composition, print the covariate benchmark with the §14.4 caution and the overlap factor beside it.
* **Proof items** (proposed **R28**):
  * (a) identification of $\Delta^{(-k)}$ under CAR(Z) with positivity;
  * (b) validity of the bridge bound and monotonicity in m, with its LP (δ) and SOCP (κ) forms;
  * (c) the group-rate error bound in terms of κ and the overlap operator norm in R×C, generalising §14.5's 2×2 factor;
  * (d) a statement that $m^*$ is not identified, so the bridge is untestable except through emptiness.
* **Checks.**
  * N22, extended to several covariates and to R×C.
  * Coverage of the bridge bound at m ≥ m* and failure at m < m*, by simulation with a known truth.
  * The bridge bound collapses to CAR(Z) at m = 0 and to DD as m → ∞.

### 3.6 Step 7: M-EXCEPT, how many exceptions can a conclusion survive?

**What it does (plain language).** "Veterans were more supportive in every precinct" is a strong claim, and one odd precinct can break it. M-EXCEPT replaces it with "the ordering held except in precincts holding at most a share π of voters". It then reports, for each π, whether the conclusion still holds, and the smallest π at which it fails. Reviewers ask how robust a conclusion is to exceptions, and this answers them directly.

**Information or assumption used.** The ordering belief from Step 1, weakened by allowing exceptions in at most a share π of voters (and, optionally, by a tolerated size *s* of violation).

**The question it answers.** "How many exceptional precincts would it take to overturn my conclusion?"

**Worked example.** The two precincts again, with the floor of the veterans' yes rate:
* No exceptions: 0.45.
* Either precinct may be the exception (π = 50%): 0.20, because the worst case is P2 misbehaving.
* Anything goes: 0.125 (N13).

So the claim "the veterans' rate exceeds 40%" survives only if no precinct holding half the voters is exceptional. A uniform five-point relaxation cannot represent "one badly reversed precinct among many orderly ones".

**Why here.**
* It needs Step 1's ordering rows and a MILP.
* It gives the most useful robustness display for the realised beliefs.
* It is not needed before those beliefs exist. It applies equally to the sign beliefs on changes in Step 2 ("veterans' support rose more than others' except in townships holding at most π of the men").
* **Cost:** M (MILP with indicators, plus a frontier grid). Proof: small.

**Technical notes.**
* Rows: $g_i\ge -s-z_i$ (big-M on the realised gap), $z_i\in\{0,1\}$, $\sum_ip_iz_i\le\pi$.
* The frontier is monotone in $(\pi,s)$, because the feasible set grows in both. It is a step function in π for finite units.
* A convex alternative bounds $\sum_ip_i(-g_i)_+$, a different belief (average severity).
* Proposed **R23**: monotonicity; exactness of the MILP; the minimum exception mass implied by the margins (a data fact, not the true prevalence).
* **Inference**: deterministic with exact counts. With estimated margins, a joint region projected over the grid; pointwise intervals are not a frontier.
* **Interface**: `frontier(fit, target, threshold = .4, pi = seq(0, 1, .05), s = c(0, .05, .1), weight = c("voters","group","precincts"))`.
* **Tests**: N13, monotonicity, brute force over subsets in small designs.

---

### 3.7 Step 8: M-SHARP, narrowed to a benchmark and a diagnostic

**What it does (plain language).** The straight-line model is currently fitted only to the *average* vote share at each composition. M-SHARP also makes the line agree with what each precinct's arithmetic allows, *on average across precincts with the same make-up*. Two kinds of place give it something to work with:
* places where similar precincts vote very differently;
* places where small and large precincts differ.

For the actual voters it adds nothing: precinct-by-precinct arithmetic already uses every return. It can only narrow the *model-based* answer.

**Information or assumption used.** No new assumption: the whole distribution of returns at each composition (and size), instead of only its average.

**The question it answers.** "Does the full spread of returns tell the straight-line model anything that the averages do not?"

**Worked example.** At a 50/50 composition, suppose half the precincts give "yes" 10% and half give 90%.
* The average (50%) only says group 1's mean rate there is between 0% and 100%.
* Precinct by precinct, group 1's rate is at most 20% in the first kind of precinct and at least 80% in the second. Its mean is therefore between **40% and 60%** (N5).
* Size can matter in the same way. Small precincts that split 0%/100% at a 50/50 composition force the mean to 50%. If rates do not depend on size, that carries over to large precincts, whereas pooling sizes would give [25%, 75%] (N4).

**Why it comes last of the numbered steps.** (Its diagnostic can be demonstrated on the NZ data, which have many small booths.)
* In an oracle study that knows the law exactly, with precinct rates spread *unimodally* around their means, it narrowed the model's group rates by a median of 0–1% and at most 5%:
  * 2×2: 36 designs, median ratio 1.00, minimum 0.95;
  * 2×3/2×4: 24 designs, minimum 0.97–0.98 with 40,000 draws per type.
* In those designs the exact subset rows added nothing over averaged cell limits.
* **High-kink-mass designs** (added after consultation, §4.13). With *bimodal* precinct rates (Beta concentration 0.5–1), the ratio is 0.63–0.92 (24 cells; median 0.84). Small precincts (10–50 voters) with unimodal rates gain only 1–8%. So the gain is real exactly when precincts are polarised within groups, which is what the kink-mass diagnostic detects.
* With 8,000 draws per type the same study showed *spurious* gains of up to 8%. These come from sampling noise in the right-hand sides being read as information, which previews the inward bias any plug-in version will have.
* So the estimator costs a lot (smoothing, size bins and a new inference argument), and its payoff depends on the data. Build the **diagnostic** now and the **estimator** only if the data call for it. The diagnostic is the "kink mass" of NEW_AVENUES §6, a bimodality check on returns at similar compositions, and needs only Step 1's per-unit limits. **Run it on NZ and Wisconsin as soon as Step 1 exists.** If either shows substantial kink mass, promote the estimator (with $(X,N)$ conditioning, R9′) ahead of Step 9.
* **Cost:** M for the oracle and the diagnostic; L for a released estimator, which is not recommended now.

**Technical notes.**
* **Subset identity** (review Part II.1). For 2×C at a fixed $v$, the conditional mean masses $\bar y$ satisfy
  $$\bar y(A)\le E[\min\{x,T(A)\}\mid v]\ \text{for every proper subset }A,\qquad \bar y(\mathcal C)=x.$$
  This is exactly the Aumann expectation of the unit sets.
  * **Proof** (proposed **R24**): $S(x,t)$ is the base polytope of the submodular function $A\mapsto\min\{x,t(A)\}$ (a concave function of a modular one). The support function of a base polytope is the Lovász extension of its function, which is linear in the function. So the expected support function is the Lovász extension of $E\min\{x,T(A)\}$, a submodular function whose base polytope is the set above. R×2 is the same with $\min\{x(A),T\}$. Verified: N1, N2 (1,600 support-function comparisons, max error 4e-16).
  * Averaged cell limits are exact for C ≤ 3 and can be strictly looser for C ≥ 4 (N3a/b). The review's 2×4 example is confirmed (N3c).
  * For general R×C, use the lifted LP with one table per atom (N1b).
* **Conditioning on size**: condition on $(x,\text{size bin})$; pooling is a labelled outer relaxation (N4).
* **Inference, if ever built.** The review's averaged-instrument Hoeffding construction gives a finite-sample outer procedure. Its union-bound validity is standard and taken on trust here; it is not re-proved. R9 should be restated as (a) outer validity of that procedure and (b) a gain study. Drop the series-consistency target.
* **Interface**: `eiSharpCheck(fit, by = c("x_n", "x"), bins = ..., oracle = FALSE)` (conditioning on $(X,N)$ by default, per Batch B comment 4). It also reports the binomial-slope check (does the variance of returns scale like 1/n at fixed composition?), a diagnostic suggested by the NEW_AVENUES author (§4.13), not a rung reports the kink mass and the population (or binned) Θ_X/Θ_O width ratio per target. There is no change to `eiBounds()`'s default.
* **Tests**: N1–N5; recovery of the oracle on replicated designs; the Prop. 3 equality cases; the 2×4 coordinate failure.

---

### 3.7b Step 9: M-BERN, the polynomial robustness ladder (with L1′)

**What it does (plain language).** Every LCE answer assumes each group's rate changes along a straight line as composition changes. M-BERN lets the lines bend (quadratic, cubic) and shows how far the answer moves. It never narrows anything. It tells readers how much of a result rests on straightness. The package computes two bounds that bracket the true curved answer: an outer one (rates checked on a fine grid) and an inner one (exact Bernstein certificates).

**Information or assumption used.** None beyond a chosen maximum degree, which must be fixed in advance (2c §6).

**The question it answers.** "How much of my LCE result comes from assuming straight lines?"

**Worked example** (N16). Five kinds of precinct; target the group-1 rate (truth 0.632).
* *With straight lines only* (today): [0.436, 0.779].
* *With the ladder*:
  * quadratic: [0.419, 0.796];
  * cubic: [0.359, 0.820];
  * Duncan–Davis: [0.324, 0.834].

A reader sees that a cubic bend gives back most of the distance to Duncan–Davis on the lower end.

**Why here.** It is the credibility check for every LCE number, but it adds no tightening. It therefore comes after the steps that add information or beliefs. It shares its machinery with L1′ (quadratic contextual effects), which is on PHASE4's second-release list. **Cost:** M (the outer/inner LP pair, the dimension formula and degree elevation); no new proof beyond NEW_AVENUES §4.2's design.

**Technical notes.**
* Degree-$J$ homogeneous effects leave $d_J=(C-1)[R\binom{R+J-1}{J}-\binom{R+J}{J+1}]$ free directions (NA §4.2, proved by counting).
* Outer bound: validity on a grid or at subdivision vertices. Inner bound: non-negative Bernstein coefficients after degree elevation; Pólya convergence holds only for strictly positive polynomials. Report both, and the gap.
* Inference: the rows enter the CLR and projection stacks as for LCE.
* Interface: `ladder(..., rungs = c("M0", "BERN2", "BERN3", "DD"), bern = list(grid = 81, elevate = 4))`.
* Tests:
  * N16;
  * nesting in $J$;
  * outer ⊇ inner;
  * the gap shrinks with grid refinement and elevation;
  * reproduce LCE at $J=1$.

### 3.7d Step 10: C2, contests as districts (several contests on the same ballot)

**What it does (plain language).** On one ballot, voters often decide several contests: governor and a referendum, president and senate, several propositions. Each contest is a separate EI problem with the same groups. On its own, knowing how contests combine does not help (NEW_AVENUES C1, proved). C2 adds one assumption: **context works the same way in every contest, up to a level and a scale**. If high-veteran townships pull veterans' support up in one contest, they do so in the others, more or less strongly. Contests with lopsided results, where some group votes almost unanimously, then pin down the hidden pattern, and that information carries over to the close contests. The assumption is **testable**: it implies that the identified parts of the contests line up, which the data can contradict.

**Information or assumption used.** Returns for two or more contests in the same units, plus the assumption that contextual structure is shared up to level and scale.

**The question it answers.** "If context shapes these contests alike, what do the other contests tell me about this one?"

**Worked example** (N24; two contests, two groups, five kinds of precinct; numbers illustrative). Contest 1 alone gives group 1's rate **[0.498, 0.840]**. Contest 2 has the same pattern scaled by 1.1 and shifted down by 0.12. Under C2, the scale is identified from the two contests' identified parts, and contest 2's own limits carry over. Group 1's rate in contest 1 narrows to **[0.539, 0.798]** (truth 0.715). The gain grows with more contests and with more lopsided ones.

**Why here (Step 10, just before E4).** It is general (any ballot with several contests), and it reuses the M1 machinery with contests in place of districts. It comes late because it rests on a substantive cross-contest assumption, and the case for it varies by application: for example, governor and suffrage questions may not share contextual structure. It precedes E4 by your decision. **Cost:** M; proof item R29.

**Technical notes.**
* **Model** (NEW_AVENUES §3.2 C2): $H^{(j)}=a_jJ+b_jH^{(1)}$ per chosen column. Then $Q^{(j)}=a_jJ+b_jQ^{(1)}$ (checked, N24) and $K^{(j)}=b_jK^{(1)}$. The bound for $K^{(1)}$ is $\bigcap_j\{K:Q^{(j)}+b_jK\ \text{valid}\}$: LP rows once $a_j$ and $b_j$ are fixed by the identified $Q$'s.
* **Columns.** Contests with different outcome sets need a declared column correspondence, for example "Republican or yes". By default C2 is per declared column pair.
* **Testability.**
  * Each extra contest adds the restrictions $Q^{(j)}=a_jJ+b_jQ^{(1)}$ beyond its two parameters. The count of overidentifying restrictions per contest is the number of free entries of $Q$ per column minus 2, summed over columns.
  * With three or more contests, a pooling-type Wald test across contests (wild bootstrap, as for M1) has more power, and the emptiness test applies to the intersection.
  * A rejection rejects C2 (with LCE), not the contests' own LCE fits.
* **Identification caveats.** $b_j$ is weakly identified when $Q^{(1)}$ has little curvature beyond its level. A contest with no contextual structure ($b_j\approx0$) adds nothing. A negative $b_j$ (context reverses across contests) is allowed but should be reported.
* **Relations.**
  * **Nesting (Step 2).** Nesting ("yes ⇒ Republican, except ε") is a count-layer restriction on the *joint* table of contests in each unit. C2 is a model-layer restriction on contests' *mean* structures. They combine by intersection (Part 3A) and do not double-count.
  * **M-CVR (G2).** Cast-vote records give observed joint patterns and bound group-specific cross-contest *differences* exactly, without C2. Where CVRs exist, C2's implied group-specific differences, $m^{(j)}_r-m^{(1)}_r=a_j+(b_j-1)m^{(1)}_r$, can be compared with the CVR bounds: a further check of C2.
* **Proof item R29**:
  * (a) identification of $a_j$ and $b_j$ from the $Q$'s, with a rank condition;
  * (b) sharpness for the regression, by the 2b Thm 2 intersection argument (claimed in NEW_AVENUES, to write);
  * (c) the overidentification count and the test's calibration (wild bootstrap; simulation-supported, as for the M1 test);
  * (d) inference for a fixed number of contests (fixed dimension, so projection applies).
* **Status**: open until R29 is reviewed. Grouped versions follow R16.
* **Interface**:
  * `eiBounds(list(gov = f1, ref = f2), data = d, contests = list(columns = c(gov = "R", ref = "yes"), shared = "level_scale"))`;
  * `tests$contests` reports the pooling test.
* **Tests**:
  * N24;
  * recovery on noise-free designs with known $a_j$, $b_j$;
  * the test's size under C2 and its power when one contest's pattern is changed;
  * C2 with one contest equals plain LCE.
* **Reviews**: proof review (R29), then code review.

### 3.7c Step 11: known tables (M-DATA): E4 transfer (validated on NZ); exact uses (E5) as general-user features

**What it does (plain language).** Sometimes more is published than precinct totals:
* a county-wide count of veterans who voted yes;
* veterans and votes both broken down by age group within each precinct;
* an electorate-level cross-tab (New Zealand).

The **exact** uses add these facts as equalities on the actual votes, with no model, so they never make a bound wrong if the published numbers are right. A separate **transfer** use assumes that the hidden pattern behind a known table also holds in areas without a table, and uses that to narrow bounds there. Transfer is an assumption. It comes with a dial for how different the areas may be, and New Zealand lets us test it.

**Worked examples.**
* **A published total.** The county reports that 45 of its 80 veterans voted yes (P1 and P2 together). P2's own veteran yes rate narrows from [0.167, 0.833] to **[0.417, 0.750]** (N8b).
* **Stratified margins.** P2 reports veterans and yes votes by age:

  | P2 band | Men | Veterans | Yes votes |
  |---|---|---|---|
  | Young | 30 | 25 | 10 |
  | Old | 70 | 35 | 40 |

  Arithmetic within each band gives P2's veteran yes count in **[10, 45]** instead of [10, 50] (N9). This is covariate adjustment with no assumption at all. It serves the covariate wish on your list at almost no cost.
* **Transfer (E4), demonstrated on NZ-type data.** Two electorates under the straight-line model (two groups): electorate A has precincts with 10% and 30% veterans, electorate B has precincts with 60% and 90%. *Without the module* the model-implied veteran rate in B is in **[0.529, 0.751]**. *With A's published table*, and the assumption that the hidden pattern is the same in both, B's rate is pinned at **0.656** (with two groups one table fixes the pattern completely, NEW_AVENUES E1). *Allowing the pattern to differ by an amount worth 2 points on A's rate*, B's rate is in **[0.650, 0.662]** (N15). How large a difference was needed in past NZ elections is exactly what leave-one-electorate-out measures; NEW_AVENUES §5.2 suggests it may be large, in which case the gain shrinks accordingly.

**Why now.**
* **For you**, the value is E4 on NZ: NZ publishes every electorate's table, so leave-one-electorate-out tests whether transfer of the hidden pattern works, and gives the model layer its first real validation. It needs a proof (rank) and the actual-vs-model allowance, hence its place after the cheaper steps.
* **For other users**, the exact parts (published totals, stratified and crossing margins) are almost free once Step 1 exists: a `total` row is just an equality, and strata are extra unit-level polytopes. You have no such tables, so ship them as general-user features with Step 1's generic rows, not as a separate step.
* E4 is the repository's own top data priority (NEW_AVENUES §13.8).
* **Cost:** S for E5 (mostly in Step 1), M for E4 (proof plus the leave-one-out benchmark).

**Technical notes.**
* **E5 exact.**
  * Totals over unit sets are equality rows on cells.
  * Stratified margins (NEW_AVENUES §13.2(a)) become extra unit-strata with their own transportation polytopes, linked by the target.
  * Crossing partitions (§13.2(b)) become a latent table over atoms with both partitions' margins.
  * All of these are LPs in the count layer. Totals spanning several units couple them, so the R19(c) integer note applies.
* **E4 transfer.**
  * In the *model* layer, a table is an equality on $B$, not on $\tilde B$, so it needs a residual allowance (NEW_AVENUES §13.1, caveat 1). Use P3.5(c1) with the projection interval.
  * Rank: NEW_AVENUES §13.1 table, checked numerically. It needs a proof stated in terms of the composition second-moment matrices: proposed **R21**.
  * Show which target directions $a$ lie in the row space of the stacked table map, because those targets are point-identified even when $K$ is not (review Part I.6).
  * Dial: $\lVert K_j-K\rVert\le\delta$ reuses M1δ.
* **Inference.**
  * E5: none (exact).
  * E4: equality rows enter the existing CLR and projection stacks. Coverage follows R10 for pooled M0. For grouped transfer, R16-type caveats apply.
* **Interface.**
  ```r
  eiCounts(..., restrictions = list(total = list(units = d$county == "Dane", row = "vet", col = "yes", value = 45)),
           strata = "age")                            # long data: unit x stratum
  eiBounds(..., restrictions = list(table = list(units = ..., B = known_table, share_K = TRUE, delta = 0.05)))
  eiTransfer(tables, holdout = "Wellington Central", delta = seq(0, .2, .02))  # NZ leave-one-out
  ```
* **Tests.**
  * N8b and N9.
  * Brute force on small strata.
  * NZ leave-one-electorate-out: coverage, width and the δ needed for coverage, by election.
  * The overidentification test from shared $K$ (NEW_AVENUES §13.1).

---

### 3.8 General-user track: modules for data you do not have

*You have told us there are no poll books, linked records, CVRs or sub-tables for your projects. The modules below need such data. They are worth keeping because other users have it (archives with poll books, US cast-vote records, linked historical records), and each is small once Step 1 exists. Build them after Step 7, or earlier if a user with the data asks; G1 before G2 by effort. G3 moved to your track as Step 1b. The exact E5 rows (totals, strata) belong here too, but ship inside Step 1 (§3.7c).*

#### G1. M-AUDIT, choosing which precinct to validate (general-user feature)

**What it does (plain language).** Archives sometimes let you recover the full table for a few precincts, from poll books linked to service records or from a validation study. That is expensive. M-AUDIT tells you **in advance** how much each precinct's table would narrow your answer, and it guarantees that narrowing whatever the table turns out to show. You can then spend the budget where it counts. Big precincts dominated by the group you study are often a poor choice: their arithmetic already pins the answer.

**Information or assumption used.** The margins, plus the cost of recovering each precinct's table. The planner assumes the recovered table is complete and correct for the same voters.

**The question it answers.** "Which precincts should I pay to validate, and what will I get?"

**Worked example.** The same two precincts:
* P1 can move the veterans' yes count by 20 votes, out of a total range of 60.
* P2 can move it by 40.

Recovering P2's table removes exactly 40 of the 60 votes of uncertainty (the rate interval shrinks from 0.75 wide to 0.25 wide), whatever P2's true table is. Recovering P1 removes 20. With equal costs, audit P2 (N8a: in 200 random designs the gain equals the precinct's own range exactly).

There is one complication. If a published total links the precincts, the gains stop adding up. Take a published county total of 45 veteran yes votes and the target "P1's rate". Auditing P2 then pins P1 completely: it removes all 20 votes of P1's range, though P2 is not in the target. Without the total it would remove nothing (N8c).

**Why here (in the general-user track, not Step 2).**
* **It serves other users' data, not yours**: no poll books or linked records exist for Wisconsin/Iowa.
* For users who can recover tables it changes what they do in the archive, and it is a thin layer on Step 1 whose main fact is elementary.
* **Cost:** small code, no new data beyond costs, one small proof item.

**Technical notes.**
* **Separable case** (margins plus per-unit rows): the width equals $\sum_id_i$, with $d_i$ the unit's own range for the target. Auditing $i$ removes exactly $d_i$. Equal costs: rank by $d_i$. Unequal costs: a 0–1 knapsack.
* **Coupled case** (totals, pooled realised signs): use the worst-case remaining width. For each candidate unit $j$, the remaining width is the difference between a max-LP value (concave in unit $j$'s revealed count) and a min-LP value (convex in it), so it is concave in that count. The worst case is then a one-dimensional concave maximisation per unit. Checked in N8d: 0 violations in 1,640 checks.
* An expected-width criterion needs a prior over the truth. Offer it only as an option, with the prior stated.
* **Proposed R20**: (a) separability; (b) concavity of the remaining width in the revealed count; (c) audits of several units at once: greedy is not optimal in the coupled case, so state the exhaustive or MILP version.
* **Inference**: deterministic.
* **Interface.**
  ```r
  eiAudit(fit, target = "rate(vet, yes)", cost = d$cost, budget = 3,
          criterion = c("guaranteed", "worst_remaining", "expected"), prior = NULL)
  ```
* **Tests.**
  * N8a–N8d.
  * Hide-and-restore on the NC data that has known truth (`senc`).
  * Compare with largest-precinct and random audit plans on simulated designs.

---

#### G2. M-CVR, ballot patterns across contests (general-user feature)

**What it does (plain language).** Many US jurisdictions publish cast-vote records: for each ballot, the choice in every contest. Knowing that no ballot split its ticket tells you that **every group voted the same way in both contests**, even though the records still cannot say how any group voted in either one. M-CVR bounds group-specific *differences* across contests: split-ticket voting, roll-off and down-ballot defection.

**Worked example.** One precinct of 100 voters, 50 in each of groups A and B, two yes/no contests.
* If 50 ballots are "no–no" and 50 are "yes–yes", group A's support in contest 1 can be anything from 0% to 100%. But its support in contest 1 minus its support in contest 2 is **exactly 0**.
* With 4 ballots "yes–no", the difference is in **[0, 0.08]** (N11).
* Without the module, the package would treat each contest separately and say nothing about the difference.

**Why here (general-user track).**
* With Step 1 in place it is a change of columns and nothing else.
* It is a strict gain, and it leaves the proved negative result (C1: joint margins do not help a single contest's cell) intact.
* **It serves other users' data**: you have no CVR application, but public CVRs exist (Kuriwaki et al. 2024) and split-ticket questions are common.
* **Cost:** S; data are public (Kuriwaki et al. 2024).

**Technical notes.**
* The columns are ballot patterns. Targets are linear contrasts of pattern cells (for example $n_{A,10}-n_{A,01}$).
* It is a transportation LP, so it has integer vertices.
* A partially recorded precinct gets an "unknown pattern" block.
* Proofs: none beyond R19.
* Interface: `ballotPatterns(d, contests = c("pres","sen"))` builds the columns, then `eiCounts(..., targets = list(diff_contests("A", "pres:yes", "sen:yes")))`.
* Tests: N11, brute force, and reconciliation of CVR totals with the official returns.

---

#### G3. M-MISCLASS: moved to your track as Step 1b

Your Wisconsin veteran counts come from probabilistic matching, so uncertain group labels are your problem, not only other users'. The module is now Step 1b (§3.1b). The audit-bounds variant, error counts bounded by an audit rather than given by match probabilities, stays here as a general-user option of the same code (N14 is its example).

* **R16/R15 (M-SIGN out of experimental).** Unchanged. They are on the theory side, they run in parallel, and they are not a build step.

---

## Part 3A. Which modules combine, and how (the layering matrix)

*Your question: "Can M-CONTRAST and the context approaches be layered with the LCE / M-SIGN approaches at all?" Short answer: **yes, almost everything combines**, in one of two ways, and a few combinations are incoherent or double-count. The rule underneath comes from the PACKAGE_DESIGN author (§4.13).*

**Plain language.** The package has two layers.
* The **count layer** reasons about the actual voters' table: realised beliefs, exact arithmetic.
* The **model layer** reasons about expected rates under a model: LCE and everything built on it.

Within a layer, beliefs go into **one optimisation** and are combined exactly. Across layers they are **never mixed inside one optimisation**: a belief about the actual table says nothing per precinct about the model's averages, and vice versa. Instead, each layer gives its own range for the actual rate, and the package reports where the two ranges **overlap**. That overlap is valid for the actual rate, with the confidence level of the model-layer interval. If they do not overlap at all, the model (with its beliefs) is rejected at that level.

**Codes used in the matrix.**
* **L**: combine in the same LP (same layer).
* **∩**: intersect the two answers for the actual target (cross-layer). The confidence level is that of the model-layer interval: projection proved for pooled M0 and few-large-groups designs, CLR algorithmic.
* **↑**: one is built on top of the other (post-processing).
* **alt**: alternative bases. Choose one, or report side by side; never intersect as if both held.
* **R**: redundant. The first implies the second, so the second adds nothing for any target (the redundancy diagnostic shows zero gain).
* **×**: incoherent or double-counting as usually stated (footnote).

Target: **a** = actual rate (exact, or with the stated coverage); **m** = model-implied rate. Status: **P** = sharpness and inference proved; **P0** = proved for pooled M0 only; **O** = open; **E** = exact (no inference needed); **—** = not meaningful.

**Count-layer block** (all exact for the actual target, given exact or set-valued margins):

| | M-COUNT (realised order/sign/CR, odds ratio) | CONTRAST-count (levels / changes) | M-EXCEPT | DiD descriptive (Δ0, Δ4) | E5 rows | estimated margins (1b) |
|---|---|---|---|---|---|---|
| **M-COUNT** | L·a·E (MILP for coupled rows) | L·a·E | L·a·E; **×** if a strict ordering and its exception version are imposed on the same pair¹ | L·a·E (two-period product) | L·a·E | L·a·coverage 1−α (R25) |
| **CONTRAST-count** | | L·a·E | L·a·E | L·a·E (changes version) | L·a·E | L·a·R25 |
| **M-EXCEPT** | | | L·a·E | L·a·E (exceptions to change orderings) | L·a·E | L·a·R25 |
| **DiD descriptive** | | | | L·a·E; DiD causal is **↑** on top | L·a·E | L·a·R25 (both periods) |

**Model-layer block** (model-implied target; actual target via P3's allowance):

| | M0 | M1, M1(δ), M2, L4 | L5/L6 | M-SIGN (signs, CR) | CONTRAST-LCE (pairs / hull = M-LOC) | B1(M)+κ | CAR(Z)-δ | M-BERN / L1′ | E4 | FD/DiD-S (changes) |
|---|---|---|---|---|---|---|---|---|---|---|
| **M0** | · | alt (ladder) | alt | L·m·**P0** (R3, R10) | L·m·**P0** (R3: linear rows in k) | L·m·O (R8; R3 for the B0 model unproved) | L·m·O | alt (M-BERN ⊃ LCE) | L·m·P0 + allowance | — (FD is the changes analogue) |
| **M1…L4** | | · | alt | L·m·**O** (R16) | L·m·O (R16-type) | **×** if z is constant within groups³; else L·m·O | L·m·O (district dummies are a Z)⁴ | alt | L·m·O | — |
| **L5/L6** | | | · | **—** (point-identified: a sign is a check, not a restriction) | — | — | — | — | — | — |
| **M-SIGN** | | | | L (one LP; CR a union) | pairs: L·m·P0 (non-nested). **hull: R**, hull ⇒ δ_B and, in 2×2 with C1, δ_W⁵ | L·m·O (state the sign on the full mean incl. county term)⁶ | L·m·O (same caveat) | L·m·O (outer only) | L·m·P0 | Δ signs: L (model DiD) |
| **CONTRAST-LCE** | | | | | L | L·m·O | L·m·O | L·m·O (the midpoint identity N12c is LCE-only; secant rows stay valid) | L·m·O | DiD-S **is** the changes version: L·DiD·validity P, sharpness O (R11a) |
| **B1(M)+κ** | | | | | | · | L·m·O. **×** as "two independent checks" if z ∈ Z: then CAR(Z) = B1(M=0)⁷ | L·m·O (unexamined) | L·m·O (needs Γ shared across areas) | M-CTX2-Δ: L·DiD·O |
| **CAR(Z)-δ** | | | | | | | · | L·m·O | L·m·O | conditional parallel trends: O |
| **M-BERN** | | | | | | | | · | O | O |

**The two additions (Steps 5 and 10) in the matrix.**
* *M1(δ) neighbours* is a model-layer rung option: L with M1(δ)'s centre form (they nest within a factor of the graph diameter), L with M-SIGN and CONTRAST-LCE (status O, R16-type). With B1(M) it is × if the districts are B1's areas; otherwise L·m·O. With E4 it is L·m·O. It is not meaningful with pooled M0.
* *C2 contests* is model-layer: L with M0 and with M-SIGN (status O until R29); ∩ with nesting (count layer) and with M-CVR bounds (count layer). Stacking contests and districts at once is L but untested (O).

**Cross-layer rule** (every count-layer module × every model-layer module): **∩·a**. The count-layer range intersects the model-layer interval for the actual rate (model interval plus P3's allowance). Coverage is at least 1−α with identity count rows, and at least 1−α−P(belief rows false) otherwise. Disjointness tests the model and its rows together with the count rows, with the count rows maintained (not conditioned on): level α with identity rows only (PHASE6a R27). Status: inherited from the model-layer interval (projection: P0 or few-large-groups; CLR: algorithmic). Signs imposed in both layers: footnote 2.

**Footnotes: the incoherent or double-counting cases.**
1. **Strict ordering + M-EXCEPT on the same pair**: the strict row removes every exception, so the frontier collapses to π = 0. Impose one or the other.
2. **CR in both layers.** Realised CR (2c M2, count layer) and model CR (M5, model layer) are different beliefs, about actual voters vs expected rates. Imposing both is allowed, but it is **two** assumptions. The `object =` tag forces the user to say which, and `print()` lists them separately. Calling them one belief double-counts it. The same applies to every sign.
3. **B1(M) with district rungs.** The county share is constant within a county, so county dummies (M1/M2 with county groups) absorb it and the county effect Γ is unidentified. Use the pooled model, or groups coarser than counties, or a z that varies within groups (NEW_AVENUES §2.5).
4. **CAR(Z) with district dummies.** The dummies are a covariate. Including both a district rung and a district Z double-enters the same variable; treat districts as part of Z.
5. **CONTRAST "hull" + M-SIGN W1.** Hull contrasts imply the model between sign: δ_B(r,s) is a weighted average of the hull gaps with non-negative weights (checked, N18). In 2×2 they also imply δ_W ≥ 0, if Conjecture C1 holds (P-N). Adding the signs is redundant, not wrong. The "pairs" version is not nested with the signs, so that combination is a genuine intersection.
6. **M-SIGN with B1(M).** Once a county term enters the mean, "δ_B ≥ 0" must say whether it signs the full mean (own plus county term) or only the own-composition part. Both are linear, but they are different beliefs.
7. **B1(M) + CAR(Z) when the county share is a covariate.** Then CAR(Z) is B1(M = 0). Report one model with two dials (M for the county-benchmarked own effect, δ for the residual), not two independent robustness checks.

**The main combinations you would use, with small examples.**
* **Realised ordering + the LCE interval (∩).** In the two-precinct example, the realised ordering gives the veteran yes rate [0.45, 0.875] (N7).
  * Suppose an LCE fit on the full data gives an actual-rate interval of [0.30, 0.60]. The combined answer is **[0.45, 0.60]**, valid for the actual rate at the LCE interval's confidence level. (Illustrative LCE numbers.)
  * If the LCE interval were [0.20, 0.40], the two would not overlap. That rejects the conjunction of "LCE plus its rows" and the ordering: a false ordering can also cause it (PHASE6a R27(ii)). Only with margins alone does a rejection fall on the model.
* **M-CONTRAST + an M-SIGN sign (L, model layer).** Under pooled M0 both are linear rows in k, so they share one LP and the result is sharp for the model target (R3).
  * With **pairs** (chosen comparisons), the contrast caps the gap at each comparison's midpoint, and the sign constrains a weighted average across all precincts. Neither implies the other, so both can narrow.
  * With **hull** (every comparison), the between sign is implied (N18), and `eiRedundancy()` shows zero gain.
  * In the count layer, the contrast's bins of actual voters combine with a *realised* sign in one LP. With a *model* sign they combine by ∩.
* **B1(M) + M-SIGN (L, model layer, open).** Example: Wisconsin 1865, where B1(M = 1) puts the gap at [0.049, 0.297] (plug-in). A positive between sign on the same mean then adds nothing at the lower end, because the gap is already positive. It could still trim the upper end. The redundancy diagnostic reports which. Sharpness for the augmented model is open (R3 for B0; R8). Under county rungs, footnote 3 applies.
* **DiD + CR on changes (L, count layer, exact).** Realised CR on changes (2d Δ4) confines the descriptive DiD to the interval between the neighbourhood model and Goodman, both computed on changes. Its sign is the sign of Goodman on changes.
  * Example: the arithmetic range for the DiD is [−0.30, 0.60]. Goodman on changes is 0.25 and the neighbourhood model on changes is 0.10. With CR on changes the DiD is in **[0.10, 0.25]**. (Illustrative numbers.)
  * The causal layer then shifts this by the trend dial (↑). A model-layer FD interval for the same DiD combines with it by ∩, with the actual-DiD allowance (2d §11).

**What this means for the API.**
* Every restriction carries `object = "actual" | "model"`.
* `eiBounds()` builds one LP per layer, intersects the layers' actual-rate answers, and prints the cross-layer test.
* Cells marked × are refused with a message.
* Cells marked R run, and the redundancy line explains them.
* Cells marked O print "experimental: sharpness/inference not proved for this combination".

**Proof items this adds.**
* **R27**: coverage of the cross-layer intersection, and the level of the disjointness test. This is a one-line union-bound argument, but it needs writing out.
* R3 extended to the B0 model (LCE + Γz), and to LCE + η(Z).
* N18's implication written as a lemma.

"""
---

## Part 4. NEW_AVENUES.md, avenue by avenue

*This part works through `NEW_AVENUES.md` section by section. Each entry has four parts:*
1. *what the avenue is, in plain language;*
2. *a small worked example, with and without it;*
3. *my verdict: **build now** (in the numbered order), **build later**, **research only**, or **drop**. The verdict says whose data it serves: **yours** (Wisconsin/Iowa aggregates without poll books, NZ with electorate tables) or **other users'**;*
4. *how it fits the plan.*

*Numbers marked [NA] come from NEW_AVENUES' own scripts in `new_avenues/`.*
* *Re-run here and **reproduced byte for byte**:*
  * `check_structural_claims`, `check_review_items` (parts (i)–(iii); part (iv), the Wisconsin covariance, needs data not in the repository), `check_car_sensitivity`, `check_single_index_dim`, `check_nesting`, `check_nesting_lce`;
  * `sim_second_moments`, `sim_two_scale_context`;
  * `sim_local_near_constancy` and `sim_size_grid` (log in `new_avenues_rerun_log.txt`).
* *Taken on trust (the input data are not in the repository):*
  * the Wisconsin/Iowa figures (`apply_B_suffrage`, `apply_kappa_suffrage`, and the §13.5 covariance figures);
  * the NZ figures (`nz_diagnostics`, which needs an R export from `ei.Datasets`).
* *M-15 (the cluster-based figures).*
  * The NZ variance slope (§1.6) has clustered SEs, and so does the Wisconsin 1865 correlation interval (§13.5).
  * Neither verdict below depends on them. The NZ slope's *point estimates* (0.69–0.81) do not depend on how the SE was computed, and the §13.5 verdict rests on the 1857 point estimate not being PSD.
  * Both still need re-deriving under decision D3 before anything quotes their intervals. I could not re-derive them here: the container has no R and the unit-level data are not committed.

### 4.1 Avenue A: unit size and the second moments of returns (§1, D1 feasibility check)

**In plain language.** Small polling places bounce around more than big ones. *How much* more they bounce depends on how different the groups inside them are. If voters within a group vote independently (the coin-flip assumption every likelihood-based EI model already makes), the size–noise relationship reveals how far apart the groups' rates are. That is exactly what the regression cannot see.
* **A1** uses independence plus "rates do not depend on size".
* **A2** adds "shocks do not push the groups in opposite directions".
* **A0** (no assumption) is already implied by the sharp set and never binds.

**Worked example** [NA, `sim_second_moments`, `sim_size_grid`]. A 2×2 design with true loyalty 0.840.
* *Without A*, the current package's set is [0.521, 0.941].
* *With A1* it is [0.521, 0.841].
* *With A1, A2 and ordering* it is [0.839, 0.841], almost a point.

That is the population picture. In samples:
* the plug-in version covers the truth only 4–17% of the time;
* the valid version gains nothing at 64 or 500 units;
* with **25,000 units of median size 300**, the A1+A2 width is 0.044 against 0.363 without it, with coverage 0.93 (`sim_size_grid`).

NZ has 31,359 booths with median size 262, and its variance slope on the binomial term is 0.69–0.81 [NA, D1; point estimates; SEs under M-15]. That is consistent with coin-flip voting plus heterogeneity.

**Verdict: research only** (not in the build order). It applies to **NZ-type data**, but only if the national file is pooled. That requires a cross-electorate model (a common mean structure), which is itself a strong assumption. A2, the source of most of the gain, is the least defensible piece (§1.4).
* *Depends on*:
  * moment-inequality inference for conditional variance functions (Andrews–Shi type);
  * the R×C version of A2 (a conjecture);
  * a decision on how to pool NZ across electorates.
* A0: **drop** (implied by Θ_X, never binds).
* If Avenue A research proceeds, the covariance route of §13.5 (estimating the cross-group correlation from Var(T | x) under homoskedastic deviations) is the candidate replacement for A2's sign assumption. It is otherwise dropped from the plan.

**Fit.**
* LCE layer: A1/A2 are SOCP or LMI rows, or LP rows with ordering in 2×2.
* No change to the order. A good self-contained research project on NZ once the numbered steps are done.
* A0's relation to Step 8: it is a weaker form of the same full-law information.

### 4.2 Avenue B: context measured at a second scale, B1(M) (§2)

**In plain language.** A town's own veteran share has effects the data cannot assign to veterans or to non-veterans. The *county's* veteran share is different: its effect on each group **can** be read off the data. B1(M) assumes the town's own make-up moves each group's rate at most *M* times as much as the county's make-up does. *M* is a dial, benchmarked by a measured contextual effect of the same kind as the threat. §14.5 argues that this makes it better than benchmarking against ordinary covariates.

**Worked examples.**
* *Simulation* [NA, `sim_two_scale_context`]. Truth 0.727.
  * The bound without B is [0.495, 0.805].
  * With B1(M=1) it is [0.693, 0.805]; with M=0.5, [0.722, 0.778].
  * When the truth violates the assumption (own effects 1.5× the area's), M ≤ 1 excludes the truth, as it should.
* *Your data* [NA, `apply_B_suffrage`; trusted; plug-in, no confidence intervals]. Wisconsin 1865, trimmed units, veteran–non-veteran gap in the yes rate:

  | Rung | Gap |
  |---|---|
  | Duncan–Davis | [−0.273, 0.594] |
  | Pooled LCE (M0) | [−0.066, 0.452] |
  | B1(M=1), yes column | [0.049, 0.297] (the sign of the gap is settled) |
  | B1(M=0.5), yes column | [0.111, 0.235] |

  Two warnings:
  * Imposed on **all three columns** at once, B1 is empty unless M ≥ 3.85 (1857) or 6.65 (1865).
  * With all units used as constraint points, even B0 is empty in 1865, and in Iowa 1868 M0 itself is empty.

  Plug-in emptiness is not a test (the package's own rule), so these need the emptiness test before they are read as evidence.

**Verdict: build now, as Step 4** (with κ annotation, §4.8(a); after E4 on NZ, which calibrates it). It applies to **Wisconsin/Iowa-type data** (aggregates with a second-scale composition): the county share is the second scale, and the effects are identified. It is NEW_AVENUES' first choice for near-term use, and it gives the DiD its dial. M-CTX2-Δ replaces FD3's free δ with a measured yardstick (PACKAGE_DESIGN §7).
* *Depends on*:
  * R8, the B1 rank condition (NA Prop. B1, a sketch). Proving it is the first task of the step;
  * the parametrisation-invariant M (§2.3; the Wisconsin numbers above use the entry-wise form);
  * the existing inference machinery: few rows, so CLR and projection apply under pooled M0 by R3/R10. Grouped versions follow R16.

**Fit.**
* LCE layer only.
* It changes the order: it enters at Step 4, ahead of CAR(Z)-δ (Step 6, which shares its model) and M-EXCEPT.
* Per the NEW_AVENUES author (§4.13): do it **per column, in the range form** of §2.3, mainly for the yes column or the yes/no gap. The entry-wise all-columns figures (M ≥ 3.85/6.65) partly reflect tiny county effects in the no/abstain columns. Taken at face value, they say county share is a weak benchmark for those columns, which is worth saying in the paper.
* Report it per column *and* jointly, and show M_min (the smallest M the data allow) as a diagnostic beside it.

### 4.3 Avenue C: several contests sharing the same rows (§3; relation to M-CVR)

**In plain language.** Two contests on the same ballot (governor and suffrage referendum) are two EI problems with the same groups.
* **C1** (proved): knowing how the two contests' votes combine, without any assumption linking them, does not narrow either contest's group rates.
* **C2** treats contests like districts: context works the same way in each, up to level and scale. Contests with lopsided results then pin the pattern for the others.
* **C3/C4** (a single index or a tensor model) are stronger.
* **M-CVR** (general-user track G2) is different from all of these. It uses *observed* ballot patterns to bound *differences* across contests and needs no linking assumption. C1 and M-CVR are consistent: patterns pin cross-contest differences, not single-contest rates.

**Worked example** [NA, `check_structural_claims`, reproduced]. In 300 random tables, adding joint contest margins changed the DD bound for a contest-1 cell by at most 2e−16. That is C1: no gain without a link. With C2, a lopsided referendum in which some group votes almost unanimously pins that group's contextual pattern. Under the linking assumption it then also pins the pattern for the closer contest. No numbers exist in the repository for C2.

**Verdict.**
* **C1**: proved, nothing to build.
* **C2**: **Step 10** (added by your decision; §3.7d). It reuses the M1 intersection machinery, with contests in place of districts.

  For your data, Wisconsin has the governor race and the suffrage referendum on the same men. But context working "the same up to scale" for a partisan race and a suffrage question is a strong claim. 2d §9 deliberately uses party only through the count-layer results P2/P3.
* **C3/C4**: research only / drop. C4's local independence is implausible for partisan contests (§3.2).
* *Nesting with slack* (§12.3(d), "suffrage yes ⇒ Republican, except ε"). This is the count-layer form of a contest link that you can use, so fold it into **Step 2's party layer**. On a single outcome it changes nothing (verified, `check_nesting*`). On cross-outcome contrasts it narrows the width to a median of 0.39 [NA].

**Fit.**
* C2: LCE layer.
* Nesting: count layer, in Step 2.
* M-CVR: count layer (G2).
* No change to the order beyond adding nesting to Step 2.

### 4.4 Avenue D: the polynomial ladder with LP certificates (§4)

**In plain language.** The straight-line (LCE) assumption can be relaxed to gentle curves: quadratic, cubic. This *widens* the bounds on purpose. Its use is as a robustness check: "how much of my answer comes from assuming straight lines?" The repository's contribution is an exact pair of LPs, an outer bound (validity checked on a grid) and an inner one (Bernstein certificates), so the true curved bound is always bracketed.

**Worked example** (new, N16). Five kinds of precinct, true rates linear, target the group-1 rate (truth 0.632):

| Rung | Bound |
|---|---|
| Straight lines | [0.436, 0.779] |
| Quadratic | [0.419, 0.796] |
| Cubic | [0.359, 0.820] |
| No model (Duncan–Davis) | [0.324, 0.834] |

The ladder shows how quickly the LCE answer gives way as curvature is allowed.

**Verdict: build, as Step 9** (after Step 8). It applies to **any data**: it is the natural robustness display for every LCE result in the NZ and Wisconsin papers.
* *Depends on*: the outer/inner pair (designed, NA §4.2). The degree must be prespecified (2c §6 counterexample).
* L1′ (quadratic contextual effects, on the repository's list) is its J=2 case. Build them together.
* It does not tighten anything, so it does not compete with the steps above for user value. It does compete for credibility, so it should not slip far.

**Fit.** LCE layer. Step 9 in the revised order (§3.7b).

### 4.5 Avenue E: known tables at another level (§5, §13.1–13.3)

**In plain language.**
* **E1/E4**: a published cross-tab for a larger area pins part of the hidden pattern if the pattern is shared. Two tables for differently composed areas can pin all of it and test the sharing.
* **E2**: transfer of the pattern across elections.
* **E3** (proved): coarser margins of the same classification add nothing.
* **E5**: richer margins of the same people (age-stratified, crossing boundaries) give gains with no assumption.
* **E6**: validation microdata (exit polls, linked records), with a selection-odds dial Γ.

**Worked examples.**
* E4 (N15): one electorate's table pins a neighbour's rate from [0.529, 0.751] to 0.656 if the pattern is shared, and to [0.650, 0.662] if it may differ a little.
* Rank [NA, `check_review_items`, reproduced]: in 3×3, one table pins 4 of 6 unknown directions and two tables pin all 6, with 2 testable overidentifying restrictions.
* E5 [NA, reproduced]: age-stratified margins cut the DD width for an aggregate rate to 0.68 of its width. Crossing partitions moved the bound from [0.098, 0.939] (intersection of the two DD bounds) to [0.118, 0.930] (joint LP), a modest gain.
* E2 [NA, D2; trusted]: NZ's table-implied K correlates only 0.03–0.51 across consecutive elections, so transfer across elections needs a large δ.

**Verdict.**
* **E4: build last, as Step 11** (yours, NZ; moved last by your decision, since most applications will lack such tables).
* **E1**: part of E4.
* **E2**: drop as an identifying rung, and keep as a dial within E4. The NZ evidence is weak.
* **E3**: proved, nothing to build.
* **E5**: other users'. Its rows ship inside Step 1.
* **E6**: other users' (no validation data for your projects). Build later as part of the general-user track, alongside G1 (it is the same "recovered information" logic plus the Γ dial).

**Fit.** E4: LCE layer with the actual-vs-model allowance. E5/E6: count layer. No change to the order.

### 4.6 Avenue F: computing the sharp set (§6; relation to the M-SHARP oracle study)

**In plain language.** Use the whole spread of returns at each composition, not only the average. (Same as Step 8.)

**Worked example.** N5 (Step 8): a group's mean at a 50/50 composition goes from [0, 1] to [0.4, 0.6] when half the precincts vote 10% and half 90%. With unimodal heterogeneity the gain is a median of 0–1%; with bimodal (polarised) precinct rates it is 8–37% (my oracle study, extended after consultation).

**Verdict: build later, as a diagnostic (Step 8).** It applies to **NZ-type data** as the place to test the size channel. Two pieces of §6 are now settled or redirected:
* The "finite core-determining class" conjecture is **proved** for 2×C and R×2 by the review's subset identity (my R24). The class is the proper subsets of outcomes (or of groups).
* The "kink-mass" diagnostic is the right trigger for running the oracle at all.

**Fit.** LCE layer. Unchanged position.

### 4.7 The lower-ranked avenues (§7)

| Avenue | Plain language | Small example | Verdict, and whose data |
|---|---|---|---|
| **G, label symmetry** | Assume the model treats, say, Labour and National symmetrically. This roughly halves the unknown directions | In a symmetric yes/no referendum the free directions halve. In NZ's asymmetric party system the assumption is hard to defend | **Drop** as an identifying rung. At most a labelled sensitivity check for symmetric referenda (other users) |
| **H, deconvolution** | Recover the law of precinct shocks from the spread of returns | Proved: the shock law does not involve K. It trims only the support (m ≥ −min e), a small effect | **Drop** |
| **I, spatial structure** | Neighbouring polls as replicates for Avenue A's variance estimates | Neighbours with similar composition give local estimates of Var(T \| x, n) without a global polynomial | **Research only**, inside Avenue A (yours, NZ, if A proceeds). Neighbour-only *similarity* gives no gain (§12.3(f), `sim_local_near_constancy`) |
| **J, entropy / OT / Choquet** | Pick the "most random" table inside the bounds | Maximum entropy at fixed margins *is* "no within-precinct association": it selects a point, it does not narrow | **Drop** as identification. The honest form, distance to an explicit reference table (KL ball, §12.3(e)), is in §4.8(e). Choquet capacities are the computational form of R24 (Step 8) |

### 4.8 New rungs, modifiers and combinations (§12.3–12.5)

*(a) Variance-capped contextual effects, κ, and κ-annotated endpoints.*
* **Plain language.** Rates may vary with composition in any pattern, but not by much. The reporting device needs no commitment: label each end of a bound with "this end needs contextual variation of at least κ".
* **Example** [NA, `apply_kappa_suffrage`; trusted; plug-in]. Wisconsin 1865 veteran–non-veteran gap:

  | κ | Gap |
  |---|---|
  | none (M0) | [−0.065, 0.452] |
  | 0.25 | [0.019, 0.452] |
  | 0.1 | [0.140, 0.396] |

  The data floor is κ = 0.015.
* **Verdict: build now, inside Step 4**, with B1(M). The annotation costs little. The cap is LP rows in 2×2 and an ellipsoid in R×C. Yours (Wisconsin/Iowa) and others'. *Depends on*: R8 (validity and monotonicity); the SOCP form in R×C.

*(b) Single-index contextual effects.*
* **Plain language.** Context matters through one dimension of a place's make-up, for example a left–right index of the electorate.
* **Example** [NA, `check_single_index_dim`, reproduced]. In 3×3 the unknown directions fall from 6 to 2.
* **Verdict: research only for now.** It could be large for **3×3 tables such as NZ's**, and E4's leave-one-out is the natural test. But R14 is a conjecture and the computation is nonconvex (an LP per direction on a grid). Revisit after Step 11.

*(c) Odds-ratio dial.*
* **Plain language.** In each precinct, the group-by-vote odds ratio lies in a range. "At least 1" is the ordering belief of Step 1. An upper limit adds "but not extremely polarised".
* **Example** (new, N17). Precinct P2's veteran yes rate:
  * DD: [0.167, 0.833];
  * odds ratio ≥ 1: [0.500, 0.833];
  * odds ratio between 1 and 4: [0.500, 0.633].
* **Verdict: build now, inside Step 1**, as a unit-level restriction type, for 2×2 and pairwise 2×C. It serves yours and others' data. Label it as untestable when the dial contains 1.

*(d) Outcome nesting with slack.* Count layer, in **Step 2's** party work (§4.3).

*(e) Distance to an explicit reference (KL ball).*
* **Plain language.** "The truth is within a stated distance of a defensible reference table": a survey, a prior election, or another jurisdiction's published table.
* **Example.** No numbers in the repository. With the reference set to the independence table it would pull bounds toward "no association", which is why §12.3(e) forbids that reference.
* **Verdict: build later**, as a dial in M-DATA, when a user has a reference. Your only candidate is a prior NZ election, and E2's weak persistence argues against it. Convex, not LP.

*(f) All-pairs similarity.* **Drop** at the unit level. *Not the same as Step 5:* Step 5 bounds differences between *named neighbouring districts'* contextual patterns under LCE, a district-level restriction on K deviations, not unit-level similarity of rates.
* At best it reaches 0.5–0.6 of DD width, with large units and low heterogeneity [NA, `sim_local_near_constancy`].
* The data floor is a refutation threshold, not a calibration.
* The neighbour-only version gains nothing.

*(g) Data-informed dispersion ceiling.* **Research only.**
* A variance ceiling gives no hard unit bounds.
* In Wisconsin the ceilings are weak for veterans (SD ceiling 0.42/0.37) [NA; trusted].

*(h) Finite-sample bands.* Already Theorem P3.8. **Nothing new to build.**

*Combinations (§12.4) that change the plan.*
* **"No offsetting" + ordering ⇒ the gap ≤ the observed slope.** This is the local form of M-CONTRAST (Step 3). M-CONTRAST's two signs (the mix term and the rate term are both ≥ 0) cap the gap by the observed change between two kinds of place. That is the finite, derivative-free version of §12.5's first-ranked item, so Step 3 delivers it.
* **M1 + B1(M) with z varying within groups.** This is Step 4 under grouped rungs. Wisconsin's county share varies within counties only through the unit's own share, so B1 must use the pooled model or groups coarser than counties.
* **Validation tables calibrate the dials.** Step 11 (E4 on NZ) should output the empirical distribution of κ, δ and M that the NZ truths require. This turns E4 into a calibration source for Step 4's dials.
* **Coin-flip (A) + accounting** and the **dispersion ceiling + sup-norm cap**: research, with Avenue A.

*§12.5 list.* Items:
1. No-offsetting ceiling → Step 3.
2. κ and κ annotation → Step 4.
3. Single index → research.
4. Dispersion ceiling → research.
5. Odds-ratio dial → Step 1.

### 4.9 CAR as a rung with a dial (§13.4) and §14

**In plain language.** CAR says that, among places alike on covariates (wealth, nativity, denomination), the group mix tells you nothing about how each group behaves: no contextual effects once the covariates are fixed. CAR-δ relaxes this by letting rates still depend on composition by up to δ. §14 adds three findings:
* no test can detect a violation of the kind EI fears most;
* the Cinelli–Hazlett-style sensitivity bounds implemented in `seine` are systematically optimistic about contextual effects when benchmarked on ordinary covariates;
* B1(M) is the benchmark of the right kind.

**Worked examples** [NA, `check_car_sensitivity`, reproduced].
* A pure contextual violation leaves the regression *exactly* unchanged (max deviation 1e−16) but biases the rate by 0.080.
* A sensitivity bound benchmarked on an ordinary covariate gives 0.006–0.032 for that same bias of 0.080 (understated by up to 12×).
* Overlap: with compositions confined to [0.3, 0.7], a δ = 0.05 relaxation leaves a half-width of 2.5 × 0.05 = 0.125 at the vertex. With [0.4, 0.6] the factor is 5.

**Verdict: build now, as Step 6** (revised after your decision D5: covariates are commonly available in applications).
* **Your data and others'.** It needs township-level covariates with composition varying within covariate cells (positivity); which covariates is decision D10. With county-level Z only, CAR collapses into Avenue B (§14.7), so Step 6 shares Step 4's model and is reported beside it.
* *Depends on*: R8 (validity, monotonicity), the overlap factor beyond 2×2 (a conjecture), and the SOCP form of the κ/r² dials.
* **Cinelli–Hazlett / Chernozhukov et al. contours**: **do not build**. Point users to `seine`, and print §14.4's caution when they do.
* Step 1's strata option already gives the assumption-free version for users with stratified tables.

**Fit.** LCE layer (CAR-δ is an LP with a dial, centred on Goodman-within-Z). It is not in the numbered order.

### 4.10 Covariance of unit-level rates (§13.5, §14.6)

**In plain language.** If how much group rates vary around their means does not depend on composition, the variance of returns reveals how the groups' rates co-vary across places, including their correlation. It does *not* reveal the hidden pattern.

**Worked example** [NA; trusted; the interval is under M-15].
* Wisconsin 1857: the estimated covariance matrix is not even valid (implied "correlation" 3.35).
* Wisconsin 1865: the correlation is −0.16, with an interval of [−0.49, 0.70] (cluster bootstrap, pre-D3).

So the data cannot sign the correlation. Without the module, A2 must *assume* it is non-negative. With the module, it would be estimated, but here it is not estimable.

**Verdict: dropped from the plan** (your decision, 21 Nov). It does not move aggregate bounds (§13.7 D), and the Wisconsin illustration was negative (1857 not PSD) or unsigned (1865). Its only remaining role, replacing A2's sign assumption, is recorded as one line inside Avenue A research (§4.1).

### 4.11 Excluded composition shifters (§13.6)

**In plain language.** Something that changes *who lives where* without changing *how anyone votes* would identify the group rates, like an instrument. In EI, almost anything that changes the mix of people around a voter also changes their context, so the exclusion is hard to defend. The clean case is re-drawing *the same people* into different reporting units. That is crossing partitions (E5(b)), which needs no exclusion.

**Worked example** [NA, `check_review_items`, reproduced]. Crossing partitions moved a bound from [0.098, 0.939] to [0.118, 0.930] with homogeneous atoms, and barely moved it otherwise. A true shifter has no example in the repository.

**Verdict: drop** as a module. The clean case lives in E5 (other users', count layer). Boundary changes over time add turnover (A1).

### 4.12 Reconciliation with NEW_AVENUES' own ranking

| NEW_AVENUES ranking (§9; revised in §13.8 and §14.7) | Their rank | Mine | Why we differ |
|---|---|---|---|
| E4 several validation tables (+δ) | 1 | Step 11 | Same direction. I put the count layer and the exact DiD first because your Wisconsin/Iowa work needs no tables, and the E4 rank result needs a proof (R21) plus the allowance |
| E5 stratified / crossing / multi-way margins | 2 | Inside Step 1 (other users) | You have no such tables. The rows are nearly free once the count layer exists, so they ship with it rather than as a step |
| B second-scale context, B1(M) | 3 (1st in §9) | **Step 4** | Agree that it is the most practical new rung for data with a second-scale composition (county share in the Wisconsin/Iowa illustration). The plug-in emptiness on all columns makes inference (CLR or projection plus the emptiness test) part of the step, not an afterthought |
| CAR(Z)-δ | 4 | **Step 6** (after your decision D5) | Shares B1(M)'s model and comes right after it. With county-level Z it reduces to B. B1's benchmark is like-for-like for contextual effects; covariate benchmarks are optimistic (§14.4) |
| E6 validation microdata with Γ | 5 | General-user track | No validation data for your projects |
| §13.5 covariance → dispersion ceilings | 6 | Dropped (one line inside Avenue A research) | Wisconsin estimate not PSD (1857) or unsigned (1865); local only |
| No-offsetting ceiling; single index; κ (§12) | 7 | No-offsetting → Step 3 (M-CONTRAST); κ → Step 4; single index → research | M-CONTRAST is the finite, derivative-free form of the no-offsetting ceiling. κ annotation is cheap and pairs naturally with B1(M). The single index rests on a conjecture (R14) and is nonconvex |
| A coin-flip rungs | 8 | Research | Agree: strong identification, weak estimability. `sim_size_grid` suggests NZ's national file is the one place it could work, so I would name that as the research project |
| Excluded shifters | 9 | Drop (clean case in E5) | Agree |
| D Bernstein ladder | 5th (§9) | Step 9 (robustness) | Agree that it is infrastructure. I give it a firmer place because it is the credibility check for every LCE number in your papers |
| F sharp set | 6th (§9), "pinned" | Step 8, diagnostic only | The oracle study quantifies §9's "power 2": median gain 0–1% with unimodal precinct heterogeneity, 8–37% with bimodal (polarised) precincts; the kink-mass diagnostic decides |
| G–J | low | Drop (G as a sensitivity check; I's replicates inside A) | Agree |
| §12.3(c) odds-ratio dial | §12.5 item 5 | **Step 1** | It is a count-layer row, and ordering is its special case, so it costs nothing extra |
| §12.3(d) nesting | not ranked | **Step 2** (party layer) | It gives gains for the cross-outcome contrasts your party analysis uses |

**Where NEW_AVENUES is now superseded or needs a note.**
* §6's finite core-determining class is proved for 2×C/R×2 (R24).
* §12.5 item 1 is delivered by M-CONTRAST in finite form.
* The M-15 intervals still need re-deriving.
* The Wisconsin B1 results should be recomputed with the package's inference once Step 4 exists, before any of them is quoted.

### 4.13 Consultation with the authors of NEW_AVENUES and PACKAGE_DESIGN

*At the owner's suggestion I briefed both sessions on everything since their last work and asked for their response to this plan. Session `session_01S2ZdrGzqoc3nxJzqVGzaWE` wrote NEW_AVENUES (last active 2 Nov). Session `session_01S2ZdrLS99i21KeJ5ubyJq1` wrote PACKAGE_DESIGN and the Batch F brief (last active 18 Nov; the owner approved its reply). Their answers are recorded as information. Verdicts changed only where their arguments persuaded me, and each change is marked.*

**NEW_AVENUES author.**

| Point | Their view | Agree? | Effect on the plan |
|---|---|---|---|
| Avenue A research only; A0 dropped | Accepts both. NZ's binomial slope (0.69–0.81, not 1) already argues against the pure coin-flip rung. Pooling 31k booths would add size invariance and a common $A(x)$, $B(x)$ across electorates, which are assumptions, not free data. Keep the binomial-slope check as a cheap diagnostic | Agree | **Changed:** the "does variance scale like 1/n?" check joins Step 7's diagnostics as a diagnostic, not a rung |
| The oracle result | Fits a pattern it found throughout: local or full-law information adds little to aggregate bounds unless something imposed (a ceiling or dial) does the work | Agree for unimodal heterogeneity. The extended oracle (below) shows the exception: bimodal precincts | none beyond Step 7's revision |
| B1(M)+κ for Wisconsin/Iowa; the emptiness reading | Reading confirmed: all-columns plug-in empty unless M ≥ 3.85/6.65; yes column M_min 0.05/0.3; B0 empty in 1865 with all units but not with trimmed points. Its interpretation: (a) plug-in inward bias and small or extreme units, so report M_min with an interval; (b) the entry-wise M blows up where no/abstain county effects are tiny, so use the range form per column; (c) at face value, county share is a weak benchmark for no/abstain. Iowa 1868 is empty at M0, a data or base problem | Agree | **Changed:** Step 4 is per column, in range form, mainly for yes or the yes/no gap, with M_min and its interval as a diagnostic |
| CAR-δ and D7 | Later. With county Z it is B1(M=0). Its offered §14.8 κ benchmark is B1(M) in variance-share form, so fold it into the B1 step. Ordinary-covariate benchmarks are optimistic (§14.4) | Agree | **Changed:** the §14.8 κ benchmark folds into Step 4 |
| What to move earlier | **E4 on NZ, before or alongside B1(M).** It is the only place the truth is known, and its κ/δ/M output is the one empirical calibration for reading the Wisconsin/Iowa dials | **Persuaded** at the time | Moved E4 to Step 4 on 20 Nov. **Overruled by your decision (21 Nov):** E4 is now last (Step 11 after later additions), because most applications will lack the tables. Its calibration role is kept as a reference distribution printed beside Step 4's dials once Step 9 exists |
| Drop list; general-user track for E5/E6; covariance research only | Agrees with all | Agree | none |
| M-15 flags in §1.6 and §13.5 | Offers to add them to NEW_AVENUES | — | Not mine to authorise: **owner decision D8** |

**PACKAGE_DESIGN author.**

| Point | Their view | Agree? | Effect on the plan |
|---|---|---|---|
| Moving M-SHARP down | Yes, with conditions. N5 shows polarised precincts can move a bound a lot. Test high-kink-mass designs (bimodal returns, many small units) and run the NZ kink-mass diagnostic before demoting it for good. Restating R9 is right. Make $(X,N)$ conditioning the default | **Persuaded** on the test, which I ran. Bimodal precinct rates give 8–37% narrowing (ratio 0.63–0.92); small precincts alone give 1–8% | **Changed:** Part 1 and Step 7 now say the gain depends on polarised precincts. The kink-mass diagnostic runs on NZ and Wisconsin as soon as Step 1 exists, and promotes the estimator if it finds bimodality. `eiSharpCheck()` conditions on $(X,N)$ by default |
| Folding M-LOC into M-CONTRAST | Yes. It sidesteps Batch B comment 3 (the unidentified local derivative). Keep the P-S/R12 references so `at = "hull"` inherits the proved sandwich. For R ≥ 3, attributing a contrast to one pair is a belief. Data-defined bins must be prespecified or reported as a family, otherwise it is belief-shopping | Agree | **Changed:** Step 3 requires prespecified bins (or a reported family) and cites P-S/R12 |
| Count layer and the combination rules | Fits. The object attribute splits "unit" into counts vs model mean. Rule: rows on counts and on $k$ never share one LP. The actual-rate answer is the intersection of the exact count set with the LCE interval plus allowance. That covers $B$ at $1-\alpha$, and disjointness is a free test of LCE (add to §4.3). Footnote f becomes this rule | Agree, and it improves my plan | **Changed:** Step 1 states the rule and the free test |
| Caution (a): one vocabulary | "Veterans more supportive" can be realised (2c M2) or model (M5). Make the user choose with an `object =` tag; never impose both under one name | Agree | **Changed:** Step 1's interface uses `sign`/`cr` with `object =`, not separate `*_realised` names |
| Caution (b): integrality | Pooled realised signs and CR break total unimodularity. Realised CR is two MILPs when integrality matters; "never refutable" must say "continuous tables"; the R2 face test carries over to the count LP | Agree (consistent with N10c); later narrowed by PHASE6a: one-vs-rest W1 rows are totals and stay integral, so only pairwise, W2 and coupled rows need MILPs | **Changed:** Step 1 notes added; `eiRedundancy()` extends to the count layer |
| Review gates | The count layer's realised rows rest on 2c M2 (Batch B comment 1 still open). Step 2 rests on Δ0/Δ4 (Batch C). Step 4 rests on R8. All carry "experimental" until resolved, and the EXTERNAL_REVIEW_PLAN gate on DiD wrappers applies to the descriptive DiD too | **Partly.** "Experimental" labels until review: yes. But I keep building Step 2's descriptive layer before Batch C returns. Its arguments are short and the count layer's brute-force tests check them directly, so the gate should apply to *removing the label*, not to building | Labels added. The disagreement on the build gate is noted here for the owner (decision D4) |
| Proof register | Add R19, R22, R24; reword R1; restate R9; update "Blocks" now that R5 and R6 are deferred | Agree | **Changed:** register rows added (Part 8) |
| Reversing §6's order | Agrees: proof readiness applies at least as well to the count layer | — | none |

**Net effect.**
* The order changed: E4 on NZ moved earlier. Your later decision moved it to last (Step 9).
* M-SHARP's verdict is now conditional on the data rather than "little gain".
* Step 1 gains an explicit combination rule, a free test of LCE, and one restriction vocabulary.
* Steps 3 and 5 gain the guardrails listed above.


---

## Part 5. Not now, and why

| Item | Why not now | What would change that |
|---|---|---|
| **M-SHARP series estimator** (NEXT_BUILD item 4, R9 as stated) | The oracle gain is median 0–1% (at most 5%) under unimodal heterogeneity, 8–37% under bimodal precinct rates (`sharp_gain_study*.py`); build the estimator only if the kink-mass diagnostic finds bimodality. Plug-in versions are biased inward. "Tightens every rung" is not a sensible release criterion | Step 8's diagnostic finds data with large kink mass, or the NZ benchmark shows gains |
| **M-LOC as its own module** | It is the "everywhere" option of M-CONTRAST under LCE (N12c). C1 (R1) is needed only for the implication to $\delta_W$ | none; build it with Step 3 |
| **2×C joint/ordinal signs (R6), R×2 spanning trees (R5)** | Evidence of user demand is thin, and R5(c) and R6 are open or only checked. These are LCE-only refinements | an application with ordered outcomes or ≥3 groups that asks for them |
| **M-DISP, Avenue A, single index, tensor contests** (see Part 4 for each; the covariance diagnostic is dropped) | A and the dispersion route need thousands of units and a pooling model; the single index rests on conjecture R14 and is nonconvex; tensor contests need implausible local independence | a decision to pool NZ nationally (A); a proof of R14 |
| **Causal DiD before the descriptive layer and Batch C** | PHASE1b and PHASE2d are unreviewed; the causal labels need explicit trend assumptions | Batch C returned |
| **Switching the default interval to projection** | **This would reverse your decision.** I do not recommend it now. 0.3.1 already prints projection beside CLR with a per-fit status, and Step 0 adds the R2-regime line | Conjecture P3.5(b3) fails in a design users meet, or a referee insists |
| **General R×C full-law estimation** | The cost is high and the gain unproven (above) | Step 8 shows gains |
| **Covariate adjustment beyond stratified margins** | Low priority by your decision. Step 1's strata option delivers the assumption-free version (stratified DD) for users who have the tables | your call |
| **Superpopulation target, clustering** | Excluded by your decisions | — |

---

## Part 6. Where I agree and disagree with the external review

**Agree:**
1. **Lead with the quantity the question asks about.** Step 1 and decision D1. Part of this is done already (`ci_actual` ∩ DD, `proj_actual` on P3.5(c1)).
2. **A count-table interface beside the LCE solver, keeping the skew-symmetric LP for LCE.** This is Step 1. I would also adopt the review's request that each module declare which variables it restricts (actual cells or model means). The `meta$object` field already exists, so the remaining work is to make combination rules check it.
3. **The subset-inequality identity is right.** Verified (N1–N3c), with a short submodularity proof (R24). Averaged cell limits suffice for C ≤ 3.
4. **Targeted validation.** General-user feature G1 (no audit data exist for your projects). I add the concavity fact (N8d), which keeps the coupled case cheap.
5. **Finite-composition contrasts instead of local derivatives.** Step 3. I go further: under LCE the contrast rows *are* the P-S objects at the midpoint (N12c), so M-LOC need not be a separate module.
6. **Integrality**: an integer check where rows break total unimodularity. N10c gives a concrete LP-feasible, integer-infeasible system.
7. **Descriptive DiD before causal**, and the "MTR" mislabel.
8. **Exception-mass frontiers** and **group misclassification** as useful later modules.
9. **Ballot patterns** as a strict gain from richer data, with C1 left intact.
10. **B1(M) where a second scale has adequate variation.** The review places it at its step 5, and so do I. It applies because Wisconsin/Iowa's county share is that second scale, and it is the DiD's yardstick. Once Step 11 (E4 on NZ) exists, its truths provide a reference for which M values are plausible.

**Disagree, or agree with qualifications:**
1. **On M-SHARP's place.** The review puts M-SHARP second, before M-LOC and after the count layer. I put it last of the numbered steps and narrow it to a diagnostic. The review's argument rests on population information. My oracle study finds that information small under unimodal precinct heterogeneity (≤5% width; median about 0), but 8–37% when precinct rates are bimodal. So the diagnostic decides. The review itself notes that finite-sample bands can erase it. I agree with the review against the repository on narrowing ("remove 'tightens every rung'") and on conditioning on size. I disagree with the repository's NEXT_BUILD item 4 on priority.
2. **On the default interval.** The review would make the coverage claim pick the default, with projection where proved. You decided to keep CLR as the default, labelled algorithmic, with projection beside it. 0.3.1 implements that transparently. I would not reverse it now. What remains of the review's point is cheap and does not change the default: the R2-regime line (Step 0). Revisit if Conjecture P3.5(b3) fails in a realistic design.
3. **On E4 before E5.** The review's build order promotes E5 and E4 together. I would ship the exact uses (E5, published totals, strata) in the count layer first. They need no allowance and no rank proof. Transfer (E4) is an assumption and should follow with its δ dial and NZ evidence. NEW_AVENUES §5.2 E2 already found weak persistence of $K$ in NZ.
4. **On "short proved chains" replacing the ladder.** Agreed for `ladder()`'s printout (Step 0). The plausibility tiers can stay in the documents as an *elicitation aid*, labelled as such, because users find them helpful.
5. **On the review's ordering of M-LOC after the validation work.** I agree it is late. I would not build it separately at all (point 5 under Agree).
6. **On the priority of targeted validation and external data.** The review ranks targeted validation "high" and builds the audit calculation into its first step. For your projects no audit, linkage, CVR or sub-table data exist, so these become a general-user track (§3.8) behind the modules that work with aggregate margins only. The review's other high-priority idea, finite comparisons, needs no data and moves up to Step 3.

**Where I disagree with the repository's own plan:**
* PACKAGE_DESIGN §6 orders LCE refinements (M-LOC, 2×C, R×2) before data modules on grounds of proof readiness. The count-layer modules are at least as proof-ready: standard LP facts plus 2c M2, which is proved. They also answer the users' first question. So that rationale no longer favours the LCE refinements.
* NEXT_BUILD item 4 says M-SHARP "tightens every rung with no new assumption". True in principle, but negligible in the designs tested.

---

## Part 7. Decisions: settled and remaining

**Settled (21 Nov):**

| Decision | Your answer | Where it is reflected |
|---|---|---|
| D1 print order | actual-rate block first, DD first; `target = "model"` restores the old order | Step 1 |
| D2 function | separate `eiCounts()`, embedded and printed by `eiBounds()` | Step 1 |
| D3 re-ordering | accepted; E4 moved last | Part 1, Step 11 |
| D4 DiD | changing populations deferred; general theory and proofs first (descriptive and causal), with gaps named; data later | Step 2 (theory track, gap register) |
| D5 covariates | potentially available, in general | Step 6 (CAR(Z)-δ) |
| D6 integers | "auto"; but matched veteran counts are non-integer and estimated | Step 1 (auto-off for non-integer margins) and Step 1b (margin sets, R25/R26) |
| D7 reviews | two separate reviews per module (proofs, code) | review table, Part 3 |
| D8 M-15 flags | the NEW_AVENUES author adds them | — |

**Answered (21 Nov, second round):**
* **D9 (match calibration) and D10 (covariate set): deferred** to your research project. In the plan, the calibration slack and the covariate set are user inputs of general methods (Steps 1b and 5).
* **D11 (DiD scope)**: T ≥ 2 elections with pre-periods. All three causal layers (trend dial, relative magnitudes, MTR/ND/PS), with the set-valued pre-trend theory and the degeneration result in Step 2.
* **D12**: CAR(Z) before M-EXCEPT. Confirmed: Step 6, then Step 7.

**Added (21 Nov, third round):**
* C2, contests as districts: Step 10, just before E4 (now Step 11, last).
* Bounded differences between named neighbours: an option of M1(δ), Step 5, with neighbours defined by an adjacency matrix or a distance matrix (Lipschitz default).

**Final choices (21 Nov, fourth round):**
* The covariance diagnostic is **dropped** from the plan. One line inside Avenue A (§4.1) records its only research role.
* Covariate-benchmarked κ is kept as written in Step 6.

**No decisions remain open.**

## Part 8. What I verified, and what I take on trust

**Verified numerically here** (`verify_next_steps.py`, 34/34 PASS):

| Label | Claim | Result |
|---|---|---|
| N1, N1b | Subset identity 2×C equals the Aumann expectation (lifted LP) | max diff 3.7e-16 over 800 comparisons |
| N2 | The same for R×2 | 4.4e-16 over 800 comparisons |
| N3a–c | Cell limits exact for C=3; strictly looser in 103/150 (C=4) and 126/150 (C=5) random laws; the review's 2×4 example | confirmed |
| N4 | Size-conditioning example (1/2 vs [1/4, 3/4]) | confirmed |
| N5–N7, N9 | Worked examples (DD, realised ordering, stratified margins) | confirmed |
| N8a–d | Audit gain equals $d_i$ (separable); not separable with a total; remaining width concave in the revealed count | 200 designs; 1,640 concavity checks, 0 violations |
| N10a–c | Integer endpoints for margins only; a fractional endpoint with an ordering row; LP-feasible but integer-infeasible with coupled rows | confirmed |
| N11 | Ballot-pattern example (difference exactly 0; ≤ 0.08) | confirmed |
| N12a–c | Kitagawa identity; contrast example [0, 0.20]; under LCE the secant terms equal the P-S objects at the midpoint | max diff 2.5e-16 |
| N13, N14 | Exception-mass and misclassification examples | confirmed |
| N16 | Polynomial-ladder example: straight lines [0.436, 0.779] ⊂ quadratic [0.419, 0.796] ⊂ cubic [0.359, 0.820] ⊂ DD [0.324, 0.834] | confirmed |
| N17 | Odds-ratio dial example: OR ≥ 1 gives [0.500, 0.833]; 1 ≤ OR ≤ 4 gives [0.500, 0.633] | confirmed |
| N18 | Hull contrasts (O_rs ≥ 0 at every support point) imply the model between sign δ_B(r,s) ≥ 0 | 0 violations |
| N19 | Berkson correction for matched composition: E[T \| q] = q̄'Qq̄ + Var(u)·(Q_vv − 2Q_vn + Q_nn) | max error 3e−5 against simulation |
| N20 | Estimated margin: expected count 60 treated as exact gives [0.167, 0.833]; true count in its 95% Poisson-binomial set [56, 64] gives [0.107, 0.893] | confirmed |
| N23 | Neighbour-bounded deviations: A in [−0.10, 0.10] (M2 and centre form δ = .05), [−0.05, 0.05] (adjacency δ = .03; Lipschitz λ = .03), [−0.10, 0.10] if A is isolated | confirmed |
| N24 | Contests as districts: Q-relation holds exactly; group-1 rate [0.498, 0.840] → [0.539, 0.798] (truth 0.715) | confirmed |
| N22 | Covariate-benchmarked κ with a known contextual violation: covariates tied to composition need m* ≈ 0; covariates explaining the outcome but not composition need m* ≈ 2,800; covariates absorbing most composition (R²(X~Z) = .92) give κ_true = .006 but a group-rate error of −.123 (overlap factor 5.0) | confirmed |
| N21 | Set-valued pre-trends: pre-period DiD in [−1.23, 1.27]; relative-magnitude ATT (M̄ = 1) [−2.35, 2.65] jointly vs [−2.35, 2.69] endpoint-wise (no information); trend dial M = .10 gives [−1.18, 1.52]; with pre-period totals pinned to ±2 votes, ATT [−0.60, 0.90] | confirmed (grid over the joint set) |
| N15 | E4 transfer example (B's rate [0.529, 0.751] → 0.656; ±2 points on A → [0.650, 0.662]) | confirmed |

**Population studies** (`sharp_gain_study.py`, `sharp_gain_study_2xC.py`). Pooled LCE, oracle law at replicated composition types, Beta or Dirichlet heterogeneity.
* Unimodal heterogeneity, large precincts: the full law narrows the model-implied group rates by little (median 0–1%, at most 5%).
* High-kink-mass block (added after consultation): bimodal precinct rates (φ = 0.5, 1) give width ratios 0.63–0.92. Small precincts (n = 10–50) with unimodal rates give 0.92–0.99.

*Limits*:
* only 2-row designs for the high-kink block;
* the size channel is tested only through binomial noise, not through size-dependent laws (the N4 mechanism);
* linear truths.

Bimodal precinct rates are where the gains appear; the size channel (N4) is a second route. Step 8's diagnostic should look for both in NZ and Wisconsin.

**NEW_AVENUES re-runs** (`new_avenues_rerun_log.txt`). Ten of its scripts were re-run under `python3 -I`: `check_structural_claims`, `check_review_items`, `check_car_sensitivity`, `check_single_index_dim`, `check_nesting`, `check_nesting_lce`, `sim_second_moments`, `sim_two_scale_context`, `sim_local_near_constancy`, `sim_size_grid`. All reproduce their committed outputs byte for byte. The one exception is `check_review_items`: parts (i)–(iii) are identical, and part (iv) (the Wisconsin covariance) is skipped because its data are not in the repository. That script also rewrites its own output file when run. My re-run overwrote it and dropped part (iv); I restored the committed file. **Not re-run** (input data not in the repository; taken on trust): `apply_B_suffrage`, `apply_kappa_suffrage`, the §13.5 Wisconsin covariance figures (`check_review_items` part (iv)), and `nz_diagnostics`. **M-15** (clustered SEs for the NZ slope; cluster-bootstrap interval for Wisconsin 1865): not re-derived (no R, no unit-level data here); no verdict in Part 4 depends on them.

**Proved here (sketch, for review):** R24, the subset identity, via submodular base polytopes (Part 3.7).

**Taken on trust (not re-derived):**
* the Hoeffding/union-bound outer procedure for M-SHARP (standard);
* Rambachan–Roth's relative-magnitude framework;
* Molinari (2008) for misclassification;
* the concavity of LP value functions in their right-hand side. This is standard; I checked it numerically (N8d) but did not write out a proof.
* the repository's statuses for 2c M2, PHASE2d Δ0/Δ4/CL and P3.5, as given in the technical report §12;
* the review's literature references, which I have not checked.

**Proposed proof-register additions:**

| ID | Claim | Status here | Blocks |
|---|---|---|---|
| R19 | Count layer: LP over a product of transportation polytopes; integral polytope (after rounding rows) for margins, known cells, per-unit boxes and three classes of totals; otherwise the LP is outer and a MILP decides | proved and proof reviewed (PHASE6a) | Step 1 |
| R20 | Audit: separability; concavity of the remaining width; a joint audit by MILP | (a) checked N8a; (b) checked N8d (standard LP sensitivity) | G1 |
| R21 | Rank of the stacked table map in terms of the composition second moments (E4) | checked (NEW_AVENUES §13.1); proof open | Step 11 (E4) |
| R22 | Finite-composition identity; LCE midpoint form; realised bin rows | checked N12a–c; proofs are algebra | Step 3 |
| R23 | Exception-mass frontier: monotone, MILP exact | elementary | Step 7 |
| R24 | Subset identity (2×C, R×2) = Aumann expectation | proof sketch here; checked N1–N3 | Step 8 |
| R25 | Count layer with set-valued margins: coverage ≥ 1−α if P(V ∈ 𝒱) ≥ 1−α; joint margin sets from match probabilities; monotone in a calibration slack | (a) one-line argument; (b), (c) to write | Step 1b |
| R26 | LCE layer with matched (Berkson) composition: known regression correction tr(QΣ_u); allowance for the target's weights; nondifferential matching assumption | correction checked (N19); rest open | Step 1b |
| R29 | C2 contests as districts: identification of a_j, b_j; sharpness by the 2b Thm 2 argument; overidentification count and test; fixed-contest inference | identity checked (N24); proofs to write | Step 10 |
| R30 | M1(δ) neighbour option: nesting and monotonicity; diameter relation with the centre form; components; sharpness via 2b Thm 1′ (R16-type, open) | example checked (N23); (a)–(c) elementary | Step 5 |
| R28 | Covariate-benchmarked κ/δ: identification of the drop-one violation under CAR(Z); validity and monotonicity of the bridge bound; group-rate error bound via the overlap operator in R×C; m* not identified | checked by simulation (N22); proofs to write | Step 6 |
| R27 | Cross-layer intersection: coverage ≥ 1−α−P(belief rows false); disjointness tests model and count rows jointly (level α with identity rows) | proved and proof reviewed (PHASE6a) | Steps 1, 3A |
| R1 (reword) | C1 is needed only for the implication weighted-local ⇒ $\delta_W\ge0$, not to compute M-CONTRAST/M-LOC | J | the M-LOC ⇒ SIGN-W cell only |
| R5, R6 (Blocks column) | now block nothing in the numbered order (deferred) | as before | — |
| R8 (existing) | B1(M): the rank condition (NA Prop. B1), monotonicity in M and κ, the invariant-M LP form | partly (NA, sketch) | Step 4 |
| R9′ | (replaces R9) outer validity of the instrument-averaged M-SHARP procedure, plus a gain study | open | Step 8 (only if an estimator is built) |

---

## Part 9. For applications (outside this repository's plan)

*Your scope rule keeps these out of the plan. They are listed so that they are not lost.*
* **Wisconsin/Iowa.**
  * Assembling the referendum and governor returns.
  * The common base of eligible men in each period, and the changing-population choices (entry, exit, interval bases). The theory's gap register (Step 2) names where each enters.
  * Township selection.
  * Calibration of the veteran match probabilities and the slack (D9; Step 1b's inputs).
  * The township covariate set (D10; Step 6's input).
  * Dial ranges for B1(M), κ and δ.
  * Re-running `apply_B_suffrage`/`apply_kappa_suffrage` with the package's inference once Steps 4–6 exist.
* **NZ.** The leave-one-electorate-out benchmark (Step 11's validation), and the kink-mass diagnostic run (Step 8). These double as package benchmarks, so their *scripts* belong in `benchmarks/`, and their substantive reading belongs in the application.
* **M-15.** The NEW_AVENUES author adds the flags (D8). Re-deriving the two intervals needs the unit-level data and R.
