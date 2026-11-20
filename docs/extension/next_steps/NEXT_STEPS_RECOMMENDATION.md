# What to build next in eiPartialID: a recommendation

*Planning note, 20 November 2026. Written against `claude/wizardly-volta-lbh4bm` at 7077dfd (v0.3.1 on `master`, the round-2 corrections, the technical report and the methods paper). Nothing in the package, the theory documents, the report or the paper is changed. Checks are in this folder: `verify_next_steps.py` (labels N1–N17, 27/27 PASS, output in `verify_next_steps_output.txt`) and two population studies, `sharp_gain_study.py` and `sharp_gain_study_2xC.py` (outputs alongside).*

*How to read it. Part 1 is the one-page summary. Part 3 has one section per recommended step. Each section starts with three plain-language parts (what it is, a worked example, why it comes at this point), which you can read without the technical notes that follow them. Part 4 works through NEW_AVENUES.md avenue by avenue. Parts 5–8 cover what to leave for later, the external review, the decisions needed from you, and what was verified.*

---

## Part 1. One-page summary

**The recommendation in one sentence.** Build next a **count-table (“realised electorate”) layer** beside the existing LCE solver, and put on top of it, in this order, the descriptive Wisconsin/Iowa DiD, finite comparisons, NEW_AVENUES' second-scale context dial (B1(M), with κ-annotated endpoints), exception frontiers, and NZ table transfer. Keep audits, ballot patterns and misclassification as a general-user track for data you do not have. Narrow **M-SHARP** to a diagnostic and a population benchmark. Fold **M-LOC** into a module of finite comparisons.

**Why.** Three findings drive the order.

1. **Most applied questions are about the people who actually voted**, for example “what share of these veterans voted yes?”. Today the package answers that question only through the Duncan–Davis column and `ci_actual`, which come last in `print()`. Every optional restriction acts on the *model-implied* rate (`eiModules()`: zeros, order, sign, slope, linear, external, assoc and cr all restrict the mean). The paper (§7) lists sign rows on realised tables as not implemented. The realised-table results are already proved (2c Prop. M2; DD; PHASE2d Δ0/Δ4), so the gap is in the code, not the theory.
2. **M-SHARP buys little at the population level.** In an oracle study that knows the law of returns exactly (`sharp_gain_study*.py`), using the whole law instead of the regression narrowed the model-implied group rates by a median of **0–1%** and by at most **5%**. That held across 2×2, 2×3 and 2×4 designs, with precinct heterogeneity up to an SD of about .25. In those designs the exact subset rows of the review added nothing beyond the averaged cell limits. In samples, smoothing error will eat some or all of this. The review's subset identity is correct (verified, N1–N3), and it is the right tool for M-SHARP when M-SHARP is built. It does not justify a second-release slot.
3. **Your applications have aggregate margins only.** The count layer serves them directly: the exact DiD rungs are realised-table results, and comparisons and exception frontiers need no extra data. The modules with the largest gains for *other* users bring in information (an audited precinct, a published total, a cast-vote record); in the count layer each is an exact equality or extra columns, so they are cheap to add later (verified, N6–N11).

**What changed after working through NEW_AVENUES (Part 4).** Avenue B's B1(M) enters as Step 4: it serves your Wisconsin/Iowa data directly (the county share is the second scale), and it is the DiD's yardstick. κ annotation joins it. The odds-ratio dial joins Step 1, outcome nesting joins Step 2, and the Bernstein robustness ladder becomes Step 8. Avenue A, the single index and the covariance route stay research; CAR-δ waits on township covariates (D7).

**What changed after the owner's data answers (20 Nov).** No poll books, linked records, CVRs or sub-tables exist for your projects; your applications use aggregate margins only (NZ electorate tables; the Wisconsin/Iowa DiD). The order below therefore ranks modules by their value for that work first, and moves the modules that need external, audit or ballot-level data into a **general-user track** (§3.8). They stay in the plan, because other users have such data, but they are labelled as serving other users.

**Recommended order**

| Step | Module (working name) | What the user gets | Serves | New proofs | Data | Effort |
|---|---|---|---|---|---|---|
| 0 | Finish in-flight work | tag v0.3.1; schedule batches D and E (gates met); cheap print fixes (§3.0) | all | none | none | S |
| 1 | **M-COUNT**: count-table layer, `eiCounts()` | exact bounds for the actual rates, gaps and sub-area totals; realised signs, orderings, odds-ratio dial and CR (2c M2; NA §12.3(c)); known cells and totals; integer check with a witness table. `print()` leads with the actual rate | you + all | R19 (small) | margins | M |
| 2 | **DiD on the count layer, descriptive first** (Δ0, Δ4, composition polytope, party P2/P3 with outcome nesting); causal layer after Batch C | change bounds on the stayer population; ATT only with an explicit trend dial | **you** (Wisconsin/Iowa) | descriptive: tests only; causal: Batch C, R11 | two elections' margins | M, then M |
| 3 | **M-CONTRAST**: finite-composition comparisons, in levels and in changes (M-LOC as its "everywhere" option under LCE; the finite form of NA's no-offsetting ceiling) | signs on *observed* comparisons (low vs high veteran share) with an upper cap on the gap | you + all | R22 (identity; checked) | margins | S–M |
| 4 | **M-CTX2**: context at a second scale, B1(M), with κ-annotated endpoints (NEW_AVENUES B, §12.3(a)) | a dial benchmarked by the *identified* county-level contextual effect; each endpoint labelled with the contextual variation it needs; the DiD's yardstick | **you** (Wisconsin/Iowa county share) | R8 (B1 rank; monotonicity) | a second-scale composition | M |
| 5 | **M-EXCEPT**: exception-mass frontiers | "how many exceptional precincts before my conclusion fails?" | you + all | R23 (small) | margins | M |
| 6 | **M-DATA E4**: shared-*K* transfer, NZ leave-one-electorate-out; also the empirical distribution of κ, δ, M that NZ truths need | the model layer's first real validation; calibration for Step 4's dials | **you** (NZ) | R21 (rank) + allowance | NZ tables | M |
| 7 | **M-SHARP, narrowed**: subset rows, size-conditioned oracle, "kink mass" diagnostic on NZ | when (and by how much) the full law could tighten an LCE bound | you (diagnostic) | R24 (proved here), R9′ | margins | M |
| 8 | **M-BERN**: polynomial outer/inner ladder (with L1′ quadratic) | how much of an LCE answer rests on straight lines | you + all | designed (NA §4) | margins | M |
| G1–G3 | **General-user track**: M-AUDIT (which precinct to validate), M-CVR (ballot patterns), M-MISCLASS (uncertain labels); E5 exact rows ship inside Step 1 | features for users with audits, CVRs or linkage data | **other users** | R20; none; small | external data | S, S, M |
| later | CAR(Z)-δ if township covariates exist (D7); KL-ball reference dial; R16/R15 to take M-SIGN out of experimental (theory track); 2×C joint signs (R6), R×2 trees (R5). Research: Avenue A on pooled NZ, single index (R14), covariance/dispersion | — | — | — | — | — |

**Not now** (Part 5): a series-based M-SHARP estimator; M-LOC as a separate module; CAR-δ without township covariates; dispersion/covariance and Avenue A (research); causal DiD before the descriptive layer and Batch C; switching the default interval (this would reverse a decision of yours); general R×C full-law estimation.

**Decisions needed from you before building starts** (Part 7):
* D1. Which target should lead `print()`: the actual rate or the model-implied rate? (Open; no earlier decision. My recommendation: the actual rate. It reverses nothing.)
* D2. Should the count layer be a new function or a `target =` mode of `eiBounds()`?
* D3. Do you accept moving M-SHARP and M-LOC down?
* D4. The remaining Wisconsin/Iowa DiD data decisions (common base, composition and stayers, township selection), now that it is known no poll books exist. They gate batch C and the causal layer of Step 2.
* D5. How should integers be handled by default?
* D7. Do Wisconsin/Iowa have township-level covariates (e.g. 1860 census aggregates) with veteran share varying within covariate cells? If yes, CAR(Z)-δ becomes a candidate; if no, Step 4 covers the CAR role.

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
| M-SHARP specification | **Partly.** NEXT_BUILD now says "tighter than Θ_O, not sharp" when size is informative | Step 7: a narrower scope and a different place in the order |
| Short proved chains instead of one master ladder | **Not done** | `ladder()` prints proved implications only (Part 3.0) |
| ORD-Δ labelled "MTR" | **Not done** (PACKAGE_DESIGN §7 still says "M-ORD-Δ (MTR)") | rename to "ordering of changes" (Part 3.2) |

---

## Part 3. The recommended steps

### 3.0 Step 0: finish what is in flight (S; no new theory)

* **Plain language.** Release v0.3.1 as a tag and get the code read by someone outside the project. Then make three small changes to what `print()` says, so that what the package claims matches what is proved.
* **Contents.**
  1. Tag `v0.3.1` (NEXT_BUILD item 3 says the tag is still to follow).
  2. Schedule batch E (code review; no reviewer chosen yet) and batch D (design; its gate is met). E is the cheapest batch, and the count layer will reuse the LP code it covers.
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
* For your own work it is the base of the Wisconsin/Iowa DiD: the exact DiD rungs (Δ0, Δ4) are realised-table results.
* Its theory is done: DD, 2c Prop. M2 for realised signs and CR, and PHASE2d Δ0/Δ4.
* It answers the question applied users ask first.
* It needs only aggregate margins, which is all your applications have.
* **Cost:** moderate code (a new LP builder, targets, witnesses and an optional MILP), one small proof item, no data.

**Technical notes.**
* **Object.** Unknown counts $n_{irc}\ge0$ with the observed row and column margins in each unit. The target is any $\sum_i c_i^\top n_i$, divided by a known denominator (or by a sum of known denominators). Restrictions are linear rows on the $n$:
  * structural zeros and known cells;
  * published totals over unit sets (§3.6; a general-user feature);
  * the realised unit ordering $n_{irc}N_{is}\ge n_{isc}N_{ir}$ (N1/B5);
  * the realised association signs $\delta_B,\delta_W$ of 2c M1 (rows in $n$);
  * CR as a union of two LPs;
  * slack dials, as in 2c M2(d);
  * the odds-ratio dial of NEW_AVENUES §12.3(c) for 2×2 and pairwise 2×C: each unit's group-by-outcome odds ratio in $[\psi_L,\psi_U]$, a unit-level interval because the cell is monotone in the odds ratio. Ordering is $\psi_L=1$. In P2 of the example, an odds ratio between 1 and 4 gives the veteran rate [0.500, 0.633] instead of [0.500, 0.833] (N17). Label it untestable when the dial contains 1.
* **Proof items.**
  * Proposed **R19**:
    * (a) the bound is the LP over the product of unit transportation polytopes, intersected with the rows. This is standard.
    * (b) For margins, known cells and per-unit box rows the constraint matrix is totally unimodular, so the LP endpoints of integer-coefficient count targets are integers. This is standard; checked in N10a (300 tables, 0 fractional optima).
    * (c) Rows that couple cells across units or within a unit (pooled realised signs, ordering with fractional coefficients) can break this. Then the LP endpoint is an outer bound, and an integer program decides. N10b gives a one-unit fractional endpoint. N10c gives a three-unit system that is LP-feasible but has **no** integer table, which matches Batch B comment 1.
  * Realised signs and CR for continuous tables: 2c M2 (proved). The integer version: R15, which becomes a computed check rather than a theorem.
* **Inference status.** None is needed. The margins are treated as exact counts. If the margins are themselves estimates (survey margins), the module does not apply. Say so in `print()`.
* **Interface sketch.**
  ```r
  fit <- eiCounts(formula = cbind(yes, no) ~ cbind(vet, other), data = d, unit = "precinct",
                  targets = list(rate("vet", "yes"), gap("vet", "other", "yes"),
                                 rate("vet", "yes", where = county == "Dane")),
                  restrictions = list(order_unit = list(row = "vet", over = "other", col = "yes"),
                                      known = list(unit = "P7", row = "vet", col = "yes", value = 12),
                                      total = list(units = c("P1","P2"), row = "vet", col = "yes", value = 45),
                                      sign_realised = list(type = "dB", sign = ">=0"),
                                      cr_realised = TRUE),
                  integer = c("auto", "lp", "exact"))   # auto: MILP only when rows break TU
  print(fit)        # target, lower, upper, width, integer status, binding rows
  witness(fit)      # unit tables attaining each endpoint
  ```
  Also: `eiBounds(..., target = "actual")` leads its printout with the `eiCounts()` block (decision D1).
* **Tests.**
  * Reproduce `eiBounds()`'s `dd_lower`/`dd_upper` to 1e-12 on every reference dataset.
  * Brute-force enumeration of all integer tables in tiny examples.
  * Closed forms of 2c M2 and M3 in 2×2 (reuse the C7b/C8 checks in `verify_phase2c.py`).
  * Every witness satisfies every row.
  * The N10c system is reported as integer-infeasible.
  * Restrictions only narrow the bound.

---

### 3.2 Step 2: DiD on the count layer, descriptive first, then causal

**What it does (plain language).** For the same men voting in two elections, the descriptive part bounds **how much more** the veterans' support changed than the others' support. Arithmetic does this in each township, optionally with signs on how changes relate to composition. A second, clearly separated step asks the causal question: what would veterans have done had they not served? That needs an explicit belief about how far the two groups' trends could have differed without the war. The package reports the descriptive range first. It reports the causal range only with that belief written down as a dial.

**Information or assumption used.** Descriptive part: the two elections' returns for the same population, plus optional signs on changes. Causal part: a bound on the gap in untreated trends.

**The question it answers.** First, "how much more did veterans' support change than other men's?". Then, separately, "how much of that was caused by service?".

**Worked example** (arithmetic only).
* Suppose the descriptive layer finds that veterans' support rose between 5 and 25 points more than others'.
* If, absent the war, the two groups' trends could have differed by up to 10 points either way, the effect of service is between −5 and +35 points.
* If they could have differed by at most 2 points, it is between +3 and +27.

Without the module the user must assemble this by hand. PACKAGE_DESIGN §7 also labels an ordering of *changes* "MTR", which reads a descriptive fact as a causal one.

**Why here.**
* It is your own application (Wisconsin/Iowa), and it uses aggregate margins only, so nothing in it waits for external data.
* The exact rungs (Δ0, Δ4, the composition polytope, the party results P2/P3) are realised-table results, so they belong on the count layer, and Step 1 makes them cheap.
* PHASE1b and PHASE2d have not been reviewed. The descriptive exact rungs rest on short arguments (product sets, 2c M3 applied to changes) that the count layer's brute-force tests check directly, so they can be built before Batch C. The causal layer and the FD rungs should wait for Batch C, which is gated on the open DiD data decisions (PHASE1b §D; decision D4).
* **Cost:** M for the descriptive layer (after Step 1); M for the causal layer (after Batch C).

**Technical notes.**
* **Descriptive.**
  * Δ0: product of the two periods' unit sets.
  * Δ4: 2d Prop. Δ4, rows on changes.
  * Composition: 2d §7, interval bases in count coordinates.
  * Party: 2d P2/P3, plus outcome nesting with slack ε ("suffrage yes ⇒ Republican, except ε"; NEW_AVENUES §12.3(d)). It changes nothing for one outcome but narrows cross-outcome contrasts to a median 0.39 of their width [NA, `check_nesting*`].
* **Causal.** $ATT=D-b$, with $b$ the untreated trend gap, and $\lvert b\rvert\le M$ or relative-magnitude constraints on one feasible sequence of pre-trends (Rambachan–Roth). 2d Prop. CL (MTR+ND+PS) is an alternative layer. Label it causal, and keep the "ordering of changes" rows separate.
* **Proofs.** Batch C, R11. DiD-S validity is proved but not reviewed. The integer versions of Δ4 fall under R15, which the count layer's integer check covers computationally.
* **Inference.**
  * Exact rungs: none.
  * FD rungs: the existing machinery (CLR algorithmic, projection).
  * Causal dial: deterministic in $M$ given the descriptive set.
* **Interface**: `eiChange(d0, d1, unit = "township", stayers = ..., restrictions = list(cr_change = TRUE))`; `eiATT(change_fit, trend_gap = 0.1, pre = NULL)`.
* **Tests**:
  * 15/15 of `verify_phase2d.py` ported;
  * brute force on small two-period tables;
  * placebo on pre-periods, reported with its width (a wide placebo that contains zero is weak evidence).

---

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
* It needs no data beyond the margins, so it serves your NZ and Wisconsin work directly. Its version on *changes* (comparing low- and high-veteran townships' changes) is the finite form of DiD-S (PACKAGE_DESIGN §7) and needs no linearity in either period.
* It comes after Steps 1–2 because those are exact; this one adds beliefs.
* **Cost:** S–M, one small proof item.

**Technical notes.**
* **Identity** (Kitagawa 1955):
  $$\tau_z-\tau_x=(z-x)^\top\bar m+\bar x^\top(m_z-m_x).$$
  For $z-x=t(e_r-e_s)$, the two signs confine each term to $[0,\tau_z-\tau_x]$.
* **Count layer.** Bins are user-chosen sets of units. Bin-level group rates are linear-fractional in counts with known denominators, so they are linear. The rows are exact for those bins' actual voters. In R≥3, use only comparisons that the user names. A general displacement has a well-defined mix term, but attributing it to one pair is an assumption.
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

**Worked example** (your data; [NA, `apply_B_suffrage`, `apply_kappa_suffrage`]; trusted; plug-in, no confidence intervals). Wisconsin 1865, trimmed units, veteran minus non-veteran yes rate.
* *Without the module*: pooled LCE gives **[−0.066, 0.452]**, which does not settle the sign.
* *With B1(M = 1)* on the yes column: **[0.049, 0.297]**, a positive gap. With M = 0.5: [0.111, 0.235].
* *κ annotation* (same data, pooled LCE; see the κ table in §4.8(a) for the cap's effect): the lower end of −0.065 needs κ above 0.25. Capping κ at 0.25 lifts the lower end to 0.019. The data's own floor is κ = 0.015.
* *A warning that belongs in the output*: imposed on all three columns jointly, B1 is empty at the plug-in unless M ≥ 3.85 (1857) or 6.65 (1865). Iowa 1868 is empty even under plain M0. This is why the step includes the emptiness test and proper intervals, not plug-in sets.

**Why here.**
* It is NEW_AVENUES' first choice for near-term use (§9), and it applies directly to **your Wisconsin/Iowa data**, where the county share is the second scale.
* It gives the DiD its yardstick. M-CTX2-Δ replaces FD3's free δ with a measured contextual effect (PACKAGE_DESIGN §7). It is also the benchmark of the right kind for any CAR-style sensitivity (NEW_AVENUES §14.5).
* It comes after Steps 1–3 because those are exact, or use only directional beliefs. This step adds a magnitude belief on the LCE layer, and it needs a proof first.
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

### 3.5 Step 5: M-EXCEPT, how many exceptions can a conclusion survive?

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

### 3.6 Step 6: known tables (M-DATA): E4 transfer on NZ for you; exact uses (E5) as general-user features

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
* **Transfer (E4), the part that serves your NZ work.** Two electorates under the straight-line model (two groups): electorate A has precincts with 10% and 30% veterans, electorate B has precincts with 60% and 90%. *Without the module* the model-implied veteran rate in B is in **[0.529, 0.751]**. *With A's published table*, and the assumption that the hidden pattern is the same in both, B's rate is pinned at **0.656** (with two groups one table fixes the pattern completely, NEW_AVENUES E1). *Allowing the pattern to differ by an amount worth 2 points on A's rate*, B's rate is in **[0.650, 0.662]** (N15). How large a difference was needed in past NZ elections is exactly what leave-one-electorate-out measures; NEW_AVENUES §5.2 suggests it may be large, in which case the gain shrinks accordingly.

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

### 3.7 Step 7: M-SHARP, narrowed to a benchmark and a diagnostic

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

**Why it comes last of the numbered steps.** (Its diagnostic can run on your NZ data; that is the one place to look for the size channel.)
* In an oracle study that knows the law exactly, it narrowed the model's group rates by a median of 0–1% and at most 5%:
  * 2×2: 36 designs, median ratio 1.00, minimum 0.95;
  * 2×3/2×4: 24 designs, minimum 0.97–0.98 with 40,000 draws per type.
* In those designs the exact subset rows added nothing over averaged cell limits.
* With 8,000 draws per type the same study showed *spurious* gains of up to 8%. These come from sampling noise in the right-hand sides being read as information, which previews the inward bias any plug-in version will have.
* So it costs a lot (smoothing, size bins and a new inference argument) for little expected width. Its value is as a **benchmark** (how much is being left on the table?) and a **diagnostic** (the "kink mass" of NEW_AVENUES §6) that flags the rare data where it matters.
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
* **Interface**: `eiSharpCheck(fit, by = c("x", "x_n"), bins = ..., oracle = FALSE)` reports the kink mass and the population (or binned) Θ_X/Θ_O width ratio per target. There is no change to `eiBounds()`'s default.
* **Tests**: N1–N5; recovery of the oracle on replicated designs; the Prop. 3 equality cases; the 2×4 coordinate failure.

---

### 3.7b Step 8: M-BERN, the polynomial robustness ladder (with L1′)

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

**Why here.** It is the credibility check for every LCE number in your NZ and Wisconsin papers, but it adds no tightening. It therefore comes after the steps that add information or beliefs. It shares its machinery with L1′ (quadratic contextual effects), which is on PHASE4's second-release list. **Cost:** M (the outer/inner LP pair, the dimension formula and degree elevation); no new proof beyond NEW_AVENUES §4.2's design.

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

### 3.8 General-user track: modules for data you do not have

*You have told us there are no poll books, linked records, CVRs or sub-tables for your projects. The modules below need such data. They are worth keeping because other users have it (archives with poll books, US cast-vote records, linked historical records), and each is small once Step 1 exists. Build them after Step 5, or earlier if a user with the data asks; their order among themselves follows effort: G1, then G2, then G3. The exact E5 rows (totals, strata) belong here too but ship inside Step 1 (§3.6).*

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

#### G3. M-MISCLASS, uncertain group labels (general-user feature)

**What it does (plain language).** When "veteran" comes from linking names across records, some men are labelled wrongly: a few veterans are missed and a few non-veterans are counted as veterans. If an audit says how many labels can be wrong, the module widens the bounds just enough to cover those errors. The answer then stays honest about the group's definition instead of looking more precise than the data are.

**Information or assumption used.** An upper bound on the number of wrong labels (false and missed matches), from a clerical audit or a second register. Errors may depend on how men voted unless evidence says otherwise.

**The question it answers.** "How much could my result move if some group labels are wrong?"

**Worked example.** Precinct P2 records 60 veterans among 100 men, with 50 yes votes. *Without the module*, treating every label as correct, the veterans' yes rate is in [0.167, 0.833]. *With the module* and up to 2 mislabelled men, it is in **[0.138, 0.862]** (N14). The module widens rather than narrows. Its value is that a narrow model-based interval is not taken at face value when the group itself is measured with error.

**Why here (general-user track, last).** No linkage audit exists for your projects, so it serves other users only. It is the most expensive of the three general-user modules (ratios with unknown denominators), and it matters only when a user has an audit.

**Technical notes.**
* Counts by true group × recorded label × outcome, with the observed recorded-label and outcome totals and audit bounds on false and missed matches.
* A rate whose denominator (the true group size) is unknown is linear-fractional: profile over integer group totals, or use a Charnes–Cooper transformation for a single ratio. A difference of two such ratios needs joint profiling.
* Start at the accounting rung (Step 1). Plugging an uncertain composition into the LCE polynomial breaks the simple LP.
* Framework: Molinari (2008). Proof item: small (validity of the profiled LP; monotone in the error budget).
* **Inference**: deterministic given exact audit counts; a design-based region if the audit is a sample.
* **Interface**: `eiCounts(..., misclass = list(row = "vet", false_pos = 2, false_neg = 2, outcome_dependent = TRUE))`.
* **Tests**: N14; brute force on small tables; bounds widen monotonically in the error budget and reduce to Step 1 at zero.

* **R16/R15 (M-SIGN out of experimental).** Unchanged. They are on the theory side, they run in parallel, and they are not a build step.

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

**Verdict: research only** (not in the build order). It serves **your NZ data**, but only if the national file is pooled. That requires a cross-electorate model (a common mean structure), which is itself a strong assumption. A2, the source of most of the gain, is the least defensible piece (§1.4).
* *Depends on*:
  * moment-inequality inference for conditional variance functions (Andrews–Shi type);
  * the R×C version of A2 (a conjecture);
  * a decision on how to pool NZ across electorates.
* A0: **drop** (implied by Θ_X, never binds).

**Fit.**
* LCE layer: A1/A2 are SOCP or LMI rows, or LP rows with ordering in 2×2.
* No change to the order. A good self-contained research project on NZ once Steps 1–8 are done.
* A0's relation to Step 7: it is a weaker form of the same full-law information.

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

**Verdict: build now, as new Step 4** (with κ annotation, §4.8(a)). It serves **your Wisconsin/Iowa data**: the county share is the second scale, and the effects are identified. It is NEW_AVENUES' first choice for near-term use, and it gives the DiD its dial. M-CTX2-Δ replaces FD3's free δ with a measured yardstick (PACKAGE_DESIGN §7).
* *Depends on*:
  * R8, the B1 rank condition (NA Prop. B1, a sketch). Proving it is the first task of the step;
  * the parametrisation-invariant M (§2.3; the Wisconsin numbers above use the entry-wise form);
  * the existing inference machinery: few rows, so CLR and projection apply under pooled M0 by R3/R10. Grouped versions follow R16.

**Fit.**
* LCE layer only.
* It changes the order: it enters at Step 4, ahead of M-EXCEPT.
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
* **C2**: build later. It reuses the M1 intersection machinery, with contests in place of districts.

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

**Verdict: build, as Step 8** (after Step 7). It serves **your data and others'**: it is the natural robustness display for every LCE result in the NZ and Wisconsin papers.
* *Depends on*: the outer/inner pair (designed, NA §4.2). The degree must be prespecified (2c §6 counterexample).
* L1′ (quadratic contextual effects, on the repository's list) is its J=2 case. Build them together.
* It does not tighten anything, so it does not compete with the steps above for user value. It does compete for credibility, so it should not slip far.

**Fit.** LCE layer. Step 8 in the revised order (§3.7b).

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
* **E4: build now as Step 6** (yours, NZ).
* **E1**: part of E4.
* **E2**: drop as an identifying rung, and keep as a dial within E4. The NZ evidence is weak.
* **E3**: proved, nothing to build.
* **E5**: other users'. Its rows ship inside Step 1.
* **E6**: other users' (no validation data for your projects). Build later as part of the general-user track, alongside G1 (it is the same "recovered information" logic plus the Γ dial).

**Fit.** E4: LCE layer with the actual-vs-model allowance. E5/E6: count layer. No change to the order.

### 4.6 Avenue F: computing the sharp set (§6; relation to the M-SHARP oracle study)

**In plain language.** Use the whole spread of returns at each composition, not only the average. (Same as Step 7.)

**Worked example.** N5 (Step 7): a group's mean at a 50/50 composition goes from [0, 1] to [0.4, 0.6] when half the precincts vote 10% and half 90%. Across realistic LCE designs, though, the gain is a median of 0–1% (my oracle study).

**Verdict: build later, as a diagnostic (Step 7).** It serves **your NZ data** as the place to test the size channel. Two pieces of §6 are now settled or redirected:
* The "finite core-determining class" conjecture is **proved** for 2×C and R×2 by the review's subset identity (my R24). The class is the proper subsets of outcomes (or of groups).
* The "kink-mass" diagnostic is the right trigger for running the oracle at all.

**Fit.** LCE layer. Unchanged position.

### 4.7 The lower-ranked avenues (§7)

| Avenue | Plain language | Small example | Verdict, and whose data |
|---|---|---|---|
| **G, label symmetry** | Assume the model treats, say, Labour and National symmetrically. This roughly halves the unknown directions | In a symmetric yes/no referendum the free directions halve. In NZ's asymmetric party system the assumption is hard to defend | **Drop** as an identifying rung. At most a labelled sensitivity check for symmetric referenda (other users) |
| **H, deconvolution** | Recover the law of precinct shocks from the spread of returns | Proved: the shock law does not involve K. It trims only the support (m ≥ −min e), a small effect | **Drop** |
| **I, spatial structure** | Neighbouring polls as replicates for Avenue A's variance estimates | Neighbours with similar composition give local estimates of Var(T \| x, n) without a global polynomial | **Research only**, inside Avenue A (yours, NZ, if A proceeds). Neighbour-only *similarity* gives no gain (§12.3(f), `sim_local_near_constancy`) |
| **J, entropy / OT / Choquet** | Pick the "most random" table inside the bounds | Maximum entropy at fixed margins *is* "no within-precinct association": it selects a point, it does not narrow | **Drop** as identification. The honest form, distance to an explicit reference table (KL ball, §12.3(e)), is in §4.8(e). Choquet capacities are the computational form of R24 (Step 7) |

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
* **Verdict: research only for now.** It could be large for **your NZ 3×3 tables**, and E4's leave-one-out is the natural test. But R14 is a conjecture and the computation is nonconvex (an LP per direction on a grid). Revisit after Step 6.

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

*(f) All-pairs similarity.* **Drop.**
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
* **Validation tables calibrate the dials.** Step 6 (E4 on NZ) should output the empirical distribution of κ, δ and M that the NZ truths require. This turns E4 into a calibration source for Step 4's dials.
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

**Verdict: build later, conditional on covariates.**
* **Your data.** It needs township-level covariates with composition varying within covariate cells (positivity). Whether Wisconsin/Iowa have them (for example 1860 census township aggregates) is a new question for you (decision D7). With county-level Z only, CAR collapses into Avenue B (§14.7), which Step 4 already delivers.
* **Other users'**: build after Step 6 if they have rich covariates.
* *Depends on*: R8 (validity, monotonicity), the overlap factor beyond 2×2 (a conjecture), and the SOCP form of the κ/r² dials.
* **Cinelli–Hazlett / Chernozhukov et al. contours**: **do not build**. Point users to `seine`, and print §14.4's caution when they do.
* Covariate adjustment is low priority by your decision. Step 1's strata option already gives the assumption-free version for users with stratified tables.

**Fit.** LCE layer (CAR-δ is an LP with a dial, centred on Goodman-within-Z). It is not in the numbered order.

### 4.10 Covariance of unit-level rates (§13.5, §14.6)

**In plain language.** If how much group rates vary around their means does not depend on composition, the variance of returns reveals how the groups' rates co-vary across places, including their correlation. It does *not* reveal the hidden pattern.

**Worked example** [NA; trusted; the interval is under M-15].
* Wisconsin 1857: the estimated covariance matrix is not even valid (implied "correlation" 3.35).
* Wisconsin 1865: the correlation is −0.16, with an interval of [−0.49, 0.70] (cluster bootstrap, pre-D3).

So the data cannot sign the correlation. Without the module, A2 must *assume* it is non-negative. With the module, it would be estimated, but here it is not estimable.

**Verdict: research only** (local inference and an r² dial need many units with good overlap). For your data: Wisconsin evidence negative; NZ untested. It does not move aggregate bounds (§13.7 D). Its two eventual uses:
* replacing A2's sign assumption inside Avenue A research;
* a printed diagnostic: the PSD test, reported before any clipping.

**Fit.** Unit or second-moment level (local). It belongs with Avenue A.

### 4.11 Excluded composition shifters (§13.6)

**In plain language.** Something that changes *who lives where* without changing *how anyone votes* would identify the group rates, like an instrument. In EI, almost anything that changes the mix of people around a voter also changes their context, so the exclusion is hard to defend. The clean case is re-drawing *the same people* into different reporting units. That is crossing partitions (E5(b)), which needs no exclusion.

**Worked example** [NA, `check_review_items`, reproduced]. Crossing partitions moved a bound from [0.098, 0.939] to [0.118, 0.930] with homogeneous atoms, and barely moved it otherwise. A true shifter has no example in the repository.

**Verdict: drop** as a module. The clean case lives in E5 (other users', count layer). Boundary changes over time add turnover (A1).

### 4.12 Reconciliation with NEW_AVENUES' own ranking

| NEW_AVENUES ranking (§9; revised in §13.8 and §14.7) | Their rank | Mine | Why we differ |
|---|---|---|---|
| E4 several validation tables (+δ) | 1 | Step 6 | Same direction. I put the count layer and the exact DiD first because your Wisconsin/Iowa work needs no tables, and the E4 rank result needs a proof (R21) plus the allowance |
| E5 stratified / crossing / multi-way margins | 2 | Inside Step 1 (other users) | You have no such tables. The rows are nearly free once the count layer exists, so they ship with it rather than as a step |
| B second-scale context, B1(M) | 3 (1st in §9) | **Step 4** | Agree that it is the most practical new rung for your data (county share in Wisconsin/Iowa). The plug-in emptiness on all columns makes inference (CLR or projection plus the emptiness test) part of the step, not an afterthought |
| CAR(Z)-δ | 4 | Later, conditional on township covariates (D7) | With county-level Z it reduces to B. Covariate adjustment is low priority by your decision |
| E6 validation microdata with Γ | 5 | General-user track | No validation data for your projects |
| §13.5 covariance → dispersion ceilings | 6 | Research (with A) | Wisconsin estimate not PSD (1857) or unsigned (1865); local only |
| No-offsetting ceiling; single index; κ (§12) | 7 | No-offsetting → Step 3 (M-CONTRAST); κ → Step 4; single index → research | M-CONTRAST is the finite, derivative-free form of the no-offsetting ceiling. κ annotation is cheap and pairs naturally with B1(M). The single index rests on a conjecture (R14) and is nonconvex |
| A coin-flip rungs | 8 | Research | Agree: strong identification, weak estimability. `sim_size_grid` suggests NZ's national file is the one place it could work, so I would name that as the research project |
| Excluded shifters | 9 | Drop (clean case in E5) | Agree |
| D Bernstein ladder | 5th (§9) | Step 8 (robustness) | Agree that it is infrastructure. I give it a firmer place because it is the credibility check for every LCE number in your papers |
| F sharp set | 6th (§9), "pinned" | Step 7, diagnostic only | The oracle study quantifies §9's "power 2": median gain 0–1% in LCE designs |
| G–J | low | Drop (G as a sensitivity check; I's replicates inside A) | Agree |
| §12.3(c) odds-ratio dial | §12.5 item 5 | **Step 1** | It is a count-layer row, and ordering is its special case, so it costs nothing extra |
| §12.3(d) nesting | not ranked | **Step 2** (party layer) | It gives gains for the cross-outcome contrasts your party analysis uses |

**Where NEW_AVENUES is now superseded or needs a note.**
* §6's finite core-determining class is proved for 2×C/R×2 (R24).
* §12.5 item 1 is delivered by M-CONTRAST in finite form.
* The M-15 intervals still need re-deriving.
* The Wisconsin B1 results should be recomputed with the package's inference once Step 4 exists, before any of them is quoted.

---

## Part 5. Not now, and why

| Item | Why not now | What would change that |
|---|---|---|
| **M-SHARP series estimator** (NEXT_BUILD item 4, R9 as stated) | The oracle gain is median 0–1%, and at most 5% (`sharp_gain_study*.py`). Plug-in versions are biased inward. "Tightens every rung" is not a sensible release criterion | Step 7's diagnostic finds data with large kink mass, or the NZ benchmark shows gains |
| **M-LOC as its own module** | It is the "everywhere" option of M-CONTRAST under LCE (N12c). C1 (R1) is needed only for the implication to $\delta_W$ | none; build it with Step 3 |
| **2×C joint/ordinal signs (R6), R×2 spanning trees (R5)** | Evidence of user demand is thin, and R5(c) and R6 are open or only checked. These are LCE-only refinements | an application with ordered outcomes or ≥3 groups that asks for them |
| **CAR-δ, M-DISP/covariance, Avenue A, single index, tensor contests** (see Part 4 for each) | CAR-δ needs township covariates with overlap (D7); A and the covariance route need thousands of units and a pooling model; the single index rests on conjecture R14 and is nonconvex; tensor contests need implausible local independence | township covariates (CAR); a decision to pool NZ nationally (A); a proof of R14 |
| **Causal DiD before the descriptive layer and Batch C** | PHASE1b and PHASE2d are unreviewed; the causal labels need explicit trend assumptions | Batch C returned |
| **Switching the default interval to projection** | **This would reverse your decision.** I do not recommend it now. 0.3.1 already prints projection beside CLR with a per-fit status, and Step 0 adds the R2-regime line | Conjecture P3.5(b3) fails in a design users meet, or a referee insists |
| **General R×C full-law estimation** | The cost is high and the gain unproven (above) | Step 7 shows gains |
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
10. **B1(M) where a second scale has adequate variation.** The review places it at its step 5. I place it at Step 4, because Wisconsin/Iowa's county share is that second scale and it is the DiD's yardstick.

**Disagree, or agree with qualifications:**
1. **On M-SHARP's place.** The review puts M-SHARP second, before M-LOC and after the count layer. I put it last of the numbered steps and narrow it to a diagnostic. The review's argument rests on population information. My oracle study finds that information small in LCE designs (≤5% width; median about 0). The review itself notes that finite-sample bands can erase it. I agree with the review against the repository on narrowing ("remove 'tightens every rung'") and on conditioning on size. I disagree with the repository's NEXT_BUILD item 4 on priority.
2. **On the default interval.** The review would make the coverage claim pick the default, with projection where proved. You decided to keep CLR as the default, labelled algorithmic, with projection beside it. 0.3.1 implements that transparently. I would not reverse it now. What remains of the review's point is cheap and does not change the default: the R2-regime line (Step 0). Revisit if Conjecture P3.5(b3) fails in a realistic design.
3. **On E4 before E5.** The review's build order promotes E5 and E4 together. I would ship the exact uses (E5, published totals, strata) in the count layer first. They need no allowance and no rank proof. Transfer (E4) is an assumption and should follow with its δ dial and NZ evidence. NEW_AVENUES §5.2 E2 already found weak persistence of $K$ in NZ.
4. **On "short proved chains" replacing the ladder.** Agreed for `ladder()`'s printout (Step 0). The plausibility tiers can stay in the documents as an *elicitation aid*, labelled as such, because users find them helpful.
5. **On the review's ordering of M-LOC after the validation work.** I agree it is late. I would not build it separately at all (point 5 under Agree).
6. **On the priority of targeted validation and external data.** The review ranks targeted validation "high" and builds the audit calculation into its first step. For your projects no audit, linkage, CVR or sub-table data exist, so these become a general-user track (§3.8) behind the modules that work with aggregate margins only. The review's other high-priority idea, finite comparisons, needs no data and moves up to Step 3.

**Where I disagree with the repository's own plan:**
* PACKAGE_DESIGN §6 orders LCE refinements (M-LOC, 2×C, R×2) before data modules on grounds of proof readiness. The count-layer modules are at least as proof-ready: standard LP facts plus 2c M2, which is proved. They also answer the users' first question. So that rationale no longer favours the LCE refinements.
* NEXT_BUILD item 4 says M-SHARP "tightens every rung with no new assumption". True in principle, but negligible in the designs tested.

---

## Part 7. Decisions you need to make before building starts

* **D1. Print order (open: you have made no earlier decision on this).** Should `print.eiRxC()` lead with the actual-rate block (DD, then any count-layer restrictions, then `ci_actual`/`proj_actual`) and put the model-implied block second? *My recommendation: yes, as the default, with `target = "model"` restoring the current order.* This does not touch your CLR decision.
* **D2. One function or two.** Either a new `eiCounts()` that `eiBounds()` calls for its actual-rate block, or a single `eiBounds(target = )` with two engines. *Recommendation: a separate `eiCounts()` (clean semantics: exact, no inference) that `eiBounds()` embeds.*
* **D3. Re-ordering.** Do you accept moving M-SHARP (to a diagnostic) and M-LOC (into M-CONTRAST) below the count-layer modules? This changes `NEXT_BUILD.md` item 4 and PACKAGE_DESIGN §6 steps 1b–4.
* **D4. The DiD data decisions for Wisconsin/Iowa** (PHASE1b §D), restated for what is now known. Your own driver is Wisconsin/Iowa with **aggregate margins only** (no poll books, no linked service records), plus NZ. So the DiD is built entirely on the exact rungs (Δ0, Δ4, the composition polytope, party P2/P3) and on FD, with no validation data, and audit-based ideas do not apply to it. The decisions still open are: the common base of eligible men in each period; how to treat aging-in and exits (stayers, interval bases from the 1860/1870 censuses); and which townships enter (the selection point in PHASE1b §D.6). Once these are fixed, batch C can be sent, and the causal layer of Step 2 can follow it. The descriptive layer does not wait for them.
* **D5. Integers.** Should the default be `integer = "auto"` (MILP only when rows break total unimodularity, otherwise the LP with a certificate)? Or always LP with a flag? *Recommendation: "auto", with a time limit that falls back to the LP bound plus a labelled gap.*
* **D7. Township covariates for Wisconsin/Iowa.** Are township-level covariates available, with veteran share varying within covariate cells (positivity)? If yes, CAR(Z)-δ (NEW_AVENUES §13.4, §14) becomes a candidate after Step 6, with B1(M) as its benchmark. If no, Step 4 already covers the CAR role at the county scale, and CAR-δ stays a general-user feature.
* **D6 (minor).** Batch D's gate is met but it is unscheduled, and batch E has no reviewer. Send batch D with this recommendation attached, so the reviewer sees the revised order; choose a reviewer for batch E.

---

## Part 8. What I verified, and what I take on trust

**Verified numerically here** (`verify_next_steps.py`, 27/27 PASS):

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
| N15 | E4 transfer example (B's rate [0.529, 0.751] → 0.656; ±2 points on A → [0.650, 0.662]) | confirmed |

**Population studies** (`sharp_gain_study.py`, `sharp_gain_study_2xC.py`). Pooled LCE, oracle law at replicated composition types, large precincts, Beta or Dirichlet heterogeneity. The result is that the full law narrows the model-implied group rates by little.

*Limits of these studies*:
* only 2-row designs;
* no binomial noise;
* no size channel (all precincts large);
* linear truths.

The size channel (N4) is the main way a larger gain could appear. Step 7's diagnostic should look for it in NZ.

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
| R19 | Count layer: LP over a product of transportation polytopes; TU ⇒ integer endpoints for margins, known cells and per-unit boxes; otherwise the LP is outer and a MILP decides | (a, b) standard, checked N10a; (c) checked N10b–c | Step 1 |
| R20 | Audit: separability; concavity of the remaining width; a joint audit by MILP | (a) checked N8a; (b) checked N8d (standard LP sensitivity) | G1 |
| R21 | Rank of the stacked table map in terms of the composition second moments (E4) | checked (NEW_AVENUES §13.1); proof open | Step 6 (E4) |
| R22 | Finite-composition identity; LCE midpoint form; realised bin rows | checked N12a–c; proofs are algebra | Step 3 |
| R23 | Exception-mass frontier: monotone, MILP exact | elementary | Step 5 |
| R24 | Subset identity (2×C, R×2) = Aumann expectation | proof sketch here; checked N1–N3 | Step 7 |
| R8 (existing) | B1(M): the rank condition (NA Prop. B1), monotonicity in M and κ, the invariant-M LP form | partly (NA, sketch) | Step 4 |
| R9′ | (replaces R9) outer validity of the instrument-averaged M-SHARP procedure, plus a gain study | open | Step 7 (only if an estimator is built) |
