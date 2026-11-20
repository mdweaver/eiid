# What to build next in eiPartialID: a recommendation

*Planning note, 20 November 2026. Written against `claude/wizardly-volta-lbh4bm` at 7077dfd (v0.3.1 on `master`, the round-2 corrections, the technical report and the methods paper). Nothing in the package, the theory documents, the report or the paper is changed. Checks are in this folder: `verify_next_steps.py` (labels N1–N14, 24/24 PASS, output in `verify_next_steps_output.txt`) and two population studies, `sharp_gain_study.py` and `sharp_gain_study_2xC.py` (outputs alongside).*

*How to read it. Part 1 is the one-page summary. Part 3 has one section per recommended step. Each section starts with three plain-language parts (what it is, a worked example, why it comes at this point), which you can read without the technical notes that follow them. Parts 4–7 cover what to leave for later, the external review, the decisions needed from you, and what was verified.*

---

## Part 1. One-page summary

**The recommendation in one sentence.** Build next a **count-table (“realised electorate”) layer** beside the existing LCE solver, and put the information-based modules on top of it: auditing, known tables, finite comparisons, exceptions, and descriptive DiD. Narrow **M-SHARP** to a diagnostic and a population benchmark. Fold **M-LOC** into a module of finite comparisons.

**Why.** Three findings drive the order.

1. **Most applied questions are about the people who actually voted**, for example “what share of these veterans voted yes?”. Today the package answers that question only through the Duncan–Davis column and `ci_actual`, which come last in `print()`. Every optional restriction acts on the *model-implied* rate (`eiModules()`: zeros, order, sign, slope, linear, external, assoc and cr all restrict the mean). The paper (§7) lists sign rows on realised tables as not implemented. The realised-table results are already proved (2c Prop. M2; DD; PHASE2d Δ0/Δ4), so the gap is in the code, not the theory.
2. **M-SHARP buys little at the population level.** In an oracle study that knows the law of returns exactly (`sharp_gain_study*.py`), using the whole law instead of the regression narrowed the model-implied group rates by a median of **0–1%** and by at most **5%**. That held across 2×2, 2×3 and 2×4 designs, with precinct heterogeneity up to an SD of about .25. In those designs the exact subset rows of the review added nothing beyond the averaged cell limits. In samples, smoothing error will eat some or all of this. The review's subset identity is correct (verified, N1–N3), and it is the right tool for M-SHARP when M-SHARP is built. It does not justify a second-release slot.
3. **The new modules with the largest gains bring in information, not assumptions**: an audited precinct, a published total, an age-stratified table, a cast-vote record. In the count layer each enters as an exact equality or extra columns. None needs an allowance or a new inference theorem, and most need only standard LP facts (verified, N6–N11).

**Recommended order**

| Step | Module (working name) | What the user gets | New proofs | Data | Effort |
|---|---|---|---|---|---|
| 0 | Finish in-flight work | tag v0.3.1; send review batch E (code); cheap print fixes (Part 3.0) | none | none | S |
| 1 | **M-COUNT**: count-table layer, `eiCounts()` | exact bounds for the actual rates, gaps and sub-area totals; realised signs, orderings and CR (2c M2); known cells; integer check with a witness table. `print()` leads with the actual rate | R19 (small) | none | M |
| 2 | **M-AUDIT**: `eiAudit()` | which precinct to validate, with the guaranteed narrowing computed in advance | R20 (small; checked) | costs | S |
| 3 | **M-DATA, exact first**: published totals, stratified and crossing margins (E5); then shared-*K* transfer (E4) with NZ leave-one-out | sub-area bounds from coarser tables; covariate-stratified DD with no assumption | E5: none. E4: R21 (rank proof) + allowance | published tables | M |
| 3b | **M-CVR**: ballot patterns (if there is an application) | group-specific split-ticket and roll-off bounds | none | CVRs | S |
| 4 | **M-CONTRAST**: finite-composition comparisons (M-LOC as its "everywhere" option under LCE) | signs on *observed* comparisons (low vs high veteran share) with an upper cap on the gap | R22 (identity; checked) | none | S–M |
| 5 | **M-EXCEPT**: exception-mass frontiers | "how many exceptional precincts before my conclusion fails?" | R23 (small) | none | M |
| 6 | **DiD, descriptive first** (Δ0, Δ4 and the composition polytope on the count layer); causal layer afterwards | change bounds on the stayer population; ATT only with an explicit trend dial | Batch C review; R11 | two elections | M, then M |
| 7 | **M-SHARP, narrowed**: subset rows, a size-conditioned oracle, a "kink mass" diagnostic | when (and by how much) the full law could tighten an LCE bound | R24 (identity proved here), R9 restated | none | M |
| later | M-MISCLASS (when an audit exists); R16/R15 to take M-SIGN out of experimental; 2×C joint signs (R6), R×2 trees (R5) | — | — | — | — |

**Not now** (Part 4): a series-based M-SHARP estimator; M-LOC as a separate module; B1(M), κ, CAR-δ and dispersion; causal DiD before the descriptive layer and Batch C; switching the default interval (this would reverse a decision of yours); general R×C full-law estimation.

**Decisions needed from you before building starts** (Part 6):
* D1. Should `print()` lead with the actual rate?
* D2. Should the count layer be a new function or a `target =` mode of `eiBounds()`?
* D3. Do you accept moving M-SHARP and M-LOC down?
* D4. Which application should be the first test case: Wisconsin with an audit, NZ tables, or CVRs?
* D5. How should integers be handled by default?

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

  So the build order in PACKAGE_DESIGN §6 has not had an outside check.
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
| ORD-Δ labelled "MTR" | **Not done** (PACKAGE_DESIGN §7 still says "M-ORD-Δ (MTR)") | rename to "ordering of changes" (Part 3.6) |

---

## Part 3. The recommended steps

### 3.0 Step 0: finish what is in flight (S; no new theory)

* **Plain language.** Release v0.3.1 as a tag and get the code read by someone outside the project. Then make three small changes to what `print()` says, so that what the package claims matches what is proved.
* **Contents.**
  1. Tag `v0.3.1` (NEXT_BUILD item 3 says the tag is still to follow).
  2. Send batch E (code review) now. It is the cheapest batch, and the count layer will reuse the LP code it covers.
  3. In regime R2 (many small groups), print the effective number of units for the focal group, $1/\sum_iw_{ir}^2$, and the width of the finite-sample interval next to the DD width. This adds information without changing the default.
  4. Have `ladder()` print only proved implications (for example M1 ⊂ M1δ ⊂ M2 and "signed slopes ⇒ weighted-local") and list the other rungs as side by side.
  5. Rename "M-ORD-Δ (MTR)" in PACKAGE_DESIGN §7 to "ordering of changes". MTR is a causal assumption.
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
* Most later modules (auditing, known tables, ballot patterns, exceptions, descriptive DiD) are a few rows or columns added to this layer.
* Its theory is done: DD, 2c Prop. M2 for realised signs and CR, and PHASE2d Δ0/Δ4.
* It answers the question applied users ask first.
* It needs no data beyond what users already have.
* **Cost:** moderate code (a new LP builder, targets, witnesses and an optional MILP), one small proof item, no data.

**Technical notes.**
* **Object.** Unknown counts $n_{irc}\ge0$ with the observed row and column margins in each unit. The target is any $\sum_i c_i^\top n_i$, divided by a known denominator (or by a sum of known denominators). Restrictions are linear rows on the $n$:
  * structural zeros and known cells;
  * published totals over unit sets (Step 3);
  * the realised unit ordering $n_{irc}N_{is}\ge n_{isc}N_{ir}$ (N1/B5);
  * the realised association signs $\delta_B,\delta_W$ of 2c M1 (rows in $n$);
  * CR as a union of two LPs;
  * slack dials, as in 2c M2(d).
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

### 3.2 Step 2: M-AUDIT, choosing which precinct to validate

**What it does (plain language).** Archives sometimes let you recover the full table for a few precincts, from poll books linked to service records or from a validation study. That is expensive. M-AUDIT tells you **in advance** how much each precinct's table would narrow your answer, and it guarantees that narrowing whatever the table turns out to show. You can then spend the budget where it counts. Big precincts dominated by the group you study are often a poor choice: their arithmetic already pins the answer.

**Information or assumption used.** The margins, plus the cost of recovering each precinct's table. The planner assumes the recovered table is complete and correct for the same voters.

**The question it answers.** "Which precincts should I pay to validate, and what will I get?"

**Worked example.** The same two precincts:
* P1 can move the veterans' yes count by 20 votes, out of a total range of 60.
* P2 can move it by 40.

Recovering P2's table removes exactly 40 of the 60 votes of uncertainty (the rate interval shrinks from 0.75 wide to 0.25 wide), whatever P2's true table is. Recovering P1 removes 20. With equal costs, audit P2 (N8a: in 200 random designs the gain equals the precinct's own range exactly).

There is one complication. If a published total links the precincts, the gains stop adding up. Take a published county total of 45 veteran yes votes and the target "P1's rate". Auditing P2 then pins P1 completely: it removes all 20 votes of P1's range, though P2 is not in the target. Without the total it would remove nothing (N8c).

**Why now.**
* It is a thin layer on Step 1.
* It changes what researchers do in the archive.
* It sets up the validation work in Step 3.
* Its main fact is elementary.
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

### 3.3 Step 3: known tables and richer margins (M-DATA: exact uses first, then transfer)

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
* **Transfer (E4).** Neighbouring electorates' tables are used to bound a held-out electorate's table under "the same hidden pattern ± δ". Without the module you have the LCE bounds. With it you have a narrower interval whose validity depends on δ. NZ shows which δ was needed historically.

**Why now.**
* The exact parts reuse Step 1 directly and need no proof.
* E4 is the repository's own top data priority (NEW_AVENUES §13.8).
* NZ has the tables.
* E4 needs a proof (rank) and the actual-vs-model allowance. Doing the exact parts first keeps that cost separate.
* **Cost:** S for E5, M for E4 (proof plus the leave-one-out benchmark).

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

### 3.3b Step 3b: M-CVR, ballot patterns across contests (order flexible; small)

**What it does (plain language).** Many US jurisdictions publish cast-vote records: for each ballot, the choice in every contest. Knowing that no ballot split its ticket tells you that **every group voted the same way in both contests**, even though the records still cannot say how any group voted in either one. M-CVR bounds group-specific *differences* across contests: split-ticket voting, roll-off and down-ballot defection.

**Worked example.** One precinct of 100 voters, 50 in each of groups A and B, two yes/no contests.
* If 50 ballots are "no–no" and 50 are "yes–yes", group A's support in contest 1 can be anything from 0% to 100%. But its support in contest 1 minus its support in contest 2 is **exactly 0**.
* With 4 ballots "yes–no", the difference is in **[0, 0.08]** (N11).
* Without the module, the package would treat each contest separately and say nothing about the difference.

**Why here.**
* With Step 1 in place it is a change of columns and nothing else.
* It is a strict gain, and it leaves the proved negative result (C1: joint margins do not help a single contest's cell) intact.
* Run it earlier if you have a split-ticket application (decision D4); otherwise it can wait.
* **Cost:** S; data are public (Kuriwaki et al. 2024).

**Technical notes.**
* The columns are ballot patterns. Targets are linear contrasts of pattern cells (for example $n_{A,10}-n_{A,01}$).
* It is a transportation LP, so it has integer vertices.
* A partially recorded precinct gets an "unknown pattern" block.
* Proofs: none beyond R19.
* Interface: `ballotPatterns(d, contests = c("pres","sen"))` builds the columns, then `eiCounts(..., targets = list(diff_contests("A", "pres:yes", "sen:yes")))`.
* Tests: N11, brute force, and reconciliation of CVR totals with the official returns.

---

### 3.4 Step 4: M-CONTRAST, finite comparisons between observed kinds of places (M-LOC folded in)

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
* It comes after Steps 1–3 because those add information, which is safer than adding beliefs.
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
* It is not needed before those beliefs exist.
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

### 3.6 Step 6: DiD on the count layer, descriptive first, then causal

**What it does (plain language).** For the same men voting in two elections, the descriptive part bounds **how much more** the veterans' support changed than the others' support. Arithmetic does this in each township, optionally with signs on how changes relate to composition. A second, clearly separated step asks the causal question: what would veterans have done had they not served? That needs an explicit belief about how far the two groups' trends could have differed without the war. The package reports the descriptive range first. It reports the causal range only with that belief written down as a dial.

**Information or assumption used.** Descriptive part: the two elections' returns for the same population, plus optional signs on changes. Causal part: a bound on the gap in untreated trends.

**The question it answers.** First, "how much more did veterans' support change than other men's?". Then, separately, "how much of that was caused by service?".

**Worked example** (arithmetic only).
* Suppose the descriptive layer finds that veterans' support rose between 5 and 25 points more than others'.
* If, absent the war, the two groups' trends could have differed by up to 10 points either way, the effect of service is between −5 and +35 points.
* If they could have differed by at most 2 points, it is between +3 and +27.

Without the module the user must assemble this by hand. PACKAGE_DESIGN §7 also labels an ordering of *changes* "MTR", which reads a descriptive fact as a causal one.

**Why here.**
* The exact rungs (Δ0, Δ4, the composition polytope) are realised-table results, so they belong on the count layer, and Step 1 makes them cheap.
* PHASE1b and PHASE2d have not been reviewed: Batch C is the gate for the causal layer and for the FD rungs.
* The Wisconsin application is the owner's. Doing it after Steps 1–3 lets it use audits and stratified margins too.
* **Cost:** M for the descriptive layer (after Step 1); M for the causal layer (after Batch C).

**Technical notes.**
* **Descriptive.**
  * Δ0: product of the two periods' unit sets.
  * Δ4: 2d Prop. Δ4, rows on changes.
  * Composition: 2d §7, interval bases in count coordinates.
  * Party: 2d P2/P3.
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

**Why it comes last of the numbered steps.**
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

### 3.8 Conditional modules (when the data exist)

* **M-MISCLASS, uncertain group labels.**
  * *Plain language*: when "veteran" comes from record linkage, some labels are wrong. Bounding the number of wrong labels with an audit gives honest (wider) bounds.
  * *Example*: P2 records 60 veterans. If up to 2 men could be mislabelled, the veterans' yes rate widens from [0.167, 0.833] to [0.138, 0.862] (N14).
  * *Technical*: counts by true group × recorded label × outcome, with audit bounds on false and missed matches. A ratio with an unknown denominator is linear-fractional, so use profiling over integer group totals. Build it on Step 1. Molinari (2008) supplies the framework.
  * *When*: once a linkage audit is in hand (decision D4).
* **R16/R15 (M-SIGN out of experimental).** Unchanged. They are on the theory side, they run in parallel, and they are not a build step.

---

## Part 4. Not now, and why

| Item | Why not now | What would change that |
|---|---|---|
| **M-SHARP series estimator** (NEXT_BUILD item 4, R9 as stated) | The oracle gain is median 0–1%, and at most 5% (`sharp_gain_study*.py`). Plug-in versions are biased inward. "Tightens every rung" is not a sensible release criterion | Step 7's diagnostic finds data with large kink mass, or the NZ benchmark shows gains |
| **M-LOC as its own module** | It is the "everywhere" option of M-CONTRAST under LCE (N12c). C1 (R1) is needed only for the implication to $\delta_W$ | none; build it with Step 4 |
| **2×C joint/ordinal signs (R6), R×2 spanning trees (R5)** | Evidence of user demand is thin, and R5(c) and R6 are open or only checked. These are LCE-only refinements | an application with ordered outcomes or ≥3 groups that asks for them |
| **B1(M), κ, CAR-δ, M-DISP, single index, Bernstein, tensor contests** | Each needs either a new modelling commitment or rich data, and none relieves the current bottlenecks (credible information, honest inference with sparse designs, validation) | a dataset with a second-scale composition (B1) or covariates with overlap (CAR) |
| **Causal DiD before the descriptive layer and Batch C** | PHASE1b and PHASE2d are unreviewed; the causal labels need explicit trend assumptions | Batch C returned |
| **Switching the default interval to projection** | **This would reverse your decision.** I do not recommend it now. 0.3.1 already prints projection beside CLR with a per-fit status, and Step 0 adds the R2-regime line | Conjecture P3.5(b3) fails in a design users meet, or a referee insists |
| **General R×C full-law estimation** | The cost is high and the gain unproven (above) | Step 7 shows gains |
| **Covariate adjustment beyond stratified margins** | Low priority by your decision. Step 3 delivers the assumption-free version (stratified DD) | your call |
| **Superpopulation target, clustering** | Excluded by your decisions | — |

---

## Part 5. Where I agree and disagree with the external review

**Agree:**
1. **Lead with the quantity the question asks about.** Step 1 and decision D1. Part of this is done already (`ci_actual` ∩ DD, `proj_actual` on P3.5(c1)).
2. **A count-table interface beside the LCE solver, keeping the skew-symmetric LP for LCE.** This is Step 1. I would also adopt the review's request that each module declare which variables it restricts (actual cells or model means). The `meta$object` field already exists, so the remaining work is to make combination rules check it.
3. **The subset-inequality identity is right.** Verified (N1–N3c), with a short submodularity proof (R24). Averaged cell limits suffice for C ≤ 3.
4. **Targeted validation.** Step 2. I add the concavity fact (N8d), which keeps the coupled case cheap.
5. **Finite-composition contrasts instead of local derivatives.** Step 4. I go further: under LCE the contrast rows *are* the P-S objects at the midpoint (N12c), so M-LOC need not be a separate module.
6. **Integrality**: an integer check where rows break total unimodularity. N10c gives a concrete LP-feasible, integer-infeasible system.
7. **Descriptive DiD before causal**, and the "MTR" mislabel.
8. **Exception-mass frontiers** and **group misclassification** as useful later modules.
9. **Ballot patterns** as a strict gain from richer data, with C1 left intact.

**Disagree, or agree with qualifications:**
1. **On M-SHARP's place.** The review puts M-SHARP second, before M-LOC and after the count layer. I put it last of the numbered steps and narrow it to a diagnostic. The review's argument rests on population information. My oracle study finds that information small in LCE designs (≤5% width; median about 0). The review itself notes that finite-sample bands can erase it. I agree with the review against the repository on narrowing ("remove 'tightens every rung'") and on conditioning on size. I disagree with the repository's NEXT_BUILD item 4 on priority.
2. **On the default interval.** The review would make the coverage claim pick the default, with projection where proved. You decided to keep CLR as the default, labelled algorithmic, with projection beside it. 0.3.1 implements that transparently. I would not reverse it now. What remains of the review's point is cheap and does not change the default: the R2-regime line (Step 0). Revisit if Conjecture P3.5(b3) fails in a realistic design.
3. **On E4 before E5.** The review's build order promotes E5 and E4 together. I would ship the exact uses (E5, published totals, strata) in the count layer first. They need no allowance and no rank proof. Transfer (E4) is an assumption and should follow with its δ dial and NZ evidence. NEW_AVENUES §5.2 E2 already found weak persistence of $K$ in NZ.
4. **On "short proved chains" replacing the ladder.** Agreed for `ladder()`'s printout (Step 0). The plausibility tiers can stay in the documents as an *elicitation aid*, labelled as such, because users find them helpful.
5. **On the review's ordering of M-LOC after the validation work.** I agree it is late. I would not build it separately at all (point 5 under Agree).

**Where I disagree with the repository's own plan:**
* PACKAGE_DESIGN §6 orders LCE refinements (M-LOC, 2×C, R×2) before data modules on grounds of proof readiness. The count-layer modules are at least as proof-ready: standard LP facts plus 2c M2, which is proved. They also answer the users' first question. So that rationale no longer favours the LCE refinements.
* NEXT_BUILD item 4 says M-SHARP "tightens every rung with no new assumption". True in principle, but negligible in the designs tested.

---

## Part 6. Decisions you need to make before building starts

* **D1. Print order.** Should `print.eiRxC()` lead with the actual-rate block (DD, then any count-layer restrictions, then `ci_actual`/`proj_actual`) and put the model-implied block second? *My recommendation: yes, as the default, with `target = "model"` restoring the current order.* This does not touch your CLR decision.
* **D2. One function or two.** Either a new `eiCounts()` that `eiBounds()` calls for its actual-rate block, or a single `eiBounds(target = )` with two engines. *Recommendation: a separate `eiCounts()` (clean semantics: exact, no inference) that `eiBounds()` embeds.*
* **D3. Re-ordering.** Do you accept moving M-SHARP (to a diagnostic) and M-LOC (into M-CONTRAST) below the count-layer modules? This changes `NEXT_BUILD.md` item 4 and PACKAGE_DESIGN §6 steps 1b–4.
* **D4. The first application to drive the build.**
  * Wisconsin/Iowa with poll-book audits: Steps 1, 2, 6 and M-MISCLASS.
  * NZ tables: Step 3, E4.
  * A US CVR split-ticket application: Step 3b early.
* **D5. Integers.** Should the default be `integer = "auto"` (MILP only when rows break total unimodularity, otherwise the LP with a certificate)? Or always LP with a flag? *Recommendation: "auto", with a time limit that falls back to the LP bound plus a labelled gap.*
* **D6 (minor).** Whether to send batch D (PACKAGE_DESIGN/NEW_AVENUES) with this recommendation attached, so that the outside reviewer sees the revised order rather than the old one.

---

## Part 7. What I verified, and what I take on trust

**Verified numerically here** (`verify_next_steps.py`, 24/24 PASS):

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

**Population studies** (`sharp_gain_study.py`, `sharp_gain_study_2xC.py`). Pooled LCE, oracle law at replicated composition types, large precincts, Beta or Dirichlet heterogeneity. The result is that the full law narrows the model-implied group rates by little.

*Limits of these studies*:
* only 2-row designs;
* no binomial noise;
* no size channel (all precincts large);
* linear truths.

The size channel (N4) is the main way a larger gain could appear. Step 7's diagnostic should look for it in NZ.

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
| R20 | Audit: separability; concavity of the remaining width; a joint audit by MILP | (a) checked N8a; (b) checked N8d (standard LP sensitivity) | Step 2 |
| R21 | Rank of the stacked table map in terms of the composition second moments (E4) | checked (NEW_AVENUES §13.1); proof open | Step 3 (E4) |
| R22 | Finite-composition identity; LCE midpoint form; realised bin rows | checked N12a–c; proofs are algebra | Step 4 |
| R23 | Exception-mass frontier: monotone, MILP exact | elementary | Step 5 |
| R24 | Subset identity (2×C, R×2) = Aumann expectation | proof sketch here; checked N1–N3 | Step 7 |
| R9′ | (replaces R9) outer validity of the instrument-averaged M-SHARP procedure, plus a gain study | open | Step 7 (only if an estimator is built) |
