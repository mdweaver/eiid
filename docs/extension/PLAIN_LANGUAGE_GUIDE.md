# A plain-language guide to the bounds

*A companion to `PHASE2a_rxc.md` and `PHASE2b_districts.md`, written for using the methods rather than proving them. Each entry says what the method is, what you are assuming, what goes wrong if the assumption is false, and what you can check. Formal statements are referenced in brackets, e.g. [2b Thm 2].*

*The examples come from voter transitions (BC) and from other ecological inference problems (literacy by race in counties, veterans and suffrage). Nothing here is specific to elections.*

---

## 1. The one thing the data cannot tell you

You observe, for each unit (a poll, a county), the share of each group (rows: e.g. last election's in-party voters, out-party voters, non-voters) and the share of each outcome (columns: e.g. this election's choices). You want the group-by-outcome rates, e.g. "what share of last time's in-party voters stayed loyal".

Plotting outcomes against composition reveals **patterns**, e.g. "P does worse in polls with many out-party voters". What the data can never reveal is **whose behaviour produces a pattern**. It could be in-party voters in those polls defecting more. It could be out-party voters there crossing over less. It could be any mix of the two.

Every method below is a different answer to "how much are we willing to assume about who is behind the patterns?" Assume little and you get a range (a *bound*). Assume enough and you get a single number (a *point*). The more you assume, the narrower the answer, and the more it depends on the assumption being right.

In the formal documents, the unknown "who is behind it" numbers are called $K$. A 2×2 table has 1, a 2×3 table has 2, and a 3×3 table has 6 [2a Lemma 1].

---

## 2. Foundations shared by every method

These apply to every rung below. If one fails, *all* the bounds can be off.

| Assumption | Plain meaning | If it fails | How to check / fix |
|---|---|---|---|
| **Same people in both margins** (A1) | The row shares and column shares describe the same population. For transitions, the same electorate in both elections; for veterans, the same eligible men. | The arithmetic that links rows to columns is wrong, so the bounds can miss the truth. Turnover between elections, people moving, deaths, new voters and boundary changes all cause this. | Use one base for both margins (e.g. registered voters). Add rows or columns for "entered" and "left". Use an abstain column instead of dropping non-voters. |
| **Rates don't depend on unit size**, once composition is accounted for (A3) | Big polls and small polls with the same mix behave alike. | The province-wide rate is biased, because big units carry more weight. | Add (log) size as a covariate (2b-style grouping or the appendix of the original paper). |
| **Enough variation in composition** (A5) | Units differ enough in their group shares to fit the patterns. | Bounds are wide and unstable. This was the original paper's problem with gender, which is about 50/50 everywhere. | Look at the spread of compositions and the condition number of the regression. |

---

## 3. The rungs: from no assumption to full constancy

Listed from **most** to **least** plausible. Each rung's bound sits inside the one above it [2b Thm 2(c)]. That nesting is the point of the ladder: you can report all of them and let readers see how much each extra assumption buys.

### Rung 0: Duncan–Davis (no model)

* **What it is:** the range of rates that is arithmetically possible, unit by unit, given only the margins. Averaged over units, it gives the overall bound.
* **Key assumption:** none beyond the foundations.
* **If wrong:** it cannot be wrong unless a foundation fails.
* **Note:** grouping units (ridings, states) adds *nothing* at this rung. Any gain from groups comes from an assumption about them [2b Prop D0].
* **Use:** always report it as the baseline that everyone can agree on.

### Rung 1: district-specific linear contextual effects (M2)

* **What it is:** within each group (riding, state), rates may change with composition, but in a straight line. Each group has its own lines, and nothing is shared across groups.
* **Key assumption:** linearity within each group. For example, in-party loyalty rises steadily, not in jumps, as a poll's in-party share rises.
* **If wrong:** with threshold or tipping effects (a party machine that mobilises only above 50%), or U-shapes (competitive polls mobilising both sides), the bound can miss the truth. More data does not fix this. Some such violations are invisible in the data [2a Prop 9].
* **Check:** add cubic terms to the regression; an empty bound means a violation.
* **Note:** the group bounds simply add up, so the widths add too [2b Prop D5]. This rung is safe but often not much tighter than Rung 0. It needs enough units in each group.

### Rung 1b: the pooled model (M0); not on the ladder, reported alongside

* **What it is:** the original paper's method (and the package), extended to R×C. One set of straight-line contextual effects for *all* units, ignoring groups. Its bound is the "outer" bound, computed from the regression [2a Prop 2]. A tighter "sharp" version also uses how much similar units vary [2a Thm 1]; it matters only when similar units sometimes hit an arithmetic corner [2a Prop 3].
* **Key assumption:** linearity across all units. That implies groups don't differ in baseline in a way related to composition.
* **If wrong:** when strongholds also have different baselines (e.g. a popular local candidate in the strongest ridings), the pooled fit mistakes those baseline differences for a contextual effect. The bound can miss or come out empty [2b Prop D4; P1 §7.3].
* **Check:** compare with Rungs 1–3. An empty bound tells you something is wrong.
* **Use:** when you have no groups (e.g. stage 1 with 93 ridings as the units).

### Rung 2: shared contextual effects, with bounded differences across groups (M1(δ))

* **What it is:** as Rung 3, but each group's "who is behind the pattern" numbers may differ from the shared ones by up to δ. δ is in probability units: δ = 0.05 allows the gap between, say, "in-party voters in out-party areas" and "out-party voters in in-party areas" to differ by up to 10 points across ridings.
* **Key assumption:** groups are *similar*, though not identical, in how composition matters.
* **If wrong:** the same failure as Rung 3, but softened. A larger δ protects more but narrows less.
* **Use:** as a sensitivity dial. Report how the bound moves as δ goes from 0 (Rung 3) to large (Rung 1) [2b Prop D6].

### Rung 3: shared contextual effects, group baselines (M1, "riding dummies")

* **What it is:** each group may be more or less loyal across the board (its own baseline), but **the way rates respond to composition is the same in every group**. Two groups' lines are parallel.
* **Key assumption:** no interaction between group and composition. For example, a popular candidate raises loyalty by the same number of points in strongholds and in hostile polls.
* **Why it helps:** each group's baseline is estimated from its own units. A group with an extreme baseline *and* lopsided units nearly reveals one group's rate, and because the contextual effects are shared, that information tightens the bound for every group [2b Prop D3].
* **If wrong:**
  * The effect of the group differs by composition: an incumbent matters most in competitive polls.
  * Ceiling or floor effects: a +5-point shift is impossible where loyalty is already 97%. This is most likely in the most lopsided groups, which are exactly the ones that tighten the bound.
  * The bound can then exclude the truth and still look narrow.
* **Check:** one part is testable. The curvature terms of the regression must line up across groups [2b Lemma D1(b)], and an empty bound refutes the model. The shared "who is behind it" part cannot be tested.
* **Caution:** with few units per group (about 10 polls per riding), the simple plug-in version is too narrow. It covered the truth only 72% of the time in simulation; Phase 3 provides the corrected version.
* **Also covers time:** with periods as the "groups", this is the assumption "contextual effects unchanged between elections". It gives point-identified changes for difference-in-differences designs [2b Remark 6.2].

### Rung 4: one group has no contextual effect (L4(r))

* **What it is:** one chosen group, r, behaves the same way in every kind of unit. Its rates don't depend on composition. Every pattern that involves group r is then attributed to the *other* group, which removes those unknowns.
* **Example (2×2):** turnout falls as the Black share rises. If white turnout is assumed the same everywhere, the whole decline is due to Black residents voting less, and the answer is a single number. In 3×3, assuming previous non-voters are unaffected by the poll's mix pins 4 of the 6 unknowns. What is left is only the in-party versus out-party split.
* **Key assumption:** that particular group is insensitive to context.
* **If wrong:** the patterns that really belong to group r are given to the others, which biases them. With non-voters this is a real risk, because mobilisation often depends on local campaign intensity [2b Prop D7(b)].
* **Check:** largely untestable, except through an empty bound. Which group you pick changes the answer, so report each choice as a sensitivity check.
* **History:** the 2×2 version is the old suggestion "set the white contextual slope to zero" (Wakefield's w₁ = 0, discussed in the original paper).

### Rung 5: constancy within each group (L5)

* **What it is:** within a group, rates don't depend on composition at all. Groups may differ in level. This is Goodman's regression run separately inside each group, and it gives one number [2b Prop D7(a)]. It is Cho & Manski's instrumental-variable assumption with the group as the cell.
* **Key assumption:** no contextual effects for *any* group within a group.
* **If wrong:** a confident, wrong point. In the 3×3 example it gave 0.870 against a truth of 0.836.
* **Check:** partly testable: the within-group fit must have a specific structure [2b Lemma D1(c)].
* **Use:** report it as a diagnostic beside the bounds. If it sits outside a weaker rung's bound, constancy is contradicted by the data.

### Rung 6: Goodman (constancy everywhere)

* **What it is:** each group's rates are the same in every unit.
* **Key assumption:** no contextual effects and no group differences.
* **If wrong:** the classic failure: estimates outside the possible range.
* **Use:** a historical reference point only.

---

## 4. Optional extras that apply to any rung

Pinned in `ASSUMPTION_MENU.md` and not yet derived: structural zeros (a party not running), survey benchmarks, group ordering (in-party voters at least as likely to vote P as out-party voters), sign restrictions on contextual effects, stability across elections, similarity of neighbouring units, and sensitivity dials for each.

## 5. Reading results

* **A bound covers "everything consistent with the data and the assumption".** The truth can be anywhere inside it. The midpoint is not a better guess.
* **An empty bound is information:** the assumption, or a foundation, is contradicted by the data (or, in small samples, it may be noise; Phase 3 separates the two).
* **"Trusting the line to the extremes":** using the whole range 0–1 rather than the observed range narrows bounds, but it assumes the straight lines remain valid where no units exist (e.g. a pure out-party poll) [2a Prop 5]. Use the observed range unless you can defend extrapolation.
* **Suggested table:** rungs as rows, quantities of interest as columns. Each cell shows the bound, its confidence interval (Phase 3), and whether the bound was empty.

---

## 6. From bounds to confidence intervals (Phase 3, in plain terms)

* **The estimated bound is not a confidence interval.** It is computed from estimated regression coefficients, so it moves with the sample.
* **The estimated bound tends to be too narrow, not too wide.** Each end of the bound is the tightest of many candidate limits. With noisy data, the tightest-looking candidate is often tight partly by chance. So the plug-in range is biased inward, and badly so when there are many ridings with few polls each. In simulation with 93 ridings × 10 polls, it almost never contained the truth and was often empty [P3 Thm 4, §6].
* **The confidence interval fixes this.** It widens each candidate limit by its own uncertainty before taking the tightest. Noisy ridings are therefore unlikely to set the limit [P3 Thm 5]. It contained the truth in every simulation design.
* **The price.** With few units per riding, the valid interval can be about as wide as Duncan–Davis. That is honest: there is too little information per riding to exploit the riding assumption. More polls per riding, coarser groups, or the δ dial help [P3 §6].
* **Two quantities, two intervals.**
  * The *actual* rate in this electorate includes poll-level luck that can never be estimated, so its interval carries an extra allowance.
  * The *model-implied* rate at the observed compositions (the "behavioural" rate) does not, so its interval is narrower.
  * Say which one you report.
* **An empty confidence set is a real rejection.** If even the widened constraints cannot be satisfied, the assumptions are rejected at the stated level [P3 Prop 6].
* **Sensitivity in one number.** For a headline claim (e.g. "loyalty above 30%"), report the *breakdown value*: how far average rates could depart from straight lines, undetectably, before the claim fails [P3 §7].
