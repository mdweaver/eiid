# Difference-in-differences with ecological data: a plain-language version of Phase 2d

*Companion to `PHASE2d_did.md`. No formulas. Bracketed references point to the formal results. The running example is veterans and the Wisconsin 1857/1865 suffrage referenda.*

---

## The question

Did veterans' support for Black suffrage change more between 1857 and 1865 than non-veterans' support did? That is the **difference-in-differences (DiD)**: the veterans' change minus the non-veterans' change.

The difficulty is that we never see how veterans voted, in either year. We see only each town's veteran share and its yes votes. So both years are ecological-inference problems, and the DiD inherits the uncertainty of both.

## 1. With no link between the years [D0]

If we assume nothing connecting 1857 to 1865, the DiD range is simply the 1865 range for the veteran–non-veteran gap minus the 1857 range. The widths add up.

- **Proved here:** this is the best possible bound under that information. Nothing tighter is justified without a link between the years. Phase 1 only conjectured this.
- **Your appendix §C.3 bounds** are exactly this calculation, with JKST straight-line bounds in each year.

## 2. Only the changes need to follow straight lines [FD]

**The assumption.** Each town may have its own persistent political culture, of any shape. The *change* from 1857 to 1865 must vary in a straight line with the town's veteran share. This is much weaker than requiring straight lines in each year.

**One extra rule turns the DiD into a single number [FD2].** Each rule says something about how the war's effect varied across towns:

| Rule | Plain meaning | Main way it could fail |
|---|---|---|
| **R1** | Non-veterans changed by the same amount whether many or few veterans lived nearby: no *extra* persuasion of neighbours where veterans were numerous. | Veterans persuaded or provoked neighbours in proportion to their numbers. |
| **R2** | A veteran changed by the same amount whether he came home to many or few fellow veterans. | Veterans' groups at home reinforced each other. |
| **R4** | On average, the two groups' contextual changes cancel out. This is what your regression of the change in yes share on veteran share assumes. | Both groups moved the same way with veteran share. |

What follows from the rules:
- **R4's answer always lies between R1's and R2's.** If R1 and R2 give similar numbers, the choice of rule hardly matters.
- **R1 and R2 together** is the old "contextual effects unchanged" assumption. It is testable: the curvature of the town pattern must be the same in both years.

**How much would it take to overturn the result? [FD3]** Allow R1 to fail by up to some amount. The **breakdown value** is the smallest extra spillover onto non-veterans that would erase the effect. This is your appendix Figure C1, made formal.

**With no rule at all [FD4],** straight-line changes still give a bound, even though each town's levels are arbitrary.
- **What it uses:**
  1. The regression of the change in yes share on veteran share. It pins down everything about the change except **one shared number per outcome**: how the curvature of the change splits between veterans and non-veterans.
  2. Both years' average yes shares at each observed veteran share, which give each year's possible rates there.
- **How it works:** a value of that shared number is kept only if, at every composition, the implied change for each group can be written as (a possible 1865 rate) minus (a possible 1857 rate).
- **Where it tightens:** compositions where a group's rate is near 0 or 1 in either year (for example, towns with almost no 1857 yes votes). There the change cannot go beyond what the floor or ceiling allows, and because the number is shared, that limit applies to all towns.
- **Compared with "no link" (section 1):** it is not guaranteed to be narrower.
  - It gains from the shared number, which links all compositions.
  - It loses because it works with averages of towns at the same composition rather than each town's own range.
  - Report both.
- **Compared with your appendix bounds:** it is always at least as wide as the straight-lines-in-levels bounds, because it assumes less.

## 3. No straight lines at all: monotone EI on changes [D4]

Apply the Monotone-EI idea to the *changes* themselves. This needs only direction assumptions:

- **If spillovers ran toward suffrage** (towns with more veterans moved more, for both groups), then your Goodman-on-changes estimate (0.265) is an **upper bound** on the DiD.
- **If there was backlash** (non-veterans reacting against veterans), the same estimate is a **lower bound**.
- **If the spillover and the veteran–non-veteran gap point the same way** (contextual reinforcement), the DiD has the **same sign** as Goodman on changes. This assumption can never be rejected by the data.
- **The direction of the within-town veteran gap** is not a usable assumption, because it amounts to assuming the answer.

## 4. Assuming the contextual effects are identical in both years [§5]

This is the "periods as groups" model. It assumes R1 and R2 *and* straight lines in each year. In exchange, it also narrows the bounds on each year's levels.

## 5. When the electorate changes [D6, C1, C2]

**Population counts change.** Even R1 then gives a short **interval**, not a single number. The interval is narrow when population growth was unrelated to the veteran share.

**Aging-in.** Using the 1865 eligible men as the base for 1857 counts young men who could not yet vote in 1857 as "did not vote yes". Veterans were younger, so this can make their change look larger.
- How much it matters depends on pre-war support. If young men vote like older men of their group in 1865, the extra amount is roughly the veterans' newly-of-age share times their 1857 yes-rate, minus the same for non-veterans.
- Pre-war support was near zero in most towns, so the extra amount should be small. The new tools can also remove it entirely.

**Deaths and men who did not return.** The 1857 vote includes men who later died. Your earlier assumption is that they voted like the survivors. A dial lets them differ by a stated amount.

**Exact bounds for the "same men in both years" DiD [C1].** Splitting each town's men into stayers, newly eligible men and departures gives a calculation that remains exact.
- The census counts from 1860 and 1870 bracket each town's eligible men.
- The number of veterans who were too young in 1857 is known only within a range, unless ages at enlistment are available.

## 6. From DiD to the effect of service [CL]

Your appendix §C.4 assumptions bound what veterans would have done without serving:
- **(a)** service did not lower support;
- **(b)** without the war they would not have become less supportive than in 1857;
- **(c)** without service they would have been at least as supportive as the men who stayed home.

The effect is then between **0** and an upper limit. Your script's upper limit is exactly right: it is the best possible bound, not a conservative one.

The three assumptions are **jointly testable**. If no configuration consistent with the data satisfies all three, they are rejected. Assumption (c) does most of the work.

Alternatively, **parallel trends** (future soldiers would have changed like everyone else) makes the effect equal to the DiD.

## 7. Party [P2, P3]

- **"Democrats never voted yes" does not, by itself, narrow the veteran bounds.** Yes votes never exceed Republican votes plus men who skipped the governor's race, so the restriction rules out nothing.
- **Adding "before the war, within each party, future soldiers voted like their neighbours"** does narrow them. Each town's pre-war veteran yes-rate then depends only on the veterans' party mix, which is itself bounded by EI. This is often strictly narrower, and it targets the 1857 end, which is the widest part of the DiD.
- **Why the assumption fits:** it is the formal version of your argument that selection into the army ran through party.
- **How it could fail:** within the Republican party, the most anti-slavery men enlisted more.

## How the DiD tools line up (most to least plausible)

1. **No link between years:** always valid.
2. **Monotone EI on changes:** directions only.
3. **Straight-line changes:** towns may keep their own cultures.
4. **Plus a bounded spillover:** a dial, with a breakdown value.
5. **Plus R1, R2 or R4:** one number, or a short interval if counts changed.
6. **Identical contextual effects in both years:** also tightens each year's levels.

The **composition corrections** (stayers, newly eligible men, departures) and the **party assumption** can be added to any rung. The **causal step** (assumptions (a)–(c), or parallel trends) is a separate layer on top.

## What still has to be done (Phase 3)

- Uncertainty statements for the rungs that use regression (straight-line changes, R1/R2/R4), including a confidence interval for the breakdown value.
- Slack dials for the exact rungs.
- Treating the census-based population ranges as given.
