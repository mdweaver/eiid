# Monotone ecological inference: a plain-language version of Phase 2c

*Companion to `PHASE2c_monotone.md`. No formulas. Bracketed references point to the formal results, e.g. [M3] is Proposition M3. Updated after the Batch B referee reports (Nov 2026) to match the corrected Phase 2c.*

---

## The problem in one paragraph

For each town we know how many men are veterans and how many voted yes. We want the yes-rate *among veterans* and *among everyone else*. The town totals cannot tell us who cast the yes votes. Duncan and Davis's method of bounds gives the range that is arithmetically possible. It is always valid but usually wide. Every other method narrows that range by assuming something about behaviour.

## What "monotone EI" assumes

The paper by Elzayn, Goldin, Guage, Ho and Morton narrows the range with assumptions about the **direction** of two patterns, never their size:

1. **The "between" pattern:** within the same town, does one group vote yes more often than the other, on average across towns?
2. **The "within" (contextual) pattern:** does each group's yes-rate rise or fall with the town's make-up? Example: are veterans *and* non-veterans more pro-suffrage in towns with many veterans?

A third option, **contextual reinforcement**, assumes only that the two patterns point the same way, without saying which way.

These assumptions connect to two familiar estimates:
- Goodman's regression is right exactly when the within pattern is zero.
- The neighbourhood model (each group gets its towns' average yes-rate, weighted by its numbers) is right exactly when the between pattern is zero.

Knowing the *direction* of a pattern tells you which side of the truth each estimate falls on. That is what cuts the range.

## What Phase 2c shows

**1. The bounds are sharp for the actual population, treating rates as continuous [M2].** Combining the method of bounds with any set of direction assumptions is a single linear program (a standard optimisation). The answer contains every value that some town-by-town table, consistent with the vote totals and the assumptions, could produce, and nothing else, **if a table may split people into fractions**. No sampling theory is needed: these are statements about the towns you have.

*With whole people the answer can differ [Remark M2.1].* In a town of four voters, one vote cannot be split between two groups. Two small examples from the referees show that the continuous calculation can accept "contextual reinforcement" when no whole-ballot table satisfies it, and can move an endpoint (0 continuously, 2/15 with whole ballots). For two groups and margins only, the difference is always less than one "step" of one person; at realistic sizes that step is about six parts in a million, so it matters only for tiny tables or near-ties. Beyond that case it has not been bounded.

**2. What the town totals alone can reject depends on what else you assume [M3].** With two groups, continuous rates, and no other restrictions:
- an assumed **direction of the within pattern** can turn out impossible. When that happens, the bound comes out empty and the assumption is refuted.
- an assumed **direction of the between pattern**, and **contextual reinforcement**, can never be rejected by town totals *in this setting*. That is not a general guarantee:
  - added to the straight-line model (JKST), reinforcement can be rejected even when the straight-line model alone fits;
  - with whole-person counts, reinforcement can be impossible outright (point 1);
  - known cells, structural zeros or a unit-level ordering can also rule it out.

  Where it survives, its credibility rests on substantive argument, or on individual-level data for a subsample, as the original paper used.

Under contextual reinforcement, the sign of the group difference equals the sign of Goodman's slope, provided at least one town is mixed. If every town is all one group, the difference is known exactly and both directions give the same answer.

**3. Imposing directions on several outcomes at once helps [M4].** With yes / no / abstain columns, assuming the contextual direction for two columns at once gave a strictly narrower range in about half of the test cases. That is why the package will solve everything jointly rather than one column at a time.

**4. With more than two groups there are two ways to define "the within pattern" [§1].**
- **One group versus everyone else** (the paper's version, our default). Simple, but "everyone else" is a different mix in each town, and the predictor for both sides is the focal group's share.
- **Within each original group, pooled.** This conditions on the original groups, but it is still **one sum** over the groups, not a direction for each group: some groups can go the other way. It enters Goodman's multiple-regression errors, but signing it does not by itself sign the error of a particular group's Goodman estimate. A direction for each named group needs separate restrictions.

The two can disagree: in a two-town example where every group's rate rises by two points, the one-vs-rest version is negative and the pooled version positive.

The "between" pattern is an average over mixed towns only: towns of one group get no weight, and one large mixed town can dominate. It is not an ordering in every town or against every rival.

**5. These directions fit our existing straight-line models [M5].** If you already assume straight-line contextual effects (JKST), each direction assumption is one extra line in the same calculation, also with riding (district) effects. Our existing "ordering" option (group A at least as supportive as group B everywhere, in every riding) implies the between direction.

*One shortcut holds only without riding effects.* With one set of lines for the whole map, "both groups' rates rise with the share" implies a positive within pattern. With riding-specific intercepts it does not: in a six-town example with two ridings, both groups' rates rise within each riding, but the riding with more of the focal group has lower rates overall, and the pooled within pattern is negative. With riding effects the sign is one pooled restriction across ridings, and the sharpness results for that case are not yet proved.

**6. A new rung with no straight lines: "no reversals" (N2) [§4].**
- *What it is:* each group's average yes-rate moves in one direction as the town's make-up changes, in any shape: flat, curved, or step-like. It is defined on the town compositions actually observed.
- *How it relates:* it is weaker than "straight lines with the same signs", so its range always contains that one [N2.2(c)]. With two groups, it implies the within direction above [N2.2(d)]. With three or more groups that link breaks, and "in one direction" must say which other share adjusts [N2.2(e)].
- *What it is not:* not sharp. Using a coarser grid than the observed compositions loses information; and town size can carry information that a composition-only calculation ignores.
- *The local version is dropped.* The paper's town-by-town results use the slope of the overall town pattern at each composition. With a finite set of observed compositions that slope is not pinned down by the data: two curves can agree at every observed town and have different slopes, and using the wrong one can exclude the truth [N2.3]. Even with smooth curves, "rising at every composition" and "rising on average" do not imply each other.
- *If wrong:* a reversal, such as backlash where a group is a large minority, can push the truth outside the range.
- *Check:* an empty range refutes it, together with whatever else is assumed.

## How the pieces line up (most to least plausible)

1. **Method of bounds:** always valid.
2. **Plus a between direction, or contextual reinforcement:** not rejectable from totals alone with continuous rates, so judge it on substance; rejectable once other assumptions or whole-person counts enter.
3. **Plus a within direction:** rejectable.
4. **"No reversals" (N2):** needs an estimated town curve.
5. **Straight lines with signed slopes:** the narrowest, and the most assumed.

Straight-line models and these direction rungs are different kinds of assumption, so neither contains the other. Report them side by side.

## What still has to be done (Phase 3)

- A "slack" dial for each direction assumption, and how far the pattern could go the "wrong" way before a conclusion fails. With several dials the answer is a frontier, not one number; a single number needs a stated path (for example, all slacks equal). The M-SIGN module already has a per-restriction slack dial for the straight-line version.
- For N2: estimating the town curve with a polynomial whose degree is fixed in advance (choosing it from the same data invalidates the usual intervals), and, with many towns, grouping them into bins weighted by people, not by towns.
