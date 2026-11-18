# Monotone ecological inference: a plain-language version of Phase 2c

*Companion to `PHASE2c_monotone.md`. No formulas. Bracketed references point to the formal results, e.g. [M3] is Proposition M3.*

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

**1. The bounds are exact for the actual population [M2].** Combining the method of bounds with any set of direction assumptions is a single linear program (a standard optimisation). The answer contains every value that some town-by-town table, consistent with the vote totals and the assumptions, could produce, and nothing else. No sampling theory is needed: these are statements about the towns you have.

**2. The data can reject only one of these assumptions [M3].** With two groups:
- an assumed **direction of the within pattern** can turn out impossible. When that happens, the bound comes out empty and the assumption is refuted.
- an assumed **direction of the between pattern**, and **contextual reinforcement**, can *never* be rejected by town totals. Their credibility rests entirely on substantive argument, or on individual-level data for a subsample, as the original paper used.

Under contextual reinforcement, the sign of the group difference equals the sign of Goodman's slope.

**3. Imposing directions on several outcomes at once helps [M4].** With yes / no / abstain columns, assuming the contextual direction for two columns at once gave a strictly narrower range in about half of the test cases. That is why the package will solve everything jointly rather than one column at a time.

**4. With more than two groups there are two ways to define "the within pattern" [§1].**
- **One group versus everyone else** (the paper's version, our default). Simple, but "everyone else" is a different mix in each town.
- **Each group separately.** Easier to argue about substantively, and the version that explains Goodman's multiple-regression errors.

**5. These directions fit our existing straight-line models [M5].** If you already assume straight-line contextual effects (JKST), each direction assumption is one extra line in the same calculation. Our existing "ordering" option (group A at least as supportive as group B everywhere) implies the between direction.

**6. A new rung with no straight lines: "no reversals" (N2) [§4].**
- *What it is:* each group's average yes-rate moves in one direction as the town's make-up changes, in any shape: flat, curved, or step-like.
- *How it relates:* it is weaker than "straight lines with the same signs", so its range always contains that one [N2.2(c)]. With two groups, it implies the within direction above [N2.2(d)]. With three or more groups that link breaks [N2.2(e)].
- *A local version* (the paper's town-by-town results) also works. It uses the slope of the overall town pattern at each composition.
- *If wrong:* a reversal, such as backlash where a group is a large minority, can push the truth outside the range.
- *Check:* an empty range refutes it.

## How the pieces line up (most to least plausible)

1. **Method of bounds:** always valid.
2. **Plus a between direction, or contextual reinforcement:** never rejectable, so judge it on substance.
3. **Plus a within direction:** rejectable.
4. **"No reversals" (N2):** needs an estimated town curve.
5. **Straight lines with signed slopes:** the narrowest, and the most assumed.

Straight-line models and these direction rungs are different kinds of assumption, so neither contains the other. Report them side by side.

## What still has to be done (Phase 3)

- A "slack" dial for each direction assumption, and the breakdown value: how far the pattern could go the "wrong" way before a conclusion fails.
- For N2: estimating the town curve, with an honest statement of what that adds to the uncertainty.
