# The package plan in plain language

*Companion to `PACKAGE_DESIGN.md`. No formulas. Bracketed tags point to the formal results.*

## What we are building
A toolbox where every assumption is a separate, labelled piece. You start from what the vote totals alone allow (Duncan–Davis bounds) and add pieces one at a time. Each piece says:
* what extra information or belief it uses;
* whether the data could ever show it is wrong;
* how much it narrowed the answer;
* whether it was already implied by what you had.

## What the outside (GPT) conversation got right
* **Turning the straight-line (JKST) bounds into ranges for the two Monotone-EI patterns adds no information.** It is the same uncertainty described differently [G1].
* **Adding a direction for either pattern on top of JKST usually helps.** In simulations it narrowed the range in roughly 80–90% of cases, typically to about two-thirds of its JKST width [G2, V2].
* **When it does not help, the data already implied that direction.** The package will show this before you assume anything.
* **One caution.** "The two assumptions conflict" counts as evidence only against the *pair*. A wrong direction usually does *not* produce a conflict: it just gives a narrow, wrong range [G3].

## Two new findings
1. **The "sandwich" [P-S].** Take two groups, A and B, at any town make-up. The overall vote pattern tells us how the total yes-rate changes when an A resident replaces a B resident. That change is the sum of two parts:
   * the gap between A and B in that town;
   * how everyone's rates shift with the make-up.

   If you assume both parts are non-negative, each is squeezed between zero and the observed change. This holds for any number of groups. If the observed change is negative somewhere, the data reject the pair of assumptions.
2. **Why more groups help [P-R].** With two groups there is one unknown, and every assumption pushes on the same unknown. With more groups there are more unknowns. Each pairwise comparison (A vs B) reaches only some of them, but a connected set of comparisons (A vs B, B vs C, …) reaches all of them. In 3-group simulations these comparisons narrowed about two-thirds of the group-rate ranges [V5].

## The two ladders
* **Assumptions** (from none to strong), in three tracks that are reported side by side because none contains the others:
  * *direction* (which group is more supportive; which way rates move);
  * *shape* (straight lines, curves, constancy);
  * *size* (how large contextual effects or noise may be).
* **Build order:** first the pieces whose mathematics is already proved and that use the existing machinery (direction rows, the sandwich), then the data-combination pieces (known tables), then the pieces needing new data or new theory. No piece goes into the package until its proofs pass an independent check.

## Difference-in-differences
Every piece has a version for *changes* between two elections. The sandwich applied to changes gives a new closed-form DiD range from two direction beliefs:
* veterans' support did not move less than non-veterans' in any town;
* replacing a non-veteran with a veteran did not lower the town's change.

It needs no straight-line assumption about either election's levels, only about the changes [DiD-S].
