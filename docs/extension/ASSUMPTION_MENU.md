# Assumption menu for partial identification in EI (pinned for later)

*A running list of restrictions that can be layered between Duncan–Davis (no assumptions) and Goodman-type constancy, in the spirit of Cho & Manski (2008): report what the data alone give, then add assumptions in decreasing order of plausibility. Nothing here is derived yet. Each item records what it restricts, whether it keeps the bounds a linear program (LP), and what to watch for. Notation follows PHASE1/PHASE2a.*

*The tools are meant to be general: any R×C EI with units nested in groups and/or observed at several times. Voter transitions (BC) and the veterans/suffrage difference-in-differences (section F) are two applications, not the scope.*

## How the pieces fit

* **Two kinds of layer.**
  * A *ladder*: nested model classes, where each rung implies the one above it, so the bounds shrink monotonically.
  * *Modifiers*: restrictions that can be added to any rung. Their bounds intersect with the rung's.
* **Plausibility order is a judgment; nesting is a mathematical fact.** Only nested layers give bounds that are guaranteed to shrink. Non-nested assumptions (e.g. pooled LCE vs. riding-specific LCE) are reported side by side, not as a ladder.
* **Where a restriction acts matters.**
  * *Unit level* (on each poll's $\beta_i$): shrinks each poll's DD polytope $S(X_i,T_i)$ and so tightens even the sharp set $\Theta_I$.
  * *Mean level* (on $E[\beta\mid X]$, i.e. on $H$): adds rows to the LP of Prop. 7.
  * *Cross-unit* (links polls, ridings or elections): adds intersections.
* **Every layer is refutable when its set is empty.** That is a free specification test at each rung (Prop. 8 generalised).

## The ladder (units nested in groups, e.g. polls in ridings; decreasing plausibility)

| Rung | Assumption | Unidentified dimension (3×3) | Status |
|---|---|---|---|
| L0 | Duncan–Davis: accounting identity only | — (no model) | done |
| L1 | Riding-specific linear contextual effects (M2) | 6 per riding | 2b |
| L1′ | Riding-specific *quadratic* contextual effects (weaker than L1, needs more polls per riding) | more; to derive | idea |
| L2 | Shared contextual effects with bounded riding deviations, $\lVert K^{(d)}-K\rVert\le\delta$ | 6 + slack | 2b/3 |
| L3 | Shared contextual effects, riding intercepts (M1, $\delta=0$) | 6 | 2b |
| L4 | **Row-specific constancy:** one group's rates do not depend on composition (e.g. non-voters) | fewer: each constant row pins entries of $K$ | idea |
| L5 | Constancy within riding (Cho–Manski with $W$ = riding) | 0 (point, per riding) | idea |
| L6 | Goodman constancy (province-wide) | 0 (point) | baseline |

L4 is a useful intermediate step. Assuming no contextual effect for one row $r$ means $H_{c,rs}$ is the same for every $s$. That fixes $K_{c,rs}$ for every $s$ in terms of identified quantities, removing $(C-1)(R-1)$ of the unknowns while leaving the other rows free. In BC, a candidate is "previous non-voters' rates do not depend on the poll's partisan mix".

## Modifiers

### A. Data-side information (most credible; acts at the unit level)

1. **Structural zeros and known cells.** A party not running in a riding, or ineligible or new registrants, who can only be in the "didn't vote last time" row. Such cells are fixed at 0, which shrinks $S(x,t)$. LP: yes.
2. **Integer counts.** Use the exact count polytope for small polls instead of shares. The gain is small except for tiny polls. LP: integer, though usually relaxable.
3. **Participation records.** If poll-level counts of "voted in both / only one / neither" exist (some jurisdictions publish participation lists), the turnout margins of the 3×3 become observed. This could tighten the table substantially. *Worth checking with Elections BC.*
4. **Survey moments (data combination).** A BC Election Study or panel estimate of province-wide loyalty or defection, with its CI, enters as an extra constraint on $\varphi(K)$. It intersects with any rung; inference via moment inequalities. LP: yes. Watch for survey misreporting, which is itself a partial-identification problem.

### B. Ordering and shape (plausible; mean or unit level)

5. **Group ordering (Manski–Pepper style monotone restrictions).** E.g. in-party voters are at least as likely to vote P as out-party voters in the same poll: $\beta_{P,P}\ge\beta_{O,P}$.
   * At the *unit* level it shrinks every poll's polytope and the sharp set. At the *mean* level it adds LP rows.
   * Very plausible for loyalty versus crossover.
6. **Sign or monotone contextual effects.** E.g. P-loyalty is non-decreasing in the P share (social context or local organisation). Linear in $H$, so LP rows. The risk is direction errors: a wrong sign excludes the truth.
7. **Bounded contextual effects (Lipschitz).** $\lvert\partial E[\beta_{rc}\mid x]/\partial x_s\rvert\le\lambda$. A sensitivity parameter; the bounds widen gracefully as $\lambda$ grows. LP: yes.
8. **Nonparametric shape-restricted contextual effects** (monotone and/or concave, not linear). Weaker than LCE, so it sits between L0 and L1 in plausibility. It needs the sharp-set machinery (Theorem 1) with shape constraints. Harder; LP after discretisation.

### C. Cross-unit stability (moderately plausible; cross-unit)

9. **Temporal stability.** The contextual effects $K$ are the same across election pairs (e.g. 2017→2020 and 2020→2024). Each pair identifies its own $Q$, so a common $K$ gives an intersection, like ridings but across time. Watch for realignments.
10. **Spatial smoothness.** Neighbouring polls have similar rates, $\lVert\beta_i-\beta_j\rVert\le L\cdot\mathrm{dist}(i,j)$. A unit-level restriction, so it tightens poll-level bounds too. LP: yes (pairwise rows; sparse with a neighbour graph).
11. **Latent homogeneous groups.** Constancy within clusters of similar polls (Cho 2001), between L3 and L5. The clusters must be chosen without looking at outcomes.
12. **Riding effects linked to riding covariates (M1′).** Contextual effects vary with riding characteristics (urban/rural, incumbency) instead of being common.

### D. Dispersion (less plausible; tightens the sharp set)

13. **Bounded residual dispersion.** Polls deviate from the mean model by at most $\epsilon$: $\lVert\beta_i-M_H(X_i)\rVert_\infty\le\epsilon$. As $\epsilon\to0$ this approaches deterministic contextual effects; large $\epsilon$ gives back LCE. It also supplies poll-level bounds, which LCE alone cannot.

### E. Sensitivity and robustness tools (apply to any layer)

14. **Breakdown frontier.** For a substantive conclusion (e.g. "loyalty > 0.8"), report the smallest violation of each assumption ($\delta$, $\lambda$, $\epsilon$, $\sup\lVert\Delta\rVert$ from Prop. 9) that would overturn it (Masten & Poirier 2020).
15. **Emptiness tests at each rung,** with sampling-noise correction (Phase 3).
16. **Inference with many constraints.** Intersection-bounds methods (Chernozhukov, Lee & Rosen 2013) are needed whenever constraints come from many ridings, elections or neighbours.

### F. Repeated cross-sections and difference-in-differences

Setting: the same units observed in a pre- and a post-period, e.g. counties voting on Black suffrage before and after the Civil War, with rows veteran / non-veteran. Each period $t$ has its own identified $Q_t$ and unidentified $K_t$. A typical target is
$$\text{DiD}=(B_{v,\text{post}}-B_{v,\text{pre}})-(B_{n,\text{post}}-B_{n,\text{pre}}),$$
which is affine in $(K_\text{pre},K_\text{post})$ (Lemma 2), so every item below keeps the LP structure.

17. **No cross-period restriction.** The set for $(K_\text{pre},K_\text{post})$ is the product of the two periods' sets, and the DiD bound is a difference of intervals, so the widths add. This is the baseline. Within a period the accounting identity links the groups: with $\omega$ the group-$v$ weight, $B_n=(\bar\tau-\omega B_v)/(1-\omega)$. So each period contributes one unknown in 2×2, and the gap $B_v-B_n=(B_v-\bar\tau)/(1-\omega)$ amplifies uncertainty when one group is small.
18. **Differencing out the unidentified part.** Suppose the units, $X_i$ and $N_i$ are the same in both periods, and the unidentified part of the contextual effects does not change, in whichever parametrisation is defended. One natural choice in 2×2: *the non-veteran contextual slope is unchanged*, so the war shifted non-veterans' support by the same amount everywhere, regardless of local veteran share. Then the unidentified term cancels and **each group's change and the DiD are point-identified**, even though the levels are only bounded. This is the EI analogue of DiD removing time-invariant confounding (checked in `illustrate_did.py`: exact without noise; centred on the truth with noise, Monte Carlo SD 0.02).
    * Which object is held fixed matters. "Non-veteran slope fixed", "veteran slope fixed" and "$K$ fixed" are different assumptions with different estimates. Each must be justified substantively.
    * In the veteran case, a change in the non-veteran slope *is* a spillover (veterans persuading neighbours in proportion to local veteran share). So this assumption rules out differential spillovers.
    * **How many restrictions are needed (2×2).** With $B_v=\dots-r\,w_1$ and $B_n=w_0+m_n w_1$, where $r=E[NX(1-X)]/E[NX]$ and $m_n=E[NX(1-X)]/E[N(1-X)]$:
      $$\text{DiD}=\text{(identified)}-(r+m_n)\,\Delta w_1.$$
      So **one** restriction pinning $\Delta w_1$ suffices.
      * "Both groups' contextual slopes unchanged" ($\Delta w_1=0$ and $\Delta b_1=0$) is stronger than needed. Since $b_1=d_1+w_1$, its second half reduces to $\Delta d_1=0$, which is **testable** (equal quadratic coefficients across periods).
      * One-group version: $\Delta w_1=0$ alone. Non-veterans' trend is context-free; the context-heterogeneity of veterans' trend ($\Delta b_1=\Delta d_1$) is identified from the data.
      * In R×C, $(C-1)R(R-1)/2$ restrictions on $\Delta K$ are needed (2 in 2×3).
    * **Beyond item 17 without a link?** No. Absent a cross-period restriction, the product set is expected to be sharp (to prove). Weak links that do help:
      * a sign on $\Delta w_1$ (e.g. non-negative spillover), which gives a one-sided cut;
      * $\lvert\Delta w_1\rvert\le\delta$ (item 19), which tightens the DiD and the levels;
      * common $w_1$ (item 20), which point-identifies the DiD and tightens the levels.
19. **Bounded change (sensitivity).** $\lVert K_\text{post}-K_\text{pre}\rVert\le\delta$ (or $\lvert\Delta w_1\rvert\le\delta$). The DiD bound has width proportional to $\delta$, scaled by the sensitivity weights $\mu_{rs}/\mu_r$, and collapses to item 18 at $\delta=0$. This is the EI analogue of bounded-variation DiD (Manski & Pepper 2018) and of bounded departures from parallel trends (Rambachan & Roth 2023). With several pre-war elections, $\delta$ can be calibrated by how much the unidentified part could have moved between pre-periods ("relative magnitudes").
20. **Stable contextual structure gives intersection.** A common $K$ across periods must satisfy both periods' constraint sets. This tightens the *levels* in both periods, the same mechanism as riding dummies with time playing the role of riding, and an empty intersection refutes stability.
21. **First-difference contextual effects (robust to unit heterogeneity).** Let $\beta_{i,t}=\zeta_i+M_{H_t}(X_i)+e_{i,t}$ with an *arbitrary* time-invariant unit effect $\zeta_i$. Linearity in levels fails, but linearity in changes holds: $E[\Delta\beta\mid X]=M_{\Delta H}(X)$ and $E[\Delta T\mid X]=X^\top\Delta Q\,X$. The constraints come from both periods' unit DD polytopes, $\Delta\beta_i\in S(X_i,T_{i,\text{post}})-S(X_i,T_{i,\text{pre}})$ (a Minkowski difference). This is weaker than LCE in levels and fits "unit fixed effects" logic. The bounds are wider, but valid under persistent local political cultures.
22. **Parallel trends for the causal claim (separate from EI).** Interpreting the DiD of group rates as the effect of service needs group-level parallel trends: absent the war, future veterans' and others' support would have moved together. Partial-identification version: the causal effect lies in DiD $\oplus[-M,M]$ (Rambachan & Roth 2023). Multiple pre-periods allow a placebo: if the bound on a pre-period DiD excludes 0, parallel trends are refuted.
23. **Monotone response (Manski MTR).** E.g. service did not lower support, $\Delta B_v\ge\Delta B_n$. Sign restrictions as LP rows; plausibility is case-specific.
24. **Composition change between periods (an A1 issue that is central here).**
    * Pre-war "veterans" are *future* veterans, so pre-war $X_i$ is unobserved. Many were under voting age pre-war, and war deaths and migration change the electorate.
    * The two cross-sections then cover different people, and the DiD of group rates differs from the average within-person change.
    * Options: (a) interval-valued $X_\text{pre}$ from age and mortality data, making EI with interval regressors (more partial identification); (b) explicit entry/exit rows; (c) a within-person target, which is a 3-way (veteran × pre-vote × post-vote) table with one-way margins, a multi-way Fréchet / transition problem.
    * Also: $T$ must use the same base as $X$ (yes votes / eligible men, with an abstain column), so the 2×2 becomes 2×3. R×C is needed here, not optional.

## Notes for implementation (Phase 4)

* Items 1, 4–7, 9, 10, 12 and 17–23 keep everything as one LP: they add rows or variables to $Gk\le g$ (Prop. 7). That makes a single `constraints=` interface practical.
* Unit-level items (1, 3, 5, 10, 13) matter most for the sharp set and for poll-level bounds. Mean-level items matter for aggregate bounds.
* Suggested reporting format: one table with the rungs as rows and the targets (loyalty, defection, demobilisation, mobilisation) as columns. Each cell gives the bound, its CI, whether the set was empty, and the breakdown value of the assumption that rung adds.
