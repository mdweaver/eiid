# Referee report on M-SIGN for ecological inference

**Recommendation — major revision.** Keep the option to impose a substantively motivated sign on a linear-contextual-effects model. The one-vs-rest identities explain exactly what such a belief contributes, and the package's distinction between removing a possible model and moving the target is valuable. Three issues change what an applied researcher could publish. The guide describes W2 as a direction for each group, although the software specification signs one sum. The “exact” N1 result uses divisible population mass, which can give a different verdict on reinforcement from the one possible with counted ballots. The simulation table overstates the length of a union, and its headline example can leave the impression that narrow intervals provide reliable information in settings where reinforcement is doubtful. In a stress test with correct linear means and 600 precincts, every restricted interval was nonempty and all 400 missed the true gap. I would revise the claims, examples and default reporting before recommending this module for election analysis.

**Scope.** Commit `8b98589653db71849fc7cdca500312754c26e95d`. Line numbers below refer to that commit. **2c** is `PHASE2c_monotone.md`, **5a** is `PHASE5a_sign_module.md`, **Guide** is `PLAIN_LANGUAGE_GUIDE.md`, and **Vig** is `eiPartialID_rxc.Rmd`. Implementation statements describe 5a §5. I reviewed the prescribed documentation and Python checks. A package-code audit was outside the brief.

## Assessment by belief and feature

| Belief or feature | Meaning correct? | Credible or useful where | Main risk of misuse | Location, section and lines |
|---|---|---|---|---|
| Between sign, δB | Qualified | Loyalty to a previous party, or racial bloc voting with independent evidence on mixed precincts | Treating a weighted average as an ordering in every precinct or against each rival | 2c §1, 42–47, 71–82. Guide §9, 211. Vig 153–154 |
| Within sign, δW | Qualified | Independent evidence that within-group support tends to track the focal group's prevalence | Inferring an individual response to neighbours, a common slope, or a sign from the aggregate scatterplot | 2c §§1, 3, 45–47, 144–153. Guide §9, 212. Vig 155–156 |
| W1 | Yes | The substantive comparison really is A versus a changing mix of everyone else | A change in the mix outside A can reverse the contextual sign | 2c §1, 44–47, 76–83. 5a §1, 52–59 |
| W2 | No in the guide | Analyst can defend an average after conditioning on each original row | Reading one pooled covariance as a restriction on each group's slope or on a Goodman coefficient | 2c §1, 45, 55–68. Guide §9, 214. 5a §5, 353 |
| Pairwise | Qualified | A versus a specified B in precincts containing both | In pairwise CR, interpreting the contextual term as specific to just A and B | 2c §1, 42–45. 5a §1, 92–94. Guide §9, 214 |
| Contextual reinforcement (CR) | Qualified | Independent reasons to expect contextual and between-group directions to agree under the chosen weighting | Treating an unknown common direction as a weak assumption, or feasibility as validation | 2c §2, 93–97, 113–133. 5a §1, 52–61 |
| Exact bounds under realized signs (N1) | Qualified | Known margins, an actual-rate estimand, and a continuous-rate interpretation | Claiming all values or endpoints are attainable with whole ballots, or requesting this from the shipped linear module | 2c §2, 90–111. 5a §5, 381–385 |
| Monotone conditional means (N2) | Qualified | A direction can be defended for each designated row over the chosen two-row composition range | Treating a global sign as monotonicity, treating a smooth derivative as observed at isolated support points, or reporting an outer bound as sharp | 2c §4, 165–186, 201–212, §6, 229–236 |
| Redundancy diagnostic | Qualified | Comparing assumptions on the same fitted base and target | Treating estimated implication as empirical confirmation, or an unchanged hull as an unchanged set | 5a §2, 125–145, 157–161, 200–210, §5, 357, 384 |
| Disjoint ranges | Yes for one CR | W2 or pairwise reinforcement with both sign branches feasible | Filling the gap, averaging the branches, or counting overlap twice when computing length | 5a §1, 47–50, 92–117, §5, 355–359. Guide §9, 220 |
| Descriptive Goodman and neighbourhood values | Qualified | Summarizing the aggregate pattern and showing the endpoints used in the W1 argument | Reporting either as the true gap, reversing endpoints for negative slopes, or clipping a confidence set | 2c §1, 49–60. 5a §3, 258, 309–311, §5, 362. Vig 160–165, 211–218 |

## Major comments

### 1. State whose behaviour is averaged, and with which weights

Let the focal row be A, the outcome be a vote for P, and the precinct share of A be $x_i$. Write $p_i=N_i/\sum_j N_j$, where $N_i$ counts the population represented in **both** margins. In the implemented module the row rates below are fitted conditional means $m_{is}$. For N1 they are realized rates.

* **Between.** The one-vs-rest belief says that the A-minus-rest gap is nonnegative on average, weighting precincts by $p_i x_i(1-x_i)$. Pure precincts get zero weight in this comparison. Pairwise A versus B uses $p_i x_{iA}x_{iB}$ and compares only their two rates. The sign of the sum over rivals need not hold against any named rival in particular. A large precinct with both groups can matter more than many small precincts where the ordering reverses.
* **Within, W1.** Compare precincts with different A shares, first among A voters and then among everyone else. Combine the two covariances using the population shares of A and everyone else. Either covariance can be negative when their sum is positive. The predictor in **both** covariances is A's share. “Own share” in the vignette gives the wrong predictor for the complement. Each covariance also reflects sorting, campaign targeting and other differences between precincts. A causal interpretation needs a separate design.
* **Within, W2.** For focal share $X_A$, the signed scalar is $\Delta_A=\sum_s\mu_s\operatorname{Cov}_{w^s}(X_A,m_s)$ with $w_i^s=p_iX_{is}/\mu_s$. Conditioning on the original rows removes the changing-composition effect inside the complement. The restriction still allows individual summands to have either sign. Also, $g-B=M^{-1}\Delta$ for multiple Goodman regression. Signing $\Delta_A$ generally does not sign the bias of A's fitted coefficient.
* **Pairwise CR.** In the pairing described in 5a line 94, the between term compares A with B and the within term is $\Delta_A$, which includes all rows. An argument based entirely on A and B can justify the former while leaving the latter open. Print the two terms separately alongside every CR result. At zero the weak product condition admits either direction for the other term.

A small example shows why the choice matters. The rows are A, B and C in that order. The units have equal populations. From unit 1 to unit 2, **every original group's P rate rises by two percentage points**.

| Unit | Shares of A, B, C | P rates for A, B, C | P rate outside A |
|---|---|---|---:|
| 1 | .20, .70, .10 | .40, .80, .10 | .7125 |
| 2 | .60, .05, .35 | .42, .82, .12 | .2075 |

W1 is **−.02633** and W2 is **+.00110** for A's share. As A becomes more prevalent, low-support C voters replace high-support B voters within the complement. These are association calculations and require no regression or rank condition. I would put this example next to the W1/W2 selector and describe W2 as “average within the original groups.” If a user wants each named group to satisfy the direction, that takes separate restrictions.

### 2. An argument for a between-group sign rarely supplies an argument for reinforcement

| Application and denominator | Plausible argument for a between sign | What must also be credible for a within sign or CR |
|---|---|---|
| Racial bloc voting, all eligible people | A representative individual survey or validated record linkage can support higher A support for P within mixed precincts. A statewide racial gap alone does not answer this question. | Mobilization in high-A precincts could reinforce the gap. Countermobilization, different candidate slates or turnout barriers could reverse the contextual association. In a multiracial jurisdiction, name the comparator and show which mixed precincts carry weight. |
| Partisan transitions, previous and current choices including abstention | Prior P voters may be likelier to vote P again than prior rival voters. The case against all other rows also depends on prior abstainers, who may enter in large numbers. | Safe-seat demobilization could lower P turnout precisely where prior P voters are concentrated. A sign for casting a P ballot among eligible people can differ from a sign for supporting P among votes cast. New entrants, movers and deaths need reconciled margins before either comparison. |
| Veteran status and a suffrage referendum, yes votes per eligible man | A service-linked survey, poll book with individual choice where legitimately available, or comparable historical evidence might support a veteran/non-veteran difference among otherwise comparable towns. Enlistment alone supplies no sign. | Veterans' organizations could mobilize neighbours. They might instead organize most intensely where veterans are scarce, or high-veteran towns might differ in party composition. A veteran loyalty argument alone says little about non-veterans' covariance with the veteran share. Distinguish eligibility, turnout and choice. |

**Suggested reporting practice.** Before examining the tightened result, specify the row, column, population, geography, sign or pairing, and independent reason for each restriction. For auxiliary records, estimate the same weighted associations and describe geographic coverage, linkage and transport to the target electorate. Show results without the belief and, when its direction is contested, under the alternative. With district-specific regressions, the implementation describes **one pooled restriction** (5a §5, line 353). It neither enforces the sign in every district nor removes cross-district composition differences.

### 3. Distinguish three kinds of sharpness and two targets

**EGGHM.** The v1 definitions are the expected conditional covariances. Lemma 1 gives the variance-ratio link between the neighbourhood and regression gaps. Propositions 2–3 supply the accounting bounds, 4 and 6 the bias identities, 5 and 7 the sign bounds, and Definition 1, Proposition 8 and Corollary 3 the reinforcement bracket and direction. Propositions 9–10 concern group means, and Theorems 1–3 distinguish the gap and the two means. Section 4 explicitly suggests the W1 collapse. These are faithful mappings when applied to a bounded column indicator and the population masses used in the paper. W2, the proposed pair weighting and reinforcement with those pairings are extensions. EGGHM's mechanism is a reason one might believe an association sign. Its identification result makes no causal claim. For local Propositions 12–13, distinguish a neighbourhood-specific gap under their Assumption 1 from a gap averaged over neighbourhoods at the same prevalence (appendix Propositions 14–16). [Elzayn et al., arXiv v1](https://arxiv.org/pdf/2504.14752v1)

**JKST plus a sign.** The quadratic regression identifies the symmetric part $Q$. Each admissible skew part $K$ gives a different array of linear group means with the same aggregate regression. M5 substitutes those means in each covariance, so a signed moment adds a linear inequality in the free coordinates of $K$. This is a coherent intersection under the common design, weighting, support and conditional-mean assumptions. No independent gain follows from feeding a range for a sign statistic implied by the base model back into the same fit. Tightening requires an additional belief that excludes previously feasible values. Inferring the realized-population sign from the model-implied sign, demanding it separately within districts, or treating the sign as an order at every composition would each add a further assumption. The base model also assumes the conditional means do not vary with precinct size given composition (2a §0.2, lines 21–25). Population weighting by itself cannot repair a violation. [Jiang, King, Schmaltz and Tanner, *Ecological Regression with Partial Identification*](https://gking.harvard.edu/publication/ecological-regression-with-partial-identification/)

**What R3 delivers.** Call the point output an *estimated bound for the model-implied gap, using the aggregate mean regression and the design*. The population bound can be wider than the sharp set using the full distribution of $T$ given composition and size. For example, suppose all precincts at each of several mixed compositions have $T$ equal to 0 or 1 with equal probabilities. Both groups then have rate $T$ in every precinct and mean ½ at every composition. The mean regression $\tau(x)=½$ alone also fits unequal, linear group means $m_1=½+a(1-x)$ and $m_0=½-ax$ for small $a$. Those means could arise with different returns having the same regression. They cannot reproduce these extreme returns. This illustrates exactly what the word “regression” costs in Guide line 215. A realized election gap additionally includes unit residuals. Report whether the CI targets $\tilde D$ or the actual $D$ with its separate residual allowance. R10's conditions and coverage attach to the specified scalar target.

**What M2 delivers.** For continuous cell masses with known margins, linear restrictions and a scalar target, the LP is sharp. With counted ballots the assertion at 2c lines 99–109 that every value is attained by a unit table fails. Consider two four-voter precincts with respectively one and two A voters, and **one P vote in each**. If $k$ of the two P votes came from A, the actual gap is $k/3-(2-k)/5$. Its possibilities are $-2/5$, $2/15$, and $2/3$. Both precincts have $T=1/4$, so $ER=NM=0$ and CR requires $D=0$. The continuous formulation is feasible at zero. **No integer table satisfies CR.** Even endpoint sharpness fails under the positive between sign: its continuous lower endpoint is zero, while its smallest attainable integer gap is $2/15$. EGGHM bound probability distributions with fixed neighbourhood masses. Their setup does not establish whole-ballot attainability. Label N1 “sharp for continuous rates given the margins.” If ballot-exact inference is intended, use integer count variables and report discrete attainable values or explicitly label an endpoint envelope. This affects the promised N1 interpretation, not the implemented mean-level LP.

### 4. Say exactly what a refutation establishes

For a single column of a two-row continuous table with nonconstant composition, assigning both rows the observed $T_i$ in every precinct gives $\delta_B=0$. It satisfies either weak between sign and CR, whatever $\delta_W$ is. In fact the equal-rates witness works for all rows and columns without extra cell restrictions. For W1 it puts the gap at $NM$. In contrast, $\delta_W\ge0$ requires $D\le ER$, so even the model-free positive within sign fails when $ER$ is below the Duncan–Davis lower bound for $D$. Distinguish the existence of this model-free check from a rejection returned by the shipped LCE fit.

With just the margins and continuous rates the claimed non-refutability is sound. The phrase “only with straight lines” is too narrow outside that comparison: known cells, structural zeros, integer counts or other credible restrictions can also rule out the witness. An empty **population** outer set rejects the stated conjunction. An empty plug-in estimate and rejection by an uncertainty-adjusted test are different events (5a §3, lines 218–256, §5, line 360). Show the base-only result alongside the joint result. If both accept, the sign still needs its original substantive justification. If just the joint test rejects, at least one element of the stated conjunction fails, subject to sampling error.

Correct M3(d)'s branch count. When all units are pure, $\gamma=1$ and a nonzero $ER$ is possible. Both latent associations vanish, so both weak-sign branches contain the same known gap. The proof at 2c line 127 recognizes this case, while the theorem at line 117 says exactly one branch survives. Require a mixed positive-weight unit as in EGGHM v1 for that assertion. With $ER=0$ the singleton describes the **gap's image**, and can correspond to many latent tables.

### 5. Report disconnected sets and uncertainty in a form a reader can use

* For a single W1 restriction the **population identification** region adds the interval between $ER^*$ and $\gamma ER^*$ in either order. Fix Vig line 163: the displayed ordering is reversed for the negative gap in 5a line 373. W1's single-LP result does not license selecting an inferential branch from the sign of the observed regression. For a single W2 or pairwise CR, a scalar projection may have two pieces, one piece, or none. A failed common-point certificate leaves overlap undecided. Several CR restrictions can generate more than two branches and potentially more than two projected pieces. The two-piece theorem concerns one CR. Clarify what the summary columns preserve when the full piece list is longer.
* Print the estimated identification set and the confidence set on separate lines. Preserve every connected component of each, allow point components, and label the population and units (percentage points for a gap). State that values inside an estimated gap are excluded by the fitted assumptions. The confidence set can still cover those values. Branches have no probability weights. A hull includes every internal gap, and marginally attainable endpoints from different cells need not form a jointly attainable transition table.
* R10's full-level branch union is sensible because one true branch suffices for coverage. No branchwise Bonferroni split is required for that scalar target. It gives no automatic simultaneous claim for all cells. The conservative W1 method inflates the estimated pieces before taking extrema. The default simulated-candidate CLR method has conjectural coverage even under a correct model. State this on the result line and show the conservative result by default or alongside CLR. The 0.359 width and 0.992 coverage in 5a line 300 refer to the **conservative** method. The separate default-CLR experiment at line 377 has different results. Uniformity near point identification remains open.
* Mark the descriptive values “weighted aggregate regression slope” and “neighbourhood gap,” name the focal group, give $\gamma$, and state whether the calculation uses fitted $\tau$ or observed $T$. Agreement under size-weighted least squares follows from its normal equations. Under unit-weighted fitting, equality need not hold and the size of the discrepancy is empirical. The voter-vaccination replication is a check of arithmetic. It cannot establish voter-level validity in another election. Intersecting an estimated bracket with a confidence interval treats a noisy endpoint as known: at a binding edge, downward estimation error can remove the truth even as both errors shrink. The observed 8–10% miss rates concern particular designs and procedures. They are not a universal error rate.
* R2's three verdicts work on a nonempty convex base with a feasible added convex restriction. “Implied” concerns all feasible parameters. “Cuts the set only” is relative to one reported target. Label each verdict with its base, target, and first/last ordering. A zero last-added width gain does not by itself establish logical redundancy. With unions, compare membership or merged components: removing the middle of $[0,1]$ leaves its endpoints unchanged. Report both total occupied length and hull span, and label the width-only fallback at 5a line 357 as inconclusive about logical implication. Near ties the implementation lacks a CI for the width difference (line 384). Show the fitted gain and base association interval, and mark the verdict as estimated.
* Include DD for the realized target, the unrestricted LCE result, and LCE plus the stated belief. Identify the target of each column. Show sign direction, comparator, weighting, baseline and joint rejection flags, number of units, composition range and size concentration. If a point-identified rung ignores an entered sign (5a line 363), check whether the fitted point violates it or mark the combination unevaluated. A precisely estimated point can still contradict a belief.

### 6. Correct the simulation table and show the price of the narrower result

The reference script sums the lengths of the two branch intervals before merging overlaps (`verify_R10.py`, lines 181–195). This preserves its membership calculation but double-counts length. In the near-zero design the two intervals overlap in **591 of 600** replications. The median length of their union is **0.242**, compared with **0.430** in 5a §3.4. The slab is also 0.242 to three decimals. The large apparent length penalty for the branch union in that row disappears. The headline comparison of 0.359 with 1.106 uses the slab and survives this correction. Replace the union column as detailed below. Distinguish total component length, hull span and sum of unmerged branch lengths throughout.

The source design holds precinct composition and size fixed within a scenario. It draws interior shares from a rescaled Beta(2,2), sizes nearly uniformly from 300–3,000, and two independent beta-distributed rates at each precinct around correctly specified linear means. The boundary design deliberately selects cancelling slopes. This isolates a binding constraint and sampling uncertainty. It gives little basis for predicting gains in segregated precincts, elections with uneven turnout or spatially dependent campaigns, historical towns dominated by a few large units, or three-party transitions. Independent beta draws also produce continuous rates rather than ballot counts.

In my matched-regression experiment, the observable mean of P is the same linear function in every scenario, the between sign is positive throughout, and only the truth of CR changes as the slopes move below zero. Baseline coverage stays high. With 600 units, the conservative CR procedure covers **398/400** times on the edge, **110/400** under a moderate departure, and **0/400** under a larger departure. It returns a nonempty interval every time. Precision can make the wrong substantive answer more reproducible. Put a misspecification comparison next to the existing successful example. Expose one-belief-at-a-time removal now and, if a slack parameter is later implemented, calibrate it in interpretable units. A sensitivity parameter for departures from linearity answers a different question.

### 7. N2 requires an explicit support and conditioning argument

The two-row shape restriction is stronger than a pooled covariance sign. Both row means increasing implies the global mean-level within sign; an average sign permits reversals. On a finite enumerated support, adjacent inequalities represent monotonicity there. On a continuous population support the distinct sampled shares are only a grid. Multiple extensions agreeing on sampled points may violate the accounting restriction elsewhere. On genuinely isolated support, differentiability of an extension does not identify $\tau'(x_j)$: adding $a\prod_j(x-x_j)$ to a smooth fit preserves every fitted support value while changing the derivative at each node. Separate the finite-support value LP, a discretized outer approximation and a local derivative model with support and smoothing conditions. Label derivative-based results as conditional on those conditions (2c §§4, 6, lines 177–190, 201–212, 233–236).

Two further details matter if N2 is developed. First, its defining model conditions on $V=(X,N)$. A sharp conditional random-set calculation would generally need $m(X)\in E[S(X,T)\mid X,N]$ almost surely, including compatibility across sizes at a given composition. Conditioning on $X$ alone as in 2c line 182 can discard restrictions. Second, for binned observations the identity involves $E[X m_1(X)\mid\mathrm{bin}]$. It generally differs from $E[X\mid\mathrm{bin}]E[m_1(X)\mid\mathrm{bin}]$. Define the binned row mean using the row's own population weights, prove the resulting nesting for the target, and give conditions for fitting $\tau$. A common averaging operator and convexity alone do not prove the proposed binning result. Cross-validation of polynomial degree likewise needs post-selection justification for inference. These concerns apply to a proposed rung, which §5 explicitly excludes from the current build.

## Minor comments and missing applied checks

1. Add a hand-worked turnout example. A prior P voter's probability of casting a P ballot is the product of their turnout probability and their conditional P support among participants. Positive loyalty among participants can coexist with a negative eligible-population gap. Use an exhaustive abstention column. Where turnout is measured only among registered voters, say who enters that denominator. A prior nonvoter and a newly eligible person need different treatment when their populations can be distinguished.
2. Report the largest precinct weight and the effective number $1/\sum_i(w_i^r)^2$ for each important row. For between comparisons also report total overlap and the largest normalized overlap weight. Nearly pure precincts can identify much of a district rate while contributing little evidence about a between sign in mixed precincts. A leave-largest-unit-out display would reveal dependence on one city. The root-number-of-units residual claim at 2c line 160 needs conditions on weights and interunit dependence.
3. Position this module beside King's EI and hierarchical R×C estimators. Those approaches can produce distributions over latent tables through additional distributional and pooling assumptions. M-SIGN reports feasible scalar ranges and sampling uncertainty under specified moments. A posterior mean inside a sign bound gives no independent support for the sign. A particularly informative benchmark would aggregate a held-out individual panel with both election choices and abstention observed, lock sign choices using other units, then score actual-gap inclusion, width, overlap and rejection by precinct size. Where linked ballots are unavailable, use a documented mock election or linked survey benchmark and report selection and measurement error. [Rosen, Jiang, King and Tanner, *Bayesian and Frequentist Inference for Ecological Inference: The R×C Case*](https://gking.harvard.edu/publication/bayesian-and-frequentist-inference-for-ecological-inference-the-rxc-case/)
4. In a 2×3 transition, signs on multiple columns can sharpen another column via the row totals. In a 3×2 transition there is also a substantive choice of complement, comparator and contextual statistic. Neither dimension alone predicts the size of the gain. State that restrictions on all outcome columns must respect adding-up, and that signs for complementary binary outcomes reverse. Mark the ladder's “both signs” under CR at 2c line 218 as *both weakly positive or both weakly negative*. Opposite signed restrictions generally do not imply reinforcement.
5. Use one status key. The current build documents mean-level association rows and CR. N1, N2, local monotonicity and sign slack are outside the implementation (5a line 385). The planned API at 2c line 7 is different from the shipped `restrictions` syntax. Make sure readers can tell which worked example they can run today.

## Numerical explorations

The two supplied identity-check scripts reproduced their reported totals (2c, 12/12; Phase 1, 10/10). I used the pinned R10 simulator as a Python reference for conservative intervals. These experiments do not exercise the package's CLR implementation. The independent checks above on W1/W2 and on integer attainability use the displayed precinct tables.

**A. Correcting the union lengths.** Same seed, design, 600 repetitions per scenario, 4,000 Gaussian draws per critical value, and the original call order. Widths are in probability units. Nonempty branch intervals are sorted and their intersection is counted once.

| Reference design | Published median, sum of branch lengths | Median length of the union | Runs with positive overlap | Published median, slab |
|---|---:|---:|---:|---:|
| Edge, $ER^*=.15$ | .386 | .357 | 134/600 | .359 |
| Edge, $ER^*=.005$ | .430 | .242 | 591/600 | .242 |
| Interior | .480 | .480 | 0/600 | .480 |

**B. Deliberately wrong CR.** Using the reference design generator and sampling code, set $m_0(x)=.40-sx$ and $m_1(x)=.50+s-sx$. All rates are valid throughout $[0,1]$. For every $s$, the aggregate regression is exactly $\tau(x)=.40+.10x$ and the within-precinct difference is $.10+s>0$. For $s>0$, both contextual slopes are negative and CR is false. In fact the true model-implied gap is $.10+(1-\gamma)s$, strictly above the CR upper endpoint of $.10$ for mixed precincts. There is always an alternative CR-compatible mean model with the identical regression, namely $s=0$. The simulation fixes the regression across scenarios. It need not fix the entire observable distribution, since beta variances can change. I held the design fixed across the three scenarios at each sample size.

| Precincts | $s$ | True gap | LCE covered / 400 | CR slab covered / 400 | Median LCE width | Median CR width |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 0 | .100 | 400 | 399 | 1.262 | .306 |
| 150 | .08 | .170 | 400 | 336 | 1.256 | .304 |
| 150 | .25 | .318 | 400 | 2 | 1.254 | .303 |
| 600 | 0 | .100 | 400 | 398 | 1.074 | .194 |
| 600 | .08 | .169 | 400 | 110 | 1.069 | .196 |
| 600 | .25 | .317 | 400 | 0 | 1.066 | .195 |

No baseline or conservative slab interval was empty in these 2,400 runs. The binomial Monte Carlo standard error for each coverage fraction is at most .025. Zero of 400 is a finite simulation finding. The analytic exclusion of the true gap from the limiting bracket explains the deterioration as the number of precincts grows. A false belief that happens to leave a different target's true value inside the restricted projection need not produce this pattern for that target.

## Outside scope — background documents

* **1b §A2, lines 91–98.** The preliminary 2×2 discussion says CR is nonconvex and requires two LPs. The later R3(d) shows that the matched W1 projection reduces to a slab, including in two rows. Mark the earlier statement superseded. The same paragraph's “different moments” does not itself prove non-nesting. What matters is which models and which targets are restricted.
* **1b §D.1, lines 305–311.** The account gives each W2 sign to one identifiable group. The formula in its own §A3 (line 138) averages over all original groups. Carry the correction through the owner-facing description. In §A1, lines 66–69, attribute Table 1 to Theorem 1 and Tables 2–3 to Theorems 2–3. In lines 84–87, describe the linked vaccination evidence as support for sign transport, with selection and party-measure caveats acknowledged by EGGHM §6.
* **Guide §2, lines 23–29.** A3 (no size dependence in mean behaviour) and A5 (composition rank) are labelled shared requirements. The realized, continuous DD accounting bounds need neither, and a finite N1 LP does not need a quadratic regression. State the respective denominators and rank requirements at each rung.
* **1b §D.2, lines 313–318.** A reversal in aggregate yes-share is not by itself evidence that either group's conditional mean reverses. Changing mixture shares can create it. Plotting returns provides context, while feasibility under a specified composition curve and monotone cell means addresses the actual restriction.

## Comparison with the earlier internal review

I wrote the recommendation, assessment, major comments and minor comments before reading 5a §4. I agree with its separation of conservative slab inference from CLR branch inference, its distinction between one empty branch and two disjoint pieces, and its correction to conditional coverage. For a fixed procedure, coverage implies a nonempty confidence set. Conditioning on the event that this particular set is nonempty cannot reduce its coverage probability when that event has positive probability. This argument supplies no protection when the researcher selects a belief because it gives an attractive result.

My recommendation differs from the internal “accept with fixes” for the narrower R3, R2 and R10 review. Reproducing the printed R10 table, as §4 line 342 reports, does not check whether its length statistic measures the union. The overlap calculation, the counted-ballot counterexample, the W1/W2 sign reversal and the false-CR experiment all change how I would communicate the results. None requires discarding the central mean-level LP construction.

## Reproducibility appendix

Run the following as one Python script from the pinned commit's `docs/extension/` directory with NumPy installed. The reference script is loaded fresh for each experiment so the first simulation's random stream does not carry into the second. The integer and W1/W2 checks at the end are independent of the reference simulator.

```python
from itertools import product
from pathlib import Path
import runpy
import numpy as np

reference = Path("verify_R10.py")

def union_experiment():
    r = runpy.run_path(str(reference))
    for kind in ("boundary, ER far from 0", "boundary, ER near 0", "interior"):
        x, p = r["design"](150)
        theta, ab = r["make_truth"](x, p, kind)
        parts = r["linear_parts"](x, p)
        rows = []
        for _ in range(600):
            intervals = r["one_rep"](x, p, theta, ab, 15, parts)["CR-union"]
            if intervals is None:
                continue
            intervals = sorted(intervals)
            raw = sum(hi - lo for lo, hi in intervals)
            overlap = 0
            if len(intervals) == 2:
                (lo1, hi1), (lo2, hi2) = intervals
                overlap = max(0, min(hi1, hi2) - lo2)
            rows.append((raw, raw - overlap, overlap))
        a = np.asarray(rows)
        print(kind, np.median(a[:, :2], axis=0), sum(a[:, 2] > 0))

def misspecification_experiment():
    r = runpy.run_path(str(reference))
    # All reference functions share this global generator.
    r["one_rep"].__globals__["RNG"] = np.random.default_rng(1739)
    for n in (150, 600):
        x, p = r["design"](n)
        parts = r["linear_parts"](x, p)
        for s in (0, .08, .25):
            ab = (.40, -s, .50+s, -s)
            theta = np.array([.40, .45, .50, -.05-s])
            truth = parts[3](theta)["D"]
            results = {name: [] for name in ("JKST", "CR-slab")}
            for _ in range(400):
                out = r["one_rep"](x, p, theta, ab, 15, parts)
                for name, records in results.items():
                    ci = out[name]
                    if ci is None:
                        records.append((1, 0, np.nan))
                    else:
                        lo, hi = ci
                        records.append((0, lo <= truth <= hi, hi-lo))
            for name, records in results.items():
                a = np.asarray(records)
                print(n, s, name, truth, a[:, :2].sum(axis=0),
                      np.nanmedian(a[:, 2]))

def small_tables():
    x = np.array([[.2, .7, .1], [.6, .05, .35]])
    m = np.array([[.4, .8, .1], [.42, .82, .12]])
    p = np.array([.5, .5])
    focal = x[:, 0]
    def contribution(share, rate):
        weights = p * share
        center = np.average(focal, weights=weights)
        return weights @ ((focal - center) * rate)
    rest = np.sum(x[:, 1:] * m[:, 1:], axis=1) / (1-focal)
    w1 = contribution(focal, m[:, 0]) + contribution(1-focal, rest)
    w2 = sum(contribution(x[:, j], m[:, j]) for j in range(3))
    print("Rest rates, W1, W2", rest, w1, w2)
    # Each four-voter precinct has one P vote, which can go to either row.
    gaps = {k/3 - (2-k)/5 for votes in product((0, 1), repeat=2)
            for k in [sum(votes)]}
    print("Integer-feasible gaps", sorted(gaps))

union_experiment()
misspecification_experiment()
small_tables()
```
