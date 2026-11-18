# Referee report B · Partial identification and inference

**Recommendation — major revision.** The main geometric results survive scrutiny. In the stated pooled linear contextual effects (LCE) model, R3 distinguishes sharpness for the full observable law from sharpness for the regression and design. R2 gives a necessary and sufficient test for a restriction to move an endpoint, including when the optimum is an entire facet. R10 gives a valid way to handle a union of restrictions, conditional on a valid interval for the true branch. Four changes matter for users of the bounds. Integer electorates can reject contextual reinforcement (CR) when the “exact” LP accepts it. District intercepts can reverse the pooled within-association even when all contextual slopes are positive. The observed compositions do not generally identify the derivative used by the local monotonicity LP. Finally, the inference guarantees depend on the candidate construction and asymptotic regime. The report gives counterexamples, replacement statements and a fixed-design proof of persistent undercoverage for the naive CR interval.

**Reviewed version.** All line numbers refer to commit `8b98589653db71849fc7cdca500312754c26e95d`. In the tables, **2c** denotes `PHASE2c_monotone.md`, **5a** denotes `PHASE5a_sign_module.md`, and **P3** denotes `PHASE3_inference.md`, all under `docs/extension/`. I recorded the verdicts before opening 5a §4. I reviewed the listed theory and verification scripts. The implementation was outside this review. [Pinned 2c](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE2c_monotone.md) · [Pinned 5a](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE5a_sign_module.md)

## 1. Verdicts

Severity refers to the consequence of the qualification or error for the result as advertised. CR means the union of the two weak same-sign branches, so a zero association matters.

| Result | Verdict | Severity | Location |
|---|---|---|---|
| M1 | **Established.** W1 is the collapsed within-group association. W2 is the vector entering the multivariate normal equations. Both follow from accounting. | — | 2c §1, 42–84 |
| M2 | **Established under the continuous-rate interpretation of the stated simplex.** Integer-count sharpness and integer refutation are incorrect. Slack nesting is componentwise. | Major for actual count tables | 2c §§0.1, 2, 15–19, 90–111 |
| M3 | **Incorrect** as a whole. Parts (a)–(c), the CR interval and its sign conclusion hold for the unrestricted continuous margin model. The branch-count assertion in (d) fails at an allowed boundary. | Moderate for branch count, major for integer refutation | 2c §2, 113–133 |
| M4 | **Established.** Joint feasibility implies each individual restriction. Strict improvement is possible and emptiness is possible. | — | 2c §2, 135–138 |
| M5 | **Incorrect** across all advertised rungs. Parts (a)–(b) hold using the full district parameterization and district-specific ordering. Part (c) holds for a common pair of lines or within a district. | Major | 2c §3, 144–160 |
| N2.1 | **Established.** Both conditional identities hold almost surely. | — | 2c §4, 167–175 |
| N2.2 | **Established under additional or stated conditions** for the finite-support outer LP, nesting and weighted covariance results. The proposed sharp replacement in (a) generally needs conditioning on all of $V$. Part (b) cannot fill missing support constraints. | Major for the proposed sharp extension and general-support reading | 2c §4, 177–198 |
| N2.3 | **Not established** as an identified local restriction under the stated support assumptions. The product rule and interval equivalence hold for compatible differentiable functions. The non-implication in (d) is true even for smooth functions. | Major | 2c §4, 201–212 |
| R3(a) | **Established** for the model-implied affine target, continuous latent rates, full observable law, and (C1)–(C2) in the stated pooled model. | — | 5a §1, 18–40, 65–70 |
| R3(b) | **Established under stated conditions** for the union over regression-equivalent laws. Calling its computation an LP requires a polyhedral restriction or an equivalent finite LP representation. | Minor terminology | 5a §1, 24–32, 42–45, 72–75 |
| R3(c) | **Established.** Compare the two scalar images, including empty and singleton images. A common parameter witness is sufficient. | — | 5a §1, 47–50, 77–79 |
| R3(d) | **Established** for one W1 pairing, positive composition variance, and a specified column. The slab can have empty intersection with the LCE set. | — | 5a §1, 52–59, 81–87 |
| R3(e) | **Established** as the joint refutation implication. The word “never” needs qualification when a belief also has a model-free observable implication. | Minor wording, consequential application at 5a line 309 | 5a §§1, 3, 61, 89–90, 118, 309 |
| R2(a) | **Established** given attainment on both nonempty sets. The displayed feasibility LP in proof step 4 omits the target intercept. | Moderate for a literal implementation | 5a §2, 125–138, 165–169 |
| R2(b) | **Established.** The statement concerns all optimal multipliers of the added row for the endpoint under examination. | — | 5a §2, 140–143, 171–176 |
| R2(c) | **Established** for the family of all linear targets on the free parameters and a closed convex restriction. | — | 5a §2, 145, 178–181 |
| R2(d) | **Established under stated conditions.** Positive variance, $\gamma<1$, a nonconstant gap target, and nonempty intersections for numerical gains. | — | 5a §2, 147–155, 183–186 |
| R2(e) | **Established under a target-specific meaning of redundancy.** The final sentence is incorrect if it means the row is implied by the others. | Moderate | 5a §2, 157–161, 188 |
| R10(a) | **Established under the augmented-system assumptions** for conservative inference, full-candidate CLR under its actual conditions, and simultaneous projection. Simulated-candidate CLR remains conjectural. Actual targets require a valid residual allowance. | Major if generalized across regimes | 5a §3, 218–235, 262–266 |
| R10(b) | **Established under stated conditions** for conservative inference. Simultaneous inflation must include all pieces and the target constant. | — | 5a §3, 237–243, 268–272 |
| R10(c) | **Established under stated conditions** for point coverage when at least one true branch has valid inference conditional on the design. No split across branches. | — | 5a §3, 245–249, 274–279 |
| R10(d) | **Established under stated conditions** pointwise, with the specified single emptiness rule and error allocation for the chosen version. No uniform near-boundary result follows. | — | 5a §3, 249–256, 281–289 |
| R10(e) | **Established.** The Monte Carlo demonstrates finite-sample failure. Its sketch gives a sufficient variance-ratio condition. Comment 7 supplies a design that verifies the condition analytically. | Proof completion | 5a §3, 258, 287–311 |

## 2. Major comments

### 1. A continuous LP and an integer electorate can reach opposite refutation decisions

In M2 the decisive step is 2c line 106. The equality $n_{irc}=N_iX_{ir}\beta_{irc}$ gives nonnegative **real** cell masses with the prescribed margins. It does not give integer counts. For real masses the proof is sound. Every feasible point is a table in that model. Convexity gives every intermediate target value, and compactness gives both endpoints. Closed polyhedral unit restrictions preserve this argument. Each CR branch has its own interval. Right-hand-side LP duality gives continuous polyhedral endpoint functions on the nonempty feasibility region.

Consider three units of size four. Their group-1 counts are $1,2,3$, and their numbers of positive outcomes are $0,1,0$. Each group contains six people. The single positive outcome belongs to one group or the other, so the possible integer gaps are exactly $\{-1/6,1/6\}$. Symmetry gives $D_{ER}=D_{NM}=0$ and $\gamma=1/6$. The continuous CR interval is $\{0\}$, attained by allocating half of the positive outcome to each group. **Both integer allocations violate CR.** This changes the endpoints and the emptiness decision, even with mixed units and without any auxiliary unit restrictions.

With integer margins alone, the DD extrema for one column are attainable by integer allocations. The sign identities remain exact, so each single weak between sign still has at least one integer witness on its side of NM, and the single within-sign emptiness test still works by comparing ER with the integer DD extremum. What can fail is attaining the cut point itself, filling the interval, or satisfying both weak signs together.

For continuous tables, M3(a)–(c) follow by intersecting the unrestricted DD gap interval with a half-line. Setting both row distributions equal to the observed column margins in each unit attains the NM gap jointly with the other column margins. Extra orderings, zeros or previously imposed signs can remove that witness. M3(d) has a separate error even for continuous rates. With two equally weighted pure units, $(X,T)=(0,0),(1,1)$, one has $D=ER=NM=1$, $\gamma=1$, and both associations zero. Both weak branches contain the solution. The proof at line 127 identifies the exception, then the proposition omits it.

**Correction.** Label N1 “sharp for continuous unit tables conditional on the observed margins.” If integer-table inference is intended, impose integrality on the cells and allow a discrete target set, with separate integer feasibility checks for the branches. In M3(d), exactly one branch is feasible when $ER\ne0$ and $\gamma<1$. At $\gamma=1$, both give the same identified gap. At $ER=0$, both give gap zero in the unrestricted continuous model. For M2(d), a scalar breakdown value additionally requires a specified conclusion, treatment of equality, and a scalar relaxation path. A vector of slacks generally gives a breakdown frontier.

### 2. The district rungs preserve affine rows but can reverse the sign predicted by M5(c)

For fixed design write $m_i=m_{\rm id}(V_i)+L_i k$. Under M2 the vector $k$ stacks district-specific free coordinates. Under M1 the intercept shifts enter the identified part. Under M1($\delta$) the common center may be an additional variable with zero coefficient in the association row. Every association in M1(d) is a design-weighted linear functional of the stacked means. Composition of these two maps proves M5(a) for each rung. District-specific ordering on the corresponding district support proves M5(b) by multiplying by $p_iX_{ir}X_{is}\ge0$ and summing. The single formula $M_{Q+K}(X_i)$ at 2c line 150 proves only the pooled specialization.

For M5(c), take six equally weighted units with $X=.1,.2,.3$ in district 1 and $X=.7,.8,.9$ in district 2. In both groups set the positive-outcome mean to $.8+.1X$ in district 1 and $.2+.1X$ in district 2. The other column takes the complement. This satisfies M1 with common strictly positive slopes, valid probabilities, and three distinct compositions per district. Direct computation gives pooled $\delta_W^*=-0.0492711$. In the within-group weighted covariance, the district intercept is correlated with $X$. The proof at line 153 treats that intercept as constant. For M1 the missing terms are $\sum_g \mu_g\operatorname{Cov}_{w^g}(X,\zeta_{g,D})$. Under M2 the slopes themselves can also vary by district.

**Correction.** Keep the affine-row conclusion for every rung. State the slope shortcut for M0, or compute it separately under weights conditional on each district. A pooled shortcut requires accounting for the between-district covariance. The same distinction applies if R3 or R10 is advertised for districts. Their extension must retain the district design, identified district coefficients, the appropriate $A_d$, and all free coordinates.

### 3. N2 needs a support convention, conditioning on size, and an identified derivative

**Finite support.** For N2.1, conditioning accounting on $X$ gives $m(X)^\top X=\tau(X)$ almost surely, and the row simplexes are closed and convex. On a complete finite support every atom has positive mass. Adjacent order constraints suffice, and a common piecewise-constant extension of the *whole matrix* preserves all row sums and designated directions. This proves N2.2(a)–(b) for the outer restrictions on that support. The signed-slope LCE functions satisfy those same constraints at each atom, proving (c). For (d), under each row's tilted distribution the covariance of $X$ and a nondecreasing function of $X$ is nonnegative, by the independent-copy covariance identity.

“Distinct observed compositions” and “support points” at line 177 can be identified if the design law is explicitly the finite empirical design. For a subset of a larger support, extendibility subject to the **other** accounting equations can fail. Take true support $\{1/4,1/2,3/4\}$, row-1 means $(.4,.4,.4)$, row-2 means $(8/15,.7,.8)$. Both rows are increasing weakly and $\tau=(.5,.55,.5)$. On a grid keeping only the extremes, assigning both rows the value .5 at both extremes passes the LP. Any monotone extension assigns .5 to both at the middle and cannot give .55. Dropping support constraints can yield a relaxation. It does not show that nothing is lost for the full model.

**Full-law sharpness.** Line 182 proposes $E[S(X,T)\mid X=x]$ for the sharp version. The model imposes $E[\beta\mid X,N]=m(X)$. Fix $X=1/2$ and give two sizes equal probability. For the first size let $T$ equal 0 or 1 equiprobably. For the second let $T=1/2$ surely. Both sizes have the same aggregate conditional mean. At the first size both row rates must be $T$, hence both row means must be .5. Conditioning only on $X$ allows the pair $(.75,.25)$ by selecting $(1,0)$ at the second size. The $S(x,\tau)$ LP also admits it. With $V=(X,N)$, the required sharp condition is
$$
m(X)\in E[S(X,T)\mid V]\quad P_V\text{-a.s.}
$$
together with the shape restrictions. Equivalently, enforce every conditional support-function inequality for the same $m(x)$ at almost every size given $x$. The selection construction of 2a Theorem 1 then applies to each admissible measurable $m$. Conditioning on $X$ alone suffices if size adds no information relevant to those selection sets, or if the model conditions only on $X$.

**Local slope.** N2.3 differentiates an almost-sure identity at a specified point. On an isolated support point even infinitely many repetitions give no information about a derivative there. Take $x\in\{1/4,1/2,3/4\}$, $m_1=.7,m_0=.3$. A second smooth extension adds $h(x)=(x-1/4)(x-1/2)(x-3/4)$ to both functions. Both are valid probabilities on $[0,1]$, agree on the complete support, and induce the same law there. Their structural aggregate derivatives at $1/2$ differ by $h'(1/2)=-1/16$. The true first model satisfies weak local CR with gap .4 and zero local within effect. Using the other compatible derivative gives an upper threshold below .4. An LP using that threshold need not cover the first model's gap.

**Correction.** State the calculus identity for the structural version $\tau_m(z)=z m_1(z)+(1-z)m_0(z)$. For an observable restriction, require the derivative to be identified, for instance at an interior point of interval support with compatible smooth versions, or under an identified, maintained regression family. Then $\delta_W(x)=\tau'(x)-D(x)$ gives exactly $D(\tau'-D)\ge0$, or $D$ between 0 and $\tau'$, including a zero derivative. An outer LP can omit smooth extendibility once the correct derivative is supplied. An estimated derivative adds sampling uncertainty. At a pure composition, EGGHM's local between association is $x(1-x)D(x)=0$; replacing its sign by the sign of $D(x)$ changes the belief there.

The last assertion of N2.3(d) is true. To avoid the kinks in the cited V6 example, take $X$ uniform on $[1/5,2/5]$, $m_1=x-x^2/2$ and $m_0=1-x^2/2$. Both are probabilities and $xm_1'+(1-x)m_0'=0$ everywhere. Integration under the two row weights gives $\delta_W^*=-1/47250$. A finite-design version with equal mass at the two endpoints gives $-1/5250$ for the same smooth structural functions. For the reverse non-implication, complement both functions in the continuous example, shrink them slightly into the interior, then subtract a sufficiently small positive multiple of $x$ from each. The aggregate within association stays strictly positive by continuity and the local within effect becomes strictly negative.

### 4. Binning and cross-validation in 2c §6 require different arguments

At line 236 the row simplex is fixed, but the accounting constraint changes with $x$. With equal bin weights at $x=.25,.75$, set $m_1=x,m_0=0$. Then $E[T\mid{\rm bin}]=E[X^2]=.3125$. Evaluating the vector of ordinary bin-average row rates at $E[X]=.5$ gives .25. Applying conditional convexity to the varying sets $S(x,\tau(x))$ does not prove the proposed step.

There is a repair for the **global two-row monotonicity** LP. Use person-weighted bin expectations and define
$$
\bar x_b=E_bX,\qquad \bar m_{1b}=E_b[Xm_1]/\bar x_b,\qquad
\bar m_{0b}=E_b[(1-X)m_0]/(1-\bar x_b).
$$
These rates reproduce $E_bT$, preserve full-design row targets with the correct bin row masses, and inherit monotone ordering across ordered, nonoverlapping bins. Treat absent groups by a compatible zero-weight extension. In the example the correct first-row bin mean is .625. This argument establishes an outer relaxation of the global mean constraints. Local derivative rows need a separate justification.

Line 235 also conditions on a polynomial degree chosen using the outcomes. Fixed-degree inference does not give the conditional law after such selection. For example, compare intercept-only and intercept-plus-slope fits by leave-one-out squared error in a balanced two-level design with independent equal-variance errors and true slope zero. Both candidate mean models contain the truth. Their leverages are respectively $1/n$ and $2/n$. Comparing ${\rm RSS}/(1-h)^2$ selects the extra slope asymptotically when its standardized estimate has $Z^2>2$. The nominal 95% slope interval then covers conditionally with probability $1-.05/P(|Z|>\sqrt 2)$, about 68%. **Correction.** Prespecify an identifiable degree or justify a selection-aware or independent-sample procedure and any approximation bias. If $J$ denotes the number of distinct points at line 177, degree $J$ at line 234 has $J+1$ coefficients and is not identified from those values. Use separate notation for support size and degree.

### 5. R3 is sharp for specified information. The refutability claim changes with that information

**Full law and regression.** In (a), a compatible latent law for $k$ exists by the definition of $\Theta_I(P;A)$. Imposing a restriction evaluated at the same $P$ simply discards incompatible $k$'s. Closedness of $C$ preserves compactness, convexity preserves intermediate affine target values. In (b) the information is exactly the design law $P_{N,X}$, the fixed $A$, and the quadratic regression $\tau(x)=x^\top Qx$ column by column. The comparison class has the same design and satisfies $E_{P'}[T'\mid X,N]=\tau(X)$. For every feasible $k$, the deterministic construction $\beta=M_{Q+K}(X),\ T'=\tau(X)$ witnesses attainability in this class. (C2) holds fixed the restriction when this law is substituted. Dependence on a conditional variance or the conditional distribution of $T$ would invalidate that step. (C2) is stronger than needed in (a). A ball is an admissible $C$ under (C1)–(C2), so sharpness alone does not supply a finite LP.

The dimension condition is sufficient. $R$ affinely independent points of $A$ are linearly independent in the unit-sum hyperplane. Bounded values of $K_c x$ there bound every matrix entry. A5 implies this dimension, since confinement to a proper affine hyperplane would give a nonzero quadratic vanishing on the support. $P(X\in A)=1$ ensures both the deterministic witness and the full target use the region on which the mean is specified. Trimming requires separate outside-$A$ terms. Fixing a design permits the auxiliary convention of assigning equal probability to its units to represent the finite-design target via a population law. Conditioning alone does not equate that target to an arbitrary superpopulation rate.

For a concrete distinction, let $X=1/4,1/2,3/4$ equiprobably, same size, and $T$ a fair binary total at every composition. The full law forces both row means to .5 everywhere on the support. The regression law with $T'=.5$ surely admits $m_1=.5+k(1-x),m_0=.5-kx$. For $C=[.2,.4]$, the full-law identified set is empty and the regression-only gap set is $[1/6,1/3]$. This is exactly the distinction (a)–(b) promise.

**Union and slab.** Images commute with unions. Two nonempty compact scalar intervals fill their hull exactly when they overlap or touch. One empty branch can be dropped, and two empty branches mean rejection. In the square $[-1,1]^2$, the affine pair $g_B=u+.5,g_W=u-.5$ gives two disjoint vertical strips. Target $u$ has a gap of length one. Target $v$ has identical branch images despite an infeasible common-witness LP. For W1 specifically, the accounting identities give positive multiples of $D-\gamma ER$ and $ER-D$. Their product gives the stated slab for every row. Ties can leave both branches feasible.

**Refutation.** A population-empty LP rejects the conjunction. Emptiness alone does not attribute the failure. A belief with a meaning outside LCE can sometimes also be tested separately. For example on $x=1/4,1/2,3/4$ let both row means be $1.5-2x$, and let the aggregate equal this function. Then LCE holds on $A=[1/4,3/4]$, while $ER=-2$ and any feasible gap is at least $-1$. Positive within association is impossible even in the model-free continuous table problem. For restrictions defined solely on a skew LCE parameter, interpretation of $C$ already presumes that model.

More importantly, the “unrefutable belief” language in 5a line 309 imports M3 into the wrong base model. Let $X=0,1/4,1$ equiprobably, $A=[0,1]$, $m_0=1-x,\ m_1=(49/55)x$, with deterministic $T$. The regression is $(1-x)^2+(49/55)x^2$. Its Bernstein cross coefficient is zero, so the two nonnegative off-diagonal entries of $H$ must both be zero. Thus LCE fixes both lines. Direct substitution gives $ER=0$ and $\tilde D=-261/1925$. CR requires $\tilde D=0$, so **both LCE branches are empty**, although LCE alone fits. The unrestricted continuous tables in M3 can still satisfy CR. Say explicitly which base system the non-refutation claim uses.

### 6. R2 certifies changes to a specified endpoint

Write $U=\max_F(\varphi_0+h^\top k)$. If any point of its optimal face survives, the augmented optimum is $U$. Otherwise the attained new maximizer has strictly smaller value. This proof works for every face dimension, even a constant target. For the feasibility LP, impose $h^\top k\ge U-\varphi_0$. The printed $h^\top k\ge U$ at line 169 can reject the entire old face. On $[0,1]^2$ with target $2+x$, the printed row is $x\ge3$.

For one added row let $b^*=\min_{F_U^*}a^\top k$. An upper loss occurs exactly at feasible $b<b^*$. The right derivative of the concave LP value is the smallest optimal multiplier. It is strictly positive to the left of $b^*$ by the chord bound and zero at or to the right. Equivalently, a zero multiplier lets the augmented dual certify the base optimum after deleting the new row. To see the one-positive-multiplier trap, maximize $x$ on $[0,1]$ and duplicate $x\le1$. Every split $(1-\lambda,\lambda)$, $0\le\lambda\le1$, of the dual unit cost is optimal. In particular the new row can receive multiplier one without moving the endpoint. A positive multiplier **can** certify tightening if its uniqueness or a positive lower bound for all optimal multipliers has separately been established.

For (c), separation of a point in $F\setminus N$ from closed convex $N$ provides a witness linear objective. This conclusion quantifies over all linear parameter objectives. An arbitrary separating direction need not be a reported cell-rate target. For (d), the slope of $\tilde D$ in $w_1$ is $-E[X(1-X)]/[\mu(1-\mu)]=-(1-\gamma)$. The four formulas are intersections of the base gap interval with the appropriate half-line, interpreting reversed endpoints as emptiness.

For (e), take the same square $[0,1]^2$, target $x$, and rows $x+y\le.8$ and $x\le.5$. Adding the first row alone lowers the upper endpoint to .8. Adding it last gives zero gain for $x$, yet still removes $(.5,1)$ and lowers the upper endpoint for $x+y$ from 1.5 to .8. **Correction.** Report “zero incremental gain for this target.” Use the inclusion test in (c) to label a constraint redundant for the feasible set.

### 7. R10 inherits a coverage theorem only when the augmented procedure has one

**Added rows and the slab.** In the conservative procedure there is one simultaneous event: every population row intercept is at most its inflated estimate, and the fitted target constant is within its two-sided allowance. All row normals and the target normal are fixed conditional on the design. On this event the true $k_0$ is feasible for the inflated LP and the projected interval contains the model-implied target. Adding finitely many affine rows changes the joint Gaussian quantile. Deterministic intercepts enter exactly. For a parallel disjunction, at least one piece holds at $(\beta_0,k_0)$. Inflating all pieces simultaneously and then maximizing dominates that piece, even at a tie. The common normal makes the resulting sample constraint a single LP row. This proof does not require strict feasibility. Taking the gradient of a winning estimated piece at a tie would be a different procedure.

**Union and emptiness.** Conditional on the design, the branch rows evaluated at the population regression and $k_0$ are fixed. Those evaluated at $\hat\beta$ are random. Choose the first true population branch in a fixed ordering if both hold. Its interval is contained in the reported union on every sample. This proves point coverage with each branch at the full error budget. It does not prove simultaneous coverage of both identified branch images or validate choosing a branch using $\widehat{ER}$. If designs form a triangular array and the true label switches with $n$, one needs branch validity along each subsequence on which that branch is true. There are only two such labels.

For (a), an empty inflated feasible set implies failure of the simultaneous event. For a candidate-based interval using the crossing rule, $\hat L>\hat U$ requires a one-sided failure about the same true **model-implied** target. A union bound gives $\limsup P({\rm empty})\le\alpha_U+\alpha_L$. A simultaneous projection at total coefficient-error level $\alpha_U+\alpha_L$ gives this bound directly. Widening for the realized-rate residual cannot create a crossing. Both branches empty implies the true branch empty. Conditional coverage given nonemptiness is $P({\rm cover})/P({\rm nonempty})$ whenever the denominator is positive. These arguments also work at a fixed point-identified law. Pointwise statements do not imply a uniform guarantee over sequences of nearly infeasible laws.

| Interval version | What R10 can claim |
|---|---|
| P3.5(a), simultaneous row inflation | (a)–(d), including the direct parallel slab, under a valid joint row-and-constant approximation for the design sequence. |
| P3.5(b1), all dual vertices, CLR | (a), (c), (d) under the full augmented CLR assumptions, including one-sided validity. The new rows change the dual feasible set, its vertices, gradients, selection and quantile. Simply appending row standard errors to the old CLR maximum is insufficient. |
| P3.5(b2)(i), a fixed simultaneous superset | The same branch-union and crossing arguments work with the fixed-superset bounds and allocated endpoint levels. |
| P3.5(b2)(ii), projection | The joint coefficient event controls both ends together, including any data-dependent collection of dual-feasible candidates. Use the coefficient-error budget once. Account separately for residual error in actual targets. |
| P3.5(b3), CLR calibrated on simulated candidates | The true-branch inclusion and crossing implications are deterministic. Their asserted nominal probabilities remain conjectural. Simulation coverage cannot promote them to theorems. |

The status notes at 5a lines 228–233 mostly make the right version distinctions. Two limits need to appear alongside them. First, P3 §0.5 lines 64–75 lists growing-district M1/M2 in regime R3 as conjectural even before candidate collection. The statement “R1/R3” at 5a line 218 does not settle increasing dimension. Second, the full vertex set being finite at each $n$ does not verify CLR Condition P along the sequence. The cited condition requires, among other things, a suitable fixed compact index representation, uniformly regular gradients bounded away from zero, well-conditioned covariance, and a covariance-estimation rate. Nonzero gradients for one design are insufficient for that citation. A zero-variance **candidate endpoint** is a separate issue from an exactly known row intercept. Check or separately treat deterministic candidate bounds. [CLR, Condition P, Lemma 4 and Theorem 4](https://arxiv.org/pdf/0907.3503v5)

For actual $B_{rc}$, carry over P3-S1a, the design weight condition and the residual allowance. For an actual gap, use its own residual calculation. For example the unit term $w_{i1}e_{i1}-w_{i0}e_{i0}$ has range length at most $w_{i1}+w_{i0}$. Thus $\tfrac12\{\sum_i(w_{i1}+w_{i0})^2\}^{1/2}$ is a valid standard-deviation envelope under conditional unit independence. The single-row $S_1$ need not bound this contrast. No allowance is needed for the model-implied gap in the supplied Monte Carlo.

**An analytic counterexample for (e).** Repeat each of $X=1/4,1/2,3/4$ equally often, with equal unit sizes. Take $m_0=.4,\ m_1=.55$ constant. Independently across units and groups, draw bounded group rates of these means with positive variance, for example beta draws of fixed positive concentration. All true base probability rows have strict slack, OLS has a full-rank fixed design and positive definite asymptotic covariance, and HC3 is consistent. The positive CR branch has $\tilde D_0=ER=.15$ on its upper edge; the other bracket end is $.025$.

Use the script's basis $((1-x)^2,2x(1-x),x^2)$ and skew coordinate $k$. Symmetry gives, for **any** coefficient vector $q$,
$$
ER(q)=q_{11}-q_{00},\qquad c_D(q)=\tfrac{7}{12}(q_{11}-q_{00}).
$$
Indeed $E[X^2]/E[X]=7/12$ and the middle coefficient cancels. Thus ${\rm se}(c_D)/{\rm se}(\widehat{ER})=7/12$ exactly. With probability tending to one the plug-in sign is positive and the base rows allow the $k$ attaining the plug-in upper edge, by strict slack and continuity. The naive upper endpoint is then $\widehat{ER}+k_{\rm crit}\,{\rm se}(c_D)$. Its calibration uses eight one-sided base rows and two target tails. The union bound gives the population Gaussian critical value $k_{\rm crit}\le z_{.995}$ at level .95, irrespective of their correlations. Consequently
$$
\limsup_n P(\tilde D_0\in{\rm CI}_{\rm naive})
\ \le\ \Phi\!\left(\tfrac7{12}z_{.995}\right)=0.933525<0.95 .
$$
This proves failure for a single fixed repeating design sequence, without extrapolating the three reported simulation sizes. It also survives the script's 4,000 critical-value draws. Each draw exceeds $z_{.996}$ with conditional probability at most .04 by the same union bound. An empirical .95 quantile above that threshold needs at least 200 exceedances. The limiting coverage is therefore at most $\Phi(7z_{.996}/12)+P\{\mathrm{Bin}(4000,.04)\ge200\}=.940092$.

The supplied far-edge design family also satisfies the sufficient condition in its large-sample design limit. Its unweighted composition is $X=.1+.8U$ for $U$ beta$(2,2)$, and the independent bounded sizes give the same limiting person-weighted distribution. Symmetry gives $ER(q)=q_{11}-q_{00}$ and $c_D(q)=(.5+2\operatorname{Var}X)(q_{11}-q_{00})$ in that limit. Here $\operatorname{Var}X=.032$, hence the limiting standard-error ratio is .564. The script calibrates the true gap to .15 at every design, its slopes tend to $+.2,-.2$ and its intercepts to $.4,.75$ for groups 0 and 1. The base rows therefore have limiting strict slack and the variance of the ER estimator is positive. Along convergent typical designs the same ten-tail bound gives limiting coverage at most $\Phi(.564z_{.995})$, about .927 with population critical values. Thus the persistence argument also applies to the documented design-calibrated experiment. The repeating three-point example avoids even that dependence of the structural coefficients on the design.

For the narrower procedure that directly intersects a confidence interval with the raw estimated bracket, the upper endpoint cannot exceed $\widehat{ER}$ with probability tending to one in this example. Its limiting coverage is at most one half. This direct argument avoids an unproved containment at 5a line 311. In general, dropping an empty pre-expansion intersection and expanding endpoints do not commute. With raw base interval $[0,.1]$, bracket $[.15,.2]$ and target allowance .1, the former intersection is empty although the expanded base intersects the bracket. Also distinguish the two sample ERs when the regression and target use different weights.

## 3. Minor comments and mapping to EGGHM

- **W1 versus W2.** Expanding $E[XT_c]$ group by group gives $MB+\Delta$, and $\mathbf1^\top\Delta=0$ since the composition sums to one. Subtracting $E[X_r\beta_r]$ from $E[X_rT]$ gives minus the sum of the pairwise between associations. The individual indicator of row $r$ has variance $\mu_r(1-\mu_r)$; the across-unit regressor has variance ${\rm Var}_p(X_r)$. Collapsing the other rows with their within-unit shares gives both W1 identities with these distinct denominators. Invertibility of $M$ is needed for $g$, only positive regressor variance for $ER^r$. Choose a bounded value for the absent rest group at $X_r=1$. W2 signs concern components of a **sum** of row-specific covariances. They need not sign every summand.
- **Residual order.** At 2c line 160, (S2) is a mean condition. It implies neither conditional independence of residuals given the entire design nor a $p^{-1/2}$ rate by itself. Even with independence the natural scale is the square root of the sum of squared design coefficients. A unit retaining positive weight can prevent concentration. State the sampling and weight conditions or delete the rate.
- **The ladder.** At 2c line 218, “both signs” implies CR when their specified directions agree. Two opposite signs generally do not. At lines 219–225, N2 admits steps whereas a classical local derivative is undefined at a jump. Specify a differentiable subclass or a precise almost-everywhere interpretation with its support conditions. State comparisons between realized and model-implied associations only after choosing the same population of tables or means. At $R\ge3$ specify a simplex chart and the share that adjusts. For example $m_s=X_3$ for every group is constant in $X_1$ holding $X_3$ fixed. Under a symmetric Dirichlet design $E[X_3\mid X_1]=(1-X_1)/2$; tilting by $X_1$ or $1-X_1$ preserves that conditional mean, giving negative covariance under both W1 weights. This supplies the association counterexample, beyond the unweighted covariance reported in C5(c). A globally increasing *univariate* function of $X_1$ still obeys the usual covariance inequality.
- **Identification and sample labels.** R3(e) is a population statement. Infeasibility of a raw fitted LP is a sample diagnostic. A level-controlled rejection uses R10's calibrated rule. R2 §2.5 line 208 needs its nonconstant-target qualification. Width-gain inference at line 210 needs joint treatment of the base and augmented endpoints, including the estimated threshold if there is one. An interval for either endpoint alone is insufficient.

| EGGHM v1 reference | Assessment of the mapping |
|---|---|
| §2, Lemma 1, Props. 1, 4, 6, 9–10 | Correct accounting and bias identities for bounded group means with individual population weights. ER is the slope of a population linear projection with intercept. EGGHM assumes positive mass in at least one mixed neighborhood (p. 8), giving $\gamma<1$. Positive composition variance is also needed for ER. Extending the algebra to purely segregated, varying designs is legitimate; the strict branch distinction at nonzero ER requires a mixed unit. |
| Props. 2–3, 5, 7, 8 and Cor. 3 | Correct MOB, signed and CR formulas on the continuous mean-probability model. Sharpness in v1 pp. 12–13 ranges over compatible joint distributions. That proof does not license integer-table interpolation at fixed known unit counts. The relation ${\rm sign}(ER-NM)={\rm sign}(ER)$ uses $\gamma<1$. |
| §4 ending, Theorems 1–3 | The one-vs-rest collapse is exactly the extension described on p. 20. W2 and simultaneous multiple-column LPs need their own proofs, which M1 and M4 supply. 1b line 69 attributes all three tables to Theorem 1; the gap is Theorem 1, the group levels are Theorems 2 and 3. |
| §5, Props. 11–13, Appendix Props. 14–16 | The derivative decomposition and same-sign bracket have the cited numbering. Neighborhood-specific interpretations in Props. 12–13 use equal group means at a given prevalence (Assumption 1) and meaningful derivatives. N2 imposes a common **conditional mean** and targets a behavioral contrast. It permits realized unit rates to differ at equal prevalence. Without unit-level homogeneity the appendix's prevalence-average interpretation is closer. Neither citation identifies a nonparametric derivative on an isolated atom. |

Source checked against the specified version. [Elzayn, Goldin, Guage, Ho and Morton, *Monotone Ecological Inference*, arXiv v1](https://arxiv.org/pdf/2504.14752v1)

## 4. Numerical checks run

I ran each of the five requested scripts from the pinned `docs/extension/` using `python3 -I`, with single-threaded BLAS. Runtime versions were Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. I compared the complete output of each run with its pinned text file.

| Script | Rerun and interpretation |
|---|---|
| `verify_phase2c.py` | **12/12 PASS.** C1 identity error $6.2\times10^{-15}$. C2 closed-form error $5.4\times10^{-15}$ against $5.6\times10^{-15}$ saved. Two positive-within restrictions were empty. C3 found 75 strict improvements among 139 retained feasible joint comparisons. C4 checks affinity for pairwise between terms and W2 in a pooled model. Its ordering construction changes one column without rebalancing the other columns, so it tests the scalar sign implication without always constructing a valid full transition matrix. C5(b) compares one projected target across 47 cases. C6 differentiates a chosen extension. None of these tests checks integer feasibility, district intercepts or derivative identification. |
| `verify_phase1_mei_did.py` | **10/10 printed checks PASS** (V7 has three parts). V1 and V4 check the covariance identities; V2 handles 1,093 nonempty and 107 empty cases with closed-form error $8.9\times10^{-16}$. V2 builds its sign rows from the bias identities themselves. V3 tests pooled lines, V5 the row-weighted monotonicity implication, and V6 finds a piecewise-linear witness, minimum about $-.0030$. Its segment constraints check one-sided derivatives at knots. V7–V8 concern the companion DiD algebra. Only V7b's rounding changed, $2.5\times10^{-16}$ versus $2.2\times10^{-16}$. |
| `verify_R3.py` | Pinned output reproduced. S1 attained 1,000 target values with worst violation $1.8\times10^{-15}$ using deterministic totals, evidence for **(b)**. S2 reported one positive-width branch in 2,000 designs. Its code labels a branch nonempty only if its width exceeds $10^{-9}$, and labels zero positive-width branches “both degenerate,” also allowing two empty sets. S3 found no both-feasible W1 cases. S4 found 357/1,000 and 614/1,000 both-feasible cases for the alternative pairings. S5 found 3 gaps in 3,600 target problems, 951 both-feasible and 945 common-witness certificates. Cases reuse designs across objectives. A frequency on these draws is descriptive. |
| `verify_R2.py` | T1 agreed in 800/800 cases, with a positive solver multiplier in 133 of 400 active unchanged-endpoint cases. T1b matched on all tested cuts over 2,907 positive-width faces. T2 matched 333/333 feasible cases and T3 matched all 7,409 eligible comparisons. **T4 differed:** this run gave 153 implied, 494 set-only and 1,153 target-tightening out of 1,800, versus 147, 494 and 1,153 out of 1,794 saved. T4 skips an unsuccessful base solve without distinguishing status; the old output does not record enough information to determine why six comparisons are absent. The substantive fraction stays near 28%. Tests use zero target intercept in the face LP and do not detect the line 169 error. |
| `verify_R10.py` | Pinned output reproduced. For the far positive edge at $n=150$, 600 replications gave JKST 1.000, slab .992, branch union .992, naive .917 (Monte Carlo SE about .011 for naive). The other designs gave 1.000 for all four and no empties occurred. Each design is held fixed across replications; the target is $\tilde D$, and both branch intervals use the conservative procedure. These runs do not evaluate CLR, residual coverage, the crossing rule, or a growing-dimensional regime. Lines 189–190 sum the two interval lengths even when they overlap, so the reported “union length” can exceed its Lebesgue length. Merge the intervals first. |

**Independent checks.** The appendix is a standalone script. Running it verified the boundary accounting identities, the integer and pure-design examples, strict joint-column improvement, the district reversal, the missing support-point and smooth derivative counterexamples, full-law versus regression sharpness, target-dependent holes, LCE refutation of CR, and the facet and redundancy examples. It uses an exact Gaussian experiment to stress R10's event inclusions at a point-identified target. Across 200,000 draws, the joint event had frequency .94943, coverage .99766, and emptiness .00094. A two-branch tie gave coverage .95078, and a positive-edge union .97557. The draws illustrate the event proofs; the persistent-failure bounds .933525 and .940092 in comment 7 use analytic inequalities.

## 5. Outside scope · background issues encountered

1. **2a §0.3 line 57, Lemma 1 lines 144–150, and Lemma 3(c) lines 310–314.** Equality of continuous functions almost surely *on* $A$ extends to the support of the measure restricted to $A$. This can be smaller than ${\rm supp}(P_X)\cap A$. For $X$ uniform on $[0,1]$ and $A=\{1/2\}$, any two functions agree almost surely on $A$ in that sense. The constants zero and one differ at $1/2$. This also exposes the tension between Lemma 3(c) and Remark 5.1. Define support using the restricted measure when trimming and impose A5 on that support. R3's full-mass hypothesis avoids the problem.
2. **P3 §5 lines 179–189.** The rank-$q$ pseudoinverse projection proof bounds gradients in the covariance range. A singular CLT alone need not control smaller errors in null directions selected as objectives. For instance, an estimation error $(Z/\sqrt n,-1/n)$ has limiting covariance ${\rm diag}(1,0)$; an upper bound for the second coordinate using zero standard error misses surely. Work in a full-rank reduced parameterization, require exact null-coordinate identities, or control the omitted components. This caveat concerns the conditions for projection, independently of CR.
3. **P3 §0.3 line 40 and §5 lines 201–203.** Small maximum squared weight relative to the sum of squared weights need not give a variance-standardized CLT if only one summand has positive variance. The residual tail bound with the conservative envelope $S_1$ can still be proved. Along subsequences on which ${\rm sd}(E)/S_1\to0$, use Chebyshev. If the ratio is bounded away from zero, the maximum-weight condition and bounded summands give Lindeberg relative to the actual standard deviation. This supplies the missing argument without strengthening the conclusion.
4. **1b §§A1–A2, lines 59–80 and 116–121.** Correct the theorem-to-table mapping noted above. For pure units, unit-level row ordering cannot be equated with the observable implication of a local between sign, since the between covariance is identically zero there. These points affect attribution and boundary cases.

## 6. Comparison with the earlier internal review

I agree with §4 on using one emptiness rule per version, anchoring the crossing proof at the model-implied target, conditioning on nonemptiness, inflating the slab piece by piece, and using a branch union for CLR. I also agree that S1 checks regression sharpness, that a common witness is only a sufficient certificate, that the general R2(c) argument needs separation, and that the larger-face test fills a genuine gap in T1. Combining two tests would give an upper bound from a union bound; their actual rejection probability need not double.

I differ from the overall “accept with fixes” assessment for the requested batch. The internal review covered 5a §§1–3 and three scripts. Its scope did not address the integer-table, district-slope, derivative-identification or binning problems in 2c. Within 5a, it missed the target intercept in the face LP, the distinction between set and target redundancy in (e), the empty LCE–CR example, and double-counting in the Monte Carlo union lengths. Its R10(e) sketch states the right variance-ratio condition but leaves its verification to simulation. Comment 7 verifies it in a fixed-design model. The inference labels should also preserve P3's growing-district qualification separately from its candidate-simulation conjecture.

## Appendix · Reproducible independent checks

Save the following block as `referee_checks.py` and run `python3 -I referee_checks.py`. It reads no repository or data files.

```python
"""Independent checks. Run with python3 -I referee_checks.py (NumPy and SciPy)."""
from fractions import Fraction as F
from itertools import product

import numpy as np
from numpy.polynomial import Polynomial
from scipy.optimize import linprog
from scipy.stats import norm

def cov(w, a, b):
    w = w / w.sum()
    return (w * (a - w @ a)) @ (b - w @ b)

def within(x, m1, m0):
    p = np.ones(len(x)) / len(x)
    mu = p @ x
    return mu*cov(p*x, x, m1) + (1-mu)*cov(p*(1-x), x, m0)

# M1. Direct weights with absent cells, including a pure first row.
X = np.array([[1,0,0],[0,.6,.4],[.3,0,.7],[.2,.5,.3]])
p = np.array([1,2,3,4])/10
m = np.array([[.8,.3,.1],[.2,.9,.5],[.3,.1,.6],[.1,.5,.9]])
T, mu = (X*m).sum(axis=1), p @ X
B = (p @ (X*m))/mu
M = X.T @ (p[:, None]*X)
W2 = sum(mu[s]*np.array([cov(p*X[:,s], X[:,r], m[:,s])
                         for r in range(3)]) for s in range(3))
assert np.allclose(X.T @ (p*T)-M@B, W2) and abs(W2.sum())<1e-12
for r in range(3):
    x = X[:,r]
    d = B[r] - (p@T-mu[r]*B[r])/(1-mu[r])
    other = np.divide(T-x*m[:,r], 1-x, out=np.zeros(len(x)), where=x!=1)
    W1 = mu[r]*cov(p*x,x,m[:,r])+(1-mu[r])*cov(p*(1-x),x,other)
    between = p @ (x*(m[:,r]-T))
    assert np.isclose(W1, cov(p,x,T)-cov(p,x,x)*d)
    assert np.isclose(between, mu[r]*(1-mu[r])*d-cov(p,x,T))
print("M1: all boundary-design identities hold")

# M2/M3. Enumerate individual counts for the stated three-unit example.
rows, yes = [1,2,3], [0,1,0]
legal = [range(max(0,y-(4-r)),min(r,y)+1) for r,y in zip(rows,yes)]
gaps = {F(sum(z),6)-F(sum(yes)-sum(z),6) for z in product(*legal)}
assert gaps == {F(-1,6),F(1,6)}
assert np.isclose(cov(np.ones(3),np.array(rows)/4,np.array(yes)/4),0)
x=np.array([0.,1.])
assert np.isclose(within(x,x,x),0) and np.all(x*(1-x)==0)
print("M2/M3: integer gaps -1/6,+1/6; CR empty. Pure units: both branches.")

# M4. Two 2x2 units, complements are the second column. Opposite column
# within signs jointly force D=ER. At T=.5 both individual cuts are feasible.
x = np.array([.25,.75])
# m0=(.5-x*m1)/(1-x); at these designs D = weights@m1 - 1.
weights = 2*x
def ends(A,b,bounds,c):
    opt = lambda s: linprog(s*c,A_ub=A,b_ub=b,bounds=bounds,method="highs")
    lo,hi=opt(1),opt(-1)
    assert lo.success and hi.success
    return np.array([lo.fun,-hi.fun])
bounds=[(0,1),(1/3,2/3)]
one=ends([weights],[1],bounds,weights)-1
both=ends([weights,-weights],[1,-1],bounds,weights)-1
assert np.allclose(one,[-.5,0]) and np.allclose(both,[0,0])
print("M4: joint within signs change [-.5,0] to {0}")

# M5. Two district intercepts and positive slopes in both rows.
x=np.array([.1,.2,.3,.7,.8,.9])
m_base=np.where(x<.5,.8,.2)+.1*x
aff=lambda k: within(x,m_base+k*(1-x),m_base-k*x)
assert aff(0)<0 and np.isclose(aff(.02),(aff(0)+aff(.04))/2)
print("M5: district-pooled within =",aff(0),"; affine in k")

# N2.1/N2.2. Extreme totals force both rates to be the total.
b=np.array([[0,0],[1,1]])
assert np.allclose(b.mean(0),[.5,.5]) and np.isclose(np.mean([.75,.25]),.5)
x=np.array([.25,.5,.75])
tau=x*.4+(1-x)*np.array([8/15,.7,.8])
assert np.allclose(tau,[.5,.55,.5])
print("N2: middle composition prevents constant endpoint extension")

# N2.3. A polynomial vanishing on the design changes the derivative.
z=Polynomial([0,1])
h=(z-.25)*(z-.5)*(z-.75)
assert np.allclose(h(x),0) and np.isclose(h.deriv()(.5),-1/16)
# Uniform X on [1/5,2/5]. Integrate each size-weighted covariance.
def E(f):
    integ=f.integ()
    return 5*(integ(.4)-integ(.2))
def weighted_term(w,m):
    return E(w*z*m)-E(w*z)*E(w*m)/E(w)
m1=z-z*z/2
m0=1-z*z/2
assert np.allclose((z*m1.deriv()+(1-z)*m0.deriv()).coef,0)
delta=weighted_term(z,m1)+weighted_term(1-z,m0)
assert np.isclose(delta,-1/47250,atol=1e-14)
print("N2.3: same support values, derivative shift -1/16; smooth local zero, aggregate",delta)
# Section 6. Check row-specific bin averages.
x=np.array([.25,.75])
assert np.isclose(np.mean(x*x),.3125) and np.isclose(x.mean()**2,.25)
assert np.isclose((x@x)/x.sum(),.625)

# R3(a,b). Full-law T in {0,1} fixes both means; T'=1/2 allows skew k.
x=np.array([.25,.5,.75])
for k in [.2,.3,.4]:
    m1,m0=.5+k*(1-x),.5-k*x
    assert np.allclose(x*m1+(1-x)*m0,.5)
    assert np.isclose(x@m1/x.sum()-(1-x)@m0/(1-x).sum(),5*k/6)
# R3(c). Two separated vertical strips, CR pair (u+.5,u-.5).
for target, expected in [(np.array([1,0]),[[-1,-.5],[.5,1]]),
                         (np.array([0,1]),[[-1,1],[-1,1]])]:
    pieces=[ends([[s,0]],[-.5],[(-1,1)]*2,target) for s in [1,-1]]
    assert np.allclose(pieces,expected)
# R3(d,e). Forced skew because off-diagonal probabilities sum to zero.
x=np.array([0.,.25,1.])
m1,m0=(49/55)*x,1-x
tau=x*m1+(1-x)*m0
d=x@m1/x.sum()-(1-x)@m0/(1-x).sum()
assert abs(cov(np.ones(3),x,tau))<1e-14 and np.isclose(d,-261/1925)
print("R3: sharp empty vs outer [1/6,1/3]; CR gap depends on target; LCE excludes CR")

# R2. Optimal facet with target intercept 2 and cuts that keep/remove it.
x_target=np.array([1,0])
for rhs,expect in [(1.2,3),(.8,2.8)]:
    assert np.isclose(2+ends([[1,1]],[rhs],[(0,1)]*2,x_target)[1],expect)
# Zero last-added gain for x, positive last-added gain for x+y.
for c in [x_target,np.ones(2)]:
    u0=ends([[1,0]],[.5],[(0,1)]*2,c)[1]
    u1=ends([[1,0],[1,1]],[.5,.8],[(0,1)]*2,c)[1]
    assert np.isclose(u0-u1,0 if c[1]==0 else .7)
print("R2: partial facet survives; full facet removed; target-specific redundancy")

# R10(a-d). Simulate an exact Gaussian experiment at ties, including emptiness.
rng=np.random.default_rng(831)
z1,z2,z3=rng.standard_normal((3,200000))
from scipy.optimize import brentq
q=brentq(lambda u:norm.cdf(u)**2*(2*norm.cdf(u)-1)-.95,1,4)
# k >= beta1 and k <= beta2; beta1=beta2=0, target k+beta3.
lo,hi=z1-q,z2+q
Egood=(z1<=q)&(-z2<=q)&(np.abs(z3)<=q)
cover=(lo<=hi)&(lo+z3-q<=0)&(0<=hi+z3+q)
assert np.all(cover[Egood])
# Parallel row k <= max(beta1,beta2); add k>=0, true k=0.
slab=np.maximum(z1+q,z2+q)
assert np.all(slab[(z1>=-q)&(z2>=-q)]>=0)
# W1 toy slab at beta=0; plus and minus branches calibrated separately
# with +/-Z directions, critical value z_.975. Empty if the ends cross.
t=norm.ppf(.975)
gamma=.25
pl,pu=gamma*(z1-t),z1+t
ml,mu=z1-t,gamma*(z1+t)
union=((pl<=0)&(0<=pu))|((ml<=0)&(0<=mu))
assert np.array_equal(union,abs(z1)<=t)
# Away from the tie, scaled true beta=4 and k=4 is on the true plus edge.
v=z1+4
cp=(gamma*(v-t)<=4)&(4<=v+t)
cm=(v-t<=4)&(4<=gamma*(v+t))
assert np.all((cp|cm)[cp])
print("R10: simultaneous, cover, empty, tie union, positive-edge union:",
      *[round(float(v),5) for v in [Egood.mean(),cover.mean(),(lo>hi).mean(),
                                     union.mean(),(cp|cm).mean()]])
# R10(e). Symmetric fixed design, fixed population edge, eight base rows
# and two target tails. The marginal-quantile union bound is deterministic.
x=np.array([.25,.5,.75])
design=np.array([(1-x)**2,2*x*(1-x),x*x]).T
p=np.ones(3)/3
er=np.array([cov(p,x,f)/cov(p,x,x) for f in design.T])
g=np.array([-((1-x)**2).sum()/(1-x).sum(),0,(x*x).sum()/x.sum()])
assert np.allclose(er,[-1,0,1]) and np.allclose(g,(7/12)*er)
print("R10(e): limiting naive coverage at most",norm.cdf((7/12)*norm.ppf(.995)))
from scipy.stats import binom
finite_draw_upper = norm.cdf((7/12)*norm.ppf(.996)) + binom.sf(199,4000,.04)
assert finite_draw_upper < .95
print("R10(e): bound with 4000 critical-value draws",finite_draw_upper)
```
