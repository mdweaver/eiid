# New avenues: what else in the data can tighten the bounds?

*A Phase-1 research proposal. It surveys ways to extract more identifying information from ecological data, states what each one assumes, and ranks them. Nothing here is a finished derivation. Every claim is labelled **proved** (with a short argument that Phase 2 will write out with tags), **conjecture**, or **heuristic**. Verification scripts and their outputs are in `docs/extension/new_avenues/`. Notation follows P1 and 2a: $H_c=Q_c+K_c$; $Q$ is identified; $K\in\mathcal K$ is not; $m_{rc}(x)=e_r^\top H_cx$ is the conditional-mean rate.*

*Out of scope, because other sessions cover them: monotone EI, difference-in-differences and other cross-period restrictions, and links to Bayesian EI (Imai–Lu–Strauss, `eco`, Stan). Where an avenue touches one of these, the boundary is noted.*

---

## 0. Summary

**The diagnosis.** The data enter the current machinery in two ways:

* through the quadratic regression $E[T\mid X]$, which identifies $Q$; and
* (sharp set only) through the conditional law of $T$ given $X$, via Duncan–Davis kinks.

Every avenue below finds information that **neither of these uses**, or imports some from outside. There are four kinds:

1. **Unit size.** Sampling noise in a poll of $n$ voters shrinks like $1/n$. How fast it shrinks reveals how different the rows are *within* polls. (Avenue A)
2. **Context measured at another scale.** Any contextual variable other than the poll's own composition has *row-specific* effects that are identified. Only own-composition effects suffer from the "who is behind the pattern" problem. (Avenue B)
3. **Known tables at another level.** In New Zealand, for example, the electorate cross-tabulation is published. (Avenue E)
4. **Several contests sharing the same rows.** (Avenue C)

Three structural facts frame the search:

* **(proved)** joint ballot-level margins across contests add nothing without a cross-contest assumption;
* **(proved)** a known aggregate table pins only $(R-1)(C-1)$ of the $(C-1)R(R-1)/2$ unknown directions;
* **(proved)** independent poll-level shocks, however well deconvolved, never reveal $K$.

**Ranking** (details in §9):

| # | Avenue | New information used | Extra assumption | Identifying power | Estimability today |
|---|---|---|---|---|---|
| A | Unit size and second moments | $\operatorname{Var}(T\mid X,N)$ and its $1/N$ slope | binomial voting within poll-row cells (A1-rung); poll shocks not negatively correlated across rows (A2-rung) | **large**: one endpoint moves to the truth (A1), near-point (A2) | **weak**: needs many small units; plug-in is badly biased |
| B | Context at a second scale | a neighbourhood composition $z\neq x$ | own-poll effects bounded relative to the identified neighbourhood effects (dial $M$) | **large**, smooth in $M$ | good (few constraints) |
| E | Known tables at another level | published cross-tabs (e.g. NZ electorate tables) | the model holds within the area whose table is known | large for sub-area targets; point-identifies $K$ when $R=2$ | good |
| C | Several contests | other contests on the same ballot | cross-contest structure (contests as "districts", or single index) | medium–large | data-dependent |
| D | Polynomial (Bernstein) ladder | — (robustness, not tightening) | weaker than LCE | negative: widens, by design | good, and an exact LP |
| F | Computing the sharp set | $P_{T\mid X}$ (already pinned) | none | small–medium | good |
| G–J | Symmetry, deconvolution, entropy/OT, Choquet | — | various | small or none | — |

**Main messages.**

* **A is the most interesting avenue scientifically, and it is not a new assumption to most of the EI literature.** Every likelihood-based EI model (King, Wakefield, Brown–Payne, Greiner–Quinn, Imai–Lu–Strauss) already assumes binomial or multinomial voting within poll-row cells.
* **A works in New Zealand.** In 31,359 NZ polling places the $1/n$ term is clearly present and binomial-sized: the slope on its binomial prediction is 0.69–0.81, with SE about 0.10 (§1.6).
* **But A is statistically expensive.** The information sits in the smallest units, and naive plug-in bounds are badly biased inward.
* **B is the most practical new rung.** It uses one more variable that analysts usually have, and its dial $M$ fits the existing sensitivity machinery.

---

## 1. Avenue A: unit size and the second moments of $T$ given $X$

### 1.1 Plain language

* **What it is.** In a poll of 30 voters the outcome shares bounce around much more than in a poll of 3,000. The amount of bouncing depends on how different the voting groups in the poll are from each other.
  * If every group voted P with probability 0.5, small polls would be very noisy.
  * If in-party voters almost always vote P and others almost never do, a small poll's outcome is nearly fixed by its composition, so there is little noise.

  So **comparing how noisy small and large polls are reveals how far apart the groups' rates are.** That spread is exactly what the regression cannot see.
* **Key assumption.** Within a poll, voters of the same group vote independently given that poll's group-level propensities. Formally: binomial or multinomial voting, and the poll-level propensities do not depend on poll size.
* **If it is wrong.** Correlated voting within a poll, e.g. households or a local get-out-the-vote drive, makes small polls *noisier* than binomial. The method then reads that noise as "groups are similar", and the bound becomes too narrow, possibly excluding the truth. Small polls that differ systematically (hospitals, rural booths) also bias it.
* **What can be checked.**
  * The binomial slope can be compared with its binomial prediction; slopes above 1 signal correlated voting (§1.6).
  * The constraint must leave a non-empty set.
  * The size–variance relationship can be inspected directly.

### 1.2 The information we currently discard

The regression uses $E[T\mid X]$. The sharp set $\Theta_I$ uses all of $P_{T\mid X,N}$, but only through the realised unit polytopes $S(X,T)$, i.e. under the *realised-share* model. Neither uses the fact that $T_i$ is an average of $n_i$ individual choices. Under the propensity model that fact gives a variance decomposition with a known dependence on $n$.

### 1.3 Three rungs (decreasing plausibility)

Fix a column $c$ and drop it from the notation. Write $\tau(x)=E[T\mid x]$ and $W(m;x)=\sum_rx_r\big(m_r(x)-\tau(x)\big)^2$, the $x$-weighted *between-row* spread of the mean rates. $W$ is a convex quadratic function of $K$.

**A0. No new assumption (realised shares, A3 only).** For every $(x,n)$ in the support:
$$\textstyle\sum_rx_r\,m_r(x)^2\ \le\ \tau(x)-\operatorname{Var}(T\mid X=x,N=n). \tag{A0}$$

* *Status:* **proved**. Sketch:
  1. $W_i:=\sum_rx_r\beta_{ir}^2-T_i^2\le\sum_rx_r\beta_{ir}-T_i^2=T_i(1-T_i)$, since $\beta^2\le\beta$ on $[0,1]$.
  2. $E[W_i\mid x,n]=W(m;x)+\sum_rx_r\operatorname{Var}(\beta_r)-\operatorname{Var}(T)$.
  3. $\operatorname{Var}(T)=\operatorname{Var}(\sum_r x_r\beta_r)\le\sum_rx_r\operatorname{Var}(\beta_r)$, by Jensen with weights $x_r$.
  4. Combining: $W(m;x)\le E[T(1-T)\mid x,n]$, which rearranges to (A0).
* It is implied by $\Theta_X$, so it is a cheap *outer approximation* of the sharp set, not new identification.
* In the simulation it never bound (§1.5): realised-share variance in polls of 20 or more voters is too small.

**A1. Binomial voting.** Assume:
* (i) given poll-row propensities $\pi_{ir}$, the $N_{ir}$ voters choose independently;
* (ii) the law of $\pi_i$ given $(X,N)$ does not depend on $N$. This is A3 strengthened from the mean to the law; only the first two moments are needed.

Then $\operatorname{Var}(T\mid x,n)=A(x)+B(x)/n$ exactly, with

* $A(x)=\operatorname{Var}(\sum_rx_r\pi_r\mid x)$, and
* $B(x)=\sum_rx_rE[\pi_r(1-\pi_r)\mid x]$.

Both are identified from variation in $n$ at fixed $x$. Extrapolating to $n=1$ gives
$$W(m;x)\ \le\ \omega(x):=\tau(1-\tau)-A(x)-B(x). \tag{A1}$$

* *Status:* **proved**. $\omega(x)=E[\sum_rx_r(\pi_r-\tau_i)^2\mid x]$, the expected within-poll between-row variance of propensities, and $W(m)\le\omega$ by the same Jensen step as A0.
* **Reading.** A0 evaluated at a poll of one voter. Since $\omega-W(m)=\sum_rx_r\operatorname{Var}(\pi_r)-A(x)$ is the poll-level heterogeneity, the constraint is **tight at the truth when polls are homogeneous given $x$**.
* It caps how *different* the rows can be. In 2×2 it reads $|m_1(x)-m_0(x)|\le\sqrt{\omega(x)/(x(1-x))}$.
* With the ordering $m_1\ge m_0$ (a pinned modifier), this is a **linear row in $\kappa$**, so the LP machinery applies unchanged.
* In R×C it is a convex quadratic (an SOCP row) per $(x,c)$. The cross-column version is a matrix inequality,
  $$\textstyle\sum_rx_r(m_{r\cdot}-\tau)(m_{r\cdot}-\tau)^\top\preceq\Omega(x),$$
  an LMI, with $\Omega$ identified from the $1/n$ slope of $\operatorname{Cov}(T\mid x,n)$. Polyhedral outer approximations by tangent cuts keep everything an LP.

**A2. Poll shocks are not negatively correlated across rows.** Assume additionally $\operatorname{Cov}(\pi_r,\pi_s\mid x)\ge0$. Poll-level shocks such as candidate popularity or local swing move the groups together. This is the same flavour as the "non-negative correlation" option agreed for the unit-level term (P4, second release).

Then the heterogeneity is bounded by the identified $A(x)$. In 2×2,
$$(m_1-m_0)^2\ \ge\ \frac{\omega(x)}{x(1-x)}-\frac{A(x)}{\min(x,1-x)^2},$$
a **lower** bound on the spread, so both sides are now bounded. With ordering it is again a linear row.

* *Status:* **proved** in 2×2. The steps are $\operatorname{Var}(\pi_1-\pi_0)\le V_1+V_0$ and $A\ge x^2V_1+(1-x)^2V_0$. The R×C version is a **conjecture**: a matrix analogue via $\operatorname{diag}$ dominance.

**A3 (extension). Higher factorial moments.**

* At fixed $x$, the $k$-th factorial moment of the counts is a polynomial in $n$ of degree $k$. Its $k$ non-constant coefficients identify $k$ linear functionals of the $k+1$ degree-$k$ mixed moments $E[\pi_1^j\pi_0^{k-j}\mid x]$. One direction per degree stays free, mirroring the single free $\kappa$ (**proved** for 2×2 by direct expansion; this is the $k$-th-order version of Plackett's ambiguity).
* Requiring these to be moments of *some* law on $[0,1]^2$ is a Hausdorff moment problem: positive semidefinite moment and localising matrices, an SDP (Lasserre hierarchy).
* **Conjecture:** in the no-heterogeneity limit, A1 plus x-variation point-identifies $K$ up to a finite set for every $R$. For $R=2$ this is **proved**: $(g(x)+\kappa)^2=D(x)$ for all $x$ has a unique solution when the regression's $x^2$ coefficient is non-zero, and a sign ambiguity otherwise. The ambiguity is resolved by ordering.

### 1.4 Ladder position

* A0 is implied by the existing model.
* A1 and A2 are **modifiers that apply to any rung** (M0, M1, …). They are nested: $\Theta\cap A2\subseteq\Theta\cap A1\subseteq\Theta$.
* They restrict the *within-poll* behaviour of individuals, not contextual effects, so they are orthogonal to the existing ladder.
* They are not refutable from $E[T\mid X]$, but they are partly testable (§1.6).

### 1.5 Verification (`new_avenues/sim_second_moments.py`)

**Design.**
* 2×2, with true loyalty $0.75+0.15x$ and crossover $0.05+0.20x$.
* Poll propensities are Beta with dispersion $\phi$ (heterogeneity), correlated across rows through a Gaussian copula with correlation $\rho$.
* Votes are binomial, and poll sizes are NZ-like (median 300, log-sd 0.9).
* Target: model-implied loyalty, with truth 0.840.

**Identification (population moments).** These are bounds on loyalty, and the pattern is the same in all four designs, including large heterogeneity:

| constraint set | bound | width / width(O) |
|---|---|---|
| O (outer set, current package) | [0.521, 0.941] | 1.00 |
| O + A0 (polls of ≥ 20) | [0.521, 0.941] | 1.00 |
| O + A1 | [0.521, **0.841**] | 0.76 |
| O + ordering | [0.578, 0.941] | 0.86 |
| O + A1 + ordering | [0.578, **0.841**] | 0.62 |
| O + A1 + A2 + ordering | [0.839, 0.841] | ≈ 0.00 |
| high heterogeneity ($\phi=8$), O + A1 + A2 + ordering | [0.832, 0.845] | 0.03 |

**Estimation.**

* **The pointwise plug-in fails.** It takes a minimum over a grid of noisy constraints, which is the Phase-3 inward-bias phenomenon at its worst. It covers the truth in only 4–17% of replications and is often empty.
* **A valid integrated version** averages the constraint over $x$ and moves its right-hand side by 1.645 standard errors. It covers (93–98%), but buys almost nothing until $p$ is in the thousands:

  * **O + A1 + ord (valid):** the O + ordering width (0.36) at 64 and 500 units in every design. At 5,000 units, widths are 0.30–0.31 in the two low-heterogeneity, NZ-size designs and there is no gain otherwise.
  * **O + A1 + A2 + ord (valid):** in the two low-heterogeneity, NZ-size designs, widths of 0.18 at 500 units and 0.08–0.09 at 5,000. Coverage is only 0.74–0.80 at 500 and 0.88–0.91 at 5,000.
  * **High heterogeneity or large units (median 3,000):** little or no gain even at 5,000 units.

**Honest summary.** The identifying power is large. It is concentrated in the smallest units, and exploiting it needs far more units than one electorate provides: pooling across many tables, or national precinct files. It also needs a proper moment-inequality treatment. The A2 interval's coverage below 95% at moderate $p$ shows that the simple correction is not enough.

### 1.6 Real-data feasibility check (`new_avenues/nz_diagnostics.py`, D1)

**Data.** 450 NZ 3×3 tables with 31,359 polling places: median size 262, and 29% have fewer than 100 voters.

**Method.** Leverage-corrected squared residuals of the within-table quadratic fit are regressed on:
* table fixed effects;
* a quadratic in $x$ (for $A$); and
* $b_i/n_i$, where $b_i=\sum_rx_{ir}B_{rc}(1-B_{rc})$ is built from the true aggregate rates.

**Result.** The slope is **0.69 (LAB), 0.81 (NAT), 0.69 (OTH)**, with cluster SE about 0.10.

**Reading.** The $1/n$ term is unmistakably present, and smaller than the binomial ceiling. A slope below 1 is what binomial voting with some poll heterogeneity predicts. There is **no sign of over-dispersion** from correlated voting, the main threat to A1. This is evidence of feasibility, not a test of A1 itself.

### 1.7 References

* **Plackett (1977, *Biometrika*), "The marginal totals of a 2×2 table".** Identification from the distribution of a sum of two binomials. *Citation to verify.*
* **Wakefield (2004, *JRSS A* 167), "Ecological inference for 2×2 tables".** The convolution likelihood, and information from small areas.
* **Brown & Payne (1986, *JASA* 81), aggregated-data model for voting transitions.** Dirichlet-multinomial over-dispersion used in estimation. The closest precedent; *verify volume and pages*.
* Lasserre (2001, *SIAM J. Optim.*) and Schmüdgen / Putinar for the moment-problem SDP.
* Bhatia & Davis (2000, *Amer. Math. Monthly*).

**What is new here** relative to those references:
* an inequality (partial-identification) version that needs no Dirichlet or other parametric law;
* its integration with the LP ladder; and
* the two-sided A2 rung.

**Boundary with the Bayesian-link session:** these are frequentist moment inequalities, not a likelihood.

---

## 2. Avenue B: context measured at a second scale

### 2.1 Plain language

* **What it is.** Suppose loyalty depends on the political make-up of the *neighbourhood* (from adjacent polls, the census or the riding) as well as on the poll itself. The data can tell apart *which group* responds to the neighbourhood: in-party loyalty or out-party crossover. They can never do this for the poll's own make-up.
* **Using it.** We then assume that the poll's own make-up matters no more than $M$ times as much as the neighbourhood's, for each group. $M$ is a dial.
* **Key assumption.** Own-poll context effects are bounded by a multiple of neighbourhood context effects. At the extremes:
  * $M=0$: the poll's own composition doesn't matter beyond the neighbourhood (Goodman with neighbourhood context);
  * proportional effects: context works the same way at both scales, only more or less strongly.
* **If it is wrong.** If own-poll effects are larger than $M$ times the neighbourhood effects, the bound can exclude the truth. The simulation shows exactly that failure at $\lambda=1.5$, $M=1$.
* **What can be checked.**
  * "No own effect" is testable: in 2×2 the $x^2$ coefficient must be zero; in the simulation it was rejected (empty set), as it should be.
  * Proportionality is testable in R×C.
  * The *identified* half of the own effect can be compared with the identified neighbourhood effect, which is evidence about the size of $M$.

### 2.2 Mathematical idea

Let $z_i\in\Delta^{R-1}$ be a second composition: the neighbourhood, the riding, or census blocks. Consider
$$E[\beta_{rc}\mid x,z]=e_r^\top H_cx+e_r^\top\Gamma_cz\quad\Longrightarrow\quad E[T_c\mid x,z]=x^\top Q_cx+x^\top\Gamma_cz.$$

**Proposition B1 (proved; sketch).** If the support of $(x,z)$ is rich enough, then:

* the bilinear form identifies $\Gamma_c$ up to $\Gamma_c\mapsto\Gamma_c-v\mathbf 1^\top$ together with $H_c\mapsto H_c+v\mathbf 1^\top$. This is the same function, because $\mathbf 1^\top x=\mathbf 1^\top z=1$, so it is just an intercept normalisation;
* hence the row-specific neighbourhood effects $\partial_zE[\beta_{rc}\mid x,z]$ are **identified**;
* only $K$, the own-composition skew part, is not.

This generalises Lemma D1 (district intercepts identified) and JKST App. C (covariate coefficients identified). **The "who is behind the pattern" problem is confined to the poll's own composition.**

The precise rank condition is "no non-zero $(S,G)$ with $S$ symmetric, $x^\top Sx+x^\top Gz\equiv0$ on the support, other than the normalisation". **Conjecture:** it holds whenever $z$ is not an affine function of $x$ on the support.

### 2.3 Rungs (decreasing plausibility)

* **B0.** $M_0+\Gamma z$: weaker than M0, which is the case $\Gamma=0$. The bound may widen.
* **B1(M), relative magnitudes.** $|H$-effect of row $r|\le M\,|\Gamma$-effect of row $r|$, entry-wise in a chosen parametrisation (e.g. slopes in 2×2). These are LP rows, nested in $M$. The idea is the EI analogue of Rambachan & Roth's (2023) relative magnitudes: the unidentified component is calibrated by an identified one.
* **B2, proportionality.** $H_c^{\text{(slopes)}}=\lambda\,\Gamma_c^{\text{(slopes)}}$.
  * The symmetric part identifies $\lambda$, and then $K=\lambda\,\operatorname{skew}(\Gamma)$ is **point-identified** (**proved** in 2×2: $\lambda=d/(g_1-g_0)$).
  * In R×C, one scalar $\lambda$ must reconcile several identified symmetric entries per column, so B2 is overidentified and testable. The exact count, after the intercept normalisation, is for Phase 2.
* **B3, no own effect.** A point, refutable.

### 2.4 Verification (`new_avenues/sim_two_scale_context.py`)

**Design.** 2×2; 93 neighbourhoods × 10 polls; $z$ = the neighbourhood's mean poll composition. Truth: own slopes are $\lambda\times$ neighbourhood slopes.

**Case $\lambda=0.4$** (assumption $M\ge0.4$ true). Truth is 0.727; DD is [0.31, 0.88].

| rung | identification | width | estimated: width / coverage (plug-in) |
|---|---|---|---|
| B0 (O) | [0.495, 0.805] | 0.31 | 0.33 / 1.00 |
| B1(2) | [0.637, 0.805] | 0.17 | 0.18 / 0.97 |
| B1(1) | [0.693, 0.805] | 0.11 | 0.10 / 0.84 |
| B1(0.5) | [0.722, 0.778] | 0.06 | 0.05 / 0.30 |
| B3 | empty (refuted) | — | always empty |
| B2 | point 0.727 | 0 | median abs. error 0.04, but 5–95% error range [−0.54, 0.39] |

**Case $\lambda=1.5$.** B1(1) and B1(0.5) **exclude the truth**, as they should: the assumption is false. B1(2) covers it.

**Reading.**
* The gain is large and degrades smoothly in $M$.
* Plug-in coverage drops as constraints become nearly binding, which is expected; Phase 3's CLR applies directly, since there are few constraints.
* The proportional point is a ratio of estimated curvatures, so it is **weakly identified with heavy tails**. It should be reported only as a diagnostic or with weak-IV-robust inference.

### 2.5 Data needs and relation to existing rungs

* $z$ must vary *within* the groups used for M1. Otherwise the riding dummies absorb it.
* Candidates: k-nearest-poll composition, census dissemination areas (BC), meshblocks (NZ).
* $z$ that includes poll $i$ itself is fine, as long as $z\neq x$.
* Not the same as M1′ (riding covariates moderating $H$). That is a cross-group restriction; this is an exclusion or magnitude restriction at the unit level.

### 2.6 References

* Achen & Shively (1995), *Cross-Level Inference*, Chicago. Contextual models with distinct regressors. *Verify the chapter.*
* King (1997), *A Solution to the EI Problem*, ch. 9: covariates in EI.
* Rambachan & Roth (2023, *REStud*).

---

## 3. Avenue C: several contests sharing the same rows

### 3.1 Plain language

On one ballot people vote in several contests: president, senate, referenda. Every contest is a separate EI problem with the same rows.

* **On its own,** this helps not at all. Even ballot-level records showing how each contest's votes combine do not help without an assumption tying the contests together.
* **With such an assumption** it can help a lot. For example:
  * "context shapes all contests in the same way, up to scale", which treats contests like the districts of M1; or
  * "each group's support across many referenda lines up on one left–right dimension".
* **If wrong:** a realignment specific to one contest (e.g. a candidate with cross-party appeal) is read as structure.
* **Checks:**
  * the symmetric parts must satisfy the same cross-contest relation (testable);
  * the intersection must be non-empty;
  * with many referenda, a unit-level overidentification test.

### 3.2 Results and ideas

**C1 (proved; verified in `check_structural_claims.py`, 300 random tables, max difference 2e−16).** Joint margins (cast-vote records, CVRs) of contests 1 and 2 do not tighten the Duncan–Davis bound for a contest-1 cell.

*Reason:* fix the $(r,y_1)$ margin. For each $y_1$, a non-negative $R\times C_2$ array with row sums $n_{\cdot y_1}$ and column sums $t_{12}(y_1,\cdot)$ exists, because the totals match: a transportation polytope is non-empty iff the totals agree.

The same argument shows that a third election's one-way margin does not tighten a two-election transition table. So **gains need a cross-contest assumption**, and with only margins the data are a product.

**C2. Contests as districts (linear).** Suppose $H^{(j)}=a_jJ+b_jH^{(1)}$, i.e. context works the same way up to level and scale.

* Then $Q^{(j)}=a_jJ+b_jQ^{(1)}$, which is **testable**, and $K^{(j)}=b_jK^{(1)}$.
* The bound is $\bigcap_j\{K:\ Q^{(j)}+b_jK\ \text{feasible}\}$.
* This is the M1 intersection mechanism with contests in place of ridings. Contests with extreme baselines (a lopsided referendum) pin $K$ for all.
* LP; **proved** by the same argument as 2b Thm 2.
* *Caveat:* a common *linear* factor across contests cannot create information about $K$ on its own, because each contest reveals only its own symmetric part. The gain comes solely from the intersection of feasibility regions.

**C3. Single index (nonlinear).** Suppose $\beta^{(j)}_{ir}=F_j(\alpha_{ir})$, with $F_j$ monotone and known up to identified parameters: one latent position per poll-group and contest-specific cut-points.

* With $J\ge R$ contests, each unit gives $J$ equations in $R$ unknowns.
* **Conjecture:** generic unit-level point identification when $J\ge R$ (Jacobian rank), and $J-R$ overidentifying restrictions per unit when $J>R$.
* A partial-identification weakening, "the ordering of groups is the same in every contest", is a disjunction of LP rows given an ordering.

**C4. Tensor / Kruskal.** With $\ge3$ contests, CVRs and local independence of contest choices given (poll, group), the poll-level $\beta$'s are identified up to label switching.

* Known $x$ resolves the labels (Kruskal 1977; Allman, Matias & Rhodes 2009, *Ann. Statist.*).
* **Implausible for partisan contests:** partisanship correlates choices within a group. It is listed for completeness.

**Verification.** C1 is numeric (above). C2–C4 need data with several contests per ballot; none is in the repository. A simulation is proposed for Phase 2.

**Data.** US CVRs (e.g. the cast-vote-record collections assembled by Kuriwaki and co-authors; *citation to verify*), Swiss multi-referendum days, and California propositions. BC has a single contest per ballot. NZ's 2020 referenda were reported at polling-place level (*to verify*).

---

## 4. Avenue D: a polynomial ladder with exact LP certificates (weaker rungs)

### 4.1 Plain language

LCE assumes straight lines. We can allow gentle curves (quadratic, cubic) instead. This *widens* the bound: it is a robustness rung above L1, not a tightening. What is new is a cheap and exact way to compute it.

* **Key assumption:** contextual effects are polynomials of degree $k$ on the hull.
* **If wrong:** the same as L1, but less likely.
* **Check:** the degree-$(k+1)$ regression fit.

### 4.2 Content

**Dimension (proved, by counting).** With degree-$k$ homogeneous contextual effects, the unidentified dimension is
$$d_k=(C-1)\Big[R\binom{R+k-1}{k}-\binom{R+k}{k+1}\Big].$$

* Check: $d_0=0$ (Goodman) and $d_1=(C-1)R(R-1)/2$ (LCE). For $k=2$ and $R=2$, $d_2=2(C-1)$.
* The argument: the map $(m_r)_r\mapsto\sum_rx_rm_r$ from $R$-tuples of degree-$k$ forms onto degree-$(k+1)$ forms is surjective, so the kernel dimension is the difference.

**Exact LP feasibility.** "$m(x)\in$ simplex on $A$" becomes polynomial non-negativity on a polytope.

* On the simplex, non-negative Bernstein coefficients certify it, and degree elevation makes the certificate exact in the limit for strictly positive polynomials (Pólya 1928; Powers & Reznick 2001).
* On a general polytope $A$, Handelman's theorem plays the same role.
* All of these give **LP rows**, so the package needs no new solver.
* Monotonicity and convexity are sign conditions on Bernstein differences, also LP rows. This makes pinned item 8 (shape-restricted effects) computable.

**Conjecture.** As $k\to\infty$ with no smoothness bound, the bound tends to the model-free (DD) bound. With a bound on the $(k+1)$-th derivative, it interpolates. This gives a **curvature dial** that is more interpretable than P3's $\sup\lVert\Delta\rVert$.

**Already pinned (L1′, item 8).** What is new here is the dimension formula, exactness via Bernstein/Handelman LPs, and the curvature dial.

---

## 5. Avenue E: known tables at another level (data combination)

### 5.1 Plain language

* **What it is.** Sometimes the cross-tabulation we want is published for a *larger* area. NZ publishes each electorate's party-by-candidate table, and US states publish some exit-poll or registration tables. We then want rates for *sub*-areas (urban vs rural booths) or for areas without a table.
* **Using it.** If the same contextual model holds across the larger area, its known table pins part of the unknown "who is behind the pattern", and sub-area bounds tighten a lot.
* **Key assumption:** the model, with shared $K$, holds within the area whose table is known.
* **If wrong:** the sub-area bounds can exclude the truth.
* **Checks:**
  * the known table must lie inside the model's bound for the larger area (testable, a free specification test);
  * emptiness.

### 5.2 Content

**E1 (proved; verified in `check_structural_claims.py` for 7 table shapes).** $K\mapsto B^*(K)$ has rank $(R-1)(C-1)$.

*Reason:* in each column, the $\mu_r$-weighted sum of the $K$-terms cancels, giving $C-1$ independent columns times $R-1$ independent rows.

A known aggregate table therefore pins:
* **all of $K$ when $R=2$**, so any sub-area target is point-identified under M0 within the area;
* only $(R-1)(C-1)$ of $(C-1)R(R-1)/2$ directions when $R\ge3$. For 3×3, 4 of 6 directions are pinned and 2 remain free.

These are **equality rows** in the LP.

**E2. Transfer.** $\lVert K_{\text{target}}-K_{\text{known}}\rVert\le\delta$ across areas, a dial.

* This overlaps with cross-period stability, which is the DiD session's topic. Only the *cross-area* version is proposed here.
* **Exploratory evidence** (`nz_diagnostics.py`, D2): the part of each NZ table that M0 must attribute to $K$ is large (median absolute value 0.11–0.15 for the main cells). Its correlation across consecutive elections in the same electorate is only 0.03–0.51.
* So a transfer dial would need a large $\delta$. Strong stability of $K$ is not supported in NZ. *This is shared with the DiD session as a data point.*

**E3. Nested geography with margins only adds nothing (proved, trivial).** Coarse margins are sums of fine margins, so they carry no new information. Gains need either cross-tabs at some level, or categories cut differently across levels. The latter is a multi-way Fréchet / data-fusion problem; see Ridder & Moffitt (2007, *Handbook of Econometrics* 6B) and Fan, Sherman & Shum (2014, *Econometrica*).

**Verification for the benchmark:** NZ has the electorate truth but no sub-electorate truth. Phase 2 can validate E1 by splitting each electorate into pseudo-sub-areas in simulation, and on NZ by "leave-one-electorate-out transfer".

---

## 6. Avenue F: better use of $P_{T\mid X}$ (sharp set; already pinned)

The series-based $\Theta_X$ is agreed for the second release. This proposal adds:

* **A0 (§1.3)** as a closed-form, low-dimensional outer approximation of $\Theta_X$. In our designs it was not binding.
* **A finite core-determining class (conjecture).** For fixed $x$, the transportation polytopes $S(x,t)$ take finitely many combinatorial types as $t$ varies. Hence (T1) needs only the finitely many directions $u$ in the common refinement of their normal fans, and $\Theta_X$ with $x$ discretised is an exact finite LP.
  * Related: core-determining classes (Galichon & Henry 2011, *REStud*); Beresteanu, Molchanov & Molinari (2011).
* **Diagnostic:** a "kink mass" statistic, the share of units whose $T$ lies on both sides of a DD kink within an $x$-neighbourhood. It predicts when $\Theta_X$ will beat $\Theta_O$ (P2a Prop. 3). In NZ, removing small booths from the composition range (`M0 trim`) lowers the share of empty plug-in sets from 46% to 21%. *Heuristic reading:* small units carry much of the extreme-composition and kink information, and also much of the noise.

---

## 7. Lower-ranked avenues

**G. Label symmetry (equivariance).** When rows and columns share labels (NZ: party × candidate), assume the model is invariant under swapping Labour and National. Then $K$ lies in the invariant subspace, roughly halving $d$ (**proved** as linear algebra).
* The $Q$-side analogue is testable.
* Plausibility is low in asymmetric party systems; it is more natural for yes/no referenda framed symmetrically.

**H. Deconvolution (proved negative).** If $\beta_r=m_r(x)+e_r$ with $e$ independent of $x$, then the law of $T\mid x$ is the law of $\tau(x)+\sum_rx_re_r$, which does not involve $K$.
* Kotlarski-type arguments (row-independent $e$) identify the shock laws.
* $K$ enters only through the support condition $m_r(x)+\operatorname{supp}e_r\subseteq[0,1]$. That tightens Θ_O's rows from $m\ge0$ to $m\ge-\min e_r$. It is small, and row-independence is implausible.

**I. Spatial structure (pinned item 10).** The new use is as *replicates* for Avenue A: neighbouring polls with similar $x$ give local estimates of $\operatorname{Var}(T\mid x,n)$ without a global polynomial.

**J. Entropy, optimal transport and Choquet (reporting, no identifying power).**
* The maximum-entropy point of the identified set is a *selection*, not an estimate. It is the information-projection idea of Judge, Miller & Cho (2004, in King, Rosen & Tanner (eds.), *Ecological Inference: New Methodological Strategies*; *verify chapter*). In OT language it is the entropic (Sinkhorn) projection onto the Aumann set.
* It may be offered as a labelled default point with its distance to each endpoint.
* Capacity (Choquet) characterisations (Artstein 1983) are equivalent to Theorem 1, so they are useful only computationally (F).

---

## 8. What would fail, and why (cross-cutting)

* **A's inward bias** is the Phase-3 problem squared: there are many nearly binding, noisy constraints that depend on second moments. The integrated (averaged) constraint is the only version shown valid so far. Phase 3 for A needs moment-inequality inference for conditional variance functions:
  * Andrews & Shi (2013, *Econometrica*), conditional moment inequalities;
  * CLR (2013), already used.
* **A is invalid under within-poll dependence.** A design-effect allowance $\text{deff}\in[1,\bar d]$ replaces $B$ by $B/\bar d$ in (A1). That is a sensitivity dial, and NZ's slope below 1 suggests $\bar d$ near 1 there.
* **B, C and E are all assumptions about the unidentified $K$.** They are calibrated, or partially tested, by identified quantities, but never fully testable (Prop. 9 still applies).

## 9. Ranking (criteria from the brief)

Scores run 1–5, where 5 is best.

| Avenue | Power | Plausibility (no new assumption = 5) | Computability / LP fit | Testability | Overall |
|---|---|---|---|---|---|
| **A** unit size (A1; A2) | 5 | 3 (A1), 2 (A2) | 4 (2×2 LP with ordering; SOCP/LMI or cuts in R×C) | 4 (overdispersion slope, emptiness) | **1st for research, 3rd for near-term use** |
| **B** two-scale context | 4 | 3 | 5 (LP rows, few constraints) | 3 (B3 and B2 testable; B1 calibrated) | **1st for near-term use** |
| **E** known tables / transfer | 4 (sub-areas) | 4 (E1 within area) / 2 (E2) | 5 (equality rows) | 4 | **2nd** |
| **C** several contests | 3–5 (data-dependent) | 3 (C2) / 2 (C3) / 1 (C4) | 5 (C2) / 3 | 4 | 4th |
| **D** Bernstein ladder | — (robustness) | 5 | 5 | 3 | 5th (useful infrastructure) |
| **F** sharp-set computation | 2 | 5 | 3 | — | 6th (already pinned) |
| G / H / I / J | 1 | 2–5 | 3–5 | 2 | low |

**Suggested Phase 2 order:**
1. **B**: Prop. B1 with the R×C rank condition, the B1(M) LP, and B2's overidentification test.
2. **A**: the A1/A2 rungs in R×C, an LMI or cut formulation, and the identification of $A(x)$, $B(x)$ with a minimal functional form.
3. **E1**: equality rows and sub-area targets.
4. **C2.**

D can be done at any point, as infrastructure.

## 10. Clarifying questions

1. **Binomial voting as a rung (Avenue A).**
   * Are you willing to put "independent voting within poll × group cells" on the ladder, given that every likelihood-based EI model assumes it?
   * Would you rather present A only with the design-effect dial?
2. **Data for Avenue B.** For BC stage 2, is a sub-riding composition $z$ available? Options are census dissemination areas or the composition of adjacent polls, which would need poll geography. For NZ, meshblock-level data?
3. **Pooling for Avenue A.** A needs thousands of units. Is a "variance structure shared across electorates (or ridings)" assumption acceptable? It is much weaker than sharing $K$: it shares only the over-dispersion and heterogeneity ratio.
4. **Several contests (C).** Do you have, or want, a multi-contest dataset? Candidates are US CVRs, Swiss referendum days, NZ 2020 referenda by booth, or BC's 2018 electoral-reform referendum by poll, if published. Without one, C stays theoretical.
5. **Solvers.** In R×C, A1 is a second-order-cone or semidefinite constraint. Is a dependency such as `CVXR`, `ROI` or `scs` acceptable, or should everything stay in `lpSolve` via tangent-cut outer approximations (valid but slightly wider)?
6. **Scope boundary with the DiD session.** E2 (transfer of $K$) overlaps with temporal stability. Should I hand the NZ persistence evidence (D2) to that session and keep only cross-area transfer here?

## 11. Files

| File | What it checks | Result |
|---|---|---|
| `new_avenues/sim_second_moments.py` → `_output.md` | Avenue A, identification and estimation, four designs | large identification gain; valid estimation needs many units |
| `new_avenues/sim_two_scale_context.py` → `_output.md` | Avenue B, rungs B0–B3, true and false $M$ | smooth gain in $M$; correct failure when $M<\lambda$ |
| `new_avenues/check_structural_claims.py` → `_output.md` | C1 (CVRs add nothing), E1 (rank of $K\mapsto B^*$) | PASS / PASS |
| `new_avenues/export_nz_units.R`, `new_avenues/nz_diagnostics.py` → `_output.md` | D1: binomial slope in NZ; D2: persistence of the K-part | slope 0.69–0.81 (SE ≈ 0.10); persistence 0.03–0.51 |

To reproduce the NZ diagnostics:

```sh
Rscript docs/extension/new_avenues/export_nz_units.R <dir>   # needs ei.Datasets
python3 docs/extension/new_avenues/nz_diagnostics.py <dir>
```
