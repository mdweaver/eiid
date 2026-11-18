# Response to the external review of Phases 1, 2a, 2b and 3

*Two reports: the first on Phases 1 and 2a; the revised report on Phases 1, 2a, 2b and 3. The revised documents are on branch `claude/wizardly-volta-lbh4bm` of https://github.com/mdweaver/eiid. The R package is `eiPartialID` 0.3.0. Our working ledger is `docs/extension/REVIEW_RESPONSE_PLAN.md`.*

## Cover note

We thank the referee for two careful reports. We accepted every point. Each mathematical correction comes with a numerical check, in the `verify_*.py` scripts, that **fails on the old statement and passes on the new one**. The revision went beyond the corrections in three ways, described below:
* the inference framework now states where its randomness comes from;
* the regime coverage of each inference result is tabulated;
* a later review of the M-SIGN module found three further corrections to the core documents, which we report here.

**Errors corrected.**
* **S3 and the hull (R-04).** Probability validity follows automatically only when $A$ is covered by the hull of the support inside $A$ (2a Lemma 3(c), with the referee's singleton counterexample as Remark 5.1). The later review sharpened this to the support of the measure restricted to $A$ (see below).
* **Compactness (R-07).** The condition is $\operatorname{affdim}(A)=R-1$. A5 implies it, so Theorem 1 itself was unaffected (Remark 5.2).
* **Value-function regularity (R-06).** It is concave in $Q$ for fixed $\mu$, convex in the objective for fixed $Q$, and in general neither jointly (the bilinear example).
* **D4(b) (R-23).** It now conditions on $(X,N)$; $D\perp X$ alone is not enough (Remark D4.1).
* **Smaller corrections.** The M1 restriction count (R-30), turnover marked open (R-35), Theorem P3.4(ii)'s condition on distinct projected gradients (R-27), and "only $\tilde B$ is affine in $K$" (R-33).

**Scope and assumptions.**
* Phase 1 has a new §1.7 separating the three identified sets ($\Theta_I\subseteq\Theta_X\subseteq\Theta_O$, with the information each uses), the three hulls (population, empirical, sample hull as an estimate) and the three targets ($B^*$, $B$, $\tilde B$).
* Every result that uses the affine form of the estimand now states $P(X\in A)=1$. A new Proposition 7′ bounds units outside $A$ by Duncan–Davis.
* The M1 rank condition is now a stacked rank condition (2b Lemma D1(d)).
* Theorem 1's functional-analytic steps are written out (F1–F7).
* The citation table records which theorem numbers were checked in the source. The CLR result we use is their Theorem 4, not Theorem 2.

**Clustering (R-28).** We accepted that the earlier Theorem P3.5(a) assumed independence that the sampling assumption did not grant. The fix goes further than adding a cluster allowance (P3 §0.4).
* Following Abadie, Athey, Imbens and Wooldridge (2020, 2023), clustering is justified by clustered sampling or clustered assignment. Ecological inference has neither: every precinct and riding is observed, and nothing is assigned.
* The uncertainty is therefore model-based, conditional on the design. The contextual-effects assumption concerns $E[\beta\mid V]$, which needs the precinct rates to be random given the design.
* The targets $B$ and $\tilde B$ refer to this electorate, so a riding-level shock is a fixed feature of it:
  * under M1/M2 it is a district effect, part of the mean;
  * under M0 it is a misspecification, which a new pooling test detects (P3.6(c), conditional on the design).
* Units are assumed independent given the design and the district effects (P3-S1a). Dependence within a riding beyond its riding effect, such as neighbourhood spillovers, is a stated limitation.
* A superpopulation target that averages over re-drawn riding shocks is not pursued. The package's `cluster` option was removed.

**Inference: what is proved, and where.**
* **The candidate collection (R-29).**
  * (b1) CLR with the full dual-vertex set is a cited theorem.
  * (b2) is proved directly from LP weak duality: any candidate set, with a critical value taken over a fixed superset or by projection.
  * (b3) The package's simulated candidate collection is a conjecture, with a proved reduction to (b1).
  * The package reports the CLR interval, labelled algorithmic, and the projection interval beside it.
* **Regimes (R-31).** P3 §0.5 tabulates the regime coverage of each result. M0 is covered in every regime; M1/M2 with a fixed number of districts. With many small districts, coverage of the CLR interval rests on simulation. A new finite-sample interval (Theorem P3.8) is valid in every regime but, in every design we tried, as wide as the logical range. We say so.
* **Sample hull (R-32).** Coverage of the true parameter does not depend on the sample hull, because the true parameter satisfies the constraints at the sample-hull vertices. The hull affects only width and the endpoint's distribution, so P3.4 is stated for the sample-hull endpoint.

**Monte Carlo (R-37).**
* 1,000 replications per design, with Monte Carlo standard errors.
* Designs: riding shocks, a design where Slater fails, and up to 200 districts with 10 polls each.
* Results:
  * No proved interval under-covered beyond Monte Carlo error.
  * The CLR interval covered at least 0.994 except where Slater fails (0.944, at the nominal level).
  * With many small districts only the CLR interval is informative.
  * The pooling test over-rejected with many restrictions. Both tests are now calibrated by a restricted wild bootstrap, with sizes 0.043 and 0.054 at nominal 5%.
* The NZ benchmark found that M0 and M1 by electorate are both rejected, and the M1 CLR interval missed in 5 of 63 non-empty cells. This is reported as a known limit.

**Further core corrections from the later review (Batch B).** The review of the M-SIGN module included three issues in these documents. We report them here because they correct results you reviewed.
1. **Support of the restricted measure (2a).** Almost-sure equality on $A$ extends only to the support of $P_X$ restricted to $A$, which can be smaller than $\operatorname{supp}P_X\cap A$. The notation $\operatorname{supp}_AP_X$ is now used throughout (2a §0.1, Lemma 1, Lemma 3(c), Remarks 5.1–5.2), and the tension between Lemma 3(c) and Remark 5.1 is resolved.
2. **Projection under a singular covariance (P3.5(b2)(ii)).** The rank-$q$ argument does not control errors in null directions. The proposition now requires a positive definite $\Sigma$, or a reduced full-rank parameterisation with exactly zero null-direction error. The package drops candidates whose gradients use a numerical null direction. In the Monte Carlo the covariance had full rank in every replication.
3. **Residual CLT (P3.5(a), step 5).** A small maximum weight does not by itself give a variance-standardised CLT. The step is completed by a Chebyshev/Lindeberg subsequence argument, and the allowance is unchanged.

**Not done, by decision.** Covariate adjustment (R-15) is a planned extension, at low priority, because it needs stronger assumptions. Bounded turnover remains open (R-35).

## Point-by-point
| Point | Reviewer § | Type | Location | Change | Status |
|---|---|---|---|---|---|
| R-01 | 2.1 | X | P1 §1.3 | Replace the "open subset of the simplex" argument for $\sum_cH_c=J$ with the A5 argument from 2a Lemma 1. | Done |
| R-02 | 2.1 | S | P1 A5; 2a Lemma 1 | State A5 in matrix form ($x^\top Dx=0$ on $\operatorname{supp}P_X\cap A$ ⇒ $D=0$ for symmetric $D$). Add a remark that for $R=3$, three distinct points are not enough: the design must identify all $R(R+1)/2$ coefficients. | Done |
| R-03 | 2.2 | X | 2a Lemma 2 | Replace $\partial B^*/\partial K_{c,rs}$ with "the coefficient on $k_{c,rs}$ in the free parametrisation of Prop. 7 ($c<C$, $r<s$)". | Done |
| R-04 | 2.4 | **E** | 2a §0.2 (S3), §10 item 2, Lemma 3(c) | Restate: S3 follows from S1–S2 when $A=\operatorname{conv}\operatorname{supp}P_X$, or more generally when $A\subseteq\operatorname{conv}(\operatorname{supp}P_X\cap A)$. Otherwise S3 is an extra assumption. Add the reviewer's singleton counterexample as a remark. | Done |
| R-05 | 2.5 | S | 2a Prop. 4–7, Cor. 1 | Add $P(X\in A)=1$ to the hypotheses of every result that uses Lemma 2 (Prop. 5 satisfies it automatically). Add Prop. 7′ (D2). Keep Remark 2.1 as the bridge to 2b. | Done |
| R-06 | 2.6 | **E** | 2a Prop. 7, closing paragraph | Replace with: (i) for fixed $\mu$, the value is concave and piecewise linear in $Q$; (ii) for fixed $Q$, it is convex and piecewise linear in $\gamma$; (iii) jointly it is neither in general. Add the bilinear example $V=\mu Q$. Point to P3 Lemma P3.1 for the fixed-design case. | Done |
| R-07 | 2.7 | **E** | 2a Theorem 1 (compactness), Lemma 3(b), "polytope" wording | Replace "nonempty relative interior" with $\operatorname{affdim}(A)=R-1$. Add the singleton counterexample with $Kx_0=0$. Call $\Theta_O$ a polyhedron until boundedness is shown. Add a remark that A5 plus $A=\operatorname{conv}\operatorname{supp}P_X$ implies full dimension. | Done |
| R-08 | 2.3 | P | 2a Theorem 1 proof | Expand the sufficiency proof in six numbered facts: $\mathcal S\subset L^1(P;\mathbb R^{R\times C})$; $L^\infty$-bounded ⇒ uniformly integrable; convex and norm-closed ⇒ weakly closed; Dunford–Pettis ⇒ weakly compact; $E[\cdot\mid V]$ is weak-to-weak continuous; the support function is attained by a measurable maximiser. Add a paragraph on why the inequality must hold for all $u$ outside one common null set. | Done |
| R-09 | 2.3, 6 | P | 2a §0.3 († citations) | Check every † citation against a named edition: Hiai–Umegaki, Aumann/Kuratowski–Ryll-Nardzewski, measurable maximum (Aliprantis–Border 18.19), Dunford–Pettis, Fang–Santos, CLR Thm 2. Record edition and page in a new table, §0.4. | Done for the papers we could read (CLR 2013, Fang–Santos 2019; one wrong theorem number corrected). Textbook theorem numbers are marked 'standard' and still to be checked against the printed editions |
| R-10 | 3.2 | S | 2a Prop. 2, 3, 4, 7 | Wherever $\Theta_O$ is called sharp, add "for the regression function". Prop. 4 (JKST equivalence) and Prop. 7 bound the **outer** set. Keep "sufficient" for Prop. 3. | Done |
| R-11 | 3.3 | S | 2a §0.2; P1 §1.2 | Add a short subsection "Three hulls": population $A=\operatorname{conv}\operatorname{supp}P_X$; the empirical-law hull (finite-sample identification); the sample hull $\hat A$ as a plug-in. Point forward to P3 Prop. P3.2. | Done |
| R-12 | 3.4 | S | 2a §4 (before Theorem 1) | Add a remark that sharpness is relative to the continuous latent-rate model, not to the finite set of integer cell counts. | Done |
| R-13 | 4.2 | X | 2a Prop. 5 | Add: the interval is empty unless $0\le Q_{c,rs}$ and $\sum_cQ_{c,rs}=1$ hold, as the probability restrictions require. | Done |
| R-14 | 4.4 | X | 2a Cor. 6.1(c) | Replace "strong contextual effects" with the exact condition: part (b)'s common-sign condition fails, which (given the sum restrictions) means at least one pure-poll extrapolated rate is negative. | Done |
| R-15 | 3.1, rev. 5.1 | S | P1 §1.5, §1.6, A3 | Apply D1. Rewrite A3 as $E[\beta\mid X,N]=M_H(X)$ (with district shifts under M1/M2) and remove $Z$ from $\mathcal G$. Move M1′ and continuous covariates to a "Planned extension" paragraph that names the extra assumptions. No text may say the covariate case is proved. | Done for this revision (no covariates in the mean; A3 rewritten). Covariate adjustment is a planned extension (Phase 2c), not claimed |
| R-16 | 5 | X | 2a Prop. 9 | State the domain: the restrictions on $\Delta$ hold on all of $\Delta^{R-1}$ for the "linear $\Delta=M_{K'}$" characterisation. Move a short version of Prop. 9 into the introduction or §2 as motivation. | Done |
| R-17 | 6 | X | 2a new §4.3 | Add "What Theorem 1 and Prop. 7 each need": the sharp set needs the law of $T$ given $(X,N)$, which in practice means smoothing or pooling because every $(X,N)$ is unique. The outer set needs only the regression. Point to the series-based $\Theta_X$ planned for the second release. | Done |
| R-18 | 3.2 | X | 2b Prop. D0 | Separate the sharp interval for each cell from the joint set for $B^*$, which is not a Cartesian product. | Done |
| R-19 | 3.3 | S | 2b Lemma D1, Remark D1.1; P1 A9 | Replace "$R$ distinct compositions per district" with: the stacked design (pooled quadratic monomials plus district × $X_r$ shifts) has full column rank. Under M1 this requires each district's support to be affinely spanning ($R$ affinely independent points). Raise Remark D1.1 (M1 needs only linear variation within a district) into the main text. | Done |
| R-20 | 3.5 | P | 2b Theorem 1′ | List the hypotheses explicitly: $v\mapsto m_\theta(v)$ Borel; uniformly bounded; $\theta\mapsto m_\theta$ affine and continuous; selection set nonempty and compact. | Done |
| R-21 | 3.6 | X | 2b Theorem 2(c), ladder figure | Label M1(0) as a *sensitivity benchmark*: it keeps the common skew part but drops the testable symmetric restrictions. | Done |
| R-22 | 3.7 | S | 2b Prop. D3 | Make $0<l_d\le u_d<1$ a visible condition. Derive the limiting constraint for districts with pure-composition units ($t\in\{0,1\}$) as a separate case. | Done |
| R-23 | 3.8 corr. | **E** | 2b Prop. D4(b); P1 §1.6 (D4) | Condition on $(X,N)$: M0 holds iff $E[Z_D\mid X,N]$ is affine in $X$ and does not depend on $N$. The sufficient condition becomes $D\perp(X,N)$. Add a counterexample to "$D\perp X$ suffices" (district depends on $N$). Qualify "exactly when": the conclusion concerns the mean model, and the outer bounds also need S3 on the pooled region. | Done |
| R-24 | 3.9 | X | 2b Remark 6.2 → Proposition | Promote the DiD result to a proposition. State in words that the $NX_rX_s$ and $NX_r$ moments must be equal across periods, not just similar share distributions. | Done: now 2b Proposition D8, with the exact moment condition |
| R-25 | 4.2 | S | P3 Lemma P3.1 | Add: this holds for a fixed constraint system and objective. It says nothing about joint variation in $(Q,\mu)$ (see R-06). | Done |
| R-26 | 4.3 | S/P | P3 Prop. P3.2 | (i) State convergence along a nested sampling sequence, or add an argument for non-nested hulls. (ii) Call Slater a strict-interiority condition and say it fails when an implied rate is 0 at a vertex. Note that P3.5's coverage does not need Slater (step 2 uses only feasibility). (iii) Split the error into estimation error ($O_p(n^{-1/2})$) and hull error (≥ 0, rate not shown). | Done |
| R-27 | 4.4 | **E** (minor) | P3 Theorem P3.4(ii) | The non-Gaussian limit requires at least two active gradients with $(g_\lambda-g_{\lambda'})^\top\Sigma(g_\lambda-g_{\lambda'})>0$. Parametrise only $C-1$ outcome columns, so $\Sigma$ is nonsingular. | Done |
| R-28 | 4.5 | **E** (assumption mismatch) | P3 P3-S1, Theorem P3.5(a); P1 A6; P2b D4 | Apply D3. **Docs:** (1) Replace P3-S1 with *S1a: given the design $\{(N_i,X_i,D_i)\}$ and the district effects, the unit residuals $e_i$ are independent.* Delete "arbitrary dependence within districts for M0, cluster-robust variances." (2) New P3 §0.4, "Where the randomness comes from": every precinct and riding is observed, so the uncertainty is neither sampling-based nor design-based in the Abadie et al. sense. It is model-based: the contextual-effects assumption is a statement about $E[\beta\mid V]$, which needs $\beta_i$ to be random given $V_i$. The targets $B$ and $\tilde B$ condition on this electorate, so riding shocks are fixed and belong in the mean (M1/M2), not in the residual. Under M0, a riding shock that is non-affine in composition is the omitted-variable failure of Prop. D4, which clustered SEs would not repair. (3) State the remaining assumption plainly: no dependence within a riding beyond its riding effect (e.g. neighbourhood spillovers). Name it as a limitation. **Response:** accept the reviewer's point that P3.5(a) needed independence, and say that the framework now assumes it explicitly, with the reasons above. Use the reviewer's common-shock example: under our target it is a riding effect, absorbed by M1 and detected by the new pooling test. | Done, and changed in substance: no clustered variance (see 'Clustering' below) |
| R-29 | 4.6 | **E** (claim) | P3 Theorem P3.5(b) | Split it. **(b1) Theorem:** CLR with the full vertex set $\Lambda$, which is finite for a fixed design. **(b2) Algorithm:** collect candidates at simulated draws. **Proof route for a valid (b2):** by weak duality, $\theta_\lambda(\beta_0)\ge\varphi_U(\beta_0)$ for *every* dual-feasible $\lambda$. So $\min_{\lambda\in\hat V}(\hat\theta_\lambda+\hat k\hat\sigma_\lambda)$ covers $\varphi_U$ whenever $\hat k$ controls the maximum over a *fixed* set containing $\hat V$. A data-dependent subset is valid if its critical value is computed over the fixed superset (or the collected set is independent of $\hat\beta$ given the draws, to be checked). Write this as Prop. P3.5(b2) if it closes, otherwise as a Conjecture. List rarely selected vertices, degenerate faces, singular $\Sigma$ and a growing $\lvert\Lambda\rvert$ as the conditions to address. | Done: (b1) theorem, (b2) proved, (b3) conjecture; the package reports CLR (algorithmic) with the projection interval beside it |
| R-30 | 4.7 | **E** (P1) | P1 §1.6; P3 Prop. P3.6(b) | Fix P1's count to $(D^\star-1)(C-1)R(R-1)/2$, with the derivation $(C-1)R(R+1)/2-(C-1)R$. State the Wald/F test for fixed dimension. Under a growing number of restrictions it becomes a conjecture with simulation support. | Done |
| R-31 | 4.8–4.9 | **E** (scope) | P3 §0.1, §6, Prop. P3.7 | Apply D5. Replace "results are for R1/R3" with a table of which result covers each model and regime (D5 (i)–(iv)). Add **Theorem P3.8 (finite-sample CI)**. P3.7(a), "dropping constraints enlarges the set", stays a proposition. P3.7(b) becomes a conditional statement: *if* the procedure is valid on the kept districts, choosing them from the design alone keeps conditional coverage. "CLR restores coverage" in R2 becomes a *Simulation finding*. | Done; the finite-sample interval (Theorem P3.8) is valid but uninformative in the designs tried |
| R-32 | 4.3, 9 | S | P3 §3–4 | Central point to add: **coverage of the true parameter does not depend on the hull.** The true $k_0$ satisfies the population constraints, which include those at sample-hull vertices ($\hat A_d\subseteq A_d$). So P3.5 step 2 holds with $\hat A$, and the hull affects only width and the distribution of the endpoint estimator. State P3.4 as a fixed-design result for the sample-hull endpoint $\varphi(\beta_0;\hat A)$. | Done |
| R-33 | 4.10 | S/X | P3 §0.2 | Three targets: $B^*$ (superpopulation), $B$ (realised), $\tilde B$ (conditional mean given the design, $=E[B\mid\text{design}]$). Inference is conditional on the design. Covering $B^*$ would also need design-sampling variation, which is left open. Fix "both linear in $K$": only $\tilde B$ is affine in $K$. | Done |
| R-34 | 4.11 | S | P3 Prop. S1 | Label it an outer bound under bounded deviations, not the sharp set. | Done |
| R-35 | rev. 5.2 | **E** (inconsistent) | P1 A1; P3 §7 | Change P1 A1 to say turnover is an open extension, as P3 §7 says. Note that it changes both the DD sets and the regression identity. | Done |
| R-36 | rev. 5 | S | P1 A6, A10 | Align A6 with S1a (D3) and A10 with the regime table (R-31). | Done |
| R-37 | rev. 6 | X | P3 §8 | Report Monte Carlo standard errors, and use at least 1000 replications for coverage cells. Add designs: (i) riding shocks correlated with composition, showing that M0 fails, M1 covers, and the pooling test detects the shocks; (ii) truth on the boundary, so Slater fails; (iii) R2 with growing $D^\star$, comparing the CLR (simulation finding) and finite-sample (theorem) intervals. | Done (1,000 replications per design, MC standard errors, new designs; P3 §8) |

## What we now regard as established (update of the referee's table)

| Result | Referee's assessment | Status now |
|---|---|---|
| R×C identification decomposition; dimension of $\mathcal K$ | Established | Unchanged; A5 in matrix form, on $\operatorname{supp}_AP_X$ |
| Conditional-selection characterisation of $\Theta_I$ | Established, subject to conditions | Unchanged; proof expanded (F1–F7) |
| Compactness | Required a precise dimensionality condition | $\operatorname{affdim}(A)=R-1$; implied by A5 |
| Joint LP value-function concavity | Incorrect as stated | Corrected (three separate statements) |
| D4(b) | Incorrect conditioning | Corrected (conditions on $(X,N)$) |
| Sample-hull one-sided conservatism | Established | Unchanged; coverage of the parameter shown not to depend on the hull |
| Root-n plug-in under Slater | Established for fixed dimension | Unchanged; Slater stated as strict interiority; not needed for coverage |
| Directional delta method | Established for fixed design and hull | Stated for the sample-hull endpoint; distinct projected gradients |
| Conservative CI for $\tilde B$; for $B$ with independent residuals | Established | Unchanged; residual step completed (Chebyshev/Lindeberg) |
| CI for $B$ with within-district dependence | Not established | Withdrawn: no clustering (P3 §0.4) |
| CLR with the full candidate set | Appropriate, subject to conditions | Cited theorem (CLR Thm 4, Condition P) |
| Simulated CLR candidate generation | Outstanding | Conjecture P3.5(b3), with a proved reduction; the package labels it algorithmic |
| Projection interval | — (new) | Proved (P3.5(b2)), with the nonsingularity condition |
| Inference in R2/R3 | Requires additional theory | Regime table (P3 §0.5); finite-sample interval valid but uninformative; CLR in R2 a simulation finding |
| Bounded nonlinear sensitivity | Established as an outer bound | Labelled as an outer bound |
| Bounded-turnover sensitivity | Open | Open |
