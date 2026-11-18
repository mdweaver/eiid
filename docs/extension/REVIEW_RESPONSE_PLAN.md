# Plan: responding to the two reviewer reports

*Covers the first report (Phases 1 and 2a) and the revised report (Phases 1, 2a, 2b and 3). Each comment gets an ID (R-xx), a type, the files it touches, the action, any package impact, and how the fix will be checked. Work is grouped into workstreams. They run in the order given because later ones cite the corrected statements from earlier ones.*

**Types:**
* **E**: a mathematical claim that is wrong as stated
* **S**: a correct result whose scope or assumptions must be stated
* **P**: a proof that is correct but needs more detail or checked citations
* **X**: exposition or emphasis only
* **C**: the package (`eiPartialID`) must change, not only the docs

**Status labels used after the revision (P1 §5 rule 1, made stricter):**
* *Theorem/Proposition*: proved, with its full hypothesis list
* *Proposition (fixed design)*: proved conditional on the design and a fixed hull
* *Conjecture*: a proof route exists but is not finished
* *Algorithm*: a computational procedure with no coverage claim of its own
* *Simulation finding*: Monte Carlo evidence only

---

## 0. Decisions needed before starting

These change what gets written, so they come first. The recommended choice is listed first.

| # | Decision | Recommended | Alternative |
|---|---|---|---|
| D1 | Covariates $Z$ (R-15) | **Decided (Nov 2026): Option B.** The package gets covariate adjustment. Develop a new Phase 2c with two rungs: (i) $Z$ shifts baselines only, $E[\beta\mid X,N,Z]=M_H(X)+\Gamma(Z)$ with $\Gamma(Z)$ linear in $Z$ (this generalises M1 and JKST App. C, and adds no unidentified dimensions); (ii) $Z$ also modifies contextual effects, $H_c(Z)=H_{c0}+\sum_\ell Z_\ell G_{c\ell}$ (M1′; the skew parts of the $G_{c\ell}$ add $L(C-1)R(R-1)/2$ unidentified dimensions). $\log N$ as a covariate relaxes A3. | — |
| D2 | Units outside $A$ (R-05) | **Carry $P(X\in A)=1$ through §6 of 2a.** Also add a short Prop. 7′ that writes down what `range = "trim"` already does: the LP over kept units plus Duncan–Davis for dropped units. The code already does this, so the proof is cheap. | Assumption only, with no Prop. 7′. |
| D3 | Unit-level term under clustering (R-28) | **Both.** State P3.5(a) for independent units (Option A). Add P3.5(c), a cluster allowance (Option B). The package already accepts `cluster`, so Option A alone would leave a published option without a valid interval. | Option A only, and stop `ci_actual` when `cluster` is given. |
| D4 | Default CI in the package (R-29) | **Decided (Nov 2026).** **Keep CLR-type as the default, but compute the critical value over a fixed superset of candidates (see WS3, R-29).** This makes it provably valid. Report the simulated-collection version as an *Algorithm* option. | Make the conservative P3.5(a) interval the default, and keep CLR as an option. |
| D5 | Formal asymptotic scope (R-31) | **R1 only (fixed $D^\star$, fixed dimension).** State R2/R3 results as conjectures backed by simulation. | Also prove R3 for M0 with growing independent clusters. This is feasible: standard cluster asymptotics. |

---

## 1. Ledger: every comment, mapped to an action

### First report (Phases 1, 2a)

| ID | Reviewer § | Type | Location | Action | Package |
|---|---|---|---|---|---|
| R-01 | 2.1 | X | P1 §1.3 | Replace the "open subset of the simplex" argument for $\sum_cH_c=J$ with the A5 argument from 2a Lemma 1. | — |
| R-02 | 2.1 | S | P1 A5; 2a Lemma 1 | State A5 in matrix form ($x^\top Dx=0$ on $\operatorname{supp}P_X\cap A$ ⇒ $D=0$ for symmetric $D$). Add a remark that for $R=3$, three distinct points are not enough: the design must identify all $R(R+1)/2$ coefficients. | Check the rank of the design matrix with a clear error (it may already exist; verify) |
| R-03 | 2.2 | X | 2a Lemma 2 | Replace $\partial B^*/\partial K_{c,rs}$ with "the coefficient on $k_{c,rs}$ in the free parametrisation of Prop. 7 ($c<C$, $r<s$)". | — |
| R-04 | 2.4 | **E** | 2a §0.2 (S3), §10 item 2, Lemma 3(c) | Restate: S3 follows from S1–S2 when $A=\operatorname{conv}\operatorname{supp}P_X$, or more generally when $A\subseteq\operatorname{conv}(\operatorname{supp}P_X\cap A)$. Otherwise S3 is an extra assumption. Add the reviewer's singleton counterexample as a remark. | — |
| R-05 | 2.5 | S | 2a Prop. 4–7, Cor. 1 | Add $P(X\in A)=1$ to the hypotheses of every result that uses Lemma 2 (Prop. 5 satisfies it automatically). Add Prop. 7′ (D2). Keep Remark 2.1 as the bridge to 2b. | Document that `trim` is Prop. 7′ |
| R-06 | 2.6 | **E** | 2a Prop. 7, closing paragraph | Replace with: (i) for fixed $\mu$, the value is concave and piecewise linear in $Q$; (ii) for fixed $Q$, it is convex and piecewise linear in $\gamma$; (iii) jointly it is neither in general. Add the bilinear example $V=\mu Q$. Point to P3 Lemma P3.1 for the fixed-design case. | — |
| R-07 | 2.7 | **E** | 2a Theorem 1 (compactness), Lemma 3(b), "polytope" wording | Replace "nonempty relative interior" with $\operatorname{affdim}(A)=R-1$. Add the singleton counterexample with $Kx_0=0$. Call $\Theta_O$ a polyhedron until boundedness is shown. Add a remark that A5 plus $A=\operatorname{conv}\operatorname{supp}P_X$ implies full dimension. | Warn when a group's hull is lower-dimensional (unbounded LP) |
| R-08 | 2.3 | P | 2a Theorem 1 proof | Expand the sufficiency proof in six numbered facts: $\mathcal S\subset L^1(P;\mathbb R^{R\times C})$; $L^\infty$-bounded ⇒ uniformly integrable; convex and norm-closed ⇒ weakly closed; Dunford–Pettis ⇒ weakly compact; $E[\cdot\mid V]$ is weak-to-weak continuous; the support function is attained by a measurable maximiser. Add a paragraph on why the inequality must hold for all $u$ outside one common null set. | — |
| R-09 | 2.3, 6 | P | 2a §0.3 († citations) | Check every † citation against a named edition: Hiai–Umegaki, Aumann/Kuratowski–Ryll-Nardzewski, measurable maximum (Aliprantis–Border 18.19), Dunford–Pettis, Fang–Santos, CLR Thm 2. Record edition and page in a new table, §0.4. | — |
| R-10 | 3.2 | S | 2a Prop. 2, 3, 4, 7 | Wherever $\Theta_O$ is called sharp, add "for the regression function". Prop. 4 (JKST equivalence) and Prop. 7 bound the **outer** set. Keep "sufficient" for Prop. 3. | Print and vignette wording: "outer bounds" |
| R-11 | 3.3 | S | 2a §0.2; P1 §1.2 | Add a short subsection "Three hulls": population $A=\operatorname{conv}\operatorname{supp}P_X$; the empirical-law hull (finite-sample identification); the sample hull $\hat A$ as a plug-in. Point forward to P3 Prop. P3.2. | — |
| R-12 | 3.4 | S | 2a §4 (before Theorem 1) | Add a remark that sharpness is relative to the continuous latent-rate model, not to the finite set of integer cell counts. | — |
| R-13 | 4.2 | X | 2a Prop. 5 | Add: the interval is empty unless $0\le Q_{c,rs}$ and $\sum_cQ_{c,rs}=1$ hold, as the probability restrictions require. | — |
| R-14 | 4.4 | X | 2a Cor. 6.1(c) | Replace "strong contextual effects" with the exact condition: part (b)'s common-sign condition fails, which (given the sum restrictions) means at least one pure-poll extrapolated rate is negative. | — |
| R-15 | 3.1, rev. 5.1 | S | P1 §1.5, §1.6, A3; new Phase 2c | Apply D1 (Option B): A3 becomes $E[\beta\mid X,N,Z]=E[\beta\mid X,Z]$, matched by the model $E[\beta\mid X,N,Z]=M_{H(Z)}(X)+\Gamma(Z)$. Phases 2a/2b are the special cases with no $Z$ and with $Z$ = district dummies. Until Phase 2c is written, no text may say the covariate case is proved. | New `covariates =` argument (WS2c) |
| R-16 | 5 | X | 2a Prop. 9 | State the domain: the restrictions on $\Delta$ hold on all of $\Delta^{R-1}$ for the "linear $\Delta=M_{K'}$" characterisation. Move a short version of Prop. 9 into the introduction or §2 as motivation. | — |
| R-17 | 6 | X | 2a new §4.3 | Add "What Theorem 1 and Prop. 7 each need": the sharp set needs the law of $T$ given $(X,N)$, which in practice means smoothing or pooling because every $(X,N)$ is unique. The outer set needs only the regression. Point to the series-based $\Theta_X$ planned for the second release. | — |

### Revised report (adds Phases 2b and 3)

| ID | Reviewer § | Type | Location | Action | Package |
|---|---|---|---|---|---|
| R-18 | 3.2 | X | 2b Prop. D0 | Separate the sharp interval for each cell from the joint set for $B^*$, which is not a Cartesian product. | — |
| R-19 | 3.3 | S | 2b Lemma D1, Remark D1.1; P1 A9 | Replace "$R$ distinct compositions per district" with: the stacked design (pooled quadratic monomials plus district × $X_r$ shifts) has full column rank. Under M1 this requires each district's support to be affinely spanning ($R$ affinely independent points). Raise Remark D1.1 (M1 needs only linear variation within a district) into the main text. | Check rank per group with a message naming the district |
| R-20 | 3.5 | P | 2b Theorem 1′ | List the hypotheses explicitly: $v\mapsto m_\theta(v)$ Borel; uniformly bounded; $\theta\mapsto m_\theta$ affine and continuous; selection set nonempty and compact. | — |
| R-21 | 3.6 | X | 2b Theorem 2(c), ladder figure | Label M1(0) as a *sensitivity benchmark*: it keeps the common skew part but drops the testable symmetric restrictions. | `ladder()` print annotation |
| R-22 | 3.7 | S | 2b Prop. D3 | Make $0<l_d\le u_d<1$ a visible condition. Derive the limiting constraint for districts with pure-composition units ($t\in\{0,1\}$) as a separate case. | Handle $t\in\{0,1\}$ vertices in the R=2 fast path, if one exists (verify) |
| R-23 | 3.8 corr. | **E** | 2b Prop. D4(b); P1 §1.6 (D4) | Condition on $(X,N)$: M0 holds iff $E[Z_D\mid X,N]$ is affine in $X$ and does not depend on $N$. The sufficient condition becomes $D\perp(X,N)$. Add a counterexample to "$D\perp X$ suffices" (district depends on $N$). Qualify "exactly when": the conclusion concerns the mean model, and the outer bounds also need S3 on the pooled region. | — |
| R-24 | 3.9 | X | 2b Remark 6.2 → Proposition | Promote the DiD result to a proposition. State in words that the $NX_rX_s$ and $NX_r$ moments must be equal across periods, not just similar share distributions. | DiD wrappers (second release) cite it |
| R-25 | 4.2 | S | P3 Lemma P3.1 | Add: this holds for a fixed constraint system and objective. It says nothing about joint variation in $(Q,\mu)$ (see R-06). | — |
| R-26 | 4.3 | S/P | P3 Prop. P3.2 | (i) State convergence along a nested sampling sequence, or add an argument for non-nested hulls. (ii) Call Slater a strict-interiority condition and say it fails when an implied rate is 0 at a vertex. Note that P3.5's coverage does not need Slater (step 2 uses only feasibility). (iii) Split the error into estimation error ($O_p(n^{-1/2})$) and hull error (≥ 0, rate not shown). | Optional: report the minimum LP slack as a boundary flag |
| R-27 | 4.4 | **E** (minor) | P3 Theorem P3.4(ii) | The non-Gaussian limit requires at least two active gradients with $(g_\lambda-g_{\lambda'})^\top\Sigma(g_\lambda-g_{\lambda'})>0$. Parametrise only $C-1$ outcome columns, so $\Sigma$ is nonsingular. | Check that $\hat\Sigma$ is built on $C-1$ columns (verify `rxc_core.R`) |
| R-28 | 4.5 | **E** | P3 P3-S1, Theorem P3.5(a) steps 4–5 | Apply D3. Split P3-S1 into (S1a) units independent given the design, and (S1b) districts independent with arbitrary dependence within. P3.5(a) assumes S1a. New **P3.5(c)** under S1b: allowance $S_1^{cl}=\tfrac12(\sum_dW_{dr}^2)^{1/2}$, $W_{dr}=\sum_{i\in d}w_{ir}$; CLT over districts with $D^\star\to\infty$ and $\max_dW_{dr}^2/\sum_dW_{dr}^2\to0$. Cluster Bhatia–Davis version: $\sum_{i\in d}w_ie_i$ lies in an interval of length $W_d$ with mean 0, so its variance is at most $W_d^2\bar m_d(1-\bar m_d)$. Add the reviewer's common-shock example. | **Yes.** With `cluster`, `ci_actual` must use $S_1^{cl}$ (it now uses $\sum_iw_i^2$, see `R/rxc_infer.R:291`). Warn when there are few clusters. |
| R-29 | 4.6 | **E** (claim) | P3 Theorem P3.5(b) | Split it. **(b1) Theorem:** CLR with the full vertex set $\Lambda$, which is finite for a fixed design. **(b2) Algorithm:** collect candidates at simulated draws. **Proof route for a valid (b2):** by weak duality, $\theta_\lambda(\beta_0)\ge\varphi_U(\beta_0)$ for *every* dual-feasible $\lambda$. So $\min_{\lambda\in\hat V}(\hat\theta_\lambda+\hat k\hat\sigma_\lambda)$ covers $\varphi_U$ whenever $\hat k$ controls the maximum over a *fixed* set containing $\hat V$. A data-dependent subset is valid if its critical value is computed over the fixed superset (or the collected set is independent of $\hat\beta$ given the draws, to be checked). Write this as Prop. P3.5(b2) if it closes, otherwise as a Conjecture. List rarely selected vertices, degenerate faces, singular $\Sigma$ and a growing $\lvert\Lambda\rvert$ as the conditions to address. | Per D4: compute the critical value over the union of collected vertices *and* the plug-in vertex; print "algorithmic" when the theorem's conditions are not checked |
| R-30 | 4.7 | **E** (P1) | P1 §1.6; P3 Prop. P3.6(b) | Fix P1's count to $(D^\star-1)(C-1)R(R-1)/2$, with the derivation $(C-1)R(R+1)/2-(C-1)R$. State the Wald/F test for fixed dimension. Under a growing number of restrictions it becomes a conjecture with simulation support. | Check the df in the Wald test code |
| R-31 | 4.8–4.9 | **E** (scope) | P3 §0.1, §6, Prop. P3.7 | Apply D5. Replace "results are for R1/R3" with the reviewer's six-row table of what each model and regime supports. P3.7(a), "dropping constraints enlarges the set", stays a proposition. P3.7(b) becomes a conditional statement: *if* the procedure is valid on the kept districts, choosing them from the design alone keeps conditional coverage. "CLR restores coverage" in R2 becomes a *Simulation finding*. | `eiBounds()` note when the median group size is small (R2 warning) |
| R-32 | 4.3, 9 | S | P3 §3–4 | Central point to add: **coverage of the true parameter does not depend on the hull.** The true $k_0$ satisfies the population constraints, which include those at sample-hull vertices ($\hat A_d\subseteq A_d$). So P3.5 step 2 holds with $\hat A$, and the hull affects only width and the distribution of the endpoint estimator. State P3.4 as a fixed-design result for the sample-hull endpoint $\varphi(\beta_0;\hat A)$. | — |
| R-33 | 4.10 | S/X | P3 §0.2 | Three targets: $B^*$ (superpopulation), $B$ (realised), $\tilde B$ (conditional mean given the design, $=E[B\mid\text{design}]$). Inference is conditional on the design. Covering $B^*$ would also need design-sampling variation, which is left open. Fix "both linear in $K$": only $\tilde B$ is affine in $K$. | Help text for `ci_model`/`ci_actual` |
| R-34 | 4.11 | S | P3 Prop. S1 | Label it an outer bound under bounded deviations, not the sharp set. | `sensitivity()` doc |
| R-35 | rev. 5.2 | **E** (inconsistent) | P1 A1; P3 §7 | Change P1 A1 to say turnover is an open extension, as P3 §7 says. Note that it changes both the DD sets and the regression identity. | — |
| R-36 | rev. 5 | S | P1 A6, A10 | Align A6 with the split S1a/S1b, and A10 with the regime table (R-31). | — |
| R-37 | rev. 6 | X | P3 §8 | Report Monte Carlo standard errors, and use at least 1000 replications for coverage cells. Add designs: (i) a common shock within districts, showing that P3.5(a) under-covers and P3.5(c) does not; (ii) truth on the boundary, so Slater fails; (iii) R2 with growing $D^\star$, reported as simulation findings. | Benchmark/MC scripts |

---

## 2. Workstreams and order

### WS1: Correct outright errors (docs only; first)
R-04, R-06, R-07, R-23, R-27, R-30, R-35, plus the wording fix in R-33.
* For each, write the corrected statement, a one-line counterexample to the old statement, and a verification check.
* **Done when:** every row marked **E** above has a corrected statement, and `verify_*.py` contains a check that fails on the old statement (§3).

### WS2: Scope and assumptions (docs; needs D1, D2)
R-01, R-02, R-05, R-10, R-11, R-12, R-15, R-19, R-22, R-25, R-33, R-36.
* Add a new **P1 §1.7, "Objects and targets"**, cited everywhere: $\Theta_I\subseteq\Theta_X\subseteq\Theta_O$ with the information each uses; the three hulls; the three targets; continuous-rate sharpness.
* Phases 2a, 2b and 3 then point to it, instead of repeating it.

### WS2c: Covariate adjustment, new Phase 2c (needs D1; may run alongside WS3)
Derive for both rungs:
* the identification decomposition: which parts of $\Gamma$ and $G$ are identified, and the unidentified dimension;
* the estimand (the Lemma 2 analogue, with weights that now involve moments of $(N,X,Z)$);
* the outer set. Because $M_{H(Z)}(x)$ is bilinear in $(x,z)$, the vertex reduction of 2a Lemma 3(a) does not carry over to the joint hull of $(X,Z)$. The constraints must be imposed at every observed $(X_i,Z_i)$, or at the vertices of $A\times\operatorname{conv}(\operatorname{supp}Z)$ if a product region is assumed;
* sharpness, via 2b Theorem 1′ with $V=(X,N,Z)$;
* the pooling result (the D4 analogue, conditioning on $(X,N)$): when omitting $Z$ leaves M0 valid.

Also: a rank condition on the design; inference, which reuses P3 unchanged because the target stays affine in the unknowns; a Monte Carlo design with a confounding covariate; and a `covariates =` argument in `eiBounds()`/`ladder()`, adding the covariate rungs to the ladder.

### WS3: Proof hardening (docs; may run alongside WS2)
R-03, R-08, R-09, R-13, R-14, R-16, R-18, R-20, R-26, R-29 (b1 and the weak-duality route for b2), R-34.
* R-09 produces a citation table, with an edition and page for every † tag. Remove the † marks once checked.

### WS4: Restructure inference in Phase 3 (needs WS1, D3–D5)
R-28, R-29, R-31, R-32, and the regime table.

New layout for P3 §5–6:
1. **Theorem P3.5(a):** conservative CI; S1a; covers $\tilde B$ and $B$.
2. **Theorem P3.5(b1):** CLR with the full vertex set.
3. **Proposition/Conjecture P3.5(b2):** the collected candidate set (weak-duality argument).
4. **Theorem P3.5(c):** cluster allowance, under S1b with growing $D^\star$.
5. **Proposition P3.7(a)**, plus the conditional statement P3.7(b).
6. **Table:** what is proved for each model and regime.
7. **§10 Open items:** add R2/R3 theory, design-sampling inference for $B^*$, the hull-error rate, and boundary cases without Slater.

### WS5: Package changes (`eiPartialID`; needs WS4)

| Change | Files | Item |
|---|---|---|
| Cluster allowance ($\tfrac14\sum_dW_d^2$, plus a cluster Bhatia–Davis version) when `cluster` is given. Warn when there are fewer than about 30 clusters. | `R/rxc_infer.R` (`.rxc_fpc_sd`, line ≈291), `R/rxc_api.R` | R-28 |
| Critical value over the fixed superset of candidates. Add an `ndraw` note in the output. Add `ci_method = c("clr", "conservative")`. | `R/rxc_infer.R` (`.rxc_collect`, `.rxc_ci_clr`) | R-29, D4 |
| Messages naming the failing group for rank and dimension problems (quadric-free support for M0/M2; affinely spanning support for M1; full-dimensional hull). | `R/rxc_core.R` | R-02, R-07, R-19 |
| Boundary flag (minimum slack at the optimum) and a small-group (R2) note. | `R/rxc_infer.R`, `R/rxc_api.R` | R-26, R-31 |
| Wording: "outer bounds"; M1(0) as a benchmark; S1 as an outer bound; $\tilde B$ vs $B$ in help. | `R/rxc_api.R` roxygen, vignettes, README | R-10, R-21, R-33, R-34 |
| Check that the Wald test df equals $(D^\star-1)(C-1)R(R-1)/2$, and that $\hat\Sigma$ uses $C-1$ columns. | `R/rxc_infer.R`, `R/rxc_core.R` | R-27, R-30 |
| Update the theory column and inference defaults in `PHASE4_implementation.md`. | docs | all |

**Tests to add:**
* the cluster allowance is at least the unit allowance, with equality for singleton clusters;
* the critical value over the superset is at least the one over the selected set;
* a rank error for degenerate designs;
* the Wald df.

Then: `R CMD check`, and rerun the JKST replication test to confirm it still matches at $10^{-10}$ (the fixes must not change the 2×2 path).

### WS6: Monte Carlo and benchmarks (needs WS5)
R-37.
* Rerun `phase3_montecarlo.py` at 1000 or more replications and report Monte Carlo SEs.
* Add the common-shock, boundary and growing-$D^\star$ designs.
* Rerun the NZ benchmark with the clustered allowance, using electorates as clusters.
* Update `benchmarks/results/ei_datasets_nz_summary.md`.

### WS7: Exposition and the response letter (last)
* Exposition items: R-16 (Prop. 9 earlier), R-17 (Theorem 1 vs Prop. 7 information), R-21, R-24 (DiD proposition).
* Write **`REVIEW_RESPONSE.md`**: for each R-xx, quote the comment, give the response, and point to the new location (file § and result number).
* Update both "What I regard as established" tables with the new status labels. Rows the reviewer marked *Incorrect* or *Not established* must show either the corrected result or a downgrade to Conjecture or Open.

---

## 3. Verification additions

Add these to the existing `verify_*.py` scripts and log them in each phase's verification section:

| Check | Script | Item |
|---|---|---|
| Singleton $A$ with $Kx_0=0$: the LP is unbounded in $t$ | `verify_phase2a.py` | R-07 |
| Singleton $A$ inside the hull: $M_H(x_0)\notin\mathfrak B$ is possible while S1–S2 hold | `verify_phase2a.py` | R-04 |
| $V(Q,\mu)=\mu Q$ is not jointly concave (finite-difference Hessian) | `verify_phase2a.py` | R-06 |
| D4(b): $D\perp X$ but $D$ depends on $N$ ⇒ pooled $E[\beta\mid X,N]$ is not of the form $M_{H'}(X)$ | `verify_phase2b.py` | R-23 |
| M1 df: the rank difference between the M2 and M1 designs equals $(D^\star-1)(C-1)R(R-1)/2$ | `verify_phase2b.py` | R-30 |
| Common shock: $\operatorname{Var}(E)$ exceeds $S_1^2$ and is at most $(S_1^{cl})^2$ | `verify_phase3.py` | R-28 |
| Any dual-feasible $\lambda$ gives $\theta_\lambda(\beta_0)\ge\varphi_U(\beta_0)$ (weak duality, random LPs) | `verify_phase3.py` | R-29 |
| Equal active gradients that differ only in null directions of $\Sigma$ give a Gaussian limit | `verify_phase3.py` | R-27 |

Each new check should **fail on the old statement and pass on the new one**, so the log shows the correction and not only the fix.

---

## 4. Milestones

| Milestone | Contents | Gate |
|---|---|---|
| M1 | Decisions D1–D5 | — |
| M2 | WS1 and its verification checks | All **E** items corrected; verify scripts pass |
| M3 | WS2 and WS3 | The P1 §1.7 objects section exists; citation table complete |
| M4 | WS4 (Phase 3 restructure) | Regime table; every inference result carries a status label |
| M5 | WS5 (package 0.2.1) | `R CMD check` clean; new tests pass; JKST replication unchanged |
| M6 | WS6 | Monte Carlo with SEs; the common-shock design shows the fix |
| M7 | WS7 | `REVIEW_RESPONSE.md` covers R-01…R-37 |

The reviewer's precondition for Phase 3, the link between $\operatorname{supp}P_X$, $A$ and $\hat A$, is closed by R-11 (definitions) and R-32 (the hull does not affect coverage of the parameter). Both are in M3/M4, before any code change that depends on them.
