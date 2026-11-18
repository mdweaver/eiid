# Phase 3: estimation, inference and sensitivity

*Builds on P1, 2a and 2b. Proof steps carry justification tags (2a §0.3, 2b §0.3, plus §0.3 below). Theorems marked "(cited)" apply a published result whose conditions are verified here but whose proof is not repeated. Checks: `verify_phase3.py` (11/11 PASS; `verify_phase3_output.txt`). Monte Carlo: `phase3_montecarlo.py` (`phase3_montecarlo_output.txt`).*

---

## 0. Setting

### 0.1 Sampling

* There are $n$ units with observables $W_i=(N_i,X_i,T_i,D_i)$.
* **(P3-S1)** Units are independent across districts. Within a district they are independent conditionally on the district (rungs M1/M2 carry district parameters). For the pooled model M0, allow arbitrary dependence within districts and use cluster-robust variances.
* **(P3-S2)** The design is conditioned on: inference is conditional on $\{(N_i,X_i,D_i)\}$, so the vertices of the sample hulls $\hat A_d$ are fixed.
* **Regimes** (P1 A10): (R1) $D^\star$ fixed, $n_d\to\infty$; (R2) $D^\star\to\infty$, $n_d$ bounded; (R3) both grow. Results below are for R1/R3 unless stated. R2 is treated in §6.

### 0.2 Targets

There are three related objects (P1 §1.7):
* **Superpopulation rate** $B^*_{rc}=E[NX_r\beta_{rc}]/E[NX_r]$, the object of Phases 2a–2b.
* **Finite-population rate** (JKST's target): $B_{rc}=\sum_iN_iX_{ir}\beta_{irc}\big/\sum_iN_iX_{ir}$, the realised rate for this electorate.
* **Finite-design conditional-mean rate:** $\tilde B_{rc}=\sum_iw_{ir}\,m_{rc}(V_i)$, with $w_{ir}=N_iX_{ir}/\sum_jN_jX_{jr}$ and $m=E[\beta\mid V]$. This is the conditional expectation of $B$ given the observed design (under P3-S1). It is a sample-weighted analogue of $B^*$, not equal to it in general. Under suitable design-sampling conditions it estimates $B^*$, but inference for $B^*$ would also need the variation of the design, which is not pursued here.

All inference below is conditional on the design and covers $B$ or $\tilde B$.

They differ by the residual term
$$E_{rc}=B_{rc}-\tilde B_{rc}=\sum_iw_{ir}e_{irc},\qquad e=\beta-m.$$
$\tilde B$ is the right target when the question is about behaviour rather than the particular electorate. Only $\tilde B$ is affine in $K$, with the *sample* weights $w$, so $\mu$ is never estimated. $B$ is $\tilde B$ plus the residual term $E$, which does not depend on $K$ and is handled separately (§5). *(Corrected after review: an earlier version said both are linear in $K$.)*

### 0.3 Additional tags

| Tag | Meaning |
|---|---|
| [LPdual] | LP strong duality: a feasible LP with bounded value has an optimal dual solution with equal value; the optimum of a pointed polyhedron is attained at a vertex; a polyhedron has finitely many vertices |
| [Slater-b] | under strict feasibility, the optimal dual vectors are bounded: $\lVert\lambda\rVert_1\le(\text{value}-\text{value at the strictly feasible point})/\varepsilon$ |
| [Hausdorff] | Hausdorff/Kuratowski convergence of compact convex sets implies convergence of support functions |
| [dDM] | delta method for Hadamard directionally differentiable maps (Shapiro 1991; Fang & Santos 2019, Thm 2.1†) |
| [minLin] | the directional derivative of $\beta\mapsto\min_{j\in J}g_j^\top\beta$ ($J$ finite) at $\beta_0$, in direction $h$, is $\min_{j\in J^*(\beta_0)}g_j^\top h$, with $J^*$ the active set |
| [CLR] | Chernozhukov, Lee & Rosen (2013, *Econometrica* 81:667–737): parametric bounding functions (Condition P, Lemma 4, Thm 4, via Thms 1–2); finite index set as in their Example 1. *Corrected after review: earlier cited as Thm 2; see 2a §0.4.* |
| [Pop] | Popoviciu: $\operatorname{Var}(Y)\le(b-a)^2/4$ for $Y\in[a,b]$; [BD] Bhatia–Davis: $\operatorname{Var}(Y)\le(M-\mu)(\mu-m)$ |
| [Lyap] | Lyapunov CLT for independent bounded summands with $\max_iw_i^2/\sum_iw_i^2\to0$ |
| [union] | $P(A\cup B)\le P(A)+P(B)$ |
| [OLS] | consistency and asymptotic normality of OLS with heteroskedasticity- or cluster-robust variance, under full rank and finite fourth moments |

---

## 1. Estimators

1. **Regression.** Fit each column $T_c$ by OLS on the homogeneous quadratic monomials $X_rX_s$ ($r\le s$) for M0. For M1, also include the district shifts $1\{D=d\}X_r$, $d\ge2$; this imposes the overidentifying restrictions of 2b Lemma D1(b). For M2, fit one regression per district. The coefficients are $\hat\beta$, and the estimated columns are coherent automatically (2a Remark 1.2). Use an HC3 covariance $\hat\Sigma$ (cluster-robust for M0 with clustered units).
2. **Constraint terms.** For each district $d$, vertex $v\in\operatorname{vert}(\hat A_d)$, row $r$ and column $c$, the constraint is $s_j(\beta)+a_j^\top k\ge0$, with $s_j$ linear in $\beta$ and $a_j$ fixed (2a Prop. 7, 2b Thm 2(b)).
3. **Target.** $\tilde B(k)=c(\beta)+\gamma^\top k$, with $c$ linear in $\beta$ and $\gamma$ computed from the sample weights (2a Lemma 2, 2b Lemma D2).
4. **Plug-in bound.** $[\hat\varphi_L,\hat\varphi_U]=c(\hat\beta)+[\min,\max]\{\gamma^\top k:s(\hat\beta)+Ak\ge0\}$.

---

## 2. Structure of the bound

**Lemma P3.1 (the endpoint is a minimum of finitely many linear functions).** Suppose the LP is feasible with bounded value at $\beta$. Then
$$\varphi_U(\beta)=c(\beta)+\min_{\lambda\in\Lambda}\lambda^\top s(\beta),\qquad\Lambda=\operatorname{vert}\{\lambda\ge0:A^\top\lambda=-\gamma\},$$
a finite set. $\varphi_U$ is therefore concave and piecewise linear in $\beta$, and likewise $\varphi_L(\beta)=c(\beta)-\min_{\lambda\in\Lambda'}\lambda^\top s(\beta)$, with $\Lambda'=\operatorname{vert}\{\lambda\ge0:A^\top\lambda=\gamma\}$.

*Proof.*
1. $\max\{\gamma^\top k:s+Ak\ge0\}$ has Lagrangian $\gamma^\top k+\lambda^\top(s+Ak)$. Stationarity in the free $k$ requires $\gamma+A^\top\lambda=0$, and the dual is $\min\{\lambda^\top s:\lambda\ge0,\ A^\top\lambda=-\gamma\}$. [LPdual]
2. Primal feasibility and boundedness give equal optimal values. [LPdual; 2a:L3(b)]
3. The dual feasible set lies in $\{\lambda\ge0\}$, so it contains no line and is pointed. A linear objective bounded below on it attains its minimum at a vertex, and there are finitely many vertices. [LPdual]
4. The dual feasible set does not depend on $\beta$ ($A$ and $\gamma$ are fixed given the design). So $\varphi_U(\beta)-c(\beta)=\min_{\lambda\in\Lambda}\lambda^\top s(\beta)$ for every $\beta$ at which the primal is feasible and bounded. [3; Def]
5. A minimum of finitely many linear functions is concave and piecewise linear. [alg] The lower endpoint follows by applying 1–5 to $-\gamma$. ∎

Checked: concavity and piecewise linearity along random directions (`verify_phase3.py`).

*Scope.* The lemma is about variation in the regression coefficients $\beta$ for a **fixed** constraint matrix $A$ and objective $\gamma$. Conditioning on the design fixes both, because the target uses the sample weights. It is the first case of 2a Prop. 7's regularity paragraph. It says nothing about joint variation in the identified coefficients, the population moments and the constraint set, where the value function is in general neither concave nor piecewise linear.

**Reading.** JKST's Proposition 4 is this lemma with $|\Lambda|=4$ (their $g_{lj},g_{uj}$). Its "tie-breaking" assumption says the minimising $\lambda$ is unique. In R×C, and especially with many districts, $\Lambda$ is large and near-ties are the norm. That is why the inference below treats the endpoint as an intersection bound.

---

## 3. Consistency and rate

**Proposition P3.2.** Assume (P3-S1–S2), [A1]–[A5] for the rung used, $\sqrt n(\hat\beta-\beta_0)=O_p(1)$, and **(Slater)**: there are $\tilde k$ and $\varepsilon>0$ with $s_j(\beta_0)+a_j^\top\tilde k\ge\varepsilon$ for all constraints defined by the population hulls. Then:
(a) $\hat\varphi_U-\varphi_U(\beta_0;\hat A)=O_p(n^{-1/2})$, and likewise for $\varphi_L$;
(b) $\varphi_U(\beta_0;\hat A)\to\varphi_U(\beta_0;A)$ when $\hat A\to A$ in Hausdorff distance, either along a nested sequence ($\hat A_n\subseteq\hat A_{n+1}$, as when observations are added), or for any Hausdorff-convergent sequence whose feasible sets $\Theta(\hat A_n)$ eventually lie in a common bounded set;
(c) since $\hat A_d\subseteq A_d$, $\varphi_U(\beta_0;\hat A)\ge\varphi_U(\beta_0;A)$ and $\varphi_L(\beta_0;\hat A)\le\varphi_L(\beta_0;A)$: **estimating the hull only ever widens the bound.**

*Proof.*
1. (c): a vertex of $\hat A_d$ is an observed composition, hence in $A_d$. So $\hat A_d\subseteq A_d$, the constraints over $\hat A_d$ are a subset of those over $A_d$, and the feasible set is weakly larger. [Def; 2a:L3(a)]
2. (a): by Lemma P3.1, $\varphi_U(\beta)-\varphi_U(\beta')=c(\beta-\beta')+\big[\min_\lambda\lambda^\top s(\beta)-\min_\lambda\lambda^\top s(\beta')\big]$. [P3.1]
3. The bracket is bounded by $\max\lVert S^\top\lambda\rVert\,\lVert\beta-\beta'\rVert$, the max taken over optimal duals at $\beta$ or $\beta'$. [P3.1; $s$ linear, $s(\beta)=S\beta$]
4. For $\beta$ in a neighbourhood of $\beta_0$, $\tilde k$ remains $\varepsilon/2$-strictly feasible, so the optimal duals are uniformly bounded. [Slater; continuity of $s$; Slater-b]
5. Hence $\lvert\hat\varphi_U-\varphi_U(\beta_0;\hat A)\rvert\le L\lVert\hat\beta-\beta_0\rVert=O_p(n^{-1/2})$, with $L$ uniform in the hull because the constraint set over $\hat A$ is a subset of the one over $A$. [2–4; 1]
6. (b): the sets $\Theta(\hat A)$ are compact and convex, decrease as $\hat A$ grows, and contain $\Theta(A)$. [1; 2a:L3(b)]
7. If $k\in\Theta(\hat A_p)$ for all $p$ and $v\in A$, take $v_p\in\hat A_p$ with $v_p\to v$. The constraints are continuous in $v$, so $k$ satisfies the constraint at $v$. Hence $\bigcap_p\Theta(\hat A_p)=\Theta(A)$. [Def; continuity]
8. A decreasing sequence of compact convex sets converges in Hausdorff distance to its intersection, so the support functions, and hence the LP values, converge. [Hausdorff] ∎(nested)
9. Non-nested case: $\Theta(A)\subseteq\Theta(\hat A_n)$ for every $n$ by 1. If $k_n\in\Theta(\hat A_n)$ and $k_{n_j}\to k$, then for each $v\in A$ pick $v_j\in\hat A_{n_j}$ with $v_j\to v$. Continuity of the constraints in $(k,v)$ gives the constraint at $v$, so $k\in\Theta(A)$. The upper Kuratowski limit is therefore $\Theta(A)$, which is also contained in every term. In a common bounded set this is Hausdorff convergence. [1; Def; continuity; Hausdorff] ∎(b)

**Remarks.**
* **Error decomposition.** $\hat\varphi_U-\varphi_U(\beta_0;A)=\underbrace{[\hat\varphi_U-\varphi_U(\beta_0;\hat A)]}_{O_p(n^{-1/2})\text{ by (a)}}+\underbrace{[\varphi_U(\beta_0;\hat A)-\varphi_U(\beta_0;A)]}_{\ge0,\ \to0\text{ by (b), rate not shown}}$. The second term is not shown to be $O_p(n^{-1/2})$, and can be larger. So Theorem P3.4's limit is stated for the **sample-hull** endpoint $\varphi_U(\beta_0;\hat A)$, conditional on the design, and not for the population-hull endpoint.
* **Slater is a strict-interiority condition.** It requires a $k$ at which every implied rate at every hull vertex is at least $\varepsilon>0$. It fails when every feasible $k$ makes some implied rate exactly 0 at some vertex, e.g. a true structural zero, or the boundary design C of §8, which is not exotic in EI. The Lipschitz bound in (a) is proved only under Slater; boundary cases are open (§10). The coverage results of §5 do **not** use Slater: P3.5(a) step 2 uses only that the true $k_0$ is feasible.
* The hull effect in (c) is a one-sided (conservative) bias. Its size depends on how much mass $P_X$ puts near the boundary of $A_d$; it can vanish slowly. Report the bound's sensitivity to trimming the hull (P1 A4).
* The $O_p(n^{-1/2})$ rate is for the estimation part. The target $B$ adds $E=O_p(n^{-1/2})$ as well (§5). Bounds do not shrink to a point: their width converges to the population width, which is $O(1)$.

---

## 4. Asymptotic distribution

**Theorem P3.4 (cited structure: [dDM]).** Under Proposition P3.2's assumptions and $\sqrt n(\hat\beta-\beta_0)\Rightarrow\mathcal Z\sim N(0,\Sigma)$, with the design and the sample hull $\hat A$ fixed, so that $\varphi_U(\beta_0)$ means $\varphi_U(\beta_0;\hat A)$ (see the error decomposition after P3.2):
$$\sqrt n\big(\hat\varphi_U-\varphi_U(\beta_0)\big)\Rightarrow\min_{\lambda\in\Lambda^*}\big(c+S^\top\lambda\big)^\top\mathcal Z,$$
where $\Lambda^*\subseteq\Lambda$ is the set of optimal dual vertices at $\beta_0$.
* (i) If $\Lambda^*$ is a singleton (JKST's tie-breaking), the limit is $N(0,g^{*\top}\Sigma g^*)$ with $g^*=c+S^\top\lambda^*$.
* (ii) If two active gradients differ in a direction in which the estimator varies, i.e. $(g_\lambda-g_{\lambda'})^\top\Sigma(g_\lambda-g_{\lambda'})>0$ for some $\lambda,\lambda'\in\Lambda^*$, the limit is a minimum of correlated normals that are not a.s. equal, and it has strictly negative mean. The plug-in upper endpoint is then biased *downward*, and the lower endpoint upward, both by order $n^{-1/2}$. That is, the plug-in bound is too narrow. If instead all active gradients agree on the support of $\mathcal Z$, all the $g_\lambda^\top\mathcal Z$ are a.s. equal and the limit is normal, as in (i), even though $\Lambda^*$ is not a singleton. *(Corrected after review: an earlier version claimed a negative mean whenever $\Lambda^*$ is not a singleton.)*
* (iii) In case (ii) the standard bootstrap is inconsistent (Fang & Santos 2019, Cor. 3.1†).

*Proof.*
1. $\varphi_U(\beta)=\min_{\lambda\in\Lambda}g_\lambda^\top\beta$ with $g_\lambda=c+S^\top\lambda$: a finite minimum of linear maps. [P3.1]
2. Such a map is Lipschitz, and its directional derivative is $h\mapsto\min_{\lambda\in\Lambda^*}g_\lambda^\top h$. A Lipschitz map with directional derivatives in finite dimensions is Hadamard directionally differentiable. [minLin]
3. Apply the directional delta method. [dDM; 1–2]
4. (i): the minimum over a singleton is linear, giving a normal limit. [3]
5. (ii): let $Y=g_\lambda^\top\mathcal Z$ and $Y'=g_{\lambda'}^\top\mathcal Z$ with $\operatorname{Var}(Y-Y')=(g_\lambda-g_{\lambda'})^\top\Sigma(g_\lambda-g_{\lambda'})>0$. Then $E\min_{\Lambda^*}g^\top\mathcal Z\le E\min(Y,Y')=\tfrac12E(Y+Y')-\tfrac12E\lvert Y-Y'\rvert=-\tfrac12E\lvert Y-Y'\rvert<0$. If $\operatorname{Var}(g_\lambda^\top\mathcal Z-g_{\lambda'}^\top\mathcal Z)=0$ for all active pairs, the minimum equals each $g_\lambda^\top\mathcal Z$ a.s., which is normal. [3; $\min(a,b)=\frac{a+b}2-\frac{\lvert a-b\rvert}2$; $Y-Y'$ is a nondegenerate normal]
6. (iii): in case (ii) the derivative is nonlinear on the support of $\mathcal Z$, so the condition for bootstrap consistency fails. [Fang & Santos 2019†] ∎

**Remark P3.4.1 (parametrise $C-1$ columns).** The shares satisfy $\sum_cT_c=1$, so the coefficient vectors of the $C$ column regressions sum to a known constant, and the joint covariance of all $C$ columns is singular. State $\hat\beta$, $\Sigma$ and the gradients $g_\lambda$ for the first $C-1$ columns and recover column $C$ by adding up. Then $\Sigma$ can be nonsingular. When it is positive definite, the condition in (ii) reduces to $g_\lambda\neq g_{\lambda'}$ in the reduced parametrisation.

Checked: with two tied constraints and unit noise, the mean plug-in error is $-0.564\approx-1/\sqrt\pi$; with a unique active constraint it is 0.007. With two active gradients that differ only in a null direction of $\Sigma$, the mean is 0 (normal limit), as the corrected (ii) states.

**Boundary cases.**
* **Truth at a vertex.** When the true $K$ lies on the boundary of the identified set, the side of the interval that touches it carries the full error probability, so coverage is close to nominal. When the true $K$ is interior, coverage tends to 1 (JKST App. A Remark). This is the Imbens–Manski distinction between covering the parameter and covering the set.
* **Width near zero.** When the bound's width is near zero (point identification: L4 in 2×2, L5, the DiD with unchanged contextual effects), one-sided intervals per endpoint under-cover unless the critical value adapts (Stoye 2009). The two-sided Bonferroni split used below is uniformly valid, at some cost in width.

---

## 5. Confidence intervals

**Theorem P3.5 (validity).** Fix $\alpha_U,\alpha_L,\alpha_E\in(0,1)$. Write $\theta_\lambda(\beta)=g_\lambda^\top\beta$ for the candidate upper endpoints ($\lambda\in\Lambda$), with standard errors $\hat\sigma_\lambda$ from $\hat\Sigma$.

**(a) Conservative version.** Let $\hat k$ be the $(1-\alpha_U-\alpha_L)$-quantile of $\max_j$ over all constraint terms and the constant of the standardised Gaussian errors (one-sided for the constraints, two-sided for the constant). Put
$$\hat\Theta^+=\{k:s(\hat\beta)+\hat k\hat\sigma_s+Ak\ge0\},$$
$$\mathrm{CI}=\Big[c(\hat\beta)-\hat k\hat\sigma_c+\min_{\hat\Theta^+}\gamma^\top k-z_{\alpha_E}S_1,\ \ c(\hat\beta)+\hat k\hat\sigma_c+\max_{\hat\Theta^+}\gamma^\top k+z_{\alpha_E}S_1\Big],$$
with $S_1=\tfrac12\big(\sum_iw_{ir}^2\big)^{1/2}$ and $z_{\alpha_E}=\Phi^{-1}(1-\alpha_E/2)$. Then
$$\liminf_nP(B_{rc}\in\mathrm{CI})\ge1-\alpha_U-\alpha_L-\alpha_E.$$

**(b1) CLR with the full candidate set (cited: [CLR], Thm 4).** Let $\Lambda$ be the full set of dual vertices of Lemma P3.1, which is finite for a fixed design. Assume CLR's Condition P: $\Sigma$ is positive definite (use $C-1$ columns, Remark P3.4.1, and a full-rank design), and $g_\lambda\ne0$ for every $\lambda\in\Lambda$. For the upper end:
1. Select $\hat V=\{\lambda:\hat\theta_\lambda\le\min_\mu(\hat\theta_\mu+k_p\hat\sigma_\mu)+2k_p\hat\sigma_\lambda\}$, where $k_p$ is the $\gamma_n=1-0.1/\log n$ quantile of $\max_\lambda Z_\lambda$, and $Z$ is the standardised Gaussian vector with correlation from $\hat\Sigma$.
2. Set $\hat U=\min_{\lambda\in\hat V}(\hat\theta_\lambda+\hat k\hat\sigma_\lambda)$, with $\hat k$ the $(1-\alpha_U)$-quantile of $\max_{\lambda\in\hat V}Z_\lambda$.
3. Do the same for the lower end, and add $z_{\alpha_E}S_1$ on each side (each side's E-allowance at level $\alpha_E/2$).

The same coverage bound holds. With $\alpha_E=0$ and no $S_1$ term, either version covers $\tilde B_{rc}$ with probability at least $1-\alpha_U-\alpha_L$.

**(b2) Any candidate set, with a simultaneous critical value (proved here).** Let $\hat\Lambda$ be any data-dependent set of dual-feasible points (vertices or not; for example collected by simulation), and $\hat U=\min_{\lambda\in\hat\Lambda}(\hat\theta_\lambda+\hat k\hat\sigma_\lambda)$, with $\hat U=+\infty$ if $\hat\Lambda=\emptyset$. Then $\liminf_nP(\varphi_U(\beta_0)\le\hat U)\ge1-\alpha_U$ in either of two cases:
* (i) $\hat k$ is the $(1-\alpha_U)$-quantile of $\max_{\lambda\in\bar\Lambda}Z_\lambda$ for a fixed finite set $\bar\Lambda\supseteq\hat\Lambda$ a.s. (e.g. all of $\Lambda$), with $\sigma_\lambda>0$ on $\bar\Lambda$;
* (ii) projection: $\hat k=\sqrt{\chi^2_{q,1-\alpha}}$ with $q=\operatorname{rank}\Sigma$. This one event covers both ends, any $\hat\Lambda$, at total level $\alpha$.

*Proof.*
1. Under the model the primal LP at $\beta_0$ is feasible ($k_0$ is feasible). For any dual-feasible $\lambda$ and any primal-feasible $k$: $\gamma^\top k=-\lambda^\top Ak\le\lambda^\top s(\beta_0)$. So $\theta_\lambda(\beta_0)\ge\varphi_U(\beta_0)$ for **every** dual-feasible $\lambda$, not only optimal ones. [weak duality; $\lambda\ge0$, $s+Ak\ge0$]
2. On $\mathcal E=\{\theta_\lambda(\beta_0)\le\hat\theta_\lambda+\hat k\hat\sigma_\lambda\ \forall\lambda\in\bar\Lambda\}$: $\hat U\ge\min_{\lambda\in\hat\Lambda}\theta_\lambda(\beta_0)\ge\varphi_U(\beta_0)$. [1; $\hat\Lambda\subseteq\bar\Lambda$]
3. (i): $(\theta_\lambda(\beta_0)-\hat\theta_\lambda)/\hat\sigma_\lambda$, $\lambda\in\bar\Lambda$, is a fixed finite collection of standardised linear functionals of $\hat\beta$, jointly asymptotically $N(0,\text{corr})$. So $P(\mathcal E)\to1-\alpha_U$. [OLS; $\hat\Sigma$ consistent; Def $\hat k$]
4. (ii): by Cauchy–Schwarz in the metric of $\Sigma$, $\lvert g^\top(\hat\beta-\beta_0)\rvert\le\lVert\Sigma^{+1/2}(\hat\beta-\beta_0)\rVert\,(g^\top\Sigma g)^{1/2}$ for every $g$ in the range of $\Sigma$. The first factor, scaled by $\sqrt n$, is asymptotically $\chi_q$. So $\mathcal E$ holds for all $\lambda$ at once, and for both ends, with probability $\to1-\alpha$. [alg; OLS] ∎

(ii) needs no enumeration of vertices and is valid for any collection rule. It is conservative, increasingly so as $q$ grows (many districts). (i) is sharper, but needs a fixed superset, and enumerating all dual vertices is feasible only in small problems. Checked (`verify_phase3.py`): weak duality on 500 random dual-feasible points. In a small LP, the projection bound built from the plug-in dual alone covered in 100% of 3000 replications, against 94.5% for a one-candidate $z$-bound.

**Conjecture P3.5(b3) (the implemented CLR procedure).** The package collects $\hat\Lambda$ as the optimal dual vertices at $\hat\beta$ and at `ndraw` Gaussian draws, and computes both $k_p$ and $\hat k$ over $\hat\Lambda$. The conjecture is that it has asymptotic coverage $\ge1-\alpha_U$ under (b1)'s conditions as `ndraw` $\to\infty$.
* *Proved reduction.* Suppose $k_p$ is computed over $\Lambda$, and $\hat\Lambda$ contains every vertex selected by CLR's rule over $\Lambda$ and a minimiser of $\hat\theta_\mu+k_p\hat\sigma_\mu$ over $\Lambda$. Then the threshold, the selected set, $\hat k$ and $\hat U$ all coincide with (b1), so (b1)'s coverage carries over. (The threshold is the same because its minimiser is in $\hat\Lambda$. The selected set over $\hat\Lambda$ is CLR's selected set intersected with $\hat\Lambda$, which is CLR's selected set.)
* *Open.* (1) Recovery: vertices that are optimal only for rare Gaussian directions, degenerate optimal faces, and singular $\Sigma$. (2) Computing $k_p$ over $\hat\Lambda\subseteq\Lambda$ can only lower it, so this error goes in the anti-conservative direction. (3) $\lvert\Lambda\rvert$ grows with the number of districts, outside CLR's fixed-index setting.

Until this is settled, the CLR-type interval as implemented is an *Algorithm* with simulation support (§8). The intervals of (a) and (b2) carry the formal guarantee.

*Proof of (a).*
1. Let $\mathcal E$ be the event that $s_j(\beta_0)\le s_j(\hat\beta)+\hat k\hat\sigma_{s_j}$ for all $j$ and $\lvert c(\beta_0)-c(\hat\beta)\rvert\le\hat k\hat\sigma_c$. Then $P(\mathcal E)\to\ge1-\alpha_U-\alpha_L$. [OLS; linear maps of an asymptotically normal $\hat\beta$; Def of $\hat k$; $\hat\Sigma$ consistent]
2. On $\mathcal E$, the true $k_0$ satisfies $s(\beta_0)+Ak_0\ge0$, hence $s(\hat\beta)+\hat k\hat\sigma_s+Ak_0\ge0$, i.e. $k_0\in\hat\Theta^+$. [2b:Thm 2 / 2a:P1, the true value is feasible; $\mathcal E$]
3. On $\mathcal E$, $\tilde B=c(\beta_0)+\gamma^\top k_0$ lies in $[c(\hat\beta)-\hat k\hat\sigma_c+\min_{\hat\Theta^+}\gamma^\top k,\ c(\hat\beta)+\hat k\hat\sigma_c+\max_{\hat\Theta^+}\gamma^\top k]$. [2; $\mathcal E$; 2b:D2 with sample weights]
4. $E=\sum_iw_ie_i$ with $e_i$ independent given the design, mean 0, and $\operatorname{Var}(e_i\mid V_i)\le\tfrac14$. So $\operatorname{Var}(E\mid\text{design})\le S_1^2$. [P3-S1; CE; Pop]
5. $P(\lvert E\rvert>z_{\alpha_E}S_1)\to\le\alpha_E$, provided $\max_iw_i^2/\sum_iw_i^2\to0$. [Lyap; 4]
6. If $B\notin\mathrm{CI}$, then either $\tilde B$ lies outside the interval in 3, or $\lvert E\rvert>z_{\alpha_E}S_1$. [$B=\tilde B+E$; triangle inequality]
7. So $P(B\notin\mathrm{CI})\le P(\mathcal E^c)+P(\lvert E\rvert>z_{\alpha_E}S_1)\to\le\alpha_U+\alpha_L+\alpha_E$. [union; 1, 3, 5, 6] ∎

*(b1)* replaces step 1 by CLR's coverage of each intersection bound, which needs the finite index set $\Lambda$ (Lemma P3.1) and their Condition P [CLR]. Steps 4–7 are unchanged. *(b2)* is proved above. The simulation-based collection of candidates is the subject of Conjecture P3.5(b3); it is not claimed as a theorem.

**Remarks.**
* **Better variance bound for the residual term.** $\operatorname{Var}(e_{irc}\mid V)\le m(1-m)$ [BD], with $m=E[\beta_{irc}\mid V]$, which is bounded through the identified set. Taking the maximum of $m(1-m)$ over that set's range at $X_i$ gives a smaller, still valid $S_1$.
* **Optional assumption for a much smaller residual term.** If the poll-level deviations of different groups are not negatively correlated, then $\operatorname{Var}(X_{ir}e_{irc})\le\operatorname{Var}(T_{ic}\mid V_i)$, which is estimable. In simulation this would remove most of the residual allowance. It is an added assumption and belongs on the menu, not in the default.
* **JKST's $\mathrm{CI}_x$** is version (b) with a single candidate per side (tie-breaking), no selection step, and $S_L=S_1+\hat\sigma_{\lambda^*}$, which combines the two error sources by adding standard deviations instead of by the union bound.

**Proposition P3.6 (tests).**
(a) **Emptiness test.** Reject the model (rung, foundations) at level $\alpha_U+\alpha_L$ if $\hat\Theta^+$ is empty.
(b) **Test of M1's testable part.** Run a Wald/F test of the M1 regression against M2 (the overidentifying restrictions of 2b D1(b)), with $(D^\star-1)(C-1)R(R-1)/2$ restrictions. In 2×2: equal curvature across districts.
*Proof of (a).*
1. Under the model, $k_0\in\hat\Theta^+$ on $\mathcal E$. [P3.5(a) step 2]
2. So $P(\hat\Theta^+=\emptyset)\le P(\mathcal E^c)\to\le\alpha_U+\alpha_L$. [1]

(b) is the standard Wald test [OLS]. Under R2 the number of restrictions grows with $D^\star$, so use the F-form with HC3 or a wild bootstrap. ∎

---

## 6. Many districts with few units (regime R2)

**Proposition P3.7 (dropping constraints is conservative; design-based selection is valid).** For any subset $\mathcal D'$ of districts, $\Theta^{M1}_O(\mathcal D')\supseteq\Theta^{M1}_O$. If $\mathcal D'$ is chosen as a function of the design $\{(N_i,X_i,D_i)\}$ only, every result of §§3–5 holds conditionally on $\mathcal D'$.
*Proof.*
1. Fewer constraints give a weakly larger feasible set. [Def; 2b:Thm 2(b)]
2. Conditioning on the design fixes $\mathcal D'$, so the proofs apply with $\mathcal D'$ in place of all districts. [P3-S2] ∎

**What happens under R2.**
* Each district's shifts $\hat\zeta_d$ come from about $n_d$ units, and their sampling error does not shrink.
* The number of candidates grows with $D^\star$.
* The plug-in is badly too narrow: in the Monte Carlo (§8) it covers almost never, and is empty about 40% of the time. This is the Theorem P3.4(ii) bias, compounded across many districts.
* The CLR version restores coverage by inflating each candidate in proportion to its own standard error. Noisy districts therefore rarely bind, which is precision weighting built into the method.
* The price is width. With 10 units per district, the valid interval is about as wide as Duncan–Davis, so the identification gain from district dummies (2b) is consumed by sampling noise.

Practical options, all valid:
* (i) Keep the constraints of districts with enough units or composition spread (P3.7).
* (ii) Use M1($\delta$), $\delta>0$, which bounds how far any one district can pull the shared parameter.
* (iii) Pool districts into larger groups.
* (iv) Target $\tilde B$ to drop the residual term.

A hierarchical restriction on the district shifts (they are draws from a distribution) would change the model and is left open.

---

## 7. Sensitivity analysis

**Proposition S1 (bounded nonlinearity).** Replace (S2) by $E[\beta\mid V]=M_H(X)+\Delta(X)$, with $\Delta$ in the undetectable class of 2a Prop. 9 ($x^\top\Delta(x)=0$, $\Delta(x)\mathbf 1=0$) and $\lVert\Delta(x)\rVert_\infty\le\eta$. Then $Q$ is still identified as before, and every value of $\tilde B_{rc}$ consistent with the data lies in
$$\Big[\min_{\Theta_O(\eta)}\varphi-\eta,\ \max_{\Theta_O(\eta)}\varphi+\eta\Big],\qquad\Theta_O(\eta)=\{k:s+Ak\ge-\eta\mathbf 1\}.$$
For a general $\varphi=\sum a_{rc}B_{rc}$, widen by $\eta\sum\lvert a_{rc}\rvert$ instead. This is an **outer** bound under bounded deviations: it is valid, but it is not claimed to be the sharp set of the relaxed model.

*Proof.*
1. $E[T\mid X]=X^\top Q X+X^\top\Delta(X)=X^\top QX$, so $Q$ is unchanged. [2a:L1(a); $x^\top\Delta=0$]
2. $M_H(x)+\Delta(x)\in\mathfrak B$ implies $M_H(x)\ge-\Delta(x)\ge-\eta$ elementwise on $A$. $M_H$ is linear, so this is equivalent to the same inequality at the vertices: $s+Ak\ge-\eta$. [S3; 2a:L3(a)]
3. $\tilde B_{rc}=\sum_iw_i[M_H(X_i)+\Delta(X_i)]_{rc}=\varphi(k)+\sum_iw_i\Delta_{rc}(X_i)$, and the last term lies in $[-\eta,\eta]$ because the weights sum to 1. [2b:D2; Def] ∎

Checked: a nonlinear undetectable violation with $\eta=0.02$ stays inside the bound, and the bounds are nested in $\eta$.

**Breakdown value.** For a conclusion $\mathcal C$, e.g. "loyalty > 0.30", define $\eta^*(\mathcal C)=\sup\{\eta:\mathcal C\text{ holds for every value in the }\eta\text{-bound}\}$. The bounds are nested in $\eta$, so the conclusion holds for all $\eta<\eta^*$ and fails above it, and $\eta^*$ can be found by bisection. In the 3×3 example, $\eta^*(\text{loyalty}>0.30)=0.060$: the conclusion survives any undetectable nonlinearity of up to 6 percentage points. Apply the CI version of the bound to account for sampling error.

**Other dials, same logic** (nested sets, bisection for breakdown):
* $\delta$ for district heterogeneity (2b Prop. D6);
* $\delta$ for change between periods (menu item 19);
* population turnover $\epsilon$ (A1). Turnover requires relaxing the Duncan–Davis sets and the regression identity together. It is left as an open item, with interim advice to add explicit entry/exit rows.

---

## 8. Monte Carlo

`phase3_montecarlo.py`, 200 replications per design. Unit rates are Beta/Dirichlet draws whose conditional means are exactly linear, so the model holds. Each design reports the plug-in bound, the conservative CI (Thm P3.5(a)), the CLR CI (P3.5(b)) for $B$, and the CLR CI for $\tilde B$, all at nominal 95%.

Coverage is the share of 200 replications containing the truth; width is the median. The DD width is shown for reference; every CI may be intersected with DD.

| Design | DD width | Plug-in: cover / empty / width | Conservative CI for $B$ | CLR CI for $B$ (Popoviciu) | CLR CI for $B$ (Bhatia–Davis) | CLR CI for $\tilde B$ |
|---|---|---|---|---|---|---|
| A: 2×2, M1, 10 districts × 93 | 0.62 | 0.970 / 0.00 / 0.282 | 1.000 / 0.681 | 1.000 / 0.522 | 1.000 / 0.531 | 1.000 / 0.424 |
| A: 2×2, M1, 93 × 30 | 0.56 | 0.665 / 0.00 / 0.158 | 1.000 / 0.702 | 1.000 / 0.476 | 1.000 / 0.481 | 1.000 / 0.419 |
| A: 2×2, M1, 93 × 10 | 0.56 | **0.055 / 0.44** / 0.049 | 1.000 / 1.048 | 1.000 / 0.571 | 1.000 / 0.580 | 1.000 / 0.471 |
| C: 2×2 pooled, truth on the boundary | 0.33 | 0.660 / 0.11 / 0.031 | 1.000 / 0.140 | 1.000 / 0.135 | 1.000 / **0.077** | **0.990** / 0.042 |
| B: 3×3 pooled, loyalty P→P | 0.84 | 1.000 / 0.00 / 0.536 | 1.000 / 0.689 | 1.000 / 0.653 | 1.000 / 0.661 | 1.000 / 0.560 |
| B: 3×3 pooled, defection P→O | 0.64 | 1.000 / 0.00 / 0.352 | 1.000 / 0.498 | 1.000 / 0.466 | 1.000 / 0.471 | 1.000 / 0.373 |
| B: 3×3 pooled, demobilisation | 0.53 | 1.000 / 0.00 / 0.283 | 1.000 / 0.428 | 1.000 / 0.397 | 1.000 / 0.396 | 1.000 / 0.303 |
| B: 3×3 pooled, mobilisation N→P | 0.89 | 1.000 / 0.00 / 0.530 | 1.000 / 0.720 | 1.000 / 0.660 | 1.000 / 0.668 | 1.000 / 0.561 |

**Reading the results.**

1. **The plug-in is not a confidence interval and fails exactly where Theorem P3.4 predicts.**
   * With many districts and few units it almost never covers: 5.5% at 93 × 10, where it is also empty 44% of the time. At 93 × 30 it covers 66.5%.
   * It covers 66% in the boundary design.
   * It is fine only with a few large districts, or pooled data with the truth in the interior.
2. **Both CIs are valid in every design**, including the boundary case. That case is the informative one: with the truth at an endpoint, the CLR CI for $\tilde B$ covers 99%, close to its nominal level. The 100% elsewhere reflects truths in the interior of the identified set (§4).
3. **CLR is much narrower than the conservative version**: 0.571 vs 1.048 at 93 × 10, for example. It should be the default.
4. **The residual term is a large part of the width.** Targeting $\tilde B$ instead of $B$ narrows the 3×3 intervals by 0.09–0.10 and the 2×2 M1 intervals by 0.06–0.10.
5. **With 10 units per district the valid CI for $B$ is no narrower than DD** (0.571 vs 0.56), and for $\tilde B$ it is 0.471. The identification gain of M1 is mostly consumed by sampling noise in this regime (§6). With 93 units per district, the CI for $\tilde B$ is 0.424 against DD's 0.62.
6. **Bhatia–Davis versus Popoviciu.** Bhatia–Davis narrows the actual-rate CI sharply when the plausible range of rates excludes 0.5: in the boundary design it is 0.077 against 0.135. Elsewhere that range includes 0.5, so the allowance equals Popoviciu's and the interval is about 0.01 wider, the price of the α share spent bounding the variance. It remains the package default: it is never worse by more than that small share and can be much better.
7. **In the pooled 3×3 design the CIs are 0.30–0.56 wide, against DD's 0.53–0.89.** A real gain, but these are wide intervals. Narrowing them further is the role of the stronger rungs and the menu restrictions.

---

## 9. Practical recipe (input to Phase 4)

1. Fit the rung's regression (M0 / M1 / M2) with HC3 or cluster-robust covariance, and the sample hulls.
2. Build the LP (2a Prop. 7 / 2b Thm 2) and report the plug-in bound, labelled as an estimate of the identified set, not a confidence interval.
3. Report the CLR confidence interval for $\tilde B$ (default) and for $B$, intersected with Duncan–Davis (valid because $B\in DD$ always).
4. Report the emptiness test and, for M1, the test of its testable part.
5. Report the sensitivity curves in $\eta$ (and $\delta$), with the breakdown value of the headline conclusion.
6. With few units per district, apply P3.7 selection, or prefer M1($\delta$) or coarser groups, and say so.

## 10. Open items

* Conjecture P3.5(b3): validity of collecting dual vertices by simulation, with $k_p$ computed over the collected set. A proved reduction to (b1) is given; recovery and the $k_p$ effect are open.
* Boundary cases where Slater fails (implied rates exactly 0 at a vertex): rate and limit of the endpoint estimator.
* The rate of the hull error $\varphi(\beta_0;\hat A)-\varphi(\beta_0;A)$, and inference for the population-hull endpoint.
* Inference for the superpopulation target $B^*$, which needs design-sampling variation.
* Uniform validity near point identification for the CLR version (Stoye-type adaptation).
* Inference under R2 with a hierarchical model for the district shifts.
* Joint relaxation for population turnover.
* Coverage of the plug-in under misspecification, with the menu's dials.
