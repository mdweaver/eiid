# Phase 1: Structural foundations and assumptions

*Extending Jiang, King, Schmaltz & Tanner (2019), "Ecological Regression with Partial Identification" (JKST below), and the `eiPartialID` package to (1) R×C tables and (2) covariates, especially district fixed effects. This is a draft for review. Results marked **(to prove)** are stated precisely here and get proved in Phase 2. Nothing in this phase is a derivation.*

---

## 0. What the source materials actually do (and what that means for "sharpness")

I read the paper, the supplementary appendix (`appendix.tex`, Appendices A–C) and `package/eiPartialID/R/main.R`. Six facts shape everything below.

1. **Model.** JKST Assumption 1: $(\beta^b_i,\beta^w_i,X_i,N_i)$ are iid, with $E(\beta^w_i\mid X_i,N_i)=w_0+w_1X_i$ and $E(\beta^b_i\mid X_i,N_i)=b_0+b_1X_i$. The accounting identity then implies $E(T_i\mid X_i)=w_0+c_1X_i+d_1X_i^2$. That identifies $\theta=(w_0,c_1,d_1)$ and leaves one free parameter, $w_1$.
2. **Prop. 1 (JKST)** bounds $w_1$ through $E(L_i\mid X_i)\le w_0+w_1X_i\le E(U_i\mid X_i)$ on a set $A$. Prop. 2 relaxes this with Jensen's inequality and the quadratic regression. **The package implements only Prop. 2** with $A=[\min X_i+10^{-5},\ \max X_i-10^{-5}]$ and $\lambda=0$, so the package's bound is an *outer* set, not the sharp one.
3. **Inference (Prop. 4, App. A).** $\mathrm{CI}_x=[\hat B_L-xS_L,\ \hat B_U+xS_U]\cap DD$, with HC1 sandwich variance and the Popoviciu bound $\operatorname{Var}(\beta)\le 1/4$ for the finite-population term ($S_1$). It needs tie-breaking (a unique active constraint) and covers the *parameter* $B$ with probability at least $\Phi(x)$, in the spirit of Imbens–Manski.
4. **Selection heuristic.** The package falls back to DD when $\mathrm{CI}_0$ intersected with DD is empty. App. B justifies this: under the model, $\mathrm{CI}_0\neq\emptyset$ in the large-$p$ limit.
5. **Covariates (App. C).** These enter additively in the intercept: $E(\beta^w\mid X,Z)=\tilde w_0+\tilde w_1X+\tilde w_2Z$, and likewise for $\beta^b$, with a common slope. Five parameters are identified and only $\tilde w_1$ is not. The appendix bounds it by restricting $A$ to a line segment $Z=a_0+a_1X$ (Prop. 3\*). That choice throws away most of what a covariate can contribute. Apart from the $\log N$ example, it does not analyse *discrete* covariates such as district dummies. None of App. C is in the package.
6. **Verified equivalence.** I rewrote Prop. 2 in the homogeneous form used below (§1.4) and checked numerically that it reproduces `generateBounds` to $10^{-16}$ (`docs/extension/check_2x2_equivalence.py`). The R×C generalisation therefore contains the existing package as the special case $R=C=2$.

**Consequence for the brief.** You asked for a proof that "the bounds" are sharp. Under JKST's assumptions the bound the package computes is **not** sharp. There are three nested sets, $\Theta_I\subseteq\Theta_X\subseteq\Theta_O$ (§1.5), and Phase 2 will prove sharpness for $\Theta_I$, prove validity for $\Theta_O$, and characterise when they coincide. I will not claim sharpness for $\Theta_O$.

---

## 1. Identification regime

### 1.1 Sampling scheme

There are units (precincts) $i=1,\dots,p$ nested in districts $d=1,\dots,D$. $D_i\in\{1,\dots,D\}$ is unit $i$'s district, $\mathcal I_d=\{i:D_i=d\}$, and $p_d=|\mathcal I_d|$. Each unit gives
$$
W_i=(N_i,\;X_i,\;T_i,\;Z_i,\;D_i),\qquad X_i\in\Delta^{R-1},\quad T_i\in\Delta^{C-1},
$$
where $\Delta^{k-1}=\{v\in[0,1]^k:\mathbf 1^\top v=1\}$, $X_{ir}=N_{ir}/N_i$ is the row margin (e.g. vote shares in election 1, including abstention), and $T_{ic}$ is the column margin (e.g. shares in election 2: abstain, party A, party B, …). $Z_i\in\mathcal Z\subseteq\mathbb R^{k_z}$ are other covariates, such as $\log N_i$ or census variables.

**Latent quantities.** These are the unit cell counts $N_{irc}$ (the unobserved $R\times C$ table) and the transition probabilities $\beta_{irc}=N_{irc}/N_{ir}$, defined when $N_{ir}>0$. Write $\beta_i=[\beta_{irc}]\in\mathfrak B:=(\Delta^{C-1})^R$ (each row is a point in the simplex).

**Accounting identity (definitional for a closed population, see A1):**
$$
T_{ic}=\sum_{r=1}^R X_{ir}\beta_{irc}\quad(c=1,\dots,C),\qquad\text{i.e.}\quad T_i=\beta_i^\top X_i. \tag{1.1}
$$

### 1.2 Estimands

* **Finite-population (district- or region-wide) transition rates** for any target set of units $\mathcal S$ (all units, or $\mathcal I_d$):
  $$B_{rc}(\mathcal S)=\frac{\sum_{i\in\mathcal S}N_iX_{ir}\beta_{irc}}{\sum_{i\in\mathcal S}N_iX_{ir}}.\tag{1.2}$$
  The JKST $B$ is $B_{\text{black},\text{vote}}(\{1..p\})$.
* **Superpopulation analogue** (the object that is *identified*; $B-B^*=O_p(p^{-1/2})$, Phase 3):
  $$B^*_{rc}=\frac{E[N X_r\beta_{rc}]}{E[NX_r]}.\tag{1.3}$$
* **Transition functionals** for your application ($R=2$ rows A, B; $C=3$ columns $0$ = abstain, A, B): loyalty $\lambda_r=B_{r,r}$, defection $\delta_r=B_{r,r'}$ ($r'\neq r$, $r'\neq 0$), demobilisation $\mu_r=B_{r,0}$. "Same/different party" is a row-dependent relabelling of columns, so every quantity of interest is a **linear functional** $\varphi=\sum_{rc}a_{rc}B_{rc}$. Since $\sum_c B_{rc}=1$, the three are linked by $\lambda_r+\delta_r+\mu_r=1$.

### 1.3 Structural model: linear contextual effects in homogeneous form

Because $\mathbf 1^\top X=1$, an affine function of $X$ is a *homogeneous linear* function of $X$: $a+\sum_{s<R}g_sX_s=\sum_s(a+g_s)X_s$ with $g_R:=0$. This removes the intercept/slope bookkeeping and makes the R×C algebra symmetric.

**Model M0 (pooled LCE, generalising JKST Assumption 1).** For each column $c$ there is a matrix $H_c\in\mathbb R^{R\times R}$ such that
$$
E[\beta_{irc}\mid X_i,N_i]=e_r^\top H_cX_i=\sum_{s=1}^R H_{c,rs}X_{is}\qquad\text{a.s. on }\{X_i\in A\}. \tag{1.4}
$$
Rows of $\beta$ sum to 1, so $\sum_c e_r^\top H_c x=1=\mathbf 1^\top x$ on an open subset of the simplex, which forces
$$\textstyle\sum_{c=1}^C H_c=\mathbf 1\mathbf 1^\top=:J. \tag{1.5}$$
(2×2 dictionary, with row 1 = black and $X_1=X$: $w_0=H_{22}$, $w_1=H_{21}-H_{22}$, $b_0=H_{12}$, $b_1=H_{11}-H_{12}$.)

**Observable implication.** Taking $E[\cdot\mid X]$ of (1.1),
$$
\tau_c(x):=E[T_c\mid X=x]=x^\top H_cx=x^\top Q_cx,\qquad Q_c:=\tfrac12(H_c+H_c^\top). \tag{1.6}
$$
Write $H_c=Q_c+K_c$ with $K_c:=\tfrac12(H_c-H_c^\top)$ skew-symmetric. A skew part never shows up in a quadratic form, since $x^\top Kx=0$.

**Lemma 1 (identification decomposition) (to prove).** Assume A5 below (the support of $X$ is not contained in a quadric). Then:
(a) $Q=(Q_1,\dots,Q_C)$ is point-identified from $\tau(\cdot)$;
(b) $K=(K_1,\dots,K_C)$ is *completely* unidentified by $\tau$ and is constrained only by $\sum_cK_c=0$ (from 1.5 together with $\sum_cQ_c=J$);
(c) so the unidentified parameter lives in the linear space
$$\mathcal K=\Big\{(K_1,\dots,K_C):K_c^\top=-K_c,\ \textstyle\sum_cK_c=0\Big\},\qquad \dim\mathcal K=(C-1)\,\tfrac{R(R-1)}{2}. \tag{1.7}$$

Check: $2\times2\Rightarrow1$ (JKST's $w_1$, with $K_{1,12}=c_1/2-w_1$). **$2\times3\Rightarrow 2$.** $3\times3\Rightarrow 6$. Interpretation: for each pair of rows $(r,s)$ and each column $c$, the data identify the *sum* of "the effect of $s$-share on row $r$" and "the effect of $r$-share on row $s$", but not how that sum splits between them. This is the R×C version of "is turnout low because blacks vote less, or because whites in black areas vote less?" (JKST Fig. 1).

**Estimand as a function of the parameters.** Under A3 (mean independence from $N$),
$$
B^*_{rc}(K)=\frac{1}{\mu_r}\Big[\sum_{s}Q_{c,rs}\mu_{rs}+\sum_{s\neq r}K_{c,rs}\mu_{rs}\Big],\qquad \mu_{rs}:=E[NX_rX_s],\ \mu_r:=E[NX_r]. \tag{1.8}
$$
This is *affine in $K$* with sensitivities $\partial B^*_{rc}/\partial K_{c,rs}=\mu_{rs}/\mu_r$ (for 2×2 this is JKST's $r$, eq. 18).

### 1.4 Constraint sets: the outer set (generalised JKST Prop. 2/3)

Let $A\subseteq\Delta^{R-1}$ be a compact convex polytope on which (1.4) is assumed, with vertex set $\operatorname{vert}(A)$.

$$
\boxed{\;\Theta_O(P;A)=\Big\{K\in\mathcal K:\ (Q_c(P)+K_c)\,v\ \ge 0\ \text{(elementwise)}\quad\forall c,\ \forall v\in\operatorname{vert}(A)\Big\}\;} \tag{1.9}
$$

**Proposition 2-RC (to prove).** (i) $\Theta_O$ is a convex polytope containing the true $K$ under M0. (ii) For $R=C=2$ and $A=[l,u]$ it equals JKST Prop. 2; for $A=\Delta^{R-1}$, where $\operatorname{vert}(A)=\{e_s\}$ and the constraint becomes $H_{c,rs}\ge0$, it equals JKST Prop. 3. (iii) It is exactly the set of $K$ for which every implied conditional-mean row $E[\beta_{r\cdot}\mid X=x]$ lies in $\Delta^{C-1}$ for all $x\in A$.

*Why it has this form:* the JKST Jensen step $E[U\mid X]\le\min\{1,E[T\mid X]/(1-X)\}$ generalises to "the conditional mean of a selection of the DD polytope lies in the DD polytope evaluated at $E[T\mid X]$". Because $x^\top H_c x=\tau_c(x)$ holds by construction, the only constraint left is nonnegativity, and a linear function on a polytope is minimised at a vertex. (Remark 1.4: this is the form checked numerically against the package.)

**Bounds on any linear functional** are linear programs:
$$
[\varphi_L,\varphi_U]=\Big[\min_{K\in\Theta_O}\varphi(K),\ \max_{K\in\Theta_O}\varphi(K)\Big]. \tag{1.10}
$$
**Closed form for $R=2$ (any $C$) (to prove).** Each constraint involves a single $\kappa_c:=K_{c,12}$, so $\Theta_O=\big(\prod_c[\ell_c,h_c]\big)\cap\{\sum_c\kappa_c=0\}$. Its projection is
$$\kappa_c\in\Big[\max\{\ell_c,-\textstyle\sum_{c'\neq c}h_{c'}\},\ \min\{h_c,-\sum_{c'\neq c}\ell_{c'}\}\Big].$$
**Corollary (to prove):** the joint 2×C bound is weakly tighter than running $C$ separate collapsed 2×2 analyses ("column $c$ vs the rest"), and strictly tighter whenever the sum-to-zero coupling binds. So the 2×3 transition analysis gains from being done jointly.

### 1.5 The sharp identification region

**Primitive objects.** Let $(\Omega,\mathcal F,\mathbb Q)$ carry $(W,\beta)$. The observable law is $P=\mathbb Q\circ W^{-1}$. The conditioning σ-algebra is $\mathcal G=\sigma(X,N,Z)$ under A3; it is $\sigma(X,Z)$ if A3 is dropped and a covariate stands in for $N$. A **structure** is a pair $s=(H,\mathbb Q_{\beta\mid W})$. The model $\mathcal M(A)$ is the set of laws $\mathbb Q$ such that

* (S1) $\beta\in\mathfrak B$ a.s. and $T=\beta^\top X$ a.s. (A1);
* (S2) $E_{\mathbb Q}[\beta\mid\mathcal G]=M_H(X):=[e_r^\top H_cX]_{rc}$ a.s. on $\{X\in A\}$ (A2, A3).

**Definition (sharp identified set).**
$$
\Theta_I(P;A)=\big\{K\in\mathcal K:\ \exists\,\mathbb Q\in\mathcal M(A)\text{ with }\mathbb Q\circ W^{-1}=P\text{ and }H=Q(P)+K\big\}. \tag{1.11}
$$
$\Theta_I$ is a correspondence $P\mapsto 2^{\mathcal K}$: it contains exactly the parameters that are observationally equivalent to the data. For the estimand, $\mathcal B^*_I(P)=\{B^*(K):K\in\Theta_I(P)\}$, which is an interval for a scalar $\varphi$ because $\Theta_I$ is convex and $B^*$ is affine.

**Random-set characterisation.** Let $S(x,t)=\{b\in\mathfrak B: b^\top x=t\}$ be the unit-level DD polytope (a transportation polytope in cell-share coordinates $n_{rc}=x_rb_{rc}$), with support function $h_S(u;x,t)=\max_{b\in S(x,t)}\langle u,b\rangle$, an LP that is concave in $t$.

**Theorem 1 (to prove).** Assume A1–A3 and that $\mathcal F$ is rich enough to allow auxiliary randomisation. Then
$$
K\in\Theta_I(P;A)\iff \langle u,M_{Q+K}(x)\rangle\le E\big[h_S(u;X,T)\mid X=x,N=n,Z=z\big]\quad\forall u\in\mathbb S^{RC-1},\ P\text{-a.e. }(x,n,z),\ x\in A, \tag{1.12}
$$
equivalently $M_{Q+K}(x)\in\mathbb E[S(X,T)\mid\mathcal G]$, the conditional Aumann (selection) expectation. Tools: Hiai–Umegaki conditional selection expectation, the Kuratowski–Ryll-Nardzewski measurable selection theorem, and Beresteanu, Molchanov & Molinari (2011). For 2×2 conditioned on $X$ only, (1.12) reduces *exactly* to JKST Prop. 1, which confirms their "tightest bound" claim in that case.

**Nesting (to prove).**
$$
\Theta_I(P;A)\ \subseteq\ \Theta_X(P;A)\ \subseteq\ \Theta_O(P;A). \tag{1.13}
$$
$\Theta_X$ is (1.12) conditioning on $X$ only, i.e. JKST Prop. 1; the gap $\Theta_I\subsetneq\Theta_X$ appears when the model conditions on $N$ (or $Z$) and $E[h_S\mid X,N]$ varies with $N$. $\Theta_O$ is (1.9). The second inclusion is the Jensen step. **Equality condition (conjecture, to prove or refute):** $\Theta_X=\Theta_O$ if, for a.e. $x$, $h_S(u;x,\cdot)$ is affine on the conditional support of $T\mid X=x$ for every $u$ that is active at the optimum. In 2×2 this means $T/(1-X)$ does not straddle 1 (and $T/X$ does not straddle 0 or 1) within an $X$-slice.

**Support caveat.** (1.8) uses the linear model on all of $\operatorname{supp}(X)$, so $A\supseteq\operatorname{supp}(X)$ is required for the formula. If $A\subsetneq\operatorname{supp}(X)$, the sharp set for $B^*$ must use the model-free Aumann bound on $\operatorname{supp}(X)\setminus A$ (to formalise). The package sets $A=$ the range of the data, which is consistent with this.

### 1.6 Covariates and district fixed effects

Write the district-specific coefficient matrices as $H_c^{(d)}$.

| Model | Assumption on $E[\beta_{rc}\mid X,N,D=d]$ | Unidentified parameter | Identified set |
|---|---|---|---|
| **M0** pooled | $e_r^\top H_cX$ (no district) | $K\in\mathcal K$ | $\Theta_O(P;A)$ |
| **M1** district intercepts, common contextual effects | $e_r^\top H_cX+\zeta_{c,d,r}$, $\ \sum_c\zeta_{c,d}=0$, $\ \zeta_{c,1}=0$ | common $K\in\mathcal K$ | $\bigcap_d\Theta_O^{(d)}$ |
| **M1′** slopes depend on district covariates $V_d$ | $e_r^\top(H_c+\textstyle\sum_\ell V_{d\ell}G_{c\ell})X+\zeta_{c,d,r}$ | $K$ and the skew parts of $G$ | polytope; to derive |
| **M2** saturated | $e_r^\top H^{(d)}_cX$, all district-specific | $K^{(d)}$ for each $d$ | $\prod_d\Theta_O(P_d;A_d)$ |

Here $\Theta_O^{(d)}=\{K:(Q^{(d)}_c+\operatorname{skew}(\zeta_{c,d}\mathbf1^\top)+K_c)v\ge0\ \forall v\in\operatorname{vert}(A_d),\ \forall c\}$, with $A_d$ the convex hull of the support of $X$ in district $d$. JKST App. C is M1 with a continuous $Z$.

**Structural facts (to prove in Phase 2b):**

* **(D0) Without a model, dummies add nothing.** Under A1 alone, the sharp set for $B_{rc}$ is the DD bound whatever $Z$ or $D$ is, because DD is already computed unit by unit and aggregation is additive. *Every gain from district dummies is therefore model-borne*: it comes only from cross-district restrictions such as common $K$.
* **(D1) Intercept shifts are point-identified.** Under M1, $Q_c^{(d)}=Q_c+\tfrac12(\zeta_{c,d}\mathbf1^\top+\mathbf1\zeta_{c,d}^\top)$, so $\zeta_{c,d}=\operatorname{diag}(Q^{(d)}_c-Q^{(1)}_c)$ is identified. Additive covariates add constraints but no new unidentified dimensions: $\dim$ of the unidentified part is still $(C-1)R(R-1)/2$, the same as M0. This generalises App. C, where $\tilde w_2$ and $\tilde b_2$ are identified.
* **(D2) Overidentification.** M1 implies $(Q^{(d)}_c-Q^{(1)}_c)_{rs}=\tfrac12(\zeta_{c,d,r}+\zeta_{c,d,s})$ for $r\ne s$. That is $R(R-1)/2$ testable restrictions per $(c,d\neq1)$; in 2×2 it is "the $X^2$ coefficient is equal across districts". In addition, $\bigcap_d\Theta_O^{(d)}=\emptyset$ refutes M1, which generalises the JKST revert-to-DD rule.
* **(D3) Intersection gain.** $\Theta^{M1}=\bigcap_d\Theta_O^{(d)}\subseteq\Theta_O^{(d)}$ for every $d$. For $R=2$ the bound on $\kappa_c$ is $[\max_d\ell_{c,d},\ \min_d h_{c,d}]$, so a single district whose intercept puts some row's conditional mean close to 0 or 1 at an extreme of its own $X$-range pins down the common contextual effect for *all* districts.
* **(D4) When dummies matter.** Under M1, $E[\beta\mid X]=M_H(X)+E[\zeta_D\mid X]$. So M0 holds (with a shifted intercept) iff $E[\zeta_D\mid X]$ is affine in $X$, e.g. when $D\perp X$, and **M0 is misspecified when district effects are correlated with $X$ non-affinely**. This is an ecological version of omitted-variable bias. If $\zeta\equiv0$, then $\Theta^{M1}=\Theta_O(P;\operatorname{conv}\bigcup_dA_d)=\Theta^{M0}$: dummies cost precision but neither help nor hurt identification. Corollary: dummies (i) protect against between-district aggregation bias, and (ii) can tighten bounds through (D3). M0 and M1 are **non-nested**; M2 is weaker than both.
* **(D5) M2 sums intervals.** $\mathcal B^{*,M2}_I=\sum_d\omega_d[\,B^{*}_{d,L},B^{*}_{d,U}]$, a Minkowski sum. It is valid under district-specific contextual effects and can be wider or narrower than M0.
* **By-product:** M1 and M2 give **district-level bounds** $B_{rc}(\mathcal I_d)$ as well, e.g. transition rates per constituency.

---

## 2. Assumptions

Each entry gives the statement, its mathematical role in the bounds, what it means for voter transitions, the consequence of a violation, and whether it can be tested.

**A1. Closed population / accounting identity.** $T_i=\beta_i^\top X_i$ exactly, with $N_i$ the same people in both margins.
*Role:* defines $S(x,t)$, DD and $\tau(x)=x^\top Q x$; without it nothing is bounded.
*Transitions:* both elections are measured on the *same* electorate (the same register), and rows include previous abstainers and the ineligible if they can enter.
*Violation:* turnover (deaths, moves, new registrants) or boundary changes give $T_i=\beta_i^\top X_i+\varepsilon_i$. The DD sets shift, $Q$ absorbs $E[\varepsilon\mid X]$, and the bounds can exclude the truth. Phase 3 treats $\|\varepsilon_i\|_1\le\epsilon$, which widens the bounds by $O(\epsilon)$. Partial remedy: extra "entered" row and "exited" column categories.
*Testable:* only through implied infeasibility, e.g. $T_i$ outside the range compatible with $X_i$.

**A2. Linear contextual effects (LCE) on $A$:** equation (1.4).
*Role:* the *only* identifying restriction; it turns the infinite-dimensional $E[\beta\mid X]$ into $H$, and gives (1.6), (1.8) and (1.9).
*Transitions:* A-voters' loyalty in precinct $i$ is linear in the precinct's previous shares, e.g. "A-voters in A strongholds are more loyal" (social context, local organisation).
*Violation:* threshold or tipping effects (party machines mobilise only above 50%), or U-shapes (competitive precincts mobilise both sides). The quadratic regression is then misspecified, $\Theta_O$ is computed for the wrong $Q$, and the bound can miss the truth by $O(1)$ that does not shrink with $p$. There are *undetectable* violations: in the R×C analogue of JKST Example 6, quadratic contextual effects cancel in $\tau$.
*Testable:* partly. A cubic term in the $T\sim X$ regression rejects A2; emptiness of $\Theta_O$ or of $\Theta_O\cap DD$ rejects A2 (App. B generalised); exact cancellations cannot be detected.

**A3. Mean independence from size:** $E[\beta\mid X,N,Z]=E[\beta\mid X,Z]$.
*Role:* makes the finite-population term $\sum_iN_iX_{ir}e_{irc}/\sum_iN_iX_{ir}$ mean zero, so that (1.8) holds with weights $N X_r$.
*Transitions:* large (urban) precincts do not differ in loyalty from small ones *given* their partisan mix.
*Violation:* $B^*\neq$ the regression functional; the bias is $\operatorname{Cov}(N X_r,\beta_{rc}-M(X))/E[NX_r]$. Remedy: put $\log N$ in $Z$ (App. C).
*Testable:* partly, by adding $X\cdot\log N$ interactions to the $T$ regression.

**A4. Support and extrapolation:** $\operatorname{supp}(X)\subseteq A$, compact, and A2 holds on all of $A$.
*Role:* $\Theta_O$ uses the vertices of $A$. A larger $A$ gives tighter bounds (more constraints), but only if A2 is believed there.
*Transitions:* the linear form holds from the least to the most A-dominated precinct observed. For $R\ge3$, $A$ is the convex hull of $X$ in the simplex, and its vertices are the extreme precincts.
*Violation:* extrapolating beyond the data, or outlier precincts at the hull vertices, make the constraints bind spuriously, so bounds are too narrow or empty. The vertices are estimated (App. C, Remark 4\*).
*Testable:* sensitivity to trimming $A$.

**A5. Rank / richness:** $\operatorname{supp}(X)$ is not contained in any quadric of the simplex (for $R=2$: at least 3 distinct values of $X$).
*Role:* point-identifies $Q$ (Lemma 1a).
*Transitions:* fails, or nearly fails, when shares barely vary, e.g. stable turnout across precincts for $R=3$. This is JKST's "gender" problem.
*Violation:* $Q$ is weakly identified, the bound inherits huge variance, and $\Theta_O$ is very wide or unstable.
*Testable:* yes, through the condition number of the design matrix.

**A6. Sampling: independence and identical distribution.** JKST: iid across units. Here: independent across districts with arbitrary dependence within, plus iid within district under M1/M2.
*Role:* LLN and CLT for $\hat Q$ and $\hat\mu$, and the $S_1$ term.
*Transitions:* precincts in the same constituency share candidates, campaigns and shocks.
*Violation:* plain HC1 errors are too small, so the CIs under-cover; the bounds' probability limits are unaffected. Remedy: cluster-robust variance by district.
*Testable:* yes (intraclass correlation of residuals).

**A7. Regularity for inference** (JKST tie-breaking): a unique active constraint at the optimum, $\mu_{rs}>0$, and $\Theta_O$ not a singleton.
*Role:* the LP value function is differentiable at $Q$, so the delta method and $\mathrm{CI}_x$ are valid.
*Violation:* at ties the bound is only directionally differentiable; the naive bootstrap and delta method fail, and Phase 3 needs Fang–Santos type inference. This is not a bias in identification.
*Testable:* partly (estimated gaps between constraint values).

**A8 (M1). Homogeneous contextual effects across districts:** $K$ (and the symmetric slope parts) are common to all $d$; only intercepts vary.
*Role:* this restriction is what produces the intersection $\bigcap_d\Theta_O^{(d)}$, which is the entire gain (D0, D3).
*Transitions:* how a precinct's previous A-share shifts A-voters' loyalty is the same in every constituency, while constituency-level popularity shifts (candidate quality, local swing) are free.
*Violation:* e.g. a strong incumbent changes the contextual gradient in one district. The intersection can exclude the truth, often with *falsely narrow* bounds, because the max/min across districts picks the misspecified district.
*Testable:* yes, by the overidentification (D2) and emptiness.

**A9 (M1/M2). Within-district richness:** at least $R$ distinct $X$ values per district (M1), or a quadric-free support per district (M2).
*Role:* identifies $\zeta_{c,d}$ (M1) or $Q^{(d)}$ (M2).
*Transitions:* districts with only a few, similar precincts identify almost nothing.
*Violation:* weak identification per district. Combined with $\max_d$/$\min_d$ this biases the bounds *inward*.
*Testable:* yes.

**A10. Asymptotic regime (a choice, not a fact):** (R1) $D$ fixed, $p_d\to\infty$; (R2) $D\to\infty$, $p_d$ bounded; (R3) both grow.
*Role:* under R2 the $\zeta_d$ are incidental parameters that cannot be estimated consistently, and $\max_d\hat\ell_{c,d}$ is biased upward. This needs intersection-bounds methods with precision correction (Chernozhukov, Lee & Rosen 2013) or a hierarchical restriction on $\zeta$.
*Transitions:* the answer depends on your design (see the questions at the end).

**Assumptions explicitly not made:** no distributional form for $\beta$ (unlike King's truncated bivariate normal), no independence between $\beta$ and $X$ (contextual effects are the point), no stationarity over time, and no point-identifying nonlinearity.

**Optional shape restrictions** (each is a convex constraint on $K$, so the bounds stay LPs):
* **Sign or monotone contextual effects**, e.g. "A-voters are no less loyal where A was stronger": $\partial_{x_A}E[\beta_{A,A}\mid x]\ge0$, a half-space in $K$. The bounds shrink; if the sign is wrong, the truth can be excluded.
* **Bounded heterogeneity** $|K^{(d)}-K|\le\delta$ (between M1 and M2): a sensitivity parameter that interpolates between the two (Phase 3).
* **Bhatia–Davis variance bound** $\operatorname{Var}(\beta)\le m(1-m)\le\tfrac14$: tightens $S_1$ (Phase 3); valid without extra assumptions.

---

## 3. Notation dictionary

Also available as a LaTeX `tabular` in `docs/extension/notation.tex`.

| Symbol | Meaning |
|---|---|
| $i=1,\dots,p$ | unit (precinct) index |
| $d=1,\dots,D$; $D_i$; $\mathcal I_d$; $p_d$ | district index; unit $i$'s district; units in $d$; their number |
| $r,s=1,\dots,R$ | row categories (election-1 choice / group) |
| $c=1,\dots,C$ | column categories (election-2 choice) |
| $N_i$, $N_{ir}$, $N_{irc}$ | unit size; row count; latent cell count |
| $X_i\in\Delta^{R-1}$ | observed row shares $N_{ir}/N_i$ |
| $T_i\in\Delta^{C-1}$ | observed column shares |
| $Z_i$ | additional covariates (e.g. $\log N_i$) |
| $\beta_i=[\beta_{irc}]\in\mathfrak B=(\Delta^{C-1})^R$ | latent transition matrix of unit $i$ |
| $W_i=(N_i,X_i,T_i,Z_i,D_i)$ | observables of unit $i$ |
| $P$ | law of $W$ (observable); $\mathbb Q$: law of $(W,\beta)$ |
| $(\Omega,\mathcal F)$ | underlying measurable space |
| $\mathcal G=\sigma(X,N,Z)$ (or $\sigma(X,N,Z,D)$) | conditioning σ-algebra of the contextual model |
| $\mathcal O_p=\sigma(W_1,\dots,W_p)$ | σ-algebra of the sample |
| $B_{rc}(\mathcal S)$ | finite-population transition rate over units $\mathcal S$ (1.2) |
| $B^*_{rc}$ | superpopulation analogue (1.3) |
| $\varphi=\sum a_{rc}B_{rc}$ | generic linear target (loyalty $\lambda_r$, defection $\delta_r$, demobilisation $\mu_r$) |
| $A\subseteq\Delta^{R-1}$, $\operatorname{vert}(A)$ | region where LCE is assumed; its vertices ($R=2$: $\{l,u\}$) |
| $H_c\in\mathbb R^{R\times R}$ | homogeneous LCE coefficients, column $c$ (1.4) |
| $Q_c=\operatorname{sym}(H_c)$ | identified part, from $E[T_c\mid X]=X^\top Q_cX$ |
| $K_c=\operatorname{skew}(H_c)$; $\mathcal K$ | unidentified part; its space, $\dim=(C-1)R(R-1)/2$ |
| $\kappa_c=K_{c,12}$ | the unidentified scalar per column when $R=2$ |
| $M_H(x)=[e_r^\top H_cx]_{rc}$ | implied conditional-mean transition matrix |
| $\tau_c(x)=E[T_c\mid X=x]$ | regression function (quadratic under A2) |
| $\mu_r=E[NX_r]$, $\mu_{rs}=E[NX_rX_s]$ | weighting moments; $\mu_{rs}/\mu_r$ = sensitivity of $B^*_{rc}$ to $K_{c,rs}$ |
| $S(x,t)$; $h_S(u;x,t)$ | unit DD polytope; its support function |
| $\mathbb E[S\mid\mathcal G]$ | conditional Aumann (selection) expectation |
| $\Theta_I,\ \Theta_X,\ \Theta_O$ | sharp set; JKST-Prop.-1 set; outer (Jensen / vertex) set |
| $\mathcal B^*_I$ | identified set for the estimand, $\{B^*(K):K\in\Theta_I\}$ |
| $\zeta_{c,d}\in\mathbb R^R$ | district-$d$ intercept shifts (M1), identified |
| $Q^{(d)}_c$, $\Theta^{(d)}_O$, $A_d$ | district-level identified part, constraint set, support |
| $\mathrm{CI}_x$, $S_1$, $S_L$, $S_U$ | JKST confidence interval and its standard-error terms |
| $\Phi$ | standard normal CDF |

---

## 4. Proposed refinements of the phase plan

The three-phase structure works, but goal 2 and the code need their own steps. I propose:

* **Phase 2a: R×C derivations.** Lemma 1; Prop. 2-RC with the vertex reduction and 2×2/JKST equivalence; Theorem 1 (sharpness of $\Theta_I$ via selection expectations, with both directions of observational equivalence); the nesting (1.13) and the equality condition; the closed form for $R=2$ and the joint-vs-collapsed corollary; non-emptiness (generalised App. B); the R×C version of JKST Example 6 (undetectable violation).
* **Phase 2b: covariates and districts.** D0–D5 with proofs; sharpness of $\bigcap_d\Theta^{(d)}$ under M1; a characterisation of *when the gain is strict* in terms of $\zeta_d$ and $A_d$; the misspecification bias of M0 under M1, in closed form; M1′.
* **Phase 3: inference and sensitivity.** Consistency at rate $\sqrt p$ (or $\sqrt{D}$ under clustering) of LP-value estimators; asymptotic distribution at a unique active vertex (normal) and at ties (max/min of Gaussians, directional delta method, Fang & Santos 2019); generalised $\mathrm{CI}_x$ (cover-the-parameter, Imbens–Manski/Stoye) vs cover-the-set; the many-districts problem (CLR 2013 intersection bounds); breakdown-frontier sensitivity for A1 ($\epsilon$), A2 (sup-norm deviation $\delta$ from linearity) and A8 (slope heterogeneity).
* **Phase 4 (new): R implementation.** New exported functions (suggested names `eiBoundsRxC()` and `eiBoundsDistrict()`) that leave the existing API untouched. Closed forms for $R=2$; a small LP otherwise. The package's dependencies are `sandwich`, `MASS` and `eco`; general $R$ needs an LP solver (`lpSolve` is light), which is a dependency decision for you. Cluster-robust `vcovCL`. Regression test: for $R=C=2$ the new code must reproduce `generateBounds` and `evaluateBounds` to $10^{-10}$.
* **Phase 5 (new): validation.** Monte Carlo (correctly specified, $K$ at a vertex for tightness of coverage, and each assumption violated in turn); ground-truth data (`eco::census` for 2×2; `eiPack::senc`, NC party registration by race with known cells, for a true 3×3 check); width ratio and coverage tables as in JKST Table 1.

## 5. Keeping the work rigorous

1. **Label every claim**: *Theorem / Proposition* (proved), *Conjecture*, or *Heuristic* (e.g. the revert-to-DD rule). Never call an outer set "sharp".
2. **State the measure-theoretic level explicitly**: conditional statements a.s.; the existence of measurable selections cited, not assumed; auxiliary randomisation stated.
3. **Cite exact theorems**, e.g. Artstein (1983); Hiai & Umegaki (1977); Molchanov (2017, *Theory of Random Sets*, 2nd ed.); Beresteanu, Molchanov & Molinari (2011); Molinari (2020, *Handbook of Econometrics* 7A); Cross & Manski (2002), the canonical partial-identification treatment of EI with covariates; Chernozhukov, Lee & Rosen (2013); Imbens & Manski (2004); Stoye (2009); Fang & Santos (2019). Use Rosen, Jiang, King & Tanner (2001) and Greiner & Quinn (2009) for the R×C EI context. I will check each citation against the source text before relying on it.
4. **Check algebra mechanically**: verify each reparametrisation symbolically (sympy) and numerically, and solve every LP closed form two ways (vertex enumeration and a generic solver), as already done for the 2×2 equivalence.
5. **Nested-model sanity checks**: every R×C or district result must collapse to JKST when $R=C=2$, $D=1$.
6. **Adversarial simulation for every assumption**: for each assumption, a DGP that violates it, to confirm the predicted failure mode (direction and order of the error).
7. **Separate identification from inference**: Phase 2 statements concern population objects ($P$ known). Sampling enters only in Phase 3.
8. **Justification tags**, as you require, will reference specific results ("[Lemma 1(a)]", "[Hiai–Umegaki 1977, Thm 5.1]") rather than generic labels wherever possible.

## 6. Questions that change the derivations

1. **Rows of the transition table.** Is the electorate restricted to election-1 voters of A and B (true 2×3, $\dim\mathcal K=2$), or do previous abstainers and new entrants form rows (3×3 or more, $\dim\mathcal K=6$)? Are both margins on the same register (A1)?
2. **District design.** Roughly how many districts $D$, and how many units per district $p_d$? Do units within a district vary in $X$? This decides between regimes R1 and R2 (A10), and whether Phase 3 needs intersection-bounds corrections.
3. **Targets.** Overall rates only, district-level rates $B_{rc}(\mathcal I_d)$ too, or both?
4. **Priority.** I recommend implementing $\Theta_O$ first (closed form / LP, extends the package directly), then a nonparametric estimate of the sharp $\Theta_I$ as an option. Do you agree?
5. **Dependencies.** May the package take on `lpSolve`, needed only for $R\ge3$?

---

## 7. Decisions after first review (27–28 Oct 2026)

* **General R×C** is the target. The running example is BC transitions, **3×3** (in-party voters, out-party voters, non-voters in each election), so $\dim\mathcal K=(3-1)\cdot3=6$; $X_i$ lives in the 2-simplex, and $A$ is the convex hull of the observed polls.
* **Target:** overall (province-wide) rates $B_{rc}$ only. District-level rates are a by-product, not a goal.
* **Data design.**
  * *Stage 1:* 93 ridings as the units, so $p=93$ and $D=1$. Riding dummies are impossible at this stage (one unit per riding). Only M0 applies, and with $p=93$ and six quadratic coefficients per column, A5 and A7 (precision) are the binding concerns.
  * *Stage 2:* about 10 polls per riding, so $p\approx930$, $D=93$, $p_d\approx10$. This is regime **R2** (many districts, few units each). The incidental-parameter problem of A10 is first-order here, not a footnote (see the illustration below). Phase 3 must deliver a precision-corrected intersection (Chernozhukov–Lee–Rosen) and/or partial pooling of the riding intercepts.
* **Priority:** $\Theta_O$ first, as an extension of the package; the nonparametric $\Theta_I$ comes later as an option.
* **Dependencies:** `lpSolve` is approved.

### 7.1 Illustration: what riding dummies do (`illustrate_district_dummies.py`)

The script uses a 2×2 setup for transparency: 93 ridings × 10 polls, riding baselines correlated with the riding's mean $X$, and contextual slopes common to all ridings (M1 true). Results over 300 replications:

| Bound | Covers truth | Median width |
|---|---|---|
| DD (model-free) | 1.00 | about 0.54 |
| Pooled M0 (current package) | 1.00 | 0.305 |
| Riding M1, **true** coefficients (identification only) | 1.00 | 0.244 |
| Riding M1, **estimated** from 10 polls per riding | **0.72** | 0.172 |

* **Identification gain (row 3 vs row 2):** about a 20% narrower width, because each riding constrains the common contextual effect from its own baseline.
* **Estimation problem (row 4):** taking the max/min over 93 noisy riding-level constraints makes the bounds *too narrow* and they miss the truth 28% of the time. The plain plug-in intersection must not be used at 10 polls per riding; this is what Phase 3 has to fix.
* With step-shaped riding effects (pooled model misspecified), the pooled bound still covered the truth here (0.997). The robustness argument (D4) is a theoretical guarantee, but in these designs its practical value was small compared with the tightening.

### 7.2 What $H$ and $K$ mean, and what M1 shares (refines A8)

* $E[\beta_{rc}\mid x]=\sum_sH_{c,rs}x_s$, so $H_{c,rs}$ is the rate at which row-$r$ voters go to column $c$ **in a (hypothetical) poll made up entirely of row-$s$ voters**. Diagonal entries are "own-territory" rates. Off-diagonal entries are rates of $r$-voters in $s$-territory, and the contextual effect for row $r$ is how $H_{c,r\cdot}$ varies across $s$.
* $Q_{c,rs}=\tfrac12(H_{c,rs}+H_{c,sr})$ is identified. $K_{c,rs}=\tfrac12(H_{c,rs}-H_{c,sr})$, half the difference between "$r$-voters in $s$-territory" and "$s$-voters in $r$-territory", is not. $K$ is therefore the **unidentified part of the contextual effect, not the whole of it**.
* **M1 shares all of $H$ except an additive riding shift:** $E[\beta_{rc}\mid x,d]=\zeta_{c,d,r}+\sum_sH_{c,rs}x_s$. Equivalently, there is no interaction between riding and poll composition on the probability scale. The shared symmetric part is testable (D2); the shared $K$ is not, and it is the substantive leap.
* Additivity on the probability scale is least plausible near 0 or 1 (ceiling and floor compression). A weaker option for Phase 2b is to let $Q^{(d)}$ be free and share only $K$, or to allow $\|K^{(d)}-K\|\le\delta$ as a sensitivity parameter. The free-$Q^{(d)}$ version needs 6 coefficients per column per riding in 3×3, which 10 polls cannot support.

### 7.3 Worked 2×2 example: M0 vs M1 (`example_2x2_ridings.py`)

Rows: $P$ = voted P last time (share $x$), $O$ = everyone else. Outcome: votes P now. *Loyalty* is the rate for P-voters; *crossover* is the rate for O-voters. The example is noise-free, so it shows identification only.

| | M0: one model for the province | M1: riding dummies |
|---|---|---|
| Loyalty in a poll | $a_P+g_Px$ | $a_{P,d}+g_Px$ |
| Crossover in a poll | $a_O+g_Ox$ | $a_{O,d}+g_Ox$ |
| Same everywhere | all four numbers | the slopes $g_P,g_O$ only |
| Riding matters | only through its polls' $x$ | shifts each group's baseline |
| Data reveal | $T=a_O+(a_P-a_O+g_O)x+(g_P-g_O)x^2$, pooled | the same, riding by riding, with common $x^2$ term |
| Unknown | $g_O$ (one number) | $g_O$ (one number, shared) |
| Testable | cubic term, empty bound | plus: $x^2$ term equal across ridings |
| Not testable | how the curvature splits into $g_P,g_O$ | that the split is the same in every riding |

In words, M1 says that **the gap between two ridings' loyalty is the same in every kind of poll** (parallel lines), and likewise for crossover.

Truth: $g_P=g_O=0.10$. Riding 1 (P stronghold, polls $x\in[0.5,0.9]$): loyalty $0.89+0.1x$, crossover $0.20+0.1x$. Riding 2 (P weak, $x\in[0.1,0.5]$): loyalty $0.45+0.1x$, crossover $0.02+0.1x$.

* Riding 1's polls reveal $T=0.20+0.79x$. Loyalty at $x=0.5$ cannot exceed 1, so $g\ge-0.02$. Riding 2's polls reveal $T=0.02+0.53x$. Loyalty at $x=0.1$ cannot fall below 0, so $g\le0.611$.
* **M1:** $g\in[-0.020,0.611]$, with riding 1 setting the lower limit and riding 2 the upper. This gives $B\in[0.618,0.866]$ (truth $0.819$; DD $[0.453,0.879]$).
* **M0** fits one curve through both ridings: $T=-0.095+0.936x+0.305x^2$. It reads the stronghold's higher baselines as a contextual effect, and the implied bound on the slope is empty ($[1.160,0.969]$). The package's rule then reverts to DD. Here M0's failure was *detectable*; it need not always be.
