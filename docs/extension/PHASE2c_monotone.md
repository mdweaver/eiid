# Phase 2c: monotone ecological inference in R×C, derivations and proofs

*Builds on P1, 2a, 2b, P3 and the Phase-1 proposal `PHASE1b_monotone_did.md` (cited as 1b). It follows Elzayn, Goldin, Guage, Ho & Morton, "Monotone Ecological Inference" (EGGHM; arXiv v1 numbering, see 1b §0). Every proof step ends with a justification tag: those of 2a §0.3, 2b §0.3 and P3 §0.3, plus §0.3 below. Checks: `verify_phase2c.py` (output `verify_phase2c_output.txt`, 12/12 PASS). Owner decisions applied (1b §D):*
* *exact (finite-population) bounds for the sign rungs;*
* *one-vs-rest (W1) as the default, with W2 as an option;*
* *contextual reinforcement available as an option;*
* *API `monotone=` plus `model=`.*

---

## 0. Setting

### 0.1 Finite population (the exact rungs)

* The units $i=1,\dots,p$ are fixed (the design is conditioned on, as in P3 (P3-S2)). Unit weights are $p_i=N_i/\sum_jN_j$.
* $X_i\in\Delta^{R-1}$, $T_i\in\Delta^{C-1}$, and the realised rates satisfy $\beta_i\in S(X_i,T_i)$ (2a §0.2).
* $E_p[f]:=\sum_ip_if_i$, and $\mu_r:=E_p[X_r]>0$.
* Row weights $w^r_i:=p_iX_{ir}/\mu_r$. $\operatorname{Cov}_w$ is the covariance under weights $w$.
* **Target** (P3's actual rate): $B_{rc}=E_p[X_r\beta_{rc}]/\mu_r=\sum_iw^r_i\beta_{irc}$.

### 0.2 Superpopulation (the K-rows and N2 rungs)

* $m(v)=E[\beta\mid V=v]$, and P3's model-implied rate is $\tilde B_{rc}=\sum_iw^r_im_{rc}(V_i)$.
* A superscript $^*$ on an association means it is computed with $m$ in place of $\beta$.

### 0.3 Additional tags

| Tag | Meaning |
|---|---|
| [EGGHM:x] | result x of EGGHM, arXiv v1 (2504.14752). Their identities hold for any bounded outcome; checked again in 1b V1 |
| [Cheb] | Chebyshev's association inequality: for nondecreasing $f,g:\mathbb R\to\mathbb R$ and any probability measure on $\mathbb R$, $\operatorname{Cov}(f(Z),g(Z))\ge0$ (e.g. Boucheron, Lugosi & Massart 2013, Thm 2.14†) |
| [ext] | a function that is monotone on a finite subset of $\mathbb R$ extends to a monotone function on $\mathbb R$ (piecewise-constant extension) |
| [conn] | a continuous image of a connected set is connected; a connected subset of $\mathbb R$ is an interval |
| [CE-cvx], [CE-lin], [CE-pull] | as in 2a §0.3 |

---

## 1. The associations in R×C and their identities

Fix a column $c$ and drop it from the notation.

**Definition 1.**
1. **Between, pairwise:** $\delta_B(r,s):=E_p[X_rX_s(\beta_r-\beta_s)]$ for $r\ne s$. It is antisymmetric.
2. **Between, one-vs-rest:** $\delta^r_B:=\sum_{s\ne r}\delta_B(r,s)$.
3. **Within, W2 (each group separately):** $\Delta:=\sum_s\mu_s\operatorname{Cov}_{w^s}(X,\beta_s)\in\mathbb R^R$.
4. **Within, W1 (one-vs-rest; EGGHM §4 end).** Let $\beta^{\neg r}_i:=\sum_{s\ne r}X_{is}\beta_{is}/(1-X_{ir})$, the rate of everyone outside row $r$ (any value if $X_{ir}=1$, where its weight is 0), and $w^{\neg r}_i:=p_i(1-X_{ir})/(1-\mu_r)$. Then
$$\delta^r_W:=\mu_r\operatorname{Cov}_{w^r}(\beta_r,X_r)+(1-\mu_r)\operatorname{Cov}_{w^{\neg r}}(\beta^{\neg r},X_r).$$

**Estimators** (functions of the data only):
* Goodman's no-intercept regression: $g:=M^{-1}E_p[XT_c]$, with $M:=E_p[XX^\top]$.
* The neighbourhood model: $NM_r:=E_p[X_rT_c]/\mu_r$.
* The one-vs-rest regression with intercept: $ER^r:=\operatorname{Cov}_p(T_c,X_r)/\operatorname{Var}_p(X_r)$.
* The one-vs-rest gap: $D^r:=B_r-B_{\neg r}$, with $B_{\neg r}:=\sum_{s\ne r}\mu_sB_s/(1-\mu_r)$.

**Lemma M1 (identities; accounting identity only).** Assume [A1]. Assume $M$ is invertible, i.e. the $X_i$ span $\mathbb R^R$, the linear version of A5.
(a) $\mathbf 1^\top\Delta=0$ and $g-B_{\cdot c}=M^{-1}\Delta$.
(b) $NM_r-B_r=-\delta^r_B/\mu_r$.
(c) The collapsed population $(X_r,\beta_r,\beta^{\neg r})$ is a valid two-group population for column $c$. Hence
$$ER^r-D^r=\frac{\delta^r_W}{\operatorname{Var}_p(X_r)},\qquad (NM_r-NM_{\neg r})-D^r=-\frac{\delta^r_B}{\mu_r(1-\mu_r)},$$
where $NM_{\neg r}=E_p[(1-X_r)T_c]/(1-\mu_r)$.
(d) $B$, $\delta_B(r,s)$, $\delta^r_B$, $\Delta$ and $\delta^r_W$ are affine functions of $\beta=(\beta_i)_i$, with coefficients that depend only on the design.

*Proof of (a).*
1. $T_{ic}=\sum_sX_{is}\beta_{is}$. [A1]
2. $E_p[XT_c]=\sum_sE_p[XX_s\beta_s]$. [1; lin]
3. $E_p[XX_s\beta_s]=\mu_sE_{w^s}[X\beta_s]=\mu_s\big(E_{w^s}[X]B_s+\operatorname{Cov}_{w^s}(X,\beta_s)\big)$. [Def $w^s$; Def Cov; $B_s=E_{w^s}\beta_s$]
4. $\mu_sE_{w^s}[X]=E_p[XX_s]=Me_s$. [Def]
5. So $E_p[XT_c]=MB_{\cdot c}+\Delta$, and $g=B_{\cdot c}+M^{-1}\Delta$. [2–4; $M$ invertible]
6. $\mathbf 1^\top\Delta=\sum_s\mu_s\operatorname{Cov}_{w^s}(\mathbf 1^\top X,\beta_s)=0$, since $\mathbf 1^\top X=1$ is constant. [lin; simplex] ∎(a)

*Proof of (b).*
1. $\mu_rNM_r=E_p[X_rT_c]=E_p[X_r^2\beta_r]+\sum_{s\ne r}E_p[X_rX_s\beta_s]$. [A1; lin]
2. $\mu_rB_r=E_p[X_r\beta_r]=E_p[X_r^2\beta_r]+\sum_{s\ne r}E_p[X_rX_s\beta_r]$, using $X_r=1-\sum_{s\ne r}X_s$. [Def; simplex]
3. Subtract: $\mu_r(NM_r-B_r)=\sum_{s\ne r}E_p[X_rX_s(\beta_s-\beta_r)]=-\delta^r_B$. [1–2; Def 1] ∎(b)

*Proof of (c).*
1. $\beta^{\neg r}_i$ is a convex combination of $\beta_{is}\in[0,1]$ with weights $X_{is}/(1-X_{ir})$, so it lies in $[0,1]$. [Def; simplex]
2. $T_{ic}=X_{ir}\beta_{ir}+(1-X_{ir})\beta^{\neg r}_i$. [A1; Def]
3. $E_p[(1-X_r)\beta^{\neg r}]/(1-\mu_r)=\sum_{s\ne r}\mu_sB_s/(1-\mu_r)=B_{\neg r}$. [Def; lin]
4. So $(X_r,\beta_r,\beta^{\neg r})$ is a two-group population in EGGHM's sense, with group means $B_r$ and $B_{\neg r}$. [1–3]
5. Its between association is $E_p[X_r(1-X_r)(\beta_r-\beta^{\neg r})]=E_p[X_r((1-X_r)\beta_r-\sum_{s\ne r}X_s\beta_s)]=\delta^r_B$. Its within association is $\delta^r_W$ by definition. [Def; simplex]
6. Apply [EGGHM:P4] and [EGGHM:P6] to it. Here the individual-level variance is $\operatorname{Var}(X)=\mu_r(1-\mu_r)$ and the between-unit variance is $\operatorname{Var}_p(X_r)$. [4–5] ∎(c)

*Proof of (d).* $B$ and $\delta_B$ are linear in $\beta$ with design weights. [Def] $\operatorname{Cov}_w(\beta_s,X_r)=\sum_iw_i\beta_{is}(X_{ir}-E_wX_r)$ is linear in $\beta$. [Def Cov] In $\delta^r_W$, $(1-\mu_r)w^{\neg r}_i\beta^{\neg r}_i=p_i\sum_{s\ne r}X_{is}\beta_{is}$, so the division by $1-X_{ir}$ cancels. [Def] ∎

Checked: C1 (R = C = 3, 100 populations, max error $6\times10^{-15}$).

---

## 2. Exact bounds under sign restrictions (rung N1)

**Definition 2.**
* A *sign restriction* is $g_k(\beta)\ge-\varepsilon_k$ with $g_k=\pm$ one of the associations of Definition 1, for any column, and slack $\varepsilon_k\ge0$ (default 0).
* Optional *unit-level* restrictions $\beta_i\in U_i$, with $U_i$ polyhedral: unit-level ordering, structural zeros (MENU A1, B5).
* The **N1 feasible set** is
$$F:=\{\beta:\ \beta_i\in S(X_i,T_i)\cap U_i\ \forall i;\ g_k(\beta)\ge-\varepsilon_k\ \forall k\}.$$
* **Contextual reinforcement** for a chosen pair $(g_B,g_W)$ is $F_{CR}:=F_+\cup F_-$, where $F_\pm:=F\cap\{\pm g_B\ge0,\pm g_W\ge0\}$.

**Proposition M2 (N1 is exact).** Let $\varphi=\sum a_{rc}B_{rc}$.
(a) The set of values of $\varphi$ attained by unit tables consistent with the margins and the restrictions is $\varphi(F)$. If $F\neq\emptyset$ it equals $[\min_F\varphi,\max_F\varphi]$: **sharp in the finite-population sense**.
(b) $F=\emptyset$ iff the restrictions contradict the data.
(c) Under CR the set is $\varphi(F_+)\cup\varphi(F_-)$, a union of at most two intervals.
(d) $\varepsilon\mapsto\varphi(F(\varepsilon))$ is nested, with continuous piecewise-linear endpoints where non-empty. So the breakdown value of a conclusion is well defined.

*Proof.*
1. A table $n_{irc}=N_iX_{ir}\beta_{irc}$ reproduces the margins iff $\beta_i\in S(X_i,T_i)$. [A1; Def $S$]
2. It satisfies the restrictions iff $\beta\in F$. [Def 2; L M1(d)]
3. $F$ is a bounded polyhedron, since each $S(X_i,T_i)$ is compact (2a L0(b)) and the restrictions are linear. $\varphi$ is affine in $\beta$. [L M1(d)]
4. So $\varphi(F)$ is a compact interval with both endpoints attained, which gives (a). (b) is the definition. [conn; 3]
5. (c): $F_+$ and $F_-$ are each of the form in (a). [Def 2; 4]
6. (d): $\varepsilon$ enters only the right-hand side. [LP-par; nest] ∎

**Proposition M3 (two rows: closed forms, refutability, CR).** Let $R=2$, write $X=X_1$, and fix a column. Write $D=B_1-B_2$, $\gamma=\operatorname{Var}_p(X)/(\mu(1-\mu))$, and let $[D^-,D^+]$ be the MOB interval (DD for $D$). Assume $\operatorname{Var}_p(X)>0$.
(a) $\delta_W\ge0\iff D\le D_{ER}$, and $\delta_B\ge0\iff D\ge D_{NM}$.
(b) $\{\delta_W\ge0\}$ gives $[D^-,\min\{D^+,D_{ER}\}]$. It is **empty iff $D_{ER}<D^-$**. Symmetrically, $\{\delta_W\le0\}$ is empty iff $D_{ER}>D^+$.
(c) $D_{NM}\in[D^-,D^+]$ always. So sign restrictions on $\delta_B$ **are never refuted by the data**: they give $[D_{NM},D^+]$ or $[D^-,D_{NM}]$.
(d) $D_{NM}=\gamma D_{ER}$ with $\gamma\in(0,1]$. Under CR, exactly one of $F_\pm$ is non-empty when $D_{ER}\neq0$; both give $\{0\}$ when $D_{ER}=0$. The identified set is
$$[D_{NM},\min\{D^+,D_{ER}\}]\ \text{if }D_{ER}\ge0,\qquad[\max\{D^-,D_{ER}\},D_{NM}]\ \text{if }D_{ER}<0.$$
**CR is never refuted by the data**, and under it $\operatorname{sign}D=\operatorname{sign}D_{ER}$.

*Proof.*
1. $\delta_W=(D_{ER}-D)\operatorname{Var}_p(X)$ and $\delta_B=(D-D_{NM})\mu(1-\mu)$. [L M1(c) with $R=2$, where $\neg1=2$] This gives (a). [$\operatorname{Var}_p(X)>0$]
2. $D$ is affine in $\beta$ over the product of the unit polytopes. Its range is $[D^-,D^+]$, with every value attained. [L M1(d); M2(a); EGGHM:P3]
3. (b): $F=\{\beta:D(\beta)\le D_{ER}\}$. Its image is $[D^-,D^+]\cap(-\infty,D_{ER}]$, which is empty iff $D_{ER}<D^-$. [1; 2; conn]
4. (c): put $\beta_{i1}=\beta_{i2}=T_i$ in every unit. This lies in $S(X_i,T_i)$, gives $\delta_B=0$, and hence $D=D_{NM}$. [Def $S$; 1]
5. (d): $D_{NM}=\gamma D_{ER}$ [EGGHM:L1]. $\gamma>0$ by $\operatorname{Var}_p(X)>0$. $\gamma\le1$ because the between-unit variance is at most the total variance $\mu(1-\mu)$ (law of total variance).
6. Suppose $D_{ER}>0$. Then $F_-$ requires $D\ge D_{ER}$ and $D\le D_{NM}=\gamma D_{ER}\le D_{ER}$. So $D=D_{ER}=D_{NM}$, which forces $\gamma=1$. Otherwise $F_-=\emptyset$. If $\gamma=1$, every unit has $X_i\in\{0,1\}$ and the problem is trivial. [5; 1]
7. $F_+=\{D_{NM}\le D\le D_{ER}\}\cap[D^-,D^+]=[D_{NM},\min\{D^+,D_{ER}\}]$. This is non-empty because $D^-\le D_{NM}\le D_{ER}$ and $D_{NM}\le D^+$. [1; 4; 5]
8. Then $D\ge D_{NM}>0$. The case $D_{ER}<0$ is symmetric, and $D_{ER}=0$ gives $D_{NM}=0=D$. [5–7] ∎

Checked: C2 (300 populations; LP = closed forms, $\max$ error $6\times10^{-15}$; emptiness iff $D_{ER}<D^-$; one CR branch; sign agreement in all 152 populations satisfying realised CR).

**Reading.** In two rows, only signs on the *within* association can be contradicted by the data. Between signs and CR can never be rejected from aggregates, so their credibility rests entirely on substantive argument (or auxiliary individual data, as in EGGHM §6).

**Proposition M4 (joint restrictions).** Imposing sign restrictions on several columns or rows jointly gives a subset of $F$ for any single restriction. So joint bounds are weakly tighter.
*Proof.* Adding constraints shrinks the feasible set. [Def 2; nest] ∎

Checked (C3): in 2×3, signing the within association of two columns jointly narrowed the bound for column 1 strictly in 75 of 139 random populations (max width gain 0.27). This is why the package solves one LP rather than applying EGGHM's 2×2 closed forms column by column.

---

## 3. Sign restrictions under linear contextual effects (K-rows)

**Proposition M5.** Under (S2) for any LCE rung (M0, M1, M1(δ), M2; 2a, 2b), with $m=M_H(X)$:
(a) $\delta^*_B(r,s)$, $\delta^{r*}_B$, $\Delta^*$ and $\delta^{r*}_W$ are affine in the free coordinates $k$ of $K$ (2a Prop. 7). So each sign restriction is one LP row added to $Gk\le g$, and CR is a union of two LPs.
(b) Mean-level ordering on $A$, $M_H(x)_{rc}\ge M_H(x)_{sc}$ for $x\in A\supseteq\operatorname{supp}$, implies $\delta^*_B(r,s)\ge0$.
(c) For $R=2$, $\delta^*_W=a(w_1+d_1)+b\,w_1$ with $a=\mu\operatorname{Var}_{w^1}(X)\ge0$ and $b=(1-\mu)\operatorname{Var}_{w^2}(X)\ge0$. Hence both contextual slopes $\ge0$ ⇒ $\delta^*_W\ge0$.

*Proof.*
1. $m(X_i)=M_{Q+K}(X_i)$ is affine in $k$. [2a L1; Def]
2. Each association is linear in $m$ with design coefficients. [L M1(d) with $m$ for $\beta$] So it is affine in $k$. [1] (a) follows. [Def 2]
3. (b): $\delta^*_B(r,s)=E_p[X_rX_s(m_r-m_s)]$, and the integrand is $\ge0$ at every $X_i\in A$. [Def; hypothesis]
4. (c): $\operatorname{Cov}_{w^1}(b_0+b_1X,X)=b_1\operatorname{Var}_{w^1}(X)$, and likewise for row 2. Then $b_1=w_1+d_1$ (P1 §1.3 dictionary). [Def Cov; alg] ∎

Checked: C4 (affinity in $k$ for R = C = 3; ordering ⇒ $\delta^*_B\ge0$) and 1b V3 (closed form in 2×2).

**Remark M5.1 (which δ).**
* N1 uses the realised associations of the actual population (owner decision 1).
* K-rows use $\delta^*$, the behavioural version.
* Under (S2), $\delta-\delta^*$ is a weighted sum of unit residuals, with mean zero given the design and order $p^{-1/2}$ (P3 §0.2).
* Intersecting N1 with an LCE rung therefore mixes the two readings. The package will report them separately.

---

## 4. Monotone nonparametric contextual effects (rung N2; two rows)

**Model N2(σ).**
* (S1) holds, and $E[\beta\mid V]=m(X)$ for a Borel $m:[0,1]\to\mathfrak B$ (A3). Here $X$ is the row-1 share.
* For designated pairs $(r,c)$ with directions $\sigma_{rc}\in\{+1,-1\}$, $\sigma_{rc}m_{rc}$ is nondecreasing on $\operatorname{supp}P_X$.
* No linearity is assumed. $\tau(x):=E[T\mid X=x]$ is identified.

**Lemma N2.1.** $m(x)\in S(x,\tau(x))$ for $P_X$-a.e. $x$.
*Proof.*
1. $m(X)=E[\beta\mid X]\in\mathfrak B$ a.s. [CE-tower with A3; CE-cvx]
2. $m(X)^\top X=E[\beta^\top X\mid X]=E[T\mid X]=\tau(X)$ a.s. [CE-pull; CE-lin; A1] ∎

**Outer set and LP.** Let $x_1<\dots<x_J$ be the support points (the distinct observed compositions). Then
$$\mathcal M_O:=\{(m(x_j))_j:\ m(x_j)\in S(x_j,\tau(x_j))\ \forall j;\ \sigma_{rc}\big(m_{rc}(x_{j+1})-m_{rc}(x_j)\big)\ge0\ \forall j,(r,c)\}.$$
Bounds on $\tilde B_{rc}=\sum_iw^r_im_{rc}(X_i)$ are the LP $[\min,\max]_{\mathcal M_O}$, which has $2CJ$ variables.

**Proposition N2.2.**
(a) Every $m$ allowed by N2(σ) restricts to a point of $\mathcal M_O$. So the LP bounds are an **outer** set for $\tilde B$. They are not claimed sharp: a sharp version would replace $S(x,\tau(x))$ by $\mathbb E[S(X,T)\mid X=x]$, as in 2a Thm 1.
(b) Conversely, every point of $\mathcal M_O$ extends to a monotone Borel $m$ on $[0,1]$ satisfying Lemma N2.1 on the support. So restricting to support points loses nothing.
(c) **Nesting with LCE.** If $\tau$ is the quadratic of 2a Lemma 1, then $\{M_{Q+K}:K\in\Theta_O\cap\{\sigma\text{-signed slopes}\}\}\subseteq\mathcal M_O$. Hence $\mathcal B(\text{M0}\cap\text{signed slopes})\subseteq\mathcal B(\text{N2})$, for both outer sets.
(d) **Nesting with N1** (superpopulation reading). If $m_{1c}$ and $m_{2c}$ are both nondecreasing, then $\delta^*_W\ge0$ for column $c$.
(e) For $R\ge3$, monotonicity of $m_{rc}$ in one coordinate does **not** sign $\operatorname{Cov}_w(m_{rc},X_s)$. So (d) has no general R×C analogue.

*Proof.*
1. (a): N2(σ) ⇒ Lemma N2.1 at each $x_j$, and monotonicity holds on the support. [L N2.1; Def] The target depends on $m$ only at the support points. [Def $\tilde B$]
2. (b): extend each $m_{rc}$ as a step function. [ext] The extension is Borel and monotone. Its values at the support points are unchanged, so Lemma N2.1 still holds there. [Def]
3. (c): for $K$ in the outer set, $M_{Q+K}(x)\ge0$ on $A\supseteq\operatorname{supp}$ [2a L3(a)] and its rows sum to 1 [2a L1(b)], so it lies in $\mathfrak B$. Also $x^\top M_{Q+K}(x)=x^\top Qx=\tau(x)$ [2a L1(a)], so $M_{Q+K}(x)\in S(x,\tau(x))$. A linear function with a slope of sign σ is σ-monotone. [Def]
4. (d): $\delta^*_W=\mu\operatorname{Cov}_{w^1}(m_1(X),X)+(1-\mu)\operatorname{Cov}_{w^2}(m_2(X),X)$, with both covariances $\ge0$. [Def 1 with $m$; Cheb applied under $w^1$ and under $w^2$]
5. (e): counterexample. With $R=3$, take $m_{rc}(x)=x_3$. It is constant, hence weakly nondecreasing, in $x_1$. Yet $\operatorname{Cov}(X_3,X_1)<0$ under a Dirichlet design. [C5(c)] ∎

Checked:
* C5(a): N2 contains the truth under a monotone nonlinear DGP.
* C5(b): the nesting (c) held in 47 random cases.
* C5(c): the counterexample, with Cov = −0.019.
* 1b V5: Chebyshev.

**Proposition N2.3 (local sign and local CR; EGGHM Props 12–13 in this setting).** Suppose $m_{1c}$ and $m_{2c}$ are differentiable at $x$. Write $D_c(x)=m_{1c}(x)-m_{2c}(x)$ and $\delta_{W,c}(x)=x\,m_{1c}'(x)+(1-x)m_{2c}'(x)$. Then:
(a) $\tau_c'(x)=D_c(x)+\delta_{W,c}(x)$.
(b) $\delta_{W,c}(x)\ge0\iff D_c(x)\le\tau'_c(x)$. This is linear in $m(x)$ given $\tau'$.
(c) Local CR, $D_c(x)\delta_{W,c}(x)\ge0$, holds iff $D_c(x)$ lies between 0 and $\tau'_c(x)$. This is convex given $\tau'$, because the branch is fixed by the identified sign of $\tau'_c(x)$.
(d) Both rows nondecreasing ⇒ local sign $\ge0$. The local sign everywhere does **not** imply $\delta^*_W\ge0$ (1b V6).

*Proof.*
1. Differentiate $\tau_c(x)=x\,m_{1c}(x)+(1-x)m_{2c}(x)$. [Lemma N2.1; product rule] This gives (a), and (b) follows. [alg]
2. (c), case $\tau'\ge0$: if $D\ge0$ and $\delta_W\ge0$, then $0\le D=\tau'-\delta_W\le\tau'$. If $D\le0$ and $\delta_W\le0$, then $D=\tau'-\delta_W\ge\tau'\ge0\ge D$, so $D=0$. Conversely, $D\in[0,\tau']$ gives $D\ge0$ and $\delta_W=\tau'-D\ge0$. The case $\tau'\le0$ is symmetric. [a; alg]
3. (d): $x,1-x\ge0$. [Def] ∎

Checked: C6 (500 random draws). The LP version imposes (b) or (c) at the support points with an estimated $\tau'$. It does not impose differentiability on $m$, so it is again an outer set.

---

## 5. The ladder (2 rows), with status

$$DD\ \supseteq\ \text{N1}(\delta_B\text{ sign})\ ,\ \text{N1(CR)}\ ,\ \text{N1}(\delta_W\text{ sign})\ \supseteq\ \text{N1(both signs)}\qquad\text{[M2; nest]}$$
$$\text{N2(both rows }\uparrow)\ \subseteq\ \text{N2-local(sign}\ge0)\qquad\text{[N2.3(d)]}$$
$$\mathcal B(\text{M0}\cap\text{signed slopes})\ \subseteq\ \mathcal B(\text{N2})\quad(\text{outer, }\tau\text{ quadratic})\qquad\text{[N2.2(c)]}$$
$$\text{N2(both rows }\uparrow)\ \Rightarrow\ \delta^*_W\ge0\qquad\text{[N2.2(d)]};\qquad \text{mean-level order}\ \Rightarrow\ \delta^*_B\ge0\qquad\text{[M5(b)]}.$$
**Not nested:**
* M0 with N1 or N2;
* N2-local with N1($\delta_W$) (1b V6);
* for $R\ge3$, N2 with N1($\delta_W$) (N2.2(e)).

---

## 6. What Phase 3 must supply

* **N1** needs no sampling inference (owner decision 1). Phase 3 adds the slack dial ε and the breakdown value (M2(d)).
* **K-rows:** the existing P3 machinery, with one extra row each.
* **N2 and N2-local** need $\tau$ (and $\tau'$) estimated at the support points.
  * Proposal: a weighted polynomial of degree $J$, selected by cross-validation; with $J=2$, N2(c) holds exactly.
  * The constraints are then linear in the coefficients, and P3 Thm 5 applies *conditional on $J$*: valid if $\tau$ is a polynomial of degree $J$, a heuristic otherwise.
  * With many support points, group compositions into bins. Binning keeps the outer property: averaging Lemma N2.1 within a bin uses [CE-cvx].
* **Implementation sizes.** N1 has $p\,R\,C$ variables. That is about 370×6 for Wisconsin in 2×3, which `lpSolve` handles; the equalities will be eliminated to $p(R-1)(C-1)$ free variables.
