# Phase 2b: units nested in groups, and the assumption ladder

*Builds on `PHASE1_foundations.md` (P1) and `PHASE2a_rxc.md` (2a). Every proof step ends with a justification tag; tags are those of 2a §0.3 plus the new ones in §0.3. All closed forms, inclusions and examples are checked in `verify_phase2b.py` (output: `verify_phase2b_output.txt`, 21/21 PASS).*

*"District" means any partition of the units into groups: ridings, counties within states, schools within boards, or, as Remark 6.2 shows, election periods. The results do not use anything specific to voter transitions.*

---

## 0. Setting

### 0.1 Data and notation

* Observables per unit: $W=(N,X,T,D)$ with $D\in\{1,\dots,D^\star\}$ the group, finite, and $\pi_d=P(D=d)>0$. $V:=(X,N,D)$, $\mathcal G=\sigma(V)$.
* $A_d$ is a compact convex set with $P(X\in A_d\mid D=d)=1$; the default is $A_d=\operatorname{conv}\operatorname{supp}P_{X\mid D=d}$, the district's own composition range. In samples it is a polytope with vertices among the district's polls.
* Moments: $\mu^{(d)}_{rs}=E[NX_rX_s1\{D=d\}]$, $\mu^{(d)}_r=E[NX_r1\{D=d\}]$, and $\mu_{r}=\sum_d\mu^{(d)}_r>0$, $\mu_{rs}=\sum_d\mu^{(d)}_{rs}$.
* For $Z\in\mathbb R^{R}$, $Z\mathbf 1^\top$ is the $R\times R$ matrix whose row $r$ is constant and equal to $Z_r$. Since $\mathbf 1^\top x=1$ on the simplex, $(Z\mathbf 1^\top)x=Z$: an additive shift written in homogeneous form.

### 0.2 The district model and the rungs

**District model $\mathcal M(\mathbf H;\mathbf A)$**, with $\mathbf H=(H^{(1)},\dots,H^{(D^\star)})$ and each $H^{(d)}=(H^{(d)}_1,\dots,H^{(d)}_C)$: laws $\mathbb Q$ of $(W,\beta)$ with

* (S1) $\beta\in\mathfrak B$ and $T=\beta^\top X$ a.s.;
* (S2-d) $E[\beta\mid V]=M_{H^{(D)}}(X)$ a.s.;
* (S3-d) $M_{H^{(d)}}(x)\in\mathfrak B$ for all $x\in A_d$, all $d$.

Write $Q^{(d)}=\operatorname{sym}H^{(d)}$ and $K^{(d)}=\operatorname{skew}H^{(d)}$ column by column. Let $k^{(d)}\in\mathbb R^{d_0}$ be the free coordinates of $K^{(d)}$ (2a Prop. 7), $d_0=(C-1)R(R-1)/2$, and let $c_d$ denote the coordinates of $\operatorname{skew}(\zeta_{\cdot,d}\mathbf 1^\top)$.

**Rungs** (restrictions on $\mathbf H$; decreasing plausibility; P1 §7, `ASSUMPTION_MENU.md`):

| Rung | Restriction $\mathcal R_\rho$ on $\mathbf H$ |
|---|---|
| **M2** | none: every district has its own linear contextual effects |
| **M1($\delta$)** | $\exists\bar k$: $\lVert k^{(d)}-\bar k-c_d\rVert_\infty\le\delta$ for all $d$, where $\zeta_{c,d}:=\operatorname{diag}(Q^{(d)}_c-Q^{(1)}_c)$ |
| **M1** | $H^{(d)}_c=H_c+\zeta_{c,d}\mathbf 1^\top$ with $\sum_c\zeta_{c,d}=0$ and $\zeta_{c,1}=0$: shared contextual effects, district baselines |
| **L4($r$)** | M1 and, for the designated row $r$, $H_{c,rs}=H_{c,rr}$ for all $s,c$: row $r$ has no contextual effect |
| **L5** | $H^{(d)}_c=\alpha_{c,d}\mathbf 1^\top$: constancy within each district |
| **L6** | L5 with $\alpha_{c,d}$ the same for all $d$ (Goodman) |
| *M0* | *(not nested with the others)* $H^{(d)}=H$ for all $d$ with $A=\operatorname{conv}\bigcup_dA_d$: the pooled model of 2a |

The district $\delta$-ball in M1($\delta$) is centred on the M1 prediction $\bar k+c_d$. $K$ entries are half-gaps between two hypothetical rates (2a §6.2), so $\delta$ is measured in probability units: $\delta=0.05$ lets the gap between "$r$-voters in pure-$s$ polls" and "$s$-voters in pure-$r$ polls" differ across districts by up to 10 points.

**Identified sets.** For rung $\rho$, $\Theta^\rho_I(P)$ is the set of $(k^{(1)},\dots,k^{(D^\star)})$ (and $\bar k$ for M1($\delta$)) for which some $\mathbb Q\in\mathcal M(\mathbf H;\mathbf A)$ with $\mathbf H\in\mathcal R_\rho$, $Q^{(d)}=Q^{(d)}(P)$, reproduces $P$. $\Theta^\rho_O$ replaces the sharp condition by nonnegativity at the vertices of each $A_d$. $\mathcal B^\rho$ denotes the corresponding set of values of a linear target $\varphi$.

### 0.3 Additional tags

| Tag | Meaning |
|---|---|
| [2a:X] | result X of `PHASE2a_rxc.md` |
| [cond-D] | conditioning on $\{D=d\}$, an atom of positive probability: (S1)–(S3-d) restricted to $\{D=d\}$ are the 2a model for the conditional law $P(\cdot\mid D=d)$ |
| [Gale] | Gale's supply–demand theorem (equivalently max-flow/min-cut): a bipartite transportation problem is feasible iff every set of supply nodes demands no more than the capacity of its neighbours |
| [sep] | minimising a separable sum over a product set equals summing the separate minima |
| [nest] | if model class $a\subseteq b$, then $\Theta^a_I(P)\subseteq\Theta^b_I(P)$ and the identified sets of any functional are nested |
| [LP-par] | the optimal value of a feasible bounded LP is a continuous, piecewise-linear and convex (min) / concave (max) function of its right-hand side on the region of feasibility |

---

## 1. Without a model, grouping adds nothing

**Proposition D0.** Under (S1) alone, the sharp identified set for $B^*_{rc}$ is
$$DD^*_{rc}=\Big[\tfrac{E[N\max\{0,X_r+T_c-1\}]}{\mu_r},\ \tfrac{E[N\min\{X_r,T_c\}]}{\mu_r}\Big],$$
whatever the law of $D$ (or of any covariate). Any gain from district dummies therefore comes entirely from a cross-district restriction.

*Proof.*
1. Every $\mathbb Q$ satisfying (S1) gives $B^*_{rc}\in DD^*_{rc}$. [2a:P8 steps 3–5]
2. Fix $(x,t)$ and put $m=\min\{x_r,t_c\}$. Consider the transportation problem with row margins $x$, column margins $t$ and cell $(r,c)$ fixed at $m$. The residual margins $x-me_r$ and $t-me_c$ are nonnegative with equal totals. [alg]
3. The residual problem has all edges present, so Gale's condition reduces to equal totals, and it is feasible. Hence $\max_{b\in S(x,t)}x_rb_{rc}=m$. [Gale; 2a:L0]
4. Put $m'=\max\{0,x_r+t_c-1\}$ and forbid further mass in cell $(r,c)$. The residual row $r$ needs $x_r-m'\le\sum_{c'\ne c}t_{c'}=1-t_c$, which holds because $m'\ge x_r+t_c-1$. Symmetrically, column $c$ needs $t_c-m'\le1-x_r$. [alg]
5. Any other set of rows is adjacent to every column, so Gale's condition holds for it automatically. So the residual problem is feasible and $\min_{b\in S(x,t)}x_rb_{rc}=m'$. [Gale]
6. A Borel selection attaining the maximum (or minimum) of $x_rb_{rc}$ over $S(x,t)$ exists. [2a tag MMT; 2a:L0(c)]
7. Integrating $N\times$ that selection gives the endpoints of $DD^*_{rc}$. Both are attained by valid structures (S1), with $\beta$ a function of $W$. [3, 5, 6; lin]
8. Neither the constraint $\beta\in S(X,T)$ nor the objective involves $D$, so the set does not depend on the law of $D$. [Def] ∎

**Cell intervals versus the joint set.** $DD^*_{rc}$ is the sharp set for one cell. The joint sharp set for the whole matrix $B^*$ is
$$\mathcal B^*_{DD}=\Big\{\big[E[NX_r\xi_{rc}]/\mu_r\big]_{rc}:\ \xi\text{ a measurable selection of }S(X,T)\Big\},$$
also free of $D$. It is convex and compact. It is **not** the Cartesian product $\prod_{rc}DD^*_{rc}$: the accounting identity links the cells of each unit, so in general the per-cell endpoints cannot be attained together. For a linear target $\sum a_{rc}B^*_{rc}$, the sharp interval is $\big[-E[N\,h_S(-\tilde a(X);X,T)],\ E[N\,h_S(\tilde a(X);X,T)]\big]$ with $\tilde a(x)_{rc}=a_{rc}x_r/\mu_r$ (the off-$A$ term of 2a Prop. 7′, applied to every unit), which can be strictly inside the sum of the cell intervals. Verified: per-cell extremes are attained in 300 random tables.

---

## 2. What each rung identifies

**Lemma D1.** Assume (A5-d): for each $d$, $\operatorname{supp}(P_{X\mid D=d}|_{A_d})$ lies in no quadric of the simplex. Here $\operatorname{supp}(P_{X\mid D=d}|_{A_d})$ is the support of the restricted measure (2a §0.1); since $P(X\in A_d\mid D=d)=1$ (§0.1), it equals $\operatorname{supp}P_{X\mid D=d}$. *(Corrected after review, Batch B, Nov 2026: earlier written $\operatorname{supp}P_{X\mid D=d}\cap A_d$, which differs when $A_d$ carries less than full mass; see 2a §0.1.)*
(a) Under any rung, $Q^{(d)}$ is identified from $E[T\mid X,D=d]=X^\top Q^{(d)}X$, $\sum_cH^{(d)}_c=J$, and each $K^{(d)}$ is unidentified with $K^{(d)}\in\mathcal K$.
(b) Under M1: $\zeta_{c,d}=\operatorname{diag}(Q^{(d)}_c-Q^{(1)}_c)$; the off-diagonal entries satisfy the **overidentifying restrictions**
$$(Q^{(d)}_c-Q^{(1)}_c)_{rs}=\tfrac12(\zeta_{c,d,r}+\zeta_{c,d,s})\quad(r\ne s);$$
and $K^{(d)}=K+\operatorname{skew}(\zeta_{\cdot,d}\mathbf 1^\top)$ with $K=K^{(1)}$ common. The unidentified dimension is $d_0$, the same as the pooled model.
(c) Under L5: $\alpha_{c,d}=\operatorname{diag}Q^{(d)}_c$, $K^{(d)}=\operatorname{skew}(\alpha_{\cdot,d}\mathbf 1^\top)$ is identified, and $Q^{(d)}_c=\operatorname{sym}(\alpha_{c,d}\mathbf 1^\top)$ is testable.

(d) **Under M1, less richness suffices (revised after review).** Replace (A5-d) by the stacked rank condition (A5-M1): for symmetric $D$ and vectors $z_2,\dots,z_{D^\star}$,
$$x^\top Dx+z_d^\top x=0\ \ \forall x\in\operatorname{supp}(P_{X\mid D=d}|_{A_d}),\ \forall d\quad(z_1:=0)\ \Longrightarrow\ D=0,\ z_d=0.$$
This is the population version of "the design with the pooled quadratic monomials and the district × $X_r$ shifts ($d\ge2$) has full column rank". Then $Q$, the $\zeta_{c,d}$ and the conclusions of (b) hold. A sufficient condition: A5 holds within one district, and every other district's support contains $R$ affinely independent points. Having $R$ *distinct* compositions is not enough: they must be affinely independent.

*Proof of (a).*
1. On $\{D=d\}$, (S1)–(S3-d) are the 2a model with $H=H^{(d)}$ and $A=A_d$ for the law $P(\cdot\mid D=d)$. [cond-D]
2. Apply 2a Lemma 1 to that conditional law. [2a:L1; A5-d] ∎(a)

*Proof of (d).*
1. Under M1, $E[T_c\mid X=x,D=d]=x^\top Q_cx+x^\top(\zeta_{c,d}\mathbf 1^\top)x=x^\top Q_cx+\zeta_{c,d}^\top x$ on the district's support. [2a:L1(a) steps 1–6 within $d$; Def M1; $\mathbf 1^\top x=1$]
2. If two parameter sets give the same regression in every district, their differences $(D,z_d)$ satisfy the premise of (A5-M1), so they coincide. [1; A5-M1]
3. Sufficiency: if A5 holds in district $d_0$, then for $d_0=1$ the premise gives $x^\top Dx=0$ there, so $D=0$. For $d_0\ne1$, $x^\top Dx+z_{d_0}^\top x=x^\top(D+\operatorname{sym}(z_{d_0}\mathbf 1^\top))x=0$ gives $D=-\operatorname{sym}(z_{d_0}\mathbf 1^\top)$, and district 1's equation then forces $z_{d_0}^\top x=0$ on an affinely spanning set, so $z_{d_0}=0$ and $D=0$. Then $z_d^\top x=0$ on each district's affinely spanning support. Such a support spans $\mathbb R^R$ linearly, because it lies in $\{\mathbf 1^\top x=1\}$, which misses the origin. So $z_d=0$. [alg; A5]
4. (b)'s formulas follow as before, from the identified $Q$ and $\zeta$. ∎(d)

*Proof of (b).*
1. $H^{(d)}_c=H_c+\zeta_{c,d}\mathbf 1^\top$. [Def M1]
2. $\operatorname{sym}(\zeta\mathbf 1^\top)_{rs}=\tfrac12(\zeta_r+\zeta_s)$ and $\operatorname{skew}(\zeta\mathbf 1^\top)_{rs}=\tfrac12(\zeta_r-\zeta_s)$. [alg]
3. So $Q^{(d)}_c=Q_c+\operatorname{sym}(\zeta_{c,d}\mathbf 1^\top)$ and $K^{(d)}_c=K_c+\operatorname{skew}(\zeta_{c,d}\mathbf 1^\top)$, with $Q=Q^{(1)}$ and $K=K^{(1)}$ because $\zeta_{\cdot,1}=0$. [1–2; lin]
4. Taking diagonals in 3 gives $\operatorname{diag}(Q^{(d)}_c-Q^{(1)}_c)=\zeta_{c,d}$, identified by (a). [2–3]
5. The off-diagonal entries of 3 give the overidentifying restrictions. [2–3]
6. The only unidentified object left is $K$; $\sum_c\zeta_{c,d}=0$ keeps $K^{(d)}\in\mathcal K$. [3; (a)] ∎(b)

*Proof of (c).*
1. $H^{(d)}_c=\alpha_{c,d}\mathbf 1^\top$, so $Q^{(d)}_c=\operatorname{sym}(\alpha_{c,d}\mathbf 1^\top)$, whose diagonal is $\alpha_{c,d}$. [Def L5; (b)2]
2. Hence $K^{(d)}_c=\operatorname{skew}(\alpha_{c,d}\mathbf 1^\top)$ is a known function of identified quantities. The off-diagonal entries of 1 are testable restrictions. [1; (a)] ∎(c)

**Remark D1.1 (why (d) matters in practice).** Given the common $Q$, district $d$'s regression is $E[T_c\mid X,d]=X^\top Q_cX+\zeta_{c,d}^\top X$, which is linear in $\zeta_{c,d}$. So under M1 (Lemma D1(d)), $\zeta_{c,d}$ needs only $R$ affinely independent compositions in district $d$, rather than A5-d. This is identification, not precision: with about 10 polls in a riding, $\hat\zeta_d$ is identified but noisy (Phase 3, regime R2). For example, 2 distinct polls per district in 2×2, provided $Q$ is identified elsewhere. This matters with about 10 polls per riding.

Verified on the 3×3 example: $\zeta$ recovered to $10^{-15}$; the overidentifying restrictions and the $K^{(d)}$ relation hold to $10^{-15}$.

**Lemma D2 (estimand).** Under (S1)–(S2-d) with $P(X\in A_D)=1$:
$$B^*_{rc}=\frac1{\mu_r}\sum_d\Big[\sum_sQ^{(d)}_{c,rs}\mu^{(d)}_{rs}+\sum_{s\ne r}K^{(d)}_{c,rs}\mu^{(d)}_{rs}\Big].$$
Under M1, its sensitivity to the common $K_{c,rs}$ is $\mu_{rs}/\mu_r$, exactly as in the pooled model.
*Proof.*
1. $E[NX_r\beta_{rc}]=\sum_dE[NX_r\beta_{rc}1\{D=d\}]$. [lin; $D$ finite]
2. $=\sum_dE[NX_r1\{D=d\}\,e_r^\top H^{(d)}_cX]$. [CE-tower; CE-pull ($NX_r1\{D=d\}$ is $V$-measurable); S2-d]
3. $=\sum_d\sum_sH^{(d)}_{c,rs}\mu^{(d)}_{rs}$, then split $H=Q+K$. [lin; skew]
4. Under M1, $K^{(d)}=K+\operatorname{skew}(\zeta\mathbf 1^\top)$, and $\sum_d\mu^{(d)}_{rs}=\mu_{rs}$. [D1(b); Def] ∎

---

## 3. Sharp and outer sets for each rung

**Theorem 1′ (sharpness for any affine mean model).** Let $\theta$ range over a closed subset $\Theta_0$ of $\mathbb R^p$. Let $\mathcal V_A=\{v=(x,n,d):x\in A_d\}$, and assume:
* (H1) for each $\theta$, $v\mapsto m_\theta(v)\in\mathbb R^{R\times C}$ is Borel on $\mathcal V_A$;
* (H2) $\sup_{v\in\mathcal V_A}\lVert m_\theta(v)\rVert<\infty$ for each $\theta$;
* (H3) $\theta\mapsto m_\theta(v)$ is affine for each $v$, and continuous uniformly in $v$ on bounded sets of $\theta$;
* (H4) the selection set $\mathcal S$ of 2a Theorem 1 is unchanged: $S(x,t)$ is nonempty, compact and convex, with a measurable graph (2a Lemma 0), and is not affected by $\theta$.

Define the model by (S1), $E[\beta\mid V]=m_\theta(V)$ a.s. on $\{V\in\mathcal V_A\}$, and $m_\theta(v)\in\mathfrak B$ for all $v\in\mathcal V_A$. Then $\theta\in\Theta_I$ iff the latter holds and, for $P_V$-a.e. $v$ and all $u$,
$$\langle u,m_\theta(v)\rangle\le\psi(u;v).$$
*Proof.* Replace $M_H(X)$ by $m_\theta(V)$ in the proof of 2a Theorem 1. (H4) leaves facts F1–F7 there unchanged. That proof uses only boundedness (H2) and $V$-measurability (H1) of $M_H(X)$ (steps 6, 8, 15, 16 and 21 of "(ii)⇒(i)"; steps 6–7 of "(i)⇒(ii)"), and affinity and continuity in the parameter (H3; the convexity and closedness part). Each step therefore carries over with the same tag. For the rungs of Theorem 2, $m_\theta(x,n,d)=M_{Q^{(d)}+K^{(d)}(\theta)}(x)$ is linear in $x$, so (H1)–(H2) hold on compact $A_d$, and (H3) holds because $K^{(d)}(\theta)$ is linear in $\theta$. [2a:T1] ∎

**Theorem 2 (the rungs).** Assume (A1), (A5-d) and $P(X\in A_D)=1$.
(a) **Sharp set.** For each rung $\rho$,
$$\Theta^\rho_I=\Big\{\theta\in\mathcal R_\rho:\ \text{for each }d,\ \langle u,M_{Q^{(d)}+K^{(d)}(\theta)}(x)\rangle\le\psi(u;x,n,d)\ \ \forall u,\ P_V\text{-a.e. }(x,n)\text{ given }D=d\Big\}.$$
That is the rung restriction intersected with each district's own sharp condition (2a Theorem 1 within district $d$).
(b) **Outer set.**
$$\Theta^\rho_O=\Big\{\theta\in\mathcal R_\rho:(Q^{(d)}_c+K^{(d)}_c(\theta))v\ge0\ \ \forall d,c,\ v\in\operatorname{vert}(A_d)\Big\}\supseteq\Theta^\rho_I.$$
It is a polyhedron, and every target bound is an LP. For M1 the variables are just the common $k\in\mathbb R^{d_0}$ and
$$\Theta^{M1}_O=\bigcap_d\Theta^{(d)}_O,\qquad\Theta^{(d)}_O=\{k:(Q^{(d)}_c+\operatorname{skew}(\zeta_{c,d}\mathbf 1^\top)+K_c(k))v\ge0\ \ \forall c,\ v\in\operatorname{vert}A_d\}.$$
(c) **Ladder.** For $0\le\delta\le\delta'$:
$$\mathcal B^{L6}\subseteq\mathcal B^{L5}\subseteq\mathcal B^{L4(r)}\subseteq\mathcal B^{M1}\subseteq\mathcal B^{M1(0)}\subseteq\mathcal B^{M1(\delta)}\subseteq\mathcal B^{M1(\delta')}\subseteq\mathcal B^{M2}\subseteq DD^*,$$
for both the sharp and the outer versions, whenever the smaller rung's model is non-empty.

*Proof of (a).*
1. Each rung is an affine mean model in the sense of Theorem 1′, with $\theta=(k^{(d)})_d$ (plus $\bar k$) restricted to the polyhedron $\mathcal R_\rho$, and $m_\theta(x,n,d)=M_{Q^{(d)}+K^{(d)}(\theta)}(x)$. [Def; D1(a); lin]
2. Apply Theorem 1′. The condition holds $P_V$-a.e. iff it holds a.e. given each atom $D=d$. [T1′; cond-D] ∎(a)

*Proof of (b).*
1. Within district $d$, sharp ⊆ outer by 2a Prop. 1 applied to the conditional law, with vertex reduction by 2a Lemma 3(a). [cond-D; 2a:P1; 2a:L3]
2. Intersecting over $d$ and with $\mathcal R_\rho$ preserves the inclusion. [1]
3. All constraints are linear in $\theta$, and $\mathcal R_\rho$ is polyhedral for every rung (an $\ell_\infty$ ball is a polyhedron). The target is affine in $\theta$ by Lemma D2. So the bounds are LPs. [D2; Def]
4. For M1, $K^{(d)}(k)=K(k)+\operatorname{skew}(\zeta\mathbf 1^\top)$ with $\zeta$ identified, so the constraints of district $d$ involve only the common $k$. [D1(b)] ∎(b)

*Proof of (c).*
1. **M2 ⊆ DD:** the M2 model class is contained in the (S1)-only class. [Def; nest; D0]
2. **M1($\delta$) ⊆ M1($\delta'$) ⊆ M2:** the $\delta$-ball is contained in the $\delta'$-ball, and M2 has no restriction. [Def; nest]
3. **M1 ⊆ M1(0):** under M1, $k^{(d)}=k(K)+c_d$ by D1(b), so $\bar k=k(K)$ satisfies the M1(0) condition. [D1(b); nest]
4. **L4($r$) ⊆ M1:** L4 is M1 plus a further restriction. [Def; nest]
5. **L5 ⊆ L4($r$):** under L5 every row is constant in $x$. Writing $\alpha_{c,d}=\alpha_{c,1}+(\alpha_{c,d}-\alpha_{c,1})$ shows L5 is M1 with $H_c=\alpha_{c,1}\mathbf 1^\top$ and $\zeta_{c,d}=\alpha_{c,d}-\alpha_{c,1}$, where $\sum_c\zeta_{c,d}=\mathbf 1-\mathbf 1=0$; and row $r$ is constant. [Def; alg]
6. **L6 ⊆ L5:** L6 is L5 with $\alpha$ common. [Def; nest]
7. **Outer versions:** each step adds constraints to, or restricts the variables of, the same polyhedral system, so the feasible sets and hence the LP ranges are nested. [(b); Def] ∎

Verified on the 3×3 example (§6): the chain DD ⊇ M2 ⊇ M1(.10) ⊇ M1(.05) ⊇ M1(.02) ⊇ M1 holds for all four targets, and the truth lies in every rung.

---

## 4. Closed forms and properties

**Proposition D3 (M1 with two rows).** Let $R=2$, $\operatorname{vert}(A_d)=\{(t,1-t):t\in\{l_d,u_d\}\}$ with $0<l_d\le u_d<1$ (needed: the formula divides by $t$ and $1-t$; pure-composition vertices are treated below), and let $\lambda_c,\rho_c$ be defined from the reference $Q=Q^{(1)}$ as in 2a Prop. 6. Under M1:
$$\ell_c=\max_d\max_{t\in\{l_d,u_d\}}\Big[\lambda_c(t)-\frac{\zeta_{c,d,1}}{1-t}\Big],\qquad h_c=\min_d\min_{t\in\{l_d,u_d\}}\Big[\rho_c(t)+\frac{\zeta_{c,d,2}}{t}\Big],$$
and $\Theta^{M1}_O=\{\kappa:\ell\le\kappa\le h,\ \sum_c\kappa_c=0\}$, projected by 2a Lemma 4.

*Derivation.*
1. Row 1 in district $d$ at $v=(t,1-t)$: $e_1^\top(H_c+\zeta_{c,d}\mathbf 1^\top)v=Q_{c,11}t+(Q_{c,12}+\kappa_c)(1-t)+\zeta_{c,d,1}\ge0$. [Def M1; $\mathbf 1^\top v=1$]
2. $\iff\kappa_c\ge-\frac{Q_{c,11}t+Q_{c,12}(1-t)+\zeta_{c,d,1}}{1-t}=\lambda_c(t)-\frac{\zeta_{c,d,1}}{1-t}$. [alg; 2a:P6]
3. Row 2: $(Q_{c,12}-\kappa_c)t+Q_{c,22}(1-t)+\zeta_{c,d,2}\ge0\iff\kappa_c\le\rho_c(t)+\frac{\zeta_{c,d,2}}{t}$. [alg; 2a:P6]
4. Intersecting over the vertices and the districts gives the max/min. [D3 statement; Thm 2(b)]
5. Adding $\sum_c\kappa_c=0$ gives a box intersected with a hyperplane. [Def $\mathcal K$; 2a:L4] ∎

Verified against the LP in 126 random column instances.

**Pure-composition vertices.** If $u_d=1$ (a poll with only group 1), the row-1 constraint at $t=1$ reads $Q_{c,11}+\zeta_{c,d,1}\ge0$. It does not involve $\kappa_c$: it is a testable restriction, and if it holds it adds no bound on $\kappa_c$. The term $t=1$ is therefore dropped from the maximum defining $\ell_c$. The row-2 constraint at $t=1$ is $\kappa_c\le Q_{c,12}+\zeta_{c,d,2}=\rho_c(1)+\zeta_{c,d,2}$, as the formula gives. Symmetrically, if $l_d=0$, the row-2 constraint at $t=0$ becomes the testable $Q_{c,22}+\zeta_{c,d,2}\ge0$ and is dropped from the minimum defining $h_c$, while row 1 gives $\kappa_c\ge\lambda_c(0)-\zeta_{c,d,1}=-Q_{c,12}-\zeta_{c,d,1}$. With $\zeta=0$ the same rule applies to 2a Prop. 6.
*Derivation.* Set $t=1$ (or $t=0$) in steps 1 and 3 above before dividing; the coefficient of $\kappa_c$ vanishes in exactly one of the two rows. [alg]

**Reading the formula.** A district raises the lower bound when its group-1 baseline for column $c$ is low ($\zeta_{c,d,1}<0$) and it contains polls that are mostly group 1 ($t$ near 1, where $1/(1-t)$ is large). In such polls group 1's rate is nearly observed, the shared contextual effect carries that information to every district, and so **one riding with lopsided polls and an extreme baseline can tighten the province-wide bound.** The upper bound behaves symmetrically, through group 2.

**Proposition D4 (when dummies matter).**
(a) If $\zeta\equiv0$, then $\Theta^{M1}_O=\Theta^{M0}_O(\operatorname{conv}\bigcup_dA_d)$: dummies neither help nor hurt identification.
(b) Under M1 with $P(X\in A_D)=1$, write $Z_d:=[\zeta_{c,d,r}]_{rc}$. Here M0 means the pooled model of 2a, a restriction on $E[\beta\mid X,N]$, which does not involve $D$. Then
$$E[\beta\mid X,N]=M_H(X)+E[Z_D\mid X,N],$$
and M0 holds (for some $H'$) iff there is a $G$ with $E[Z_D\mid X,N]=M_G(X)$ a.s.; that is, iff $E[Z_D\mid X,N]$ is a.s. equal to a function of $X$ alone that is affine on $\operatorname{supp}P_X$. In particular $D\perp(X,N)$ suffices, and so does $D\perp N\mid X$ together with $E[Z_D\mid X]$ affine. **$D\perp X$ alone does not suffice** (Remark D4.1). When M0 holds this way, its (S3) on $A=\operatorname{conv}\operatorname{supp}P_X$ is automatic by 2a Lemma 3(c). For any other pooled region $A$ it remains a separate condition.

*Proof of (a).*
1. $\zeta\equiv0$ gives $Q^{(d)}=Q$ for all $d$, so the pooled regression function is also $X^\top QX$. [D1(b); CE-tower]
2. The M1 constraints then require $(Q_c+K_c)v\ge0$ for $v\in\bigcup_d\operatorname{vert}(A_d)$. By linearity that is equivalent to the constraint on $\operatorname{conv}\bigcup_dA_d$, i.e. the M0 outer set. [Thm 2(b); 2a:L3(a)] ∎(a)

*Proof of (b).*
1. $E[\beta\mid X,N]=E\big[E[\beta\mid X,N,D]\mid X,N\big]=E[M_H(X)+Z_D\mid X,N]=M_H(X)+E[Z_D\mid X,N]$. [CE-tower, $\sigma(X,N)\subseteq\sigma(V)$; S2-d with M1; 0.1, $(Z\mathbf 1^\top)x=Z$; CE-pull, $M_H(X)$ is $\sigma(X,N)$-measurable]
2. On the simplex, an affine function of $x$ equals a homogeneous linear one: $a+b^\top x=(a\mathbf 1+b)^\top x$. Every linear map $\mathbb R^R\to\mathbb R^{R\times C}$ is $x\mapsto M_G(x)$ for some $G$. [alg; $\mathbf 1^\top x=1$]
3. "If": if $E[Z_D\mid X,N]=M_G(X)$ a.s., then $E[\beta\mid X,N]=M_{H+G}(X)$, which is (S2) of 2a with $H'=H+G$. [1–2]
4. "Only if": if M0 holds with $H'$, then $E[Z_D\mid X,N]=M_{H'}(X)-M_H(X)=M_{H'-H}(X)$ a.s., a function of $X$ only, linear in $X$. [1; 2a Def of (S2)]
5. If $D\perp(X,N)$, then $E[Z_D\mid X,N]=E[Z_D]$, a constant, hence affine in $X$. If $D\perp N\mid X$, then $E[Z_D\mid X,N]=E[Z_D\mid X]$. [indep.]
6. (S3) on $\operatorname{conv}\operatorname{supp}P_X$: under M0 with $H'$, (S1)–(S2) hold, so 2a Lemma 3(c) applies. [2a:L3(c)] ∎(b)

**Remark D4.1 ($D\perp X$ is not enough; corrected after review).** An earlier version conditioned on $X$ only and claimed $D\perp X$ suffices. The size-weighted estimand of 2a Lemma 2 uses $E[\beta\mid X,N]$, not $E[\beta\mid X]$. Counterexample: two districts with the same distribution of $X$, district 2 has larger precincts, and $\zeta_{\cdot,2}\neq0$. Then $E[Z_D\mid X,N]=Z_2\,P(D=2\mid N)$ varies with $N$, so (S2) of M0 fails. The $X$-only pooled model of $\Theta_X$ still holds, since $E[Z_D\mid X]=Z_2\pi_2$ is constant. But the identity $E[NX_r\beta_{rc}]=E[NX_rM_{H'}(X)_{rc}]$ behind Lemma 2 fails whenever size is correlated with the district shift, so pooled bounds on $B^*$ can be biased. The violation is generally visible in the data: $E[T\mid X,N]$ then depends on $N$ unless $x^\top Z_2=0$, which (T1) of 2a Theorem 1 rules out. It is checked numerically in `verify_phase2b.py`.

*Meaning:* for the conditional-mean model, dummies are needed exactly when the expected district baseline, given composition and size, either varies non-affinely with composition (e.g. stronghold districts with different baselines) or varies with precinct size. The pooled *outer* bounds additionally need probability validity on the pooled region, which is automatic for the population hull and an assumption otherwise. When M0 fails, the pooled quadratic fit is a projection with no bounding guarantee. In the 2×2 example (P1 §7.3) the pooled outer set is empty. In the 3×3 example it happens to cover the truth, which is not guaranteed.

**Proposition D5 (M2 adds district intervals).**
$$\mathcal B^{M2}_O=\sum_d\big[\varphi^{(d)}_L,\varphi^{(d)}_U\big],$$
where $\varphi^{(d)}$ is district $d$'s term in Lemma D2, bounded over $\Theta_O(P_d;A_d)$ (2a Prop. 7 within the district). The same holds for the sharp version.
*Proof.*
1. Under M2, the feasible set is the product $\prod_d\Theta^{(d)}$. [Thm 2; Def M2]
2. $\varphi=\sum_d\varphi^{(d)}(k^{(d)})$ is separable. [D2]
3. So the min (max) of the sum is the sum of the mins (maxes). [sep] ∎

Verified numerically. As a consequence, M2 widths add across districts. This is why M2 is usually little tighter than DD unless districts have wide, informative composition ranges.

**Proposition D6 (bounded deviation interpolates).**
(a) $\delta\mapsto\mathcal B^{M1(\delta)}_O$ is nested (shrinks as $\delta$ falls).
(b) For $\delta\ge\bar\delta:=\max_d\max_{k\in\Theta^{(d)}_O}\lVert k-c_d\rVert_\infty$ (finite), $\mathcal B^{M1(\delta)}_O=\mathcal B^{M2}_O$.
(c) The endpoints are continuous, piecewise-linear functions of $\delta$ wherever the set is non-empty.
(d) $\mathcal B^{M1(0)}_O=\mathcal B^{M1}_O$ computed with the same $Q^{(d)}$ and $\zeta$.
*Proof.*
1. (a) follows from Theorem 2(c), step 2. [Thm 2(c)]
2. (b): each $\Theta^{(d)}_O$ is bounded, so $\bar\delta<\infty$. Take $\bar k=0$; every point of $\prod_d\Theta^{(d)}_O$ then satisfies the $\delta$-condition, so the restriction is slack. [2a:L3(b); Def]
3. (c): $\delta$ enters only the right-hand side of the LP, linearly. [LP-par]
4. (d): at $\delta=0$ the M1($\delta$) rows become $k^{(d)}=\bar k+c_d$, which are the M1 equalities. The objective is identical. [Def; D2; Thm 2(b)] ∎

Note on (d): M1(0) does not impose the testable restrictions of D1(b). So it stays defined even when the data reject M1's symmetric part; it is "M1 without its testable implication". It is therefore a **sensitivity benchmark**, not a substantive rung on a par with M1: it separates a failure of M1's testable implications from uncertainty about the unidentified contextual effects. *(Emphasis added after review.)*

**Proposition D7 (the most restrictive rungs).**
(a) **L5** point-identifies the target:
$$B^*_{rc}=\sum_d\frac{\mu^{(d)}_r}{\mu_r}\,\alpha_{c,d,r},\qquad\alpha_{c,d}=\operatorname{diag}Q^{(d)}_c.$$
Equivalently, $\alpha_{\cdot,d}$ is the Goodman regression of $T$ on $X$ without intercept, within district $d$. If L5 holds, this point lies in every larger rung's bound.
(b) **L4($r$)** pins $K_{c,rs}=Q_{c,rr}-Q_{c,rs}$ (reference district) for all $s\ne r$ and all $c$, i.e. $(C-1)(R-1)$ of the $d_0$ coordinates. For $R=2$ the target is point-identified.
*Proof of (a).*
1. Under L5, $E[\beta_{rc}\mid V]=\alpha_{c,D,r}$, so $B^*_{rc}=\sum_dE[NX_r1\{D=d\}]\alpha_{c,d,r}/\mu_r$. [CE-tower; Def L5]
2. $\alpha_{c,d}=\operatorname{diag}Q^{(d)}_c$. [D1(c)]
3. Within district $d$, $E[T_c\mid X,d]=\sum_rX_r\alpha_{c,d,r}$. [2a:L1(a) proof steps 1–2; Def L5]
4. The point lies in every larger rung's bound. [Thm 2(c)] ∎(a)

*Proof of (b).*
1. $H_{c,rs}=H_{c,rr}$, and $H_{c,rr}=Q_{c,rr}$ because $K_{c,rr}=0$. [Def L4; skew]
2. $H_{c,rs}=Q_{c,rs}+K_{c,rs}$, hence $K_{c,rs}=Q_{c,rr}-Q_{c,rs}$. [Def; alg]
3. The pairs containing $r$ number $R-1$, in each of $C-1$ free columns. [alg]
4. For $R=2$ that is all $d_0=C-1$ coordinates. [alg] ∎(b)

Verified: under a DGP with zero contextual effects, the L5 point equals the truth and lies inside the M1 bound for all four 3×3 targets.

**Remark D7.1 (L5 as a diagnostic).** When contextual effects are not zero, the L5 "point" is a misspecified estimate. In the 3×3 example it is 0.870 for loyalty (truth 0.836): inside M1 but wrong. Reporting it next to the bounds shows how much the zero-contextual-effect assumption is doing.

---

## 5. Remarks on scope

**Remark 6.1 (estimation regimes).** Everything above is population identification. With about 10 polls per riding (BC stage 2), each $Q^{(d)}$ and $\zeta_d$ is noisy. The plug-in intersection in Theorem 2(b) then takes a max/min over many noisy constraints and is biased inward (P1 §7.1: 72% coverage). Remark D1.1 (fewer parameters per district under M1) and M1($\delta$) with $\delta>0$ (which softens the max) both help. A valid procedure needs Phase 3: precision-corrected intersection bounds (Chernozhukov, Lee & Rosen 2013) and inference for LP values with non-unique solutions. *(Revised: an earlier version also listed cluster-robust variances. Inference is conditional on the design and the district effects, with no clustered variance; P3 §0.4.)*

**Remark 6.2 (periods as districts: the DiD link).** Take "districts" to be election periods: the same counties observed twice, $D\in\{\text{pre},\text{post}\}$.
* **M1** then says the contextual effects are the same in both periods, with only baselines shifting, which is exactly the "contextual effects unchanged" assumption in `ASSUMPTION_MENU.md` item 18.
* Lemma D1(b) identifies the baseline shifts $\zeta$, and Proposition D8 below gives the exact condition under which the change $B^*_\text{post}-B^*_\text{pre}$ is point-identified. That recovers the owner's earlier result as a special case of M1.
* Theorem 2(b) adds that the *levels* in both periods are bounded by the intersection of the two periods' sets.
* **M1($\delta$)** is the bounded-change version (item 19).
* The overidentifying restrictions of D1(b) are the testable part, "equal quadratic coefficients across periods" (menu note to item 18).

**Proposition D8 (the DiD under M1 with periods as districts; promoted from Remark 6.2 after review).** Let $D\in\{\text{pre},\text{post}\}$, assume M1, (S1)–(S2-d), $P(X\in A_D)=1$ and (A5-M1). For period $t$ write $\mu^{[t]}_{rs}=E[NX_rX_s\mid D=t]$, $\mu^{[t]}_r=E[NX_r\mid D=t]$ and $B^*_{rc,t}=E[NX_r\beta_{rc}\mid D=t]/\mu^{[t]}_r$. Then
$$B^*_{rc,\text{post}}-B^*_{rc,\text{pre}}=\text{(identified)}+\sum_{s\ne r}K_{c,rs}\Big(\frac{\mu^{[\text{post}]}_{rs}}{\mu^{[\text{post}]}_r}-\frac{\mu^{[\text{pre}]}_{rs}}{\mu^{[\text{pre}]}_r}\Big).$$
So the change is point-identified, whatever the identified set for $K$, iff
$$\frac{\mu^{[\text{post}]}_{rs}}{\mu^{[\text{post}]}_r}=\frac{\mu^{[\text{pre}]}_{rs}}{\mu^{[\text{pre}]}_r}\qquad\text{for all }s\ne r.$$
This holds in particular when the size-and-composition moments are equal across the periods, $\mu^{[\text{post}]}=\mu^{[\text{pre}]}$. It is a condition on the $NX_rX_s$ and $NX_r$ moments, not merely on similar-looking distributions of vote shares.
*Proof.*
1. Within period $t$, Lemma D2 applied to the conditional law gives $B^*_{rc,t}=\big[\sum_sQ^{(t)}_{c,rs}\mu^{[t]}_{rs}+\sum_{s\ne r}K^{(t)}_{c,rs}\mu^{[t]}_{rs}\big]/\mu^{[t]}_r$. [D2; cond-D]
2. Under M1, $Q^{(t)}$ and $\operatorname{skew}(\zeta_{\cdot,t}\mathbf 1^\top)$ are identified, and $K^{(t)}=K+\operatorname{skew}(\zeta_{\cdot,t}\mathbf 1^\top)$ with $K$ common. [D1(b), D1(d)]
3. Subtracting the two periods, every term except $\sum_{s\ne r}K_{c,rs}(\cdot)$ is identified. [1–2; lin]
4. If every coefficient is zero, the change does not depend on $K$, so it is point-identified. If some coefficient is nonzero, the change moves with the corresponding entry of $K$, so it is point-identified only if the identified set happens to fix that linear combination. "Point-identified for every $K$" in the statement means: whatever the identified set for $K$ is. [3; alg] ∎
*Reading.* The DiD of transition rates is point-identified when the contextual effects are unchanged across periods (M1) and the electorate's size-weighted composition moments are the same in both periods. With changing composition the change is only bounded; the bound is the projection of the common identified set for $K$ (Theorem 2(b)) through the coefficients above.

**Remark 6.3 (other groupings).** Nothing above uses the nesting being geographic. Groups can be states (with counties as units), demographic strata of units, or time.

---

## 6. Worked examples (noise-free; identification only)

**2×2 (P1 §7.3, two ridings).** Province-wide P-voter loyalty, truth 0.819.

| Rung | Bound |
|---|---|
| DD | [0.453, 0.879] |
| M0 pooled | empty (reverts to DD) |
| M2 | [0.563, 0.870] |
| M1 | [0.618, 0.866] |
| L5 point | 0.858 (misspecified: slopes are 0.10) |

**3×3 (six districts with different composition ranges and baselines; common $H$ from 2a §8).** Each bound is [lower, upper].

| Target | Truth | DD | M0 pooled* | M2 | M1(.10) | M1(.05) | M1(.02) | M1 | L5 point |
|---|---|---|---|---|---|---|---|---|---|
| Loyalty P→P | 0.836 | [0.167, 0.953] | [0.423, 0.898] | [0.391, 0.916] | [0.392, 0.916] | [0.397, 0.914] | [0.412, 0.899] | [0.430, 0.880] | 0.870 |
| Defection P→O | 0.095 | [0.000, 0.610] | [0.026, 0.331] | [0.024, 0.404] | [0.024, 0.403] | [0.024, 0.390] | [0.024, 0.367] | [0.042, 0.339] | 0.061 |
| Demobilisation P→abstain | 0.068 | [0.000, 0.551] | [0.029, 0.296] | [0.024, 0.350] | [0.024, 0.339] | [0.024, 0.323] | [0.030, 0.309] | [0.044, 0.295] | 0.069 |
| Mobilisation N→P | 0.135 | [0.001, 0.685] | [0.038, 0.410] | [0.045, 0.470] | [0.045, 0.457] | [0.051, 0.444] | [0.070, 0.428] | [0.089, 0.413] | 0.122 |

*M0 is misspecified here (district baselines vary with composition). It happens to cover the truth and to be narrower than M2, but Prop. D4 gives no guarantee for it.

What the table shows:
* The ladder shrinks monotonically and always contains the truth.
* Moving from M2 to M1 narrows loyalty from width 0.525 to 0.450.
* δ traces a smooth path between the two rungs.
* The baseline shifts here are small (SD 0.04). The gain grows with more heterogeneous baselines and more lopsided polls (Prop. D3 reading).

---

## 7. Verification log

`python3 -I verify_phase2b.py` (21 PASS, 0 FAIL) checks:
* Lemma D1: recovery of $\zeta$, the overidentifying restrictions and the $K^{(d)}$ relation.
* Theorem 2: ladder nesting for four targets; the truth in every valid rung.
* Prop. D6: monotone in $\delta$; equals M2 for large $\delta$.
* Prop. D5: M2 equals the sum of district bounds.
* Prop. D4(a): $\zeta=0$ makes dummies equal to pooled.
* Prop. D7: the L5 point under a true L5.
* Prop. D3: closed form against the LP (126 instances).
* Prop. D0: per-cell Fréchet extremes attained (300 tables).
* The 2×2 example reproduces P1 §7.3.
* Corrections after review (each check fails on the earlier statement):
  * Prop. D4(b): with $D\perp X$ but $D$ correlated with $N$, $E[\beta\mid X,N]$ depends on $N$ while $E[\beta\mid X]$ is affine, and the pooled size-weighted rate is biased. With $D\perp(X,N)$ it is exact.
  * P1 (D2) and P3 Prop. P3.6(b): the rank gap between the M2 and M1 designs is $(D^\star-1)R(R-1)/2$ per column, and column $C$'s restrictions are implied. That gives $(D^\star-1)(C-1)R(R-1)/2$ in total (three table shapes).
* Scope additions after review:
  * Lemma D1(d): the stacked M1 design has full rank with A5 in one district and three affinely independent points in the others, and is rank-deficient with three distinct but collinear points;
  * Prop. D3 with pure-composition vertices: the modified closed form matches the LP (195 instances).

Not machine-checked: the measure-theoretic content of Theorem 1′, which rests on 2a Theorem 1.

## 8. Open items

* **Phase 3:** noisy-intersection inference (Remark 6.1); cluster dependence; the finite-population gap; a test of D1(b)'s overidentifying restrictions; sensitivity in $\delta$ (the breakdown value of a conclusion).
* **Not yet done:** M1′ (contextual effects depending on district covariates) is a polyhedral restriction like the others and fits Theorem 2 directly. It is left until a covariate set is chosen.
* **Strict-gain characterisation for $R\ge3$:** when $\Theta^{M1}$ is strictly smaller than the best single district's set, beyond the $R=2$ reading of Prop. D3.
