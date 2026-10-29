# Phase 2a: R×C bounds, derivations and proofs

*Builds on `PHASE1_foundations.md` (cited as P1). Every proof step ends with a justification tag in brackets; the tag legend is in §0.3. Every closed form and every inclusion claimed here is also checked numerically or symbolically in `verify_phase2a.py` (output in `verify_phase2a_output.txt`; 22/22 PASS). §10 lists the corrections this phase makes to P1.*

---

## 0. Setting, conventions and tags

### 0.1 Objects

* $\mathcal W=\mathbb N\times\Delta^{R-1}\times\Delta^{C-1}$ is the space of one unit's observables $W=(N,X,T)$; it is Polish. $V:=(X,N)$ and $\mathcal G:=\sigma(V)$. Covariates and districts come in Phase 2b.
* $\mathfrak B=(\Delta^{C-1})^R\subset\mathbb R^{R\times C}$, with $\langle u,b\rangle=\sum_{rc}u_{rc}b_{rc}$ and the Frobenius norm. Every $b\in\mathfrak B$ satisfies $\|b\|\le\sqrt R$.
* $J=\mathbf 1\mathbf 1^\top$. $\operatorname{Skew}_R=\{K\in\mathbb R^{R\times R}:K^\top=-K\}$.
* $\mathcal K=\{(K_1,\dots,K_C)\in\operatorname{Skew}_R^C:\sum_cK_c=0\}$.
* For $H=(H_1,\dots,H_C)$, $M_H(x)\in\mathbb R^{R\times C}$ has entries $M_H(x)_{rc}=e_r^\top H_cx$.
* $A\subseteq\Delta^{R-1}$ is compact and convex. The default is $A=\operatorname{conv}(\operatorname{supp}P_X)$. A *polytope* $A$ has a finite vertex set $\operatorname{vert}(A)$, and that is always the case for the sample hull.
* **Standing integrability:** $E[N]<\infty$, and $\mu_r:=E[NX_r]>0$ for all $r$. $\mu_{rs}:=E[NX_rX_s]$.

### 0.2 Model, identified sets (these refine P1 §1.5)

**Model $\mathcal M(H;A)$:** the set of laws $\mathbb Q$ of $(W,\beta)$ on $\mathcal W\times\mathfrak B$ such that

* (S1) $\beta\in\mathfrak B$ and $T=\beta^\top X$, $\mathbb Q$-a.s. *(A1, accounting identity)*;
* (S2) $E_{\mathbb Q}[\beta\mid V]=M_H(X)$ $\mathbb Q$-a.s. on $\{X\in A\}$ *(A2 + A3, linear contextual effects, no dependence on $N$)*;
* (S3) $M_H(x)\in\mathfrak B$ for every $x\in A$ *(A4: the linear model describes valid probabilities on all of $A$. Lemma 3(c) shows this is implied when $A\subseteq\operatorname{conv}\operatorname{supp}P_X$. It is a genuine extra assumption only when $A$ extrapolates beyond the data, as in JKST Prop. 3)*.

**Identified part.** $Q(P)=(Q_1,\dots,Q_C)$ denotes the unique symmetric matrices with $E_P[T_c\mid X=x]=x^\top Q_cx$ for $P_X$-a.e. $x\in A$, when they exist (Lemma 1). If they do not exist, set $\Theta_I(P;A):=\emptyset$.

**Sharp identified set:**
$$\Theta_I(P;A)=\{K\in\mathcal K:\ \exists\,\mathbb Q\in\mathcal M(Q(P)+K;A),\ \mathbb Q\circ W^{-1}=P\}.$$

**$\Theta_X(P;A)$** is the same with (S2) replaced by the weaker $E_{\mathbb Q}[\beta\mid X]=M_H(X)$; this is JKST Prop. 1.

**Outer set:**
$$\Theta_O(P;A)=\{K\in\mathcal K:\ (Q_c(P)+K_c)x\ge0\ \text{elementwise}\ \ \forall c,\ \forall x\in A\}.$$

**Unit Duncan–Davis polytope and its support function:**
$$S(x,t)=\{b\in\mathfrak B:b^\top x=t\},\qquad h_S(u;x,t)=\max_{b\in S(x,t)}\langle u,b\rangle.$$
$P_{T\mid V}$ is a regular conditional distribution of $T$ given $V$ under $P$, and
$$\psi(u;v):=\int h_S(u;x,t)\,P_{T\mid V=v}(dt).$$

**Estimand:** $B^*_{rc}=E[NX_r\beta_{rc}]/\mu_r$, and $\varphi=\sum_{rc}a_{rc}B^*_{rc}$.

### 0.3 Tag legend

Theorem numbers marked † are cited from memory and must be checked against the stated edition before publication. The mathematical content of each is standard.

| Tag | Meaning / source |
|---|---|
| [Def] | the definition just introduced or referenced |
| [S1], [S2], [S3], [A1]–[A5] | model conditions §0.2, assumptions P1 §2 |
| [lin] | linearity (of expectation, inner product, matrix products) |
| [alg] | elementary algebra, shown in the same line |
| [CE-lin], [CE-mono], [CE-pull], [CE-tower] | conditional expectation: linearity, monotonicity, "taking out what is known", tower property |
| [CE-cvx] | if $Y\in F$ a.s. with $F$ closed convex and $Y$ integrable, then $E[Y\mid\mathcal G]\in F$ a.s. (conditional Jensen / separation) |
| [RCD] | existence of regular conditional distributions on Polish spaces, and disintegration: $E[f(V,T)\mid V]=\int f(V,t)P_{T\mid V}(dt)$ a.s. for bounded jointly measurable $f$ (Kallenberg 2002, Thms 6.3–6.4†) |
| [supp] | if continuous $f,g$ agree $P_X$-a.e. on $A$, they agree on $\operatorname{supp}P_X\cap A$ |
| [Frob] | $\lvert\langle u,b\rangle\rvert\le\lVert u\rVert\,\lVert b\rVert$ (Cauchy–Schwarz) |
| [Mink] | a polytope is the convex hull of its vertices (Minkowski / Krein–Milman for polytopes) |
| [Rock] | closed convex sets are determined by their support functions (Rockafellar 1970, Cor. 13.1.1†) |
| [Berge] | Berge maximum theorem: max of a continuous function over a compact-valued u.h.c. correspondence is u.s.c. |
| [MMT] | measurable maximum theorem: a measurable argmax selector exists (Aliprantis & Border 2006, Thm 18.19†) |
| [DP] | Dunford–Pettis: uniformly integrable subsets of $L^1$ are relatively weakly compact |
| [Mazur] | convex norm-closed subsets of a normed space are weakly closed |
| [HB] | Hahn–Banach strict separation of a point from a closed convex set, with $(L^1)^*=L^\infty$ for a probability space |
| [wk-cont] | bounded linear operators are weak-to-weak continuous |
| [Riesz] | $L^1$ convergence implies a.s. convergence along a subsequence |
| [FH] | Fréchet–Hoeffding bounds for a cell of a table with given margins |
| [poly] | a polynomial in $k$ variables vanishing on a nonempty open set is the zero polynomial |
| [IVT] | intermediate value theorem / the image of a convex set under a linear map is convex |
| [exch] | exchange argument for linear objectives |
| [L0]–[L4], [T1], [P1]–[P9] | results of this document |

---

## 1. The unit Duncan–Davis polytope

**Lemma 0.** For all $x\in\Delta^{R-1}$ and $t\in\Delta^{C-1}$:
(a) $S(x,t)\neq\emptyset$;
(b) $S(x,t)$ is compact and convex;
(c) $(x,t)\mapsto S(x,t)$ has closed graph, is upper hemicontinuous and Borel measurable;
(d) $h_S(u;x,t)$ is sublinear in $u$, $\sqrt R$-Lipschitz in $u$, concave in $t$, and Borel in $(u,x,t)$;
(e) $S(x,t)=\{b:\langle u,b\rangle\le h_S(u;x,t)\ \forall u\}$.

*Proof.*
1. Put $b_{rc}:=t_c$. Each row is $t\in\Delta^{C-1}$, so $b\in\mathfrak B$. [Def]
2. $(b^\top x)_c=\sum_rx_rt_c=t_c$, using $\sum_r x_r=1$. [alg] So (a) holds. [1–2]
3. $\mathfrak B$ is a finite product of simplices, hence compact and convex. [Def]
4. $\{b:b^\top x=t\}$ is an affine subspace, hence closed and convex. [lin]
5. $S(x,t)$ is the intersection of 3 and 4, hence compact and convex. That is (b). [3–4]
6. The graph $\{(x,t,b):b\in\mathfrak B,\ b^\top x=t\}$ is the zero set of a continuous map intersected with a compact set, hence closed. [alg]
7. A closed-graph correspondence into the compact $\mathfrak B$ is upper hemicontinuous with compact values. Such correspondences are Borel measurable. That is (c). [6; Aliprantis & Border 2006, Ch. 17–18†]
8. $h_S$ is a maximum of the linear functions $u\mapsto\langle u,b\rangle$, hence convex and positively homogeneous in $u$, i.e. sublinear. [Def; lin]
9. $\lvert h_S(u)-h_S(u')\rvert\le\max_{b\in\mathfrak B}\lvert\langle u-u',b\rangle\rvert\le\sqrt R\,\lVert u-u'\rVert$. [Frob; $\lVert b\rVert\le\sqrt R$]
10. If $b\in S(x,t)$, $b'\in S(x,t')$ and $\lambda\in[0,1]$, then $\lambda b+(1-\lambda)b'\in\mathfrak B$ and $(\lambda b+(1-\lambda)b')^\top x=\lambda t+(1-\lambda)t'$. [3; lin]
11. Taking $b,b'$ to be maximisers, $h_S(u;x,\lambda t+(1-\lambda)t')\ge\lambda h_S(u;x,t)+(1-\lambda)h_S(u;x,t')$: concavity in $t$. [10; Def]
12. $h_S$ is u.s.c. in $(u,x,t)$, hence Borel. That completes (d). [Berge; 7]
13. (e) holds because $S(x,t)$ is closed and convex. [Rock; 5] ∎

---

## 2. What the data identify

**Lemma 1 (identification decomposition).** Let $\mathbb Q\in\mathcal M(H;A)$ with $\mathbb Q\circ W^{-1}=P$, and assume [A5]: for symmetric $D$, $x^\top Dx=0$ for all $x\in\operatorname{supp}P_X\cap A$ implies $D=0$. Then
(a) $E[T_c\mid X=x]=x^\top Q_cx$ with $Q_c=\tfrac12(H_c+H_c^\top)$, and $Q_c$ is unique;
(b) $\sum_cH_c=J$, so $\sum_cQ_c=J$ and $\sum_cK_c=0$ where $K_c=\tfrac12(H_c-H_c^\top)$;
(c) for every $K'\in\mathcal K$, $Q+K'$ gives the same regression function; $\dim\mathcal K=(C-1)R(R-1)/2$.

*Proof of (a).*
1. $T_c=\sum_rX_r\beta_{rc}$ a.s. [S1]
2. $E[T_c\mid X]=\sum_rX_rE[\beta_{rc}\mid X]$ a.s. [CE-lin; CE-pull, since $X$ is $\sigma(X)$-measurable]
3. $E[\beta_{rc}\mid X]=E\big[E[\beta_{rc}\mid V]\mid X\big]=e_r^\top H_cX$ a.s. on $\{X\in A\}$. [CE-tower; S2; $M_H(X)$ is $\sigma(X)$-measurable]
4. Hence $E[T_c\mid X]=\sum_{r,s}X_rH_{c,rs}X_s=X^\top H_cX$ on $\{X\in A\}$. [2–3; alg]
5. Since $X^\top K_cX$ is a scalar, $X^\top K_cX=(X^\top K_cX)^\top=X^\top K_c^\top X=-X^\top K_cX$, so $X^\top K_cX=0$. [alg; $K_c$ skew]
6. Therefore $X^\top H_cX=X^\top Q_cX$. [4–5; $H_c=Q_c+K_c$]
7. If symmetric $Q_c,Q'_c$ both satisfy 6 $P_X$-a.e. on $A$, then $x^\top(Q_c-Q'_c)x=0$ on $\operatorname{supp}P_X\cap A$. [supp; both sides continuous]
8. $Q_c-Q'_c$ is symmetric, so $Q_c=Q'_c$. [A5] ∎(a)

*Proof of (b).*
1. $\sum_c\beta_{rc}=1$ a.s. [S1; $\beta\in\mathfrak B$]
2. $\sum_ce_r^\top H_cx=1$ for $P_X$-a.e. $x\in A$. [CE-lin on 1; S2; CE-tower]
3. Let $d_r^\top:=e_r^\top(\sum_cH_c-J)$. Then $d_r^\top x=1-\mathbf 1^\top x=0$ on $\operatorname{supp}P_X\cap A$. [2; $\mathbf 1^\top x=1$; supp]
4. $x^\top\tfrac12(d_r\mathbf 1^\top+\mathbf 1d_r^\top)x=(d_r^\top x)(\mathbf 1^\top x)=0$ on the same set. [alg; 3]
5. So $d_r\mathbf 1^\top+\mathbf 1d_r^\top=0$. [A5; the matrix is symmetric]
6. Its $(s,s)$ entry is $2d_{r,s}=0$ for every $s$, so $d_r=0$. [alg]
7. Hence $\sum_cH_c=J$. Taking symmetric and skew parts, $\sum_cQ_c=J$ and $\sum_cK_c=0$. [6; $J$ symmetric; lin] ∎(b)

*Proof of (c).*
1. For $K'\in\mathcal K$: $x^\top(Q_c+K'_c)x=x^\top Q_cx$. [(a) step 5]
2. $(K_1,\dots,K_{C-1})\mapsto(K_1,\dots,K_{C-1},-\sum_{c<C}K_c)$ is a linear bijection from $\operatorname{Skew}_R^{C-1}$ onto $\mathcal K$. [Def; lin]
3. $\dim\operatorname{Skew}_R=R(R-1)/2$, so $\dim\mathcal K=(C-1)R(R-1)/2$. [2] ∎

**Remark 1.1 (A5 in practice).** For $R=2$, A5 holds iff $\operatorname{supp}P_X\cap A$ contains at least three distinct points.
*Proof.*
1. On $x=(t,1-t)$, $x^\top Dx$ is a polynomial of degree at most 2 in $t$. [alg]
2. If it vanishes at three distinct $t$, it is the zero polynomial. [poly]
3. If $D_{11}t^2+2D_{12}t(1-t)+D_{22}(1-t)^2\equiv0$, then $t=0$ gives $D_{22}=0$, $t=1$ gives $D_{11}=0$, and then $D_{12}=0$. [alg]
4. With only two support points $t_1,t_2$, the nonzero form $D$ given by $(t-t_1)(t-t_2)$ in homogenised form vanishes on both. [alg] ∎

For general $R$, A5 says that the design of the quadratic regression on the polls has full rank. It is checkable from the condition number of the regressor matrix (Phase 4).

**Remark 1.2 (sum-to-one in estimation).** Fit each column by OLS on the $R(R+1)/2$ quadratic monomials of $X$ (equivalently the homogeneous form). Then $\sum_c\hat Q_c=J$ exactly.
*Proof.*
1. $\sum_cT_{ic}=1=(\mathbf 1^\top X_i)^2=X_i^\top JX_i$ for every $i$. [Def; alg]
2. OLS is linear in the response, so $\sum_c\hat Q_c$ is the OLS fit of the response $1$. [lin]
3. By 1 that response is fitted exactly by $J$, and the fit is unique when the design has full rank (sample A5). [1; A5] ∎

So the estimated columns are automatically coherent; no constrained regression is needed.

---

## 3. The estimand in terms of $K$

**Lemma 2.** Under (S1)–(S2) with $P(X\in A)=1$:
$$B^*_{rc}=\frac1{\mu_r}\Big[\sum_sQ_{c,rs}\mu_{rs}+\sum_{s\neq r}K_{c,rs}\mu_{rs}\Big],\qquad\frac{\partial B^*_{rc}}{\partial K_{c,rs}}=\frac{\mu_{rs}}{\mu_r}\ge0.$$

*Proof.*
1. $0\le NX_r\beta_{rc}\le N$, so $E[NX_r\beta_{rc}]$ is finite. [S1; standing integrability]
2. $E[NX_r\beta_{rc}]=E\big[NX_rE[\beta_{rc}\mid V]\big]$. [CE-tower; CE-pull, since $NX_r$ is $V$-measurable]
3. $=E[NX_r\,e_r^\top H_cX]$. [S2; $P(X\in A)=1$]
4. $=\sum_sH_{c,rs}E[NX_rX_s]=\sum_sH_{c,rs}\mu_{rs}$. [lin]
5. $H_{c,rs}=Q_{c,rs}+K_{c,rs}$ and $K_{c,rr}=0$. [Def; skew]
6. Divide by $\mu_r>0$. [standing integrability]
7. The derivative is read off from the affine form. $\mu_{rs}\ge0$ because $N,X_r,X_s\ge0$. [alg] ∎

**2×2 check.** $B^*=[Q_{11}\mu_{11}+(Q_{12}+K)\mu_{12}]/\mu_1$. With the dictionary of Prop. 4 below and $\mu_{11}+\mu_{12}=\mu_1$, this is JKST eq. (16) in the population: $w_0+c_1+d_1\mu_{11}/\mu_1-w_1\mu_{12}/\mu_1$.

**Remark 2.1 ($A\subsetneq\operatorname{supp}$).** If $P(X\notin A)>0$, then
$$B^*_{rc}=\mu_r^{-1}\Big\{E\big[NX_rM_{H,rc}(X)1\{X\in A\}\big]+E\big[NX_r\beta_{rc}1\{X\notin A\}\big]\Big\},$$
and the second term is bounded only by the Duncan–Davis terms of the polls outside $A$. This is formalised in Phase 2b, where $A$ varies by district. Until then $P(X\in A)=1$ is assumed, which is the package's choice.

---

## 4. The sharp identified set

**Theorem 1 (characterisation of $\Theta_I$).** Assume (A1), (A5), standing integrability, and $H=Q(P)+K$ with $K\in\mathcal K$. The following are equivalent:
(i) $K\in\Theta_I(P;A)$;
(ii) (S3) holds, and there is a $P_V$-null set $\mathcal N$ such that for all $v=(x,n)\notin\mathcal N$ with $x\in A$,
$$\langle u,M_H(x)\rangle\le\psi(u;v)\quad\text{for all }u\in\mathbb R^{R\times C}.\tag{T1}$$
Moreover, $\Theta_I(P;A)$ is convex and closed, and it is compact when $A$ has nonempty relative interior.

*Interpretation:* the implied mean transition matrix at composition $x$ must be an average of selections from the polls' Duncan–Davis polytopes among polls with that composition and size, i.e. $M_H(x)\in\mathbb E[S(X,T)\mid V=v]$. Directions $u$ and $-u$ along the equality constraints force $x^\top M_H(x)=E[T\mid V=v]$, so (T1) also requires $E[T\mid V]=E[T\mid X]$, an observable implication of A3.

*Proof of (i)⇒(ii).*
1. Take $\mathbb Q\in\mathcal M(H;A)$ with $\mathbb Q\circ W^{-1}=P$. (S3) holds by membership. [Def $\Theta_I$; S3]
2. $\beta\in S(X,T)$ $\mathbb Q$-a.s. [S1; Def $S$]
3. For each fixed $u$, $\langle u,\beta\rangle\le h_S(u;X,T)$ a.s. [2; Def $h_S$]
4. Both sides are bounded by $\sqrt R\lVert u\rVert$, hence integrable. [Frob; L0(d)]
5. $E[\langle u,\beta\rangle\mid V]\le E[h_S(u;X,T)\mid V]$ a.s. [3–4; CE-mono]
6. The left side equals $\langle u,E[\beta\mid V]\rangle=\langle u,M_H(X)\rangle$ a.s. on $\{X\in A\}$. [CE-lin; S2]
7. The right side equals $\psi(u;V)$ a.s. The RCD of $T$ given $V$ under $\mathbb Q$ is a version of $P_{T\mid V}$, because $\mathbb Q\circ W^{-1}=P$. [RCD; L0(d)]
8. So for each $u$ there is a null set $\mathcal N_u$ outside which (T1) holds at that $u$. [5–7]
9. Let $U_0$ be a countable dense subset of $\mathbb R^{R\times C}$ and $\mathcal N=\bigcup_{u\in U_0}\mathcal N_u$. Then $\mathcal N$ is null. [countable union of null sets]
10. For $v\notin\mathcal N$, the left side of (T1) is linear in $u$, and $\lvert\psi(u;v)-\psi(u';v)\rvert\le\sqrt R\lVert u-u'\rVert$. [lin; L0(d) integrated]
11. So (T1) extends from $U_0$ to all $u$. [9–10; density and continuity] ∎

*Proof of (ii)⇒(i).* Work on the canonical space $(\mathcal W,\mathcal B(\mathcal W),P)$. Let $\mathcal S=\{\xi:\mathcal W\to\mathbb R^{R\times C}\text{ Borel}:\ \xi(w)\in S(x,t)\ P\text{-a.s.}\}$, viewed in $L^1$.

1. $\xi^0(w)_{rc}:=t_c$ is Borel, and $\xi^0\in\mathcal S$, so $\mathcal S\neq\emptyset$. [L0(a) proof]
2. $\mathcal S$ is convex. [L0(b)]
3. $\lVert\xi\rVert\le\sqrt R$ for every $\xi\in\mathcal S$, so $\mathcal S$ is bounded in $L^\infty$ and uniformly integrable. [$\lVert b\rVert\le\sqrt R$]
4. If $\xi_n\to\xi$ in $L^1$ with $\xi_n\in\mathcal S$, a subsequence converges a.s. Since $S(x,t)$ is closed, $\xi\in S$ a.s., so $\mathcal S$ is norm-closed. [Riesz; L0(b)]
5. $\mathcal S$ is weakly compact in $L^1$. [3; DP; 2, 4; Mazur]
6. $\Pi\xi:=E[\xi\mid V]$ is a bounded linear operator $L^1\to L^1(\mathcal G)$, hence weak-to-weak continuous. [CE contraction; wk-cont]
7. $\mathcal C:=\Pi(\mathcal S)$ is convex and weakly compact, hence weakly closed, in $L^1(\mathcal G)$. [5–6; lin]
8. Put $\eta:=M_H(X)1\{X\in A\}+\Pi(\xi^0)1\{X\notin A\}$. $\eta$ is $\mathcal G$-measurable and bounded, because $M_H$ is linear and $A$ is compact. [Def; alg]
9. Suppose, for contradiction, that $\eta\notin\mathcal C$. [—]
10. Then there exist $g\in L^\infty(\mathcal G;\mathbb R^{R\times C})$ and $\varepsilon>0$ with $E\langle g,\eta\rangle\ge\sup_{\zeta\in\mathcal C}E\langle g,\zeta\rangle+\varepsilon$. [HB on $(L^1(\mathcal G),\text{weak})$; 7]
11. For $\xi\in\mathcal S$: $E\langle g,\Pi\xi\rangle=E\big[E[\langle g,\xi\rangle\mid V]\big]=E\langle g,\xi\rangle$. [CE-pull, $g$ is $\mathcal G$-measurable; CE-tower]
12. $\sup_{\xi\in\mathcal S}E\langle g,\xi\rangle\le E[h_S(g(V);X,T)]$, pointwise. [Def $h_S$]
13. A Borel selector $\xi^\star(w)\in\arg\max_{b\in S(x,t)}\langle g(v),b\rangle$ exists, belongs to $\mathcal S$ and attains the bound in 12. So $\sup_{\mathcal S}E\langle g,\xi\rangle=E[h_S(g(V);X,T)]$. [MMT; L0(c); the objective is Carathéodory]
14. $E[h_S(g(V);X,T)]=E[\psi(g(V);V)]$. [RCD applied to $f(v,t)=h_S(g(v);x,t)$, bounded and Borel by L0(d); CE-tower]
15. On $\{X\in A\}\setminus\{V\in\mathcal N\}$: $\langle g(V),\eta\rangle=\langle g(V),M_H(X)\rangle\le\psi(g(V);V)$. [8; (T1) with $u=g(V(w))$, valid for all $u$ at once]
16. On $\{X\notin A\}$: $\langle g(V),\eta\rangle=E[\langle g(V),\xi^0\rangle\mid V]\le E[h_S(g(V);X,T)\mid V]=\psi(g(V);V)$ a.s. [8; CE-pull; $\xi^0\in S$; CE-mono; RCD]
17. Integrating 15–16: $E\langle g,\eta\rangle\le E\,\psi(g(V);V)=\sup_{\zeta\in\mathcal C}E\langle g,\zeta\rangle$, which contradicts 10. [15–16; CE-mono; 11, 13, 14]
18. Hence there is $\xi^\star\in\mathcal S$ with $E[\xi^\star\mid V]=\eta$ a.s. [9–17; Def $\mathcal C$]
19. Let $\beta:=\xi^\star(W)$ and $\mathbb Q:=$ the law of $(W,\beta)$ under $P$. Then $\mathbb Q\circ W^{-1}=P$. [Def]
20. $\beta\in\mathfrak B$ and $\beta^\top X=T$ a.s., which is (S1). [18; Def $\mathcal S$, $S$]
21. $E_{\mathbb Q}[\beta\mid V]=\eta=M_H(X)$ a.s. on $\{X\in A\}$, which is (S2). [18; 8]
22. (S3) holds by hypothesis (ii). So $\mathbb Q\in\mathcal M(H;A)$ and $K\in\Theta_I(P;A)$. [20–21; Def $\Theta_I$] ∎

*Proof of convexity, closedness, compactness.*
1. $K\mapsto M_{Q+K}(x)$ is affine, so the left side of (T1) is affine in $K$, and (S3) is a family of affine inequalities. [lin]
2. If $K,K'$ satisfy (ii) outside null sets $\mathcal N,\mathcal N'$, then every convex combination satisfies (ii) outside $\mathcal N\cup\mathcal N'$. [1]
3. If $K_n\to K$ with $K_n$ satisfying (ii) outside $\mathcal N_n$, then $K$ satisfies (ii) outside $\bigcup_n\mathcal N_n$, which is null. [continuity in $K$; countable union]
4. $\Theta_I\subseteq\Theta_O$, and $\Theta_O$ is bounded when $A$ has nonempty relative interior. [P1 below; L3(b)] ∎

**Remark 4.1 (no auxiliary randomisation).** P1 §1.5 said the construction needs extra randomisation. It does not: $S(x,t)$ is convex, so the selection can be chosen as a function of $W$ (step 18).

**Corollary 1 (sharp bounds for any linear target).** Assume Theorem 1 holds, $P(X\in A)=1$ and $\Theta_I(P;A)\neq\emptyset$. Then the sharp identified set for $\varphi=\sum a_{rc}B^*_{rc}$ is the interval
$$\Big[\min_{K\in\Theta_I}\varphi(K),\ \max_{K\in\Theta_I}\varphi(K)\Big],$$
where $\varphi(K)$ is given by Lemma 2.
*Proof.*
1. For every $\mathbb Q$ consistent with $P$ (i.e. in $\mathcal M(Q(P)+K;A)$ with $\mathbb Q_W=P$), $\varphi(B^*(\mathbb Q))=\varphi(K)$. [L2; lin]
2. So the identified set is the image $\varphi(\Theta_I)$, and every point of it is attained by some consistent $\mathbb Q$. [1; T1, (ii)⇒(i)]
3. $\varphi$ is affine and $\Theta_I$ is convex and compact, so the image is a compact interval with the stated endpoints. [lin; T1; IVT] ∎

---

## 5. The outer set: reduction, nesting, and when it is sharp

**Lemma 3.**
(a) If $A$ is a polytope and $\ell$ is linear, then $\ell\ge0$ on $A$ iff $\ell(v)\ge0$ for all $v\in\operatorname{vert}(A)$. Hence the constraint "$\forall x\in A$" in $\Theta_O$ can be replaced by "$\forall v\in\operatorname{vert}(A)$".
(b) If $A$ has nonempty relative interior, $\Theta_O(P;A)$ is bounded.
(c) Under (S1)–(S2), if $A=\operatorname{conv}(\operatorname{supp}P_X)$ then (S3) holds automatically.

*Proof.*
1. (a) "⇒" is immediate since $\operatorname{vert}(A)\subseteq A$. "⇐": write $x=\sum_k\lambda_kv_k$ with $\lambda\ge0$, $\sum\lambda_k=1$. [Mink]
2. Then $\ell(x)=\sum_k\lambda_k\ell(v_k)\ge0$. [lin] ∎(a)
3. (b) For $K\in\Theta_O$ and $x\in A$, every row of $M_{Q+K}(x)$ is nonnegative and sums to $e_r^\top Jx=1$, so it lies in $\Delta^{C-1}$. [Def $\Theta_O$; L1(b)]
4. Hence $\lvert e_r^\top K_cx\rvert\le1+\lvert e_r^\top Q_cx\rvert$ for all $x\in A$. [3; alg]
5. A relatively open subset of $\Delta^{R-1}$ contains $R$ affinely independent points $x^{(1)},\dots,x^{(R)}$. These are linearly independent because $0\notin\{\mathbf 1^\top x=1\}$. [alg]
6. $K_c=[\,e_r^\top K_cx^{(k)}\,]_{r,k}\,[x^{(1)}\cdots x^{(R)}]^{-1}$ is therefore bounded uniformly over $\Theta_O$. [4–5; lin] ∎(b)
7. (c) $E[\beta\mid V]\in\mathfrak B$ a.s., since $\mathfrak B$ is closed and convex. [CE-cvx]
8. So $M_H(x)\in\mathfrak B$ for $P_X$-a.e. $x$, and by continuity on $\operatorname{supp}P_X$. [S2; supp]
9. The rows of $M_H(x)$ are linear in $x$ with constant row sums, so $M_H(x)\in\mathfrak B$ on the convex hull $A$. [L1(b); lin; Mink] ∎(c)

**Proposition 1 (nesting).** $\Theta_I(P;A)\subseteq\Theta_X(P;A)\subseteq\Theta_O(P;A)$.
*Proof.*
1. Let $\mathbb Q\in\mathcal M(H;A)$ (conditioning on $V$). Then $E[\beta\mid X]=E[E[\beta\mid V]\mid X]=M_H(X)$ on $\{X\in A\}$, so $\mathbb Q$ belongs to the $X$-only model. [CE-tower; S2]
2. Every $K$ admitted by the stronger model is admitted by the weaker one, so $\Theta_I\subseteq\Theta_X$. [1; Def]
3. Let $K\in\Theta_X$ with witness $\mathbb Q$. Then $M_H(X)=E[\beta\mid X]\in\mathfrak B$ a.s. on $\{X\in A\}$. [S2 ($X$-version); CE-cvx]
4. (S3) gives $M_H(x)\in\mathfrak B$ on all of $A$; in particular $(Q_c+K_c)x\ge0$ for all $x\in A$. [S3; Def $\mathfrak B$]
5. So $K\in\Theta_O$. [4; Def $\Theta_O$] ∎

The JKST "Jensen" step (their proof of Prop. 2) is this inclusion, $\Theta_X\subseteq\Theta_O$. In R×C it needs no Jensen inequality: conditional expectations preserve the convex set $\mathfrak B$, and the accounting identity in means ($x^\top M_H(x)=\tau(x)$) holds by Lemma 1.

**Proposition 2 ($\Theta_O$ is sharp for the regression function).** Let $\mathcal P(P)$ be the set of laws $P'$ of $(N,X,T')$ with $P'_{(N,X)}=P_{(N,X)}$ and $E_{P'}[T'\mid V]=E_P[T\mid X]$ a.s. Assume $P(X\in A)=1$. Then
$$\Theta_O(P;A)=\bigcup_{P'\in\mathcal P(P)}\Theta_I(P';A).$$
*Proof.*
1. ⊇: For $P'\in\mathcal P(P)$, $Q(P')=Q(P)$, because the regression function and $\operatorname{supp}P_X$ are the same. [L1(a); Def $\mathcal P$]
2. $\Theta_O$ depends on $P$ only through $Q(P)$ and $A$, so $\Theta_O(P')=\Theta_O(P)$. [Def $\Theta_O$; 1]
3. $\Theta_I(P')\subseteq\Theta_O(P')=\Theta_O(P)$. [P1; 2]
4. ⊆: Let $K\in\Theta_O(P;A)$, $H=Q(P)+K$. Put $\beta:=M_H(X)$ and $T':=\beta^\top X$. [Def]
5. $\beta\in\mathfrak B$ a.s.: it is nonnegative on $A$ and its rows sum to 1. [Def $\Theta_O$; L1(b); $P(X\in A)=1$]
6. (S1) holds by construction of $T'$. [4–5]
7. $E[\beta\mid V]=M_H(X)$, because $\beta$ is $V$-measurable. This is (S2). [CE-pull]
8. (S3) holds by the definition of $\Theta_O$. [Def]
9. $T'_c=X^\top H_cX=X^\top Q_cX=E_P[T_c\mid X]$, so the law $P'$ of $(N,X,T')$ is in $\mathcal P(P)$. [L1(a) step 5; Def $Q(P)$]
10. The law of $(N,X,T',\beta)$ is in $\mathcal M(H;A)$ with $W$-marginal $P'$, so $K\in\Theta_I(P';A)$. [6–9; Def $\Theta_I$] ∎

*Meaning:* $\Theta_O$ is exactly what can be learned when only the regression of $T$ on $X$ is used. It is also the sharp set for any dataset in which polls with the same composition have identical outcomes.

**Corollary 2.** If $T=E[T\mid X]$ a.s., then $\Theta_I(P;A)=\Theta_O(P;A)$.
*Proof.*
1. $\Theta_I(P)\subseteq\Theta_O(P)$. [P1]
2. In the proof of P2, $T'=E_P[T\mid X]=T$ a.s., so $P'=P$. [P2 step 9; hypothesis]
3. Hence every $K\in\Theta_O(P)$ lies in $\Theta_I(P)$. [P2 step 10] ∎

**Proposition 3 (sufficient condition for $\Theta_I=\Theta_O$).** Assume $P(X\in A)=1$ and $E[T\mid V]=E[T\mid X]=:\tau(X)$ a.s. Suppose that for $P_V$-a.e. $v=(x,n)$ and every $u$, the map $t\mapsto h_S(u;x,t)$ is affine on $\operatorname{conv}\operatorname{supp}P_{T\mid V=v}$. Then $\Theta_I(P;A)=\Theta_O(P;A)$.
*Proof.*
1. Let $K\in\Theta_O$ and $H=Q+K$. Fix such a $v$. [Def]
2. $\psi(u;v)=\int h_S(u;x,t)P_{T\mid V=v}(dt)=h_S\big(u;x,\int t\,P_{T\mid V=v}(dt)\big)$: the integral of an affine function equals the function at the mean, and the mean lies in the closed convex hull of the support. [hypothesis; lin]
3. $=h_S(u;x,\tau(x))$. [RCD; hypothesis $E[T\mid V]=\tau(X)$]
4. $M_H(x)\in\mathfrak B$ and $x^\top M_H(x)=(x^\top H_cx)_c=\tau(x)$, so $M_H(x)\in S(x,\tau(x))$. [Def $\Theta_O$; L1(a)]
5. So $\langle u,M_H(x)\rangle\le h_S(u;x,\tau(x))=\psi(u;v)$ for all $u$. [Def $h_S$; 2–3]
6. Condition (ii) of Theorem 1 holds, with (S3) from the definition of $\Theta_O$. So $K\in\Theta_I$. [T1]
7. The reverse inclusion is Proposition 1. [P1] ∎

**Remark 3.1 (2×2 reading of Prop. 3).** In 2×2, $S(x,t)$ is the segment with $\beta^w\in[L,U]$, where $L=\max\{0,(t_1-x)/(1-x)\}$ and $U=\min\{1,t_1/(1-x)\}$.
1. The segment's direction depends on $x$ only. So for fixed $u$ the maximising endpoint is the same for every $t$, and $h_S(u;x,\cdot)$ is an affine image of $L$ or of $U$. [alg]
2. $L$ and $U$ are piecewise affine in $t_1$, with kinks at $t_1=x$ and $t_1=1-x$ respectively. [alg]
3. So Prop. 3 applies when, among polls with the same composition, the vote share never lies on both sides of $x$ or of $1-x$. [1–2]

In words: $\Theta_O$ loses nothing unless polls with similar composition are so different that some sit at a Duncan–Davis corner and others do not. Verified numerically: inclusion is strict in 46 of 400 random valid designs, and never reversed.

---

## 6. Closed forms

### 6.1 Equivalence with JKST (2×2)

**Proposition 4.** Let $R=C=2$, with row 1 = black, $X=X_1$, and $A=\{(t,1-t):t\in[l,u]\}$, $0<l<u<1$. Write $E[T\mid X]=w_0+c_1X+d_1X^2$ and $K:=K_{1,12}$. Then
(a) $Q_{1,22}=w_0$, $Q_{1,12}=w_0+c_1/2$, $Q_{1,11}=w_0+c_1+d_1$, and $K=c_1/2-w_1$;
(b) $\Theta_O$, written in terms of $w_1$, is exactly JKST Prop. 2;
(c) with $A=\Delta^1$ it is exactly JKST Prop. 3;
(d) $B^*(K)$ equals JKST (16).

*Derivation (a).*
1. $X^\top Q_1X=Q_{11}t^2+2Q_{12}t(1-t)+Q_{22}(1-t)^2$ with $t=X$. [Def]
2. $=Q_{22}+(2Q_{12}-2Q_{22})t+(Q_{11}-2Q_{12}+Q_{22})t^2$. [alg]
3. Matching coefficients with $w_0+c_1t+d_1t^2$: $Q_{22}=w_0$; $Q_{12}=w_0+c_1/2$; $Q_{11}=d_1+2Q_{12}-Q_{22}=w_0+c_1+d_1$. [L1(a) uniqueness; alg]
4. White row: $E[\beta^w\mid X]=H_{21}t+H_{22}(1-t)=H_{22}+(H_{21}-H_{22})t$, so $w_0=H_{22}$ and $w_1=H_{21}-H_{22}$. [Def]
5. $H_{21}=Q_{12}-K$ and $H_{22}=Q_{22}$, so $w_1=w_0+c_1/2-K-w_0$, i.e. $K=c_1/2-w_1$. [skew; 3–4]

*Derivation (b).* Column 2 is $H_2=J-H_1$ (L1(b)). At a vertex $v=(t,1-t)$:
1. Row 2, column 1: $H_{21}t+H_{22}(1-t)=w_0+w_1t\ge0\iff w_1\ge-w_0/t$. [4 above; alg]
2. Row 2, column 2: $1-(w_0+w_1t)\ge0\iff w_1\le(1-w_0)/t$. [L1(b); alg]
3. Row 1, column 1: $H_{11}=Q_{11}$ and $H_{12}=Q_{12}+K=w_0+c_1-w_1$. [(a); skew]
4. So $H_{11}t+H_{12}(1-t)=(w_0+c_1+d_1)t+(w_0+c_1-w_1)(1-t)=w_0+c_1+d_1t-w_1(1-t)$. [alg]
5. $\ge0\iff w_1\le(w_0+c_1+d_1t)/(1-t)$. [alg; $t<1$]
6. Row 1, column 2: $1-[\,4\,]\ge0\iff w_1\ge(w_0+c_1+d_1t-1)/(1-t)$. [L1(b); alg]
7. JKST write the last bound as $(w_0+c_1+d_1-1)/(1-x)-d_1=(w_0+c_1+d_1x-1)/(1-x)$, which is the same. [alg]
8. Taking $\max$ / $\min$ over $t\in\{l,u\}$ gives their $w_l,w_u$. [L3(a); 1–7] ∎(b)

This is also verified numerically against `generateBounds`: agreement to $10^{-16}$ (P1 §0).

*Derivation (c).* With $\operatorname{vert}(A)=\{(1,0),(0,1)\}$:
1. At $(1,0)$: rows 1 and 2 give $Q_{11}\in[0,1]$, which does not involve $K$, and $w_0+w_1\in[0,1]$. [Def; 1–2 of (b) at $t=1$]
2. At $(0,1)$: $H_{22}=w_0\in[0,1]$, which does not involve $K$, and $H_{12}=w_0+c_1-w_1\in[0,1]$. [Def; (b)3]
3. Together: $\max\{-w_0,\,c_1+w_0-1\}\le w_1\le\min\{1-w_0,\,c_1+w_0\}$, which is JKST Prop. 3. [alg] ∎(c)

The same result also follows from Prop. 5 below: verified symbolically in `verify_phase2a.py`.

*Derivation (d).*
1. $\mu_1B^*=Q_{11}\mu_{11}+(Q_{12}+K)\mu_{12}=(w_0+c_1+d_1)\mu_{11}+(w_0+c_1-w_1)\mu_{12}$. [L2; (a)]
2. $=(w_0+c_1)(\mu_{11}+\mu_{12})+d_1\mu_{11}-w_1\mu_{12}$. [alg]
3. $\mu_{11}+\mu_{12}=E[NX(X+1-X)]=\mu_1$. [lin]
4. So $B^*=w_0+c_1+d_1\mu_{11}/\mu_1-w_1\mu_{12}/\mu_1$, the population version of JKST (16). [1–3] ∎

### 6.2 $A$ = the whole simplex (general R×C)

**Proposition 5.** With $A=\Delta^{R-1}$:
$$\Theta_O=\Big\{K\in\mathcal K:\ \lvert K_{c,rs}\rvert\le Q_{c,rs}\ \ \forall c,\ r<s\Big\},$$
which also requires $Q_{c,rr}\ge0$ (a testable condition not involving $K$). The set separates by row pair $(r,s)$, and for each pair and column
$$K_{c,rs}\in[-m_{c,rs},\,m_{c,rs}],\qquad m_{c,rs}:=\min\{Q_{c,rs},\,1-Q_{c,rs}\}.$$
Hence
$$B^*_{rc}\in\frac1{\mu_r}\Big[\sum_sQ_{c,rs}\mu_{rs}\ \mp\ \sum_{s\neq r}m_{c,rs}\mu_{rs}\Big],$$
with the minus sign giving the lower bound.

*Derivation.*
1. $\operatorname{vert}(\Delta^{R-1})=\{e_1,\dots,e_R\}$, so $\Theta_O=\{K:(Q_c+K_c)e_s\ge0\ \forall c,s\}=\{K:Q_{c,rs}+K_{c,rs}\ge0\ \forall c,r,s\}$. [L3(a); alg]
2. For $r=s$: $K_{c,rr}=0$, so the condition is $Q_{c,rr}\ge0$. [skew]
3. For $r\ne s$, the pair of conditions $(r,s)$ and $(s,r)$ reads $Q_{c,rs}+K_{c,rs}\ge0$ and $Q_{c,rs}-K_{c,rs}\ge0$, i.e. $\lvert K_{c,rs}\rvert\le Q_{c,rs}$. [$Q$ symmetric; $K$ skew]
4. The only other restriction is $\sum_cK_{c,rs}=0$, which involves one pair at a time, so the set is a product over pairs. [Def $\mathcal K$; 3]
5. For a fixed pair, apply Lemma 4(b) with $\ell_c=-Q_{c,rs}$, $h_c=Q_{c,rs}$: the coordinate range is $[\max\{-Q_{c,rs},-\sum_{c'\ne c}Q_{c',rs}\},\ \min\{Q_{c,rs},\sum_{c'\ne c}Q_{c',rs}\}]$. [L4(b)]
6. $\sum_{c'}Q_{c',rs}=J_{rs}=1$, so $\sum_{c'\ne c}Q_{c',rs}=1-Q_{c,rs}$, and the range is $[-m_{c,rs},m_{c,rs}]$. [L1(b); alg]
7. $B^*_{rc}$ depends on $\{K_{c,rs}\}_{s\ne r}$, which belong to different pairs and so vary independently, each with coefficient $\mu_{rs}/\mu_r\ge0$. [L2; 4]
8. Its range is therefore attained coordinate-wise at $\mp m_{c,rs}$. [6–7; lin] ∎

*Meaning:* $H_{c,rs}$ and $H_{c,sr}$ ("$r$-voters in pure-$s$ polls" and "$s$-voters in pure-$r$ polls") are two probabilities whose average $Q_{c,rs}$ is identified. Splitting a known average of two numbers in $[0,1]$ leaves exactly $\pm\min(Q,1-Q)$ of freedom. This holds only if the linear model extrapolates validly to pure polls (A4 with $A=\Delta$). Verified against an LP solver in 200 random tables.

### 6.3 Box ∩ hyperplane (used by Props 5–6)

**Lemma 4.** Let $\ell\le h$ in $\mathbb R^C$ and $\mathcal P=\{\kappa:\ell\le\kappa\le h,\ \sum_c\kappa_c=0\}$.
(a) $\mathcal P\neq\emptyset\iff\sum\ell_c\le0\le\sum h_c$.
(b) $\{\kappa_c:\kappa\in\mathcal P\}=\big[\max\{\ell_c,-\sum_{c'\ne c}h_{c'}\},\ \min\{h_c,-\sum_{c'\ne c}\ell_{c'}\}\big]$.
(c) $\max_{\mathcal P}g^\top\kappa$ is attained by the greedy rule: start from $\kappa=\ell$, and raise coordinates in decreasing order of $g_c$, each up to $h_c$, until $\sum\kappa_c=0$.

*Proof.*
1. The box $[\ell,h]$ is convex, so its image under $\kappa\mapsto\sum\kappa_c$ is an interval containing $\sum\ell$ and $\sum h$, and contained in $[\sum\ell,\sum h]$ by monotonicity. So it equals $[\sum\ell,\sum h]$. [IVT; lin]
2. (a) follows: $0$ is attained iff $0\in[\sum\ell,\sum h]$. [1]
3. (b): $\kappa_c=k$ extends to a point of $\mathcal P$ iff $\ell_c\le k\le h_c$ and the box in the other coordinates contains a point with sum $-k$. [Def]
4. By 1 applied to the other coordinates, that second condition is $\sum_{c'\ne c}\ell_{c'}\le-k\le\sum_{c'\ne c}h_{c'}$. Intersecting gives (b). [1; alg]
5. (c): substitute $\kappa=\ell+y$ with $0\le y\le h-\ell$ and $\sum y=s:=-\sum\ell\ge0$. The objective becomes $g^\top\ell+g^\top y$. [alg; (a)]
6. If $y$ is optimal, $g_a>g_b$, $y_a<h_a-\ell_a$ and $y_b>0$, moving $\delta>0$ from $b$ to $a$ stays feasible and raises the objective by $\delta(g_a-g_b)>0$, a contradiction. So the optimum fills coordinates in decreasing order of $g$, which is the greedy rule. [exch] ∎

### 6.4 Two rows, any number of columns, $A=[l,u]$

**Proposition 6.** Let $R=2$ and $\operatorname{vert}(A)=\{(t_k,1-t_k)\}_{k=1,2}$ with $0<t_1=l<t_2=u<1$. Write $\kappa_c:=K_{c,12}$, so that $K_c=\kappa_c\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$. Define
$$\lambda_c(t)=-Q_{c,12}-Q_{c,11}\tfrac{t}{1-t},\qquad\rho_c(t)=Q_{c,12}+Q_{c,22}\tfrac{1-t}{t},$$
$$\ell_c=\max_k\lambda_c(t_k),\qquad h_c=\min_k\rho_c(t_k).$$
Then $\Theta_O=\{\kappa:\ell\le\kappa\le h,\ \sum_c\kappa_c=0\}$. The projected interval for each $\kappa_c$ is given by Lemma 4(b). The bounds are
* $B^*_{1c}=[Q_{c,11}\mu_{11}+(Q_{c,12}+\kappa_c)\mu_{12}]/\mu_1$, attained at the endpoints of $\kappa_c$'s projected interval (increasing in $\kappa_c$);
* $B^*_{2c}=[(Q_{c,12}-\kappa_c)\mu_{12}+Q_{c,22}\mu_{22}]/\mu_2$ (decreasing in $\kappa_c$);
* any $\varphi=\varphi_0+\sum_cg_c\kappa_c$, with $g_c=\mu_{12}(a_{1c}/\mu_1-a_{2c}/\mu_2)$, by the greedy rule of Lemma 4(c).

*Derivation.*
1. Row 1 at $v=(t,1-t)$: $e_1^\top(Q_c+K_c)v=Q_{c,11}t+(Q_{c,12}+\kappa_c)(1-t)\ge0$. [Def]
2. $\iff\kappa_c\ge-[Q_{c,11}t+Q_{c,12}(1-t)]/(1-t)=\lambda_c(t)$. [alg; $1-t>0$]
3. Row 2: $(Q_{c,12}-\kappa_c)t+Q_{c,22}(1-t)\ge0\iff\kappa_c\le\rho_c(t)$. [alg; $t>0$]
4. Over both vertices: $\ell_c\le\kappa_c\le h_c$. Adding $\sum_c\kappa_c=0$ gives the stated set. [L3(a); Def $\mathcal K$]
5. The $B^*$ formulas are Lemma 2 with $R=2$. Their monotonicity follows from $\mu_{12}\ge0$. [L2]
6. $\varphi$ is affine in $\kappa$ with the stated $g$, and Lemma 4(c) applies. [L2; lin; L4(c)] ∎

Verified against an LP solver in 300 random cases, for both the projection and the greedy rule.

**Corollary 6.1 (joint vs collapsed; corrects P1 §1.4).** The collapsed analysis of column $c$ treats the table as two columns $(c,\neg c)$ with $Q_{\neg c}=J-Q_c$. For $R=2$:
(a) the joint interval for $\kappa_c$ is contained in the collapsed one;
(b) the two are equal whenever, for each $r\in\{1,2\}$, the numbers $Q_{c',rr}$ ($c'\ne c$) are all of one sign. That holds in particular when the extrapolated pure-poll rates $H_{c',rr}=Q_{c',rr}$ are nonnegative;
(c) strict tightening therefore needs strong contextual effects, under which some pure-poll extrapolation is negative.

For $R\ge3$, rows are coupled through the shared $K_{c,rs}$, and the joint bound can be strictly tighter even when every extrapolation is valid. A characterisation for $R\ge3$ is left open.

*Proof.*
1. $\lambda_{\neg c}(t)=-(1-Q_{c,12})-(1-Q_{c,11})\tfrac t{1-t}=\sum_{c'\ne c}\lambda_{c'}(t)$, because $\sum_{c'}Q_{c'}=J$. [Def; L1(b)]
2. The collapsed upper end is $\min\{h_c,-\ell_{\neg c}\}$, and $-\ell_{\neg c}=-\max_k\sum_{c'\ne c}\lambda_{c'}(t_k)$. [L4(b) with two columns; 1]
3. $\max_k\sum_{c'}\lambda_{c'}(t_k)\le\sum_{c'}\max_k\lambda_{c'}(t_k)=\sum_{c'\ne c}\ell_{c'}$. [the max of a sum is at most the sum of the maxima]
4. So the collapsed upper end is at least $\min\{h_c,-\sum_{c'\ne c}\ell_{c'}\}$, the joint upper end. The lower ends compare the same way, using $\rho$. That proves (a). [2–3; L4(b)]
5. $t\mapsto t/(1-t)$ is increasing, so $\lambda_{c'}$ is monotone with direction set by the sign of $Q_{c',11}$. Likewise $(1-t)/t$ is decreasing and $\rho_{c'}$ is monotone with direction set by the sign of $Q_{c',22}$. [alg]
6. If the signs agree across $c'\ne c$, all $\lambda_{c'}$ are maximised at the same vertex, so equality holds in 3. The same argument applies to $\rho$. That proves (b). [5]
7. (c) is the contrapositive of (b), with $H_{c',rr}=Q_{c',rr}$. [skew] ∎

Numerical evidence (`verify_phase2a.py`):
* 2×3 with nonnegative pure-poll rates: never strict, as (b) predicts.
* 2×3 with the model valid only on $A$: strict in 46 cases.
* 3×3 with all rates valid: strict in 20 cases.
* Never looser in any case.

Joint fitting remains necessary anyway for targets that combine columns (e.g. loyalty minus defection) and for $R\ge3$.

### 6.5 General R×C: the linear program

**Proposition 7.** Parametrise $\mathcal K$ by $k=(k_{c,rs})_{c<C,\ r<s}\in\mathbb R^d$, $d=(C-1)R(R-1)/2$, with $K_c=\sum_{r<s}k_{c,rs}(e_re_s^\top-e_se_r^\top)$ for $c<C$ and $K_C=-\sum_{c<C}K_c$. For a polytope $A$:
$$\Theta_O=\{k:Gk\le g\},$$
with one row for each $(v,r,c)\in\operatorname{vert}(A)\times[R]\times[C]$:
$$-\big(e_r^\top K_c(k)v\big)\le(Q_cv)_r,\qquad e_r^\top K_c(k)v=\sum_{s>r}k_{c,rs}v_s-\sum_{s<r}k_{c,sr}v_s\ \ (c<C),\quad e_r^\top K_C(k)v=-\sum_{c<C}e_r^\top K_c(k)v.$$
For a linear target, $\varphi(k)=\varphi_0+\gamma^\top k$ with
$$\varphi_0=\sum_{rc}\frac{a_{rc}}{\mu_r}\sum_sQ_{c,rs}\mu_{rs},$$
and, for $c<C$ and $r<s$,
$$\gamma_{c,rs}=\mu_{rs}\Big[\frac{a_{rc}-a_{rC}}{\mu_r}-\frac{a_{sc}-a_{sC}}{\mu_s}\Big].$$
This follows from Lemma 2: $k_{c,rs}$ enters $B^*_{rc}$ with $+\mu_{rs}/\mu_r$ and $B^*_{sc}$ with $-\mu_{rs}/\mu_s$ (since $K_{c,sr}=-k_{c,rs}$), and enters column $C$ with the opposite signs (since $K_C=-\sum_{c<C}K_c$). The bounds are $\varphi_0+\min/\max\{\gamma^\top k:Gk\le g\}$. With $A$ the sample hull of $\{X_i\}$, imposing the rows at every observed $X_i$ is equivalent to imposing them at the hull vertices.

*Proof.*
1. The parametrisation is a bijection onto $\mathcal K$. [L1(c) step 2]
2. $(K_cv)_r=\sum_sK_{c,rs}v_s$, with $K_{c,rs}=k_{c,rs}$ for $r<s$, $K_{c,rs}=-k_{c,sr}$ for $r>s$, and $K_{c,rr}=0$. [Def; skew]
3. Each constraint $(Q_c+K_c)v\ge0$ becomes one row of $Gk\le g$. [L3(a); 2]
4. $\varphi$ is affine in $K$ by Lemma 2, hence affine in $k$ by 1. [L2; lin]
5. The vertices of $\operatorname{conv}\{X_i\}$ are among the $X_i$, and every $X_i\in A$. So constraints at all $X_i$ imply constraints at the vertices, and conversely by Lemma 3(a). [Mink; L3(a)] ∎

The value function $Q\mapsto\max\{\gamma^\top k:Gk\le g(Q)\}$ is concave and piecewise linear in $(Q,\mu)$, and smooth wherever the optimal basis is unique. That generalises JKST's "tie-breaking" condition and is the starting point of Phase 3.

---

## 7. Refutability and undetectable violations

**Proposition 8 (non-emptiness; generalises JKST App. B).** Suppose $P$ is generated by some $\mathbb Q\in\mathcal M(H_0;A)$ with $P(X\in A)=1$. Then $K_0\in\Theta_I\subseteq\Theta_O$, and $B^*_{rc}\in[\varphi_L,\varphi_U]\cap DD^*_{rc}$, where
$$DD^*_{rc}=\Big[\frac{E[N\max\{0,X_r+T_c-1\}]}{\mu_r},\ \frac{E[N\min\{X_r,T_c\}]}{\mu_r}\Big].$$
Hence an empty $\Theta_O$, or an empty intersection with $DD^*$, refutes (A1)–(A5) jointly.
*Proof.*
1. $K_0\in\Theta_I$ by definition, and $\Theta_I\subseteq\Theta_O$. [Def; P1]
2. $B^*_{rc}=\varphi(K_0)\in[\varphi_L,\varphi_U]$. [L2; Cor 1; 1]
3. Write $n_{rc}=X_r\beta_{rc}$. Then $\sum_cn_{rc}=X_r$, $\sum_rn_{rc}=T_c$ and $n\ge0$. [S1]
4. $n_{rc}\le\min\{X_r,T_c\}$, and $n_{rc}=X_r-\sum_{c'\ne c}n_{rc'}\ge X_r-\sum_{c'\ne c}T_{c'}=X_r+T_c-1$. [3; FH]
5. Multiply by $N\ge0$, take expectations, and divide by $\mu_r$. [CE-mono; standing integrability] ∎

In samples an empty set can also arise from sampling noise; Phase 3 turns this into a test.

**Proposition 9 (what the data can never rule out).** Let $\Delta:\Delta^{R-1}\to\mathbb R^{R\times C}$ be Borel with (i) $x^\top\Delta(x)=0$, (ii) $\Delta(x)\mathbf 1=0$, and (iii) $M_H(x)+\Delta(x)\in\mathfrak B$ on $\operatorname{supp}P_X$. Then the structure $\beta=M_H(X)+\Delta(X)$ and the structure $\beta'=M_H(X)$ generate the same $P$, but
$$B^*_{rc}(\beta)-B^*_{rc}(\beta')=\frac{E[NX_r\Delta_{rc}(X)]}{\mu_r}.$$
The linear $\Delta$'s are exactly $M_{K'}$ with $K'\in\mathcal K$, i.e. the unidentified part. A nonlinear $\Delta$ violates A2 yet leaves no trace in the data.

*Proof.*
1. $T=\beta^\top X=M_H(X)^\top X+\Delta(X)^\top X=\tau(X)+0$, and likewise $T'=\tau(X)$. With the same law of $(N,X)$, both give the same $P$. [S1; (i); L1(a)]
2. Both are valid structures. [(ii)–(iii); $M_H$ has unit row sums by L1(b)]
3. The difference in $B^*$ follows by the computation of Lemma 2. [L2 steps 1–4]
4. If $\Delta=M_D$ is linear, then (i) gives $x^\top D_cx=0$, so $\operatorname{sym}D_c=0$ by A5, and (ii) gives $\sum_cD_c=0$, so $D\in\mathcal K$. Conversely every $M_{K'}$ with $K'\in\mathcal K$ satisfies (i)–(ii). [L1(a) step 5; A5; L1(b)] ∎

*Example ($R=2$).* Take $\Delta_{1c}=\delta_c(x)x_2$ and $\Delta_{2c}=-\delta_c(x)x_1$ with $\sum_c\delta_c=0$. A constant $\delta_c$ is a shift of $\kappa_c$, invisible but inside the model. $\delta_c(x)=b_cx_1$ instead gives a *quadratic* contextual effect for row 1 that is exactly cancelled in row 2. It shifts $B^*_{1c}$ by $b_cE[NX_1^2X_2]/\mu_1$ and can push the truth outside $\Theta_O$; JKST Example 6 is the case $C=2$. This is the precise sense in which A2 is an untestable leap. Phase 3 bounds its consequences through a sensitivity parameter $\sup_x\lVert\Delta(x)\rVert$.

---

## 8. The BC 3×3 case

Rows and columns are both ordered $(P,O,N)$: in-party, out-party, non-voter, both last time and this time. There are $d=6$ unknowns:
$$k_{c,PO},\ k_{c,PN},\ k_{c,ON}\quad\text{for }c\in\{P,O\},$$
with the abstain column fixed by $K_N=-K_P-K_O$.

*Reading one unknown.* $k_{P,PO}=\tfrac12(H_{P,PO}-H_{P,OP})$ is half the gap between two hypothetical rates:
* the P-loyalty of P-voters in a poll otherwise all out-party; and
* the rate at which out-party voters switch to P in a poll otherwise all in-party.

The data identify the average of the two rates, not the gap.

*One constraint* (poll $x=(x_P,x_O,x_N)$, row $P$, column $P$):
$$Q_{P,PP}x_P+Q_{P,PO}x_O+Q_{P,PN}x_N+k_{P,PO}x_O+k_{P,PN}x_N\ \ge\ 0.$$
There are 9 such rows per poll vertex.

*Loyalty:*
$$B^*_{PP}=\frac{Q_{P,PP}\mu_{PP}+(Q_{P,PO}+k_{P,PO})\mu_{PO}+(Q_{P,PN}+k_{P,PN})\mu_{PN}}{\mu_P}.$$

**Illustration** (noise-free, identification only; arbitrary linear DGP; 2000 polls; $A$ = their hull; `verify_phase2a.py`):

| Target | Truth | $\Theta_O$ bound | DD bound | Width ratio |
|---|---|---|---|---|
| Loyalty P→P | 0.848 | [0.396, 0.932] | [0.131, 0.973] | 0.64 |
| Defection P→O | 0.086 | [0.022, 0.361] | [0.000, 0.628] | 0.54 |
| Demobilisation P→abstain | 0.066 | [0.022, 0.301] | [0.000, 0.532] | 0.52 |
| Mobilisation N→P | 0.169 | [0.048, 0.590] | [0.000, 0.896] | 0.60 |

All four contain the truth, and all are narrower than DD. They are still wide, because six numbers are unidentified. Narrowing them is what Phase 2b (riding dummies) and the optional shape restrictions (P1 §2) are for.

---

## 9. Verification log

`python3 -I verify_phase2a.py` (22 PASS, 0 FAIL; full output in `verify_phase2a_output.txt`) checks:
* Lemma 1: the dimension $(C-1)R(R-1)/2$ by rank computation for 2×2, 2×3, 3×3 and 4×3; $x^\top Kx=0$ and uniqueness of $Q$ symbolically for $R=3$.
* Prop. 4: the JKST Prop. 3 constraints symbolically; Prop. 2 numerically (P1 script).
* Prop. 5 against an LP solver (200 tables).
* Prop. 6 projection and greedy rule against an LP solver (300 cases).
* Cor. 6.1: inclusion, the equality case, and strict cases.
* Prop. 1: $\Theta_X\subseteq\Theta_O$ in 400 random valid 2×2 designs, strict in 46.
* Theorem 1 "⇒": a known feasible structure lies in $\Theta_X$ (227/227).
* Prop. 7: the closed-form $\gamma$ against direct computation (100 targets).
* The 3×3 example: the truth lies inside the bounds.

Not verified by computer: the measure-theoretic steps of Theorem 1, which rest on the cited theorems, and the theorem numbers marked †.

## 10. Corrections to Phase 1

1. **Joint vs collapsed (P1 §1.4).** For $R=2$, the joint bound is strictly tighter only under strong contextual effects (Cor. 6.1). The general claim of strict gain is withdrawn.
2. **A4.** The model now includes (S3). It is implied when $A$ is the hull of the data (Lemma 3(c)), and it is a real assumption only when extrapolating, as in JKST Prop. 3 and our Prop. 5.
3. **Theorem 1** needs no auxiliary randomisation (Remark 4.1).
4. **New result (Prop. 2):** $\Theta_O$ is exactly the sharp set for an analyst who uses only the regression function.

## 11. Open items carried forward

* **2b:** riding dummies (M1/M2), including $A\subsetneq\operatorname{supp}$ per district (Remark 2.1), and the 3×3 counterpart of P1 §7.3.
* **3:** inference for LP value functions (Prop. 7) with non-unique bases; the finite-population term $B-B^*$; the many-ridings intersection; sensitivity through $\sup\lVert\Delta\rVert$ (Prop. 9).
* **Open question:** characterise strict joint-vs-collapsed gains for $R\ge3$.

---

## 12. Relation to Cho & Manski (2008) and Cross & Manski (2002)

JKST's conclusion suggests extending to R×C tables "perhaps generalizing Cho and Manski 2008". That chapter surveys EI and summarises Cross & Manski (2002, *Econometrica* 70:357–368). The correspondence with this document:

| Cross–Manski / Cho–Manski | Here |
|---|---|
| outcome $Y$, group $Z$, precinct $X$ | column $c$, row $r$, unit $i$ |
| data: $P(Y\mid X)$, $P(Z\mid X)$ | $T_i$, $X_i$ |
| target: $P(Y\mid X,Z)$ | $\beta_i$; aggregated, $B$ |
| "data alone": feasible $P(Y\mid X=x,Z)$ = solutions of the finite-mixture problem $\sum_z\eta_zP(Z=z\mid x)=P(Y\mid x)$ | the unit DD polytope $S(x,t)$ (Lemma 0) |
| that region is convex with finitely many extreme points | $S(x,t)$ is a transportation polytope |
| region across precincts is a Cartesian product (no cross-precinct restriction) | DD for $B$ is a weighted sum of unit polytopes; P1 (D0) |
| identifying assumption: $V$ (part of the precinct identifier) as an instrument, $E(Y\mid V,W,Z)=E(Y\mid W,Z)$; Goodman's constancy is the case with no $W$ | the special case of LCE with **zero contextual effects**, $H_c=\alpha_c\mathbf 1^\top$. Then $K$ is pinned down and $B$ is point-identified |
| sharp region $H^*_w$ under that assumption; outer region $C^*_w$ from iterated expectations only, with $H^*_w\subseteq C^*_w$ | $\Theta_I\subseteq\Theta_O$; $\Theta_O$ is exactly what iterated expectations of the accounting identity deliver (Props 1–2) |
| rank of $\Pi$ = $\lvert Z\rvert$ gives point identification | A5 (rank of the quadratic design) |
| an empty region refutes the assumption | Prop. 8 |
| "layer assumptions in decreasing order of plausibility" | DD ⊇ LCE bounds ⊇ Goodman point (and P1 §6, M2 ⊇ M1 ⊇ within-riding constancy) |

**What is new relative to Cross–Manski:**
1. **The assumption is weaker than constancy.** Rates may vary linearly with composition. This leaves the $(C-1)R(R-1)/2$-dimensional $K$ unidentified, so the result is a set rather than a point. JKST's 2×2 bounds and our R×C bounds sit *between* their data-alone region and their Goodman point.
2. **It is a mean restriction on a superpopulation.** Cross–Manski treat each precinct's population as fixed, so constancy means $\beta_i=\beta_j$ exactly. LCE restricts only $E[\beta\mid X]$, and each poll may deviate. This is why the bounds are for aggregates only, why individual polls stay at DD, and why the finite-population gap $B-B^*$ must be handled in Phase 3.
3. **Composition is continuous.** Their instrument $V$ is a discrete precinct attribute with several precincts per cell. Here the "instrument" is the composition $X$ itself, entering through the linear model.

**Bearing on Phase 2b.** Their $W$-conditioned assumption (constancy across precincts *within* cells of $W$) is, with $W$ = riding, "Goodman within each riding". It is a stronger layer than M1, which lets rates vary linearly with composition inside each riding. It is a natural extra layer to report, and its point estimate offers a check: does it fall inside the M1 bounds?
