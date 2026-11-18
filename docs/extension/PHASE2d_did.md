# Phase 2d: difference-in-differences with partial identification, derivations and proofs

*Builds on P1, 2a, 2b, P3, 2c and the Phase-1 proposal `PHASE1b_monotone_did.md` (1b; decisions in 1b §D). Every proof step carries a justification tag: 2a §0.3, 2b §0.3, P3 §0.3, 2c §0.3, plus §0.3 below. Checks:*
* *`verify_phase2d.py` (output `verify_phase2d_output.txt`, 15/15 PASS);*
* *`verify_phase1_mei_did.py` V7–V8.*

*Application in view: Wisconsin 1857/1865 and Iowa 1857/1868 suffrage referenda. Townships are units; rows are veterans (v) and non-veterans (n); columns are yes, no, abstain on a common base of eligible men (1b §D.7).*

---

## 0. Setting

### 0.1 Fixed population (§§2–5)

* Units $i$, periods $t\in\{0,1\}$.
* The same men are present in both periods, and $N_i$ and $X_i$ (veteran share) are common to both. §§6–7 drop this.
* $\beta_{it}\in S(X_i,T_{it})$, and $\Delta\beta_i:=\beta_{i1}-\beta_{i0}$.
* Weights: $p_i$, $w^r_i$, $\mu_r$ as in 2c §0.1. Second moments: $\mu_{vv}=E_p[X^2]$, $\mu_{vn}=E_p[X(1-X)]$, $\mu_{nn}=E_p[(1-X)^2]$.
* Targets for column $c$:
  * the actual DiD, $\text{DiD}_c:=\sum_iw^v_i\Delta\beta_{ivc}-\sum_iw^n_i\Delta\beta_{inc}$;
  * the model-implied DiD $\widetilde{\text{DiD}}_c$, the same with $m$ in place of $\beta$ (P3 §0.2).

### 0.2 Homogeneous form for two rows

* $H_{c,rs}$ is the rate of $r$-men in a hypothetical unit made up entirely of $s$-men (P1 §7.2).
* $\kappa_c:=K_{c,vn}=\tfrac12(H_{c,vn}-H_{c,nv})$.
* The contextual slopes (P1 §1.3 dictionary) are:
  * non-veterans, $w_{1c}=H_{c,nv}-H_{c,nn}$;
  * veterans, $b_{1c}=H_{c,vv}-H_{c,vn}$;
  * curvature, $d_{1c}=b_{1c}-w_{1c}=Q_{c,vv}-2Q_{c,vn}+Q_{c,nn}$.

### 0.3 Additional tags

| Tag | Meaning |
|---|---|
| [prod-supp] | the support function of a Cartesian product is the sum of the support functions: $h_{S\times S'}((u,u'))=h_S(u)+h_{S'}(u')$ |
| [KRN] | Kuratowski–Ryll-Nardzewski: a measurable correspondence with nonempty closed values in a Polish space has a Borel selector (Aliprantis & Border 2006, Thm 18.13†) |
| [CC] | Charnes–Cooper (1962): a linear-fractional program over a polytope with a positive denominator is equivalent to an LP |
| [coupling] | two tables sharing a margin can be joined into a three-way table with both as two-way margins, e.g. by conditional independence given the shared margin |
| [bilin] | a bilinear function on a product of polytopes attains its extremes at a pair of vertices |
| [MTR], [MTS] | monotone treatment response (Manski 1997) and monotone treatment selection (Manski & Pepper 2000) as assumptions on counterfactual means |

---

## 1. Lemma G: 2a Theorem 1 for general random polytopes

**Lemma G.** Let $\mathcal S:\mathcal W\rightrightarrows\mathbb R^d$ be a correspondence with nonempty, compact, convex values, a closed graph and a uniform bound. Let $\theta\mapsto m_\theta(v)$ be affine with $m_\theta$ bounded and Borel. Consider the model "$\xi\in\mathcal S(W)$ a.s. and $E[\xi\mid V]=m_\theta(V)$ a.s.". Then $\theta$ is in its sharp identified set iff, for $P_V$-a.e. $v$ and all $u$,
$$\langle u,m_\theta(v)\rangle\le\int h_{\mathcal S}(u;w)\,P_{W\mid V=v}(dw).$$

*Proof (by inspection of 2a Thm 1).* That proof uses the following properties of $S(x,t)$, and nothing else:
* nonempty values and a Borel selector, step (ii)⇒(i).1. Here use [KRN] with the closed graph;
* compactness and convexity, steps 2, 4 and 5;
* the uniform bound, steps (i)⇒(ii).4 and 10 and (ii)⇒(i).3;
* a Borel, u.s.c. support function and a measurable argmax, steps 12–14. These follow from the closed graph by [Berge] and [MMT];
* the affine structure of $m_\theta$, as in 2b Thm 1′.

Each is a hypothesis of the lemma. [2a:T1; 2b:T1′] ∎

Uses below:
* the product $S(x,t_0)\times S(x,t_1)$ (§2);
* polytopes with interval row totals (§7).

---

## 2. No link between periods

**Proposition D0.**
(a) **Model-free, exact.** Let $[G^-_t,G^+_t]$ be the DD range of the gap $B_{vt}-B_{nt}$. Then the set of $\text{DiD}_c$ values consistent with the margins is $[G^-_1-G^+_0,\ G^+_1-G^-_0]$.
(b) **Per-period LCE (M2 with periods as groups).** The sharp identified set of $(K_0,K_1)$ is $\Theta_{I,0}\times\Theta_{I,1}$, where each factor is computed from its own period's law $P_t$ of $(N,X,T_t)$. The outer set is $\Theta_{O,0}\times\Theta_{O,1}$. So the DiD bound is the difference of the two periods' gap intervals, and it is **sharp** when computed with the sharp sets. *This settles the conjecture in MENU item 17 and 1b D0.*
(c) In 2×2 each period's gap is increasing in $\kappa_t$. So the owner's appendix formula $ATT_{L}=(L_{v1}-U_{v0})-(U_{n1}-L_{n0})$ (and $ATT_U$) equals (b) with JKST outer sets.

*Proof.*
1. (a): no restriction links $\beta_{i0}$ with $\beta_{i1}$, so the feasible set is $\prod_i S(X_i,T_{i0})\times\prod_iS(X_i,T_{i1})$. [Def; 2c:M2(a)]
2. $\text{DiD}_c=\text{gap}_1(\beta_{\cdot1})-\text{gap}_0(\beta_{\cdot0})$ is separable across the two factors. So its range is the difference of the ranges. [sep; conn] Checked: 1b V7a.
3. (b): put $\xi=(\beta_0,\beta_1)$, $\mathcal S(w)=S(x,t_0)\times S(x,t_1)$ and $m_\theta=(M_{H_0}(x),M_{H_1}(x))$. Lemma G applies, since products preserve its hypotheses. [L G; 2a:L0]
4. $h_{\mathcal S}((u_0,u_1);w)=h_S(u_0;x,t_0)+h_S(u_1;x,t_1)$. Integrating against $P_{W\mid V}$ gives $\psi_0(u_0;v)+\psi_1(u_1;v)$, where $\psi_t$ uses only the conditional law of $T_t$ given $V$. [prod-supp; CE-lin; RCD]
5. The joint condition with $u_1=0$ is period 0's condition (T1) (and $h(0)=0$); with $u_0=0$ it is period 1's. Conversely, summing the two separate conditions gives the joint one. [4; alg]
6. So joint sharp = product of the per-period sharp sets [5; L G; 2a:T1], and likewise for the outer sets [Def $\Theta_O$]. The DiD bound follows as in step 2.
7. (c): by 2a Lemma 2, $B_{vt}$ is increasing in $\kappa_t$ (coefficient $\mu_{vn}/\mu_v$) and $B_{nt}$ is decreasing (coefficient $-\mu_{vn}/\mu_n$). So $L_{vt}$ and $U_{nt}$ are attained at the same $\kappa_t$. [2a:L2; alg] ∎

---

## 3. First-difference linear contextual effects (rung FD)

**Model FD.** (S1) holds in each period, and
$$E[\Delta\beta\mid V]=M_{\Delta H}(X)\quad\text{a.s.},\qquad\textstyle\sum_c\Delta H_c=0.$$
**Levels are unrestricted.** Each township may have any persistent political culture, $\zeta_i$ in 1b B2 D1, related to $X$ in any way.

**Lemma FD1.** Assume A5 for the design.
(a) $E[\Delta T_c\mid X]=X^\top\Delta Q_cX$ with $\Delta Q_c=\operatorname{sym}\Delta H_c$. So $\Delta Q$ is identified, and $\Delta K\in\mathcal K$ is left completely unidentified by this regression.
(b) $\widetilde{\text{DiD}}_c=I_c+s\,\Delta\kappa_c$, where
$$I_c=\frac{\Delta Q_{c,vv}\mu_{vv}+\Delta Q_{c,vn}\mu_{vn}}{\mu_v}-\frac{\Delta Q_{c,nn}\mu_{nn}+\Delta Q_{c,vn}\mu_{vn}}{\mu_n},\qquad s=\mu_{vn}\Big(\frac1{\mu_v}+\frac1{\mu_n}\Big)>0.$$
It depends on $\Delta K$ only through $\Delta\kappa_c$.

*Proof.*
1. $\Delta T=\Delta\beta^\top X$. [A1 in both periods]
2. $E[\Delta T_c\mid X]=X^\top\Delta H_cX=X^\top\Delta Q_cX$. [CE-pull; FD; skew parts vanish in quadratic forms]
3. $\Delta Q$ is unique, and $\Delta K$ does not enter the regression. [2a:L1 applied to $\Delta T$; A5]
4. Rows of $\Delta\beta$ sum to 0, so $\sum_c\Delta H_c=0$, and therefore $\sum_c\Delta K_c=0$. [Def; 2a:L1(b) argument]
5. (b): 2a Lemma 2 with $\Delta H$ for $H$ and sample weights, using $K_{c,nv}=-\kappa_c$. [2a:L2; P3 §0.2; alg] ∎

**Proposition FD2 (single restrictions; two rows, one column).** Write $\Delta\kappa^{R1}:=\Delta Q_{vn}-\Delta Q_{nn}$ and $\Delta\kappa^{R2}:=\Delta Q_{vv}-\Delta Q_{vn}$, dropping $c$.
(a) **R1** (no differential spillover onto non-veterans), $\Delta w_1=0$: gives $\Delta\kappa=\Delta\kappa^{R1}$.
(b) **R2** (veterans' change does not depend on fellow veterans), $\Delta b_1=0$: gives $\Delta\kappa=\Delta\kappa^{R2}$.
(c) **R4** (no average contextual change), $a\Delta b_1+b\Delta w_1=0$ with $a=\mu_v\operatorname{Var}_{w^v}(X)$ and $b=\mu_n\operatorname{Var}_{w^n}(X)$: gives
$$\Delta\kappa=\frac{a\,\Delta\kappa^{R2}+b\,\Delta\kappa^{R1}}{a+b}.$$
**So $\widetilde{\text{DiD}}^{R4}$ lies between $\widetilde{\text{DiD}}^{R1}$ and $\widetilde{\text{DiD}}^{R2}$.**
(d) Under R4, the population slope of the linear regression (with intercept) of the model-implied $\Delta T$ on $X$ equals $\widetilde{\text{DiD}}$. That is the appendix's Goodman-on-changes estimator, which needs only R4.
(e) R1 and R2 together imply $\Delta d_1=0$, which is testable. Then all three estimates coincide.

*Proof.*
1. $w_1=H_{nv}-H_{nn}=(Q_{vn}-\kappa)-Q_{nn}$, and $b_1=H_{vv}-H_{vn}=Q_{vv}-(Q_{vn}+\kappa)$. [Def; $H=Q+K$]
2. Taking differences: $\Delta\kappa=\Delta\kappa^{R1}-\Delta w_1=\Delta\kappa^{R2}-\Delta b_1$. [1; alg] This gives (a) and (b).
3. (c): substitute 2 into $a\Delta b_1+b\Delta w_1=0$ and solve. $a,b\ge0$ with $a+b>0$, so $\Delta\kappa^{R4}$ is a convex combination. [alg] $\widetilde{\text{DiD}}$ is affine in $\Delta\kappa$. [L FD1(b)]
4. (d): $a\Delta b_1+b\Delta w_1=\delta^{\Delta*}_W$ by 2c Prop. M5(c) applied to $\Delta m$. [2c:M5(c)] Then EGGHM's identity for the outcome $\Delta m\in[-1,1]$ gives $ER^*_\Delta-\widetilde{\text{DiD}}=\delta^{\Delta*}_W/\operatorname{Var}_p(X)=0$. [2c:L M1(c)]
5. (e): $\Delta b_1-\Delta w_1=\Delta d_1$. [1; alg] ∎

Checked: D1. With arbitrary non-linear levels, each rule recovers the DiD exactly when it holds (max error $4\times10^{-16}$), R4 always lies between R1 and R2, and Goodman on changes equals the DiD under R4.

**Proposition FD3 (bounded change and breakdown).** Under FD with $\lvert\Delta w_1\rvert\le\delta$,
$$\widetilde{\text{DiD}}\in[\widetilde{\text{DiD}}^{R1}-s\delta,\ \widetilde{\text{DiD}}^{R1}+s\delta]\cap\mathcal B^{FD}_O.$$
The same holds with $\Delta b_1$ and R2. For the claim "DiD > 0" with $\widetilde{\text{DiD}}^{R1}>0$, the breakdown value is
$$\delta^*=\widetilde{\text{DiD}}^{R1}/s,$$
or $\infty$ if the lower end of $\mathcal B^{FD}_O$ is already positive. In words: **the smallest differential spillover onto non-veterans (in yes-rate points per unit of veteran share) that would erase the effect.**
*Proof.* $\widetilde{\text{DiD}}=I+s(\Delta\kappa^{R1}-\Delta w_1)$. [FD2 step 2; L FD1(b)] Intersect with the outer set. [Def] ∎ Checked: D1 (slope $-s$).

**Proposition FD4 (outer set without any restriction on $\Delta K$).** Let $\tau_t(x):=E[T_t\mid X=x]$ and let $x_j$ be the support points. Then
$$\Theta^{FD}_O:=\{\Delta k:\ \forall j\ \exists\,m_0\in S(x_j,\tau_0(x_j))\ \text{with}\ m_0+M_{\Delta H(\Delta k)}(x_j)\in S(x_j,\tau_1(x_j))\}$$
contains the true $\Delta k$. The DiD bound is an LP in $(\Delta k,(m_0(x_j))_j)$. Moreover:
* per-period LCE (M2 with periods) ⊆ FD, as model classes;
* when $\tau_t$ are the quadratics, the product of the per-period outer sets maps into $\Theta^{FD}_O$.

*Proof.*
1. $m_t(x):=E[\beta_t\mid X=x]\in S(x,\tau_t(x))$. [2c:L N2.1, per period]
2. $m_1-m_0=E[\Delta\beta\mid X]=M_{\Delta H}(X)$. [CE-lin; FD; A3] So the true $\Delta k$ is feasible, with $m_0$ as witness. [1–2]
3. Levels linear in each period implies changes linear. [Def] If $(\kappa_0,\kappa_1)$ lie in the product of the outer sets, then $M_{H_t}(x_j)\in S(x_j,\tau_t(x_j))$. [2c:N2.2(c) step 3] So $m_0:=M_{H_0}(x_j)$ is a witness for $\Delta k=k_1-k_0$. [Def] ∎

Checked:
* D2(a): the truth is feasible under non-linear levels.
* D2(b): the product bound lies inside the FD bound in 60 cases.

**Note.** $\tau_t$ (levels) is not quadratic under FD, so Phase 3 must estimate it nonparametrically, by a polynomial of chosen degree as in 2c §6.

---

## 4. Monotone EI on changes (exact, no linearity)

**Proposition D4.** In the fixed population, apply EGGHM to the individual change $\Delta Y\in\{-1,0,1\}$. The unit group means are $\Delta\beta_{iv},\Delta\beta_{in}$, and the associations of the realised changes are $\delta^\Delta_B$, $\delta^\Delta_W$. Let $ER_\Delta$ and $NM_\Delta$ be Goodman (with intercept) and the neighbourhood model computed on $\Delta T$.
(a) $ER_\Delta-\text{DiD}=\delta^\Delta_W/\operatorname{Var}_p(X)$; $NM_\Delta-\text{DiD}=-\delta^\Delta_B/(\mu(1-\mu))$; $NM_\Delta=\gamma\,ER_\Delta$.
(b) The exact sets intersect D0's interval $[D^-_0,D^+_0]$ with:
  * $\delta^\Delta_W\ge0$: $\text{DiD}\le ER_\Delta$ (empty iff $ER_\Delta<D^-_0$);
  * $\delta^\Delta_W\le0$: $\text{DiD}\ge ER_\Delta$;
  * $\delta^\Delta_B\ge0$: $\text{DiD}\ge NM_\Delta$ (never empty);
  * CR on changes: DiD lies between $NM_\Delta$ and $ER_\Delta$. **Never empty, and $\operatorname{sign}\text{DiD}=\operatorname{sign}ER_\Delta$.**
(c) If $\delta^\Delta_W=0$ (realised R4), then $\text{DiD}=ER_\Delta$, with no model.

*Proof.*
1. With the same men in both periods, $\Delta Y$ is defined per man, bounded, and has unit means $\Delta\beta_{ir}$. [0.1; Def]
2. (a): [EGGHM:P4, P6, L1] for the outcome $\Delta Y$. [1; 1b V7c]
3. (b): repeat 2c Prop. M3 with D0's product set in place of the MOB. The step "$NM\in$ range" uses $\beta_{ivt}=\beta_{int}=T_{it}$ in both periods, which gives $\Delta\beta_{iv}=\Delta\beta_{in}$, hence $\delta^\Delta_B=0$ and $\text{DiD}=NM_\Delta$. [2c:M3 steps 1–8; D0(a)]
4. (c): (a) with $\delta^\Delta_W=0$. ∎

Checked: D3 (one branch; sign agreement in all 83 populations with realised CR).

**Reading for the application.**
* **Which signs are usable.**
  * $\delta^\Delta_B$ is the $X(1-X)$-weighted average of the *within-township* veteran–non-veteran difference in changes, i.e. the DiD itself with different weights. Signing it presupposes the answer, so it is not recommended.
  * $\delta^\Delta_W\ge0$ means the *net* contextual change is towards suffrage in veteran-heavy townships, as with persuasion spillovers.
* **Consequence for the earlier analysis.** If spillovers are positive, $\text{DiD}\le ER_\Delta$: **the appendix's Goodman-on-changes estimate (0.265 in §C.3) is then an upper bound**. With backlash it is a lower bound.
* **Contextual reinforcement.** CR says the within-township gap in changes and the contextual pattern of changes have the same sign. It is never rejected by the data, and it signs the DiD by the sign of Goodman on changes.

---

## 5. Intersection across periods (M1 with periods as groups)

2b Thm 2 and Remark 6.2 apply unchanged. **Correction (1b §B3; checked 1b V8):** "contextual effects unchanged" ($\Delta H_c=\zeta_c\mathbf 1^\top$) means **R1 and R2 together**, not $\Delta\kappa=0$; in fact $\Delta\kappa=(\zeta_v-\zeta_n)/2$. M1-periods therefore lies inside FD+R1+R2, adds level linearity in each period, and tightens the levels by intersection. Its testable part is $\Delta d_1=0$.

---

## 6. Changing population counts

**Proposition D6.** Suppose the periods have different counts $N_{it}$ (and possibly $X_{it}$), with per-period moments $\mu^{(t)}$. Assume per-period LCE and R1 across periods, $w_{1}^{(1)}=w_{1}^{(0)}$. Then
$$\widetilde{\text{DiD}}=I'+(s_1-s_0)\,\kappa_0,\qquad s_t=\mu^{(t)}_{vn}\Big(\frac1{\mu^{(t)}_v}+\frac1{\mu^{(t)}_n}\Big),$$
with $I'$ identified and $\kappa_0\in\Theta_{O,0}\cap(\Theta_{O,1}-\Delta\kappa^{R1})$. **The DiD is therefore an interval of width $\lvert s_1-s_0\rvert\times$ (length of that set), not a point.** It is a point when $s_1=s_0$, e.g. when growth is unrelated to composition.

*Proof.*
1. Each period's rates are affine in $\kappa_t$ with weight $\mu^{(t)}_{vn}/\mu^{(t)}_r$. [2a:L2]
2. R1 gives $\kappa_1=\kappa_0+\Delta\kappa^{R1}$. [FD2 step 2, per period]
3. Substitute; $\kappa_0$ must satisfy both periods' constraints. [1–2; 2b:Thm 2(b)] ∎

Checked: D4 (formula exact; truth inside).

FD (§3) requires the *same men*. With changing populations, only level models (§§2, 5, 6) or the composition accounting of §7 apply.

---

## 7. Composition: aging-in, exits and interval bases

**Accounting polytope.** For unit $i$ and period $t$, let $\mathcal R_t$ be a set of row types, for example:
* $t=0$: veteran-stayer, non-veteran-stayer, veteran-exit (died or did not return), non-veteran-exit;
* $t=1$: veteran-stayer, non-veteran-stayer, veteran-aged-in, non-veteran-aged-in.

Row counts $\rho_{irt}\ge0$ are variables subject to linear constraints:
* known or interval totals: `elig` from 1860 vs 1870, at least the ballots cast, at least the veterans;
* known veteran totals (`vets_alt`, `_s`, `_h`);
* **links**: stayers are the same men, $\rho_{i,\text{v-stay},0}=\rho_{i,\text{v-stay},1}$;
* **structural zeros**: aged-in men have no period-0 row.

Cells $y_{irct}\ge0$ satisfy $\sum_cy_{irct}=\rho_{irt}$ and $\sum_ry_{irct}=$ the observed votes for $c\in\{$yes, no$\}$, with abstain $\ge0$ as the remainder.

**Proposition C1.** The feasible set $\mathcal T$ of $(\rho,y)$ is a polytope.
(a) A target that is a ratio of linear forms with a positive denominator (e.g. veteran stayers' 1865 yes-rate) has an interval of values, computable as an LP. [CC; conn]
(b) Fix the denominators $V_O=\sum_i\rho_{i,\text{v-stay}}$ and $N_O=\sum_i\rho_{i,\text{n-stay}}$. Then the stayers' DiD is linear, so for those values its range is an LP.
  * Over all feasible $(V_O,N_O)$ the set is an interval.
  * When total stayers are known, the scan is one-dimensional.
  * The endpoint functions are continuous, so refining the grid converges to the exact interval.
(c) The bounds are exact in the finite-population sense, given the type accounting (validity and sharpness as in 2c Prop. M2). Computing them by a grid is an approximation, which will be reported as such.

*Proof.*
1. All constraints are linear in $(\rho,y)$. [Def]
2. (a): [CC] with $\mathcal T$ connected. [conn]
3. (b): for fixed denominators the DiD is linear in $y$. [Def] Its value function is continuous in $(V_O,N_O)$ on the feasible region. [LP-par] The image of the connected set $\mathcal T\cap\{\text{denominators}>0\}$ is an interval. [conn]
4. (c): as in 2c M2(a). [2c:M2] ∎

Checked: D6 (aged-in split interval-valued; the bound contains the stayers' DiD in 15 random populations).

**Proposition C2 (what the post-war base mixes).** Put the 1857 votes over the post-war base, as in the owner's earlier analysis. Let $\omega_r=O_r/(O_r+A_r)$ (stayers' share of group $r$'s post-war men). Then
$$\text{DiD}^{\text{orig}}=\big[\omega_v\Delta B^{O}_v+(1-\omega_v)B^{A}_{v1}\big]-\big[\omega_n\Delta B^{O}_n+(1-\omega_n)B^{A}_{n1}\big].$$
If aged-in men vote like stayers of their group in 1865 ($B^A_{r1}=B^O_{r1}$), then
$$\text{DiD}^{\text{orig}}-\text{DiD}^{O}=(1-\omega_v)B^O_{v0}-(1-\omega_n)B^O_{n0}.$$
It lies between $-(1-\omega_n)B^O_{n0}$ and $(1-\omega_v)B^O_{v0}$, so it is small when pre-war yes-rates are near zero. Veterans are younger, so $\omega_v<\omega_n$.

*Proof.*
1. $B^{\text{orig}}_{r0}=\omega_rB^O_{r0}$, because aged-in men cast no 1857 votes. [Def; structural zero]
2. $B_{r1}=\omega_rB^O_{r1}+(1-\omega_r)B^A_{r1}$. [Def]
3. Subtract, then substitute the assumption. [alg] ∎

Checked: D5.

**Remark C3 (exits).** The 1857 pool of group $r$ includes men who later died or did not return. Its rate is $B^{\text{pool}}_{r0}=\pi_rB^O_{r0}+(1-\pi_r)B^E_{r0}$.
* The owner's earlier assumption is $B^E=B^O$.
* The dial $\lvert B^E_{r0}-B^O_{r0}\rvert\le\eta$ puts the stayers' pre-war rate within $(1-\pi_r)\eta$ of the pool rate, since $B^{\text{pool}}_{r0}-B^O_{r0}=(1-\pi_r)(B^E_{r0}-B^O_{r0})$. [alg]
* Model-free, the exit rows are free cells in $\mathcal T$.

---

## 8. The causal layer (owner's appendix §C.4)

Let $\theta$ be veterans' post-war counterfactual rate had they not served; $ATT=B_{v1}-\theta$. Assumptions:
* **(MTR)** $\theta\le B_{v1}$: service did not lower support [MTR];
* **(ND)** $\theta\ge B_{v0}$: absent service, no decline from the pre-war level;
* **(PS)** $\theta\ge B_{n1}$: future soldiers would have been at least as supportive as non-veterans; a positive-selection, [MTS]-type assumption.

**Proposition CL.** Let $\mathcal F$ be any rung's (convex) identified set for $(B_{v0},B_{v1},B_{n1})$.
(a) The ATT set is non-empty iff some $b\in\mathcal F$ has $b_{v1}\ge\max\{b_{v0},b_{n1}\}$. So **MTR+ND+PS are jointly refutable.**
(b) If non-empty, the ATT set is $[0,U]$ with
$$U=\max_{b\in\mathcal F}\min\{b_{v1}-b_{v0},\ b_{v1}-b_{n1}\},$$
an LP in epigraph form.
(c) With two rows and the product structure of D0, $U=B^U_{v1}-\max\{B^L_{v0},B_{n1}(B^U_{v1})\}$. This is the owner's script line 336, which is therefore attained.
(d) Under parallel trends, ATT = DiD. Bounded departures give ATT ∈ DiD ⊕ [−M, M] (Rambachan & Roth 2023).

*Proof.*
1. ATT set $=\{b_{v1}-\theta:b\in\mathcal F,\ \max\{b_{v0},b_{n1}\}\le\theta\le b_{v1}\}$. [Def]
2. (a): the inner interval is non-empty iff $b_{v1}\ge\max\{b_{v0},b_{n1}\}$. [1]
3. (b): each feasible $b$ contributes $[0,\,b_{v1}-\max\{b_{v0},b_{n1}\}]$. All contain 0, so the union is an interval. [1; conn] Its upper end is the stated max; $\min$ of affine functions is concave, so it is an LP in epigraph form. [lin]
4. (c): within period 1, $\mu_vB_{v1}+\mu_nB_{n1}=E_p[T_1]$ for both actual and model-implied rates. So $B_{n1}$ is decreasing in $B_{v1}$. [A1; Def]
5. Both terms of the min increase in $b_{v1}$. $b_{v0}$ ranges over its own interval independently. [D0(b)] So the max is at $b_{v1}=B^U_{v1}$ and $b_{v0}=B^L_{v0}$. [4; alg]
6. (d): definition of parallel trends. ∎

Checked: D7 (epigraph LP = closed form, 200 cases).

**Citation note.** The appendix attributes the bounds to Athey & Imbens (2006, pp. 151–3) "with one change". The derivation above is self-contained (MTR: Manski 1997; MTS: Manski & Pepper 2000). The Athey–Imbens page reference should be checked before publication.

---

## 9. Partisanship

**Setup** (1857, unit $i$). Groups $g\in\{v,n\}$ (future veterans or not); party $q\in\{R,D,G\}$ ($G$ = did not vote for governor); suffrage $\in\{$yes, other$\}$. Observed: the one-way margins $V$, $(R,D,G)$ and $(Y,\cdot)$.

**Proposition P2 (the "Democrats never yes" zero alone does not help).** If $Y\le R+G$, the zero $n_{\cdot,D,\text{yes}}=0$ leaves the set of feasible veteran-yes counts at Duncan–Davis, $[\max(0,V+Y-N),\min(V,Y)]$.
*Proof.*
1. Take any veteran × suffrage table with the given margins.
2. Take a party × suffrage table with $z_{D,\text{yes}}=0$. It exists by [Gale], since $Y\le R+G$.
3. Join them by $n_{gqs}=y_{gs}z_{qs}/Y_s$. This has both tables as two-way margins and keeps the zero. [coupling] ∎

Checked: D8 (300 cases, LP over three-way tables = DD).

**Proposition P3 (selection on party only).** Assume:
* **(SP)** in each township, within each party, future veterans and others had the same 1857 yes-rate $\pi_{iq}$;
* **(Z)** yes votes came only from Republican governor voters.

Then $\pi_{iR}=Y_i/R_i$, and veteran-yes $=v_{iR}Y_i/R_i$ with $v_{iR}\in[\max(0,V_i+R_i-N_i),\min(V_i,R_i)]$. The range
$$\Big[\max(0,V_i+R_i-N_i)\tfrac{Y_i}{R_i},\ \min(V_i,R_i)\tfrac{Y_i}{R_i}\Big]$$
is exact under SP+Z and lies inside DD; it is often strictly narrower (152 of 300 random cases in D8). The aggregate is the weighted sum of the unit ranges.

*Proof.*
1. SP+Z give unit veteran-yes $=v_{iR}\pi_{iR}$ with $\pi_{iR}=Y_i/R_i$. [Def; alg]
2. $v_{iR}$ ranges over the veteran × party DD interval, all values attainable. [2a:L0; [FH]]
3. Upper end: $\min(V,R)\,Y/R\le\min(VY/R,\,Y)\le\min(V,Y)$, since $Y\le R$. [alg]
4. Lower end: $(V+R-N)Y/R-(V+Y-N)=(1-Y/R)(N-V)\ge0$. [alg]
5. Aggregation is separable across units. [sep] ∎

Checked: D8 (inside DD and containing the truth in all 300 cases; strictly narrower in 152).

**Without (Z)** (e.g. some suffrage voters skipped the governor's race, or Democrats' yes-rate is $\le\varepsilon$), veteran-yes $=\sum_qv_{iq}\pi_{iq}$ is bilinear over a product of polytopes. Exact per-unit ranges then come from vertex enumeration. [bilin; conn]

**Uses:**
* P3 tightens the *pre-war* veteran rate, which is the widest part of the DiD.
* Post-war, SP is not assumed (the effect of service is the point).
* A party DiD (veteran × party in 1857 and 1865) uses §§2–4 unchanged.

---

## 10. The DiD ladder, with status

Decreasing plausibility; ⊇ means nested (the bounds shrink). Applies to the fixed population and, through §7, to stayers.

$$D0\ (\text{exact})\ \supseteq\ D0\cap\text{MEI-}\Delta\ (\delta^\Delta_W\text{ sign or CR; exact})$$
$$\mathcal B^{FD}_O\ \supseteq\ \mathcal B^{FD}_O\cap\{\lvert\Delta w_1\rvert\le\delta\}\ \supseteq\ \{\text{R1 point}\}\qquad[\text{FD3}]$$
$$\text{per-period LCE (product)}\ \subseteq\ \text{FD};\qquad \text{M1-periods}\ \subseteq\ \text{FD}+\text{R1}+\text{R2}\qquad[\text{FD4; }\S5]$$
**Not nested:** MEI-Δ with FD (no linearity vs linearity of changes). R4 lies between R1 and R2 (FD2(c)). Party (P3) and composition (§7) modify every rung. The causal layer (§8) sits on top of any rung.

## 11. What Phase 3 must supply

* **Exact rungs** (D0, MEI-Δ, the composition polytope, P3): no sampling inference. Add slack dials and breakdown values (2c M2(d)), plus the grid for C1(b).
* **FD:**
  * the regression of $\Delta T$ on the homogeneous quadratic, clustered by county;
  * $\tau_0,\tau_1$ by polynomial series, for the outer set;
  * CLR over the support-point constraints (P3 Thm 5);
  * the breakdown value δ* (FD3) with its CI.
* **Changing counts (D6)** and **interval bases** (§7): inference conditional on the census-based bases, with the 1860/1870 interval treated as known.
* **The two targets** (1b §D.3): the actual DiD for the exact rungs, and the model-implied DiD for FD and LCE rungs. The actual-DiD allowance uses $(\sqrt{m_1(1-m_1)}+\sqrt{m_0(1-m_0)})^2$.
