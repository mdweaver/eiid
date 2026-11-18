# Phase 2c: monotone ecological inference in R×C, derivations and proofs

*Builds on P1, 2a, 2b, P3 and the Phase-1 proposal `PHASE1b_monotone_did.md` (cited as 1b). It follows Elzayn, Goldin, Guage, Ho & Morton, "Monotone Ecological Inference" (EGGHM; arXiv v1 numbering, see 1b §0). Every proof step ends with a justification tag: those of 2a §0.3, 2b §0.3 and P3 §0.3, plus §0.3 below. Checks: `verify_phase2c.py` (output `verify_phase2c_output.txt`, 26/26 PASS). Owner decisions applied (1b §D):*
* *finite-population bounds for the sign rungs, **sharp for continuous unit tables given the margins** (corrected after review, Batch B: earlier "exact"; see Remark M2.1 for integer counts, decision U1(a));*
* *one-vs-rest (W1) as the default, with W2 as an option;*
* *contextual reinforcement available as an option;*
* *API: planned here as `monotone=` plus `model=`. The shipped M-SIGN module uses `restrictions = list(assoc = ..., cr = ...)` instead (5a §5). N1, N2 and local monotonicity are not implemented.*

*Batch B review (Nov 2026). Two external referee reports found errors in M2 (integer tables), M3(c)–(d), M5(c) (district rungs), N2.2 and N2.3, and §6 (binning, degree selection). All are corrected below and marked; the ledger is `reviews/B_PI_RESPONSE_PLAN.md` §2 (rows P-B1 to P-B14, P-B28, P-B29, P-B31, P-B32, P-B34, P-B35; P-E1, P-E3, P-E8, P-E11 to P-E13, P-E23, P-E27, P-E28). Every correction has a check in `verify_phase2c.py` that fails on the old statement and passes on the new one (C7–C14), unless marked otherwise. Status labels follow `REVIEW_RESPONSE_PLAN.md`.*

---

## 0. Setting

### 0.1 Finite population (the exact rungs)

* The units $i=1,\dots,p$ are fixed (the design is conditioned on, as in P3 (P3-S2)). Unit weights are $p_i=N_i/\sum_jN_j$.
* $X_i\in\Delta^{R-1}$, $T_i\in\Delta^{C-1}$, and the realised rates satisfy $\beta_i\in S(X_i,T_i)$ (2a §0.2).
* $E_p[f]:=\sum_ip_if_i$, and $\mu_r:=E_p[X_r]>0$.
* Row weights $w^r_i:=p_iX_{ir}/\mu_r$. $\operatorname{Cov}_w$ is the covariance under weights $w$.
* **Target** (P3's actual rate): $B_{rc}=E_p[X_r\beta_{rc}]/\mu_r=\sum_iw^r_i\beta_{irc}$.
* **Continuous-rate reading.** The realised rates $\beta_i$ are real numbers: the cell masses $N_iX_{ir}\beta_{irc}$ need not be whole numbers of people. All N1 results below are for this reading. With counted ballots (integer cell counts) attainability and refutation can change (Remark M2.1). *(Added after review, Batch B.)*

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
3. **Within, W2 (within each original group, pooled into one sum):** $\Delta:=\sum_s\mu_s\operatorname{Cov}_{w^s}(X,\beta_s)\in\mathbb R^R$. A sign on $\Delta_r$ signs **one sum** of row-specific covariances; the summands need not share its sign. *(Corrected after review, Batch B: earlier labelled "each group separately".)*
4. **Within, W1 (one-vs-rest; EGGHM §4 end).** Let $\beta^{\neg r}_i:=\sum_{s\ne r}X_{is}\beta_{is}/(1-X_{ir})$, the rate of everyone outside row $r$ (any value if $X_{ir}=1$, where its weight is 0), and $w^{\neg r}_i:=p_i(1-X_{ir})/(1-\mu_r)$. Then
$$\delta^r_W:=\mu_r\operatorname{Cov}_{w^r}(\beta_r,X_r)+(1-\mu_r)\operatorname{Cov}_{w^{\neg r}}(\beta^{\neg r},X_r).$$

**Remarks on Definition 1** *(added after review, Batch B)*.
* **Weights.** $\delta^r_B$ is the A-minus-rest gap averaged with weights $p_iX_{ir}(1-X_{ir})$: pure units get zero weight, and a large mixed unit can dominate. It is not an ordering in every unit, nor against each rival. The pairwise $\delta_B(r,s)$ uses $p_iX_{ir}X_{is}$.
* **W1 versus W2.** In W1 the predictor in **both** covariances is $X_r$ (row $r$'s share), and the complement's rate is that of a changing mix of the other rows; either covariance can be negative when their sum is positive. W2 conditions on the original rows. The two can have **opposite signs** while every group's rate rises. Example (EI referee): rows A, B, C, two units of equal size, shares $(.2,.7,.1)$ and $(.6,.05,.35)$, rates $(.40,.80,.10)$ and $(.42,.82,.12)$: every rate rises by 2 points, yet $\delta^A_W=-0.02633$ and $\Delta_A=+0.00110$ (check C14). A direction for each named group needs separate restrictions.
* **Goodman bias.** $g-B=M^{-1}\Delta$ (Lemma M1(a)): signing $\Delta_r$ does not in general sign the bias of row $r$'s Goodman coefficient.
* **Scope relative to EGGHM.** The one-vs-rest collapse W1 is EGGHM §4's suggestion. W2, the pairwise weighting, and CR with those pairings are extensions with their own proofs here (M1, M4); they are not EGGHM results. Covariances also reflect sorting and other between-unit differences; a causal reading needs a separate design.

**Estimators** (functions of the data only):
* Goodman's no-intercept regression: $g:=M^{-1}E_p[XT_c]$, with $M:=E_p[XX^\top]$.
* The neighbourhood model: $NM_r:=E_p[X_rT_c]/\mu_r$.
* The one-vs-rest regression with intercept: $ER^r:=\operatorname{Cov}_p(T_c,X_r)/\operatorname{Var}_p(X_r)$.
* The one-vs-rest gap: $D^r:=B_r-B_{\neg r}$, with $B_{\neg r}:=\sum_{s\ne r}\mu_sB_s/(1-\mu_r)$.

**Lemma M1 (identities; accounting identity only).** Assume [A1]. Part (a)'s formula for $g$ needs $M$ invertible (the $X_i$ span $\mathbb R^R$); $\mathbf 1^\top\Delta=0$, (b) and (d) need nothing more; (c) needs $\operatorname{Var}_p(X_r)>0$ for $ER^r$ and $0<\mu_r<1$. Where $X_{ir}=1$, give the absent rest group any bounded value (its weight is 0). *(Corrected after review, Batch B: earlier "Assume $M$ is invertible" for every part.)*
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

## 2. Bounds under sign restrictions (rung N1): sharp for continuous tables

**Definition 2.**
* A *sign restriction* is $g_k(\beta)\ge-\varepsilon_k$ with $g_k=\pm$ one of the associations of Definition 1, for any column, and slack $\varepsilon_k\ge0$ (default 0).
* Optional *unit-level* restrictions $\beta_i\in U_i$, with $U_i$ polyhedral: unit-level ordering, structural zeros (MENU A1, B5).
* The **N1 feasible set** is
$$F:=\{\beta:\ \beta_i\in S(X_i,T_i)\cap U_i\ \forall i;\ g_k(\beta)\ge-\varepsilon_k\ \forall k\}.$$
* **Contextual reinforcement** for a chosen pair $(g_B,g_W)$ is $F_{CR}:=F_+\cup F_-$, where $F_\pm:=F\cap\{\pm g_B\ge0,\pm g_W\ge0\}$.

**Proposition M2 (N1 is sharp for continuous unit tables).** Let $\varphi=\sum a_{rc}B_{rc}$, and read the unit tables as nonnegative **real** cell masses (§0.1).
(a) The set of values of $\varphi$ attained by continuous unit tables consistent with the margins and the restrictions is $\varphi(F)$. If $F\neq\emptyset$ it equals $[\min_F\varphi,\max_F\varphi]$: **sharp for continuous unit tables, conditional on the observed margins**.
(b) $F=\emptyset$ iff the restrictions contradict the margins in this continuous problem.
(c) Under CR the set is $\varphi(F_+)\cup\varphi(F_-)$, a union of at most two intervals.
(d) $\varepsilon\mapsto\varphi(F(\varepsilon))$ is nested **componentwise** ($\varepsilon\le\varepsilon'$ coordinatewise implies $\varphi(F(\varepsilon))\subseteq\varphi(F(\varepsilon'))$), with continuous piecewise-linear endpoints on the region where $F(\varepsilon)\ne\emptyset$. With several slacks, the set of $\varepsilon$ that preserves a given conclusion is in general a **breakdown frontier**, not a product of intervals (check C7b). A scalar breakdown value needs a specified conclusion, a rule for equality, and a scalar path $\varepsilon=t\,v$, $t\ge0$; along such a path it is well defined. (The package's `breakdown()` uses one scalar slack per module, i.e. a stated path.)

*(Corrected after review, Batch B, Nov 2026: (a) earlier said "attained by unit tables … sharp in the finite-population sense", without the continuous reading; (b) earlier said "contradict the data"; (d) earlier said "the breakdown value of a conclusion is well defined" for vector slacks. Checks C7, C7b.)*

*Proof.*
1. A table $n_{irc}=N_iX_{ir}\beta_{irc}$ of nonnegative real masses reproduces the margins iff $\beta_i\in S(X_i,T_i)$. [A1; Def $S$]
2. It satisfies the restrictions iff $\beta\in F$. [Def 2; L M1(d)]
3. $F$ is a bounded polyhedron, since each $S(X_i,T_i)$ is compact (2a L0(b)) and the restrictions are linear. $\varphi$ is affine in $\beta$. [L M1(d)]
4. So $\varphi(F)$ is a compact interval with both endpoints attained, which gives (a). (b) is the definition. [conn; 3]
5. (c): $F_+$ and $F_-$ are each of the form in (a). [Def 2; 4]
6. (d): $\varepsilon$ enters only the right-hand side, so the sets are nested componentwise. Right-hand-side LP duality gives continuous polyhedral (piecewise-linear) endpoint functions on the region where $F(\varepsilon)\neq\emptyset$. A conclusion such as "upper endpoint $\le c$" or "width $\le c$" holds on a sublevel set of such a function, which need not be a product set. [LP-par; nest] ∎

**Remark M2.1 (integer counts; added after review, Batch B).** Step 1 gives real cell masses, not whole numbers of people. With counted ballots, M2(a) and (b) can fail.
* *Referee example (partial identification).* Three units of size 4 with group-1 counts 1, 2, 3 and positive outcomes 0, 1, 0. Each group has six people, so the integer gaps are exactly $\{-1/6,1/6\}$. By symmetry $D_{ER}=D_{NM}=0$ and $\gamma=1/6$: the continuous CR set is $\{0\}$ (half of the positive outcome in each group), while **no integer table satisfies CR**.
* *Referee example (ecological inference).* Two four-voter precincts with 1 and 2 A voters and one P vote each. The integer gaps are $-2/5$, $2/15$, $2/3$; CR requires $D=0$, so integer CR is empty. Under the positive between sign the continuous lower endpoint is 0, the integer one $2/15$: an **endpoint** moves, not only emptiness.
* **Proposition (2×2, margins only, sign rows on one column).** Let $N_1,N_0$ be the group totals and $N=N_1+N_0$. The integer gap set is the Duncan–Davis interval $[D^-,D^+]$ intersected with the lattice $\{-Y/N_0+j(1/N_1+1/N_0)\}$, whose step is $1/(N\mu(1-\mu))$; both endpoints are attained. *Proof.* $D=Y_1(1/N_1+1/N_0)-Y/N_0$, where $Y_1$ is the number of group-1 positives and $Y$ the total. Each unit's group-1 count ranges over an integer interval, so their sum takes every integer between its extremes, which are the DD extremes. ∎ Each sign row of M3 is a half-line in $D$ (M3(a)), so the integer set under sign rows is the continuous set intersected with the lattice. Endpoints therefore differ by less than one step, and integer emptiness can occur only when the continuous set is shorter than one step. The step is about $6\times10^{-6}$ at 370 units of 500–4,000 people, so in this setting the difference matters only for tiny tables or near-ties. Checks C7, C7a.
* **Open (proof item R15).** With unit restrictions $U_i$, joint columns, or $R\ge3$, the integer–continuous gap is not bounded by this argument. Integer N1 would need integer cell variables, discrete target sets and separate integer feasibility checks per branch (decision U1: not built; continuous N1 is labelled as such).
* EGGHM's sharpness results range over compatible joint distributions on the continuous mean-probability model; they do not license whole-ballot interpolation either.

**Proposition M3 (two rows: closed forms, refutability, CR).** Let $R=2$, write $X=X_1$, and fix a column. Write $D=B_1-B_2$, $\gamma=\operatorname{Var}_p(X)/(\mu(1-\mu))$, and let $[D^-,D^+]$ be the MOB interval (DD for $D$). Assume $\operatorname{Var}_p(X)>0$, the continuous reading (§0.1), **no unit-level restrictions** ($U_i$ the whole simplex) and no sign rows other than those named. *(Corrected after review, Batch B: the last two hypotheses were missing; the proof of (c) uses them. Check C8b.)*
(a) $\delta_W\ge0\iff D\le D_{ER}$, and $\delta_B\ge0\iff D\ge D_{NM}$.
(b) $\{\delta_W\ge0\}$ gives $[D^-,\min\{D^+,D_{ER}\}]$. It is **empty iff $D_{ER}<D^-$**. Symmetrically, $\{\delta_W\le0\}$ is empty iff $D_{ER}>D^+$.
(c) $D_{NM}\in[D^-,D^+]$ always. So in this continuous, model-free problem, sign restrictions on $\delta_B$ are never refuted: they give $[D_{NM},D^+]$ or $[D^-,D_{NM}]$. Unit restrictions or other sign rows can remove the witness of step 4 (C8b). With integer counts a single weak between sign still has an integer witness on its side of $D_{NM}$, since the integer extremes are the DD extremes (Remark M2.1).
(d) $D_{NM}=\gamma D_{ER}$ with $\gamma\in(0,1]$. Under CR:
* if $\gamma<1$ (some unit with positive weight is mixed) and $D_{ER}\neq0$, exactly one of $F_\pm$ is non-empty;
* if $\gamma=1$ (every unit with positive weight is pure), $D$ is point-identified, $D^-=D^+=D_{ER}=D_{NM}$, and both branches equal $\{D_{ER}\}$;
* if $D_{ER}=0$, both branches give $\{0\}$. This singleton is the **image** of the gap; it corresponds to many latent tables.

The identified set is
$$[D_{NM},\min\{D^+,D_{ER}\}]\ \text{if }D_{ER}\ge0,\qquad[\max\{D^-,D_{ER}\},D_{NM}]\ \text{if }D_{ER}<0.$$
In this continuous, model-free problem CR is never refuted, and under it $\operatorname{sign}D=\operatorname{sign}D_{ER}$ (for $D_{ER}\neq0$). *(Corrected after review, Batch B: (d) earlier said "exactly one of $F_\pm$ is non-empty when $D_{ER}\neq0$", which fails at $\gamma=1$; the proof noted the case and the statement dropped it. "CR is never refuted by the data" now names its base system. Checks C8, C7.)*

*Proof.*
1. $\delta_W=(D_{ER}-D)\operatorname{Var}_p(X)$ and $\delta_B=(D-D_{NM})\mu(1-\mu)$. [L M1(c) with $R=2$, where $\neg1=2$] This gives (a). [$\operatorname{Var}_p(X)>0$]
2. $D$ is affine in $\beta$ over the product of the unit polytopes. Its range is $[D^-,D^+]$, with every value attained. [L M1(d); M2(a); EGGHM:P3]
3. (b): $F=\{\beta:D(\beta)\le D_{ER}\}$. Its image is $[D^-,D^+]\cap(-\infty,D_{ER}]$, which is empty iff $D_{ER}<D^-$. [1; 2; conn]
4. (c): put $\beta_{i1}=\beta_{i2}=T_i$ in every unit. This lies in $S(X_i,T_i)$ (and in $U_i$ only because $U_i$ is the whole simplex), gives $\delta_B=0$, and hence $D=D_{NM}$. [Def $S$; hypothesis; 1]
5. (d): $D_{NM}=\gamma D_{ER}$ [EGGHM:L1]. $\gamma>0$ by $\operatorname{Var}_p(X)>0$. $\gamma\le1$ because $\mu(1-\mu)-\operatorname{Var}_p(X)=E_p[X(1-X)]\ge0$, with equality iff every unit with $p_i>0$ is pure. [alg]
6. Suppose $D_{ER}>0$. Then $F_-$ requires $D\ge D_{ER}$ and $D\le D_{NM}=\gamma D_{ER}\le D_{ER}$. So $D=D_{ER}=D_{NM}$, which forces $\gamma=1$. Hence $F_-=\emptyset$ when $\gamma<1$. When $\gamma=1$ every positive-weight unit is pure, each unit's rate for its only group equals $T_i$, so $D$ is point-identified and equals $D_{ER}=D_{NM}$, and both branches equal $\{D_{ER}\}$. [5; 1]
7. $F_+=\{D_{NM}\le D\le D_{ER}\}\cap[D^-,D^+]=[D_{NM},\min\{D^+,D_{ER}\}]$. This is non-empty because $D^-\le D_{NM}\le D_{ER}$ and $D_{NM}\le D^+$. [1; 4; 5]
8. Then $D\ge D_{NM}>0$. The case $D_{ER}<0$ is symmetric, and $D_{ER}=0$ gives $D_{NM}=0=D$ on both branches. [5–7] ∎

Checked: C2 (300 populations; LP = closed forms, $\max$ error $6\times10^{-15}$; emptiness iff $D_{ER}<D^-$; one CR branch when $\gamma<1$; sign agreement in all 152 populations satisfying realised CR); C8 ($\gamma=1$: both branches); C8b (a unit ordering removes the witness); C7 (integer counts).

**Reading.** In two rows, in the **continuous, model-free** problem with no unit restrictions, only signs on the *within* association can be contradicted by the margins; a between sign and CR are always feasible there. That non-refutability is specific to this base system:
* with a linear rung (LCE), a between sign and CR can be **jointly** rejected even when LCE alone fits (5a Remark R3.2);
* with integer counts, CR can be rejected outright (Remark M2.1);
* known cells, structural zeros or unit orderings can remove the equal-rates witness (C8b).

The equal-rates witness ($\beta_{ir}=T_i$ for every row) works for **all rows and columns** in R×C too: every pairwise $\delta_B$ is then 0, so a weak $\delta_B$ sign and both weak CR branches hold, in the continuous problem without unit restrictions (EI referee). Where a belief survives, its credibility still rests on substantive argument (or auxiliary individual data, as in EGGHM §6). An empty plug-in set, an empty population set and a rejection by a calibrated test are three different events. *(Corrected after review, Batch B: earlier "between signs and CR can never be rejected from aggregates".)*

**Proposition M4 (joint restrictions).** Imposing sign restrictions on several columns or rows jointly gives a subset of $F$ for any single restriction. So joint bounds are weakly tighter.
*Proof.* Adding constraints shrinks the feasible set. [Def 2; nest] ∎

Checked (C3): in 2×3, signing the within association of two columns jointly narrowed the bound for column 1 strictly in 75 of 139 random populations (max width gain 0.27). This is why the package solves one LP rather than applying EGGHM's 2×2 closed forms column by column. The referee's two-column example makes the point by hand: two 2×2 units at $x=.25,.75$, opposite within signs on the two columns jointly force $D=ER$, shrinking $[-.5,0]$ to $\{0\}$, while each single restriction is feasible.

---

## 3. Sign restrictions under linear contextual effects (K-rows)

**Proposition M5.** Under (S2) for an LCE rung, with the conditional means written for the fixed design as
$$m_i=m_{\rm id}(V_i)+L_ik,$$
where $m_{\rm id}$ is identified and $k$ stacks the free coordinates of that rung: under M0 the coordinates of $K$ (2a Prop. 7), $m_i=M_{Q+K}(X_i)$; under M2 the district-specific coordinates $K^{(d)}$; under M1 the common $K$, with the district intercept shifts in the identified part (2b Lemma D1(b)); under M1(δ) possibly also a common centre that enters the association rows with zero coefficient.
(a) **Every rung.** $\delta^*_B(r,s)$, $\delta^{r*}_B$, $\Delta^*$ and $\delta^{r*}_W$ are affine in $k$. So each sign restriction is one LP row added to $Gk\le g$, and CR is a union of two LPs.
(b) **Ordering, per district.** If $m_{rc}(x,d)\ge m_{sc}(x,d)$ for every district $d$ and every $x$ in that district's support (or $A_d$), then $\delta^*_B(r,s)\ge0$, pooled over districts.
(c) **Slopes and the within sign, $R=2$.**
* *M0 (common lines).* $\delta^*_W=a(w_1+d_1)+b\,w_1$ with $a=\mu\operatorname{Var}_{w^1}(X)\ge0$ and $b=(1-\mu)\operatorname{Var}_{w^2}(X)\ge0$. Hence both contextual slopes $\ge0$ ⇒ $\delta^*_W\ge0$.
* *District rungs.* By the law of total covariance under each row weight,
$$\delta^*_W=\sum_{s}\mu_s\Big(\underbrace{E_{w^s}\big[\operatorname{Cov}_{w^s}(X,m_s\mid D)\big]}_{\text{within districts}}+\underbrace{\operatorname{Cov}_{w^s}\big(E_{w^s}[X\mid D],E_{w^s}[m_s\mid D]\big)}_{\text{between districts}}\Big).$$
Nonnegative district slopes make the within-district term $\ge0$ (as for M0, district by district). The between-district term has no sign: under M1 it contains $\sum_s\mu_s\operatorname{Cov}_{w^s}(X,\zeta_{s,D})$, the covariance of the district intercept shifts with composition. So the shortcut "positive slopes ⇒ $\delta^*_W\ge0$" holds for M0, and within each district (district-conditional weights), but **not** for the pooled association under M1/M1(δ)/M2.
* *Example (referee).* Six equally weighted units, $X=.1,.2,.3$ in district 1 and $.7,.8,.9$ in district 2; both rows' mean $.8+.1x$ in district 1 and $.2+.1x$ in district 2 (the other column the complement). This satisfies M1 with common strictly positive slopes and valid rates. Within each district $\delta^*_W=+0.00064$, but the pooled $\delta^*_W=-0.04927$: the district intercept falls as $X$ rises. (Check C9; the plan's check I2 also shows that the 2×2 local rows hold everywhere here, so the chain "signed slopes ⊂ weighted-local ⊂ $\delta_W\ge0$" also fails for pooled M1.)

*(Corrected after review, Batch B, Nov 2026. M5 was stated for every LCE rung with the single formula $m=M_{Q+K}(X)$, which proves only the pooled case. (a) now uses the stacked parameterisation of each rung; (b) needs district-specific ordering; (c) was false for the district rungs, where the proof treated the district intercept as constant. By decision U2(a), the pooled statements are the proved ones; sharpness and inference for M-SIGN under district rungs (R3′, R10′) are proof item **R16**, not proved. The W1 slab identity $\delta^*_W=\operatorname{Var}_p(X)(ER^*-D^*)$ still holds with the district-fitted $\tau$, so the CR slab rows remain valid rows under every rung.)*

*Proof.*
1. For each rung, $m_i=m_{\rm id}(V_i)+L_ik$ with $L_i$ fixed by the design. [2a L1; 2b Lemma D1; Def of the rung]
2. Each association is linear in $(m_i)_i$ with design coefficients. [L M1(d) with $m$ for $\beta$] Composing with 1, it is affine in $k$. (a) follows. [Def 2]
3. (b): $\delta^*_B(r,s)=\sum_ip_iX_{ir}X_{is}\big(m_r(X_i,D_i)-m_s(X_i,D_i)\big)$, and every term is $\ge0$ by hypothesis, since $p_iX_{ir}X_{is}\ge0$. [Def; hypothesis]
4. (c), M0: $\operatorname{Cov}_{w^1}(b_0+b_1X,X)=b_1\operatorname{Var}_{w^1}(X)$, and likewise for row 2. Then $b_1=w_1+d_1$ (P1 §1.3 dictionary). [Def Cov; alg]
5. (c), district rungs: apply the law of total covariance under $w^s$, conditioning on $D$. Within district $d$, $m_s(\cdot,d)$ is linear with slope $b_{s,d}$, so the conditional covariance is $b_{s,d}\operatorname{Var}_{w^s}(X\mid D=d)\ge0$ when $b_{s,d}\ge0$. [total covariance; 4] The example is computed in C9. ∎

Checked: C4 (affinity in $k$ for R = C = 3; ordering ⇒ $\delta^*_B\ge0$, with valid transition matrices; *rebuilt after review, Batch B: the old construction changed one column without rebalancing the others*), C9 (district reversal and the decomposition), and 1b V3 (closed form in 2×2).

**Remark M5.1 (which δ).**
* N1 uses the realised associations of the actual population (owner decision 1).
* K-rows use $\delta^*$, the behavioural version.
* Under (S2), $\delta-\delta^*=\sum_i\sum_rc_{ir}e_{ir}$ is a design-weighted sum of unit residuals, with mean zero given the design. Under P3-S1a (independence given the design) its standard deviation is at most $\tfrac12\{\sum_i(\sum_r\lvert c_{ir}\rvert)^2\}^{1/2}$, the root sum of squared design coefficients; it is small only when no unit keeps a non-negligible coefficient. (S2) alone, a mean condition, gives no rate. *(Corrected after review, Batch B: earlier "order $p^{-1/2}$", which needs these sampling and weight conditions.)*
* Intersecting N1 with an LCE rung therefore mixes the two readings. The package will report them separately.

---

## 4. Monotone nonparametric contextual effects (rung N2; two rows)

**Support convention** *(added after review, Batch B)*. N2 is defined on the **finite empirical design law**: $\operatorname{supp}P_X=\{x_1<\dots<x_{J_X}\}$ is the set of distinct observed compositions, each with positive mass. $J_X$ denotes the support size; polynomial degrees in §6 are written $q$. On a continuous population support the observed compositions are only a grid, and the results below do not carry over (N2.2(b)).

**Model N2(σ).**
* (S1) holds, and $E[\beta\mid V]=m(X)$ for a Borel $m:[0,1]\to\mathfrak B$ (A3). Here $X$ is the row-1 share.
* For designated pairs $(r,c)$ with directions $\sigma_{rc}\in\{+1,-1\}$, $\sigma_{rc}m_{rc}$ is nondecreasing on $\operatorname{supp}P_X$.
* No linearity is assumed. $\tau(x):=E[T\mid X=x]$ is identified at each support point.

**Lemma N2.1.** $m(x)\in S(x,\tau(x))$ for $P_X$-a.e. $x$, i.e. at every support point.
*Proof.*
1. $m(X)=E[\beta\mid X]\in\mathfrak B$ a.s. [CE-tower with A3; CE-cvx]
2. $m(X)^\top X=E[\beta^\top X\mid X]=E[T\mid X]=\tau(X)$ a.s. [CE-pull; CE-lin; A1] ∎

**Outer set and LP.** With the support points $x_1<\dots<x_{J_X}$,
$$\mathcal M_O:=\{(m(x_j))_j:\ m(x_j)\in S(x_j,\tau(x_j))\ \forall j;\ \sigma_{rc}\big(m_{rc}(x_{j+1})-m_{rc}(x_j)\big)\ge0\ \forall j,(r,c)\}.$$
Bounds on $\tilde B_{rc}=\sum_iw^r_im_{rc}(X_i)$ are the LP $[\min,\max]_{\mathcal M_O}$, which has $2CJ_X$ variables.

**Proposition N2.2.** Under the support convention:
(a) Every $m$ allowed by N2(σ) restricts to a point of $\mathcal M_O$. So the LP bounds are an **outer** set for $\tilde B$. They are not claimed sharp. Replacing $S(x,\tau(x))$ by $\mathbb E[S(X,T)\mid X=x]$ is **not** the sharp version when size is informative, because the model conditions on $V=(X,N)$ (Remark N2.2.1).
(b) Conversely, every point of $\mathcal M_O$ extends to a monotone Borel $m$ on $[0,1]$ satisfying Lemma N2.1 on the support: a common piecewise-constant extension of the whole matrix preserves every row sum and every designated direction. So, **on the complete finite support**, restricting to support points loses nothing. On a strict subset of the support (a coarser grid) this fails: extendibility subject to the **other** accounting equations can fail, and the grid LP is only a relaxation (Remark N2.2.2).
(c) **Nesting with LCE.** If $\tau$ is the quadratic of 2a Lemma 1, then $\{M_{Q+K}:K\in\Theta_O\cap\{\sigma\text{-signed slopes}\}\}\subseteq\mathcal M_O$. Hence $\mathcal B(\text{M0}\cap\text{signed slopes})\subseteq\mathcal B(\text{N2})$, for both outer sets.
(d) **Nesting with N1** (superpopulation reading). If $m_{1c}$ and $m_{2c}$ are both nondecreasing, then $\delta^*_W\ge0$ for column $c$. The two-row shape restriction is stronger than a pooled covariance sign: a sign on $\delta^*_W$ permits reversals.
(e) For $R\ge3$ a direction "in $x_1$" needs a **simplex chart**: which share adjusts when $x_1$ changes. With the chart $(x_1,x_3)$ (so $x_2$ adjusts), $m_{rc}(x)=x_3$ is constant, hence weakly nondecreasing, in $x_1$; yet under a symmetric Dirichlet design $E[X_3\mid X_1]=(1-X_1)/2$, and tilting by $X_1$ or $1-X_1$ preserves that conditional mean, so the covariance with $X_1$ is negative under **both** W1 weights. So (d) has no general R×C analogue. (A globally increasing *univariate* function of $X_1$ still obeys the covariance inequality.)

*Proof.*
1. (a): N2(σ) ⇒ Lemma N2.1 at each $x_j$, and monotonicity holds on the support. [L N2.1; Def] The target depends on $m$ only at the support points. [Def $\tilde B$]
2. (b): extend the whole matrix $m$ as a step function, constant between support points. [ext] The extension is Borel, each row stays in the simplex, each designated entry stays monotone, and its values at the support points are unchanged, so Lemma N2.1 still holds there. [Def]
3. (c): for $K$ in the outer set, $M_{Q+K}(x)\ge0$ on $A\supseteq\operatorname{supp}$ [2a L3(a)] and its rows sum to 1 [2a L1(b)], so it lies in $\mathfrak B$. Also $x^\top M_{Q+K}(x)=x^\top Qx=\tau(x)$ [2a L1(a)], so $M_{Q+K}(x)\in S(x,\tau(x))$ at each atom. A linear function with a slope of sign σ is σ-monotone. [Def]
4. (d): $\delta^*_W=\mu\operatorname{Cov}_{w^1}(m_1(X),X)+(1-\mu)\operatorname{Cov}_{w^2}(m_2(X),X)$, with both covariances $\ge0$. [Def 1 with $m$; Cheb applied under $w^1$ and under $w^2$]
5. (e): the example. [C5(c), C5(c′)] ∎

*(Corrected after review, Batch B, Nov 2026: "support points (the distinct observed compositions)" was ambiguous between the finite empirical support and a grid on a larger support, and (b) is false on a grid; (a)'s proposed sharp version was not sharp; (e) lacked a chart and the weighted covariances. Checks C10, C10b, C5(c′).)*

**Remark N2.2.1 (size conditioning; Conjecture).** Fix $X=1/2$ and give two sizes equal probability. At the first size $T$ is 0 or 1 equiprobably, at the second $T=1/2$ surely. At the first size both row rates must equal $T$, so the model forces $m(1/2)=(.5,.5)$. Conditioning only on $X$ (and the $S(x,\tau)$ LP) admits $(.75,.25)$, by selecting $(1,0)$ at the second size (check C10b). **Conjecture N2.2′.** The sharp N2 set is $\{m:\ m(X)\in\mathbb E[S(X,T)\mid V]\ P_V\text{-a.s., plus the shape restrictions}\}$, i.e. every conditional support-function inequality holds for the same $m(x)$ at almost every size given $x$. The proof route is 2a Theorem 1's selection construction applied to each admissible measurable $m$; it is not written. Conditioning on $X$ alone suffices if size carries no information relevant to the selection sets, or if the model conditions only on $X$. The same caution applies to M-SHARP (PD step 1b): its series bound conditions on $X$ only, so it is "tighter than $\Theta_O$", not sharp, when size is informative.

**Remark N2.2.2 (grids).** True support $\{1/4,1/2,3/4\}$, row-1 means $(.4,.4,.4)$, row-2 means $(8/15,.7,.8)$, both weakly increasing, $\tau=(.5,.55,.5)$. On a grid keeping only the extremes, assigning both rows $.5$ at both extremes passes the LP. Any monotone extension then assigns $.5$ to both rows at the middle and cannot give $.55$ (check C10). For binned or gridded data use Lemma N2.4.

Checked:
* C5(a): N2 contains the truth under a monotone nonlinear DGP.
* C5(b): the nesting (c) held in 47 random cases.
* C5(c), C5(c′): the R = 3 example, with negative covariance unweighted and under both W1 weights.
* C10, C10b: the grid and size examples.
* 1b V5: Chebyshev.

**Proposition N2.3 (structural local identity; EGGHM Props 12–13 in this setting).** Let $m_{1c},m_{2c}$ be differentiable at $x$ and write $\tau_{m,c}(z)=z\,m_{1c}(z)+(1-z)m_{2c}(z)$ for the **structural** aggregate. Write $D_c(x)=m_{1c}(x)-m_{2c}(x)$ and $\delta_{W,c}(x)=x\,m_{1c}'(x)+(1-x)m_{2c}'(x)$. Then:
(a) $\tau_{m,c}'(x)=D_c(x)+\delta_{W,c}(x)$.
(b) $\delta_{W,c}(x)\ge0\iff D_c(x)\le\tau'_{m,c}(x)$.
(c) Local CR, $D_c(x)\delta_{W,c}(x)\ge0$, holds iff $D_c(x)$ lies between 0 and $\tau'_{m,c}(x)$ (including a zero derivative).
(d) Both rows nondecreasing ⇒ local sign $\ge0$. The local sign everywhere does **not** imply $\delta^*_W\ge0$, even for smooth $m$: with $X$ uniform on $[1/5,2/5]$, $m_1=x-x^2/2$ and $m_0=1-x^2/2$, both are probabilities, $xm_1'+(1-x)m_0'=0$ everywhere, and $\delta^*_W=-1/47250$ ($-1/5250$ for the finite design with equal mass at the two endpoints). The reverse non-implication also holds: complement both functions, shrink them slightly into the interior, and subtract a small positive multiple of $x$ from each; the aggregate within association stays positive by continuity and the local one becomes negative.

**Not an observable restriction under the support convention.** On a finite support $\tau'$ is **not identified**: adding $h(x)=\prod_j(x-x_j)$ to both rows leaves the law on the support unchanged and moves the derivative at every node. Example: support $\{1/4,1/2,3/4\}$, $m_1=.7$, $m_0=.3$, $h=(x-\tfrac14)(x-\tfrac12)(x-\tfrac34)$; both models are valid probabilities on $[0,1]$ and induce the same law, but $\tau'(1/2)$ is $.4$ in one and $.3375$ in the other. The first model satisfies weak local CR with gap $.4$; an LP using the other compatible derivative caps the gap at $.3375$ and excludes the truth (check C11). An observable version needs an **identified** derivative: an interior point of an interval support with compatible smooth versions, or a maintained identified regression family, and an estimated derivative adds sampling error. At a pure composition EGGHM's local between association is $x(1-x)D(x)=0$; replacing its sign by the sign of $D(x)$ changes the belief.

**Decision U4(a):** the N2-local rows are **dropped** from the build; N2.3 is kept as the structural calculus identity above. *(Corrected after review, Batch B, Nov 2026: N2.3 was stated as an observable restriction "linear in $m(x)$ given $\tau'$", with an estimated $\tau'$ at the support points; and (d)'s witness was the kinked 1b V6. Checks C11, C13; C6 checks the identity.)*

*Proof.*
1. Differentiate $\tau_{m,c}(z)=z\,m_{1c}(z)+(1-z)m_{2c}(z)$. [product rule] This gives (a), and (b) follows. [alg]
2. (c), case $\tau'\ge0$: if $D\ge0$ and $\delta_W\ge0$, then $0\le D=\tau'-\delta_W\le\tau'$. If $D\le0$ and $\delta_W\le0$, then $D=\tau'-\delta_W\ge\tau'\ge0\ge D$, so $D=0$. Conversely, $D\in[0,\tau']$ gives $D\ge0$ and $\delta_W=\tau'-D\ge0$. The case $\tau'\le0$ is symmetric. [a; alg]
3. (d): $x,1-x\ge0$; the example is computed exactly in C13. [Def] ∎

**Lemma N2.4 (person-weighted bins; added after review, Batch B, following the referee's repair).** Let $b=1,\dots,K$ be ordered, nonoverlapping bins of compositions, $P(b)>0$ their person-weighted masses and $E_b$ the person-weighted ($p$-weighted) expectation within bin $b$. Put
$$\bar x_b=E_bX,\qquad\bar m_{1b}=E_b[Xm_1(X)]/\bar x_b,\qquad\bar m_{2b}=E_b[(1-X)m_2(X)]/(1-\bar x_b),\qquad\bar\tau_b=E_bT,$$
with a compatible zero-weight extension for a group absent from a bin ($\bar x_b\in\{0,1\}$: any value that keeps the bin rates monotone). Under N2(σ):
(i) $\bar m_b\in S(\bar x_b,\bar\tau_b)$;
(ii) $\tilde B_{1c}=\sum_bP(b)\bar x_b\bar m_{1bc}/\mu$, and likewise for row 2: the row targets are linear in the bin rates with known weights;
(iii) if $\sigma_{rc}m_{rc}$ is nondecreasing, so is $\sigma_{rc}\bar m_{rbc}$ across the ordered bins.
Hence the binned LP $\{\bar m_b\in S(\bar x_b,\bar\tau_b),\ \text{adjacent order}\}$ with the targets of (ii) is an **outer relaxation** of the global N2 bounds. It does not cover local derivative rows.
*Proof.* (i) $\bar x_b\bar m_{1b}+(1-\bar x_b)\bar m_{2b}=E_b[Xm_1+(1-X)m_2]=E_b\tau=E_bT$ [L N2.1; CE-tower], each $\bar m_{rb}$ is a weighted average of simplex points, so lies in the simplex. (ii) $E_p[Xm_1]=\sum_bP(b)E_b[Xm_1]$. [lin] (iii) $\bar m_{1b}$ is an average of $m_1$ over bin $b$ under the tilted weights $\propto pX$, so lies in $[\min_bm_1,\max_bm_1]$; ordered nonoverlapping bins and monotone $m_1$ give $\max_bm_1\le\min_{b+1}m_1$. Row 2 likewise with weights $\propto p(1-X)$. ∎
The ordinary bin average does **not** work: with equal bin weights at $x=.25,.75$, $m_1=x$, $m_2=0$, $E_bT=E[X^2]=.3125$, while the vector of ordinary bin-average rates evaluated at $E[X]=.5$ gives $.25$; the person-weighted rate is $.625$ and reproduces $.3125$ (check C12). The earlier argument ("averaging Lemma N2.1 within a bin uses [CE-cvx]") applied conditional convexity to the varying sets $S(x,\tau(x))$, which does not prove the step.

---

## 5. The ladder (2 rows), with status

All N1 statements are for continuous tables (Remark M2.1).
$$DD\ \supseteq\ \text{N1}(\delta_B\text{ sign})\ ,\ \text{N1(CR)}\ ,\ \text{N1}(\delta_W\text{ sign})\ \supseteq\ \text{N1(both signs, same direction)}\qquad\text{[M2; nest]}$$
"Both signs" means both weakly positive or both weakly negative; then N1(both signs) ⊆ N1(CR). Two opposite signs generally do not imply CR.
$$\mathcal B(\text{M0}\cap\text{signed slopes})\ \subseteq\ \mathcal B(\text{N2})\quad(\text{outer, }\tau\text{ quadratic, finite support})\qquad\text{[N2.2(c)]}$$
$$\text{N2(both rows }\uparrow)\ \Rightarrow\ \delta^*_W\ge0\qquad\text{[N2.2(d)]};\qquad \text{mean-level order in each district}\ \Rightarrow\ \delta^*_B\ge0\qquad\text{[M5(b)]};\qquad\text{M0 signed slopes}\ \Rightarrow\ \delta^*_W\ge0\qquad\text{[M5(c), M0 only]}.$$
* N2 admits step functions, where a local derivative is undefined; the structural local identity N2.3 applies only to a differentiable subclass and is not an observable rung (U4(a)).
* Compare a realised association ($\delta$, N1) with a model-implied one ($\delta^*$) only after fixing the same population of tables or means (Remark M5.1).
* For $R\ge3$, any monotonicity statement needs a simplex chart (N2.2(e)).

**Not nested:**
* M0 with N1 or N2;
* the local sign everywhere and the aggregate $\delta^*_W$ sign, in either direction (N2.3(d));
* for $R\ge3$, N2 with N1($\delta_W$) (N2.2(e));
* pooled $\delta^*_W$ and signed slopes under the district rungs (M5(c)).

*(Corrected after review, Batch B, Nov 2026: the line "N2(both rows ↑) ⊆ N2-local" is removed with the N2-local rung; "both signs" now requires the same direction; the M5 cells are scoped to M0.)*

---

## 6. What Phase 3 must supply

* **N1** needs no sampling inference (owner decision 1). Phase 3 adds the slack dials ε and the breakdown frontier, with a scalar breakdown value along a stated path (M2(d)). Labels: "sharp for continuous unit tables given the margins".
* **K-rows:** the existing P3 machinery, with one extra row each (pooled M0; grouped rungs: R16).
* **N2 (global only; N2-local dropped, U4(a))** needs $\tau$ at the support points.
  * A weighted polynomial of **prespecified** degree $q$ (separate notation from the support size $J_X$). Degree $q$ has $q+1$ coefficients and is identified only from at least $q+1$ distinct support points; $q=J_X$ is never identified. With $q=2$, N2.2(c) holds exactly.
  * The constraints are then linear in the coefficients, and P3.5 applies *conditional on the degree* ((a) and (b2) proved; CLR as implemented algorithmic, (b3)): valid if $\tau$ is a polynomial of degree $q$, with approximation bias otherwise.
  * **Do not select the degree by cross-validation and then use fixed-degree inference.** In a balanced two-level design with true slope zero, choosing between intercept-only and intercept-plus-slope by leave-one-out error selects the slope when $Z^2>2$, and the nominal 95% slope interval then covers with conditional probability $1-.05/P(\lvert Z\rvert>\sqrt2)\approx0.68$ (check C12b). Use a prespecified degree, sample splitting, or a selection-aware procedure, and state any approximation bias.
  * With many support points, group compositions into bins and use the **person-weighted** bin rates of Lemma N2.4; $\bar\tau_b$ is then the person-weighted bin mean of $T$. This is an outer relaxation of global N2.
  * *(Corrected after review, Batch B, Nov 2026: earlier "a weighted polynomial of degree $J$, selected by cross-validation" (with $J$ also the support size) and "binning keeps the outer property: averaging Lemma N2.1 within a bin uses [CE-cvx]". Checks C12, C12b.)*
* **Implementation sizes.** N1 has $p\,R\,C$ variables. That is about 370×6 for Wisconsin in 2×3, which `lpSolve` handles; the equalities will be eliminated to $p(R-1)(C-1)$ free variables.
