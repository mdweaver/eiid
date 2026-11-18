# Phase 5a: the sign module (M-SIGN), proofs for build step 1

*Proof-register items R3, R2 and R10-union from `PACKAGE_DESIGN.md` §8. Notation and justification tags follow 2a §0.3, 2c §0.3 and P3 §0.3. Each item ends with its numerical check. Status: R3 (§1), R2 (§2) and R10-union (§3) written and checked. An independent review pass (fresh session) found one error, in R10(d), and several gaps; all are fixed below and listed in §4. Two external referee reports (Batch B, Nov 2026: `reviews/B_partial_identification_report.md`, `reviews/B_ecological_inference_report.md`) found further errors and gaps, ledgered in `reviews/B_PI_RESPONSE_PLAN.md` (rows P-B15 to P-B27, P-B36 to P-B38, P-B43; P-E6, P-E12, P-E17, P-E19, P-E21, P-E22, P-E28). Each correction below is marked "Corrected after review, Batch B". **Scope:** every statement here is for the **pooled** LCE model (M0). The sign and CR rows are valid linear rows under the district rungs too (2c M5(a), as corrected), but R3 and R10 for M1, M1(δ), M2 and L4 are proof item **R16**, not proved.*

*Extra tags used here:*

| Tag | Meaning |
|---|---|
| [Weier] | a continuous function on a nonempty compact set attains its maximum |
| [pert] | for an LP, the perturbation function $v(\beta)=\max\{\varphi(k):k\in F,\ a^\top k\le\beta\}$ is concave and nondecreasing on its domain. Its superdifferential at $\beta$ is the set of optimal dual multipliers of the row $a^\top k\le\beta$, and for a concave function on $\mathbb R$ it is $[v'_+(\beta),v'_-(\beta)]$ (Rockafellar 1970, Thms 23.2, 24.1, 29.1†; Bertsimas & Tsitsiklis 1997, §5.2†) |

---

## 1. R3: sharpness of LCE plus restrictions on $k$

### 1.1 Setting

* The design $(N_i,X_i)_{i\le p}$ is fixed, and $A\supseteq\operatorname{supp}P_X$ is a polytope (2a §0).
* LCE rung, **pooled** (M0): $H=Q(P)+K$ with $K\in\mathcal K$, and $k\in\mathbb R^d$ its free coordinates (2a Prop. 7). The district rungs are not covered (proof item R16; 2c M5).
* $\Theta_I=\Theta_I(P;A)$ is the sharp set (2a Thm 1). $\Theta_O=\Theta_O(P;A)$ is the regression-only set that the LP computes (2a Prop. 7). Both are treated as subsets of $\mathbb R^d$ through $k$.
* **Target.** A model-implied linear target $\varphi(k)=\varphi_0+h^\top k$, for example $\tilde B_{rc}$ or $\tilde D=\tilde B_1-\tilde B_2$ (2a Prop. 7; 2c §0.2). (The letter $\gamma$ is reserved for the variance ratio $\gamma_r$ of (d).)
  * $\tilde B$ uses the design weights $N_i$. Fixing the design permits the auxiliary **convention** of giving its units equal probability, so that the finite-design target is represented as 2a's $B^*$ for that empirical law. Conditioning alone does not equate this target with an arbitrary superpopulation rate. *(Corrected after review, Batch B, Nov 2026: earlier worded as an identity.)* In the population, $ER$ computed from $T$ and from $\tau=E[T\mid X]$ coincide by iterated expectations, so restrictions written with the $\tau$-based $ER^*$ satisfy (C2).
  * The actual-population rate $B$ is *not* covered: $B-\tilde B$ is a weighted sum of unit residuals that the data never pin down, and it needs the P3 allowance.
* **Admissible restriction.** A set $C=C(P)\subseteq\mathbb R^d$ such that
  * (C1) $C$ is closed and convex;
  * (C2) $C$ depends on $P$ only through the design and $Q(P)$.

  (C2) is stronger than (a) needs; (b) uses it, because there the restriction must be the same for every law with the same regression and design. (C1)–(C2) admit non-polyhedral sets (a ball, for instance), so sharpness alone does not make the computation a finite LP. *(Corrected after review, Batch B, Nov 2026.)*

  The examples below are each a finite set of rows $\{Gk\le g\}$, so for them every computation is a finite LP:
  * any sign on $\delta^*_B(r,s)$, $\delta^{r*}_B$, $\Delta^*_r$ or $\delta^{r*}_W$ (2c M5(a), pooled form);
  * mean-level ordering on $A$ (2c M5(b));
  * the M-LOC rows $(x^\top H)_r-(x^\top H)_s\ge0$ on $A$, which reduce to the vertices because they are linear in $x$ (2a L3(a)).
* **Union restriction.** $C=C_+\cup C_-$ with $C_\pm$ admissible. Contextual reinforcement for a pair $(g_B,g_W)$ of affine functions of $k$ is $C_\pm=\{\pm g_B\ge0,\ \pm g_W\ge0\}$.

### 1.2 Statement

**Theorem R3.** Assume 2a (A1)–(A5), standing integrability, $P(X\in A)=1$, and that $\operatorname{affdim}(A)=R-1$. *(Corrected in the review response: "nonempty relative interior" is too weak, see 2a Remark 5.2. Under A5 the condition holds automatically.)* The dimension condition is sufficient, not necessary: $R$ affinely independent points of $A$ bound every entry of $K$ (2a L3(b)). $P(X\in A)=1$ is used twice: the deterministic witness of (b) and the full target both use the region on which the mean is specified. With trimming ($P(X\in A)<1$) the outside-$A$ terms need 2a Prop. 7′.

(a) **Sharp interval.** For admissible $C$, the sharp identified set for $\varphi$ under LCE $\cap\,C$ is
$$\varphi(\Theta_I\cap C)=\Big[\min_{\Theta_I\cap C}\varphi,\ \max_{\Theta_I\cap C}\varphi\Big],$$
or $\emptyset$.

(b) **What the LP delivers.** $\Theta_O\cap C=\bigcup_{P'\in\mathcal P(P)}\big(\Theta_I(P')\cap C(P')\big)$, where $\mathcal P(P)$ is as in 2a Prop. 2.
* So the interval $\varphi(\Theta_O\cap C)$ is sharp for the information in the regression $E[T\mid X]$ and the design. It is a finite LP when $C$ is polyhedral (all package restrictions are). *(Corrected after review, Batch B, Nov 2026: earlier called "the LP interval" for every admissible $C$.)*
* It equals the sharp interval of (a) whenever $\Theta_I=\Theta_O$ (2a Cor. 2, Prop. 3).
* Otherwise it is a weakly outer interval, and M-SHARP can tighten it.
* **Remark R3.1 (full law versus regression; added after review, Batch B).** Let $X=1/4,1/2,3/4$ with equal probability and equal sizes, and $T$ a fair binary total at every composition. The full law forces both row rates to equal $T$ in every unit, so both row means are $.5$ everywhere and $k=0$. The regression alone ($\tau=.5$) admits $m_1=.5+k(1-x)$, $m_0=.5-kx$. For $C=\{k\in[.2,.4]\}$ the full-law identified set (a) is **empty**, while the regression-only gap set (b) is $[1/6,1/3]$. This is exactly the distinction (a)–(b) draw. Example due to the referee; check S6.

(c) **Unions.** For $C=C_+\cup C_-$, the identified set is $\varphi(\Theta\cap C_+)\cup\varphi(\Theta\cap C_-)$ for $\Theta\in\{\Theta_I,\Theta_O\}$. That is at most two closed intervals.
* If one piece is empty, the set is the other piece, a single interval, and the hull is sharp.
* If both are non-empty, the hull $[\min,\max]$ is sharp exactly when the two intervals intersect.
* **Certificate:** if $\Theta\cap C_+\cap C_-\neq\emptyset$, they intersect. For $\Theta_O$ with polyhedral $C_\pm$ this is one feasibility LP. The certificate is sufficient only (§1.4).

(d) **One-vs-rest CR is one LP.** Let $(g_B,g_W)=(\delta^{r*}_B,\delta^{r*}_W)$ (W1), and assume $\operatorname{Var}_p(X_r)>0$. Put:
* $\gamma_r=\operatorname{Var}_p(X_r)/(\mu_r(1-\mu_r))\in(0,1]$;
* $ER^{r*}=\operatorname{Cov}_p(X_r,\tau(X))/\operatorname{Var}_p(X_r)$, with $\tau(x)=(x^\top Q_cx)_c$;
* $D^{r*}=\tilde B_r-\tilde B_{\neg r}$.

Then
$$C_+\cup C_-=\big\{k:\ D^{r*}(k)\in\operatorname{hull}\{\gamma_rER^{r*},\ ER^{r*}\}\big\},$$
which is two linear rows in $k$. For $R=2$ this is the only CR, so in $2\times2$ CR is a single LP. Under it, $\operatorname{sign}\tilde D=\operatorname{sign}ER^*$. Ties ($ER^{r*}=0$, or $\gamma_r=1$) can leave both branches feasible; the slab is unaffected. The slab can have an empty intersection with the LCE set (Remark R3.2).

(e) **Refutation is joint.** If $P$ is generated by a model in LCE $\cap\,C$ with skew part $K_0$, then $K_0\in\Theta_I\cap C\subseteq\Theta_O\cap C$. An empty population LP therefore rejects the **conjunction** LCE $\cap\,C$; emptiness alone does not attribute the failure to either part.
* A belief that also has a model-free meaning can sometimes be tested on its own. Example: on $x=1/4,1/2,3/4$ let both row means be $1.5-2x$ and the aggregate equal this function. LCE holds on $A=[1/4,3/4]$, $ER=-2$, and every feasible gap is at least $-1$, so a positive within association is impossible even in the model-free continuous table problem. For a restriction defined only on the LCE parameter $k$, the meaning of $C$ already presumes LCE.
* (e) is a population statement. Infeasibility of the LP at $\hat\beta$ is a sample diagnostic; a level-controlled rejection uses R10(d)'s rule.

*(Corrected after review, Batch B, Nov 2026: earlier "never $C$ alone". The referee's examples show both directions: a belief can be refutable without LCE, and LCE can exclude a belief that is never refutable model-free (Remark R3.2).)*

### 1.3 Proof

*(a)*
1. A value $t$ is in the identified set iff some $\mathbb Q\in\mathcal M(Q(P)+K;A)$ with $\mathbb Q\circ W^{-1}=P$ has $k\in C$ and $\varphi(k)=t$. [Def identified set; Def $\Theta_I$]
2. $C$ is evaluated at the true $P$ for every such $\mathbb Q$, because it depends on $\mathbb Q$ only through its $W$-marginal $P$, via the design and $Q(P)$. So the condition in 1 is $k\in\Theta_I\cap C$. [C2; 1]
3. The identified set is therefore the image $\varphi(\Theta_I\cap C)$. [2]
4. $\Theta_I$ is convex and compact. $C$ is closed and convex. So $\Theta_I\cap C$ is convex and compact, or empty. [2a T1; C1]
5. $\varphi$ is affine and continuous, so the image of a nonempty compact convex set is a compact interval with the stated endpoints. [lin; IVT] ∎(a)

*(b)*
6. ⊇: for $P'\in\mathcal P(P)$, $Q(P')=Q(P)$ and $P'_{(N,X)}=P_{(N,X)}$, so $C(P')=C(P)$. Also $\Theta_I(P')\subseteq\Theta_O(P')=\Theta_O(P)$. [Def $\mathcal P$; C2; 2a P1, P2 step 2]
7. ⊆: take $k\in\Theta_O\cap C$. The construction in 2a Prop. 2 (steps 4–10) gives $P'\in\mathcal P(P)$ with $k\in\Theta_I(P')$. Then $k\in C(P)=C(P')$. [2a P2; 6]
8. The statements about the LP interval follow from (a) applied to each $P'$, and from $\Theta_I=\Theta_O$ in the cases cited. [7; 2a Cor. 2, Prop. 3] ∎(b)

*(c)*
9. $\varphi(\Theta\cap(C_+\cup C_-))=\varphi(\Theta\cap C_+)\cup\varphi(\Theta\cap C_-)$. Each piece is a closed interval or empty, by (a), or by (b) for $\Theta_O$. [image of a union; (a), (b)]
10. If one piece is empty the union is the other piece. A union of two non-empty closed intervals is an interval iff they intersect. A point $k_0\in\Theta\cap C_+\cap C_-$ gives $\varphi(k_0)$ in both. [order on $\mathbb R$] ∎(c)

*(d)*
11. Apply 2c Lemma M1(c) with $m(X_i)=M_H(X_i)$ in place of $\beta_i$, and $\tau(X_i)$ in place of $T_i$. The identities there use only $T=\sum_sX_s\beta_s$ and the simplex, and $x^\top M_H(x)=\tau(x)$ holds. This gives
$$\delta^{r*}_W=\operatorname{Var}_p(X_r)\,(ER^{r*}-D^{r*}),\qquad \delta^{r*}_B=\mu_r(1-\mu_r)\,(D^{r*}-NM^{r*}),$$
where $NM^{r*}$ is the neighbourhood gap computed from $\tau$. [2c L M1(c); 2a L1(a); alg]
12. $NM^{r*}=\operatorname{Cov}_p(X_r,\tau)/(\mu_r(1-\mu_r))=\gamma_rER^{r*}$. [Def; alg; 2c M3(d) for the collapsed population]
13. Both $\ge0$ ⇔ $\gamma_rER^{r*}\le D^{r*}\le ER^{r*}$. Both $\le0$ ⇔ $ER^{r*}\le D^{r*}\le\gamma_rER^{r*}$. The union of these is $D^{r*}\in\operatorname{hull}\{\gamma_rER^{r*},ER^{r*}\}$. [11–12; positive factors]
14. $D^{r*}$ is affine in $k$, and $ER^{r*}$ and $\gamma_r$ depend only on the design and $Q$. So the set is two rows of an admissible $C$. The sign statement follows because $\gamma_r>0$. [2c M5(a); C2; 13] ∎(d)

*(e)*
15. $K_0\in\Theta_I$ by definition, and $k_0\in C$ by hypothesis. Then apply 2a P1. [Def; 2a P1] ∎

### 1.4 What fails, and where

* **Pairings other than W1.** For CR on $(\delta^{r*}_B,\Delta^*_r)$ or $(\delta^*_B(r,s),\Delta^*_r)$, (d) fails. The two sign conditions are on *different* affine functions of $k$, so both branches can be non-empty and their images can be disjoint (S5).
* **Not covered:**
  * N2 shape restrictions, which do not restrict $k$;
  * N1 signs on the realised $\delta$, which are unit-level, not mean-level;
  * any target that mixes $\tilde B$ with $B$ (marked ⚠ in the matrix).
* **The common-witness certificate is only sufficient** (example added after review, Batch B). In the square $[-1,1]^2$ with $g_B=u+.5$ and $g_W=u-.5$, the two branches are the disjoint strips $u\ge.5$ and $u\le-.5$. For the target $u$ the images are $[.5,1]$ and $[-1,-.5]$, a gap of length one. For the target $v$ both images are $[-1,1]$, identical, although the certificate LP $\{g_B=g_W=0\}$ is infeasible.
* **Remark R3.2 (LCE can exclude CR although LCE fits; added after review, Batch B).** Let $X=0,1/4,1$ with equal probability, $A=[0,1]$, $m_0=1-x$, $m_1=(49/55)x$, and deterministic $T$. The regression is $(1-x)^2+(49/55)x^2$. Its Bernstein cross coefficient is 0, so the vertex rows force $k=0$: LCE fits, with a single point. There $ER=0$ and $\tilde D=-261/1925$. CR needs $\tilde D\in\operatorname{hull}\{\gamma ER,ER\}=\{0\}$, so **both LCE branches are empty**, although LCE alone is not rejected. In the model-free continuous table problem (2c M3) CR is still feasible ($\beta_{i1}=\beta_{i2}=T_i$). Any "CR cannot be refuted" statement must name its base system. Example due to the referee; check S7.
* **Degenerate hull.** If $\operatorname{affdim}(A)<R-1$, $\Theta_O$ can be unbounded (2a L3(b), Remark 5.2). A singleton $A$ has nonempty relative interior and is still degenerate; under A5 this case cannot arise. For polyhedral $C$, (a) and (b) then hold with a possibly infinite endpoint, and the LP reports unboundedness. For non-polyhedral closed convex $C$ the image of an unbounded closed convex set need not be closed, so an endpoint may be a supremum that is not attained. Added rows can restore a finite endpoint (V5, near-1d support).

### 1.5 Numerical check (`verify_R3.py`, output `verify_R3_output.txt`)

| Check | Result |
|---|---|
| S1: constructive sharpness *for the LP set*, i.e. part (b) (it uses the 2a Prop. 2 construction with $x^\top\beta=\tau$). 3×2, LCE with the true signs of $\delta^{0*}_B$ and $\delta^{0*}_W$. For 1000 random $t\in[L,U]$, solve for $k$, build $\beta_i=M_H(X_i)$, and check simplex, accounting, rows and $\varphi=t$ | worst violation $1.8\times10^{-15}$ |
| S2: 2×2 CR branches (2000 designs). *Relabelled after review (Batch B): each branch is classified as empty, a point, or a positive-width interval; the old code called "no positive-width branch" "both degenerate", which also admitted two empty branches* | see output: no case with two positive-width branches and none with two empty branches. Frequencies on these draws are descriptive (the design generator is reused) |
| S3: 3×2 W1 CR, row 0 (1000 designs) | both branches non-empty 0 times, as (d) predicts |
| S4: 3×2 CR with the W2 or pairwise pairings | both non-empty in 357/1000 and 614/1000; no gaps |
| S6 (added after review, Batch B): full law versus regression, Remark R3.1 | regression gap set $[1/6,1/3]$; full-law set empty |
| S7 (added after review, Batch B): LCE excludes CR although LCE fits, Remark R3.2. *Fails on the old §3.5 "unrefutable" as applied there* | LCE: $k=0$; $ER=0$; $\tilde D=-261/1925$; both branches empty |
| S5: adversarial search (3 bases × 2 pairings × 3 targets; 3600 problems) | both non-empty 951; certificate held 945; **gaps 3**, all where the certificate failed; gap ≈ 2% of hull width on this seed. The reviewer's rerun on a fresh seed found 14 gaps, up to ≈ 31% of hull width, each confirmed by checking that both branch LPs are infeasible at the gap midpoint |

### 1.6 Consequences for the package

1. For the **pooled** linear rung (M0) plus any sign, ordering or local rows, `eiBounds()` may label the LP interval **sharp for the regression**. It is the sharp identified set when $\Theta_I=\Theta_O$, and outer otherwise. Under district rungs (M1, M1(δ), M2, L4) the rows are correct but sharpness is not proved (R16): the label should read "outer for the regression", and the note should say the sign is one pooled restriction across groups. *(Corrected after review, Batch B, Nov 2026.)*
2. One-vs-rest CR (the default, W1) and $2\times2$ CR are implemented as two rows in $k$, so identification needs no union. Inference is not quite ordinary: the rows' constants are a min/max of estimated quantities. The conservative CI inflates each piece before the max (R10(b)); the CLR CI takes the union of the two branch CIs (R10(c)). With the package's simulated candidates, the CLR version is algorithmic (Conjecture P3.5(b3)).
3. For other CR pairings the package:
   * runs the certificate LP. The certificate is sufficient, not necessary, for the two branch images to meet. In S5, 3 cases without a certificate had no gap, and every gap had a failing certificate. R3(c) allows a gap; the certificate is not a hypothesis;
   * if the certificate fails, solves both branches and reports both intervals when they are disjoint.
4. Emptiness messages always name the *combination* (LCE + restriction) that has been rejected, say whether the rung alone is rejected, and distinguish three events: an empty plug-in set at $\hat\beta$ (a sample diagnostic), an empty population set (R3(e)), and a rejection by the calibrated test (R10(d)). A statement that a belief "cannot be refuted" must name its base system (Remark R3.2). *(Corrected after review, Batch B, Nov 2026.)*

---

## 2. R2: when does an added restriction tighten the bound?

### 2.1 Setting

* $F=\{k:Gk\le g\}$ is the current feasible set: the rung, plus any rows already imposed. It is compact by 2a L3(b) when $\operatorname{affdim}(A)=R-1$, which A5 implies (2a Remark 5.2). R2's proof uses this compactness (step 3). Nonemptiness is *assumed*: R3 does not deliver it, and an empty $F$ is a rejection (R3(e)).
* $\varphi(k)=\varphi_0+h^\top k$ is a linear target, with $U=\max_F\varphi$, $L=\min_F\varphi$, and optimal faces $F^*_U=\arg\max_F\varphi$ and $F^*_L=\arg\min_F\varphi$.
* $N$ is the added restriction: an admissible set as in §1.1, closed and convex, with $F\cap N\ne\emptyset$. $U_N$ and $L_N$ are the endpoints after adding it.
* Three possible verdicts for $N$ relative to $F$:
  * **implied:** $F\subseteq N$;
  * **cuts the set only:** $F\not\subseteq N$, but $[L_N,U_N]=[L,U]$ for the target in question;
  * **tightens the target:** $[L_N,U_N]\subsetneq[L,U]$.

### 2.2 Statement

**Proposition R2.**

(a) **Face criterion.** $U_N<U$ iff $F^*_U\cap N=\emptyset$, and $L_N>L$ iff $F^*_L\cap N=\emptyset$ (given attainment on both nonempty sets, which compactness supplies). Each is one feasibility LP when $N$ is polyhedral; the LP must carry the target intercept (step 4). In words: the restriction moves an endpoint iff it removes *every* point that attains that endpoint.

(b) **One row $a^\top k\le b$.** $U_N<U$ iff $\min_{F^*_U}a^\top k>b$. This is also equivalent to: every optimal dual multiplier of the new row is $>0$ (that is, $v'_+(b)>0$).
* A positive multiplier returned by the solver at *one* optimum does **not** certify tightening. Example (added after review, Batch B): maximise $x$ on $[0,1]$ with the row $x\le1$ duplicated; every split $(1-\lambda,\lambda)$ of the unit dual cost is optimal, so the new row can carry multiplier 1 without moving the endpoint. A positive multiplier certifies tightening only if its uniqueness, or a positive lower bound on all optimal multipliers, has been established separately.
* A zero multiplier, or a row that is slack at an optimum, certifies that the row does not tighten.
* "Active" (the row holds with equality at the optimum) does **not** mean "tightening".

(c) **Implied = useless for every target.** $N$ is implied iff it tightens no linear target. For general closed convex $N$ the witness comes from a hyperplane separating a point of $F\setminus N$ from $N$; for polyhedral $N$, any non-implied row of $N$ serves. For a row $a^\top k\le b$, it is implied iff $\max_Fa^\top k\le b$, and when it is not implied the target $\varphi=a^\top k$ is a witness that tightens. So the "already implied" table answers "can this assumption buy anything at all?". Whether it buys something *for your target* is a separate question, answered by (a). The witness objective quantifies over all linear functions of $k$; it need not be a reported cell-rate target.

(d) **$2\times2$ under LCE.** Assume $\operatorname{Var}_p(X)>0$ and $\gamma<1$. Then $k=w_1$ is scalar and $\tilde D$ is strictly monotone in $w_1$, with slope $-(1-\gamma)$ times a positive factor.
* Each of the four sign rows is either implied or tightens $\tilde D$, and any target that is non-constant in $w_1$; for such targets the middle verdict cannot occur. A target constant in $w_1$ is never tightened. *(Corrected after review, Batch B, Nov 2026: the qualifications were implicit.)*
* The closed forms are those of G2:
  * $\delta^*_W\ge0$: $[L_J,\min(U_J,ER^*)]$; implied iff $ER^*\ge U_J$; empty iff $ER^*<L_J$.
  * $\delta^*_W\le0$: $[\max(L_J,ER^*),U_J]$.
  * $\delta^*_B\ge0$: $[\max(L_J,NM^*),U_J]$.
  * $\delta^*_B\le0$: $[L_J,\min(U_J,NM^*)]$.

  Here $ER^*$ and $NM^*$ are computed from $\tau$ (§1.2(d)), and $[L_J,U_J]$ is the JKST interval for $\tilde D$. The formulas are intersections of $[L_J,U_J]$ with a half-line; reversed endpoints mean the set is empty.

(e) **Several rows.** For $N=\bigcap_jN_j$, (a) applies to $N$ as a whole. A row's contribution depends on what is already imposed.
* The package reports, for each row:
  * its gain when added first (to the rung alone);
  * its gain when added last (dropped from the full stack).
* A row whose last-added gain is zero has **zero incremental gain for this target**, given the others, even if it tightens on its own. It need not be implied by the others, and it may still tighten other targets. Example: on $[0,1]^2$ with target $x$ and rows $x+y\le.8$, $x\le.5$, adding $x+y\le.8$ last gives zero gain for $x$, yet it still removes $(.5,1)$ and lowers the upper endpoint of $x+y$ from 1.5 to .8 (check T4b). A row is redundant **for the feasible set** only if the implication test (c) says so. *(Corrected after review, Batch B, Nov 2026: earlier "redundant given the others", which reads as "implied".)*

### 2.3 Proof

*(a)*
1. $F\cap N\subseteq F$, so $U_N\le U$. [nest]
2. If $k\in F^*_U\cap N$, then $U_N\ge\varphi(k)=U$. So $U_N=U$. [Def; 1]
3. If $F^*_U\cap N=\emptyset$: $F\cap N$ is nonempty and compact, so $\varphi$ attains $U_N$ at some $k'\in F\cap N$. Since $k'\notin F^*_U$, $\varphi(k')<U$. [C1; compact; Weier; Def $F^*_U$]
4. The lower endpoint is the same argument applied to $-\varphi$. If $N=\{Mk\le m\}$, then $F^*_U\cap N=\{Gk\le g,\ h^\top k\ge U-\varphi_0,\ Mk\le m\}$ is a feasibility LP. [sym; Def] ∎(a) *(Corrected after review, Batch B, Nov 2026: the row was printed as $h^\top k\ge U$, omitting the target intercept $\varphi_0$. On $[0,1]^2$ with $\varphi=2+x$ that row is $x\ge3$, which rejects the whole optimal face and so reports "moves" for every added row. The package's face test (`.rxc_face_moves`) subtracts the intercept and was not affected; `verify_R2.py` used zero intercepts and could not see the error. Check T1c.)*

*(b)*
5. $F^*_U$ is a nonempty compact face of $F$, so $\beta^*:=\min_{F^*_U}a^\top k$ is attained, and $F^*_U\cap\{a^\top k\le b\}=\emptyset$ iff $\beta^*>b$. The first claim follows from (a). [Weier; (a)]
6. Let $v(\beta)$ be the perturbation function. Then $v(\beta)=U$ for $\beta\ge\beta^*$, by (a) applied at $\beta$, and $v(\beta)<U$ for $\beta<\beta^*$. [(a); Def $v$]
7. If $b\ge\beta^*$, $v$ is constant on $[b,\infty)$, so $v'_+(b)=0$ and $0$ is an optimal multiplier. [6; pert]
8. If $b<\beta^*$: by concavity, $v'_+(b)\ge(v(\beta^*)-v(b))/(\beta^*-b)>0$. So every optimal multiplier, which lies in $[v'_+(b),v'_-(b)]$, is positive. At the left end of the domain, $b=\min_Fa^\top k$, the multiplier set is $[v'_+(b),\infty)$ ($v'_-=+\infty$); the equivalence is unaffected. [pert; 6]
9. A solver returns *one* multiplier from $\partial v(b)$, and when $\partial v(b)=[0,\lambda]$ with $\lambda>0$ it may return a positive one even though the row does not tighten. This is the degenerate case of a row through an optimal vertex: it is active, but it removes nothing at the optimum (T1). Suppose the row is slack, $a^\top k<b$, at an optimum $k$ of the *augmented* LP. Near $k$ the augmented set coincides with $F$, so $k$ is a local maximiser of the linear $\varphi$ on the convex $F$, hence a global one. Then $k\in F^*_U\cap N$, and by 2 the row does not tighten. [pert; local = global for concave maximisation; 2] ∎(b)

*(c)*
10. If $F\subseteq N$, then $F\cap N=F$, and no target changes. [Def]
11. If $a^\top k\le b$ is not implied, then $\max_Fa^\top k>b$. For $\varphi=a^\top k$, $U_N\le b<U$. [Def; 1]
11′. General closed convex $N$ with $F\not\subseteq N$: take $k_1\in F\setminus N$. Strict separation of the point $k_1$ from the closed convex set $N$ gives $(a,b)$ with $a^\top k\le b$ on $N$ and $a^\top k_1>b$. For $\varphi=a^\top k$: $U_N\le b<a^\top k_1\le U$. [sep. hyperplane; 11] ∎(c)

*(d)*
12. By 2c M3 / §1.2(d) with $R=2$: $\delta^*_W=\operatorname{Var}_p(X)(ER^*-\tilde D)$ and $\delta^*_B=\mu(1-\mu)(\tilde D-NM^*)$, with $\tilde D$ affine in $w_1$ with slope proportional to $-(1-\gamma)$, so non-constant iff $\gamma<1$. [§1.3 step 11; G2 step 3]
13. Each sign row is therefore a half-line in $\tilde D$, and $F$ maps to $[L_J,U_J]$. Intersecting gives the closed forms. [12; R3(b) and 2a §6.1, since $F$ is the LP set]
14. Since $\tilde D$ is a bijection between $F$ and $[L_J,U_J]$, cutting $F$ means cutting $[L_J,U_J]$. So "cuts the set only" is impossible for $\tilde D$, and also for any target that is non-constant in $w_1$. [12; bij] ∎(d)

*(e)* This is (a) applied to $N$ and to each sub-stack. Order dependence follows from (a), since $F^*_U$ changes as rows are added. ∎

### 2.4 Numerical check (`verify_R2.py`, output `verify_R2_output.txt`)

| Check | Result |
|---|---|
| T1: face test vs brute force. 3×2 LCE base, random targets; 400 random rows and 400 rows through the solver's optimal vertex | agree 800/800. Of the 400 active-but-redundant rows, the solver reported a **positive dual in 133** |
| T1b (added after review; random targets give a vertex face almost surely, so T1 never tests a larger face): target normal = a facet normal, so the optimal face is that facet; added rows cut through the facet's interior, cut all of it, or miss it | face test = brute force in every case (see output) |
| T2: implied ⇔ the witness target $a^\top k$ does not tighten (a sanity check, close to a tautology since $U_N=\min(U,b)$ for this target) | 333/333 |
| T3: 2×2, four sign rows, 2000 designs | closed forms = LP (error $1.7\times10^{-14}$); not implied ⇔ tightens $\tilde D$ in 7409/7409; "cuts set only" 0 |
| T4: 3×2, true-sign rows on $\delta^{0*}_B$ and $\delta^{0*}_W$, targets $\tilde B_0,\tilde B_1,\tilde B_2$ (1800 = 300 designs × 2 rows × 3 targets) | 147 implied (8%); **494 cuts set but not this target (28%)**; 1153 tightens (64%); 6 not classified, solver status recorded (base solve status 4, numerical difficulty). *Corrected after review, Batch B: unsuccessful base solves were skipped silently; the status is now printed. The referee's rerun on other library versions classified 153/494/1153 of 1800* |
| T1c (added after review, Batch B): target with intercept, $\varphi=2+x$ on $[0,1]^2$; rows $x+y\le1.2$ (keeps part of the optimal facet) and $x+y\le.8$ (removes it). *Fails on the old step 4* | old face LP ($h^\top k\ge U$) reports "moves" for both rows, wrong for the first; corrected LP ($h^\top k\ge U-\varphi_0$) agrees with brute force for both |
| T4b (added after review, Batch B): square, target $x$, rows $x+y\le.8$, $x\le.5$. *Fails on the old (e) sentence read as "implied"* | added-last gain 0 for $x$, 0.7 for $x+y$; not implied by the other row |

### 2.5 Consequences for the package

1. **The redundancy diagnostic uses the face test, not solver duals.** For each added module and each reported target, the package runs one feasibility LP per endpoint. Dual values are shown only as supplementary sensitivity information.
2. **Three verdicts per module, reported separately:**
   * "already implied by the data and the rung" (target-free; the "already implied" table);
   * "narrows the set but not this quantity";
   * "narrows this quantity by X".

   In $2\times2$ the middle verdict cannot occur for targets that are non-constant in $w_1$ (R2(d)); with three or more groups it is common (T4: 28%). Each verdict is relative to a nonempty convex base, a feasible added restriction, one target and one order (first or last), and the output should name all four.
3. **Added-first and added-last gains** are reported for stacks, so users can see which beliefs are doing the work. A zero added-last gain means "no gain for this target", not "implied" (R2(e)).
4. **Statistical caveat (feeds R10).** In samples, "implied" and "tightens" are estimated sharp-edged verdicts. Near a tie (e.g. $ER^*\approx U_J$), report the width change with uncertainty rather than a binary verdict. Inference on a width gain needs **joint** treatment of the base and augmented endpoints (and of any estimated threshold, such as $ER^*$); an interval for either endpoint alone is not enough. A one-sided test of $\max_Fa^\top k\le b$ is an LP-value test of the kind P3 Thm 5 already covers. *(Corrected after review, Batch B, Nov 2026: earlier "the width change with its confidence interval", without the joint requirement.)*

---

## 3. R10-union: inference with added rows, either/or restrictions, and near-empty sets

### 3.1 Setting

P3 §0, with the regime table of P3 §0.5: the asymptotic results are theorems for M0 in every regime and for M1/M2 with a fixed number of districts (R1); M1/M2 with growing districts (R3) is a conjecture **before** any question of candidate collection; R2 has only the finite-sample interval P3.8 and simulation findings. Everything below is for the pooled M0 model unless stated (R16 for the district rungs). *(Corrected after review, Batch B, Nov 2026: earlier "regimes R1/R3", which does not settle increasing dimension.)* $\hat\beta$ is the regression coefficient vector with $\sqrt n(\hat\beta-\beta_0)\Rightarrow N(0,\Sigma)$ and $\hat\Sigma$ consistent. The design is conditioned on (P3-S2).
* Base rows: $s_j(\beta)+a_j^\top k\ge0$ (P3 §1 step 2).
* Target: $\varphi=c(\beta)+h^\top k$.
* **Linear-in-$\beta$ added row:** $s_j(\beta)+a_j^\top k\ge0$, with $s_j$ linear in $\beta$ and $a_j$ fixed by the design. Every sign, ordering and local row in §1.1 has this form. The reason: $\delta^*$, $\tilde B$, $ER^*$ and the orderings are linear in $m=M_{Q+K}(X_i)$ with design weights, and $Q$ is linear in $\beta$ [2c M5 step 2; P3 §1].
* **Parallel disjunction:** a row of the form $\max_{l\le L}s_l(\beta)+a^\top k\ge0$, equivalently $\bigcup_l\{s_l(\beta)+a^\top k\ge0\}$. The pieces share the same $k$-normal $a$.
* **General union:** $C_+\cup C_-$, each branch a finite set of linear-in-$\beta$ rows.
* $\mathrm{CI}(\mathcal R)$ denotes the P3.5 interval (version (a) or (b)) computed with the base rows plus the row set $\mathcal R$, at total level $\alpha=\alpha_U+\alpha_L(+\alpha_E)$. **Emptiness convention**, one rule per version: for P3.5(a), $\mathrm{CI}=\emptyset$ iff the inflated set $\hat\Theta^+$ is empty; for P3.5(b), $\mathrm{CI}=\emptyset$ iff the computed ends cross, $\hat L>\hat U$. (Applying both rules to (b) and combining them by a union bound would double the **bound** in (d); the actual rejection probability need not double. *Reworded after review, Batch B.*)

### 3.2 Statement

**Theorem R10.** Assume P3.5's conditions *for the augmented row system* (base rows plus added rows), and that the restrictions hold at the truth, i.e. $k_0$ satisfies them. Pooled M0, or M1/M2 in regime R1 (P3 §0.5).

**Status by interval version** *(replaces the earlier status note; corrected after review, Batch B, Nov 2026, following the referee's version table)*:

| Interval version (P3.5) | What R10 claims | Status |
|---|---|---|
| (a) simultaneous row inflation (conservative) | (a)–(d), including the direct parallel slab (b), under a valid joint row-and-constant approximation for the design sequence | Theorem |
| (b1) CLR over all dual vertices | (a), (c), (d), under the **full augmented** CLR conditions, including one-sided validity. The added rows change the dual feasible set, its vertices, the gradients, the selection step and the quantile; appending row standard errors to the old CLR maximum is not enough | Theorem under those conditions (cited) |
| (b2)(i) fixed simultaneous superset | the same branch-union and crossing arguments, with the fixed-superset bounds and the allocated endpoint levels | Theorem |
| (b2)(ii) projection | the joint coefficient event controls both ends and any data-dependent collection of dual-feasible candidates; spend the coefficient-error budget once, and account separately for the residual of actual targets. Needs $\Sigma$ positive definite or a reduced full-rank parameterisation (P3.5(b2)(ii), corrected) | Theorem |
| (b3) CLR calibrated on simulated candidates (the package default for CR: the CLR branch union) | the true-branch inclusion and the crossing implications are deterministic; their nominal probabilities are conjectural | **Algorithm** (Conjecture P3.5(b3)); simulation coverage cannot promote it |

Further limits:
* **Regime.** M1/M2 with growing districts (R3) is a conjecture in P3 §0.5 independently of candidate collection.
* **CLR Condition P is not verified by a finite vertex set.** Their Condition P needs a fixed compact index representation along the sequence, gradients uniformly bounded away from zero, a well-conditioned covariance and a covariance-estimation rate. Nonzero gradients for one design are not enough.
* **Zero-variance candidates.** A candidate endpoint with $\hat\sigma=0$ is a different matter from an exactly known row intercept. Treat it as exact only if its gradient is exactly zero; otherwise check it separately.
* In the package, the default CR interval (CLR branch union) and, with `ci = "clr"`, the "already implied" intervals and the near-tie margin are therefore algorithmic. The gated Monte Carlo coverage (0.990) is a simulation finding.
* Version (a) needs no Slater-type condition. Statements of P3.2/P3.4 type do, and rows that bind at the truth (e.g. the CR edge) can break it.

(a) **Added rows.** For a finite set $\mathcal R$ of linear-in-$\beta$ rows, P3.5(a), P3.5(b) and P3.6(a) hold for the augmented system, version by version as in the table above, for $\tilde B$; for $B$ they hold with P3-S1a, the weight condition and the residual allowance carried over. The new rows join the Gaussian maximum that sets the critical value. A row whose $s_j$ does not depend on $\beta$ has $\sigma_j=0$: it is imposed exactly and left out of the maximum.
* **Actual targets that are contrasts.** For an actual gap $B_1-B_0$ the residual is $\sum_i(w_{i1}e_{i1}-w_{i0}e_{i0})$, and each unit term has range at most $w_{i1}+w_{i0}$, so $\tfrac12\{\sum_i(w_{i1}+w_{i0})^2\}^{1/2}$ is a valid standard-deviation envelope under conditional unit independence. The single-row $S_1$ need not bound it. The sum of the per-cell envelopes (which the package uses) is also valid, by Minkowski's inequality, and more conservative. *(Corrected after review, Batch B, Nov 2026.)*

(b) **Parallel disjunction.** Inflate each piece separately, then take the maximum:
$$\max_l\big(s_l(\hat\beta)+\hat k\hat\sigma_l\big)+a^\top k\ge0,$$
with all pieces included in the critical-value maximum. **For version (a) (conservative) only**, this covers at the nominal level. The row is not linear in $\beta$, so Lemma P3.1 fails (an endpoint becomes a max of minima) and CLR does not apply. For CLR, use (c): in W1/2×2 each branch is a set of linear-in-$\beta$ rows.
* W1 CR and $2\times2$ CR are of this form, by §1.2(d):
  * lower row: $h_D^\top k+\max(c_D-ER^*,\ c_D-\gamma_rER^*)\ge0$;
  * upper row: $-h_D^\top k+\max(ER^*-c_D,\ \gamma_rER^*-c_D)\ge0$.
* No branch is selected, so the sign of $\widehat{ER}$ never enters.
* The proof uses **one** simultaneous event: every population row intercept is at most its inflated estimate, and the fitted target constant is within its two-sided allowance. It needs no strict feasibility. At a tie, inflating all pieces and then maximising still dominates the true piece; taking the gradient of a winning estimated piece would be a different procedure. *(Added after review, Batch B.)*

(c) **General union.** $\mathrm{CI}_\cup:=\mathrm{CI}(C_+)\cup\mathrm{CI}(C_-)$, with *each* branch at the full level $\alpha$, satisfies
$$\liminf_nP(\theta_0\in\mathrm{CI}_\cup)\ \ge\ 1-\alpha.$$
* No Bonferroni split across branches is needed.
* The hull of $\mathrm{CI}_\cup$ is also valid; it is wider or equal.
* The emptiness test "reject if both branch sets are empty" has asymptotic level $\le\alpha_U+\alpha_L$.
* This is **point coverage** of the true target. It gives no simultaneous coverage of both branch images, and it does not license choosing a branch from $\widehat{ER}$. If designs form a triangular array and the true branch switches with $n$, branch validity is needed along each subsequence on which that branch is true (there are only two labels). *(Added after review, Batch B.)*

(d) **Near-empty sets.**
* Under the model, $P(\mathrm{CI}=\emptyset)\to\le\alpha_U+\alpha_L$, including when the identified set is a point or nearly one.
* Coverage holds unconditionally, so reporting "empty" is part of the allowed error, not an extra one.
* Conditioning on passing the emptiness test can only *raise* coverage: covering implies non-empty, so $P(\text{cover}\mid\ne\emptyset)=P(\text{cover})/P(\ne\emptyset)\ge P(\text{cover})$.
* Conditioning on other data events (a near-rejection, a CI that came out narrow) is not covered by the theorem.
* If the restrictions are false, no coverage statement holds. In particular, a non-empty, narrow CI is no evidence that they are true (G3). The EI referee's experiment makes this concrete: with correct linear means and a false CR, the conservative CR interval at 600 precincts covered 0 of 400 times under a large departure and was never empty (§3.4, V12).
* All statements are **pointwise** in the DGP; no uniform guarantee over sequences of nearly infeasible laws follows.
* Conditioning on non-emptiness offers no protection when the researcher selects a belief because it gives an attractive result. *(Added after review, Batch B.)*

(e) **Plug-in restrictions are not valid.** A rule that estimates the CR bracket $\operatorname{hull}\{\gamma_r\widehat{ER},\widehat{ER}\}$ and imposes it without inflating it can under-cover when the truth is on the bracket's edge, and the under-coverage does not vanish as $n$ grows. For a fixed repeating three-point design, $\limsup_nP(\tilde D_0\in\mathrm{CI}_{\rm naive})\le\Phi(\tfrac7{12}z_{.995})=0.933525<0.95$ (proof below). *(Corrected after review, Batch B: earlier supported by a Monte Carlo and an asymptotic sketch only.)*

### 3.3 Proof

*(a)*
1. Under the model, $k_0$ satisfies the base rows and the rows in $\mathcal R$ at $\beta_0$. [R3(e); hypothesis]
2. Each added $s_j(\hat\beta)$ is a linear map of $\hat\beta$, so it is jointly asymptotically normal with the base terms. Let $\mathcal E$ be P3.5(a)'s event, with the maximum taken over all rows. Then $P(\mathcal E)\to\ge1-\alpha_U-\alpha_L$. [OLS; lin; P3.5(a) step 1]
3. On $\mathcal E$, $k_0\in\hat\Theta^+$. P3.5(a) steps 3–7 and P3.6(a) then go through verbatim, since they use only $k_0\in\hat\Theta^+$ and the residual term. [1–2; P3.5, P3.6]
4. For (b) of P3.5, the endpoint is still a minimum of finitely many linear functions of $\beta$, because Lemma P3.1 needs only rows that are linear in $\beta$ with fixed $a_j$. CLR applies to the **augmented** system if its conditions hold there: the added rows change the dual feasible set, its vertices, the gradients, the selection and the quantile, so the conditions must be checked anew (status table). It inherits P3's open item on collecting candidates. [P3.1; CLR] ∎(a)

*(b)*
5. Let $l^*$ be a piece that holds at the truth: $s_{l^*}(\beta_0)+a^\top k_0\ge0$. [hypothesis; Def disjunction]
6. On $\mathcal E$, which now includes every piece, $s_{l^*}(\beta_0)\le s_{l^*}(\hat\beta)+\hat k\hat\sigma_{l^*}$. So $\max_l(s_l(\hat\beta)+\hat k\hat\sigma_l)+a^\top k_0\ge0$. [Def $\mathcal E$; 5]
7. For given data this is one linear row in $k$, and $k_0$ satisfies it on $\mathcal E$. Proceed as in step 3. [6; (a)]
8. For W1/2×2 CR: §1.2(d) gives $D^{r*}\ge\min(ER^*,\gamma_rER^*)\iff\max(D^{r*}-ER^*,\ D^{r*}-\gamma_rER^*)\ge0$, and similarly for the upper side. Each piece is linear in $\beta$ apart from the common $h_D^\top k$. [§1.3 steps 11–13; alg] ∎(b)

*(c)*
9. Fix a branch $j^*\in\{+,-\}$ with $k_0\in C_{j^*}$. It is non-random, because $k_0$ is a population quantity. [hypothesis]
10. The model "rung $\cap\,C_{j^*}$" holds, so by (a), $\liminf P(\theta_0\in\mathrm{CI}(C_{j^*}))\ge1-\alpha$. [9; (a)]
11. $\mathrm{CI}(C_{j^*})\subseteq\mathrm{CI}_\cup\subseteq\operatorname{hull}\mathrm{CI}_\cup$. [Def]
12. The other branch's CI is never used for coverage, so its level is irrelevant. That is why no split is needed. [10–11]
13. Both sets empty implies the $j^*$ set is empty, which has probability $\to\le\alpha_U+\alpha_L$. [P3.6(a) for version (a); the crossing argument of (d) step 15 for version (b); both applied to branch $j^*$] ∎(c)

*(d)*
14. For P3.5(a): $\mathrm{CI}=\emptyset$ requires $\hat\Theta^+=\emptyset$, which requires $\mathcal E^c$. [P3.5(a) step 2]
15. For P3.5(b), take $\theta_0=\tilde B$ (or $\tilde D$). A crossing $\hat L>\hat U$ requires $\hat U<\theta_0$ or $\hat L>\theta_0$. These have asymptotic probabilities $\le\alpha_U$ and $\le\alpha_L$ respectively, by CLR's one-sided validity, whatever the width of the identified set. For the target $B$, the E-allowance only widens the interval, so a crossing of the $B$-interval implies a crossing of the $\tilde B$-interval, and the same bound holds. [CLR; union; P3 allowance]
16. In both cases the empty event is contained in the non-coverage event whose probability is bounded. Neither step uses the width of the identified set. [14–15]
17. $\{\theta_0\in\mathrm{CI}\}\subseteq\{\mathrm{CI}\ne\emptyset\}$, so $P(\text{cover}\mid\ne\emptyset)=P(\text{cover})/P(\ne\emptyset)\ge P(\text{cover})$. If the restrictions are false, step 1 fails. [Def cond. prob.; 1] ∎(d)

*(e)* **Fixed-design proof** (the partial-identification referee's argument, adopted with credit; it replaces the earlier Monte Carlo sketch).
18. Repeat each of $X=1/4,1/2,3/4$ equally often, with equal sizes. Take $m_0=.4$ and $m_1=.55$ constant, and draw the group rates independently across units and groups from bounded laws with these means and positive variance (e.g. beta draws of fixed concentration). Every base probability row has strict slack, OLS has a full-rank fixed design and a positive definite asymptotic covariance, and HC3 is consistent. The positive CR branch has $\tilde D_0=ER=.15$ on its upper edge. [Def; OLS]
19. In the basis $((1-x)^2,2x(1-x),x^2)$, for **any** coefficient vector $q$: $ER(q)=q_{11}-q_{00}$ and $c_D(q)=\tfrac7{12}(q_{11}-q_{00})$, because $E[X^2]/E[X]=7/12$ and the middle coefficient cancels by symmetry. So $\operatorname{se}(c_D)/\operatorname{se}(\widehat{ER})=7/12$ exactly. [alg]
20. With probability tending to one, the plug-in sign is positive and the base rows allow the $k$ that attains the plug-in upper edge (strict slack and continuity). The naive upper endpoint is then $\widehat{ER}+k_{\rm crit}\operatorname{se}(c_D)$. [19; continuity]
21. Its critical value is calibrated on eight one-sided base rows and two target tails. By the union bound, the population Gaussian critical value satisfies $k_{\rm crit}\le z_{.995}$ at level .95, whatever the correlations. [union]
22. Non-coverage is $\{\widehat{ER}-ER<-k_{\rm crit}\operatorname{se}(c_D)\}$, so $\limsup_nP(\tilde D_0\in\mathrm{CI}_{\rm naive})\le\Phi(\tfrac7{12}z_{.995})=0.933525<0.95$. [20–21; OLS]
23. With 4,000 critical-value draws (as in `verify_R10.py`), each draw exceeds $z_{.996}$ with conditional probability at most .04 by the same union bound; an empirical .95 quantile above that threshold needs at least 200 exceedances. So the limiting coverage is at most $\Phi(\tfrac7{12}z_{.996})+P\{\mathrm{Bin}(4000,.04)\ge200\}=0.940092$. ∎ [union; binomial tail]

*Remarks.* (i) The far-edge design of §3.4 also satisfies the condition in its design limit: with $X=.1+.8U$, $U\sim$Beta$(2,2)$, symmetry gives $c_D(q)=(.5+2\operatorname{Var}X)(q_{11}-q_{00})$ with $\operatorname{Var}X=.032$, a ratio of .564 and a limiting bound of about .927 with population critical values (referee's derivation; the ratio is recomputed in V9). (ii) The Monte Carlo of §3.4 (coverage 0.917, Monte Carlo SE about 0.011, and the reviewer's 0.914, 0.909, 0.903 at $n=150$, 600, 2000) is now an illustration, not the argument. Check: V9.

**Uniformity.** The statements are pointwise in the DGP, as in P3. Uniform validity near point identification (a Stoye-type adaptation) remains P3's open item. (d) only says that emptiness is controlled pointwise, including at point-identified DGPs.

### 3.4 Monte Carlo (`verify_R10.py`, output `verify_R10_output.txt`)

Setup:
* 2×2, $n=150$ units, 600 replications.
* Unit rates are Beta draws with linear conditional means, so LCE and CR hold.
* HC3 covariance, conservative CI (P3.5(a)) at 95%, target $\tilde D$.

| Design | JKST | CR slab (b) | CR union (c): coverage / occupied length (hull span; unmerged sum; runs with overlap) | CR naive (plug-in bracket) |
|---|---|---|---|---|
| Truth on the CR edge ($\delta^*_W=0$), $ER^*=0.15$ | 1.000 / 1.106 | **0.992** / 0.359 | 0.992 / 0.357 (0.357; 0.386; 134/600) | **0.917** / 0.256 |
| Truth on the edge, $ER^*=0.005$ (near-point CR set, branch ambiguous) | 1.000 / 1.223 | 1.000 / 0.242 | 1.000 / 0.242 (0.242; 0.430; 591/600) | 1.000 / 0.147 |
| Interior truth | 1.000 / 1.232 | 1.000 / 0.480 | 1.000 / 0.480 (0.480; 0.480; 0/600) | 1.000 / 0.388 |

Entries are coverage / median length. No CI was empty in any design. *(Corrected after review, Batch B, Nov 2026: the union column earlier reported the sum of the two branch lengths, 0.386 / 0.430 / 0.480, which double-counts overlaps. It now reports the occupied (Lebesgue) length, with the hull span and the old unmerged sum beside it; the numbers reproduce the EI referee's correction exactly. The apparent length penalty of the branch union in the near-zero design disappears; the headline slab comparison, 0.359 against 1.106, is unaffected. The package merges overlapping pieces (`.rxc_union`) and was not affected. Check V11.)*

**What these runs do not test:** the CLR procedure (all four methods use the conservative version), coverage of the actual-rate target (residual allowance), the crossing emptiness rule, or a growing-dimension regime. The design is narrow: fixed composition and sizes, Beta rates, and cancelling slopes chosen to put the truth on the edge.

**A false CR (EI referee's experiment, reproduced as V12).** Same generator, with $m_0=.40-sx$ and $m_1=.50+s-sx$: the aggregate regression is $.40+.10x$ for every $s$, the between sign is positive throughout, and for $s>0$ both contextual slopes are negative, so CR is false. The design is fixed across $s$.

| Precincts | $s$ | true $\tilde D$ | JKST covered / 400 | CR slab covered / 400 | median JKST width | median CR width |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 0 | .100 | 400 | 399 | 1.262 | .306 |
| 150 | .08 | .170 | 400 | 336 | 1.256 | .304 |
| 150 | .25 | .318 | 400 | 2 | 1.254 | .303 |
| 600 | 0 | .100 | 400 | 398 | 1.074 | .194 |
| 600 | .08 | .169 | 400 | 110 | 1.069 | .196 |
| 600 | .25 | .317 | 400 | 0 | 1.066 | .195 |

No interval was empty. A false restriction gives narrow, nonempty and reproducibly wrong intervals; precision makes the wrong answer more reproducible, not more credible. This is the price R10(d) states, made visible. *(Added after review, Batch B.)*

### 3.5 Consequences for the package

1. **Every added row goes through the same CI machinery.** Sign, ordering and local rows join the critical-value maximum. Nothing new is needed beyond appending rows (R10(a)).
2. **CR (W1 default and 2×2) uses the slab form (b) for the conservative CI** (one row per side, each piece inflated before the maximum) **and the branch union (c) for CLR**, the default CI. That default is algorithmic (Conjecture P3.5(b3)); the conservative version and the projection interval carry the guarantee. In the edge design of §3.4 the **conservative slab** has median width 0.359 (coverage 0.992) against 1.106 for JKST alone; the default **CLR branch union**, in the package's gated run, had median width 0.230 and coverage 0.990 (§5). It is a large gain, and it rests on a belief whose refutability depends on the base system: in the model-free continuous table problem (2c M3, no linear rung) CR is always feasible; under LCE it can be rejected even when LCE alone fits (Remark R3.2); with integer counts it can be rejected outright (2c Remark M2.1). A false CR gives narrow wrong intervals (§3.4, V12). *(Corrected after review, Batch B, Nov 2026: earlier "it rests on an unrefutable belief (2c M3)", which imported M3 into the LCE base; and the quoted numbers did not name their method.)*
3. **Other CR pairings use the union of the branch CIs, each at the full level (c).** The package reports the union, which can be two intervals, and offers the hull as an option.
4. **Never impose an estimated bracket as if it were known.** The EGGHM bracket $\operatorname{hull}\{\gamma\,\widehat{ER},\widehat{ER}\}$ intersected with a JKST CI is not a confidence interval (e). Direct argument: in the design of (e), the upper endpoint of that intersection cannot exceed $\widehat{ER}$ with probability tending to one, so its limiting coverage of $\tilde D_0=ER$ is at most one half. *(Corrected after review, Batch B, Nov 2026: earlier "the MC naive interval contains bracket ∩ JKST CI, so the under-coverage transfers a fortiori". That containment is unproved and false in general: dropping an empty pre-expansion intersection and expanding endpoints do not commute. With raw base interval $[0,.1]$, bracket $[.15,.2]$ and target allowance .1, the raw intersection is empty although the expanded base meets the bracket. Check V10.)* When the regression and the target use different weights, the two sample ERs must be distinguished. (Our $\widehat{ER}^*$ is $\tau$-based. It equals the $T$-based EGGHM $\widehat{ER}$ under $N$-weighted OLS, because $X_r$ and $1$ lie in the span of the quadratic monomials. P3 uses unweighted OLS, so the two differ slightly in samples.)
5. **Emptiness is reported as a test result, not as an error.** Output includes:
   * the rejection flag;
   * the combination that was rejected (R3(e));
   * a warning that a narrow CI after a near-rejection is not reassurance (d).

---

## 4. Independent review log

A fresh-session adversarial review of §§1–3 and the three scripts. Verdicts: R3, R2 and R10 all **accept with fixes**. Every fix below has been applied.

| # | Item | Class | Finding | Fix |
|---|---|---|---|---|
| 1 | R10(d) | ERROR | It claimed coverage conditional on non-emptiness can fall below $1-\alpha$. In fact cover ⊆ non-empty, so conditioning raises coverage | (d) and step 17 rewritten |
| 2 | R10(b) | GAP | The slab row is not linear in $\beta$, so CLR is not covered | Slab restricted to version (a); CLR goes through (c) |
| 3 | R10 §3.1 | GAP | Two emptiness rules for version (b) would double the level | One rule per version |
| 4 | R10 step 15 | GAP | Wrong for $\theta_0=B$ (the E-allowance) | Argument run with $\tilde B$ |
| 5 | R10(a) | GAP | Rows with $\sigma_j=0$; conditions on the augmented system; Slater | Stated |
| 6 | R10(e) | GAP | Shown by MC only | Asymptotic sketch added |
| 7 | R10 minors | MINOR | Step-13 citation; "wider or equal"; a-fortiori remark; $\tau$- vs $T$-based $ER$ | Fixed |
| 8 | R3(c) | GAP | "Sharp iff intersect" is false when one branch is empty | Case split |
| 9 | R3 setting | GAP | Meaning of $\tilde B$ and of $ER^*$ under (C2) | Note added |
| 10 | R3 §1.6.2 | GAP | Contradicted R10(b): CR rows need inflation | Rewritten |
| 11 | R3 §1.4 | MINOR | Non-polyhedral images may be non-closed | Stated |
| 12 | R3 minors | MINOR | S1 tests (b), not (a); certificate is an LP only for $\Theta_O$; "weakly outer"; $\gamma$ used for two things (target normal now $h$); S5 gap sizes reach ≈ 31% on other seeds | Fixed |
| 14 | R2(c) | GAP | Proved for one row only | Separation step 11′ |
| 15 | R2 step 9 | GAP | Which optimum is meant | Local ⇒ global argument |
| 16–17 | R2 minors | MINOR | Domain boundary in [pert]; citations; nonemptiness of $F$ is an assumption | Fixed |
| 18 | R2 checks | GAP | T1 never tests a non-singleton face; T2 nearly tautological | T1b added; T2 relabelled |

Numerical re-checks by the reviewer: the identities in R3 step 11 hold for every $r$ to $1.5\times10^{-16}$; the S5 holes are genuine; the T4 28% "cuts set only" rate is not a tolerance artefact (shallowest cut $3.8\times10^{-6}$); `verify_R10.py` reproduces exactly.

**External review (Batch B, Nov 2026).** This internal review covered §§1–3 and three scripts. It missed the face-LP intercept (R2 step 4), the set-versus-target meaning of redundancy (R2(e)), the LCE example that excludes CR (Remark R3.2), the double-counted union lengths (§3.4), and the scope of the "unrefutable" sentence; and it left R10(e) to simulation. All are corrected above; the full ledger is `reviews/B_PI_RESPONSE_PLAN.md` §2. Process lesson: a pre-screen must also cover the scope document a proof relies on (here 2c), not only the new proofs.


---

## 5. Implementation status (build step 1, branch `claude/msign-build`)

**Implemented** (`package/eiPartialID`):

| Spec item | Where | Notes |
|---|---|---|
| Sign rows on $\delta_B$ (W1, pairwise) and $\delta_W$ (W1, W2 $=\Delta_r$), §1.6.1, R10(a) | `R/rxc_sign.R` (`.rxc_mlin`, `.rxc_assoc_W`, `.rxc_assoc_form`); rows added in `.rxc_extra_rows` (`R/rxc_core.R`) | `restrictions$assoc = list(list(type = "dB"/"dW", row, col, sign, version, other))`. One design-weighted ($N_i$) row per restriction, pooled over groups; unit $i$ uses its own group's $H$ |
| CR, W1 / 2×2 slab, §1.6.2, R10(b) | `.rxc_cr_parts` (slab = two disjunctive rows), `.rxc_fix_disj`; `.rxc_ci_conservative` (`R/rxc_infer.R`) | Point bounds: slab at $\hat\beta$. Conservative CI: every piece in the critical-value maximum, each inflated by its own $\hat k\hat\sigma_l$, then the max. `.rxc_ci_clr` refuses disjunctive rows |
| CR, W2 / pair, and CLR for every CR, §1.6.3, R10(c) | `.rxc_systems(cons, "branch")`, `.rxc_estimate` | Union of branch CIs, each at the full level; branches with empty / crossed / out-of-range intervals drop out; all empty = rejected. Several CRs give $2^m$ systems |
| Two intervals per target, hull option | `.rxc_union`, `.rxc_union_summary`; table columns `n_pieces`, `*_gap_lower/upper`; `fit$pieces`; `eiBounds(union = "pieces"/"hull")` | Printed as `[a, b] U [c, d]` |
| R2 redundancy diagnostic, §2.5 | `eiRedundancy()` (exported; stored in `fit$sign` by default when `assoc`/`cr` are given) | Face test via one feasibility LP per endpoint (`.rxc_meet`, `.rxc_face_moves`); duals never used. Verdicts implied / cuts set only / tightens, added first and added last, with width gains. If the base is itself a union (another non-W1 CR), the verdict falls back to a brute-force width comparison and says so (`method`) |
| "Already implied" table | `.rxc_implied_table` | Every $\delta$ (W1, W2, pairwise) under the rung alone: point interval and CI, signs flagged |
| CR certificate, §1.6.3 | `.rxc_cr_table` | $\Theta\cap\{g_B=0,g_W=0\}\neq\emptyset$ (feasibility LP), branch intervals, number of pieces |
| Emptiness, §1.6.4, R10(d) | `.rxc_empty_test` (conservative; reject iff every system's inflated set is empty; also reports whether the rung alone is rejected); notes in `eiBounds` | Messages name the combination "rung (M) + restriction labels"; 2-row note that $\delta_B$ / CR are never refuted on their own; "narrow interval after a near-rejection is not reassurance" |
| CLR no longer falls back to the conservative CI | `.rxc_collect(relax = TRUE)`, `.rxc_ci_clr` | When no draw is feasible, candidates come from rows relaxed by the smallest uniform slack. The dual feasible set does not depend on the row constants, so these vertices are dual-feasible for the original rows and bound the endpoint by weak duality (P3.5(b2) step 1). Coverage of the CLR critical value over the collected set is Conjecture P3.5(b3). The CLR rule (ends cross ⇒ empty) is kept, so the two emptiness rules are never mixed |
| Descriptive $\widehat{ER}$, $\widehat{NM}=\gamma\widehat{ER}$, $\gamma$, bracket, §3.5.4 | `.rxc_descriptive`, printed under "DESCRIPTIVE ONLY" | $\tau$-based; equals the EGGHM $T$-based value under `weights = "n"` (tested to 1e-10). Never imposed |
| Integration | `ladder(sign_rungs = ...)`, `print`/`summary`, `sensitivity()`, `breakdown()` (pass `union`) | Point-identified rungs (L5, L6) ignore sign/CR rows, with a note |
| Docs | roxygen (`eiBounds`, `eiRedundancy`, `ladder`), vignette `eiPartialID_rxc.Rmd` section "Signs on associations and contextual reinforcement", `PLAIN_LANGUAGE_GUIDE.md` §9 | Applied framing as in §3.5 |

**Tests** (`tests/testthat/test-msign.R`; reference numbers from `docs/extension/make_msign_reference.py`, which calls `verify_R3.py`'s `assoc_3x2` / `assoc_alt` and writes `msign_example_3x2.csv`, `msign_gap_3x2.csv`, `msign_reference_3x2.csv`):

* all association forms, ER, NM, γ, D and $\tilde B$ match Python to 1e-9 on a noise-free 3×2 example and on the gap example;
* the W1 identities hold for random $k$ to 1e-12;
* 3×2 CR point bounds (W1, W2, pair; three targets) equal the Python branch unions; W1 slab = union of branches; the gap example (near-1-d support plus one ordering row, found with the S5 search) gives two pieces for W2 and pair, certificate false, and the hull option;
* noise-free gap example: both CI versions return the two pieces;
* 2×2 closed forms of R2(d) for all four sign rows and CR (voter_vax and two simulated designs, both regression weightings); 2×2 CR slab = the single non-empty branch;
* MonotoneEI `voter_vax` (bundled as `voter_vax_MonotoneEI.csv`, MIT, attribution in `voter_vax_MonotoneEI.README`): CR bounds for the gap $[-0.4699, -0.0502]$, ER $=-0.470$, NM $=-0.050$, N-weighted;
* redundancy verdicts = brute force for 24 sign rows × 3 targets (both "tightens" and "cuts the set, not this target" occur); the middle verdict never occurs in 2×2; added-last gains = width differences;
* joint refutation and its messages; CR CIs contain the estimate and are narrower than the rung's; an empty CLR branch is dropped via crossing ends, not by fallback;
* `sensitivity`, `breakdown`, `ladder(sign_rungs)` run with sign/CR rows; input checks;
* gated coverage MC (`EIPARTIALID_SLOW=1`): truth on the CR edge, $n=150$, 200 reps: conservative slab **1.000** (median width 0.346), CLR branch union **0.990** (median width 0.230).

**Run:** R 4.3 (Ubuntu apt) in the build session: full `testthat` suite 0 failures (skips: `ei.Datasets` not installed, slow MC unless enabled); the slow MC was run once with the result above; `R CMD check --no-manual --no-vignettes`: no errors (pre-existing notes: unused imports `MASS`/`eco`, unavailable Suggests; a locale warning from the container). The vignette was rendered with `rmarkdown::render`.

**Open / not done:**
* Uniformity near point identification remains P3's open item; the CLR candidate collection is unchanged apart from the relaxation above.
* With several non-W1 CRs, systems grow as $2^m$; no pruning.
* The redundancy verdicts are point-level; the R2 §2.5.4 "width change with its CI" near ties is not yet reported (the "already implied" table does give CIs for the associations).
* Out of scope as instructed: N1, N2, M-LOC (R1), R×C rank (R5), DiD.
* *Corrected after review (Batch B, Nov 2026):* this list earlier also named "slack $\eta$ on $\delta$". That was stale: a per-restriction slack in probability units is implemented (`sensitivity(dial = "slack")`, `slack` in each restriction, `tests/testthat/test-msign-dial.R`). The stale line misled the EI referee.
