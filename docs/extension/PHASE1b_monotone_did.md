# Phase 1 (second cycle): monotone ecological inference and difference-in-differences

*A proposal for review, covering structure and assumptions only. Nothing here is a Phase-2 derivation. Status labels are used throughout:*
* ***proved*** *only for one-line facts whose proof is given here and that are also checked numerically;*
* ***checked*** *for identities verified numerically in `verify_phase1_mei_did.py` (output: `verify_phase1_mei_did_output.txt`, 10/10 PASS), with formal proofs due in Phase 2;*
* ***conjecture*** *or* ***heuristic*** *otherwise.*

*Cited as P1, 2a, 2b, P3 and P4: the earlier phase documents. Cited as MENU: `ASSUMPTION_MENU.md`.*

---

## 0. Sources read, and what I could not read

* **Elzayn, Goldin, Guage, Ho & Morton, "Monotone Ecological Inference"** (EGGHM below).
  * *Political Analysis* (2026), First View, pp. 1–20, doi:10.1017/pan.2026.10048, CC BY 4.0.
  * I read the full arXiv v1 (2504.14752, March 2025, 86 pages with appendices: proofs, local results without the same-prevalence assumption, data and matching, covariance estimation). I read the published HTML only in part, through a summarising fetch.
  * **Numbering differs between versions.** Below I cite arXiv v1 numbers. The published version appears to renumber; for example, the contextual-reinforcement result is reported there as Prop. 6 / Cor. 5, and the method of bounds as Prop. 4. *All published-version numbers must be checked against the PDF before citing.*
  * Not read:
    * the published Supplementary Material, which the abstract page says covers more than two groups and continuous outcomes;
    * the Code Ocean capsule (doi:10.24433/CO.6554779.v1).
* **R package `MonotoneEI`** (github.com/reglab/MonotoneEI, v0.1.0, MIT). I read `monotone_EI()`, `MOB()`, `ER()` and `NM()`. Two observations matter for later cross-checks:
  1. Its "monotone" interval is $[\min(NM,ER),\max(NM,ER)]\cap MOB$, i.e. the contextual-reinforcement interval.
  2. Its confidence limits take the bootstrap distribution of whichever estimator binds at the plug-in, a selection that ignores near-ties. Also, `MOB()` does not appear to use the weights `wt`, while `ER()` and `NM()` do. Cross-checks will therefore use equal weights.

---

# Part A. Monotone ecological inference

## A1. What the paper assumes and proves (arXiv v1 numbering)

**Setting.** There are individuals $(X,Y,N)$: binary group $X$, outcome $Y\in[\underline Y,\overline Y]$ (bounded, not necessarily binary), and neighbourhood $N$ from a finite set. Per neighbourhood one observes $X_n=E[X\mid N=n]$, $Y_n=E[Y\mid N=n]$ and $p_n=\Pr(N=n)$. This is a **finite population with no superpopulation model and no sampling**; sampling enters only in the application.

**Estimands.** $Y^x=E[Y\mid X=x]$ and $D=Y^1-Y^0$. Neighbourhood-level versions $Y^x_n$ and $D_n$ appear in §5 of the paper.

**The two conditional associations** (footnote 6 gives the 2-group forms):
$$\delta_B=E[\operatorname{Cov}(Y,X\mid N)]=\sum_np_nX_n(1-X_n)(Y^1_n-Y^0_n),$$
$$\delta_W=E[\operatorname{Cov}(Y,X_N\mid X)]=E[X]\operatorname{Cov}_{w^1}(Y^1_N,X_N)+(1-E[X])\operatorname{Cov}_{w^0}(Y^0_N,X_N),$$
with $w^1_n\propto p_nX_n$ and $w^0_n\propto p_n(1-X_n)$.
* $\delta_B$, the *between-group* association, is the within-neighbourhood gap between the groups.
* $\delta_W$, the *within-group* association, is how each group's mean moves with neighbourhood composition, i.e. a contextual effect.

**Identities** (all checked: V1).
* Prop. 1: $D(1-\gamma)\operatorname{Var}(X)=\delta_B+\delta_W$, where $\gamma=\operatorname{Var}(X_N)/\operatorname{Var}(X)$.
* Lemma 1: $D_{NM}=\gamma D_{ER}$.
* Prop. 4: $D_{ER}-D=\delta_W/\operatorname{Var}(X_N)$.
* Prop. 6: $D_{NM}-D=-\delta_B/\operatorname{Var}(X)$.
* Props 9–10: the group-mean versions,
  * $Y^1_{ER}-Y^1=\delta_W(1-E X)/\operatorname{Var}(X_N)$,
  * $Y^0_{ER}-Y^0=-\delta_W E X/\operatorname{Var}(X_N)$,
  * $Y^1_{NM}-Y^1=-\delta_B/E X$,
  * $Y^0_{NM}-Y^0=\delta_B/(1-E X)$.

Here:
* ER is Goodman's regression *with intercept*, i.e. the population linear projection, so linearity is not assumed;
* NM is Freedman et al.'s neighbourhood model, $Y^1_{NM}=E[X_NY_N]/E[X_N]$.

**Bounds.**
* Props 2–3 give the method of bounds (MOB), which is Duncan–Davis aggregated. These bounds are sharp in their finite-population sense.
* **Sign restrictions** (Props 5, 7):
  * $\delta_W\ge0\Rightarrow D\in[D^-_{MOB},\min\{D^+_{MOB},D_{ER}\}]$;
  * $\delta_W\le0\Rightarrow D\in[\max\{D^-_{MOB},D_{ER}\},D^+_{MOB}]$;
  * $\delta_B\ge0\Rightarrow D\in[D_{NM},D^+_{MOB}]$;
  * $\delta_B\le0\Rightarrow D\in[D^-_{MOB},D_{NM}]$.
* **Contextual reinforcement (CR)**, $\delta_B\delta_W\ge0$ with the common sign unknown (Def. 1, Prop. 8, Cor. 3):
  * $D$ lies between $D_{NM}$ and $D_{ER}$, intersected with MOB;
  * $\operatorname{sign}D=\operatorname{sign}(D_{ER}-D_{NM})$.
* Theorem 1 and Tables 1–3 cover all sign combinations, for $D$, $Y^1$ and $Y^0$. Every bound is stated as sharp.

**Neighbourhood-level results** (§5).
* Prop. 11: the local sign $\delta_{B,n}=\operatorname{Var}(X\mid n)D_n$ cuts the neighbourhood MOB at 0 (Manski 2018's bounded variation).
* Prop. 12:
  * Write $\mu(x)=E[Y\mid X_N=x]=x\mu_1(x)+(1-x)\mu_0(x)$.
  * The identity $\mu'(x)=D(x)+\delta_W(x)$ holds, with $\delta_W(x)=x\mu_1'(x)+(1-x)\mu_0'(x)$.
  * A sign on $\delta_W(x)$ therefore bounds $D(x)$ by the derivative $\mu'(x)$.
  * Assumption 1 says neighbourhoods with equal $X_n$ share group means. Props 14–16 drop it and target the average over neighbourhoods with that $X_n$.
* Prop. 13 is the local CR version.

**Several groups** (§4, end). Collapse to "group $g$ vs the rest" and apply the 2-group results with $\delta^g_W$ and $\delta^g_B$.

**Inference** (§6). County bootstrap (1000 replicates) for standard errors, then Imbens–Manski (2004) intervals. $\mu'$ is estimated by population-weighted local linear regression (Epanechnikov kernel).

**Application.** The partisan COVID-vaccination gap across 3,115 counties.
* MOB gives $[-77.4,52.8]$ points.
* CR gives $[-47.9,-5.5]$, i.e. $[D_{ER},D_{NM}]$, with 95% CI $[-51.8,-5.1]$.
* The signs of $\delta_B$ and $\delta_W$ are validated on a linked individual-level auxiliary sample.

## A2. How it relates to what we have

**(a) JKST linear contextual effects.** Under 2×2 LCE ($E[\beta_1\mid X]=b_0+b_1X$, $E[\beta_0\mid X]=w_0+w_1X$, $b_1=w_1+d_1$), both associations are **affine in the single unidentified parameter $w_1$** (checked: V3):
$$\delta_W=a\,(w_1+d_1)+b\,w_1,\qquad a=E_N[X]\operatorname{Var}_{w^1}(X),\ \ b=E_N[1-X]\operatorname{Var}_{w^0}(X),$$
$$\delta_B=E_N\big[X(1-X)(c_1-w_1+d_1X)\big].$$
Consequences:
* Each EGGHM sign restriction is **one half-space in $K$**, i.e. one more LP row in 2a Prop. 7.
* CR is a **union of two polyhedra**. It is not convex, but it is two LPs whose intervals are joined.
* $\delta_W$ is a *positively weighted combination of both groups' contextual slopes*. Signing it is weaker than signing each slope (MENU B6, the package's `sign=`).
* EGGHM's ER is the *linear* projection. JKST use the *quadratic* fit and get a set rather than a point because the quadratic separates $\delta_W$ into slope pieces. The two approaches use different moments, so neither nests the other.

**(b) Sharp and outer sets (2a).** EGGHM's information set is "data alone + a linear inequality on the unit-level rates". In R×C language:
* the feasible set is $\{\beta=(\beta_i)_i:\beta_i\in S(X_i,T_i)\ \forall i,\ G\beta\ge0\}$, the product of unit Duncan–Davis polytopes cut by aggregate half-spaces, because $\delta_B$ and $\delta_W$ are linear in $(\beta_i)$;
* the target is linear in $\beta$;
* so their bounds are an LP, and the closed forms coincide with that LP (checked: V2, 1093 non-empty and 107 empty cases).
* **Sharpness.** In the finite-population sense, the LP over this set is sharp by construction: every feasible $\beta$ is a table consistent with the margins. *Proved, one line.*
* **What is new for our taxonomy.** This is a **cross-unit, aggregate** restriction on *realised* unit rates. It is neither unit-level (it does not shrink any $S(X_i,T_i)$) nor mean-level (it does not restrict $E[\beta\mid X]$).
* **Superpopulation version.** If the sign is assumed for the conditional-mean analogue ($\delta^*$ computed from $m=E[\beta\mid V]$), the realised $\delta$ differs from $\delta^*$ by $O_p(p^{-1/2})$. The finite-population sign is then a slightly different assumption (see A4).
* **When the set is empty.** With $\delta_W\ge0$ the set is empty whenever $D_{ER}<D^-_{MOB}$, which happened in 107 of 1000 random populations in V2. EGGHM's tables then print an inverted interval. Emptiness *refutes* the sign assumption, as P3 Prop. 6 treats it for our rungs.

**(c) The group rungs of 2b.**
* EGGHM's neighbourhood is our *unit*; they have no districts.
* With districts, each sign can be imposed per district, on the pooled population, or **within districts**. The within-district version uses district-demeaned $X_N$, i.e. ER with district fixed effects, which is the natural partner of M1.
* The FE-ER bias formula is *to derive*. I expect it to be $\delta^{\text{within}}_W/\operatorname{Var}(X_N-E[X_N\mid D])$; *conjecture*.

**(d) Menu items B5–B8.**

| MENU item | EGGHM analogue | Relation |
|---|---|---|
| B5 group ordering, unit level ($\beta_{ic}^{r}\ge\beta_{ic}^{s}$ for each $i$) | local $\delta_{B,n}\ge0$ (Prop. 11; Manski 2018) | identical in 2 groups |
| B5 ordering, mean level (package `order=`: $M_H(x)_{rc}\ge M_H(x)_{sc}$ on $A$) | $\delta_B\ge0$ is its $X(1-X)$-weighted average | mean-level order ⇒ $\delta_B\ge0$ (*proved*: a weighted average of non-negative terms) |
| B6 sign of each row's contextual slope (package `sign=`) | $\delta_W\ge0$ is a positively weighted sum of the two rows' slopes | sign on both rows ⇒ $\delta_W\ge0$ (*proved* under LCE from A2(a); nonparametrically via V5) |
| B8 monotone nonparametric contextual effects | Prop. 12, local $\delta_W(x)\ge0$ | "each $m_r$ nondecreasing" ⇒ local sign ⇒ … but **local sign everywhere does not imply aggregate $\delta_W\ge0$** (V6 counterexample), so the local and aggregate rungs are non-nested |
| — (new) | contextual reinforcement | sign-free: needs no direction, only agreement |

**(e) Cross & Manski (2002), Cho & Manski (2008), Manski & Pepper (2000).**
* EGGHM's group mean is Cross–Manski's *short* regression $E[Y\mid X]$, and $Y^x_n$ is the *long* one (EGGHM fn. 12).
* MOB is the Cross–Manski "data alone" region aggregated, our P1 (D0).
* EGGHM's signs are *moment inequalities on the latent joint distribution*. They are neither Cross–Manski's instrument (constancy within cells of $V$) nor Manski–Pepper (2000) MIV.
* The closest Manski–Pepper analogue is "**each group's conditional mean $m_r(x)$ is monotone in the composition $x$**", i.e. MIV with composition as the monotone instrument. It implies $\delta_W\ge0$ (V5, Chebyshev's association inequality) and is what I propose as the nonparametric rung (A3, rung N2).
* For the layering, Cho–Manski's "decreasing credibility" ordering applies as in 2a §12. EGGHM's aggregate signs sit just inside DD.

## A3. Proposed objects and rungs

**R×C definitions** (column $c$ fixed; all linear in the unit rates):
* **Between, pairwise:** $\delta^c_B(r,s)=E_p[X_rX_s(\beta_{rc}-\beta_{sc})]$. This is antisymmetric. Then $NM_{rc}-B_{rc}=-\sum_s\delta^c_B(r,s)/\mu_r$ (checked: V4).
  * EGGHM's one-vs-rest $\delta^{r}_B=\sum_s\delta^c_B(r,s)$. So pairwise signs imply the one-vs-rest sign.
* **Within, two candidate generalisations** (your choice, Q2):
  * (W1) EGGHM's **one-vs-rest**: collapse rows to $r$ vs $\neg r$, then use their $\delta^r_W$ and simple ER of $T_c$ on $X_r$.
  * (W2) the **Goodman-native vector** $\Delta^c=\sum_s\mu_s\operatorname{Cov}_{w_s}(X,\beta_{sc})\in\mathbb R^R$, with $\mathbf 1^\top\Delta^c=0$. Goodman's multiple regression (no intercept) then satisfies $g_c-B_{\cdot c}=M^{-1}\Delta^c$ with $M=E_p[XX^\top]$ (checked: V4). A sign on a component, or on the bias $M^{-1}\Delta^c$ itself, is again linear.
* **Contextual reinforcement:** for chosen pairs, $\delta_B\cdot\delta_W\ge0$, evaluated as the union of the two same-sign LPs.

**Three ways to use them** (where each acts):

| Tool | Acts on | Needs linearity? | Computation | Identified set |
|---|---|---|---|---|
| **N1** "monotone EI on DD" (EGGHM generalised): DD ∩ chosen signs / CR | realised unit rates, aggregate (cross-unit) | no | one LP with $p(R-1)(C-1)$ unit variables (NZ: about 64×4) | sharp in the finite-population sense (proved, one line); joint across columns and rows, so weakly tighter than separate one-vs-rest closed forms (*proved*: a subset of constraints; strictness *conjecture*) |
| **K-rows**: the same signs under LCE, rows added to the LP in $K$ (any rung M0, M1, M1(δ), M2) | conditional means (superpopulation) | yes | existing LP plus 1–2 rows | outer (never called sharp) |
| **N2** monotone nonparametric contextual effects: each $m_r(x)=E[\beta_r\mid X=x]$ monotone in chosen coordinates (Manski–Pepper MIV on composition); optional local CR (EGGHM Prop. 13) | conditional means, pointwise in $x$ | no | LP on a grid of compositions: pointwise $m(x)\in S(x,\tau(x))$ (outer) plus monotonicity between grid neighbours; $\tau$ estimated | outer; a sharp version would need $E[S(X,T)\mid X=x]$ (2a Thm 1) |

**Ladder placement** (2 rows; inclusions among model classes, so the bounds nest by [nest]):
$$\underbrace{DD}_{L0}\ \supseteq\ \text{N1}(\delta_B\ \text{or}\ \delta_W\ \text{sign})\ \supseteq\ \text{N1(both signs)},\qquad DD\supseteq\text{N1(CR)}\supseteq\text{N1(both signs)},$$
$$\text{N1}(\delta_W\ge0)\ \supseteq\ \text{N2(each }m_r\uparrow)\ \supseteq\ \text{M0}\cap\{\text{both slopes}\ge0\},\qquad \text{N2}\cap\text{mean order}\ \subseteq\ \text{N1(both}\ge0).$$
* The first inclusion on the second line uses the superpopulation $\delta^*$ (V5).
* The last inclusion is *proved* (A2(d)).
* These inclusions are between model classes, so they nest the *sharp* sets. For the *outer* sets actually computed (vertex constraints for M0, grid constraints for N2) the nesting still has to be checked; Phase 2.
* **M0 (linearity) is non-nested with N1 and N2.** Intersections such as "M0 ∩ $\delta_W\ge0$" are valid when both hold and are reported as such.
* Plausibility ordering, my judgment for your review:
  * CR and $\delta_B$ signs are the most plausible additions to DD (no direction needed / loyalty-type);
  * then $\delta_W$ signs;
  * then N2 (monotone in every composition, row by row);
  * then the LCE rungs.

**Unit level vs mean level.** N1 restricts realised rates, which gives exact bounds for the *actual* rate. K-rows and N2 restrict conditional means, which gives bounds for the *model-implied* rate plus P3's allowance for the actual rate. EGGHM's local results (Props 11–13) are unit-level when Assumption 1 holds and mean-level otherwise (their Props 14–16).

## A4. Inference (proposal)

* **N1, finite-population reading.** If the sign is an assumption about the *realised* population, the bound is a deterministic function of the observed table and needs no CI. I propose reporting it with a **slack dial** $\delta_W\ge-\varepsilon$ and the **breakdown value** $\varepsilon^*$ of a conclusion. This reuses `breakdown()`.
* **Superpopulation reading** (EGGHM's own, via the bootstrap).
  * In 2×2 the endpoints are a min or max of finitely many smooth functionals: MOB±, ER, NM. That is exactly P3 Thm 5(b)'s finite-candidate CLR setting, so it is valid uniformly over near-ties.
  * EGGHM's bootstrap + Imbens–Manski selects the plug-in binding candidate, which can under-cover near ties. I will show this in a small Monte Carlo rather than assert it.
  * For the joint R×C LP, candidate collection by dual vertices is P3's device, with the same *heuristic* status as there.
* **K-rows.** Unchanged P3 machinery.
* **N2.** $\tau$ is fitted by a polynomial of degree $J$ (with $J=2$, N2 contains the LCE rungs). The grid constraints are then linear in the coefficients, so P3's CLR applies *conditional on $J$*. That is valid if $\tau$ is a degree-$J$ polynomial and a *heuristic* otherwise. A fully nonparametric version (series or kernel CLR) is possible but is Phase-3 work of its own; I would defer it.
* **CR.** It is a union of two LPs. Use a CI for each branch at level α and take the union (*valid by the union bound; proof in Phase 3*). Also report the test of $\operatorname{sign}(D_{ER}-D_{NM})$ (EGGHM Cor. 3).

## A5. Implementation and benchmark plan (after approval)

* **API** (backward compatible):
  * a new argument `monotone = list(between = ..., within = ..., reinforce = ...)` for `eiBounds()`, taking (row, column, sign) triples;
  * a new `model = "DD"` path that solves the unit-level LP (N1);
  * `model = "MIV"` for N2;
  * K-rows reuse the existing `linear=` machinery;
  * `ladder()` gains optional N1/N2 rows.
* **Tests:**
  * 2×2 closed forms against the LP (as V2);
  * equality with `MonotoneEI` on equal weights;
  * nesting;
  * emptiness reproduces $D_{ER}<D^-_{MOB}$.
* **Benchmark** (`eiBenchmark()`, NZ 2002–2020, 492 tables). Rungs:
  * DD;
  * N1(CR);
  * N1(between, own-party candidate);
  * N1(within +), and within − as a falsification;
  * M0 ∩ between;
  * N2;
  * plus the existing rungs.
* **What the benchmark can and cannot show.** Truth is known only for the electorate-wide table. The benchmark therefore checks coverage of the target, not whether the sign assumptions hold.

---

# Part B. Difference-in-differences with partial identification

## B1. Setting and targets

* Units $i$ (counties) are observed in periods $t\in\{0,1\}$ (pre, post), possibly with several pre-periods.
* Rows: $v$ (veteran), $n$ (non-veteran).
* Columns on a **common base** (eligible men): yes, no, abstain. So each period is **at least 2×3**, and $\dim\mathcal K=2$ per period (2a Lemma 1).
* Targets, linear unless noted:
  * $\text{DiD}_c=(B_{vc,1}-B_{vc,0})-(B_{nc,1}-B_{nc,0})$ for $c=$ yes, as a share of eligible men;
  * group changes $\Delta B_{rc}$;
  * the levels;
  * *if wanted*, support among voters, $B_{v,\text{yes}}/(B_{v,\text{yes}}+B_{v,\text{no}})$. This ratio is **linear-fractional** in $K$ (Charnes–Cooper or bisection; needs both columns' $K$).
* **Correction to MENU item 18.** A DiD for **one column** depends only on that column's $K_c$ (2a Lemma 2). With $R=2$ it therefore needs **one** restriction on $\Delta\kappa_c$, not $(C-1)R(R-1)/2$. The larger count is needed only to point-identify *every* column at once. *Proved: Lemma 2 makes $B_{rc}$ depend on $K$ only through $K_c$.*

## B2. The theory to develop (with status)

| # | Object | Status |
|---|---|---|
| D0 | **No cross-period link**: the identified set is the product; the DiD bound is the difference of the periods' intervals | model-free version *checked* (V7a). LCE version: the sharp set is the product of the periods' sharp sets. *Conjecture* with a short proof route: the Aumann expectation of the product $S(X,T_0)\times S(X,T_1)$ is the product of Aumann expectations (selections of a product are pairs of selections). Phase 2 |
| D1 | **First-difference LCE**: $\beta_{it}=\zeta_i+M_{H_t}(X_i)+e_{it}$, $\zeta_i$ arbitrary (MENU 21). Levels are *not* linear; $E[\Delta\beta\mid X]=M_{\Delta H}(X)$, $E[\Delta T\mid X]=X^\top\Delta Q X$ | **the owner's point-identification result holds here, under a weaker model than M1 across periods**: one restriction on $\Delta\kappa_c$ plus the quadratic fit of $\Delta T$ (*checked*: V7b, exact despite non-linear unit effects). Bounds on $\Delta K$: $\Delta\beta_i\in S(X_i,T_{i1})\oplus(-S(X_i,T_{i0}))$, a semi-infinite LP checked at the sample compositions. Phase 2 |
| D2 | **Single identifying restrictions** (table B3) | each is just-identifying, hence untestable except through emptiness; pairs over-identify, e.g. $\Delta d_1=0$ |
| D3 | **Bounded change** $\lvert\Delta\kappa_c\rvert\le\delta$, or $\lvert\Delta w_1\rvert\le\delta$ (Manski–Pepper 2018 bounded variation; Rambachan–Roth 2023) | DiD = point ± (sensitivity)·δ, intersected with D0. Calibrating δ: RR "relative magnitudes" uses *identified* pre-trends, but pre-period $\Delta\kappa$ is not identified. Scaling δ by identified pre-period changes (e.g. $\Delta d_1$ between pre-periods) is only a *heuristic* |
| D4 | **Monotone EI on changes** (no linearity at all): apply A1 to $Y=\Delta Y\in[-1,1]$ | $ER(\Delta T\sim X)-\text{DiD}=\delta^\Delta_W/\operatorname{Var}(X_N)$ (*checked*: V7c). Under CR of changes, **DiD lies between NM and ER on $\Delta T$**. The MOB for changes must use the two periods' polytopes jointly, not $[-1,1]$. Sign of $\delta^\Delta_W$ = sign of differential spillover; $\delta^\Delta_B$ = the DiD's own sign |
| D5 | **Intersection across periods** (M1 with periods as groups, 2b Remark 6.2) | already *proved* in 2b. Tightens the levels; refutable through $\Delta d_1\neq0$ or emptiness |
| D6 | **Parallel trends** for the causal claim, kept separate from EI. ATT of service = DiD if future veterans' counterfactual change equals non-veterans'. Bounded version: ATT ∈ DiD ⊕ [−M, M] | interplay with spillovers: a non-zero non-veteran change caused by returning veterans contaminates the control. D2's restriction (R1) rules out only the *differential* part. Placebo with two pre-periods: the placebo DiD is itself only partially identified, so it refutes PT only if its bound excludes 0 |
| D7 | **Composition change** (B4) | design choice needed from you |
| D8 | **Within-person table** (veteran × pre-vote × post-vote) | for the **DiD of means**, it adds nothing beyond the two two-way tables that share the veteran margin. *Proved, sketch*: any such pair extends to a 3-way table by coupling $Y_0$ and $Y_1$ independently given $V$. It matters only for *transition* targets such as "share of veterans who switched to yes" |

**Nesting** (model classes; bounds nest by [nest]):
$$\text{M1-periods}\subseteq\text{D1}+\{\text{R1 and R2}\},\qquad \text{per-period LCE (M2-periods)}\subseteq\text{D1},\qquad \text{D1}+\text{(any single restriction)}\subseteq\text{D1}+\text{D3}(\delta)\subseteq\text{D1}.$$
D4 needs no linearity, is non-nested with D1, and sits inside D0's model-free bound. *All inclusions follow from the definitions; written proofs in Phase 2.*

## B3. Single restrictions that identify the DiD (2 rows, one column)

Let $H_{rs}$ be the rate of $r$-men in a hypothetical county made up entirely of $s$-men (P1 §7.2); per period, $\kappa=(H_{vn}-H_{nv})/2$ is unidentified. Each restriction below pins $\Delta\kappa$:

| Restriction | Meaning (veterans) | Comment |
|---|---|---|
| (R1) $\Delta w_1=\Delta(H_{nv}-H_{nn})=0$ | non-veterans' change does not depend on the local veteran share: **no differential spillover onto non-veterans** | most natural; a uniform spillover is still allowed and biases the *causal* reading (D6) |
| (R2) $\Delta b_1=\Delta(H_{vv}-H_{vn})=0$ | veterans' change does not depend on how many fellow veterans live nearby (no extra network or veterans-organisation effect) | plausible if veterans' attitudes formed in the army, not at home |
| (R3) $\Delta\kappa=0$ | the change of veterans in non-veteran areas equals the change of non-veterans in veteran areas | hard to motivate on its own. *Not* M1 across periods: M1 shifts each row's baseline by $\zeta_r$, so $\Delta\kappa=(\zeta_v-\zeta_n)/2\neq0$ in general (checked: V8; corrected after review) |
| (R4) $\Delta\delta_W=a\Delta b_1+b\Delta w_1=0$ | the composition-weighted average contextual change is zero | = Goodman's linear regression on $\Delta T$ is unbiased (D4) |
| (R1)+(R2) | both slopes unchanged; this **is** M1 across periods ("contextual effects unchanged", 2b Remark 6.2) | over-identifying, testable: $\Delta d_1=0$ (equal $X^2$ coefficients) |

Signed versions turn these into one-sided cuts: $\Delta w_1\ge0$ means non-negative differential spillover. Bounded versions are D3.

## B4. Composition change (your input needed)

Pre-war "veterans" are *future* veterans. Many were too young to vote before the war, war deaths removed others, and migration changes both electorates. Options, in increasing ambition:
1. **Interval-valued pre-war composition.**
   * In cell-count coordinates, a unit's table with interval row totals and fixed column totals is still a polytope. So model-free and N1/D4-type bounds stay LPs, with linear-fractional targets.
   * LCE with an interval regressor is *not* an LP; it would be a heuristic outer bound.
2. **Explicit stayer / exit / entrant rows** on a common base, with shares from census age structure, mortality and migration. These shares are point or interval.
3. **Within-stayer target** $E[\Delta Y\mid V=1,\text{stayer}]-E[\Delta Y\mid V=0,\text{stayer}]$. By D8, this needs only the two two-way tables for stayers.

## B5. Inference

* Counties appear in both periods, so the two periods' regressions are stacked and **clustered by county**.
* Candidates from both periods enter one CLR intersection (P3 Thm 5).
* For D1, regress $\Delta T$ directly.
* **Two targets.** The *model-implied* DiD has no unit-level allowance. The *actual-rate* DiD's allowance uses $\operatorname{Var}(e_1-e_0)\le(\sqrt{m_1(1-m_1)}+\sqrt{m_0(1-m_0)})^2$ (Bhatia–Davis per period plus Cauchy–Schwarz). That can be up to 4× the single-period allowance. I propose the model-implied DiD as the headline (Q8).
* **Sensitivity.** Breakdown values in δ (D3) and in M (D6).

---

# Part C. Plan and questions

**Proposed order.**
1. Task 1 Phase 2 (derivations for N1, K-rows and N2), then Phases 3–4 and the NZ benchmark.
2. Task 2 Phase 2 (D0–D8) in parallel with collecting the data.
3. Task 2 Phases 3–4 after your answers on composition.

**Housekeeping.** The regenerated NZ benchmark results (with the `M0 trim` rung) are already in the main repository, so they are not re-run on this branch.

### Questions on Task 1
1. **Finite or superpopulation signs?** Should the sign assumptions be about the realised counties/polls (exact bounds, slack dial, no CI) or about behaviour (CLR intervals)? I propose offering both, with the finite version as the default for the actual rate.
2. **R×C generalisation of $\delta_W$:** EGGHM's one-vs-rest (W1), the Goodman-native vector (W2), or both?
3. **N2 scope:** is monotonicity of each row's conditional mean in the composition the nonparametric rung you want? Is a polynomial-τ implementation with CLR conditional on the degree acceptable for now, with fully nonparametric inference deferred? For $R\ge3$, monotone in which coordinates?
4. **Benchmark assumptions for NZ:** I would use CR (sign-free) and "own-party voters are more likely than others to vote for that party's candidate" (between) as the headline rungs, with within-signs in both directions as falsification checks. Agreed?
5. **API:** a `monotone=` modifier plus `model = "DD"` / `"MIV"`, or separate functions?

### Questions on Task 2
6. **Data.** Please share the veterans/suffrage data. Which states, units (counties? towns?) and referenda: pre-war, e.g. NY 1846/1860, WI 1849/1857, CT 1847, IA 1857; post-war 1865–1870? Is there more than one pre-war referendum per state (needed for placebo and δ calibration)?
7. **Veteran share.** Was it measured once? From what source (state census of 1865, 1890 veterans schedule, rosters)? On what base (post-war eligible men)? Is anything known about pre-war composition: age structure, war deaths, migration?
8. **Same units and base?** Are county boundaries and the population base (eligible men, with abstainers) the same in both periods? Which target is the headline: yes share of eligible men (linear), or yes share among voters (ratio)? Model-implied or actual-rate DiD?
9. **Identifying restriction.** Which of R1–R4 do you find substantively defensible for veterans? Is the causal (parallel-trends) claim in scope, or only the DiD of group rates?
10. **Regression test.** Can you share your earlier 2×2 JKST results so the new code can reproduce them?

---

## D. Decisions after first review (29 Oct 2026)

| # | Decision |
|---|---|
| 1 | **Exact (finite-population) bounds** for the sign rungs. They apply to the realised units, need no CI, and come with a slack dial and breakdown value. *Caveat:* "exact" is possible only for rungs without a fitted regression (N1, D0, D4). Rungs that use a regression (K-rows, N2, D1) still need P3 inference. |
| 2 | Open: W1 vs W2; explanation in D.1. |
| 3 | N2 accepted in principle; plain-language explanation in D.2 and in the guide. |
| 4 | No NZ benchmark for this work. What is needed is the **option** to impose contextual reinforcement. |
| 5 | API: `monotone=` modifier plus `model=` values (`"DD"`, `"MIV"`). |
| 6 | Data to come: **Wisconsin 1857/1865** and **Iowa 1857/1868** suffrage referenda. One pre-period per state, so no pre-period placebo; D3's δ cannot be calibrated from pre-trends and stays a sensitivity dial with breakdown values. |
| 7 | Veteran share is measured once. Soldiers who died can be removed. The 1860 census gives age structure; migration is hard. → B4 option 1 (interval pre-war composition) with deaths removed exactly and age eligibility as an interval; migration as a sensitivity bound. |
| 8 | Units are standardised to stable geographic/electoral units. **Population counts change between periods.** Headline target: **yes share of eligible men** (linear). Model-implied vs actual explained in D.3. |
| 9 | Open: plain-language R1–R4 in D.4. |
| 10 | Data and code to come. The earlier analysis applied JKST to the WI towns with near-zero 1857 YES votes, which narrowed the time-1 bounds (see D.6). |

### D.1 W1 vs W2 (matters only with three or more rows)

With two rows (veteran / non-veteran) the two definitions coincide. They differ once a third row appears, which happens here if we add an "entered" or "exited" row (B4 option 2).
* **W1 (one-vs-rest, the paper's).** For veterans: "do veterans' yes-rates, and the yes-rates of *everyone else taken together*, rise with the town's veteran share?" Simple, but "everyone else" is a mix whose make-up changes from town to town. Its sign is therefore harder to argue from substance.
* **W2 (each group separately).** For each group and each composition share: "does this group's yes-rate rise with that share?" Each sign is a statement about one identifiable group. It is also exactly what biases Goodman's multiple regression.

**Recommendation:** W1 as the default (it reproduces the paper), with W2 as an option once rows ≥ 3.

### D.2 N2 in plain language

"How a group votes may depend on who its neighbours are. N2 assumes this dependence runs **one way only**: for example, the more veterans in a town, the higher (never lower) the veterans' yes-rate, and likewise for non-veterans. It may be any shape: flat and then rising, a step, a curve. The linear rung instead assumes a straight line."
* **If wrong:** with a reversal (rates highest in middling towns, or backlash where a group is a large minority) the bound can miss the truth.
* **Check:** an empty bound refutes it. Plot town yes-share against veteran share and look for reversals.
* **Practical form:** the overall town pattern is fitted with a polynomial curve of data-chosen degree, and the uncertainty statement treats that curve's shape as correct.

### D.3 Model-implied vs actual DiD

* The **actual** DiD is the change in yes-rates among the actual men in these towns. It includes town-level idiosyncrasies (a local editor, a local controversy) that no model explains.
* The **model-implied** DiD is what the *systematic* part of behaviour (average behaviour of towns with a given veteran share) implies for these towns.
* For the exact, model-free rungs (decision 1) the target is automatically the actual DiD, and there is no sampling error.
* For rungs that fit a regression, the actual DiD's interval adds an allowance for town-level luck in **both** elections, which can be large. The model-implied interval has none.
* **Proposal:** report the actual DiD for model-free rungs, and the model-implied DiD as the headline for regression rungs (with the actual one alongside).

### D.4 R1–R4 in plain language (Wisconsin)

A group's **slope** is how its yes-rate differs between towns with many and few veterans. The war changed each group's yes-rate; the question is whether that *change* differed between veteran-heavy and veteran-light towns. Each rule fixes one aspect of this and turns the DiD into a single number. The rules generally give *different* numbers.

* **R1 (no differential spillover onto non-veterans).** Non-veterans changed by the same amount whether they lived among many or few veterans.
  * Rules out returning veterans persuading their neighbours *more where there were more of them*.
  * Still allows non-veterans everywhere to shift equally.
  * Fails if veterans persuaded or provoked neighbours in proportion to their numbers.
* **R2 (veterans' change does not depend on fellow veterans).** A veteran changed by the same amount whether he came home to a town of many or few veterans.
  * Plausible if the change happened in the army, not through veterans' networks at home (Grand Army posts, local veterans' politics).
  * Fails if veteran communities reinforced each other.
* **R3 (a technical symmetry).** The change among veterans in an imaginary town of non-veterans equals the change among non-veterans in an imaginary town of veterans.
  * It has no natural substantive reading. **I propose to drop it**, or to report it only as a technical benchmark.
  * *Correction to §B3:* R3 is not "contextual effects unchanged". That assumption is R1 and R2 together.
* **R4 (on average, no contextual change).** The weighted average of the two groups' slope changes is zero. If non-veterans moved more in veteran-heavy towns, veterans moved correspondingly less there.
  * Equivalent to a straight-line regression of the town's change in yes share on its veteran share (Goodman on changes).
  * Fails if both groups moved the same way with veteran share (reinforcement), which seems the most likely violation.
* **R1 + R2 together ("contextual effects unchanged", your earlier assumption).** Testable: the curvature of the town pattern must be the same in both elections.

**Suggested reporting:**
* the model-free bound, then each of R1, R2 and R4 as points;
* the bounded-change version of R1 (how big a differential spillover would be needed to overturn the conclusion);
* the R1+R2 test.

### D.5 New consequence: changing population counts (proved from 2a Lemma 2)

* Each period's rate depends on the unidentified part through weights built from that period's counts, $s_t=\mu_{vn,t}/\mu_{v,t}$ (and likewise for non-veterans).
* When counts change, the DiD's unidentified part is $s_1\kappa_1-s_0\kappa_0=s_1\Delta\kappa+(s_1-s_0)\kappa_0$.
* So **even R1–R4 no longer give a single number**. The DiD is bounded by a short interval whose width is proportional to $\lvert s_1-s_0\rvert$ times the width of the pre-war set for $\kappa_0$. That is small when the town mix of population growth is unrelated to veteran share.
* This also interacts with composition (B4): different counts mean partly different men.
* Phase 2 will give the exact form and a within-person (stayer) target for which the cancellation is exact.

### D.6 New point: selecting towns by their votes

* Restricting the analysis to towns with near-zero 1857 YES votes changes the population of the target: the rates refer to those towns.
* It also **selects on the outcome**. Even if contextual effects are linear in all towns, they need not be linear among towns chosen for low YES votes, so a model-based bound on the subset can be biased (*heuristic concern; to quantify in Phase 2*).
* The new tools get the same information without selection. In those towns the unit-level Duncan–Davis range is already nearly a point. The exact rungs (N1, D0, D4) use that automatically, and the `trim` option lets such towns enter through their own unit bounds.
* The earlier result will be reproduced both ways once the data arrive.
