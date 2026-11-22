# Phase 2e: a general theory of difference-in-differences with ecological data (descriptive and causal)

*Build step 2 of `NEXT_BUILD.md`. Proof-register items **R11a–d** and **R11b′** of `next_steps/NEXT_STEPS_RECOMMENDATION.md` (cited as NS) §3.2, with NS Part 3A (the combination matrix), decision D11 (T ≥ 2 elections with pre-periods; all three causal forms) and the owner's decisions of 22 Nov 2026 (changing populations, and uncertain margins and compositions, are registered, not solved). Builds on and generalises `PHASE2d_did.md` (2d) and `PHASE1b_monotone_did.md` (1b); uses `PHASE6a_count_layer.md` (6a: R19, R27, Lemma 2.1, Lemma C1, Lemma 4.1), `PHASE2c_monotone.md` (2c: Def. 1, Lemma M1, Props M2–M5), `PHASE5a_sign_module.md` (5a), `PHASE3_inference.md` (P3: §0, Thm P3.5, Thm P3.8), `PHASE2a_rxc.md`/`PHASE2b_districts.md` (2a, 2b) and `PACKAGE_DESIGN.md` (PD: §5.1 P-S, §5.2 P-R, §7, §8). Notation follows `notation.tex`; latent counts are written $n_{irct}$ as in 6a.*

*Theory only: no package code. Nothing in other documents is edited; problems found in them are listed in §11. Checks: `verify_phase2e.py` (output `verify_phase2e_output.txt`, **71/71 PASS**), which runs 2d's 15 checks unchanged and adds 56, tagged as in 2c, 5a and 6a: **[F]** fails on an earlier or tempting ("naive") statement and passes on the one made here; **[V]** verifies a statement made here; **[I]** illustrates.*

*Status in brief.* Every result marked **P; proof review pending** has a numbered-step proof here and a numerical check; none has had the independent review that PD §8's standard requires. **C** marks conjectures with numerical support; **O** marks open items. The descriptive layer (§§2–4) is complete for the fixed population: exact count-layer rungs on the T-period product (Δ0, Δ4 in R×C, contrasts on changes), the sharp identified set under FD, DiD-S sharpness, and the combination rule. The causal layer (§5) covers the three forms the owner asked for (trend dial; relative magnitudes and smoothness over the **joint** set of latent tables, with the degeneration proposition; MTR/ND/PS), a spillover dial, and the labelling rule. Inference (§6) covers unions and FD and gives the joint-region result. Integer versions (§7) cover margins, one-vs-rest rows on changes and change rows on one column. Changing populations and uncertain margins are collected in the gap register (§8). §9 says in plain language what can be said under which assumptions.

---

## 0. Tags used here

Tags of 2a §0.3, 2b §0.3, 2c §0.3, P3 §0.3, 2d §0.3 and 6a §0 apply, plus:

| Tag | Meaning |
|---|---|
| [stack] | Lemma E1: under the fixed population the T-period problem is a count-layer problem whose units are the unit-periods $(i,t)$ |
| [lin-img] | the support function of a linear image is $h_{A\mathcal S}(u)=h_{\mathcal S}(A^\top u)$; the linear image of a compact convex set is compact and convex |
| [Aumann] | for an integrably bounded random compact convex set, $E[h_{\mathcal S}(u)]=h_{E\mathcal S}(u)$ (Aumann expectation), and $E\mathcal S\subseteq\{\text{means of selections}\}$; the expectation of a Minkowski sum is the Minkowski sum of expectations |
| [param] | if a model is "model $\mathcal M$ and $\theta\in\mathcal C$" with $\mathcal C$ a restriction on the parameter only, its identified set is $\Theta_I(\mathcal M)\cap\mathcal C$ |
| [sep] | a target that is a sum of functions of independent factors of a product set has range equal to the Minkowski sum of the factor ranges |
| [GH-depth] | 6a R19(b) proof step 2: in a laminar family, signing each chosen set by $(-1)^{\text{depth}}$ gives every element a total coefficient in $\{0,1\}$ |
| [PO] | potential outcomes: $Y_{jt}(0)$ is person $j$'s outcome at election $t$ had nobody been treated |

---

## 1. Setting

### 1.1 Units, groups, outcomes and elections

* **Units** $i=1,\dots,p$ (precincts, townships), with $N_i>0$ people.
* **Groups** (rows) $r=1,\dots,R$: a fixed characteristic of each person (veteran or not; a party registration recorded once; a census category).
* **Outcomes** (columns) $c=1,\dots,C$ on a **common base** (for example yes, no, abstain among the eligible). A share among voters is a ratio target, handled by Charnes–Cooper as in 2d C1.
* **Elections** $t\in\mathcal T=\{-L,\dots,-1,0,1,\dots,P\}$ in event time, $T=L+1+P\ge2$ (decision D11).
  * Election $0$ is the **base**: the last election before treatment.
  * Elections $t<0$ are earlier **pre-periods** ($L\ge0$; relative-magnitude and smoothness bounds need $L\ge1$).
  * Elections $t\ge1$ are **post-periods**.

### 1.2 The fixed population (the base case)

**(FP)** The same $N_i$ people live in unit $i$ at every election, and each keeps the same group. So the group counts $a_{ir}$ ($\sum_ra_{ir}=N_i$) are the same in every period, while the outcome counts $b_{ict}$ ($\sum_cb_{ict}=N_i$) vary.

Every departure from (FP) (entry, exit, aging-in, interval bases, changing boundaries, changing group membership) is an item of the gap register (§8), not a result of this note.

**Unknowns.** The realised counts $n_{irct}\ge0$ with $\sum_cn_{irct}=a_{ir}$ and $\sum_rn_{irct}=b_{ict}$. Write $\mathcal P_{it}$ for unit $i$'s transportation polytope at election $t$ (6a §1.1), $\mathcal P^{\mathcal T}:=\prod_t\prod_i\mathcal P_{it}$, and $n=(n_{irct})$. Rates: $\beta_{irct}=n_{irct}/a_{ir}$ for $a_{ir}>0$.

**Lemma E1 (stacking).** Under (FP), $\mathcal P^{\mathcal T}$ is the count-layer polytope of 6a §1.1 for the unit set $\{(i,t)\}$, with row margins $a_{ir}$ (the same for every $t$) and column margins $b_{ict}$. Hence 6a R19(a) and R19(d) apply verbatim to every linear target and every finite set of linear rows on $n$, **including rows that couple periods**: the set of target values is $\varphi(F)$, an interval when $F\ne\emptyset$, valid, and sharp for continuous tables given the margins and rows; CR-type unions give unions of intervals.

*Proof.*
1. Election $t$'s margins constrain only $n_{\cdot\cdot\cdot t}$; (FP) makes the row margins of every $(i,t)$ equal to $a_{ir}$. [Def; FP]
2. So the feasible set is the product of the unit-period transportation polytopes, which is 6a's $\mathcal P$ for the unit set $\{(i,t)\}$. [1; 6a §1.1]
3. 6a R19(a) items 1–4 and R19(d) use only that $\mathcal P$ is a product of transportation polytopes and that the rows are linear. [6a:R19(a),(d)] ∎

*Check.* E1 (Minkowski sums over periods; 40 designs, R = C = 3, T = 3). **Status: P; proof review pending.**

**Lemma E2 (within-person tables add nothing for period-margin targets).** Let $\nu_{irc_{-L}\dots c_P}\ge0$ be the panel table (people in unit $i$, group $r$, with outcome sequence $(c_t)_t$). Every family of period tables $n\in\mathcal P^{\mathcal T}$ is the set of one-period marginals of some panel table, and of an **integer** panel table when $n$ is integer. Hence, for targets and rows that depend on the panel only through the period tables (all targets of this note), the sets are the same whether the unknown is $\nu$ or $n$. Transition targets ("share of group $r$ that switched to $c$") are not of this kind and need $\nu$.

*Proof.*
1. Real case: $\nu_{irc_{-L}\dots c_P}:=a_{ir}\prod_t(n_{irc_tt}/a_{ir})$ for $a_{ir}>0$ (independent coupling given $(i,r)$) has the period tables as marginals. [coupling]
2. Integer case: list the $a_{ir}$ people of $(i,r)$; at each election give the first $n_{ir1t}$ of them outcome 1, the next $n_{ir2t}$ outcome 2, and so on (in any order per election). Counting outcome sequences gives an integer $\nu$ with the right marginals. [Def; alg]
3. A target or row that depends on $\nu$ only through its marginals takes the same values on $\{\nu\}$ and on $\mathcal P^{\mathcal T}$. [1–2] ∎

*Check.* E2 (200 integer instances, T = 3). Generalises 1b B2 Δ8 (two periods, a sketch there). **Status: P; proof review pending.**

### 1.3 Targets

For a unit set $\mathcal S$ (default: all units), with $Y_{rct}(\mathcal S)=\sum_{i\in\mathcal S}n_{irct}$ and $A_r(\mathcal S)=\sum_{i\in\mathcal S}a_{ir}$:
* **group rates** $B_{rct}(\mathcal S)=Y_{rct}(\mathcal S)/A_r(\mathcal S)$ (P3 §0.2's actual rate);
* **gaps** $g_{rs,ct}=B_{rct}-B_{sct}$; one-vs-rest gaps $g_{r\neg r,ct}$ with $B_{\neg r,ct}=(Y_{ct}-Y_{rct})/(N-A_r)$;
* **DiDs** $D_c(r,s;u,t)=g_{rs,ct}-g_{rs,cu}$; the **event-study vector** $D_t:=g_{rs,ct}-g_{rs,c0}$ ($t\ne0$), relative to the base; **consecutive DiDs** $g_{rs,ct}-g_{rs,c,t-1}$.

All are linear count targets with coefficients that are functions of the design ($1/A_r(\mathcal S)$), as 6a Lemma 4.1 and R27 require. The **model-implied** versions $\tilde B$, $\tilde D$ replace $\beta_{irct}$ by $m_{rct}(V_i)=E[\beta_{irct}\mid V_i]$ (P3 §0.2; 2d §0.1). Under (FP) each person's **individual change** $\Delta Y_j(u,t)=\mathbf 1\{Y_{jt}=c\}-\mathbf 1\{Y_{ju}=c\}\in\{-1,0,1\}$ is defined, and $\Delta\beta_{irc}(u,t)=\beta_{irct}-\beta_{ircu}$ is the mean of the individual changes of group $r$ in unit $i$.

### 1.4 What is descriptive and what is causal

* A **descriptive** statement is about the realised tables of the actual groups (count layer) or about the model-implied rates of the actual groups (model layer): "veterans' support rose 5 to 25 points more than others'".
* A **causal** statement is about potential outcomes (§5): "service raised veterans' support by …". It is always **on top of** (↑) a descriptive set, never mixed into it.
* **Rule (Definition 5.1 below).** A causal label requires a stated treated group and comparison group, a stated trend assumption, a stated spillover assumption and no anticipation. Without a trend assumption the descriptive DiD says nothing causal (Prop. E15).

---

## 2. Descriptive rungs, count layer (exact)

### 2.1 Δ0 for T elections

**Proposition E3 (Δ0: no link between elections).** Margins only, (FP).
(a) For a target that is a sum of per-election targets, $\varphi(n)=\sum_t\varphi_t(n_{\cdot\cdot\cdot t})$, the set of values is the Minkowski sum $\sum_t\varphi_t(\mathcal P_t)$, every value attained (continuous tables).
(b) For a two-row, or one-vs-rest, one-column DiD, the set is $[G^-_t-G^+_u,\ G^+_t-G^-_u]$, where $[G^-_t,G^+_t]$ is election $t$'s exact gap range: $[B^L_{rct}-B^U_{sct},\ B^U_{rct}-B^L_{sct}]$ with cellwise Duncan–Davis ends for a named pair (6a Lemma 2.1(b): both ends are attained jointly), and the gap evaluated at $Y_{rct}=Y^L_{rct},Y^U_{rct}$ for one-vs-rest (the gap is increasing in $Y_{rct}$). This holds for any $R$, $C$ and $T$.
(c) **The joint set of the event-study vector is not the box of the separate ranges.** For a named pair and one column, the joint set of $(D_t)_{t\ne0}$ is $\{(g_t-g_0)_{t\ne0}:\ g_t\in[G^-_t,G^+_t]\ \forall t\}$, the image of a box under a linear map. Whenever $[G^-_0,G^+_0]$ has positive length and at least two other elections have ranges of positive length, some corner of the box of separate ranges is infeasible (for example $D_1$ at its maximum, which needs $g_0=G^-_0$, with $D_{-1}$ at its minimum, which needs $g_0=G^+_0$). Consecutive DiDs share a period with opposite signs, so the same holds for them.
(d) **Per-period LCE.** Under per-period LCE (M2 with elections as groups) the sharp identified set of $(K_t)_{t\in\mathcal T}$ is $\prod_t\Theta_{I,t}$, each factor computed from its own election's law; the outer set is $\prod_t\Theta_{O,t}$.

*Proof.*
1. (a): $\mathcal P^{\mathcal T}$ is a product over elections [E1], and $\varphi$ is a sum of functions of the factors. [sep] Each factor's range is an interval with all values attained [6a:R19(a)], so the sum's range is the Minkowski sum, all values attained. [conn]
2. (b): $D=g_t-g_u$ is of the form (a). Election $t$'s gap range is exact by 6a Lemma 2.1(b) (named pair; the separate optima are attained by one table) or by monotonicity in $Y_{rct}$ (one-vs-rest, $B_{\neg r}$ decreasing in $Y_{rct}$) together with 6a Lemma 2.1(a). [6a:L2.1; 1]
3. (c): $(g_t)_t$ ranges over the box $\prod_t[G^-_t,G^+_t]$ [1–2], and $(D_t)$ is its image under $g\mapsto(g_t-g_0)_{t\ne0}$. [Def] At the named corner $D_1=G^+_1-G^-_0$ forces $g_0=G^-_0$ and $D_{-1}=G^-_{-1}-G^+_0$ forces $g_0=G^+_0$, incompatible when $G^-_0<G^+_0$. [alg]
4. (d): 2d Δ0(b)'s proof with the $T$-fold product $\mathcal S(w)=\prod_tS(x,t_t)$ in place of the twofold one: Lemma G (2d §1) applies because products preserve its hypotheses; the support function of a product is the sum [prod-supp]; setting all but one $u_t$ to 0 gives each election's condition, and summing gives the joint one. [2d:L G; 2d:Δ0(b)] ∎

*Checks.* E1 [V] (a), (b) (40 designs, R = C = 3, T = 3); **E1 [F]** (c) (the box corner infeasible in all 40 designs). [P2d] Δ0's own checks run unchanged. **Status: P; proof review pending.** (a)–(b) generalise 2d Δ0(a) to R×C and T elections; (d) generalises 2d Δ0(b).

### 2.2 Rows on changes

**Lemma E4 (rows on changes are rows on the product).** (FP). Fix a column $c$, elections $u<t$, a unit set $\mathcal S$ with $\operatorname{Var}_p(X_r)>0$ over $\mathcal S$ (weights $p_i=N_i/N(\mathcal S)$), and write $\Delta$ for "election $t$ minus election $u$".
(a) $\Delta T_{ic}=\sum_rX_{ir}\,\Delta\beta_{irc}$, and the realised associations of 2c Def. 1 computed on $\Delta\beta$ (pairwise $\delta^\Delta_B(r,s)$, one-vs-rest $\delta^{\Delta,r}_B$, W2 $\Delta^\Delta$, W1 $\delta^{\Delta,r}_W$) are affine functions of the two-election table with design coefficients; each equals election $t$'s association minus election $u$'s.
(b) **One-vs-rest rows depend on one number.** $\delta^{\Delta,r}_B(\mathcal S)$, $\delta^{\Delta,r}_W(\mathcal S)$ and the one-vs-rest DiD $D^r_c(u,t;\mathcal S)$ are affine in the single aggregate
$$Z:=Y_{rct}(\mathcal S)-Y_{rcu}(\mathcal S),$$
with slopes $1/N(\mathcal S)$, $-\operatorname{Var}_p(X_r)\kappa_r$ and $\kappa_r:=1/A_r(\mathcal S)+1/(N(\mathcal S)-A_r(\mathcal S))$ respectively, and intercepts that depend on both elections' outcome margins.
(c) **Identities on changes.** With $ER^r_\Delta:=\operatorname{Cov}_p(\Delta T_c,X_r)/\operatorname{Var}_p(X_r)$ and $NM^r_\Delta:=\operatorname{Cov}_p(\Delta T_c,X_r)/(\mu_r(1-\mu_r))$ (the collapsed neighbourhood-model difference $NM_r-NM_{\neg r}$ computed on $\Delta T$),
$$ER^r_\Delta-D^r_c=\frac{\delta^{\Delta,r}_W}{\operatorname{Var}_p(X_r)},\qquad NM^r_\Delta-D^r_c=-\frac{\delta^{\Delta,r}_B}{\mu_r(1-\mu_r)},\qquad NM^r_\Delta=\gamma_r\,ER^r_\Delta,$$
with $\gamma_r=\operatorname{Var}_p(X_r)/(\mu_r(1-\mu_r))\in(0,1]$.

*Proof.*
1. (a): A1 holds in both elections with the same $X_i$ [FP], so $\Delta T_{ic}=\sum_rX_{ir}\Delta\beta_{irc}$. 2c Lemma M1(d) shows the associations are affine in the unit rates with design coefficients, for any bounded rates; apply it to each election and subtract. $\beta_{irct}=n_{irct}/a_{ir}$ is linear in $n$. [A1; 2c:M1(d); lin]
2. (b): $\delta^r_B=\sum_{s\ne r}E_p[X_rX_s(\beta_r-\beta_s)]=E_p[X_r(\beta_r-T_c)]$, using $\sum_{s\ne r}X_s\beta_s=T_c-X_r\beta_r$; in counts this is $(Y_{rc}(\mathcal S)-\sum_{i\in\mathcal S}X_{ir}b_{ic})/N(\mathcal S)$. [Def; simplex] $D^r$ is $\kappa_rY_{rc}-Y_c/(N-A_r)$ per election. $\delta^r_W=\operatorname{Var}_p(X_r)(ER^r-D^r)$ by 2c Lemma M1(c), with $ER^r$ a function of the data only. Subtract the two elections. [2c:M1(c); alg]
3. (c): the collapsed population $(X_r,\Delta\beta_r,\Delta\beta^{\neg r})$, with $\Delta\beta^{\neg r}_i=\sum_{s\ne r}X_{is}\Delta\beta_{is}/(1-X_{ir})\in[-1,1]$, is a two-group population for the bounded outcome $\Delta Y$ in EGGHM's sense (2c M1(c) proof steps 1–5 with $\Delta$; only A1 and boundedness are used). Apply [EGGHM:P4, P6]. For the last identity: $E_p[X_r\Delta T]/\mu_r-E_p[(1-X_r)\Delta T]/(1-\mu_r)=\operatorname{Cov}_p(X_r,\Delta T)/(\mu_r(1-\mu_r))$. [2c:M1(c); EGGHM:P4, P6; alg] ∎

*Check.* E4 (60 designs, R = C = 3: the identities from the definitions of 2c Def. 1 on random latent tables, error < $10^{-10}$). **Status: P; proof review pending.** Generalises 2d Δ4(a) (2×2) and the one-vs-rest part of 6a Lemma C1 to changes.

### 2.3 Δ4 beyond 2×2

**Proposition E5 (sign rows and CR on changes, R×C, T elections).** (FP), continuous tables, margins only unless rows are named. Fix $(r,c,u,t,\mathcal S)$ with $\operatorname{Var}_p(X_r)>0$; let $[Z^-,Z^+]=[Y^L_{rct}(\mathcal S)-Y^U_{rcu}(\mathcal S),\ Y^U_{rct}(\mathcal S)-Y^L_{rcu}(\mathcal S)]$ (Duncan–Davis extremes of each election's aggregate, 6a Lemma 2.1(a)) and $D^\pm=D^r_c(Z^\pm)$.
(a) The set of $D^r_c(u,t;\mathcal S)$ is $[D^-,D^+]$, every value attained.
(b) **W1 signs on changes.** $\delta^{\Delta,r}_W\ge0$ gives $[D^-,\min\{D^+,ER^r_\Delta\}]$, empty iff $ER^r_\Delta<D^-$; $\delta^{\Delta,r}_W\le0$ gives $[\max\{D^-,ER^r_\Delta\},D^+]$, empty iff $ER^r_\Delta>D^+$; $\delta^{\Delta,r}_B\ge0$ gives $[\max\{D^-,NM^r_\Delta\},D^+]$; $\delta^{\Delta,r}_B\le0$ gives $[D^-,\min\{D^+,NM^r_\Delta\}]$.
(c) $NM^r_\Delta\in[D^-,D^+]$ always, so between signs on changes are never refuted in this problem.
(d) **CR on changes** ($\delta^{\Delta,r}_W\,\delta^{\Delta,r}_B\ge0$): the set is $[NM^r_\Delta,\min\{D^+,ER^r_\Delta\}]$ if $ER^r_\Delta\ge0$ and $[\max\{D^-,ER^r_\Delta\},NM^r_\Delta]$ if $ER^r_\Delta<0$. It is never empty, and if $\gamma_r<1$ and $ER^r_\Delta\ne0$ exactly one branch is non-empty and $\operatorname{sign}D^r_c=\operatorname{sign}ER^r_\Delta$.
(e) **Other rows.** Pairwise $\delta^\Delta_B(r,s)$, W2, signs in several columns, several unit sets, unit-level orderings of changes ($\Delta\beta_{irc}\ge\Delta\beta_{isc}$, a two-election row in unit $i$), known cells and published totals in any election, and within-election rows: the set is the LP of Lemma E1 (sharp for continuous tables); CR is a union of two LPs.
(f) **Several period pairs.** Rows for pairs $(u_k,t_k)$, $k=1,\dots,K$, are rows on the same product: one LP; CR in each of $K$ pairs is a union of at most $2^K$ branch LPs, and the set is the union of their intervals.

*Proof.*
1. (a): $Z=Y_{rct}(\mathcal S)-Y_{rcu}(\mathcal S)$, the two terms in different factors of the product [E1], each taking every value in its Duncan–Davis interval [6a:L2.1(a)], so $Z$ takes every value in $[Z^-,Z^+]$. [sep; conn] $D^r_c$ is affine increasing in $Z$. [E4(b)]
2. (b): by E4(c), with $\operatorname{Var}_p(X_r)>0$: $\delta^{\Delta,r}_W\ge0\iff D^r_c\le ER^r_\Delta$ and $\delta^{\Delta,r}_B\ge0\iff D^r_c\ge NM^r_\Delta$. Intersect with (a). [E4(c); conn]
3. (c): put $\beta_{isct'}=T_{ict'}$ for every row $s$ and both elections $t'\in\{u,t\}$. This lies in each $S(X_i,T_{it'})$; then $\Delta\beta_{is}=\Delta T_i$ for every $s$, so $\delta^{\Delta,r}_B=E_p[X_r(\Delta\beta_r-\Delta T)]=0$ and $D^r_c=NM^r_\Delta$ by E4(c). [Def S; E4]
4. (d): $NM^r_\Delta=\gamma_rER^r_\Delta$ with $\gamma_r\in(0,1]$ [E4(c); 2c:M3 step 5]. If $ER^r_\Delta\ge0$ then $NM^r_\Delta\le ER^r_\Delta$ and the positive branch is $[NM^r_\Delta,\min\{D^+,ER^r_\Delta\}]$, non-empty since $D^-\le NM^r_\Delta\le D^+$ [3]. The negative branch needs $ER^r_\Delta\le D\le NM^r_\Delta=\gamma_rER^r_\Delta$, which for $ER^r_\Delta>0$ forces $\gamma_r=1$. The case $ER^r_\Delta<0$ is symmetric; on the non-empty branch $D$ has the sign of $ER^r_\Delta$. [2c:M3 steps 6–8 with $Z$]
5. (e), (f): every row named is affine in $n$ [E4(a); Def]; apply E1 [stack; 6a:R19(a) item 4 for unions]. ∎

*Checks.* E5 [V] (a)–(b): closed forms = LP over the R×C two-election product (360 LPs, R = C = 3); (c)–(d) (NM in range, $NM=\gamma ER$, $\gamma<1$); (e) pairwise and joint-column rows valid and able to tighten (pairwise sign tightened 40/40; a second column's sign 2/40); (f) T = 3, CR in two pairs, union of four branch LPs contains the truth (60 designs). E3 [V]: brute force over integer two-election 2×2 tables for the margins, each W1 sign and CR. **Status: P; proof review pending.** This settles PD §8's R11 item "Δ4 beyond 2×2" for the one-vs-rest rows (closed form) and for all other rows (LP).

**Reading.** 2d's reading carries over to R×C. $\delta^{\Delta,r}_B$ is the $X_r(1-X_r)$-weighted average of within-unit DiDs, so signing it presupposes the answer; $\delta^{\Delta,r}_W\ge0$ says the net contextual change runs towards $c$ where group $r$ is common. Realised CR on changes is a belief about the actual voters (`object = "actual"`); model CR on changes (Prop. E11) is a different belief (NS Part 3A footnote 2).

### 2.4 Finite contrasts on changes (the count-layer form of DiD-S)

**Proposition E6 (contrast on changes, two bins).** (FP), $R=2$ (rows 1, 2), a column $c$, elections $u<t$, and two prespecified disjoint unit sets ("bins") $\mathcal L,\mathcal H$ with pooled compositions $\bar x_b=A_1(b)/N(b)$ and $\lambda:=\bar x_{\mathcal H}-\bar x_{\mathcal L}>0$. Write $\Delta m_{r,b}=\Delta Y_{rc}(b)/A_r(b)$ (bin change rates), $\Delta\tau_b=\Delta Y_c(b)/N(b)$ (observed), and the **bin-averaged DiD** $\bar D:=\tfrac12\sum_{b\in\{\mathcal L,\mathcal H\}}(\Delta m_{1,b}-\Delta m_{2,b})$.
(a) **Identity** (Kitagawa on changes): $\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L}=\lambda\bar D+\rho$, with the "rate" term $\rho:=\bar{\bar x}(\Delta m_{1,\mathcal H}-\Delta m_{1,\mathcal L})+(1-\bar{\bar x})(\Delta m_{2,\mathcal H}-\Delta m_{2,\mathcal L})$, $\bar{\bar x}=(\bar x_{\mathcal H}+\bar x_{\mathcal L})/2$.
(b) The two beliefs "mix $\ge0$" ($\lambda\bar D\ge0$) and "rate $\ge0$" ($\rho\ge0$) are linear count rows, and together they are equivalent to $\bar D\in[0,(\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L})/\lambda]$. The set is
$$\big[\max\{0,\bar D^-\},\ \min\{\bar D^+,(\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L})/\lambda\}\big],$$
with $[\bar D^-,\bar D^+]$ the margins-only range of $\bar D$; it is empty iff $\Delta\tau_{\mathcal H}<\Delta\tau_{\mathcal L}$ (the two beliefs are jointly refuted) or the cap is below $\bar D^-$. Sharp for continuous tables.
(c) The rate belief alone gives $\bar D\le(\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L})/\lambda$ and does **not** presuppose the sign of the DiD; the mix belief alone is a sign on the bin-averaged DiD itself.

*Proof.*
1. (a): $\Delta\tau_b=\bar x_b\Delta m_{1,b}+(1-\bar x_b)\Delta m_{2,b}$ (A1 in both elections, (FP), aggregated over the bin). Apply $a_{\mathcal H}v_{\mathcal H}-a_{\mathcal L}v_{\mathcal L}=(a_{\mathcal H}-a_{\mathcal L})\frac{v_{\mathcal H}+v_{\mathcal L}}2+\frac{a_{\mathcal H}+a_{\mathcal L}}2(v_{\mathcal H}-v_{\mathcal L})$ to $(\bar x,\Delta m_1)$ and $(1-\bar x,\Delta m_2)$ and add; the first terms give $\lambda\bar D$. [A1; FP; alg]
2. (b): $\Delta\tau_b$ is observed, so $\rho=\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L}-\lambda\bar D$ and both rows are rows on $\bar D$ alone. [1] $\bar D$ is affine with positive coefficients in $(\Delta Y_{1c}(\mathcal H),\Delta Y_{1c}(\mathcal L))$ (using $\Delta Y_{2c}(b)=\Delta Y_c(b)-\Delta Y_{1c}(b)$), whose margins-only joint range is a product of intervals with all values attained [E1; 6a:L2.1(a); sep]. So $\bar D$'s range is an interval with all values attained; intersect. [conn]
3. (c): (b) with one row. ∎

*Checks.* **E3 [V]** brute force (integer 2×2 tables, the identity exact, the integer set as stated); **E6 [V]** LP = closed form in 2×3 (59 designs; empty in 32, where the joint test rejects). **Status: P; proof review pending.** This is the count-layer version of NS §3.3's contrast "on changes" and of PD §7's DiD-S; it needs no linearity in either election. With $R\ge3$ only named pairs should be used (NS §3.3).

---

## 3. Descriptive rungs, model layer

The model layer targets the **model-implied** DiD $\tilde D$; the actual DiD is reached through the residual allowance (Lemma E13) and the combination rule (§4). Model **FD** (2d §3): (S1) in each election, and
$$E[\Delta\beta\mid V]=M_{\Delta H}(X)\quad\text{a.s.},\qquad \textstyle\sum_c\Delta H_c=0,$$
for one pair of elections (or each pair used); **levels are unrestricted**. Write $\Delta H_c=\Delta Q_c+\Delta K_c$ (symmetric plus antisymmetric part).

### 3.1 FD in R×C

**Lemma E7 (FD in R×C).** Assume FD, (FP), and A5 for the design (2a).
(a) $E[\Delta T_c\mid X]=X^\top\Delta Q_cX$; so $\Delta Q_c$ is identified and $\Delta K_c$ is not restricted by this regression.
(b) For a pair $(r,s)$ (or $r$ against the rest) and column $c$, $\tilde D_c=I_c+\ell_c(\Delta K_c)$, where $I_c$ is identified and $\ell_c$ is linear with design coefficients and involves only column $c$'s $\Delta K_c$.
(c) For $R=2$: $\ell_c=s\,\Delta\kappa_c$ with $s=\mu_{vn}(1/\mu_v+1/\mu_n)>0$ (2d Lemma FD1(b)).
(d) A pair DiD loads on the directions $K_{c,rt}-K_{c,st}$; pairwise restrictions on a spanning tree of groups reach all directions (PD Prop. P-R: proved for star trees; chains checked, R5).

*Proof.*
1. (a): $\Delta T_c=\sum_rX_r\Delta\beta_{rc}$ [A1 both elections; FP], so $E[\Delta T_c\mid X]=X^\top\Delta H_cX=X^\top\Delta Q_cX$, since antisymmetric parts vanish in quadratic forms. $\Delta Q_c$ is unique under A5. [CE-pull; 2a:L1]
2. (b): $\tilde D_c=\sum_iw_{ir}(\Delta H_cX_i)_r-\sum_iw_{is}(\Delta H_cX_i)_s$ with sample weights; substitute $\Delta H_c=\Delta Q_c+\Delta K_c$. [2a:L2; P3 §0.2; lin]
3. (c): 2d FD1(b). (d): PD §5.2. ∎

*Check.* E7 (R = 3: $\Delta Q$ recovered exactly with non-linear levels; adding an antisymmetric $\Delta K'$ leaves $E[\Delta T\mid X]$ unchanged and moves $\tilde D$ by the predicted linear amount). 2d's single-restriction results (FD2: R1, R2, R4; FD3: bounded change and breakdown) apply unchanged to each column in 2×C (checks [P2d] Δ1). **Status: P; proof review pending** for (a)–(c); (d) as in PD §8 R5.

### 3.2 The sharp identified set under FD

2d FD4 gave an **outer** set. The sharp set needs Lemma G for a linear image.

**Lemma G′ (2d Lemma G for a linear image).** Let $\mathcal S$ and $m_\theta$ be as in 2d Lemma G and $A$ a fixed linear map. In the model "$\xi\in\mathcal S(W)$ a.s. and $E[A\xi\mid V]=m_\theta(V)$ a.s.", $\theta$ is in the sharp identified set iff, for $P_V$-a.e. $v$ and all $u$,
$$\langle u,m_\theta(v)\rangle\le\int h_{\mathcal S}(A^\top u;w)\,P_{W\mid V=v}(dw).$$

*Proof.*
1. $A\mathcal S(w)$ is nonempty, compact and convex [lin-img]; it is uniformly bounded; its graph is closed: if $w_k\to w$ and $A\xi_k\to\zeta$ with $\xi_k\in\mathcal S(w_k)$, a subsequence $\xi_k\to\xi\in\mathcal S(w)$ (uniform bound; closed graph), so $\zeta=A\xi$. Its support function is $h_{\mathcal S}(A^\top u;w)$. [lin-img; Def]
2. The model is equivalent to "$\zeta\in A\mathcal S(W)$ a.s. and $E[\zeta\mid V]=m_\theta(V)$". ⇒: take $\zeta=A\xi$. ⇐: the correspondence $w\mapsto\{\xi\in\mathcal S(w):A\xi=\zeta(w)\}$ has nonempty closed values and is measurable (the intersection of a closed-graph correspondence with the preimage of a measurable function), so it has a Borel selector $\xi$. [KRN]
3. Apply 2d Lemma G to $\zeta$ and $A\mathcal S$. [1–2; 2d:L G] ∎

**Proposition E8 (sharp set under FD; R11a).** Under FD for a pair $(u,t)$, take $\xi=(\beta_u,\beta_t)$, $\mathcal S(w)=S(x,T_u)\times S(x,T_t)$, $A(\beta_u,\beta_t)=\beta_t-\beta_u$ and $m_\theta=M_{\Delta Q+\Delta K(\Delta k)}$, with $\Delta Q$ identified (E7).
(a) The sharp identified set of $\Delta k$ is
$$\Theta^{FD}_I=\Big\{\Delta k:\ \langle u,M_{\Delta H}(x)\rangle\le E\big[h_S(u;X,T_t)+h_S(-u;X,T_u)\,\big|\,V=v\big]\ \ \forall u,\ P_V\text{-a.e. }v\Big\},$$
convex; the sharp set of $\tilde D_c$ is its image, an interval.
(b) $\Theta^{FD}_I\subseteq\Theta^{FD}_O$ (2d FD4).
(c) Equality holds when $T_u$ and $T_t$ are degenerate given $V$ (no dispersion of outcomes among units with the same composition). Strict inclusion is the rule when outcomes vary within a composition and cross the kinks of the unit bounds (check: strictly narrower in 100 of 100 designs).
(d) **2×2 closed form.** With $0<x<1$, $L(t;x)=\max\{0,(t-1+x)/x\}$, $U(t;x)=\min\{1,t/x\}$ and $\Delta m_v(x;\Delta\kappa)=\Delta Q_{vv}x+(\Delta Q_{vn}+\Delta\kappa)(1-x)$:
$$\Theta^{FD}_I=\Big\{\Delta\kappa:\ \Delta m_v(x;\Delta\kappa)\in\big[E(L(T_t;x)\mid x)-E(U(T_u;x)\mid x),\ E(U(T_t;x)\mid x)-E(L(T_u;x)\mid x)\big]\ \text{for }P_X\text{-a.e. }x\Big\},$$
an intersection of intervals in $\Delta\kappa$ (the coefficient $1-x$ is positive). FD4 is the same with $L,U$ evaluated at $\tau_t(x)=E[T_t\mid x]$.

*Proof.*
1. (a): Lemma G′ with $h_{\mathcal S}(A^\top u;w)=h_S(u;x,t_t)+h_S(-u;x,t_u)$, since $A^\top u=(-u,u)$. [G′; prod-supp] Linear inequalities in $\Delta k$ define a convex set; $\tilde D_c$ is affine in $\Delta k$ [E7(b)]. [conn]
2. (b): the constraint set $S(x,t)$ is defined by linear equalities $\beta^\top x=t$ and convex constraints, so the mean of any selection of $S(x,T)$ given $x$ lies in $S(x,\tau(x))$; hence $E[h_S(u;x,T)\mid x]\le h_S(u;x,\tau(x))$ for all $u$, for both elections. [Aumann; Def S] So every $\Delta k$ in $\Theta^{FD}_I$ satisfies FD4's pointwise condition $M_{\Delta H}(x)\in S(x,\tau_t)\oplus(-S(x,\tau_u))$. [2d:FD4]
3. (c): with $T$ degenerate given $V$ the conditional expectation of $h_S$ is $h_S$ at the conditional mean. [Def]
4. (d): $S(x,t)$ is the segment of $\{x\beta_v+(1-x)\beta_n=t\}\cap[0,1]^2$ with $\beta_v\in[L,U]$. $S(x,t_t)\oplus(-S(x,t_u))$ is a segment on the line $x\Delta_v+(1-x)\Delta_n=t_t-t_u$ with $\Delta_v\in[L(t_t)-U(t_u),U(t_t)-L(t_u)]$. For segments $[p,q]$ on parallel lines with direction $d$, $q-p=\ell d$, $\ell\ge0$, the support function is $\langle u,p\rangle+\ell\max\{0,\langle u,d\rangle\}$; its conditional expectation is the support function of $[Ep,Ep+(E\ell)d]$. $M_{\Delta H}(x)$ lies on the mean line (its $v$- and $n$-coordinates satisfy $x\Delta m_v+(1-x)\Delta m_n=E[\Delta T\mid x]$ for every $\Delta\kappa$, E7(a)), so the condition of (a) is the stated interval for $\Delta m_v$. [1; alg] ∎

*Checks.* E8 [V] (discrete design, 2×2: the sharp set contains the truth, lies inside FD4, strictly in 100/100 designs); E8 [I] (degenerate case). **Status: P; proof review pending.** Computing $\Theta^{FD}_I$ from data needs conditional expectations of the kinked functions $L(T),U(T)$ given $X$ (§6, Prop. E23(c)).

### 3.3 DiD-S within FD (R11a)

**Proposition E9 (DiD-S: closed form and sharpness).** Under FD, on a prespecified set $A$ of compositions (the hull, or chosen points), impose the PD §5.1 objects on changes: the ordering of changes $O^\Delta_{rs}(x)=(\Delta Hx)_r-(\Delta Hx)_s\ge0$ and the composition-weighted spillover sign $L^\Delta_{rs}(x)=(x^\top\Delta H)_r-(x^\top\Delta H)_s\ge0$ for $x\in A$.
(a) The rows are equivalent to $\lvert g^\Delta_{rs}(x;\Delta k)\rvert\le q^\Delta_{rs}(x)$ for $x\in A$ (P-S(c)), non-empty only if $q^\Delta_{rs}\ge0$ on $A$ (a data test). For $R=2$, $g^\Delta\equiv\Delta\kappa$ (constant in $x$), so the rows are $\lvert\Delta\kappa\rvert\le q_{\min}:=\min\{q^\Delta(x_L),q^\Delta(x_U)\}$ and
$$\tilde D\in[I-s\,q_{\min},\ I+s\,q_{\min}].$$
(b) The sharp identified set of $\tilde D$ under FD + DiD-S is the image of $\Theta^{FD}_I\cap\{\text{rows}\}$.
(c) The closed form (a) is the sharp set **for the regression** (the regression of $\Delta T$ on the homogeneous quadratic does not restrict $\Delta k$; 5a's sense), always. It is the sharp set for the **full law** iff the closed-form set lies in $\Theta^{FD}_I$ (2×2: iff $[-q_{\min},q_{\min}]\subseteq\Theta^{FD}_I$); otherwise (b) is strictly narrower and is computed as one LP.

*Proof.*
1. (a): P-S(b) applied to $\Delta H$ gives $O^\Delta=q^\Delta+g^\Delta$ and $L^\Delta=q^\Delta-g^\Delta$; (c) of P-S gives the equivalence. [PD:P-S] For $R=2$, $K=\kappa\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$ gives $(Kx)_v-(Kx)_n=\kappa(1-x)+\kappa x=\kappa$. $q^\Delta$ is affine in $x$, so its minimum over an interval is at an end. $\tilde D=I+s\Delta\kappa$. [alg; E7(c)]
2. (b): the rows restrict the parameter only. [param; E8(a)]
3. (c): the regression's identified set for $\Delta k$ is everything [E7(a)], so with the rows it is the row set. Under the full law it is (b), which equals the row set iff the row set lies in $\Theta^{FD}_I$. [2] ∎

*Checks.* E9 [V] (a): rows on a fine grid of the hull = closed form; (b)–(c): the truth lies in the sharp set; the closed form is sharp in most designs and strictly wider in a constructed one ($\Delta Q=(.03,0,-.03)$ and every unit at $x=.5$ starting at 0: sharp set $\{-.03\}$ against $[-.03,.03]$). **Status: P; proof review pending.** This settles PD §8's R11 item "DiD-S sharp within FD": sharp for the regression, sharp for the full law under the stated condition, and computable exactly otherwise.

### 3.4 Contrasts on changes in the model layer

**Proposition E10 (LCE-layer contrast on changes).** Under FD, for prespecified compositions $x_L,x_H$ with $x_H-x_L=\lambda(e_r-e_s)$, $\lambda>0$, and $\Delta m(x)=\Delta Hx$, $\Delta\tau(x)=x^\top\Delta Hx$:
$$\Delta\tau(x_H)-\Delta\tau(x_L)=\underbrace{\lambda\,O^\Delta_{rs}(\bar x)}_{\text{mix}}+\underbrace{\lambda\,L^\Delta_{rs}(\bar x)}_{\text{rate}},\qquad\bar x=(x_L+x_H)/2.$$
So the two contrast rows at a pair of compositions are the DiD-S rows at its midpoint; at every pair of the hull they are DiD-S on the hull. They are linear rows in $\Delta k$; Prop. E9(b)–(c) applies.

*Proof.* $\Delta\tau(x_H)-\Delta\tau(x_L)=(x_H-x_L)^\top\frac{\Delta m(x_H)+\Delta m(x_L)}2+\bar x^\top(\Delta m(x_H)-\Delta m(x_L))$ (Kitagawa). By linearity of $\Delta m$ the first term is $(x_H-x_L)^\top\Delta H\bar x=\lambda O^\Delta_{rs}(\bar x)$ and the second is $\bar x^\top\Delta H(x_H-x_L)=\lambda L^\Delta_{rs}(\bar x)$. [alg; Def] ∎

*Check.* E10 (200 random R = 3 cases, error < $10^{-12}$). **Status: P; proof review pending.** This is the changes version of NS N12c (R22(b)). It concerns the mean function at compositions; bin aggregates of model means are evaluated at row-specific weighted compositions and are not covered by this identity.

### 3.5 Levels models with model signs on changes

**Proposition E11 (per-period LCE ∩ model rows on changes).** Under per-period pooled LCE (M0 in each election; sharp set $\prod_t\Theta_{I,t}$ by E3(d)), let $\mathcal C$ be a closed convex restriction on $(k_t)_t$, for example a model sign on changes $\pm\delta^{*\Delta}_W\ge0$ or $\pm\delta^{*\Delta}_B(r,s)\ge0$ (affine in $(k_u,k_t)$ by 2c M5(a) and subtraction). Then the sharp identified set is $\big(\prod_t\Theta_{I,t}\big)\cap\mathcal C$, and the model-implied DiD's sharp set is its image. In 2×2,
$$\delta^{*\Delta}_W=a\big((Q_{vv}-Q_{vn})_t-(Q_{vv}-Q_{vn})_u\big)+b\big((Q_{vn}-Q_{nn})_t-(Q_{vn}-Q_{nn})_u\big)-(a+b)\Delta\kappa,$$
with 2c M5(c)'s $a,b$ (the same in both elections under (FP)). CR on changes is a union of two such sets.

*Proof.* [param] with E3(d). The 2×2 formula: $\delta^*_W=ab_1+bw_1$ (2c M5(c), M0), $b_1=Q_{vv}-Q_{vn}-\kappa$, $w_1=Q_{vn}-\kappa-Q_{nn}$ (2d FD2 proof step 1); subtract the two elections. ∎

*Check.* E11 (100 designs). **Status: P; proof review pending** (pooled rungs). Under district rungs the corresponding statement is R16 (open), as for M-SIGN. M1 with elections as groups (2d §5) is FD + R1 + R2 plus level linearity; it remains an alternative base to FD (alt), not an intersection.

---

## 4. How the rungs combine

### 4.1 The rule

Within a layer, rows combine in one LP (L). Across layers (count rows on the realised tables; model rows on $k$) they never share an LP; their answers for the **actual** target are intersected (∩), per 6a R27. The causal layer (§5) sits on top (↑).

**Lemma E13 (residual allowance for multi-election targets).** Let $\varphi(n)=\varphi_0+\sum_{i,r,c,t}v_{irct}n_{irct}$ with coefficients that are functions of the design, and $E_\varphi=\varphi(n^\circ)-\tilde\varphi$. Assume P3-S1a for the stacked unit vectors $(\beta_{it})_{t\in\mathcal T}$: independence **across units** given the design, with any dependence **across elections within a unit**. Then 6a Lemma 4.1 (a)–(b) hold with
$$\rho_i:=\sum_ra_{ir}\sum_t\big(\max_cv_{irct}-\min_cv_{irct}\big).$$
For the DiD $D_c(r,s;u,t)$, $\rho_i=2(w_{ir}+w_{is})$ with $w_{ir}=a_{ir}/A_r$, so $S_D=\big(\sum_i(w_{ir}+w_{is})^2\big)^{1/2}$, twice the single-election gap allowance; the bound is attained (opposite $\pm1$ changes of the two groups). The Hoeffding radius is $t_E=(\tfrac12\log(2/\alpha_E)\sum_i\rho_i^2)^{1/2}$.

*Proof.* 6a Lemma 4.1 steps 1–5 with $Z_i=\sum_ra_{ir}\sum_t\sum_cv_{irct}\beta_{irct}$: for each $(r,t)$, $\beta_{ir\cdot t}$ lies in the simplex, so $\sum_cv_{irct}\beta_{irct}\in[\min_cv_{irct},\max_cv_{irct}]$; summing over $t$ and $r$, $Z_i$ lies in an interval of length $\rho_i$, whatever the dependence across elections and groups within the unit. The $Z_i$ are independent across units. [6a:L4.1; P3-S1a; simplex] For the DiD the four non-zero coefficients per unit are $\pm1/A_r$, $\pm1/A_s$. ∎

*Check.* E12 (random dependent laws: ratio to the bound at most 0.40; the extreme law attains 1). **Status: P; proof review pending.** The Bhatia–Davis version proposed in 2d §11 and 1b B5 is **not** covered (6a Lemma 4.1: Popoviciu only for multi-cell targets; P3.5(c1) is for single cells); see §11.

**Theorem E12 (cross-layer rule on the T-election product).** Let $\mathcal C_\varphi$ be the count-layer set for an actual target $\varphi$ on the T-election tables (Lemma E1: margins and any count rows of §2), and $\mathrm{CI}_\varphi$ an actual-target model interval at level $1-\alpha$ built from any model-layer rung of §3 (FD with its rows, DiD-S, per-period LCE with model rows on changes), using Lemma E13's allowance. Then 6a R27 holds verbatim: the union bounds $(*)$; coverage at least the model interval's with **identity rows only** (margins, published totals, audited cells, in any election); with **belief rows** (realised signs, CR or contrasts on changes, orderings, boxes) coverage at least $1-\alpha-P(n^\circ\notin F)$, and only $P(\text{reject while the rows hold})\le\alpha$; the intersection is not in general the smallest valid set.

*Proof.* 6a R27's proof uses only (1) that on $\{n^\circ\in F\}$ the realised target lies in $\mathcal C_\varphi$ (R19(a) item 2), (2) an actual-target interval, and (3) the union bound. (1) holds by Lemma E1; (2) by P3.5 with Lemma E13's allowance; (3) is unchanged. [E1; E13; 6a:R27] ∎

*Check.* E13 (FD + R1, 400 units, 300 replications: coverage of the intersection equals the model interval's in every replication; disjointness rate 0.000 ≤ α). **Status: P; proof review pending**, with the model interval's status inherited (6a §4.3; §6 below). The model interval must be for the actual target: a plug-in model set is not (6a check X3).

This resolves PD §7's open items "FD4 ∩ Δ4 (mixes mean-level and realised quantities)" and "LCE-levels ∩ sign-on-changes": the first is ∩ by Theorem E12; the second is L within the model layer (E11) when the sign is a model sign, and ∩ when it is a realised sign.

### 4.2 The DiD block of the combination matrix

Codes as in NS Part 3A (L, ∩, ↑, alt, R, ×; status P, O, E).

| | Δ0 (margins) | realised signs / CR on changes (Δ4) | count contrasts on changes | FD (+ R1, R2, R4, δ-dial) | DiD-S / LCE contrasts on changes | per-period LCE + model signs on changes |
|---|---|---|---|---|---|---|
| **Δ0** | L·a·E | L·a·E (E5) | L·a·E (E6) | ∩·a (E12) | ∩·a | ∩·a |
| **realised Δ4** | | L·a·E; CR a union (E5(d),(f)) | L·a·E | ∩·a (belief rows: R27 second bullet) | ∩·a | ∩·a; **×** if called the same belief as a model sign on changes (two beliefs; NS 3A fn. 2) |
| **count contrasts** | | | L·a·E | ∩·a | ∩·a | ∩·a |
| **FD** | | | | L·m·P (E8 sharp; E23 inference) | L·m: sharp per E9(b), closed form sharp for the regression | **R**: per-period LCE ⊆ FD (2d FD4), so imposing both is per-period LCE; **alt** if offered as a robustness pair |
| **DiD-S / LCE contrasts** | | | | | L·m (E9, E10) | L·m·P (pooled; E11); O grouped (R16) |
| **per-period LCE** | | | | | | L·m·P pooled (E11); O grouped |

Other entries:
* **Causal layer**: ↑ on every cell; joint LP on the count layer (E16), joint region on the model layer (E25).
* **Estimated margins** (gap register G7): L·a with coverage $1-\alpha$ (R25(a)), one margin set for a composition shared by both elections.
* **×** FD with a changing population (FD needs the same people: G1–G3).
* **×** PD §7's "M-ORD-Δ (MTR)" as labelled: the ordering of changes is a descriptive row, not MTR (§11 item 1).

---

## 5. The causal layer

### 5.1 Potential outcomes, the estimand and the identity

* A **treated group** $r^*$ (a row; for example veterans) receives a treatment between elections 0 and 1 (service). A **comparison group** $s^*$ is another row, or all other rows ("rest").
* [PO] $Y_{jt}(0)$ is person $j$'s outcome at election $t$ had **nobody** been treated (the global no-treatment counterfactual). $B^{(0)}_{rct}=\frac1{A_r}\sum_{j\in r}\mathbf 1\{Y_{jt}(0)=c\}$.
* **(NA) No anticipation:** $Y_{jt}=Y_{jt}(0)$ for every person and every $t\le0$.
* **Estimand:** $\mathrm{ATT}_t:=B_{r^*ct}-B^{(0)}_{r^*ct}$ for $t\ge1$: the total effect of the treatment programme on the treated group's rate (it includes any effect on treated people of others' treatment).
* **Spillover onto the comparison group:** $S_t:=B_{s^*ct}-B^{(0)}_{s^*ct}$ ($S_t=0$ for $t\le0$ by NA).
* **Counterfactual trend gap:** $b_t:=(B^{(0)}_{r^*ct}-B^{(0)}_{r^*c0})-(B^{(0)}_{s^*ct}-B^{(0)}_{s^*c0})$, so $b_0=0$.
* **Parallel-trends imputation:** $\theta_t:=B^{(0)}_{r^*ct}$ and $\theta^{PT}_t(n):=B_{r^*c0}+B_{s^*ct}-B_{s^*c0}$.

**Lemma E14 (the identity).** Under (FP) and NA, with $D_t=D_c(r^*,s^*;0,t)$:
(a) $\mathrm{ATT}_t=D_t-b_t+S_t$ for $t\ge1$;
(b) $D_t=b_t$ for $t<0$: the pre-period DiDs **are** the counterfactual trend gaps, but they are set-valued, as functions of the latent tables;
(c) $\theta_t=\theta^{PT}_t+b_t-S_t$. The only logical constraints are $\theta_t\in[0,1]$ and $B_{s^*ct}-S_t\in[0,1]$.

*Proof.* NA gives $B_{r^*c0}=B^{(0)}_{r^*c0}$ and $B_{s^*c0}=B^{(0)}_{s^*c0}$. Then $D_t=(B_{r^*ct}-B^{(0)}_{r^*c0})-(B_{s^*ct}-B^{(0)}_{s^*c0})=\mathrm{ATT}_t+(B^{(0)}_{r^*ct}-B^{(0)}_{r^*c0})-(S_t+B^{(0)}_{s^*ct}-B^{(0)}_{s^*c0})=\mathrm{ATT}_t+b_t-S_t$. For $t<0$ also $B_{r^*ct}=B^{(0)}_{r^*ct}$ and $S_t=0$. (c) rearranges (a) with $\mathrm{ATT}_t=B_{r^*ct}-\theta_t$. The counterfactual rates are rates, and nothing else links them to the data under (FP). [PO; NA; alg] ∎

*Check.* E14. **Status: P; proof review pending.**

### 5.2 The labelling rule

**Definition 5.1 (causal label).** A statement about $\mathrm{ATT}_t$ may be labelled causal only if it states:
1. the treated group $r^*$ and the comparison group $s^*$;
2. a trend assumption: a dial $\lvert b_t\rvert\le M$; relative magnitudes $\bar M$ with the pre-periods used; smoothness $M$; or MTR/ND/PS;
3. a spillover assumption $\lvert S_t\rvert\le\Lambda$, with $\Lambda=0$ stated explicitly if assumed;
4. no anticipation;
5. the descriptive information it rests on (the count rows and/or model rung).

Anything else is descriptive and is printed separately.

**Proposition E15 (without a trend assumption nothing causal follows).** Under (FP), NA and any count information $F$, with no restriction on $(\theta_t,S_t)$ beyond the logical ones, the set of $\mathrm{ATT}_t$ is
$$\{B_{r^*ct}(n)-\theta:\ n\in F,\ \theta\in[0,1]\}=[B^L_{r^*ct}-1,\ B^U_{r^*ct}]$$
(ends over $F$). It always contains 0, does not depend on $s^*$, and does not use the pre-periods. Without a stated comparison group, $D_t$ and every trend assumption are undefined.

*Proof.* $\theta$ is free in $[0,1]$ [E14(c)] and $S_t$ does not enter $\mathrm{ATT}_t=B_{r^*ct}-\theta_t$ once $\theta_t$ is free. The set is $\{B-\theta\}$ over a product of an interval and $[0,1]$ [conn]. $B^U\ge0\ge B^L-1$. [alg] ∎

*Check.* E15 (LP = closed form, 30 designs). **Status: P; proof review pending.**

### 5.3 The joint causal LP

**Theorem E16 (validity, sharpness, refutability; R11b).** Let the causal assumption be a finite union $\Gamma=\bigcup_k\Gamma_k$ of polyhedra in $(n,\theta_t,S_t)$ (all forms of §§5.4–5.7 are of this kind), with $\lvert S_t\rvert\le\Lambda$. Put
$$\mathcal A=\big\{B_{r^*ct}(n)-\theta:\ n\in F,\ (n,\theta,S)\in\Gamma,\ \theta\in[0,1],\ B_{s^*ct}(n)-S\in[0,1],\ \lvert S\rvert\le\Lambda\big\}.$$
(a) **Valid:** if (FP), NA, the count information ($n^\circ\in F$) and the causal assumption hold, $\mathrm{ATT}_t\in\mathcal A$.
(b) **Sharp** for continuous tables: every value in $\mathcal A$ is the ATT of some realised table and potential-outcome configuration consistent with all the stated information. With MTR in $\Gamma$, the configuration satisfies MTR person by person.
(c) $\mathcal A$ is a union of at most $K$ intervals, each the projection of an LP (one interval when $\Gamma$ is a single polyhedron).
(d) **Refutable:** $\mathcal A=\emptyset$ refutes the conjunction of (FP), NA, the count rows and the causal assumption. In particular parallel trends with no spillover ($M=0$, $\Lambda=0$) is refuted iff $\theta^{PT}_t(n)\notin[0,1]$ for every $n\in F$.

*Proof.*
1. (a): the true $(n^\circ,\theta_t,S_t)$ satisfies every constraint: $n^\circ\in F$; $\Gamma$ by assumption; the rate limits because $\theta_t$ and $B^{(0)}_{s^*ct}$ are rates. [E14; Def]
2. (b): take a feasible $(n,\theta,S)$. Use $n$ as the realised table. For the treated group at $t$, choose unit counts $z_i\in[0,a_{ir^*}]$ with $\sum_iz_i=\theta A_{r^*}$ (possible as $\theta\in[0,1]$); under MTR take $z_i=n_{ir^*ct}\,\theta/B_{r^*ct}(n)\le n_{ir^*ct}$ (possible as $\theta\le B_{r^*ct}(n)$), and give $Y_{jt}(0)=c$ to $z_i$ of the $n_{ir^*ct}$ people of $(i,r^*)$ observed with outcome $c$, so $\mathbf 1\{Y_{jt}(1)=c\}\ge\mathbf 1\{Y_{jt}(0)=c\}$ for each of them. For the comparison group choose counterfactual counts with mean $B_{s^*ct}(n)-S\in[0,1]$. At $t\le0$ set $Y(0)=Y$ (NA). The forms of §§5.4–5.7 constrain only $\theta$, $S$ and functions of $n$, all matched. [Def; alg]
3. (c): each $\Gamma_k$ with the other constraints is a polyhedron in $(n,\theta,S)$; the ATT is linear; its image is an interval [conn]; take the union over $k$. [Def]
4. (d): (a). With $M=\Lambda=0$, $\theta_t=\theta^{PT}_t(n)$ [E14(c)], which must lie in $[0,1]$. [E14] ∎

*Checks.* E16, E18, E19, E21 (LPs). **Status: P; proof review pending.** For integer counts, $\theta A_{r^*}$ is an integer (§7, remark).

### 5.4 Form 1: the trend dial

**Corollary E17 (trend dial).** With $\Gamma=\{\lvert b_t\rvert\le M\}$ ($b_t=D_t(n)+S-(B_{r^*ct}(n)-\theta)$, linear):
(a) $\mathcal A\subseteq\mathcal D_t\oplus[-(M+\Lambda),\,M+\Lambda]$, where $\mathcal D_t$ is the descriptive set of $D_t$;
(b) equality holds if, at a maximiser $n^+$ of $D_t$ over $F$, $\theta^{PT}_t(n^+)\ge M+\Lambda$ and $B_{s^*ct}(n^+)\ge\Lambda$, and at a minimiser $n^-$, $\theta^{PT}_t(n^-)\le1-M-\Lambda$ and $B_{s^*ct}(n^-)\le1-\Lambda$ (the rate limits are slack);
(c) for the conclusion "$\mathrm{ATT}_t>0$" on the Minkowski form, the breakdown value is $M^*=\min\mathcal D_t-\Lambda$ (when positive).

*Proof.* (a): $\mathrm{ATT}=D-b+S$ [E14(a)]. (b): at $n^+$ take $b=-M$, $S=\Lambda$, so $\theta=\theta^{PT}-M-\Lambda\in[0,1]$ and $B_{s^*}-\Lambda\in[0,1]$, which gives $\max\mathcal D_t+M+\Lambda$; symmetric at $n^-$. (c): the lower end is $\min\mathcal D_t-M-\Lambda$. [E14; E16] ∎

*Checks.* E16 [V] (the LP lies in the Minkowski set; equal in 40 cases, rate limits bind in 80); **E17 [I]** NS's arithmetic: $\mathcal D=[.05,.25]$, $M=.10$ gives $[-.05,.35]$, $M=.02$ gives $[.03,.27]$, $M^*=.05$. **Status: P; proof review pending.**

### 5.5 Form 2: bounds relative to pre-period trends (R11b′)

With $L\ge1$ pre-periods, write the violations as $\delta_t:=b_t$ ($\delta_0=0$; $\delta_t=D_t(n)$ for $t<0$ by E14(b)).

**Definition 5.2 (Rambachan–Roth relative magnitudes, RM($\bar M$)).** For every post-period $t\ge1$,
$$\lvert\delta_t-\delta_{t-1}\rvert\le\bar M\,A(n),\qquad A(n):=\max_{-L<s\le0}\lvert\delta_s-\delta_{s-1}\rvert=\max_{-L<s\le0}\lvert g_s(n)-g_{s-1}(n)\rvert,$$
the largest consecutive pre-period DiD in absolute value. For one post-period this is $\lvert b_1\rvert\le\bar M A(n)$. NS §3.2 writes $\max_{t\le0}\lvert D_t\rvert$ (levels relative to the base); the two coincide for $L=1$ and differ for $L\ge2$ (check E19 [I]).

**Definition 5.3 (smoothness, SD($M$)).** $\lvert(\delta_{t+1}-\delta_t)-(\delta_t-\delta_{t-1})\rvert\le M$ for every $t$ with $t\pm1\in\mathcal T$ (pre-periods included, which makes SD testable).

**Proposition E18 (the joint form).**
(a) Under RM the ATT set (for one post-period, or the vector of post-periods) is Theorem E16's $\mathcal A$ with
$$\Gamma_{RM}=\bigcup_{-L<s\le0}\ \bigcup_{\sigma=\pm1}\Big\{\lvert\delta_t-\delta_{t-1}\rvert\le\bar M\,\sigma\big(g_s(n)-g_{s-1}(n)\big)\ \ \forall t\ge1\Big\},$$
a union of $2L$ polyhedra, so a union of $2L$ LP projections. Under SD it is one polyhedron.
(b) The **endpoint-wise** set $\mathcal A_{EW}$, which combines the separate range of $D_t$ with $\bar M\sup_FA$, contains the joint set; only the joint set is sharp (E16(b)).

*Proof.* (a): $x\le\bar M\max_ky_k$ iff $x\le\bar My_k$ for some $k$; the maximum is the same for every post-period, so one $k=(s,\sigma)$ serves all $t$. Each branch is linear in $(n,\theta,S)$ ($\delta_t$ for $t\ge1$ is linear by E14, and $\lvert\cdot\rvert\le$ linear is two linear rows). [alg; E14] (b): every feasible $(n,b)$ has $D_t(n)$ in its range and $\lvert b\rvert\le\bar MA(n)\le\bar M\sup_FA$. [Def] ∎

*Checks.* E18 [V] (N21 reproduced exactly; the union of LPs equals a grid brute force). **Status: P; proof review pending.**

**Proposition E19 (degeneration; R11b′(a)–(d)).** One post-period ($t=1$), (FP), NA.
(a) **Contained in a trend dial.** $\mathcal A_{RM}\subseteq\mathcal A_{TD}(M=\bar M\sup_FA)$ (both with the rate limits). Without the rate limits, the upper end of $\mathcal A_{RM}$ is $\max_F(D_1+\bar MA)\le\max_FD_1+\bar M\sup_FA$, with equality **iff** some $n\in F$ attains both $\max_FD_1$ and $\sup_FA$ (compatible maximisers); likewise for the lower end with $\min_FD_1$.
(b) **When it is uninformative.** With the rate limits and $\Lambda=0$, the upper end of $\mathcal A_{RM}$ equals the no-assumption upper end $B^U_{r^*c1}$ (Prop. E15) iff some $n\in F$ has $B_{r^*c1}(n)=B^U_{r^*c1}$ and $\lvert\theta^{PT}_1(n)\rvert\le\bar MA(n)$; the lower end equals $B^L_{r^*c1}-1$ iff some $n$ has $B_{r^*c1}(n)=B^L_{r^*c1}$ and $\lvert1-\theta^{PT}_1(n)\rvert\le\bar MA(n)$. Each condition is the feasibility of $2L$ LPs. A sufficient condition for "no information at all" is $\bar M\inf_FA\ge\sup_F\max\{\lvert\theta^{PT}_1\rvert,\lvert1-\theta^{PT}_1\rvert\}$. A positive spillover allowance $\Lambda$ only enlarges the set.
(c) **When it can be informative.**
  * (c1) Monotone: $F'\subseteq F$ implies $\mathcal A_{RM}(F')\subseteq\mathcal A_{RM}(F)$; so every rung that narrows the pre-period DiDs narrows the RM set: model rungs (FD/R1 on the pre-pair, through ∩), Δ4 rows, contrasts on changes, known pre-period cells.
  * (c2) **Margins only** (product structure): the pre-period tables other than the base are free in their own boxes, so for a base gap $g_0$ the allowance reaches $\bar M\bar A(g_0)$, with $\bar A(g_0)=\max\{\lvert g_0-G^-_{-1}\rvert,\lvert g_0-G^+_{-1}\rvert,\ \max_{s\le-1}\max(\lvert G^+_s-G^-_{s-1}\rvert,\lvert G^-_s-G^+_{s-1}\rvert)\}\ge\tfrac12w_{-1}$, where $w_{-1}=G^+_{-1}-G^-_{-1}$ is election $-1$'s gap width. So $\mathcal A_{RM}=\bigcup_{g_0}\mathcal A_{TD}\big(\bar M\bar A(g_0)\,\big|\,g_0\big)$: a trend dial whose $M$ is at least $\bar M w_{-1}/2$, a pure artefact of the pre-period's Duncan–Davis width.
  * (c3) With **point-identified** pre-trends ($g_s$ constant on $F$ for $s\le0$), $\mathcal A_{RM}=\{D_1(n)-b:\ n\in F,\ \lvert b\rvert\le\bar MA\}$ (with the rate limits): standard Rambachan–Roth applied to the descriptive post-period set.
(d) **Sampling.** With model-based pre-trends, sampling uncertainty enters through one joint confidence region for all elections; a plug-in for the pre-period quantities under-covers. This is Prop. E25 (§6).

*Proof.*
1. (a): every feasible $(n,b)$ has $\lvert b\rvert\le\bar MA(n)\le\bar M\sup_FA$ [Def], so it is feasible for the dial. Without rate limits the ATT range at $n$ is $[D_1(n)-\bar MA(n),D_1(n)+\bar MA(n)]$, so the upper end is $\max_F(D_1+\bar MA)$; the inequality is $\max(f+g)\le\max f+\max g$, with equality iff $f$ and $g$ have a common maximiser ($F$ compact; $D_1$ and $A$ continuous). [E18(a); alg]
2. (b): at $n$, $\mathrm{ATT}=B_{r^*c1}(n)-\theta$ with $\theta=\theta^{PT}_1(n)+b$ [E14(c)] and $\lvert b\rvert\le\bar MA(n)$, $\theta\in[0,1]$. The value $B_{r^*c1}(n)$ (i.e. $\theta=0$) is allowed iff $\lvert\theta^{PT}_1(n)\rvert\le\bar MA(n)$. The no-assumption upper end is $B^U$ [E15], attained by RM iff $\theta=0$ is allowed at some $n$ with $B_{r^*c1}(n)=B^U$. Symmetric for $\theta=1$. "$\lvert x\rvert\le\bar M\max_{s,\sigma}y$" is a union of $2L$ polyhedra [E18(a)]. Under the sufficient condition, $\theta=0$ and $\theta=1$ are allowed at every $n$, so $\mathcal A_{RM}$ is the set of E15. [E15]
3. (c1): fewer tables, fewer feasible points. [Def] (c2): under margins only $F$ is a product over elections [E1]; for fixed base and post tables, the ATT range is centred at $D_1$ with radius $\bar MA(n)$ (cut by the rate limits), increasing in $A$, so the union over the free pre-period tables is the range with the largest $A$. The largest consecutive pre DiD touching the base is $\max\lvert g_0-g_{-1}\rvert$ over $g_{-1}\in[G^-_{-1},G^+_{-1}]$, at an end of that interval, and $\max(\lvert g_0-G^-\rvert,\lvert g_0-G^+\rvert)\ge(G^+-G^-)/2$; pairs not touching the base are free at both ends. [E3(b); sep; alg] (c3): $A(n)$ is constant on $F$. [Def]
4. (d): Prop. E25. ∎

*Checks.* E19 [V] (a) (40 designs: equal in 20, strictly inside in 20, and the "iff" agrees with an LP test for a common maximiser); (b) (the vacuity LP predicts the upper end in every design: vacuous 37, informative 2, refuted 1); (c) (margins only: the allowance formula given $g_0$, and point-identified pre-trends give standard RR); E24 for (d). **Status: P; proof review pending.**

**N21 (NS §3.2) reproduced, and corrected.** Two precincts of 100 men with 20 and 60 veterans; yes votes $(30,45)$, $(32,48)$, $(30,60)$ at elections $-1,0,1$; $\bar M=1$.
* By LP (check E18): the pre-period DiD is in $[-1.229,1.271]$ and the post-period DiD in $[-1.083,1.417]$. Without rate limits the RM set is $[-2.354,2.646]$ jointly and $[-2.354,2.688]$ endpoint-wise: NS's numbers, now exact.
* **Correction.** The ATT of a rate lies in $[B^L_{v,1}-1,B^U_{v,1}]=[-0.75,1.00]$ (Prop. E15), not in the DiD's logical range $[-2,2]$. With the rate limits the joint RM set is exactly $[-0.75,1.00]$: the no-assumption set. NS's conclusion ("no information") stands, but the comparison must be with the ATT's own logical range. By (b): at the tables that put $B_{v,1}$ at its ends, the pre-period DiD can be large enough ($\lvert\theta^{PT}\rvert\le A$ and $\lvert1-\theta^{PT}\rvert\le A$).
* With the pre-period totals pinned to $\pm2$ votes (NS's illustration), the pre-period DiD is in $[-0.06,0.10]$ and the RM set is $[-0.60,0.90]$ without rate limits, $[-0.604,0.896]$ with them (check E18 [I]).

**Proposition E20 (smoothness with set-valued pre-trends).** One pre-pair ($L=1$) and one post-period. SD($M$) reads $\lvert b_1+D_{-1}\rvert\le M$ with $D_{-1}=g_{-1}-g_0$, so $\mathrm{ATT}_1\in\{g_1-2g_0+g_{-1}+S\pm M\}$. Under margins only, without rate limits,
$$\mathcal A_{SD}=\big[G^-_1-2G^+_0+G^-_{-1}-M-\Lambda,\ \ G^+_1-2G^-_0+G^+_{-1}+M+\Lambda\big],$$
whose width exceeds the parallel-trends width $w_1+w_0$ by $w_0+w_{-1}+2M+2\Lambda$: the base election's Duncan–Davis width counts twice. With $L\ge2$, SD also restricts the pre-period tables (second differences of the pre-period gaps at most $M$), a testable restriction.

*Proof.* $\delta_1=b_1$, $\delta_0=0$, $\delta_{-1}=D_{-1}$ [E14(b)], so the constraint is $\lvert b_1+D_{-1}\rvert\le M$, and $\mathrm{ATT}=D_1-b_1+S=g_1-g_0+D_{-1}+S\mp M$. The coefficients of $(g_1,g_0,g_{-1})$ are $(1,-2,1)$; the elections are independent factors with interval ranges [E3(a),(b); sep]. [alg] ∎

*Check.* E20 (closed form = LP; width identity; 30 designs). **Status: P; proof review pending.** In ecological data smoothness is therefore **more** fragile than relative magnitudes: linear extrapolation of a set-valued pre-trend doubles the base period's uncertainty.

### 5.6 Form 3: MTR, ND and PS (2d Prop. CL generalised)

At a post-period $t$, for the indicator of outcome $c$:
* **(MTR)** monotone treatment response for the treated, person by person: $\mathbf 1\{Y_{jt}(1)=c\}\ge\mathbf 1\{Y_{jt}(0)=c\}$; it implies $\theta_t\le B_{r^*ct}$;
* **(ND)** no decline: $\theta_t\ge B_{r^*c0}$;
* **(PS)** positive selection on the **untreated** comparison rate: $\theta_t\ge B^{(0)}_{s^*ct}=B_{s^*ct}-S_t$.

**Proposition E21 (MTR + ND + PS, R×C, any T, with spillovers).**
(a) The set is non-empty iff some $n\in F$ and $\lvert S\rvert\le\Lambda$ (with $B_{s^*ct}(n)-S\in[0,1]$) satisfy $B_{r^*ct}(n)\ge\max\{B_{r^*c0}(n),B_{s^*ct}(n)-S\}$. So the three assumptions are jointly refutable.
(b) If non-empty, the set is $[0,U]$ with $U=\max_{n\in F,\,\lvert S\rvert\le\Lambda}\min\{B_{r^*ct}-B_{r^*c0},\ B_{r^*ct}-B_{s^*ct}+S\}$, an LP in epigraph form.
(c) Margins only: $U=B^U_{r^*ct}-\max\{B^L_{r^*c0},\ \underline B_{s^*ct}-\Lambda\}$, where $\underline B_{s^*ct}=B^L_{s^*ct}$ for a named row (the maximum of $B_{r^*ct}$ and the minimum of $B_{s^*ct}$ are attained by one table, 6a Lemma 2.1(b)) and $\underline B_{s^*ct}=(Y_{ct}-Y^U_{r^*ct})/(N-A_{r^*})$ for the rest.
(d) With further rows the ends in (c) need not be jointly attainable; then use (b).

*Proof.*
1. (a), (b): the ATT set at $(n,S)$ is $\{B_{r^*ct}(n)-\theta:\max\{B_{r^*c0},B_{s^*ct}-S\}\le\theta\le B_{r^*ct}\}=[0,\,B_{r^*ct}-\max\{\cdot\}]$, non-empty iff the stated inequality holds; all contain 0, so the union is $[0,U]$; the minimum of affine functions is concave, so $U$ is an LP in epigraph form. [E16; 2d:CL steps 1–3]
2. (c): the base table is a free factor, so $B_{r^*c0}$ is at its minimum; at $t$ both terms of the minimum increase in $B_{r^*ct}$ and the second decreases in $B_{s^*ct}$; the maximum of the first and the minimum of the second are attained by one table (6a Lemma 2.1(b) for a named row; for the rest $B_{\neg r^*,ct}$ is a decreasing function of $Y_{r^*ct}$); $S=\Lambda$ is optimal. [E3; 6a:L2.1] ∎

*Checks.* E21 [V] (closed form = LP for the rest and for a named row, R = C = 3, with $\Lambda\in\{0,.03\}$); [P2d] Δ7 (2×2). **Status: P; proof review pending.** 2d's PS uses the **observed** comparison rate; that is (c) with $\Lambda=0$, i.e. it assumes no spillover (§11 item 2). These three are bounds on $\theta$, hence on $b_t$: they are trend assumptions in a different parameterisation, and Definition 5.1 counts them as such.

### 5.7 Spillovers onto the comparison group

**Proposition E22 (spillovers need their own dial).**
(a) **Pre-periods carry no information on $S$.** For any $S,S'$ in $[-\Lambda,\Lambda]$ there are configurations with identical realised tables at every election (hence identical data and identical pre-period DiDs), the same $b_t$, and ATTs that differ by $S'-S$ (provided the rate limits are slack).
(b) Consequently RM and SD calibrate $b_t$ only. The ATT set must allow $\lvert S_t\rvert\le\Lambda$ separately (Theorem E16 with $S$ free in $[-\Lambda,\Lambda]$); with $\Lambda=0$ when $S\ne0$ the RM set can exclude the true ATT.
(c) **What an ecological model can see.** If the spillover onto a comparison person in unit $i$ is $\sigma_0+\sigma_1X_{ir^*}$, then $S_t=\sigma_0+\sigma_1\sum_iw_{is^*}X_{ir^*}$. The composition-linked part $\sigma_1$ moves the comparison group's change with the local treated share, so it is part of the comparison group's contextual change (2×2: $\Delta w_1$), which R1 sets to zero and FD3's δ-dial bounds; it is confounded there with untreated contextual trends. The uniform part $\sigma_0$ shifts every unit's comparison change equally and is invisible to every comparison across units. R1 therefore does not rule out spillovers.

*Proof.* (a): fix $n$ and $b$; set $\theta=\theta^{PT}+b-S$ and $B^{(0)}_{s^*}=B_{s^*}-S$, and the same with $S'$ [E14(c)]; NA leaves all $t\le0$ unchanged; $\mathrm{ATT}=D-b+S$ [E14(a)]. (b): (a) with Theorem E16. (c): average the unit spillovers with the comparison group's weights; a term constant across units adds the same amount to every unit's comparison change. [Def; alg] ∎

*Checks.* E22 [V] (a); **E22 [F]** (b): with pure units (every rate identified) and a true spillover of .08, the RM set with $\Lambda=0$ excludes the true ATT and the set with $\Lambda=.08$ contains it. **Status: P; proof review pending.**

---

## 6. Inference (R11c)

**Count-layer rungs** (Δ0, Δ4, contrasts on changes, and the causal LP on top of them) are exact given the margins and rows: no sampling inference. They come with slack dials and breakdown values (2c M2(d); 6a R19(a) item 4).

**Proposition E23 (FD and its rows).** Assume (FP), P3-S1a for the stacked unit vectors (Lemma E13), P3-S2, FD, A5, and a full-rank design for the regression of $\Delta T_c$ on the homogeneous quadratic ($C-1$ free columns).
(a) $\widehat{\Delta Q}$ is a fixed-dimensional OLS estimator with the usual asymptotics [OLS].
(b) FD with any finite set of rows on $\Delta k$ whose coefficients are linear in $\Delta Q$ and the design (R1, R2, R4 as equalities; $\lvert\Delta w_1\rvert\le\delta$; DiD-S rows at prespecified compositions; LCE contrasts on changes; model signs on changes) is a constraint system of P3's form $s(\beta)+Ak\ge0$ with $\beta=\Delta Q$ and $k=\Delta k$. P3.5(a) and (b2)(ii) give model-implied intervals, and with Lemma E13's allowance actual-target intervals: **P** (asymptotic), in any regime for the pooled FD rung (fixed dimension). The CLR version is algorithmic (P3.5(b3)). P3.8 (finite sample) applies with Lemma E13's Hoeffding radius.
(c) The FD4 outer set is also of P3's form, with $k=(\Delta k,(m_u(x_j))_j)$ and $\beta$ including the coefficients of $\tau_u,\tau_t$, **if** $\tau_t(x)=E[T_t\mid X=x]$ is a polynomial of a prespecified degree; otherwise it is a heuristic, or an outer relaxation with person-weighted bins (2c Lemma N2.4). Inference for the **sharp** FD set (E8) needs conditional expectations of the kinked functions $L(T),U(T)$ given $X$: **O**.

*Proof.* (a) [OLS] under P3-S1a (independent units, bounded outcomes). (b): each row is $c(\Delta Q)+a^\top\Delta k\ge0$ with $c$ linear; P3.5(a) and (b2)(ii) use only this form, a fixed-dimensional asymptotically normal $\hat\beta$ and a consistent variance; the allowance step uses only Lemma E13's range bound. [P3.5; E13] (c): the FD4 constraints $\sum_rx_{jr}m_{u,rc}(x_j)=\tau_{uc}(x_j)$, $0\le m_u(x_j)\le1$, $0\le m_u(x_j)+M_{\Delta H}(x_j)\le1$ are linear in $(m_u,\Delta k)$ with right-hand sides linear in $(\tau$-coefficients$,\Delta Q)$. [2d:FD4; Def] ∎

*Check.* E13 (FD + R1: coverage of the actual-target interval and of its intersection with Δ0). **Status: P; proof review pending** for (a), (b) and for (c) under its polynomial condition; the sharp-set inference is **O**.

**Proposition E24 (unions across branches and elections).** Let a model be a union of branches $k\in\mathcal K$ (sign patterns across period pairs; CR branches; RM's arg-max branches), and suppose that whenever branch $k$ holds, $P(\varphi\in\mathrm{CI}_k)\ge1-\alpha$ (asymptotically). Then $\mathrm{CI}_\cup:=\bigcup_k\mathrm{CI}_k$ has $P(\varphi\in\mathrm{CI}_\cup)\ge1-\alpha$ under the union model, with no multiplicity correction, and "$\mathrm{CI}_\cup=\emptyset$" is a level-α test of the union model. Reporting $\mathrm{CI}_{\hat k}$ for a data-selected branch $\hat k$ is not covered.

*Proof.* Some branch $k^\circ$ holds; $\{\varphi\in\mathrm{CI}_{k^\circ}\}\subseteq\{\varphi\in\mathrm{CI}_\cup\}$, and $\{\mathrm{CI}_\cup=\emptyset\}\subseteq\{\varphi\notin\mathrm{CI}_{k^\circ}\}$. [union; Def] ∎

*Check.* **E23 [F]** (true branch $\theta\ge0$, $\theta=.5$: the union covers at .949 ≈ .95; the plug-in branch's interval covers at .669). **Status: P; proof review pending.** This is PD §8 R10's "unions" item for DiD.

**Proposition E25 (one joint region for model-based pre-trends; E19(d)).** Suppose the targets entering a causal form ($D_t$ for pre- and post-periods, $B_{r^*ct}$, $\theta^{PT}_t$) are model-implied linear functionals of one stacked, fixed-dimensional parameter $\beta$ (the regressions of all elections, or of all period-pair changes) and of LP variables $k$, as in P3. Let the projection event of P3.5(b2)(ii) be taken **once** for the stacked $\beta$ ($q=\dim\beta$).
(a) Computing each branch of E18(a) (or SD, or the dial) with P3.5(b2)(ii) on that one event, and taking the union over branches, gives a set that covers the model-implied ATT with probability $\ge1-\alpha$ asymptotically, simultaneously for all branches and all post-periods. For the actual ATT add the allowance events of each actual target used, with a union bound (6a Remark 4.2).
(b) Separate confidence statements per election (each at $1-\alpha$) can only be combined by the union bound, at level $1-\sum_t\alpha_t$, and only into the endpoint-wise set: they do not identify which pre- and post-period tables go together, so they cannot deliver the joint set.
(c) Plugging in point estimates of the pre-period quantities (an interval for the post-period DiD only) is not covered and under-covers.

*Proof.* (a): on the projection event all linear functionals of $\hat\beta-\beta_0$ are within $\hat k\hat\sigma$ at once [P3.5(b2)(ii) step 4], so each branch's interval contains its population value whenever the branch's constraints hold at $\beta_0$ [P3.5(b2) steps 1–2]; the RM branch that holds at the truth exists [E18(a)]; E24's union argument applies [E24]. Allowance events: [union; E13]. (b): [union]; joint pairing needs the joint law [E18(b)]. (c): check E24. ∎

*Check.* **E24 [F]** (three Gaussian gap estimates, $\bar M=1$, ATT at the RM lower end: joint region .998, plug-in .873 at nominal .95). **Status: P; proof review pending** for (a)–(b); (c) by the computed counterexample.

---

## 7. Integer versions (R11d, with R15)

**Proposition E26 (integer counts: margins and one-vs-rest rows on changes).** (FP), integer margins, $\operatorname{Var}_p(X_r)>0$.
(a) Margins only: the integer values of $D^r_c(u,t;\mathcal S)$ are $D^r_c(Z)$ for all integers $Z\in[Z^-,Z^+]$ (E5): a lattice of step $\kappa_r=1/A_r(\mathcal S)+1/(N(\mathcal S)-A_r(\mathcal S))$, every point attained, both ends equal to the continuous ends. For a named pair with $R=2$ this is the same lattice.
(b) With W1 sign rows on changes for the same $(r,c,u,t,\mathcal S)$ (each a half-line in $Z$), the integer set is the lattice intersected with the continuous set of E5(b); its ends differ from the continuous ends by less than one step; integer emptiness occurs only when the continuous set is shorter than one step.
(c) CR on changes: the integer set is the lattice points between $NM^r_\Delta$ and $ER^r_\Delta$ (within $[D^-,D^+]$); it can be empty when the continuous set is not.

*Proof.* (a): each election's $Y_{rct}(\mathcal S)$ takes every integer between its integer Duncan–Davis extremes (a sum of unit integer intervals) [6a:L2.1(a); 2c Remark M2.1 proposition], independently across elections [E1], so $Z$ takes every integer in $[Z^-,Z^+]$; $D^r_c$ is affine in $Z$ with slope $\kappa_r$ [E4(b)]. (b): the rows are half-lines in $Z$ [E4(b),(c)]. (c): the union of the two branches of (b). ∎

*Checks.* E3 [V] (brute force over all integer two-election 2×2 tables in 59 designs: margins, each W1 sign, CR — integer CR empty with a non-empty continuous set in one design — and the contrast of E6). **Status: P; proof review pending.**

**Theorem E27 (total unimodularity for change rows; R = 2).** Integer margins, $R=2$. Rows: boxes on single cells (any election), and **change rows** $\sum_{i\in U}(n_{i1ct}-n_{i1cu})\in[h^-,h^+]$ on **one** column $c$, such that (i) the period pairs $(u,t)$ used, seen as edges on the set of elections, form a graph of maximum degree 2 without cycles (a union of paths: one pair; consecutive pairs; two pairs sharing a base), and (ii) the unit sets $U$ of each pair form a laminar family. After rounding the rows (6a §1.4), $F_{\rm rd}$ is integral: the integer ends of every linear target are those of the rounded LP, attained by integer tables. (Rows on row 2 are rows on row 1 via $n_{i2ct}=b_{ict}-n_{i1ct}$.)

*Proof.*
1. Substitute $y_{itc'}=n_{i1c't}$, so $n_{i2c't}=b_{ic't}-y_{itc'}$ (an integral affine bijection). The reduced system has block rows $\sum_{c'}y_{itc'}=a_{i1}$ (each variable in exactly one block), bounds, and the change rows, which involve only column-$c$ variables $y_{itc}$. [6a:R19(b) step 5]
2. Ghouila-Houri: take any subset of rows. **Change rows.** Give each pair $\pi$ a sign $\varepsilon_\pi=\pm1$ and sign the chosen sets of $\pi$ by depth [GH-depth], times $\varepsilon_\pi$. Since a row of $\pi=(u,t)$ has $+1$ on $y_{itc}$ and $-1$ on $y_{iuc}$ for $i\in U$, pair $\pi$ contributes $\varepsilon_\pi\{0,1\}$ to each later-election variable and $-\varepsilon_\pi\{0,1\}$ to each earlier-election variable. Choose the $\varepsilon_\pi$ along each path so that at every election lying in two pairs the two contributions have opposite signs (each new pair's sign is fixed by its neighbour's; no cycles, so no conflict). Then every column-$c$ variable has partial sum $\sigma\in\{-1,0,1\}$. [GH-depth; alg]
3. **Block rows.** For a chosen block $(i,t)$ give it the sign $-\sigma_{itc}$ if $\sigma_{itc}\ne0$, else $+1$. Then $y_{itc}$ sums to 0 or $\pm1$, and every other variable of the block, which is in no change row, sums to $\pm1$ or 0. [alg]
4. So every subset of rows has a GH split; the reduced matrix is TU [GH]; with integer right-hand sides after rounding, every vertex is integral [HK]; the affine bijection maps integer vertices to integer vertices. [6a:L1.1] ∎

*Checks.* E25 [V] (LP = MILP for 540 random objectives: one pair; a consecutive path of three pairs; two pairs sharing a base; each with its own laminar family). Three pairs sharing a base (degree 3) are outside the theorem; no integrality gap was found in 180 objectives, so a wider class is conjectured (**C**). Several columns with change rows, or change rows mixed with within-election totals of 6a's classes, are outside the theorem: LP outer, MILP exact. **Status: P; proof review pending.** The W1 one-vs-rest signs on changes over one unit set, for consecutive pairs, are change rows on one column with $U=$ all units, so they are covered.

**Remark 7.1 (causal integer versions).** With integer counts the counterfactual counts are integers too, so $\theta_tA_{r^*}\in\mathbb Z$ and the ATT lies on a lattice of step $1/A_{r^*}$ inside Theorem E16's set. A general statement of how far the integer ends can lie inside the continuous ones (beyond one step for the dial with margins only) is **O**; the certified MILP of 6a C4 gives exact answers.

---

## 8. The gap register (changing populations, uncertain margins, covariates)

*By the owner's decision these items are **registered, not solved**. Changes in composition and uncertain margins are deferred to the general "uncertain totals and compositions" module at the end of the plan (`NEXT_BUILD.md`); this register is its input. For each item: what (FP) or exact margins assumed; what still holds without the item; what assumption or data-based range would close it; which rungs it modifies; where a set with confidence $1-\alpha$ would enter; and a check that the naive statement fails when the assumption is dropped.*

### 8.1 Summary table

| # | Item | Holds without it (theorem) | Closed by (assumption, or data-based range) | Rungs modified | Where a $1-\alpha$ set enters | Check |
|---|---|---|---|---|---|---|
| G1 | **Entry** (in-migrants, new arrivals) | All count-layer results for the **period-own** populations (Δ0 of period-own rates; each election separately exact, 6a R19). Within-person targets, Δ4 (individual changes) and FD are **undefined** for entrants | entrant rows with known or interval entrant counts by group (accounting polytope, 2d C1); or "entrants vote like stayers of their group" ($B^E=B^O$), or a dial $\lvert B^E_r-B^O_r\rvert\le\eta$ | Δ0 (target), Δ4 and FD (stayers only), causal (ATT of stayers vs of the period population) | entrant counts by group: a set $\mathcal E$ with $P(E\in\mathcal E)\ge1-\alpha$ | **G1 [F]**: period-own Δ0 misses the stayers' DiD in 200/200 designs; the entrant-row polytope covers it |
| G2 | **Exit** (death, out-migration) | Exact for the **stayers** if their period-$u$ table is identified separately; otherwise the period-$u$ pool mixes leavers | free exit cells (model-free); "leavers voted like stayers" ($B^E=B^O$); dial $\lvert B^E_{r}-B^O_{r}\rvert\le\eta$, which moves the stayers' rate by at most $(1-\pi_r)\eta$ (2d C3) | Δ0, Δ4, FD (base period), causal | leaver counts by group (deaths partly known: a set) | **G2 [F]**: pool rate ≠ stayers' rate (299/300); the dial identity exact |
| G3 | **Aging-in** (became eligible between elections) | On a common later base, the DiD is the stayers' DiD plus $(1-\omega_v)B^O_{v,u}-(1-\omega_n)B^O_{n,u}$ when aged-in vote like stayers at $t$ (2d C2), exact; a structural zero for the aged-in at $u$ | the structural zero (aged-in cast no $u$ votes) plus the aged-in split by group (interval, from age structure); or "aged-in vote like stayers" | Δ0 (base), all rungs through the target | aged-in counts by group | **G3 [F]**: common-base DiD ≠ stayers' DiD (300/300), gap as stated |
| G4 | **Stayers vs leavers** (selection on staying) | Every descriptive result for the **stayers** (a fixed population); causal results for the stayers' ATT | "parallel trends among stayers" stated directly (not inferred from the pool); or no differential selection on trends; or a dial between the pool's and the stayers' $b$ | causal layer (target and trend assumption) | — (an assumption, not a count) | **G4 [F]**: $b=0$ in the pool, $b=.06$ among stayers |
| G5 | **Interval bases from censuses** ($N_{it}$ between two estimates) | Exact for the **stated** base; an undercounted base can give a range that excludes the truth | $N_{it}\in[N^L_{it},N^U_{it}]$ with floors $N_{it}\ge$ ballots cast, $\ge$ group counts: cell counts stay a polytope (2d §7; 1b D.7) | every count-layer rung (margins become interval-valued); LCE layer ($X$ and weights uncertain) | the base interval as a confidence set (census error model) | **G5 [F]**: a census-projected base taken as exact excludes the truth in 300/300 designs; the interval base covers it |
| G6 | **Changing unit boundaries** | Every count-layer result on the **coarsest common units** (Lemma G6: valid, weakly wider) | a crosswalk with known population shares (an assumption); interval shares (a polytope); or coarsening | the unit set of every rung; LCE compositions | crosswalk shares | **G6 [F]**: share-allocated margins exclude the truth or are infeasible in 51/300; coarsening covers in all |
| G7 | **Estimated or uncertain margins and compositions** (matched counts, survey margins, non-integer counts), in either election | Exact for the **stated** margins (6a R19(e)); integer check off for non-integer margins | a margin set $\mathcal V$ with $P(a^\circ\in\mathcal V)\ge1-\alpha$ (R25(a)), **one** set for a composition shared by all elections; R26 (Berkson) in the LCE layer, which carries over to FD (Prop. G7(c)) | every rung | here, by definition | **G7 [F]**: expected counts as exact miss the truth (42/200); the union covers; **G7 [V]** one set for both elections is strictly narrower (20/20); **G7 [F]** R26 on changes |
| G8 | **Covariates** (low priority) | Unconditional results; with strata published as units, the count layer within strata is exact | conditional PT within strata, $\mathrm{ATT}=\sum_zw^{r^*}_z(D_z-b_z+S_z)$; CAR(Z)-δ on changes (NS step 6) | causal layer; FD (covariates in the mean) | — | **G8 [F]**: conditional PT, unconditional DiD .010 vs ATT .050 |
| G9 | **Time-varying group membership** (a group defined at each election) | Results for people whose group is constant (a fixed population of group stayers) | membership fixed by definition (e.g. "future veterans"); or switcher rows (an accounting polytope with switcher counts) | Δ4, FD (individual changes undefined for switchers), every target | switcher counts | **G9 [F]**: period-group DiD .450 vs within-person .503 |

### 8.2 Notes per item

**G1–G3 (entry, exit, aging-in) share one device: the accounting polytope** (2d §7, Prop. C1). Row types per election (stayers, leavers, entrants, aged-in, by group) with counts $\rho_{irt}\ge0$, linked across elections for stayers, structural zeros where a type cannot vote, and cells summing to the observed outcomes. Targets with known denominators are linear; with unknown denominators they are linear-fractional (2d C1(a)–(b): Charnes–Cooper or a scan). The stayers' DiD is the natural within-person target; the period-population DiD is a different descriptive target. Every rung of §§2–5 applies to the stayers' block **once its margins are known or interval-valued**; without that the stayers' margins are unknown and only the polytope bound holds.

**G4.** Selection on staying is a causal issue even when the stayers' tables are known: "parallel trends in the original population" does not imply "parallel trends among stayers" when staying differs by group and is related to trends. The trend assumption must be stated for the population whose ATT is reported.

**G5.** The interval base turns each unit's margins into $N_{it}\in[N^L_{it},N^U_{it}]$, $a_{irt}$ and $b_{ict}$ with $\sum_ra_{irt}=\sum_cb_{ict}=N_{it}$, still a polytope in cell counts with $N$ as a variable. Under (FP) the base is the same people at every election, so one interval serves all elections (as in G7).

**Lemma G6 (coarsening is valid and loses information).** Merge units into blocks that are stable across elections. Every latent table of the fine units aggregates to one of the blocks with the blocks' margins, so the set of values of any target that depends on the table only through block-level counts (group rates over the electorate or over unions of blocks) computed on the blocks contains the set computed on the fine units: valid, and weakly wider.
*Proof.* Summing a fine table over the units of a block gives a nonnegative table with the block's margins; rows stated on the blocks hold for the aggregate; the target is unchanged. So the fine feasible set maps into the coarse one. [Def; lin] ∎ **Status: P; proof review pending** (check G6).

**Proposition G7 (where a margin set enters, and what carries over).** Let the group counts be uncertain, with a set $\mathcal V$ such that $P(a^\circ\in\mathcal V)\ge1-\alpha$.
(a) Under (FP) the same $a^\circ$ enters every election. The set $\mathcal C:=\bigcup_{a\in\mathcal V}\varphi(F(a))$ covers the actual target with probability $\ge1-\alpha$ (R25(a) on the T-election product); each $F(a)$ is the LP of Lemma E1. Using a separate set per election, $\bigcup_{(a_t)\in\mathcal V^T}$, is valid but wider, and if the sets are separately random it spends $\alpha$ once per election.
(b) The causal LP (E16) and the RM union (E18) take the same union over $a\in\mathcal V$; the cross-layer rule takes $\alpha_C+\alpha_M$ (6a Remark 4.2).
(c) **R26 on changes.** With calibrated, independent, nondifferential matches, the true composition is $X=\bar q+u$ with $E[u\mid q]=0$ and $\operatorname{Var}(u\mid q)=\Sigma_u$ known, and
$$E[\Delta T_c\mid q]=\bar q^\top\Delta Q_c\bar q+\operatorname{tr}(\Delta Q_c\Sigma_u),$$
a known correction involving only the identified $\Delta Q_c$. So FD, its single restrictions and DiD-S carry over with the corrected regression; the model DiD's weights use the true $X$ and need their own allowance (as in R26).
*Proof.* (a): on $\{a^\circ\in\mathcal V\}$, $\varphi(n^\circ)\in\varphi(F(a^\circ))\subseteq\mathcal C$ [E1; 6a:R19(a)]; $\mathcal V\subseteq\mathcal V^T$ along the diagonal [Def]. (b): the same argument. (c): $E[X^\top\Delta QX\mid q]=\bar q^\top\Delta Q\bar q+E[u^\top\Delta Qu\mid q]$ and $E[u^\top\Delta Qu\mid q]=\operatorname{tr}(\Delta Q\Sigma_u)$. [CE-lin; alg] ∎ **Status: P; proof review pending** for (a)–(c) as stated; the construction of $\mathcal V$ (Poisson-binomial or joint sets, calibration slack) is R25(b)–(c) and is deferred with the module. Checks: G7 [F] (a), G7 [V] (a) one set vs separate sets, G7 [F] (c) (exact Poisson-binomial expectation; the uncorrected fit is biased by up to $8\times10^{-3}$).

**G8 (covariates; low priority).** Covariates enter at three places: (i) strata published as units in the count layer (exact within strata); (ii) the causal layer, as trend assumptions within strata, aggregated with the **treated** group's stratum shares; (iii) the model layer, as CAR(Z)-δ on changes or FD with covariates in the mean (conditional parallel trends: **O**, NS Part 3A).

**G9.** If the group is defined at each election (for example a party registration), the DiD of period-specific group rates is not a within-person comparison; switchers mix into both groups. Under (FP) with a fixed definition (veterans defined after the war are "future veterans" before it) the item does not arise.

---

## 9. What can be said, under which assumptions (plain language)

*No formulas. The question is the owner's: "EI DiD with a causal reading: what can we say, and with what assumptions?"*

### 9.1 The two questions

* **The descriptive question.** Among the same people, did one group's support change more than another's between two elections? This is the **difference-in-differences (DiD)**. We never see how each group voted, only each place's group shares and totals, so even this is a range.
* **The causal question.** Did the treatment (say, war service) change the treated group's support? That compares what the treated group did with what it **would have done** without treatment, which nobody observes. Every causal number therefore rests on an assumption about that unseen path, on top of the descriptive range.

The two must be reported separately. A causal number needs four things stated: **who is compared with whom**; **how far the two groups' untreated paths could have differed** (the trend assumption); **whether treatment could have affected the comparison group** (spillovers); and that **nobody changed before treatment** (no anticipation).

### 9.2 The ladder, with a worked example

**Example data** (check W). Four precincts of 200 men each, with 10, 50, 150 and 190 veterans. Three elections: an early pre-war election, the last pre-war election (the base) and a post-war election. Yes votes:

| Precinct (veterans of 200) | Early pre-war | Last pre-war | Post-war |
|---|---|---|---|
| 1 (10) | 60 | 62 | 74 |
| 2 (50) | 62 | 66 | 86 |
| 3 (150) | 61 | 63 | 113 |
| 4 (190) | 56 | 61 | 120 |

Support rose most where veterans were common. The question: by how much more did veterans' support rise than others', and how much of it did service cause?

| What you assume | What you can say (veterans' change minus others', or the effect of service) | Result |
|---|---|---|
| **Nothing** beyond the counts | *Descriptive*: veterans' support changed by anywhere from 41 points less to 79 points more than others'. *Causal*: the effect of service is anywhere in the logical range; the data alone never sign it (Prop. E15) | DiD in $[-0.407,0.792]$; effect in $[-0.568,0.732]$ |
| **Contextual reinforcement on the changes**, a belief about these voters: the within-precinct gap in changes and the pattern across precincts point the same way (Prop. E5(d)) | *Descriptive*: veterans' support rose 14 to 27 points more than others' | DiD in $[0.143,0.270]$ |
| + **parallel trends**, no spillovers, no anticipation | *Causal*: service raised veterans' support by 14 to 27 points | effect in $[0.143,0.270]$ |
| + instead a **trend dial**: untreated trends could differ by up to 5 points (Cor. E17) | 9 to 32 points | $[0.093,0.320]$ |
| + a **spillover allowance**: service could have moved others by up to 2 points (Prop. E22) | 7 to 34 points | $[0.073,0.340]$ |
| **Relative magnitudes** ($\bar M=1$: the post-war violation is no larger than the largest pre-war one) with the counts alone for the pre-war elections (Prop. E19) | nothing: the pre-war DiD is itself a range from $-0.60$ to $0.60$, so the allowance swamps everything | $[-0.568,0.732]$ (no information) |
| Relative magnitudes with contextual reinforcement also imposed on the pre-war change | the pre-war DiD is pinned to $[0.004,0.008]$, so the allowance is tiny: 14 to 28 points | $[0.135,0.278]$ |
| **Monotone response + no decline + positive selection**, no parallel trends (Prop. E21) | service did not lower support, and raised it by at most 48 points (31 with reinforcement) | $[0,0.482]$; with reinforcement $[0,0.311]$ |

### 9.3 What the theory says, in eight sentences

1. **The arithmetic range of the DiD is the sum of the elections' ranges, and it is the best possible without further beliefs** (Prop. E3). With more elections, the event-study DiDs are linked through the shared base election; treating them as separate boxes overstates what is possible.
2. **Beliefs about the actual voters' changes (signs, contextual reinforcement, comparisons between kinds of places) narrow the DiD exactly**, with no model and no sampling error (Props E5, E6). Contextual reinforcement on changes is never contradicted by the counts in the continuous problem and signs the DiD by Goodman's regression on changes; with whole-ballot counts it can be contradicted (Prop. E26(c)).
3. **Straight-line models of the changes (FD) allow any persistent local culture**, and give an exact range (Prop. E8) that is usually narrower than the earlier outer range; one extra rule (R1, R2 or R4) gives a single number, and the "ordering plus spillover-sign" pair (DiD-S) gives a closed-form range (Prop. E9).
4. **A model range and an arithmetic range combine by overlap** for the actual DiD, with the model's confidence level; their disagreement is a test of the model and the beliefs together (Theorem E12).
5. **Without a trend assumption the data say nothing about the effect of treatment** (Prop. E15).
6. **Relative-magnitude bounds help only when something narrows the pre-period DiD**: a model, a belief about the pre-period changes, or known pre-period counts. On the arithmetic alone they reduce to a trend dial whose size is at least the relative magnitude times half the pre-period's arithmetic width, and usually to no information (Prop. E19). Smoothness is worse: extrapolating a pre-period range counts the base election twice (Prop. E20).
7. **Spillovers onto the comparison group cannot be learned from pre-periods**, because they did not exist then; they need their own allowance (Prop. E22).
8. **People entering, leaving, aging in, censuses that disagree, and boundaries that move are not solved here.** Each is listed in the gap register (§8) with what still holds, what would close it, and where a range with stated confidence would enter.

---

## 10. Proof-register rows (for PD §8)

| ID | Claim | Status |
|---|---|---|
| R11a | Sharpness of the descriptive rungs: Δ0 for T elections (E3), Δ4 in R×C and T elections, closed forms for one-vs-rest rows, LP otherwise (E5), count contrasts on changes (E6), the sharp FD set via Lemma G′ (E8), DiD-S sharp for the regression and for the full law iff inside the sharp FD set (E9), LCE contrasts on changes (E10), per-period LCE with model rows on changes (E11) | **P; proof review pending** (2e §§2–3; `verify_phase2e.py` E1–E11). Grouped rungs with model rows: R16 (O) |
| R11b | Validity and sharpness of the causal layer for each form: the identity ATT = D − b + S (E14), no causal content without a trend assumption (E15), the joint causal LP (E16), trend dial (E17), MTR/ND/PS with spillovers (E21), spillover dial (E22), labelling rule (Def. 5.1) | **P; proof review pending** (2e §5; E14–E17, E21, E22) |
| R11b′ | Relative magnitudes and smoothness over the joint set of latent tables: union-of-LP form (E18); degeneration (a) inclusion in a trend dial with equality iff compatible maximisers, (b) the vacuity criterion, (c) monotonicity, the margins-only allowance $\ge\bar Mw_{-1}/2$, standard RR with point-identified pre-trends, (d) one joint region (E25); smoothness doubles the base width (E20); N21 exact and corrected | **P; proof review pending** (2e §5.5; E18–E20, E24) |
| R11c | Inference: FD with linear rows (E23(a)–(b)); FD4 with polynomial τ (E23(c)); unions across branches and elections (E24); one joint region for model-based pre-trends (E25); two-election allowance (E13); cross-layer rule on the product (E12) | **P; proof review pending**; inference for the sharp FD set **O**; CLR versions algorithmic (P3.5(b3)) |
| R11d | Integer versions: lattice for margins and one-vs-rest rows on changes (E26); TU for change rows on one column over paths of period pairs with laminar unit sets (E27); panel tables add nothing (E2) | **P; proof review pending**; wider classes (three pairs sharing a base; several columns; mixtures) **C**/MILP; causal integer ends **O** (Remark 7.1) |
| R31 | DiD gap register (G1–G9): for each item, the result that holds without it, the closing assumption or data range, the rungs modified, and a check that the naive statement fails; Lemma G6 (coarsening) | register **complete**; Lemma G6 **P; proof review pending**; the items themselves **O** (deferred to the uncertain-totals module) |
| R32 | Margin sets in T elections: one set for a shared composition (G7(a)), the causal and cross-layer versions (G7(b)), R26 on changes (G7(c)) | **P; proof review pending** for G7(a)–(c) as stated; R25(b)–(c) and the module **O** (deferred) |

---

## 11. Issues found in other documents (not edited)

1. **PD §7, row "M-ORD-Δ (MTR)".** The ordering of changes $\Delta m_v(x)\ge\Delta m_n(x)$ is a **descriptive** row on the model-implied changes (one of the two DiD-S rows, E9). MTR is an assumption about potential outcomes (E21). The two are linked only through a trend assumption: with parallel trends, no spillovers and no anticipation, MTR implies $\mathrm{ATT}\ge0$ and hence a non-negative DiD, but an ordering of changes implies nothing causal without such an assumption (E15). The label mixes the layers; NS §3.2 already flags this as the example of the problem. Relatedly, the row **"M-CEIL-Δ: $\lvert\Delta w_1\rvert\le s$ with $s$ benchmarked from pre-period variation (relative magnitudes)"** conflates the comparison group's contextual change $\Delta w_1$ (a descriptive model parameter) with the counterfactual trend gap $b$ that relative magnitudes bound (Def. 5.2); and the pre-period $\Delta w_1$ is not identified (1b B2 Δ3 already calls the calibration heuristic).
2. **2d Prop. CL.** (i) PS is stated with the **observed** comparison rate $B_{n1}$; that equals the untreated rate only with no spillover onto the comparison group (E21 with $\Lambda=0$). (ii) CL(d) "under parallel trends, ATT = DiD" needs no anticipation and no spillovers (E14). (iii) "ATT ∈ DiD ⊕ [−M, M] (Rambachan & Roth 2023)": a fixed bound $\lvert b\rvert\le M$ is a bounded-variation dial in the sense of Manski & Pepper (2018); Rambachan & Roth's sets are relative magnitudes and smoothness (Defs 5.2–5.3). The ⊕ form is also only outer: the sharp set intersects with the rate limits (E17(b); check E16: the limits bind in 80 of 120 cases).
3. **2d §11 and 1b B5, the actual-DiD allowance** $(\sqrt{m_1(1-m_1)}+\sqrt{m_0(1-m_0)})^2$. This Bhatia–Davis form is not proved for a DiD: P3.5(c1) covers single cells only, and 6a Lemma 4.1 says the Bhatia–Davis version must not be used for multi-cell targets. The proved allowance is Popoviciu's with $\rho_i=2(w_{ir}+w_{is})$ (Lemma E13), which also covers any dependence between the two groups' changes within a unit (2d's formula is per group).
4. **2d Prop. Δ4(b)** ("CR on changes: never empty, and sign DiD = sign $ER_\Delta$") and **PLAIN_PHASE2d §3** ("This assumption can never be rejected by the data"): true for continuous tables with no unit rows; the sign statement also needs $\gamma<1$ and $ER_\Delta\ne0$ (as 2c M3(d) was corrected after review). With whole-ballot counts CR on changes can be empty (E26(c); check E3: one design of 59), and unit-level rows can remove the equal-rates witness (2c check C8b).
5. **2d Prop. FD4** is an outer set; the sharp set (E8) is typically strictly narrower (check E8: 100 of 100 designs). 2d does not claim sharpness, but PLAIN_PHASE2d §2's comparison "not guaranteed to be narrower" than Δ0 should be read with the sharp set available as the narrower option.
6. **NS §3.2, N21.** (i) The ATT bound is compared with the DiD's logical range $[-2,2]$; an ATT of a rate lies in $[B^L_{v,1}-1,B^U_{v,1}]$, here $[-0.75,1.00]$ (E15), and the sharp RM set with the rate limits is exactly that interval (check E18 [F]). The conclusion ("no information") stands. (ii) NS writes the RM bound with $\max_{t\le0}\lvert D_t\rvert$ (levels relative to the base); Rambachan–Roth use consecutive differences; the two coincide only with one pre-period (check E19 [I]). (iii) "The ATT set is the descriptive set shifted by the trend dial's set" (NS §3.2, technical notes) is the outer form; the sharp set also imposes the rate limits, and a spillover term is needed (E16, E22).
7. **NS Part 3A, count block, "DiD descriptive × estimated margins: L·a·R25 (both periods)".** Under a fixed population the composition is shared, so one margin set must serve both elections (G7(a)); separate per-election sets are valid but wider and, if separately random, spend α once per election.
8. **NS Part 3A, "CONTRAST-LCE × FD/DiD-S: … sharpness O (R11a)"** and **PD §8 R11 "rest O"**: now characterised (E9, E10, E11, E12); status updates, not errors.
