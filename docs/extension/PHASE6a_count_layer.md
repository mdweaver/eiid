# Phase 6a: the count layer (M-COUNT), proofs for build step 1

*Proof-register items **R19** (the count-layer bound, integrality, witnesses) and **R27** (the cross-layer combination rule) from `next_steps/NEXT_STEPS_RECOMMENDATION.md` (cited as NS) §3.1 and Part 3A, with decisions D1–D6 of NS Part 7. Builds on Duncan–Davis (P1; 2a §2), `PHASE2c_monotone.md` (2c: Lemma M1, Prop. M2, Remark M2.1, Prop. M3), `PHASE3_inference.md` (P3: §0.2, Thm P3.5, Prop. P3.5(c)) and `PACKAGE_DESIGN.md` (PD: §1.2, §4.2, §8). `PHASE2d_did.md` (Δ0, Δ4) is context only: DiD is step 2. Notation follows `notation.tex`; the latent cell count written $N_{irc}$ there is written $n_{irc}$ here, as in NS.*

*Theory only. Nothing in the package or in other theory documents is changed; problems found in them are listed in §9. Checks: `verify_phase6a.py` (output `verify_phase6a_output.txt`, **36/36 PASS**), tagged as in 2c and 5a: **[F]** fails on an earlier or tempting statement and passes on the one made here; **[V]** verifies a statement made here; **[I]** illustrates.*

*Status in brief (after the proof review, `reviews/PHASE6a_proof_review.md`; its 15 issues are fixed as listed in §10). Every proved item is labelled **P; proof review pending** until the reviewer confirms the fixes. R19(a)–(e): integrality for TU classes (i)–(iii) (TU of the original matrix for (i)–(ii)); LP outer and MILP exact outside them. Prop. C2 (realised one-vs-rest signs and CR in counts, one sign set, integer versions). Lemma 4.1 (asymptotic and Hoeffding residual allowances for linear count targets with design-measurable coefficients). R27: union bounds $(*)$; exact level-α coverage and test with identity rows only (margins, published totals, audited cells); with belief rows, coverage $\ge1-\alpha-P(\text{rows false})$ and P(reject while the rows hold) ≤ α. R27(iii-a) (compatible values): P; R27(iii-b) (non-minimality): a computed example. What remains of R15 becomes a computed check (§3.4).*

---

## 0. Tags used here

Tags of 2a §0.3, 2c §0.3 and P3 §0.3 apply, plus:

| Tag | Meaning |
|---|---|
| [NW] | a transportation problem with nonnegative row sums $a$ and column sums $b$, $\sum a=\sum b$, has a nonnegative solution (north-west corner rule); integer data give an integer solution |
| [TU] | a matrix is totally unimodular (TU) if every square minor is in $\{-1,0,1\}$. Appending identity rows, negating rows and repeating rows preserve TU (Schrijver 1986, *Theory of Linear and Integer Programming*, §19.1†) |
| [GH] | Ghouila-Houri: a matrix is TU iff every subset of its rows can be split into two parts whose row-sum difference has every entry in $\{-1,0,1\}$ (Schrijver 1986, Thm 19.3(iv)†) |
| [HK] | Hoffman–Kruskal: if $A$ is TU and $b,\ell,u$ are integer vectors, every vertex of $\{x:Ax\le b,\ \ell\le x\le u\}$ (equalities written as two inequalities) is integral (Schrijver 1986, Cor. 19.2a†) |
| [vert] | a nonempty polytope has a vertex, and a linear function attains its minimum and maximum over it at vertices; a basic optimal solution returned by the simplex method is a vertex |
| [Chv] | for a row $g^\top x\ge h$ with integer $g$, every integer $x$ satisfying it satisfies $g^\top x\ge\lceil h\rceil$ (Chvátal rounding of one row) |
| [union], [Pop], [Lyap], [Chebyshev] | as in P3 §0.3 |

---

## 1. Setup

### 1.1 Data and unknowns

* Units $i=1,\dots,p$ with $N_i>0$; units with $N_i=0$ are dropped (Lemma C1 divides by $N_i$). (Precincts, or sub-unit strata when a unit's margins are published by stratum, NS N9). Unit $i$ has size $N_i$, row (group) counts $a_{ir}=N_iX_{ir}$ and column (outcome) counts $b_{ic}=N_iT_{ic}$, with $\sum_ra_{ir}=\sum_cb_{ic}=N_i$.
* **Unknowns.** The realised cell counts $n_{irc}\ge0$ (people in unit $i$, group $r$, outcome $c$). They satisfy the margins
$$\textstyle\sum_cn_{irc}=a_{ir},\qquad\sum_rn_{irc}=b_{ic}.$$
The set of nonnegative real tables with unit $i$'s margins is the **transportation polytope** $\mathcal P_i$. It is nonempty [NW] and bounded ($0\le n_{irc}\le\min(a_{ir},b_{ic})$). Write $\mathcal P=\prod_i\mathcal P_i$ and $n=(n_{irc})$.
* **Dictionary with 2c.** For $a_{ir}>0$, $\beta_{irc}=n_{irc}/a_{ir}$, and $n\in\mathcal P$ iff $\beta_i\in S(X_i,T_i)$ for every $i$ (2c M2 proof step 1). Everything below can be read in either coordinate system; counts are used because the rows of §1.3 are linear in $n$ with integer coefficients, which §2.2 needs.
* Aggregates: $A_r=\sum_ia_{ir}$, $Y_c=\sum_ib_{ic}$, $N=\sum_iN_i$, $Y_{rc}(\mathcal S)=\sum_{i\in\mathcal S}n_{irc}$, $Y_{rc}=Y_{rc}(\text{all units})$.

### 1.2 Targets

A **count target** is any linear functional of the table with known coefficients,
$$\varphi(n)=\varphi_0+\sum_i\sum_{r,c}v_{irc}\,n_{irc},$$
usually a count over a known denominator. The count layer (R19) allows any known coefficients. For the cross-layer rule (R27, Lemma 4.1) the coefficients $v_{irc}$ must be functions of the design $\{(N_i,X_i,D_i)\}$ only, so a denominator built from outcome margins (e.g. $Y_c$) is excluded there. Examples:
* a rate over a set of units: $B_{rc}(\mathcal S)=Y_{rc}(\mathcal S)/A_r(\mathcal S)$ (this is notation.tex's $B_{rc}(\mathcal S)$; $\mathcal S$ may be a county, the units above some size, ...);
* a gap: $B_{rc}-B_{sc}$, or a one-vs-rest gap $D^r_c=B_{rc}-B_{\neg r,c}$ with $B_{\neg r,c}=(Y_c-Y_{rc})/(N-A_r)$;
* a merged-group rate, a turnout-type sum over columns, any $\sum a_{rc}B_{rc}$ (2c M2), and differences of rates across unit sets.

**The target is the actual (realised) quantity.** For a rate this is P3 §0.2's finite-population rate $B_{rc}=\sum_iN_iX_{ir}\beta_{irc}/\sum_iN_iX_{ir}=Y_{rc}/A_r$, a function of the realised table. It is **not** the model-implied rate $\tilde B_{rc}=\sum_iw_{ir}m_{rc}(V_i)$, which is not a function of the table at all, and it is not the superpopulation rate $B^*$. The count layer says nothing about $\tilde B$ or $B^*$ (§5). When the target has non-integer coefficients (rates, gaps), its **numerator** $\sum v_{irc}n_{irc}$ times a common denominator has integer coefficients; integrality statements below refer to that.

### 1.3 Restriction rows

All count-layer restrictions are linear rows on $n$ (or unions of row systems), with coefficients computed from the margins. They are beliefs about **the realised table**, i.e. `object = "actual"` in PD §1.2's vocabulary.

| Row | Form in counts | Source |
|---|---|---|
| Structural zero / known cell | $n_{irc}=0$, or $n_{irc}=k$ | MENU A1; E5 (audit, poll book) |
| Cell box | $\ell_{irc}\le n_{irc}\le u_{irc}$ | user bound; the next three rows reduce to it |
| Realised one-vs-rest ordering in unit $i$ (column $c$) | $n_{irc}/a_{ir}\ge(b_{ic}-n_{irc})/(N_i-a_{ir})\iff n_{irc}\ge a_{ir}b_{ic}/N_i$: **a box** ($0<a_{ir}<N_i$). For $R=2$ this is the pairwise ordering | NS N7; MENU B5 |
| Realised pairwise ordering in unit $i$ ($R\ge3$) | $a_{is}n_{irc}-a_{ir}n_{isc}\ge0$: two cells, not a box | MENU B5 |
| Odds-ratio dial (2×2 unit) | $\psi_L\le n_{i11}n_{i22}/(n_{i12}n_{i21})\le\psi_U$; the odds ratio is increasing in $n_{i11}$ given the margins, so this is **a box** on $n_{i11}$ (convention: if a margin of the unit is 0 the cell is fixed by the margins and the dial is ignored; the endpoints of the cell's range, where the ratio is 0 or ∞, belong to the box iff $\psi_L=0$ or $\psi_U=\infty$) | NA §12.3(c); NS N17 |
| Published total over a unit set | $t^-\le Y_{rc}(\mathcal S)\le t^+$ (equality if exact) | E5, NS §3.7c |
| Realised association signs (2c M1) | rows on $n$ (Lemma C1): **one-vs-rest** $\delta^r_B,\delta^r_W$ over a unit set $\mathcal S$ are rows on the single aggregate $Y_{rc}(\mathcal S)$, i.e. totals; pairwise $\delta_B(r,s)$ and W2 $\Delta_r$ ($R\ge3$) put unequal coefficients on different cells and units | 2c M1, M2 |
| Contextual reinforcement (CR) | union of two row systems $F_+\cup F_-$ | 2c Def. 2 |
| Slack dials | any row $g^\top n\ge h$ replaced by $g^\top n\ge h-\varepsilon$, $\varepsilon\ge0$ | 2c M2(d) |

Write $F$ for the set of real tables in $\mathcal P$ satisfying all non-union rows ($F_\pm$ for the two CR branches), and $F^{\mathbb Z}=F\cap\mathbb Z^{p\times R\times C}$.

### 1.4 Integer data and the rounded system

When the margins are integers, every realised table is an integer table. The **rounded system** $F_{\rm rd}$ replaces
* each bound $\ell\le n_{irc}\le u$ by $\lceil\ell\rceil\le n_{irc}\le\lfloor u\rfloor$;
* each row with rational coefficients by the same row scaled to coprime integer coefficients, with its right-hand side rounded inward (⌈·⌉ for $\ge$, ⌊·⌋ for $\le$) [Chv];
* each equality row whose scaled right-hand side is not an integer by an infeasible row.

**Lemma 1.1.** With integer margins, $F_{\rm rd}\cap\mathbb Z^{p\times R\times C}=F^{\mathbb Z}$, and $F^{\mathbb Z}\subseteq F_{\rm rd}\subseteq F$.
*Plain reading.* Rounding each limit to whole numbers loses no whole-number table and admits no new one.

*Proof.* Each rounding step replaces a row by one with the same integer solutions [Chv; alg]; for bounds this is the one-variable case. The rounded rows are implied by the original ones, so $F_{\rm rd}\subseteq F$. ∎

*Check.* B1, T5, N10b. **Status: P; proof review pending.**

---

## 2. R19: the count-layer bound

### 2.1 R19(a): the bound is an LP, sharp for continuous tables

**Proposition R19(a).** Let the margins be given (integer or not) and let $\varphi$ be a count target.
1. $F$ is a polytope, and the set of values of $\varphi$ over all real (continuous) tables consistent with the margins and the rows is $\varphi(F)$. If $F\ne\emptyset$ it equals $[\min_F\varphi,\max_F\varphi]$, both ends attained.
2. **Valid:** if the margins are the true margins and every row holds for the realised table $n^\circ$, then $\varphi(n^\circ)\in\varphi(F)$.
3. **Sharp for continuous tables given the stated information:** every value in $\varphi(F)$ is $\varphi(n)$ for a table $n$ that reproduces the margins and satisfies every row. So no proper subset of $\varphi(F)$ is valid for every realised table consistent with the information.
4. $F=\emptyset$ iff the rows contradict the margins in the continuous problem. Under CR the set is $\varphi(F_+)\cup\varphi(F_-)$, at most two intervals. With slack dials the sets are nested componentwise in $\varepsilon$ with continuous piecewise-linear endpoints where nonempty (2c M2(d)).

*Plain reading.* With only the precinct totals and your stated beliefs about these particular voters, the exact range is what a linear program over all tables fitting the totals and beliefs returns. Every number in the range is produced by some table that fits; the true answer is in it if the totals and beliefs are right.

*Proof.*
1. $\mathcal P$ is a nonempty polytope [NW; §1.1], and the rows are finitely many linear (in)equalities, so $F$ is a polytope (possibly empty). [Def]
2. $\varphi$ is linear, $F$ compact and convex, so $\varphi(F)$ is a compact interval with attained ends. [Weier; conn]
3. Item 2 of the statement: $n^\circ\in\mathcal P$ (true margins) and satisfies the rows, so $n^\circ\in F$. [Def]
4. Item 3: by definition every $n\in F$ is consistent with all the stated information. [Def]
5. Item 4: emptiness is the definition; CR applies 1–4 to each branch; slack dials enter only right-hand sides, so 2c M2 step 6 applies verbatim. [2c:M2(c),(d)] ∎

*Relation to 2c M2.* In $\beta$ coordinates (§1.1) $F$ is 2c's N1 feasible set with unit restrictions $U_i$, plus the rows that couple units (totals). R19(a) is M2(a)–(d) for general linear count targets and count rows. **Status: P; proof review pending** (standard; NS called it standard).

### 2.1.1 Duncan–Davis as the margins-only case

**Lemma 2.1.** Margins only (no rows).
(a) For a single cell over a unit set, $Y_{rc}(\mathcal S)$ ranges over $\big[\sum_{i\in\mathcal S}\max(0,a_{ir}+b_{ic}-N_i),\ \sum_{i\in\mathcal S}\min(a_{ir},b_{ic})\big]$ (Fréchet bounds summed: the Duncan–Davis interval of P1/2a), every value attained.
(b) For a one-column contrast $\varphi=Y_{rc}/A_r-Y_{sc}/A_s$ the range is obtained by cellwise interval arithmetic: $[\underline B_{rc}-\overline B_{sc},\ \overline B_{rc}-\underline B_{sc}]$.
(c) For a target with several cells of one row, or several same-sign cells of one column (including the one-vs-rest gap $D^r_c$ for $R\ge3$, written over cells), cellwise interval arithmetic (what `eiBounds()` prints as `duncan_davis` for a combination target, `.rxc_interval_lin`) is an **outer** bound that can be strictly wider than the exact LP range.

*Plain reading.* With margins only, the exact range of one group's count for one outcome, or of the gap between two groups on one outcome, is the familiar Duncan–Davis range. For targets that add cells of one group across outcomes, or merge groups, adding up the cell ranges overstates the range.

*Proof.*
1. The problem separates across units, so ranges add. [sep]
2. Projection of $\mathcal P_i$ onto column $c$: it is $\{0\le n_{\cdot c}\le a_{i\cdot},\ \sum_rn_{irc}=b_{ic}\}$. Indeed given column $c$, the remaining $R\times(C-1)$ table has nonnegative row sums $a_{ir}-n_{irc}$ and column sums $b_{ic'}$ with equal totals, hence a solution. [NW]
3. (a): on that projection, $n_{irc}\le\min(a_{ir},b_{ic})$ and $n_{irc}=b_{ic}-\sum_{s\ne r}n_{isc}\ge b_{ic}-(N_i-a_{ir})$. Both bounds are attained: put $n_{irc}$ at the bound and spread the remainder $b_{ic}-n_{irc}$ over the other rows within their capacities $a_{is}$, whose sum is $N_i-a_{ir}$. The remainder is $\max(0,b_{ic}-a_{ir})\le N_i-a_{ir}$ at the upper bound (as $b_{ic}\le N_i$) and $\min(b_{ic},N_i-a_{ir})\le N_i-a_{ir}$ at the lower bound. Intermediate values by convexity. [2; conn]
4. (b): within a unit, the maximum of $n_{irc}$ is $\min(a_{ir},b_{ic})$. Given it, the least feasible $n_{isc}$ is $\max(0,b_{ic}-n_{irc}-(N_i-a_{ir}-a_{is}))$, which equals the unconditional Fréchet minimum $\max(0,b_{ic}-N_i+a_{is})$ in both cases $a_{ir}\le b_{ic}$ and $a_{ir}>b_{ic}$. So both separate optima are attained simultaneously; the minimum is symmetric. [2; alg]
5. (c): interval arithmetic is outer because each cell's range is taken separately. One unit with rows $(3,2,0),(2,3,0),(0,0,5)$: the merged rate of rows 1+2 in column 1 is exactly $[0,0.5]$, cellwise $[0,1]$; row 1's share in columns 1+2 is exactly $[0,1]$, cellwise $[0,2]$ (D3). ∎

*Checks.* D1 (1296 single-cell targets on random integer data: LP = Fréchet, max error $2\times10^{-16}$); D1b (the package's reference datasets `rxc_example_3x3.csv`, `msign_example_3x2.csv`, non-integer counts: LP = the package's `.rxc_dd_cells` formula, max error $3\times10^{-16}$); D2 (one-column gaps, 150 designs); **D3 [F]** (multi-cell targets). **Status: P; proof review pending.**

*Consequence for D1 (DD first).* The block printed first by `eiBounds()` (decision D1) must be the margins-only **LP** value. It coincides with the current `duncan_davis` column for single cells and one-column contrasts (Lemma 2.1(a)–(b)); for other combination targets the current column is outer and the LP is the exact margins-only bound (§9, item 3).

### 2.2 R19(b): total unimodularity and integer endpoints

**Theorem R19(b).** Let the margins be integers. Suppose every non-union row is of one of the following kinds (after writing ordering and odds-ratio rows as boxes, §1.3), with **all** the totals in **one** of the three classes below (a mixture of classes is not covered and can fail; review check RC4):
* bounds $\ell\le n_{irc}\le u$ on single cells (boxes, known cells, structural zeros);
* **totals** $t^-\le Y_{\kappa}(\mathcal S)\le t^+$ on a cell $\kappa=(r,c)$ over unit sets $\mathcal S$, organised in one of three classes:
  * **(i) one cell:** all totals are on the same cell $\kappa$, and their unit sets form the union of **two laminar families** (a family is laminar if any two of its sets are nested or disjoint);
  * **(ii) one column:** totals on cells $(r,c_0)$ of one column $c_0$ (any rows $r$), with the unit sets of each cell forming **one** laminar family (by transposition: cells of one row, one laminar family each);
  * **(iii) two rows:** $R=2$, totals on any cells, where for each column $c$ the unit sets used on that column's two cells form one laminar family (by transposition: $C=2$ with rows and columns exchanged).

Then $F_{\rm rd}$ is an integral polytope. In classes (i) and (ii) the margins-plus-totals matrix (bounds are not rows) is TU. In class (iii), after the integral affine substitution $n_{i2c}=b_{ic}-y_{ic}$, the matrix of the reduced system is TU; the original matrix need not be (review check RC1: a 7×7 minor equal to 2; T6 here). Consequently, for every linear target:
1. $F^{\mathbb Z}=\emptyset$ iff $F_{\rm rd}=\emptyset$;
2. $\min_{F^{\mathbb Z}}\varphi=\min_{F_{\rm rd}}\varphi$ and $\max_{F^{\mathbb Z}}\varphi=\max_{F_{\rm rd}}\varphi$, each attained at an integer vertex of $F_{\rm rd}$; for integer-coefficient numerators the endpoints are integers;
3. for $\varphi=Y_\kappa(\mathcal S)$ with $\mathcal S$ addable to one family of the class (e.g. a single unit, or all units in class (i)), **every integer between the endpoints is attained**. For other targets the set of integer values is a finite subset of the interval whose hull is the interval.

The LP over the **unrounded** system can have fractional endpoints (N10b), and rounding that endpoint is **not** enough when several rows are fractional: the rows must be rounded first (T5).

*Plain reading.* For margins, known counts, cell limits, unit orderings, odds-ratio limits, published totals and the one-vs-rest association signs, the linear program (after rounding each row's limit to whole numbers) already returns whole-number tables at both ends. No integer program is needed, and the ends it reports are attained by real ballots-and-people tables.

*Proof.* Write the margins as rows $E^{\rm row}_{ir}$ (support $\{(i,r,c)\}_c$) and $E^{\rm col}_{ic}$ (support $\{(i,r,c)\}_r$), and a total as a row $T_{\mathcal S,\kappa}$ (support $\{(i,\kappa)\}_{i\in\mathcal S}$). Bounds are appended identity rows, so for classes (i) and (ii) it suffices to show that the matrix of $E$ and $T$ rows is TU (steps 1–4). [TU] For class (iii) step 5 shows the same for the reduced system after an integral affine substitution, which maps the vertices of one polytope onto those of the other. Integrality then follows from [HK], since every right-hand side and bound of $F_{\rm rd}$ is an integer (§1.4), and items 1–2 from [vert] (the polytopes are bounded, $0\le n\le\min(a,b)$) and Lemma 1.1.

Classes (i), (ii): fix any subset of rows. We give a signing as in [GH].
1. *Margins.* For each unit choose an orientation $\varepsilon_i\in\{\pm1\}$: chosen $E^{\rm row}_{ir}$ get sign $\varepsilon_i$, chosen $E^{\rm col}_{ic}$ get $-\varepsilon_i$. The variable $(i,r,c)$ then receives $\varepsilon_i(e^{\rm row}_{ir}-e^{\rm col}_{ic})\in\{-1,0,1\}$ from the margins, where $e=1$ if the row was chosen. [alg]
2. *A laminar family, signed by depth.* Among the chosen totals of one laminar family, give a set the sign $(-1)^d$, where $d$ is the number of chosen sets of the family that strictly contain it. For a unit $i$ the chosen sets containing $i$ form a chain $\mathcal S_1\supset\mathcal S_2\supset\dots\supset\mathcal S_k$, and the depth of $\mathcal S_j$ is $j-1$ (any chosen set strictly containing $\mathcal S_j$ contains $i$, so lies in the chain). The family therefore contributes $\tau_i=\sum_j(-1)^{j-1}\in\{0,1\}$ to the variable $(i,\kappa)$ of its cell. Reversing all its signs gives $\{0,-1\}$. [Def laminar]
3. *Class (i).* Sign family A as in 2 (contribution $\tau^A_i\in\{0,1\}$) and family B reversed ($\tau^B_i\in\{0,-1\}$), so $\tau_i=\tau^A_i+\tau^B_i\in\{-1,0,1\}$ on the variable $(i,\kappa)$ and 0 elsewhere. The margins give $\varepsilon_i\delta_i$ on $(i,\kappa)$ with $\delta_i=e^{\rm row}_{ir}-e^{\rm col}_{ic}\in\{-1,0,1\}$. If $\tau_i\ne0$ and $\delta_i\ne0$, take $\varepsilon_i=-\tau_i\delta_i$, so the total is 0; otherwise the total is in $\{-1,0,1\}$. Other variables of unit $i$ receive margin contributions only. The choice of $\varepsilon_i$ is per unit, so all choices are compatible. [1; 2]
4. *Class (ii).* Sign each cell's family as in 2, so the variable $(i,r,c_0)$ receives $\tau_{ir}\in\{0,1\}$ from the totals of cell $(r,c_0)$. Take $\varepsilon_i=+1$ if $E^{\rm col}_{ic_0}$ is chosen and $\varepsilon_i=-1$ otherwise. In the first case $(i,r,c_0)$ receives $(e^{\rm row}_{ir}-1)+\tau_{ir}\in\{-1,0\}+\{0,1\}$; in the second, $-e^{\rm row}_{ir}+\tau_{ir}\in\{-1,0\}+\{0,1\}$. Both are in $\{-1,0,1\}$ for every $r$ at once. [1; 2]
5. *Class (iii).* With $R=2$ put $y_{ic}=n_{i1c}$; then $n_{i2c}=b_{ic}-y_{ic}$, an affine bijection with integer offsets, so it maps integer points to integer points and vertices to vertices. The system becomes: unit rows $\sum_cy_{ic}=a_{i1}$; bounds $0\le y_{ic}\le b_{ic}$ (and the boxes); totals $\sum_{i\in\mathcal S}y_{ic}\in[t^-,t^+]$, a total on the cell $(2,c)$ becoming one on $y_{\cdot c}$ with right-hand side $[\sum_{\mathcal S}b_{ic}-t^+,\sum_{\mathcal S}b_{ic}-t^-]$. The unit rows are a partition of the variables (a laminar family $\mathcal L_1$); the totals of column $c$ are a laminar family on the variables $\{(i,c)\}_i$, and families for different columns have disjoint supports, so all totals form one laminar family $\mathcal L_2$. Sign $\mathcal L_1$ positively (each variable in at most one chosen unit row: contribution $\{0,1\}$) and $\mathcal L_2$ by depth, reversed (contribution $\{0,-1\}$). The sum is in $\{-1,0,1\}$, so the reduced matrix is TU. (The GH signing of steps 1–4 cannot be reused on the original matrix here: with totals in several columns one orientation $\varepsilon_i$ must serve all columns at once.) [2; GH]
6. *Item 3.* Adding the row $Y_\kappa(\mathcal S)=k$ keeps the matrix in the class (for (iii), in the reduced system) (it is one more set of a laminar family; a single unit is laminar with every family). For integer $k$ between the endpoints the continuous system is nonempty by convexity, so it has an integer vertex. [1–5; conn; HK] ∎

**Corollary 2.2 (named rows in the classes).** With integer margins the following are integer-exact by the LP over $F_{\rm rd}$, provided the whole set of non-union rows (per CR branch) lies in **one** class of R19(b); each $\delta^r_W$ row needs $\operatorname{Var}_p(X_r)>0$ over its own unit set:
* margins, structural zeros, known cells, cell boxes;
* realised one-vs-rest unit orderings (all $R$), hence all unit orderings when $R=2$; the odds-ratio dial on 2×2 units (boxes, §1.3);
* published totals on one cell over two laminar families of unit sets (for instance totals by county, nested in a state total, plus totals by a second, crossing geography);
* the realised one-vs-rest signs $\delta^r_B(\mathcal S)$, $\delta^r_W(\mathcal S)$ on a cell $(r,c)$ over unit sets forming at most two laminar families, together with totals on the same cell (Lemma C1: they are totals); each CR branch likewise (two totals on the same set);
* one-vs-rest signs and totals on **several rows of one column** (or several columns of one row), with one laminar family of unit sets per cell (class (ii));
* for $R=2$, these signs in **several columns jointly** (class (iii));
* slack dials on any of the above (they change right-hand sides only).

*Checks.* B1 (60 tiny designs, all integer tables enumerated: integer endpoints = LP over the rounded system, 0 fractional witnesses); T6 (the class-(iii) example of review check RC1: the original matrix has a minor 2, yet rounded LP = MILP on 200 random right-hand sides); B2 (integer values of an aggregate form an integer interval); T1–T3 (classes (i)–(iii), 400 random designs each with fractional boxes and a known cell: LP = MILP at 800 endpoints each, 0 fractional dual-simplex vertices); T4 (exhaustive minors: a class-(i) matrix is TU); **T5 [F]** (two units with ordering boxes $n\ge0.2$: unrounded LP total 0.4, its ceiling 1, integer answer 2); N10a, N10b (reproduced from NS).

**Status: P; proof review pending.** (Integrality is proved for classes (i)–(iii); TU of the original matrix for (i)–(ii) only. (ii) and (iii) go beyond NS's "margins, known cells and per-unit boxes".)

### 2.3 R19(c): rows that break total unimodularity

**Proposition R19(c).** For any finite set of rows with rational coefficients and integer margins:
1. $F^{\mathbb Z}\subseteq F_{\rm rd}\subseteq F$, so the LP over $F_{\rm rd}$ (or over $F$) gives an **outer** interval for the integer values of $\varphi$;
2. the MILP over the integer variables $n$ decides integer feasibility and returns the exact integer endpoints;
3. outside the classes of R19(b) the inclusions can be strict: fractional vertices; integer endpoints beyond the rounded LP endpoints (pairwise orderings, two-cell totals, below); and $F^{\mathbb Z}=\emptyset$ while $F_{\rm rd}\ne\emptyset$ (crossing totals, N10c).

*Plain reading.* Beliefs that tie together different cells of a precinct in unequal proportions, or that weight precincts unequally, can leave the linear program reporting a range that no whole-number table attains, or even a non-empty range when no whole-number table fits at all. The package must then solve an integer program.

*Proof.* 1 is Lemma 1.1. 2 is the definition of the MILP. 3: the examples below are computed exactly in `verify_phase6a.py`. ∎

Examples (check C1 and N10c):
* **Pairwise unit ordering, $R=3$** ($a_{is}n_{irc}-a_{ir}n_{isc}\ge0$): integer endpoint beyond the rounded LP endpoint (LP 8.714, integer 7, for a maximisation).
* **Pairwise realised $\delta_B(r,s)$, $R=3$** (unequal unit coefficients, Lemma C1(b)): a fractional vertex of $F_{\rm rd}$ after Chvátal rounding of the row (LP 45.098, integer 45). Here the integer endpoint equals the rounded LP endpoint; a search of 20,000 further designs with one such row found no gap beyond rounding, so this example shows only that TU fails.
* **Totals on cells in different rows and columns over overlapping units, $R=3$, $C=3$** (with $C=2$ the transposed class (iii) applies; review check RC10): integer endpoint beyond rounding (LP −47.5, integer −46).
* **Three pairwise-crossing totals** (e.g. totals published by three geographies): three units with $N=2$, one veteran, one yes vote each, and each pair's veteran-yes total equal to 1. The LP witness puts $1/2$ in each cell; no integer table exists. The 0/1 matrix has a $3\times3$ minor equal to 2 (T4).
* **N10c** (NS, reproduced with named rows): three 2×2 units, $N=(13,10,14)$, $V=(2,8,9)$, yes $=(13,8,4)$, a published veteran-yes total of 11, and two **average-gap rows** $\sum_iw_i\big(y_i/V_i-(Y_i-y_i)/(N_i-V_i)\big)\ge0$ with $w_i=1$ (the unweighted mean of unit gaps) and $w_i=N_i-V_i$. The LP is feasible (witness veteran-yes $(2,\,6.372,\,2.628)$), it stays feasible after each row is rounded, and enumeration finds no integer table. These rows weight units unequally; they are **not** the realised associations of 2c M1, which are totals (§3, §9 item 1).

**Status: P; proof review pending.** The `integer = "auto"` rule (§6) follows: whole row set in one class of R19(b) → rounded LP; otherwise, with integer margins → MILP.

### 2.4 R19(d): the witness at each endpoint

**Proposition R19(d).** Let $\varphi$ attain an endpoint $e$ of $\varphi(F)$ (or of $\varphi(F_{\rm rd})$, or of the MILP).
1. A basic optimal solution $n^*$ of the LP (simplex, or an interior-point solution followed by crossover) is a vertex of $F$; it reproduces every margin, satisfies every row and bound, and has $\varphi(n^*)=e$.
2. Under R19(b), $n^*$ from $F_{\rm rd}$ is an **integer** table: a table of whole people attaining the endpoint.
3. The MILP solution is an integer witness of the integer endpoint.
4. Under CR, the witness of the branch attaining the endpoint satisfies that branch's rows, hence the union.
5. A witness certifies attainability (it is the content of "sharp"); it is one of many tables with the same target value and is not an estimate of the realised table.

*Plain reading.* Each end of the range comes with an actual table, precinct by precinct, that fits every total and belief and produces that end; with whole-number data it is a table of whole people.

*Proof.* 1 [vert]; feasibility of an optimal solution is its definition. 2 [HK] on $F_{\rm rd}$ with 1. 3 Def MILP. 4 $F_\pm\subseteq F_+\cup F_-$. 5 Def. ∎

*Check.* W1: all 17,693 LP and MILP solutions computed in §§2–3 of the script reproduce every margin and satisfy every row and bound; E1 prints the integer witnesses $(6,30)$ and $(20,50)$ of the realised-ordering example. **Status: P; proof review pending.**

### 2.5 R19(e): non-integer (estimated) margins

**Proposition R19(e).** Let some margin be non-integer (for example an expected count from probabilistic matching).
1. R19(a) holds verbatim **for the stated margins**: the LP range is the exact set of target values over real tables with exactly those margins and the rows, and every value is attained by such a table.
2. No integer table has a non-integer margin, so $F^{\mathbb Z}=\emptyset$ trivially and the integer check does not apply.
3. The range is a statement about the true table only if the stated margins are the true margins. For estimated margins that fails; the claim must read "exact for the stated margins, which are estimates". Rounding the margins to integers and running the integer check is not a substitute: it asserts information nobody has.

*Plain reading.* If the group sizes are estimates rather than counts, the range is exact for those estimated sizes and nothing more; the whole-number check is switched off because no whole-number table can have fractional totals.

*Proof.* 1: R19(a)'s proof never uses integrality. 2: a sum of integers is an integer. 3: R19(a) item 2 assumes the true margins. ∎

*Check.* E5 (an expected count 60 treated as exact gives $[0.167,0.833]$, NS N20; margin 59.4 gives a lower end attained by a real table with exactly those margins). Honest coverage under estimated margins is Step 1b (§7). **Status: P; proof review pending.**

---

## 3. Realised signs and CR on the count layer

### 3.1 The associations in counts

Fix a column $c$. With $p_i=N_i/N$, $X_{ir}=a_{ir}/N_i$, $\beta_{irc}=n_{irc}/a_{ir}$ and the definitions of 2c §1 (one-vs-rest weights, W1), put
$$Y^{NM}_{rc}(\mathcal S)=\sum_{i\in\mathcal S}\frac{a_{ir}b_{ic}}{N_i},\qquad\kappa_r=\frac1{A_r}+\frac1{N-A_r}=\frac N{A_r(N-A_r)},$$
and let $ER^r$ be Goodman's one-vs-rest regression with intercept (data only), all computed over the unit set $\mathcal S$ (all units by default; $N$, $A_r$, $Y_c$, $p$ then refer to $\mathcal S$).

**Lemma C1 (realised associations in counts).**
(a) $B_{rc}=Y_{rc}/A_r$ and $D^r_c=B_{rc}-B_{\neg r,c}=\kappa_rY_{rc}-Y_c/(N-A_r)$.
(b) $\delta_B(r,s)=\dfrac1N\sum_i\dfrac{a_{is}n_{irc}-a_{ir}n_{isc}}{N_i}$.
(c) $\delta^r_B=\dfrac{Y_{rc}-Y^{NM}_{rc}}N$. Hence $NM_r=Y^{NM}_{rc}/A_r$, and $\delta^r_B\ge-\varepsilon\iff Y_{rc}\ge Y^{NM}_{rc}-N\varepsilon$.
(d) If $\operatorname{Var}_p(X_r)>0$: $\delta^r_W=\operatorname{Var}_p(X_r)\,(ER^r-D^r_c)$, so $\delta^r_W\ge-\varepsilon\iff Y_{rc}\le Y^{ER}_{rc}(\varepsilon):=\big(ER^r+\varepsilon/\operatorname{Var}_p(X_r)+Y_c/(N-A_r)\big)/\kappa_r$.
(e) So every one-vs-rest sign (and each CR branch) on cell $(r,c)$ over a unit set $\mathcal S$ is a **total row** on $Y_{rc}(\mathcal S)$. Pairwise $\delta_B(r,s)$ and, for $R\ge3$, W2's $\Delta_r$ are rows with unequal coefficients across units and cells.

*Plain reading.* "Group $r$ chose $c$ more often than the rest", read as a between or a within association over the actual voters, is a floor or a ceiling on one number: how many group-$r$ people chose $c$ in total. A pairwise comparison with one other group, by contrast, weighs precincts by their composition.

*Proof.*
1. (a): Def of $B$; $B_{\neg r,c}=\sum_{s\ne r}Y_{sc}/(N-A_r)=(Y_c-Y_{rc})/(N-A_r)$. [Def; margins]
2. (b): $p_iX_{ir}X_{is}\beta_{irc}=\frac{N_i}N\frac{a_{ir}}{N_i}\frac{a_{is}}{N_i}\frac{n_{irc}}{a_{ir}}=\frac{a_{is}n_{irc}}{NN_i}$, and symmetrically. [2c Def 1]
3. (c): sum (b) over $s\ne r$ and use $\sum_{s\ne r}a_{is}=N_i-a_{ir}$, $\sum_{s\ne r}n_{isc}=b_{ic}-n_{irc}$: the numerator is $(N_i-a_{ir})n_{irc}-a_{ir}(b_{ic}-n_{irc})=N_in_{irc}-a_{ir}b_{ic}$. Then 2c M1(b) gives $NM_r$. [alg; 2c:M1(b)]
4. (d): 2c M1(c) for the collapsed population, with (a). [2c:M1(c)]
5. (e): (c), (d) are rows on $Y_{rc}(\mathcal S)$ with coefficient 1 on each $n_{irc}$, $i\in\mathcal S$; (b) has coefficients $a_{is}/N_i$ varying with $i$. [alg] ∎

*Check.* S1: 200 random R×C populations, the three count forms against 2c Definition 1 computed from shares and rates, max error $1\times10^{-16}$. **Status: P; proof review pending.** (Each $\delta^r_W$ row needs $\operatorname{Var}_p(X_r)>0$ over its own unit set.)

### 3.2 2c M2 and M3 in count coordinates

**Proposition C2.** Fix $(r,c)$, the one-vs-rest pairing (W1), and **one unit set $\mathcal S$**: every sign row is over $\mathcal S$, $Y_{rc}$ means $Y_{rc}(\mathcal S)$, and $N,A_r,Y_c,\operatorname{Var}_p(X_r)$ are computed over $\mathcal S$, with $\operatorname{Var}_p(X_r)>0$ and $0<A_r<N$. (With sign rows over several sets, or a target over another set, the closed forms below do not apply; R19(b) items 1–2 still give integer endpoints by the rounded LP when the rows are in one class. Review check RC5b.)
(a) **Continuous (2c M2, M3 restated).** Any system of one-vs-rest sign rows on $(r,c)$ (any directions, slacks) cuts $Y_{rc}$ to an interval $J$: $\delta^r_B$ rows give half-lines at $Y^{NM}_{rc}$ (shifted by $N\varepsilon$), $\delta^r_W$ rows half-lines at $Y^{ER}_{rc}(\varepsilon)$; CR (without slack) is $J_+=[Y^{NM},Y^{ER}]\cup J_-=[Y^{ER},Y^{NM}]$, one of which is a point or empty; with slack on the CR rows both branches can be nondegenerate intervals. With no unit rows, the range of $Y_{rc}$ is the DD interval $[Y^-,Y^+]$ (Lemma 2.1(a)) and 2c M3(a)–(d) hold verbatim with $D=D^r_c$, an increasing affine function of $Y_{rc}$: in particular $Y^{NM}_{rc}\in[Y^-,Y^+]$, so a between sign is never refuted in this continuous, margins-only problem, and a within sign is refuted iff $Y^{ER}<Y^-$ (resp. $>Y^+$).
(b) **Integer.** With integer margins and unit rows of R19(b)'s class (i) on the same cell (boxes, known cells, totals on $(r,c)$ over at most two laminar families, with $\mathcal S$ addable to one of them), the integer values of $Y_{rc}$ under the sign rows are **all integers** in $[\max(Y^-_{\rm rd},\lceil\text{lower thresholds}\rceil),\ \min(Y^+_{\rm rd},\lfloor\text{upper thresholds}\rfloor)]$, where $Y^\pm_{\rm rd}$ are the integer extremes without the sign rows. So the integer gap set is $\{D^r_c(k)\}$ for the integers $k$ in that interval (step $\kappa_r=1/A_r+1/(N-A_r)$). Without unit rows, or with integer boxes, this is the lattice inside the continuous interval; with fractional boxes on several units it can be strictly smaller (T5: continuous range $[0.4,2]$, integer values $\{2\}$; check S5). For any other linear target under the same rows, integer endpoints = rounded LP (R19(b)); under CR, two rounded LPs.
(c) **Integer refutation.** Without unit rows: a between sign always has an integer table ($\lceil Y^{NM}\rceil\le Y^+$, since $Y^+$ is an integer $\ge Y^{NM}$); a within sign is integer-infeasible iff $\lfloor Y^{ER}\rfloor<Y^-$ (resp. $\lceil Y^{ER}\rceil>Y^+$); CR is integer-infeasible iff no integer lies in its continuous interval of $Y_{rc}$. Since $Y^\pm$ are integers, a single sign is integer-refuted iff it is continuously refuted ($\lfloor Y^{ER}\rfloor<Y^-\iff Y^{ER}<Y^-$); only CR, or several signs together, can be integer-infeasible while continuously feasible. With unit rows a between sign can be refuted even continuously.

*Plain reading.* For "group $r$ voted $c$ more than the rest" (between or within, or both as contextual reinforcement), the belief acts only on the total number of group-$r$ people choosing $c$: it sets a floor or a ceiling on that one number. With whole people, the floor is rounded up and the ceiling down, and every whole number in between is possible. Known counts and cell limits do not change this.

*Proof.*
1. (a): by Lemma C1(c)–(d) each row is $Y_{rc}\ge$ or $\le$ a threshold, so the system is $F_0\cap\{Y_{rc}\in J\}$. CR's branches are the intersections of the two half-lines with matching signs. [L C1; 2c Def 2]
2. With no unit rows, $Y_{rc}$ ranges over $[Y^-,Y^+]$ (Lemma 2.1(a)). The projection of the R×C problem onto $Y_{rc}$ equals that of the collapsed two-row problem (both are the DD sums), and $D^r_c$ is the collapsed gap; so 2c M3, whose proof uses only this range, the identities of M1(c) and the equal-rates witness, applies. The equal-rates table $n_{irc}=a_{ir}b_{ic}/N_i$ lies in $\mathcal P$ and has $Y_{rc}=Y^{NM}_{rc}$, which proves $Y^{NM}\in[Y^-,Y^+]$ directly. [L 2.1; 2c:M1(c), M3; NW]
3. (b): the sign rows are totals on $(r,c)$ over $\mathcal S$ (L C1(e)), so the system is in class (i). R19(b) items 1–3 give integrality, integer endpoints, and (for $Y_{rc}(\mathcal S)$, addable by hypothesis) every integer between. The gap set is the image under the affine map of Lemma C1(a); it is the full lattice in the continuous interval only when $Y^\pm_{\rm rd}$ are the rounded continuous extremes, i.e. without fractional boxes on several units. [L C1; R19(b)]
4. (c): $Y^\pm$ are integers (Lemma 2.1(a) with integer margins), and $Y^{NM}\le Y^+$ (step 2), so $\lceil Y^{NM}\rceil\le Y^+$; the rest is (b) with one half-line or with $J_\pm$. With unit rows the range of $Y_{rc}$ shrinks to $[Y^-_{\rm rd},Y^+_{\rm rd}]$, which need not contain $Y^{NM}$, so a between sign can be refuted (S4). [2; 3] ∎

*Checks.* S2 **[F]**: 500 CR branches in random R×C designs (up to 3×3) with a box and a known cell: rounded LP = MILP at 1000 endpoints, 0 fractional vertices (fails on NS's statement that pooled realised signs and CR break TU, §9 item 1). S3: 2c Remark M2.1's two referee examples reproduced in counts (continuous CR image $\{0\}$ and integer CR empty; between-sign integer lower gap $2/15$). S3b: closed forms against the LP in 200 random 2×2 populations, and the integer lower end $\max(Y^-,\lceil Y^{NM}\rceil)$. S4: a known cell and a box make a between sign infeasible. S5: the T5 design (review check RC5): integer values $\{2\}$ while the continuous range $[0.4,2]$ contains 1. **Status: P; proof review pending.**

*Relation to 2c Remark M2.1.* Its proposition (2×2, margins only, sign rows on one column: integer gap set = DD ∩ lattice) is the special case $R=2$, no unit rows, of (b). Prop. C2(b) extends it to R×C one-vs-rest rows, to box/known-cell unit restrictions and published totals on the same cell, and (Corollary 2.2, class (iii)) to several columns jointly when $R=2$. 2c is not edited (§9 item 4).

### 3.3 Realised signs versus model signs

The rows above are beliefs about the realised table ($\delta$, `object = "actual"`). Model signs ($\delta^*$, 2c M5, rows on $k$, `object = "model"`) are different beliefs; imposing both is two assumptions (NS Part 3A footnote 2; 2c Remark M5.1). They never enter one LP (§4.4).

### 3.4 What is left of R15, now a computed check

R15 (PD §8) asked for integer versions of N1 and Δ4. After R19(b) and Prop. C2(b) the following are **proved**: margins, known cells, boxes (one-vs-rest unit orderings, odds-ratio dial), published totals in classes (i)–(iii), one-vs-rest signs and CR on one cell over at most two laminar families of unit sets, and for $R=2$ signs in several columns jointly.

**Not proved, decided by MILP** (exact integer answers computed per branch; no theorem bounds the integer–continuous gap, which can be large in small tables, C1 and N10c):
* pairwise $\delta_B(r,s)$ and W2 $\Delta_r$ for $R\ge3$;
* one-vs-rest signs or totals on cells in different rows **and** different columns over overlapping units ($R\ge3$ and $C\ge3$), signs on several rows of one column when a cell carries two laminar families, and any mixture of the classes;
* pairwise unit orderings for $R\ge3$;
* three or more crossing families of totals; weighted average-gap rows;
* Δ4 with integer cells (step 2).

The package can solve each of these exactly as a MILP (`integer = "auto"`), so R15 becomes a computed check: the register entry should read "integer answers computed exactly by MILP outside the TU classes; theorem for the classes of R19(b)".

---

## 4. R27: the cross-layer combination rule

### 4.1 Setting

* **Count layer** (this note): the set $\mathcal C_\varphi=\varphi(F)$ (or $\varphi(F_+)\cup\varphi(F_-)$; or integer endpoints). It is deterministic given the margins and rows.
  * Count rows are of two kinds (review issue 2). **Identity rows** have the realised value itself as their right-hand side: the margins, a published total, an audited known cell. They hold in every replication. **Belief rows** are statements about the realised table that may or may not hold in a given replication: a realised ordering, a realised sign, a box, a stated floor on a count, an odds-ratio dial.
  * Under $H_M$ the realised table $n^\circ$ is random given the design (P3-S1a). So "the count rows hold" is the **event** $\{n^\circ\in F\}$, not a hypothesis that holds or fails once and for all.
  * $H_C$ (used only for identity rows): the stated margins are the true margins and $P_{H_M}(n^\circ\in F)=1$. This holds when every count row is an identity row with exact data.
* **Model layer**: an LCE rung with optional model rows on $k$ (M-SIGN, contrasts, ...), under P3-S1a and P3-S2.
  * $H_M$: the rung and its rows hold, with the conditions of the theorem behind the interval used.
  * Under $H_M$, $B_\varphi:=\varphi(n^\circ)=\tilde B_\varphi+E_\varphi$, where $\tilde B_\varphi=\varphi_0+\sum v_{irc}a_{ir}m_{rc}(V_i)$ is affine in $k$ with known coefficients (2a Prop. 7; 2c M5(a) stacked form $m_i=m_{\rm id}(V_i)+L_ik$) and $E_\varphi=\sum v_{irc}a_{ir}e_{irc}$, $e=\beta-m$.
* An **actual-target model interval at level $1-\alpha$** is a data-dependent interval $\mathrm{CI}_\varphi$ with $\liminf_nP(B_\varphi\in\mathrm{CI}_\varphi)\ge1-\alpha$ under $H_M$ (or the finite-sample version).

### 4.2 The residual allowance for any linear count target

P3.5(a) states the allowance for a single cell $B_{rc}$. Gaps, merged groups and sums over columns need:

**Lemma 4.1.** Let the coefficients $v_{irc}$ be **functions of the design only** ($\{(N_i,X_i,D_i)\}$, P3-S2); $\varphi_0$ may depend on observed data. Let $\rho_i=\sum_ra_{ir}\big(\max_cv_{irc}-\min_cv_{irc}\big)$ and $S_\varphi=\tfrac12\big(\sum_i\rho_i^2\big)^{1/2}$. Under P3-S1a:
(a) $\operatorname{Var}(E_\varphi\mid\text{design})\le S_\varphi^2$, and if $\max_i\rho_i^2/\sum_i\rho_i^2\to0$ then $\limsup_nP(\lvert E_\varphi\rvert>z_{\alpha_E}S_\varphi)\le\alpha_E$.
(b) **Finite sample (Hoeffding).** For every $n$, $P\big(\lvert E_\varphi\rvert>t_E(\varphi)\mid\text{design}\big)\le\alpha_E$ with $t_E(\varphi)=\big(\tfrac12\log(2/\alpha_E)\sum_i\rho_i^2\big)^{1/2}$, which is P3.8's $t_E$ with $\rho_i$ for $w_{ir}$.

For $\varphi=B_{rc}$, $\rho_i=w_{ir}$ and $S_\varphi=S_1$ of P3.5(a). For a one-column gap $B_{rc}-B_{sc}$, $\rho_i=w_{ir}+w_{is}$. Targets whose coefficients use outcome margins (e.g. a share among the $c$-voters, $Y_{rc}/Y_c$) are **not** covered: then $E[E_\varphi\mid\text{design}]\ne0$ in general.

*Plain reading.* The gap between the actual and the model-implied value of any count target is a sum of independent, bounded unit terms; the usual Popoviciu allowance applies with each unit's range in place of its weight.

*Proof.*
1. $E_\varphi=\sum_iZ_i-E[Z_i\mid V_i]$ with $Z_i=\sum_ra_{ir}\sum_cv_{irc}\beta_{irc}$, since $E[\beta_i\mid V_i]=m(V_i)$. [P3-S1a; Def]
2. Each $\beta_{ir\cdot}$ lies in the simplex, so $\sum_cv_{irc}\beta_{irc}\in[\min_cv_{irc},\max_cv_{irc}]$, and $Z_i$ lies in an interval of length $\rho_i$, however the groups' rates in unit $i$ depend on each other. [simplex; lin]
3. $\operatorname{Var}(Z_i\mid V_i)\le\rho_i^2/4$, and the $Z_i$ are independent given the design, so $\operatorname{Var}(E_\varphi\mid\text{design})\le S_\varphi^2$. [Pop; P3-S1a]
4. The tail bound in (a): P3.5(a) step 5 with $\rho_i$ for $w_i$ and $S_\varphi$ for $S_1$ (Chebyshev on subsequences with $s_n/S_\varphi\to0$, Lindeberg otherwise, using $\lvert Z_i-EZ_i\rvert\le\rho_i$). [Chebyshev; Lyap]
5. (b): the $Z_i$ are independent given the design and each lies in an interval of length $\rho_i$ (steps 2–3, which use that $v$ is design-measurable), so Hoeffding's inequality applied to $\pm E_\varphi$ gives $P(\lvert E_\varphi\rvert>t)\le2\exp(-2t^2/\sum_i\rho_i^2)$, which equals $\alpha_E$ at $t=t_E(\varphi)$. [Hoeff; P3-S1a] ∎

**Remark 4.1a (the representation of a target; review issue 3).** $\tilde B_\varphi$, $E_\varphi$ and $S_\varphi$ refer to the stated representation $(\varphi_0,v)$. Adding $\lambda_{ic}\big(\sum_rn_{irc}-b_{ic}\big)$ with design-measurable $\lambda$ leaves $\varphi$ unchanged on $\mathcal P$ but changes $v$ and $\varphi_0$ (which then uses the observed $b_{ic}$), and hence $E_\varphi$ and $S_\varphi$. Example: the one-vs-rest gap written as $\kappa_rY_{rc}-Y_c/(N-A_r)$ or over cells as $Y_{rc}/A_r-\sum_{s\ne r}Y_{sc}/(N-A_r)$; the two $E_\varphi$ differ by $(Y_c-E[Y_c\mid\text{design}])/(N-A_r)$ (review check RC7: $S_\varphi$ 0.0465 vs 0.0555). Each representation is valid with its own allowance. The package may use any **fixed** representation chosen from the design alone (for example the one minimising $S_\varphi$), never one chosen after seeing the outcomes.

*Check.* G1: 200 designs with random dependent joint laws of the two group rates in each unit: $\operatorname{Var}/S_\varphi^2\le0.21$; perfectly opposite rates attain 1. G2: the Hoeffding allowance of (b) for a gap target, simulated tail probability at most $\alpha_E$. **Status: P; proof review pending** (P3.5(a) steps 4–5 and P3.8's Hoeffding step with a different range). So P3.5(a), (b2)(ii) and P3.8 give actual-target intervals for every linear count target with design-measurable coefficients, by replacing $S_1$ with $S_\varphi$ (resp. $t_E$ with $t_E(\varphi)$). The Bhatia–Davis version P3.5(c1) is proved only for single cells and must not be used for multi-cell targets.

### 4.3 Which model intervals qualify, and their status

| Model-layer interval for $B_\varphi$ | Scope | Status inherited by R27 |
|---|---|---|
| P3.5(a) conservative, allowance $S_\varphi$ | pooled M0, any regime; M1/M2 in R1 (fixed districts); with pooled M-SIGN rows: R10(a) | **P** (asymptotic) |
| P3.5(b2)(ii) projection, Popoviciu allowance $S_\varphi$ | same; $q=\dim\beta$ or an exact reduction | **P** (asymptotic) |
| projection with the Bhatia–Davis allowance P3.5(c1) | single-cell targets, $\mathrm{CI}_v$ on the rung alone, (c1)'s Lindeberg condition | **P** (asymptotic) |
| P3.8 finite-sample interval, allowance $t_E(\varphi)$ (Lemma 4.1(b)) | any rung, any regime; design-measurable targets | **P** (finite sample; currently wide) |
| CLR as implemented (`ci_actual`; P3.5(b3), (c2)) | default interval | **algorithmic** (conjecture P3.5(b3)) |
| grouped rungs (M1, M1(δ), M2, L4) **with** model sign rows | — | **no coverage theorem** (R16): the intersection is experimental |
| grouped rungs in regime R2 (many small districts) | — | simulation finding only (P3 §0.5) |

### 4.4 The rule

**Theorem R27.** Let $\mathcal C_\varphi$ be the count-layer set and $\mathrm{CI}_\varphi$ an actual-target model interval at level $1-\alpha$ for the same target $\varphi$.
Let the margins be exact. Under $H_M$ alone, for every sample size,
$$P(B_\varphi\notin\mathcal C_\varphi\cap\mathrm{CI}_\varphi)\le P(B_\varphi\notin\mathrm{CI}_\varphi)+P(n^\circ\notin F),\qquad P\big(\{\mathcal C_\varphi\cap\mathrm{CI}_\varphi=\emptyset\}\cap\{n^\circ\in F\}\big)\le P(B_\varphi\notin\mathrm{CI}_\varphi).\tag{$*$}$$
(i) **Coverage.** With identity rows only ($H_C$: $P(n^\circ\notin F)=0$), $P(B_\varphi\in\mathcal C_\varphi\cap\mathrm{CI}_\varphi)=P(B_\varphi\in\mathrm{CI}_\varphi)$, hence $\liminf P(B_\varphi\in\mathcal C_\varphi\cap\mathrm{CI}_\varphi)\ge1-\alpha$ (exactly $\ge1-\alpha$ for a finite-sample interval): the intersection covers with at least the model interval's coverage. With belief rows, $(*)$ gives coverage at least $1-\alpha-P(n^\circ\notin F)$, and $P(n^\circ\notin F)$ is unknown; the level $1-\alpha$ holds only if the belief rows are taken as certain.
(ii) **Disjointness test.** Reject when $\mathcal C_\varphi\cap\mathrm{CI}_\varphi=\emptyset$.
* With identity rows only (in particular margins only, $\mathcal C_\varphi$ = the DD range), $\limsup P(\text{reject})\le\alpha$ under $H_M$: a level-$\alpha$ test of **the model and its rows**.
* With belief rows, $(*)$ bounds the probability of **rejecting while the belief rows are true** by α. It is not a test with level α conditional on the rows: given a row that is rarely true under $H_M$, the conditional rejection rate can be near 1 (review check RC8; X6). A rejection is evidence against the model and its rows only with the belief rows maintained, and a false belief row produces rejections too (X4).
(iii) **Not jointly sharp.**
* (iii-a) **Compatible values.** Suppose every $m(V_i)$, at some $k$ in the model set, lies in the relative interior of $\mathfrak B$ (or respects the same structural zeros as $F$). Then the set of values of $B_\varphi$ compatible with the observed margins, the count rows and $H_M$ is $\mathcal C_\varphi$ itself: $H_M$ excludes no table of $F$. The information $H_M$ adds is probabilistic and enters only through $\mathrm{CI}_\varphi$. (This is a statement about the realised value, not an identified set in the population-law sense, which is not defined for one realisation; review issue 10.)
* (iii-b) The intersection is not in general the smallest valid confidence set. $\mathrm{CI}_\varphi$ is computed without the count rows, which constrain the realised residuals. Using the count information inside the model layer can give a strictly narrower set at the same level (Remark 4.4; X5).

*Plain reading.* When the count layer uses only facts (the totals, published counts, audited cells), its range always contains the actual rate, so cutting the model's interval down to it loses no coverage, and a non-overlap is a level-α rejection of the model. When it also uses beliefs about the actual voters, those beliefs can be wrong in this electorate, and the chance that they are wrong adds to the error rate; a non-overlap then says that the model or one of the beliefs is wrong. The overlap is a safe answer, not the tightest possible one.

*Proof.*
1. On the event $\{n^\circ\in F\}$, $B_\varphi=\varphi(n^\circ)\in\mathcal C_\varphi$ (R19(a) item 2; for integer endpoints, R19(b), since $n^\circ$ is an integer table of $F$). [R19(a); Def]
2. Hence $\{B_\varphi\notin\mathcal C_\varphi\cap\mathrm{CI}_\varphi\}\subseteq\{B_\varphi\notin\mathrm{CI}_\varphi\}\cup\{n^\circ\notin F\}$, and $\{\mathcal C_\varphi\cap\mathrm{CI}_\varphi=\emptyset\}\cap\{n^\circ\in F\}\subseteq\{B_\varphi\notin\mathrm{CI}_\varphi\}$. This gives $(*)$. [1; union]
3. (i), (ii) with identity rows: $P(n^\circ\notin F)=0$, so the events of step 2 coincide with $\{B_\varphi\notin\mathrm{CI}_\varphi\}$ up to a null set; take $\liminf$ / $\limsup$ under $H_M$. With belief rows, read off $(*)$. [2; Def CI]
4. (iii-a): a table $n\in F$ is consistent with the margins and rows. $H_M$ restricts only the conditional means $m(V_i)$ of the unit rates (P3-S1a). If $m(V_i)$ lies in the relative interior of $\mathfrak B$, then for any $\beta_i$ in $\mathfrak B$ the two-point law putting mass $\lambda$ on $\beta_i$ and $1-\lambda$ on $(m(V_i)-\lambda\beta_i)/(1-\lambda)$ has mean $m(V_i)$ and, for small $\lambda>0$, both points in $\mathfrak B$. So no table of $F$ is excluded by the count rows together with $H_M$; only its probability is restricted. (With structural zeros in $m$, the same holds for tables respecting them.) Hence no deterministic set smaller than $\mathcal C_\varphi$ follows. [Def; P3-S1a; convexity of $\mathfrak B$]
5. (iii-b): Remark 4.4 gives a set valid at the same level by (i)-type arguments; X5 shows it strictly narrower in every one of 1000 replications of a design. ∎

**Status:** $(*)$, (i), (ii) **P; proof review pending**, with the model interval's status inherited (table §4.3); the exact level-α statements are for identity rows only. (iii-a) **P; proof review pending**. (iii-b) a computed example; no general dominance is claimed.

**Remark 4.1 (why count rows and model rows never share one LP).** An LP with count rows on $n$ and model rows on $k$, linked unit by unit through $n_{irc}=a_{ir}m_{rc}(V_i)$, would set $E=0$: it treats realised rates as the model's means, which is false for a finite electorate (PD §4.2 footnote f; 2c Remark M5.1). Under LCE no per-unit constraint links the realised counts to $k$; the only link is the target identity $B_\varphi=\tilde B_\varphi+E_\varphi$, with $E_\varphi$ controlled in probability (Lemma 4.1). That link is used once, in the intersection. PD §4.2 footnote f becomes this rule instead of a warning.

**Remark 4.2 (several targets; level splitting).**
* (i) and (ii) are per target. For several targets at once, take a union bound over targets, or a simultaneous interval. The projection event of P3.5(b2)(ii) is simultaneous for all linear functionals, but each target's allowance event is separate.
* If the count layer is itself a confidence set at level $1-\alpha_C$ (estimated margins, Step 1b), then $P(B_\varphi\notin\mathcal C_\varphi\cap\mathrm{CI}_\varphi)\le\alpha_C+\alpha_M$, and the disjointness test has level $\alpha_C+\alpha_M$. [union] This is $(*)$ with $\alpha_C$ in the role of $P(n^\circ\notin F)$: there $\alpha_C$ is controlled by construction, whereas for belief rows it is unknown.

**Remark 4.3 (the package's existing intersection with DD).** `ci_actual` is already intersected with Duncan–Davis, and an empty intersection sets `ci_rejected` (`rxc_infer.R`). Since $\mathcal C_\varphi\subseteq$ DD, intersecting first with DD and then with $\mathcal C_\varphi$ changes nothing, provided the emptiness information is kept and reported. The margins-only case of (ii) is therefore already computed; with belief rows, the report must carry (ii)'s second bullet.

**Remark 4.4 (a narrower valid set; the content of (iii-b)).** Suppose known cells fix $n_{irc}$ for the units in $\mathcal K$ (an $H_C$ row), and let $\mathcal U$ be the other units. Then $B_{rc}=\big(\sum_{\mathcal K}n_{irc}+A_r(\mathcal U)B_{rc}(\mathcal U)\big)/A_r$. Under $H_M$, **if $\mathcal K$ is chosen from the design only** (as in P3's design-based selection, Prop. P3.7(b)) and the Lindeberg condition holds over $\mathcal U$, the model interval for $B_{rc}(\mathcal U)$, with weights and allowance over $\mathcal U$ only, covers at level $1-\alpha$. Audits aimed at close or anomalous precincts (selected on outcomes) violate this (review issue 11). So the set
$$J=\Big(\textstyle\sum_{\mathcal K}n_{irc}+A_r(\mathcal U)\,\mathrm{CI}_{\mathcal U}\Big)\Big/A_r,\quad\text{intersected with }\mathcal C_\varphi,$$
is valid at level $1-\alpha$. Its allowance is $\tfrac12(\sum_{\mathcal U}w_{ir}^2)^{1/2}<S_1$. In X5 (1000 units, LCE true, known cells in a random half of them) $J\cap\mathcal C$ was narrower than $\mathcal C\cap\mathrm{CI}$ in 1000 of 1000 replications (median widths 0.211 vs 0.258), with coverage 1.000. This is a refinement for later; it is not part of the M-COUNT contract.

*Checks.* X1 **[I]** (NS's worked example: $[0.45,0.875]\cap[0.30,0.60]=[0.45,0.60]$; with $[0.20,0.40]$, disjoint). X2 (a pooled 2×2 LCE model layer with the projection interval and the Popoviciu allowance at $\alpha=.05$, 1000 replications each: coverage of the intersection equals coverage of the model interval in every replication; disjointness rejection rate 0.000 with margins only, with 90% of cells known (80 units), and with 1000 units. The interval is conservative here, so the size is far below α; the theorem gives only "≤ α"). **X3 [F]** (truth on the boundary of the LCE set, count layer pins $B$: the plug-in model set rejects in 42.5% of replications; the P3.5 interval for $B$ in 0%): the model interval must be an actual-target interval. **X4 [F]** (model true, false belief rows: rejection rate 1.000): disjointness involves the belief rows. **X6** (belief rows, review issue 2: a belief row "$B\ge t$" with an exact interval for $B$; conditional coverage given the row falls below $1-\alpha$ and the conditional rejection rate exceeds α, while both inequalities $(*)$ hold). X5 **[I]** (Remark 4.4).

---

## 5. Sharpness relative to the information used

**What the count layer uses:** the unit margins, and the rows the user states about the realised table (known cells, totals, orderings, signs, dials). Nothing else: no regression across units, no model of how rates vary with composition, no sampling model.

**Sharp for that information (R19(a), (b), (d)).**
* For continuous tables, every value in the range is attained by a table consistent with exactly that information.
* For integer data in the classes of R19(b), both endpoints are attained by integer tables; outside them, the MILP gives the attained integer endpoints.
* The range cannot be narrowed without adding information.

**Not sharp for the model's targets.**
1. *Different estimand.* $\mathcal C_\varphi$ bounds $B_\varphi$, a function of the realised table. $\tilde B_\varphi$ (for a fixed representation of $\varphi$, Remark 4.1a) is not such a function and can lie outside $\mathcal C_\varphi$. For example, with every cell known, $\mathcal C_\varphi=\{B_\varphi\}$, while $\tilde B_\varphi=B_\varphi-E_\varphi\ne B_\varphi$ in general. The count layer gives no bound on $\tilde B$ or $B^*$.
2. *Information it ignores.* A model (LCE, M-SIGN rows) is extra information about how the realised table was generated. It can narrow the answer for $B_\varphi$ in probability (R27), never deterministically (R27(iii-a)). It also gives answers for $\tilde B$ that the count layer cannot.
3. *Within the count layer, sharpness is relative to the rows actually imposed.* A belief left out (for example, a pairwise ordering that holds but was not stated) is not used. The redundancy diagnostic (PD R2) applies unchanged: a row tightens a target iff it removes the whole optimal face.

---

## 6. The interface contract (specification for the code step; no code here)

The R implementation must meet the following; the code review checks each item.

**C1. Function and placement (D2, D1).** `eiCounts()` is a separate function. `eiBounds()` embeds it and prints its actual-rate block **first**, with the margins-only (Duncan–Davis) line first in that block. `target = "model"` restores the old order.

**C2. Targets.** Any linear count functional over a known denominator: `rate(row, col, where=)`, `gap(row1, row2, col, where=)` (one-column contrast), merged groups, sums over columns, and differences over unit sets. A cross-layer (R27) interval is offered only for targets whose coefficients are functions of the design (Lemma 4.1); targets with outcome-margin denominators get the count layer only. The margins-only line is the **LP** value (Lemma 2.1). It must equal the current `dd_lower`/`dd_upper` for single cells and two-row one-column contrasts. It is narrower for other combinations, including the one-vs-rest gap $D^r_c$ with $R\ge3$, the M-SIGN default contrast (review check RC6; §9 item 3).

**C3. Rows.** Each row carries `object = "actual"`. Every row is mapped to a linear row on $n$ in the forms of §1.3:
* orderings (one-vs-rest) and odds-ratio dials as **boxes**;
* one-vs-rest signs as **totals** on $Y_{rc}(\mathcal S)$ (Lemma C1);
* pairwise and W2 signs with their unequal coefficients;
* CR as a union of two systems; slack dials as right-hand-side shifts.

The rows of 2c M2 (signs, CR) and of DiD (step 2) print "experimental" until Batch B comment 1 and Batch C resolve (NS §3.1 review gates). Margins, known cells and totals do not.

**C4. Integer modes.**
* `"auto"` (default, D6): if any margin is non-integer, no integer check, and print C7's note. Otherwise round the rows (§1.4); if the **whole set** of non-union rows (per CR branch) lies in **one** class of R19(b), report the rounded LP as integer-exact; otherwise solve the MILP for each endpoint and each CR branch. Checking each row separately is not enough (mixtures fail, review check RC4). "Two laminar families" (class (i)) is decided by 2-colouring the crossing graph of the unit sets (two sets cross if they intersect and neither contains the other): the sets split into two laminar families iff that graph is bipartite.
* **Exact arithmetic.** Thresholds (ordering boxes, $Y^{NM}$, $Y^{ER}$, slack-shifted right-hand sides) are computed and rounded in exact rational arithmetic, or rounded with a certified tolerance that is printed: a double-precision floor of an integer-valued threshold was wrong in 1,281 of 4,008 cases (review check RC9). Margins count as integers only if they are exactly integer in the input (not $N\cdot X$ recomputed in floating point).
* `"lp"`: the continuous LP only, labelled "sharp for continuous tables".
* `"exact"`: always the MILP (integer margins required; refuse, with a message, for non-integer margins).

Never report "integer" for the ceiling of an unrounded LP endpoint (T5).

**C5. Infeasibility.** Report separately:
* continuous infeasibility: the rows contradict the margins;
* integer infeasibility with continuous feasibility (N10c; crossing totals; 2c Remark M2.1). Report it as "no whole-number table satisfies these rows", not as an empty interval with no explanation.

Under CR, report which branches are nonempty and whether the union is one interval or two.

**C6. Witness.** `witness(fit)` returns, for each target and endpoint, a table from a basic optimal solution (simplex, or interior point plus crossover) or from the MILP. Tests must check that it reproduces every margin and satisfies every row and bound, and that it is integer when C4 reports integer-exact. Print it as a certificate of attainability, not as an estimate.

**C7. `print()` must state:**
* target and denominator, lower, upper, width;
* the integer status: "integer-exact (TU class)", "integer-exact (MILP)", "continuous only" or "not applicable: non-integer margins";
* the binding rows (via the R2 face test) and the rows' `object`;
* for non-integer margins: "**exact only for the stated margins, which are not integers; if they are estimates, see estimated margins (Step 1b)**".

**C8. Cross-layer block (R27).** When a model fit is present, for each target with an actual-target model interval:
* print $\mathcal C_\varphi$, the model interval (naming its type and status from §4.3), and their intersection;
* classify every count row as an **identity row** (margins, published totals, audited known cells) or a **belief row** (orderings, signs, CR, boxes, stated floors, odds-ratio dials);
* with identity rows only, print the intersection as covering at the model interval's level and, if disjoint, "disjoint: rejects (model + its rows) at level α";
* with belief rows, print instead: "coverage at least 1 − α − P(belief rows false), which is unknown; level 1 − α only if the belief rows are taken as certain" and, if disjoint, "disjoint: the model (+ its rows) or a belief row is wrong; P(reject while the belief rows are true) ≤ α". Never print "rejects at level α given the count rows" (it reads as conditioning, which is false: R27(ii), review check RC8);
* never put count rows and model rows into one LP;
* refuse to intersect with a model **estimate** (plug-in set; X3 [F]) or with an interval for $\tilde B$ without the residual allowance (an interval for a different target, P3 §0.2);
* for targets other than single cells, use the allowance $S_\varphi$ (or $t_E(\varphi)$ for P3.8) of Lemma 4.1, with a representation of the target fixed from the design alone (Remark 4.1a); the Bhatia–Davis allowance is forbidden for multi-cell targets;
* keep the emptiness information when `ci_actual` has already been intersected with DD (Remark 4.3).

**C9. One vocabulary.** "More supportive" and the like exist once, with `object = c("actual","model")`. The same restriction may not be imposed under one label in both layers. `print()` lists realised and model beliefs separately (NS Part 3A footnote 2).

**C10. Tests (in addition to NS §3.1's list).**
* D1/D1b against the package formula;
* brute force in tiny designs;
* T1–T3 random TU-class instances (LP = MILP);
* T5;
* N10c reported as integer-infeasible;
* C1's four row types routed to MILP;
* every witness feasible;
* adding a row never widens a range;
* the cross-layer block never reports a non-empty intersection outside $\mathcal C_\varphi$;
* "auto" sends a mixture of classes (RC4) to the MILP;
* thresholds that are exact integers in rational arithmetic are rounded to themselves (RC9);
* the belief-row wording of C8 appears whenever a belief row is present.

---

## 7. Step 1b (estimated margins): pointer only

Step 1b (NS §3.1b) comes next, with proof items **R25** (coverage of the count layer over a margin set $\mathcal V$ with $P(V\in\mathcal V)\ge1-\alpha$; joint sets; a calibration slack) and **R26** (the Berkson correction in the LCE layer). It is not written here. M-COUNT must leave these hooks:
1. **Margins as variables.** The LP builder accepts unit group counts as variables $V_i$ in a polyhedral set (per-unit intervals, or joint rows such as $\sum_iV_i\in[\cdot,\cdot]$), with the cell rows $n_{iVc}\le V_i$ and $n_{iVc}\ge Y_{ic}-(N_i-V_i)$. Fixed margins are the special case of a singleton set.
2. **Linear-fractional targets.** When the denominator is $\sum_iV_i$ (unknown), the target is linear-fractional: Charnes–Cooper, or profiling over integer $\sum V$. The target interface must allow a variable denominator.
3. **Integer mode on latent counts.** The integer check applies to the true (latent) $V_i$ and $n$, not to expected counts. `"auto"` must be able to switch it on when $V$ is a variable integer vector even though the expected margins are non-integer.
4. **Level bookkeeping.** `eiCounts()` carries a level (1 for exact margins, $1-\alpha_C$ for margin sets), and R27's block splits α as in Remark 4.2.
5. **Printing.** C7's "exact only for the stated margins" line is replaced by the R25 coverage statement once margin sets are supplied, and the matching assumption (nondifferential matching) is printed.

---

## 8. Summary of results and register entries

| Result | Statement | Status | Checks |
|---|---|---|---|
| R19(a) | LP over $\prod_i\mathcal P_i\cap$ rows; valid; sharp for continuous tables given the information | P; review pending | D1, D1b, E1, E2 |
| Lemma 2.1 | DD = margins-only LP for single cells and one-column contrasts; cellwise DD outer otherwise (including $D^r_c$, $R\ge3$) | P; review pending | D1, D1b, D2, D3 [F] |
| R19(b) | integrality for classes (i)–(iii) (whole row set in one class); TU of the original matrix for (i)–(ii), of the reduced matrix for (iii); integer endpoints = rounded LP; witnesses integer; all intermediate integers for addable aggregates | P; review pending | B1, B2, T1–T4, T5 [F], T6 [F], N10a, N10b |
| R19(c) | outside the classes: LP outer, MILP exact; fractional vertices, gaps beyond rounding and integer infeasibility occur | P; review pending | C1, N10c |
| R19(d) | witness = basic optimal solution; integer under (b) | P; review pending | W1, E1 |
| R19(e) | non-integer margins: exact for the stated margins; integer check not applicable | P; review pending | E5 |
| Lemma C1 | realised associations in counts; one-vs-rest rows are totals on $Y_{rc}(\mathcal S)$ | P; review pending | S1 |
| Prop. C2 | one sign set $\mathcal S$: 2c M2/M3 in counts; integer values of $Y_{rc}(\mathcal S)$ = an integer interval (the full lattice only without fractional boxes on several units); integer refutation | P; review pending | S2 [F], S3, S3b, S4, S5 [F] |
| Lemma 4.1 | residual allowances $S_\varphi$ (asymptotic) and $t_E(\varphi)$ (Hoeffding) for linear count targets with design-measurable coefficients; representation-dependent | P; review pending | G1, G2 |
| R27 $(*)$, (i), (ii) | union bounds $(*)$; with identity rows only, coverage at the model interval's level and a level-α test of the model; with belief rows, coverage $\ge1-\alpha-P(n^\circ\notin F)$ and P(reject ∧ rows true) ≤ α | P; review pending; status inherited from the model interval | X1, X2, X3 [F], X4 [F], X6 [F] |
| R27(iii) | compatible values = $\mathcal C_\varphi$ (P; review pending); not minimal (computed example, design-based $\mathcal K$) | P / example | X5 |
| R15 remainder | rows outside the classes | computed check (MILP), no theorem | C1 |

Register rows added to PD §8: **R19** and **R27**, both "P; proof review pending" (R27 with the identity/belief-row restriction).

---

## 9. Issues found in other documents (not edited)

1. **NS §3.1 (proof items; "Realised signs and CR") and Part 3A**, and the comment in `verify_next_steps.py` N10c. NS says pooled realised signs and CR couple cells across units, so total unimodularity fails and "realised CR is two MILPs whenever integrality matters". This is false for the **one-vs-rest (W1)** signs and CR, which are the defaults. By Lemma C1 they are rows on the single aggregate $Y_{rc}(\mathcal S)$, i.e. totals, and they lie in TU class (i). They are integer-exact by the rounded LP (Prop. C2(b); check S2 [F]). MILPs are needed for pairwise $\delta_B(r,s)$, W2, rows on cells in different rows and columns ($R\ge3$ and $C\ge3$), and the other rows of §3.4. Relatedly, the N10c comment says its two weighted rows mimic "the realised $\delta_B$ and $\delta_W$ rows of 2c M1 [which] weight units differently". The realised W1 rows do not weight units differently. N10c's rows are weighted average-gap rows (with weights $w=1$ and $w=N-V$, this note's reproduction). N10c remains a valid example of coupled rows.
2. **NS Part 3A, "The main combinations you would use", first bullet:** "It does not reject the ordering, which is a belief about the actual table." Disjointness of the count set (built with the ordering) and the model interval tests the **conjunction** of the ordering and the model. A false ordering can cause it (R27(ii); X4 [F]). Only with the margins-only count layer is the rejection one of the model alone. The cross-layer rule paragraph ("Disjointness is a level-α test of the model and its rows") needs the same qualification: it holds as stated with identity rows only (margins, published totals, audited cells); with belief rows it bounds the probability of rejecting while the rows are true, i.e. it holds "with the count rows maintained" (R27 $(*)$; not "given the count rows", which reads as conditioning).
3. **Package, `eiBounds()` `duncan_davis` column** (code, not a theory document). For a combination target it is computed by cellwise interval arithmetic (`.rxc_interval_lin` in `rxc_infer.R`). That is valid but outer, and can be strictly wider than the exact margins-only range (Lemma 2.1(c); D3 [F]). This includes the one-vs-rest gap $D^r_c$ for $R\ge3$, the M-SIGN default contrast, whose exact range is the single-cell DD range of $Y_{rc}$ mapped by Lemma C1(a): strictly wider in all 200 random designs of review check RC6, by up to 0.87. NS §3.1's test "Reproduce `eiBounds()`'s `dd_lower`/`dd_upper` to 1e-12 on every reference dataset" holds for single cells and one-column contrasts only. For other combinations, `eiCounts()` will be narrower, correctly. This matters for decision D1, which prints the DD line first.
4. **2c Remark M2.1, "Open (proof item R15)":** "With unit restrictions $U_i$, joint columns, or $R\ge3$, the integer–continuous gap is not bounded by this argument." This is now partly settled. With box and known-cell unit restrictions and one-vs-rest signs in R×C, and for $R=2$ with joint columns, the integer endpoints are those of the LP after rounding the rows (Prop. C2(b), Corollary 2.2; an endpoint statement, not a lattice description). The remainder is listed in §3.4. Not an error: a scope update for 2c and for PD §8's R15 row, which this note does not edit.
5. **NS §3.1, R19(b) wording:** "For margins, known cells and per-unit box rows the constraint matrix is totally unimodular, so the LP endpoints of integer-coefficient count targets are integers." This holds only when the box bounds (and any totals) are integers. The ordering boxes $a_{ir}b_{ic}/N_i$ are fractional in general (N10b), so the rows must be rounded first, and rounding the endpoint is not enough (T5 [F]). TU is a property of the matrix; integrality of the endpoints needs integer right-hand sides as well.
6. **PD §4.2 footnote f** ("the target must be the actual rate with P3's allowance, or ... the two sets refer to different objects"). This is consistent with R27, which turns it into a rule (Remark 4.1). No error.
7. **P3.5(a)** states the residual allowance for a single cell $B_{rc}$. Count targets with several cells (gaps, merged groups) need Lemma 4.1. The package's actual-rate intervals for combination targets should be checked against $S_\varphi$ at the code review. Likewise P3.8's Hoeffding allowance $t_E$ is stated for single cells; Lemma 4.1(b) extends it. No error in P3, which never claims more.

---

## 10. Changes after the proof review

The proof review (`reviews/PHASE6a_proof_review.md`, branch `claude/mcount-proof-review` at b2dd747, merged here) accepted the core and listed 15 issues. All are fixed in this version. New checks: T6, S5, X6, G2 (36/36 PASS).

| Issue | Severity | Fix in this note |
|---|---|---|
| 1 class (iii) matrix not TU | error (statement) | R19(b) now concludes integrality; TU of the original matrix claimed for (i)–(ii) only; (iii) via the reduced system after $n_{i2c}=b_{ic}-y_{ic}$. Proof's first paragraph and step 5 say so. Check T6 [F] (minor 2; rounded LP = MILP) |
| 2 belief rows are events | gap | §4.1 splits identity and belief rows, $H_C$ = $P_{H_M}(n^\circ\in F)=1$; R27 states the union bounds $(*)$; (i)–(ii) exact for identity rows; with belief rows, coverage $\ge1-\alpha-P(n^\circ\notin F)$ and P(reject ∧ rows true) ≤ α; plain reading, proof, C8 wording, §9 item 2 and the register row changed. Check X6 [F] |
| 3 design-measurable coefficients; representation | gap | Lemma 4.1 hypothesis ($v$ functions of the design; $\varphi_0$ may use data); outcome-margin denominators excluded; Remark 4.1a on representations (fixed from the design alone); §1.2, §5 item 1, C2, C8 |
| 4 P3.8 row | gap | Lemma 4.1(b), Hoeffding $t_E(\varphi)$; §4.3 row restated. Check G2 |
| 5 C2(b) lattice sentence | error (one sentence) | integer gap set = image of the integer interval; the full lattice only without fractional boxes on several units. Check S5 [F] |
| 6 one sign set; CR slack | gap | C2 preamble fixes one set $\mathcal S$ and target $Y_{rc}(\mathcal S)$; several sets → R19(b) items 1–2 only; CR clause "without slack" |
| 7 $R\ge3$ and $C\ge3$ | wording | R19(c) example, §3.4, §9 item 1, register row |
| 8 class (ii) bullet; one class for the whole set | wording | Corollary 2.2 bullet and proviso; R19(b) statement; C4 tests the whole row set (crossing-graph 2-colouring); C10 test |
| 9 $D^r_c$ | wording | Lemma 2.1(c), C2, §9 item 3 |
| 10 "identified set" | wording | (iii-a) restated as compatible values of the realised $B_\varphi$, interiority for every unit |
| 11 design-based $\mathcal K$ | gap (minor) | Remark 4.4 hypothesis (P3.7(b)-type selection, Lindeberg over $\mathcal U$) |
| 12 exact arithmetic | wording | C4: rational arithmetic or certified tolerance; margins integer only if exactly integer; C10 test |
| 13 C2(c) | wording | single signs integer-refuted iff continuously refuted; only CR or several signs differ |
| 14 degenerate units | wording | $N_i>0$; odds-ratio convention; $\operatorname{Var}_p(X_r)>0$ over each sign set (Corollary 2.2, Lemma C1) |
| 15 $\delta_B(r,s)$ example | wording | stated as a fractional vertex only (integer endpoint = rounded LP); gaps beyond rounding rest on the other examples; a 20,000-design search found no $\delta_B$ instance beyond rounding |

