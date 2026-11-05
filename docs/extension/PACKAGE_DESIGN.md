# Package design: a modular bounding engine for ecological inference

*Design plan for a general EI bounding package with interchangeable assumption modules. It consolidates the work so far:*
* *eiPartialID 0.2.0 (PHASE1–4);*
* *the Monotone-EI and DiD cycle (PHASE1b, PHASE2c, PHASE2d);*
* *the new-avenues review (NEW_AVENUES);*
* *an outside conversation (GPT) about combining JKST with Monotone EI.*

*Nothing in the R package changes in this pass.*

*New numerical checks are in `verify_package_design.py`, with output in `verify_package_design_output.txt`, cited as [V1]–[V7]. Each result below carries a status: **proved** (with a proof or a reference), **checked** (numerically only), or **conjecture**.*

*Plain-language companion: `PLAIN_PACKAGE_DESIGN.md`.*

Notation follows PHASE2a/2c.
* $H_c=Q_c+K_c$ is the mean-level contextual-effect matrix for column $c$, with $m_{rc}(x)=(H_cx)_r$ on the composition simplex.
* $Q_c$ (the symmetric part) is identified by $E[T_c\mid X=x]=x^\top H_cx$.
* $K_c$ (antisymmetric) is unidentified, with $\dim\mathcal K=(C-1)R(R-1)/2$ (2a Lemma 1).
* In 2×2 the free coordinate is $w_1$ (the slope of row 0): $m_0=w_0+w_1x$ and $m_1=w_0+c_1-w_1+(w_1+d_1)x$.

---

## 0. Summary

1. **Architecture.** Every module turns an assumption into one of the following, acting on the same LP in $k$ (2a Prop. 7):
   * linear rows;
   * a union of LPs;
   * a grid LP;
   * an outer/inner pair.

   Each module declares: which part of the data it uses, what it restricts (a unit's rates, a conditional mean, or a link across units or periods), which target it bounds (the actual or the model-implied rate), whether aggregates can refute it, and its sensitivity dial. Combination rules then follow mechanically (§4).
2. **The GPT conversation** is mostly right, and mostly already in the repo (§2).
   * Mapping JKST bounds to $(\delta_W,\delta_B)$ adds nothing (a line segment, [V1]).
   * Sign restrictions and CR are extra rows that tighten JKST in most simulated designs (§2.2, [V2]).
   * Its "pointwise" proposal (3.B) splits into an existing module (`order`) and one genuinely new module (weighted-local contextual monotonicity, **M-LOC**).
3. **New results** (§5):
   * **(P-S, proved)** For any $R$, pairwise ordering and weighted-local monotonicity for the same pair restrict the *same* linear functional of $K$ from opposite sides. Together they confine it to an interval whose half-width is the identified directional slope of $\tau$. GPT's 3.B is the $R=2$ case.
   * **(P-R, proved for star trees; checked for others)** One pair cuts $R-1$ of the $R(R-1)/2$ unknown directions. A spanning tree of $R-1$ pairs cuts all of them [V7]. This is the precise sense in which more groups give "more opportunity". In 2×2 there is only one direction.
   * **(P-N, proved, with conjecture C1)** In 2×2 under LCE: signed slopes ⊂ weighted-local ⊂ $\delta_W\ge0$. The second inclusion holds whenever a design point $x^*$ lies in the composition hull; C1 says it always does [V3].
   * **(DiD-S, proved)** The sandwich lemma applied to *changes* gives a closed-form DiD interval under first-difference LCE. It needs two direction assumptions and no per-period linearity (§7).
4. **Two ladders.**
   * The assumption ladder (§3) has a common data-only base and three non-nested tracks: direction, smoothness and scale/dispersion.
   * The build ladder (§6) orders the work by proof readiness and gives each step a rigor checklist.
5. **Proof register** (§8): fourteen claims, with status, strategy and what each blocks.

---

## 1. Architecture

### 1.1 The core
For each column $c$, the bound on a target $\varphi=\sum a_{rc}B_{rc}$ is
$$\min/\max\ \varphi(k)\quad\text{s.t.}\quad Gk\le g\ \ (\text{rows from active modules}),$$
* possibly with auxiliary variables (unit-level cells, grid values, slack for dials);
* possibly as a union over a finite set of sign patterns (CR, unknown signs).

The package keeps this single object. Every module is a function
```
module(data, spec) -> list(rows = (G, g) | union = list(rows) | grid = ..., meta)
meta = list(source, object, target, refutable, form, sharp, dial, proof_status)
```
`eiBounds()` already has `zeros/order/sign/slope/linear/external` in this form (PHASE4). The design generalises it to a registry of modules, so `ladder()` and the comparison tables are generated from `meta`.

### 1.2 Five attributes that decide how modules combine
| Attribute | Values | Why it matters |
|---|---|---|
| Object | unit $\beta_i$ / conditional mean $M_H(x)$ / link across units or periods | Unit-level rows tighten even the sharp set. Mean-level rows bound the model-implied rate and need P3's allowance for the actual rate. |
| Target | actual rate / model-implied rate | Restrictions on different targets cannot be intersected naively (⚠ in §4.2). |
| Form | rows / union / grid / outer–inner / integer program | This sets the computational cost and the inference method (a union needs min/max over the LPs). |
| Refutability | emptiness possible / only jointly / never / Wald test | Decides which data tests apply (§4.3). |
| Dial | $\eta\ge0$ with $\eta=0$ the strict version, or none | Gives `sensitivity()` and `breakdown()` for free. |

---

## 2. Audit of the GPT conversation

### 2.1 Claims and verdicts
| # | Claim | Verdict | Basis |
|---|---|---|---|
| G1 | JKST bounds on $w_1$ map to bounds on $(\delta_W,\delta_B)$, and re-imposing those adds nothing about $D$ | **Correct.** In 2×2 the pair is a line segment, $\delta_W/\operatorname{Var}(X_N)+\delta_B/\operatorname{Var}(X)=D_{ER}-D_{NM}$. Re-imposing the marginals changes the $D$ interval by $<10^{-15}$. | M5(a); [V1] |
| G2 | Sign restrictions on $\delta_W,\delta_B$ and CR give $[L_J,\min(U_J,D_{ER})]$, $[\max(L_J,D_{NM}),U_J]$, and $[L_J,U_J]\cap[\min,\max](D_{ER},D_{NM})$ | **Correct** (proved below) | EGGHM Props 4 and 6 hold for mean-level quantities [V1] |
| G3 | "Incompatible ⇒ evidence against the assumptions" | **Correct only jointly with LCE.** Model-free, the $\delta_B$ sign and CR are never refutable (M3). Under LCE, a wrong $\delta_B$ sign gave an empty set in 223 of 2000 designs; but in 1777 it gave a non-empty interval that excludes the truth [V2]. Emptiness is a one-sided test with low power. | M3; [V2] |
| G4 | Pointwise ordering $c_1-w_1+d_1x\ge0$ on the hull | **Already implemented** as `order` (mean level, hull) | PHASE4; M5(b) |
| G5 | Pointwise contextual monotonicity $w_1+d_1x\ge0$ | **Correct algebra; new as a module.** It is EGGHM's local within term (N2.3), $x m_1'+(1-x)m_0'$, written under LCE. It is weaker than the existing `sign` rows and (given C1) stronger than global $\delta_W\ge0$. | P-N; [V3] |
| G6 | Unidentified dimension $(C-1)\binom R2$ | **Correct** | 2a Lemma 1 |
| G7 | R×C gives "more opportunity" | **Correct in a precise sense:** P-R and P-S. The caveats: W1 vs W2 for $R\ge3$; N2 no longer implies $\delta_W$ for $R\ge3$; a known table pins only $(R-1)(C-1)$ directions. | §5; [V5], [V7] |
| G8 | 2×3: adding-up gives $\sum_ct_c=0$ | **Correct but not new:** it is built into $K$. The gain comes from joint restrictions across columns, which are strictly tighter than the best single-column restriction in 18% of designs [V6]. M4 reported about half, but under a different, model-free set-up. | M4; [V6] |
| G9 | Joint sampling uncertainty | **Correct concern; already handled** if all rows are written in $k$ and enter one CLR stack (P3) | P3 Thm 5 |

### 2.2 Proof of G2 (2×2, any mean-level model)
1. EGGHM Prop. 4: $D_{ER}-D=\delta_W/\operatorname{Var}_p(X)$, and EGGHM Prop. 6: $D_{NM}-D=-\delta_B/(\mu(1-\mu))$. Both are accounting identities for any $m$ (PHASE2c Lemma M1), so they hold for the model-implied quantities. [M1; checked to $3\times10^{-16}$ in V1]
2. With $\operatorname{Var}_p(X)>0$ and $0<\mu<1$: $\delta_W\ge0\iff D\le D_{ER}$, and $\delta_B\ge0\iff D\ge D_{NM}$. [1]
3. Under LCE, $D$ is affine and non-constant in $w_1$ when $\gamma<1$, so the JKST set maps to $[L_J,U_J]$. Intersecting with the half-lines from step 2 gives G2. CR is the union of the two same-sign cases, i.e. $D$ lies between $D_{NM}$ and $D_{ER}$. ∎

**How often each adds information** [V2], over 2000 random LCE designs with true signs imposed:

| Restriction | Strictly tightens JKST | Redundant | Median width ratio |
|---|---|---|---|
| $\delta_W$ sign | 83% | 17% | 0.70 |
| $\delta_B$ sign | 89% | 11% | 0.67 |
| CR, sign unknown (union) | 80% of the 1060 designs where CR holds | — | 0.69 |

**Redundancy diagnostic.** "Redundant" happens exactly when the JKST-implied interval for that $\delta$ excludes 0. The package should print this table ("signs already implied by the data") before any sign is imposed.

---

## 3. Ladder (a): assumptions, from none to strong

Assumptions on different tracks are **non-nested**: report them side by side. Within a track, each line implies the one above it, so bounds shrink monotonically. Ticks mark refutability by aggregates.

### Base (no behavioural assumption)
| Rung | Content | Refutable | Status |
|---|---|---|---|
| L0 | Duncan–Davis / MOB, accounting identity only | no | done |
| L0+F | the same, computed sharply from $P_{T\mid X}$ (M-SHARP, avenue F) | no | pinned (2nd release) |
| L0+D | structural zeros; known tables at other levels (E4/E5); validation microdata (E6) | ✓ (inconsistent tables) | E4/E5 high priority |

### Track 1: direction (no functional form)
| Rung | Content | Refutable (2 groups) | Status |
|---|---|---|---|
| T1.1 | $\delta_B$ sign, **or** CR with unknown sign | never (M3) | proved (2c) |
| T1.2 | $\delta_W$ sign | ✓ | proved (2c) |
| T1.3 | joint signs over several columns (M4) | ✓ | proved (2c); [V6] |
| T1.4 | N2: each group's mean monotone in composition, any shape | ✓ | outer LP (2c) |
| T1.5 | local-everywhere signs (EGGHM Props 12–13, N2.3) | ✓ | non-nested with T1.2 (V6 of 1b) |
| T1.6 | unit-level ordering, $\beta_{ir}\ge\beta_{is}$ for each unit (N1/B5) | ✓ | exact (M2) |

### Track 2: smoothness (functional form)
| Rung | Content | Refutable | Status |
|---|---|---|---|
| T2.1 | Bernstein/polynomial outer–inner ladder (avenue D); weaker than LCE | ✓ | LP pair designed |
| T2.2 | L1′ quadratic contextual effects | ✓ | second release |
| T2.3 | LCE: group-specific (M2), shared with deviations ≤ δ (M1δ), shared (M1), pooled (M0) | ✓, plus a Wald test of M1 vs M2 | done |
| T2.4 | single-index contextual effects (§12.3(b)) | ✓ when $R\ge3,C\ge3$ or $R\ge4$ | conjecture |
| T2.5 | L4 row-constant (one group has no contextual effect) | ✓ | done |
| T2.6 | L5/L6 Goodman (within groups / overall) | ✓ | done |

### Track 3: scale and dispersion (how big effects or residuals may be)
| Rung | Content | Dial | Status |
|---|---|---|---|
| T3.1 | κ / slope ceilings (M-CEIL, §12.3(a)) | κ, $s$ | partly checked |
| T3.2 | second-scale context B1(M) (M-CTX2) | $M$ | simulated (NEW_AVENUES §2.4) |
| T3.3 | CAR(Z)-δ (§14) | δ | proposal |
| T3.4 | bounded residual dispersion ε; binomial second moments A1/A2 (M-DISP) | ε | simulated |

### Cross-track combinations worth building
In decreasing order of interest:
1. **LCE ∩ T1** (signs, CR, weighted-local, ordering): M5 and §5. These are cheap LP rows, and V2 shows large gains.
2. **T1.6 unit ordering ∩ LCE**: tightens the sharp set *and* the model rung. It needs the actual-vs-model allowance (⚠).
3. **LCE ∩ M-CTX2**: B1(M) is defined on top of M1 dummies (NEW_AVENUES §12.4).
4. **LCE ∩ M-CEIL**: a slope cap on M0 is an LP (2×2) or an ellipsoid (R×C).
5. **T1 ∩ M-DATA**: a known table fixes $(R-1)(C-1)$ directions; signs can cut the rest.
6. **T1.4 N2 ∩ T3.1 κ**: a variance-capped monotone shape; the outer program is convex.

---

## 4. Module catalogue

### 4.1 Comparison table

**Part A: where the bound comes from and what it assumes.**

Information source:
* *ID*: accounting identity only.
* *Q*: $E[T\mid X]$, which identifies $Q$.
* *P*: $P_{T\mid X}$, the DD kinks.
* *V*: $\operatorname{Var}(T\mid X,N)$ and its $1/N$ slope.
* *Z*: composition at a second scale.
* *Ext*: an external table or sample.
* *Con*: other contests.
* *Per*: another period.

Plausibility tier: 0 data only · 1 direction · 2 shape · 3 parametric · 4 constancy.

| Module | Information source | Assumption (plain / formal) | Restricts | Target | Tier | Refutable by aggregates | Form | Sharp |
|---|---|---|---|---|---|---|---|---|
| **M-DATA** zeros, E4/E5/E6 | ID + Ext | Known cells or tables are correct / equality rows on cells or on $\varphi(K)$ | unit / mean | both | 0 | ✓ (inconsistent tables) | rows | sharp given the table |
| **M-SHARP** (F) | P | none (a computational improvement) | unit | actual | 0 | — | series LP | sharp (goal) |
| **M-SIGN** | ID (model-free) or Q (on LCE) | One group is more supportive on average within units ($\delta_B$); rates rise with own-group share on average ($\delta_W$); CR = same sign / $\delta\gtrless0$ | mean (K-rows) or realised (N1) | model / actual | 1 | $\delta_W$ ✓; $\delta_B$, CR never (model-free), jointly with LCE ✓ | rows; CR = union | sharp (M2) |
| **M-ORD** | Q (mean) or ID (unit) | Group $r$ at least as supportive as $s$ at every composition / every unit / $m_r\ge m_s$ on the hull, $\beta_{ir}\ge\beta_{is}$ | mean or unit | model / actual | 1 | ✓ | rows | sharp |
| **M-LOC** (new) | Q | At every composition, moving residents from $s$ to $r$ does not lower the composition-weighted group rates / $(x^\top H)_r-(x^\top H)_s\ge0$ on the hull | mean | model | 1–2 | ✓, and jointly with M-ORD it is a test of the sign of the slope (§5.1) | rows | sharp under LCE |
| **M-N2** | Q + P | Each group's mean is monotone in chosen composition directions, any shape | mean on a grid | model | 2 | ✓ | grid LP | outer |
| **M-LCE** family (M0/M1/M1δ/M2, L1′, L4–L6) | Q (+ P for the sharp version) | Group means are linear (quadratic for L1′) in composition | mean | model (+ allowance) | 3–4 | ✓; Wald test M1 vs M2 | rows | sharp for the model-implied rate |
| **M-BERN** (D) | Q | Group means are polynomials of degree ≤ $J$ | mean | model | 2–3 | ✓ | outer/inner LP pair | outer/inner |
| **M-CEIL** (§12.3 a, b) | Q | Contextual variation is small: $\operatorname{Var}(E[\beta_r\mid X])\le\kappa\operatorname{Var}(Y)$, or a slope cap $s$; single index | mean | model | 2–3 | ✓ at small κ | rows (2×2) / SOCP | outer (binned) |
| **M-CTX2** (B) | Q + Z | Own-composition effect ≤ $M$ × the identified neighbourhood effect | mean | model | 3 | B2/B3 ✓; B1(M) only via $M$ | rows | sharp in the LP |
| **M-DISP** (A, §13.5) | V | Binomial voting within cells; residual dispersion ≤ ε; covariance ceilings | unit / second moment | actual | 3 | ✓ (overdispersion slope) | SOCP/LMI or cuts | outer |
| **M-CAR** (§14) | Q + covariates | Conditional on Z, composition is as-if random up to δ | mean, conditional on Z | model | 3 | partly (§14.3) | rows per Z cell | outer |
| **M-MULTI** (C) | Con | Contests share contextual structure (contests as districts; single index) | link | model | 3–4 | ✓ | rows | sharp in the LP |
| **M-PER** (DiD: FD, R1/R2/R4, D4) | Per | Changes are linear in composition; one restriction on $\Delta K$; or signs on changes | link | DiD | 1–4 | some (equal quadratics) | rows / union | D0, D4 exact; FD sharp in the LP |
| shifters (§13.6) | excluded variable | A composition shifter that does not affect rates directly | mean | model | 3 | IV-type test | rows | open |

**Part B: dials, power, practicality, extension and status.**

| Module | Dial / breakdown | Identifying power (evidence) | Estimability | R×C behaviour | DiD analogue | Data needs | Proof status | Reference |
|---|---|---|---|---|---|---|---|---|
| M-DATA | δ (table error) / yes | E4 can point-identify $K$; one table pins $(R-1)(C-1)$ directions (proved) | good | E4 with ≥2 tables of differing composition | per-period tables | published cross-tabs | E4/E5 proved in part (NA §13) | NA §5, §13 |
| M-SHARP | — | small–medium (NA §9) | series tuning | general | per period | none | consistency open (R9) | NA §6; PHASE4 |
| M-SIGN | slack $\eta$ on $\delta$ / yes | width ratio 0.67–0.70 on JKST [V2] | good | W1 vs W2 ambiguity for $R\ge3$ | D4 (exact) | none | proved (M2–M5) | 2c |
| M-ORD | slack / yes | large with M-LOC (sandwich) [V5] | good | each pair cuts $R-1$ directions (P-R) | MTR on changes | none | proved (M5(b), M2) | 2c, PHASE4 |
| M-LOC | slack / yes | median width ratio 0.75–0.90 alone, 0.47–0.77 with M-ORD [V5] | needs hull support in the chosen direction | same functional as M-ORD (P-S) | spillover-sign sandwich (DiD-S) | none | P-S proved; P-N proved apart from C1 | §5 |
| M-N2 | Lipschitz slack / yes | medium | τ estimation; grid choice | W1/W2; N2 no longer implies $\delta_W$ for $R\ge3$ | N2 on changes | none | outer validity open (R7) | 2c §4 |
| M-LCE | δ (M1δ) / yes | medium–large | good (HC3, CLR) | $\dim\mathcal K=(C-1)\binom R2$ | FD, M1-periods | none | proved (2a, 2b, P3) | PHASE1–4 |
| M-BERN | degree $J$ / — | negative by design (robustness) | good | general | yes | none | designed (NA §4) | NA §4 |
| M-CEIL | κ, $s$ / endpoint annotation | binds only at small κ in the suffrage data | good | ellipsoid in $K$ | calibrates $\delta$ in FD3 | none (κ benchmarked on covariates) | validity open (R8) | NA §12.3 |
| M-CTX2 | $M$ / yes | large, smooth in $M$ (NA §2.4) | good | needs a rank condition (NA Prop. B1) | second-scale spillovers bound $\Delta w_1$ | a second composition $z$ | B1 in progress (R8) | NA §2 |
| M-DISP | ε / yes | large in theory (A1), weak in practice | poor: needs thousands of small units; plug-in biased | LMI | per period | unit sizes | simulated | NA §1, §13.5 |
| M-CAR | δ / yes | medium | needs covariates with overlap | general | conditional parallel trends | covariates | proposal | NA §14 |
| M-MULTI | cross-contest δ / yes | medium–large (data-dependent) | data-dependent | general | — | ballots with several contests | partly | NA §3 |
| M-PER | δ (FD3) / yes | point (R1) to wide (FD4) | good | one restriction per column for a single-column DiD | — | ≥2 periods | D0, FD2–FD4, D4 proved | 2d |
| shifters | — | exploratory | needs a credible shifter | — | — | shifter | open | NA §13.6 |

**Set aside** (NEW_AVENUES §7, §9): G symmetry, H deconvolution, I entropy/OT, J Choquet.
* They add little or no power, or they assume the answer: maximum entropy *is* "no within-unit association".
* The honest form of I is the KL-ball to an explicit reference (§12.3(e)), which belongs under M-DATA with a dial.

**Structural limits** (proved in NEW_AVENUES §0). These cap what any combination can deliver:
1. Ballot margins shared across contests add nothing without a cross-contest assumption.
2. One known aggregate table pins only $(R-1)(C-1)$ of the $(C-1)\binom R2$ directions.
3. Poll-level shocks never reveal $K$, however well they are deconvolved.

### 4.2 Module × module layering matrix
Codes:
* **⇒**: the row module implies the column module, so adding the column is redundant.
* **⇐**: the column module implies the row module.
* **∩**: can be layered as one LP, and the combination can strictly tighten (evidence noted).
* **∪**: needs a union of LPs.
* **≠**: non-nested; report side by side, or intersect only if both are believed.
* **⚠**: incoherent as stated (footnote says how to fix).
* **?**: unproved.

Upper triangle only.

| | DATA | SHARP | SIGN-B | SIGN-W | CR | ORD | LOC | N2 | LCE | BERN | CEIL | CTX2 | DISP | CAR |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **DATA** | · | ∩ | ∩ | ∩ | ∪ | ∩ | ∩ | ∩ | ∩ᵃ | ∩ | ∩ | ∩ | ∩ | ∩ |
| **SHARP** | | · | ∩ᵇ | ∩ᵇ | ∪ᵇ | ∩ | ∩ | ∩ | ∩ᵇ | ∩ | ∩ | ∩ | ∩ | ∩ |
| **SIGN-B** | | | · | ∩ [V2] | ⇐ᶜ | ⇐ (M5b) | ∩ | ∩ | ∩ [V2] | ∩ | ∩ | ∩ | ∩ | ∩? |
| **SIGN-W** | | | | · | ⇐ᶜ | ∩ | ⇐ (2×2 LCE, C1)ᵈ | ⇐ (R=2), ≠ (R≥3) | ∩ [V2] | ∩ | ∩ | ∩ | ∩ | ∩? |
| **CR** | | | | | · | ∩ | ∪ | ∪ | ∪ [V2] | ∪ | ∪ | ∪ | ∪ | ∪? |
| **ORD** | | | | | | · | ∩ sandwichᵉ | ∩ | ∩ | ∩ | ∩ | ∩ | ⚠ᶠ | ∩? |
| **LOC** | | | | | | | · | ⇐ (R=2) | ∩ [V3, V5] | ∩ | ∩ | ∩? | ⚠ᶠ | ? |
| **N2** | | | | | | | | · | ≠ (LCE ∩ signed slopes ⊂ N2) | ≠? | ∩ | ? | ⚠ᶠ | ? |
| **LCE** | | | | | | | | | · | ⊂ BERN (J ≥ 1) | ∩ | ∩ (built on M1) | ∩ᶠ | ∩ |
| **BERN** | | | | | | | | | | · | ∩ | ? | ? | ? |
| **CEIL** | | | | | | | | | | | · | ∩ (NA §14.5) | ∩? | ∩ |
| **CTX2** | | | | | | | | | | | | · | ? | ∩ (§14.5) |
| **DISP** | | | | | | | | | | | | | · | ? |

Footnotes:
* a. One table fixes at most $(R-1)(C-1)$ directions, so in R×C it cannot make LCE point-identified on its own. With $R=2$ it can (NA §0).
* b. SHARP changes how a rung is computed, not which assumption it makes. The sharp version of a union is the union of the sharp versions.
* c. CR is the union of (B+, W+) and (B−, W−). Each same-sign pair implies CR.
* d. Under 2×2 LCE: weighted-local ⇒ $\delta_W\ge0$ if $x^*\in$ hull (P-N), which conjecture C1 says always holds. Model-free, local-everywhere does not imply global (PHASE1b V6). For $R\ge3$ it is open.
* e. P-S: the same functional is bounded from both sides. With a spanning tree of pairs, the pair bounds $K$ without DD (P-R).
* f. Unit-level or second-moment restrictions combined with mean-level rows: the target must be the actual rate with P3's allowance, or the unit restriction must be stated at the mean level. Otherwise the two sets refer to different objects.

**Rules used to fill the matrix** (each needs its proof reference in §8):
1. Two modules restricting the same object for the same target intersect, and the result is a subset of each.
2. ⇒/⇐ entries are theorems; an entry without a proof is shown as ?.
3. A cell where a structural limit binds gets a footnote instead of ∩.
4. Every ∩ marked with [V·] has at least one simulated design in which the combination is strictly tighter.
5. Unions are reported as min/max over their component LPs, and inference uses the union CLR (R10).

### 4.3 Data tests that narrow the assumption space
| Test | Refutes | Form | Status |
|---|---|---|---|
| Emptiness of a rung (with noise correction) | that rung, or the combination | P3 test | done |
| "Already implied" table: JKST-implied intervals for each $\delta$ and each pairwise functional $g_{rs}$ | nothing; it shows which assumptions are redundant | LP | new (G2) |
| **Sign of the identified directional slope** $q_{rs}(x)=(Q_cx)_r-(Q_cx)_s$ on the hull | the combination M-ORD ∩ M-LOC in that direction (P-S) | pointwise inequality on $\hat Q$ | new |
| Wald test M1 vs M2 | shared contextual effects | Wald | done |
| Equal quadratic coefficients across periods ($\Delta d_1=0$) | R1 + R2 | Wald | 2d |
| Placebo DiD on pre-periods | parallel trends (§8 of 2d) | interval excludes 0 | 2d |
| Overdispersion slope | binomial cells (A1) | regression on $1/N$ | NA §1.6 |
| B2/B3 overidentification | second-scale exclusion | Wald / emptiness | NA §2 |
| Hull and support diagnostics (dimension, vertices) | how much a local or pairwise module can bind (P-R needs full-dimensional support) | geometry | new |

---

## 5. New results for the extension modules

### 5.1 Sandwich lemma (P-S), any R, one column
Write $q_{rs}(x):=(Qx)_r-(Qx)_s$, which is identified, and $g_{rs}(x;k):=(Kx)_r-(Kx)_s$, which is linear in $k$. Define
$$O_{rs}(x):=m_r(x)-m_s(x)=(Hx)_r-(Hx)_s,\qquad L_{rs}(x):=\sum_tx_t\,\partial_{e_r-e_s}m_t(x)=(x^\top H)_r-(x^\top H)_s.$$
$O$ is the between-group gap. $L$ is the composition-weighted contextual derivative when a resident of group $s$ is replaced by one of group $r$.

**Proposition P-S (proved).**
* (a) $O_{rs}(x)+L_{rs}(x)=2q_{rs}(x)=\partial_{e_r-e_s}\tau(x)$. This is the local Monotone-EI decomposition.
* (b) $O_{rs}(x)=q_{rs}(x)+g_{rs}(x;k)$ and $L_{rs}(x)=q_{rs}(x)-g_{rs}(x;k)$.
* (c) Hence $O_{rs}\ge0$ and $L_{rs}\ge0$ at $x$ $\iff$ $\lvert g_{rs}(x;k)\rvert\le q_{rs}(x)$. This is non-empty only if $q_{rs}(x)\ge0$, which is a data test.

*Proof.*
1. $\tau(x)=x^\top Hx$, so $\nabla\tau=(H+H^\top)x=2Qx$. [def]
2. $x^\top H=(H^\top x)^\top$ and $H^\top=Q-K$, since $K$ is antisymmetric. [def]
3. Substitute $H=Q+K$ into $O$ and $H^\top=Q-K$ into $L$. This gives (b); adding them gives (a). [alg]
4. (c) follows from (b). ∎

Checked to $4\times10^{-16}$ [V4]. With $R=2$, $x=$ the group-1 share and $k=w_1$, this is exactly GPT's pair $c_1-w_1+d_1x\ge0$ and $w_1+d_1x\ge0$.

### 5.2 Rank lemma (P-R), binary outcome or per column
For a pair $(r,s)$, $g_{rs}(x;k)=\sum_t(K_{rt}-K_{st})x_t$. As $x$ ranges over a full-dimensional subset of the simplex, the functionals reached are $K_{rs}$ and $K_{rt}-K_{st}$ for $t\ne r,s$: that is, $R-1$ independent directions.

**Proposition P-R.**
* (a) One pair restricts an $(R-1)$-dimensional subspace of $\mathcal K_c$. **Proved** (above).
* (b) If the pairs form a star tree $\{(0,s)\}$, all $\binom R2$ directions are reached. **Proved**: $K_{0s}$ comes directly, and then $K_{st}=K_{0t}-(K_{0t}-K_{st})$.
* (c) The same holds for any spanning tree. **Checked** for chain trees with $R=3,4,5$ [V7]; general proof by induction on the tree (R5).
* (d) Combined with P-S(c), M-ORD ∩ M-LOC on a spanning tree bounds $K_c$ with no Duncan–Davis rows. **Follows from (b) and P-S** for star trees.

Rank is necessary, not sufficient, for tightening: it says which directions can be cut, and the signs and the data decide whether they are. [V5] shows the practical effect in 3×2. With true restrictions only, ordering and local rows tighten about 2/3 of group-rate targets (median width 0.77 of DD under both). With near-1-D support they tighten more often (≈94%, median 0.47), because DD alone is weaker there.

### 5.3 Nesting in 2×2 under LCE (P-N)
Let $a=\mu\operatorname{Var}_{w^1}(X)$, $b=(1-\mu)\operatorname{Var}_{w^0}(X)$, and $x^*=a/(a+b)$.
* (a) Signed slopes ($w_1\ge0$, $w_1+d_1\ge0$) ⇒ weighted-local on the hull. **Proved**: $w_1+d_1x=(1-x)w_1+x(w_1+d_1)$.
* (b) $\delta_W=(a+b)(w_1+d_1x^*)$. **Proved** (M5(c); asserted in V3). Hence weighted-local on the hull ⇒ $\delta_W\ge0$ whenever $x^*\in[x_L,x_U]$, and the implication is strict when $d_1\ne0$ and $x^*$ is interior.
* (c) **Conjecture C1:** $x^*$ always lies in the composition hull. Evidence: 2000 random designs plus an adversarial search over 50,000 discrete designs with 2–4 support points and extreme weights found no violation beyond $10^{-6}$ [V3]. A pairwise representation $\operatorname{Var}_w(X)=\tfrac12\sum w_iw_j(X_i-X_j)^2$ writes $x^*$ as a weighted mean of pair ratios, but individual ratios can leave the hull, so a proof needs a cancellation argument (R1).

---

## 6. Ladder (b): build order with rigor checklists

**Standard checklist** (every step must clear it before its module enters the API):
- [ ] **S1 Statement.** The formal assumption, its object (unit/mean/link), its target (actual/model), and a plain-language version.
- [ ] **S2 Structure.** Affinity in $k$ (or convexity), and the LP/union/grid form; for grids, outer validity.
- [ ] **S3 Sharpness.** Sharp, outer or inner, with a proof or a stated gap.
- [ ] **S4 Nesting.** Its relation to every existing rung, filled into the §4.2 matrix with proof references.
- [ ] **S5 Refutability.** Can the set be empty? Which test applies? Size of the test.
- [ ] **S6 2×2 closed form,** checked against brute-force enumeration (small integer populations) and the LP.
- [ ] **S7 R×C property tests.** A Python reference implementation agrees with R to $10^{-9}$; restrictions only narrow the bound; the truth is covered when the assumption holds.
- [ ] **S8 Inference.** Rows enter the CLR stack; coverage simulation at nominal levels, including near-empty and redundant cases.
- [ ] **S9 Dial and breakdown value,** with monotone widening in the dial proved.
- [ ] **S10 Documentation.** An entry in the plain-language guide, with `meta` filled.
- [ ] **S11 Real data.** NC `senc` (known truth), NZ, Wisconsin suffrage. Report width, emptiness and the redundancy diagnostic.
- [ ] **S12 Independent review.** A fresh session re-derives the proofs before the status changes to proved.

**Order** (the reasons are proof readiness, then power per unit of assumption, then data availability):

| Step | Work | Proof items to clear | Why at this point |
|---|---|---|---|
| 1 | Merge the branches; refactor `eiBounds` restrictions into the module registry with `meta`; add K-row $\delta$ signs, CR (union) and the redundancy table | R2, R3, R10 (union) | M5 is proved; immediate gains [V2] |
| 1b | M-SHARP series bound | R9 | tightens every rung with no new assumption |
| 2 | M-LOC in 2×2 and 2×C; data test for the sign of the directional slope | R1 (C1), R11a | P-S and P-N are mostly proved |
| 3 | 2×C joint column signs and ordinal cumulative columns | R6 | M4 proved; [V6] |
| 4 | R×2 pairwise ordering + local (sandwich) on spanning trees; hull diagnostics | R4, R5 | P-R(a, b) proved; (c) needs proof |
| 5 | M-DATA E4/E5 (shared $K$ across validation tables) | NA §13 items | high power, no behavioural assumption |
| 6 | M-CTX2 B1(M) | R8 (B1 rank condition) | the most practical new rung (NA) |
| 7 | M-CEIL κ / slope caps; annotate endpoints with their implied κ | R8 | reporting device; DiD dial calibration |
| 8 | M-N2 grid LP with estimated τ | R7 | needs an argument about estimated τ |
| 9 | M-MULTI (contests) | NA §3 | data-dependent |
| 10 | M-CAR(Z)-δ | R8 | needs covariates |
| 11 | M-DISP (A1/A2, covariance ceilings) | NA §1, §13.5 | strong in theory, weak estimability |

This differs from NEW_AVENUES §13.8, which put E4 first. Steps 1–4 come first here because their theory is already proved or nearly so, and they reuse the existing LP with no new data. E4 stays the top *data-side* priority.

---

## 7. Extension to difference-in-differences (same units, two elections)

Setting: PHASE2d. Fixed population; rows $v$ (veteran) and $n$; per-period $Q_t,K_t$; changes $\Delta H=\Delta Q+\Delta K$ under first-difference LCE (FD, Lemma FD1). In 2×2, $\text{DiD}=\text{(identified)}-(r+m_n)\Delta w_1$ (MENU item 18). Each module has a change-version obtained by applying it to $\Delta m$ in place of $m$:

| Module on changes | Restriction | Effect on DiD | Status |
|---|---|---|---|
| M-SIGN-Δ (= D4) | sign of $\delta^\Delta_W$ / $\delta^\Delta_B$ / CR on the realised changes | one-sided cuts at Goodman-on-changes and NM-on-changes | proved (D4), exact |
| **M-ORD-Δ** (MTR) | $\Delta m_v(x)\ge\Delta m_n(x)$ on the hull: $\Delta c_1-\Delta w_1+\Delta d_1x\ge0$ | upper bound on $\Delta w_1$ | rows; proved (affine) |
| **M-LOC-Δ** (weighted spillover sign) | $\Delta w_1+\Delta d_1x\ge0$ on the hull: replacing a non-veteran with a veteran does not lower the composition-weighted change | lower bound on $\Delta w_1$ | rows; proved (affine) |
| **DiD-S** = ORD-Δ ∩ LOC-Δ | P-S on changes: $\lvert\Delta g(x)\rvert\le\Delta q(x)$ | closed-form interval for $\Delta w_1$, so DiD $\in$ identified $-(r+m_n)\cdot[\max_x(-\Delta d_1x),\ \min_x(\Delta c_1+\Delta d_1x)]$ | **proved** (P-S with $\Delta$); empty iff $\Delta q<0$ somewhere (test) |
| M-CTX2-Δ | spillovers measured at a second scale bound the own-composition $\Delta w_1$ | interval around R1, of width $\propto M$ | proposal; replaces FD3's free δ with an identified yardstick |
| M-CEIL-Δ | $\lvert\Delta w_1\rvert\le s$, with $s$ benchmarked from pre-period variation (relative magnitudes, Rambachan–Roth) | FD3 with a calibrated δ | FD3 proved; calibration open |
| M-LCE levels ∩ M-SIGN-Δ | per-period LCE plus signs on changes | non-nested with FD (2d §10) | open (R11) |
| M-DATA-Δ | per-period known tables; the pre-period table pins $K_\text{pre}$ directions | tightens levels and DiD | proposal |
| 2×3 common base | DiD-S per column with $\sum_c\Delta t_c=0$ linking columns through the DD rows | joint columns may tighten (cf. V6) | R6 analogue |

**R×C DiD.** A single-column DiD needs only that column's $\Delta K_c$ (2a Lemma 2; PHASE1b B1 correction). P-R then says a spanning tree of pairwise change restrictions bounds every direction of $\Delta K_c$.

**Data tests for DiD:**
* the sign of the identified $\Delta q_{rs}(x)$, which is required for DiD-S;
* equal quadratic coefficients across periods (R1 + R2);
* a placebo DiD on pre-periods;
* emptiness at each rung.

**Composition change.** Every module above applies to the *stayer* population of C1. With interval-valued bases (2d §7, D.7), the rows are written in cell-count coordinates and remain an LP.

**Still to prove (R11):**
* sharpness of DiD-S within FD;
* FD4 ∩ D4 (the outer set intersected with the exact sign cut, which mixes mean-level and realised quantities, ⚠);
* LCE-in-levels ∩ sign-on-changes;
* inference for unions of sign patterns across periods.

---

## 8. Proof register

Status: **P** proved · **C** checked numerically · **J** conjecture · **O** open.

| ID | Claim | Status | Strategy / check | Blocks |
|---|---|---|---|---|
| R1 | C1: $x^*\in$ hull (2×2), hence weighted-local ⇒ $\delta_W\ge0$ under LCE | J (C [V3]) | pairwise variance representation plus a cancellation argument; try two-point designs first | the M-LOC ⇒ SIGN-W cell |
| R2 | An added restriction moves an endpoint iff it removes the whole optimal face (one feasibility LP); for one row, iff every optimal dual is positive (a solver's dual alone does not certify); implied ⇔ tightens no target; 2×2 closed form = G2 | **P** (5a §2; `verify_R2.py`; independently reviewed, fixes applied, 5a §4) | — | redundancy diagnostic |
| R3 | LCE ∩ any closed convex restriction on $k$ is sharp for the model-implied target; CR is a union, a single slab for W1/2×2, and may have a gap otherwise | **P** (5a §1; `verify_R3.py`; independently reviewed, fixes applied, 5a §4) | — | M-SIGN reporting |
| R4 | Refutability classification for $R\ge3$ (which signs, orderings and local rows can produce emptiness, model-free and under LCE) | O | construct feasible completions; counterexamples | §4.2 codes, tests |
| R5 | P-R(c): any spanning tree of pairs reaches all directions; behaviour with lower-dimensional support | P (star); C (chain) [V7] | induction on tree leaves; support rank | step 4 |
| R6 | When joint column restrictions in 2×C are strictly tighter than the best single one; ordinal cumulative version | C [V6] | LP geometry in the $(C-1)$-dimensional $t$-space | step 3 |
| R7 | N2 grid LP: outer validity with estimated τ; convergence as the grid is refined | O | uniform confidence band for τ, then a monotone envelope | step 8 |
| R8 | B1(M), κ, CAR-δ: each is a valid outer bound, monotone in its dial; B1 rank condition | partly (NA) | per module | steps 6, 7, 10 |
| R9 | Consistency of the M-SHARP series bound; never wider than the current $\Theta_X$ | O | sieve argument | step 1b |
| R10 | CLR validity with added rows, unions (min/max over LPs), near-emptiness, and the actual-rate allowance with new rows | **P** (5a §3: added rows; parallel disjunction = W1/2×2 CR; union of branch CIs at full level; emptiness level; plug-in bracket invalid; MC in `verify_R10.py`; independently reviewed, one error in (d) fixed, 5a §4). Still open: uniformity near point identification, and the CLR candidate collection (both inherited from P3) | — | all code |
| R11 | DiD: DiD-S sharp within FD; FD4 ∩ D4; LCE-levels ∩ sign-Δ; union inference | DiD-S validity P; rest O | as in 2d | §7 |
| R12 | P-S (sandwich) | P [V4] | — | — |
| R13 | G2 table (JKST + signs) | P (§2.2) [V1, V2] | — | — |
| R14 | Single-index dimension $C-1$ | J (NA check) | generic-rank argument | T2.4 |

**Standard for "proved":**
1. The proof is written in the numbered-step format used in PHASE2c, with a justification tag on each step.
2. It is checked against brute force or the LP, including adversarial designs.
3. An independent fresh-session review re-derives it (S12).
4. The §4 tables are then updated.

Until then, the package labels any module that depends on the claim as `experimental` in `meta$proof_status`, and `ladder()` prints that label.
