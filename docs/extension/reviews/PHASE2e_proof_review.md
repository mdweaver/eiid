# Proof review: PHASE2e (general DiD with ecological data; build step 2; register rows R11a–d, R11b′, R31, R32)

*Independent proof review of `docs/extension/PHASE2e_did_general.md` (cited as 2e) at branch `claude/did-theory`, head `7e3d8a6`. Proofs only; code is reviewed separately. The theory documents are not edited here. Checked against `PHASE6a_count_layer.md` (6a: Lemma 2.1, R19, Lemma 4.1, R27, Remark 4.2), `PHASE3_inference.md` (P3: §0, §0.5, Thm P3.5(a), (b2) as corrected, P3.8), `PHASE2c_monotone.md` (Def. 1, M1, M3, M5), `PHASE2d_did.md` (2d), `PHASE1b_monotone_did.md` (1b), `PACKAGE_DESIGN.md` §7–§8 (PD) and `next_steps/NEXT_STEPS_RECOMMENDATION.md` (NS) §3.2 and Part 3A. The standard and the severity scale are those of `reviews/PHASE6a_proof_review.md`.*

*Reviewer's checks: `phase2e_review_checks.py` (output `phase2e_review_checks_output.txt`, **15/15 PASS**, about 1 minute). Generators and LP code are the reviewer's own; nothing is imported from `verify_phase2e.py`. Tags: **[CX]** a counterexample to a statement of 2e as written (PASS = it behaves as this review says); **[V]** an independent verification; **[W]** a wording point shown numerically.*

*Reproduction of the author's script: `verify_phase2e.py` rerun with Python 3.13 and current scipy/numpy: output **identical** to `verify_phase2e_output.txt` (71/71 PASS, about 80 s).*

Severity: **error** = a statement that is false as written; **gap** = a missing hypothesis or step, where the claim is true once it is added; **wording** = imprecise or incomplete, with no false content.

---

## 1. Verdict in brief

2e is careful work and its core is sound. The stacking lemma, the margins-only rungs for T elections, the identities on changes, the closed forms for one-vs-rest rows, the sharp FD set via Lemma G′, the causal identity ATT = D − b + S, the joint causal LP, the union-of-LPs form of relative magnitudes, the degeneration proposition, the two-election allowance, the integer lattice and the TU theorem are correct in substance. I recomputed every number in the §9.2 worked example and in N21 with my own LP code and got all of them (check W, N21). Every claim of §11 about other documents is real.

Two statements are false as written, and both have one-line fixes:
* **Cor. E17(b)** (the condition under which the trend-dial set equals D ⊕ [−(M+Λ), M+Λ]) omits two of the four rate-limit conditions (issue 1; counterexample with counts, one coupling row). Under margins only the missing conditions are implied, so no number in 2e changes.
* **Prop. E6(b)**'s emptiness criterion omits the case $\bar D^+<0$ (issue 2; counterexample with four pure units).

The most substantive gap is in **Prop. E25** (one joint region for model-based pre-trends): the post-correction hypotheses of P3.5(b2)(ii), what is model-implied and what is actual inside the causal constraints, and how the allowances enter the constraints are not stated (issue 5). The labelling rule (Def. 5.1) omits the population assumption that E14 uses (issue 4), and the person-level MTR witness in E16(b) needs integer counts (issue 3).

| Claim | Verdict | Issues |
|---|---|---|
| Lemma E1 (stacking) | accept | — |
| Lemma E2 (panel tables add nothing; integer coupling) | accept | — |
| Prop. E3 (Δ0 for T elections; event-study set not a box; per-period LCE) | accept with fixes | 9 |
| Lemma E4 (rows on changes; one-number reduction; identities) | accept | — |
| Prop. E5 (Δ4 in R×C, T elections; CR; other rows by LP) | accept | 15 (labels) |
| Prop. E6 (count contrast on changes) | **accept with fixes** | **2** |
| Lemma E7 (FD in R×C) | accept ((d) inherits R5's status) | — |
| Lemma G′ | accept with fixes | 7 |
| Prop. E8 (sharp FD set; FD4 outer; 2×2 closed form) | accept with fixes | 7 |
| Prop. E9 (DiD-S closed form and sharpness) | accept with fixes | 8 |
| Prop. E10 (LCE contrasts on changes) | accept | — |
| Prop. E11 (per-period LCE ∩ model rows on changes) | accept | — |
| Lemma E13 (allowance ρ_i = 2(w_r + w_s)) | accept with fixes | 6 |
| Theorem E12 (R27 on the product) | accept (status inherited from the model interval) | 6 |
| Lemma E14 (ATT = D − b + S) | accept | 4 |
| Def. 5.1 (labelling rule) | **accept with fixes** | **4** |
| Prop. E15 (no trend assumption) | accept with fixes | 17 |
| Theorem E16 (joint causal LP: valid, sharp, refutable) | **accept with fixes** | **3**, 10 |
| Cor. E17 (trend dial, breakdown) | **accept with fixes** ((a), (c) accept; (b) false as written) | **1** |
| Def. 5.2, Prop. E18 (RR relative magnitudes, joint form) | accept with fixes | 10 |
| Prop. E19 (degeneration (a)–(d)) | accept with fixes ((d) through E25) | 5, 12, 13, 19 |
| Prop. E20 (smoothness) | accept; its "Reading" with fixes | 11 |
| Prop. E21 (MTR/ND/PS with spillover dial) | accept with fixes | 4 |
| Prop. E22 (spillovers; pre-periods uninformative) | accept | — |
| Prop. E23 (FD inference) | accept with fixes | 6 |
| Prop. E24 (unions, no multiplicity correction) | accept | — |
| Prop. E25 (one joint projection region) | **accept with fixes** ((c) as a computed counterexample) | **5** |
| Prop. E26 (integer lattice; signs; CR) | accept | — |
| Theorem E27 (TU for change rows, R = 2) | accept | 16 (remarks) |
| Remark 7.1 (causal integer ends O) | accept as stated | — |
| Gap register G1–G9 (§8.1) | accept as a register; the naive-fails checks are genuine | 14, 18 |
| Lemma G6 (coarsening) | accept | — |
| Prop. G7 (margin sets; R26 on changes) | accept with fixes | 18 |
| §9 plain language and worked example | accept with fixes (numbers all reproduced) | 11, 12 |
| §11 claimed errors in other documents | all real | — |

---

## 2. Numbered issues

### Issue 1 (error, Cor. E17(b)). The sufficient condition for D ⊕ [−(M+Λ), M+Λ] to be sharp omits two rate limits.

E17(b) says equality holds if, at a maximiser $n^+$ of $D_t$, $\theta^{PT}_t(n^+)\ge M+\Lambda$ and $B_{s^*ct}(n^+)\ge\Lambda$, and at a minimiser $n^-$, $\theta^{PT}_t(n^-)\le1-M-\Lambda$ and $B_{s^*ct}(n^-)\le1-\Lambda$. The proof sets $b=-M$, $S=\Lambda$ at $n^+$, so $\theta=\theta^{PT}(n^+)-M-\Lambda$, and asserts $\theta\in[0,1]$. That needs **also** $\theta^{PT}(n^+)\le1+M+\Lambda$; symmetrically, at $n^-$, $\theta^{PT}(n^-)\ge-(M+\Lambda)$. $\theta^{PT}=B_{r^*c0}+B_{s^*ct}-B_{s^*c0}$ ranges over $[-1,2]$, so neither is automatic.

*Counterexamples.* (i) Abstract (M = Λ = 0): $F$ the segment from $n^+=(B_{r0},B_{s0},B_{r1},B_{s1})=(.9,0,1,.2)$ to $n^-=(.5,.5,.1,.5)$. The stated conditions hold ($\theta^{PT}(n^+)=1.1\ge0$, $\theta^{PT}(n^-)=0.5\le1$), $\max D=-0.10$, but the sharp maximum of the ATT under parallel trends is $-0.15$. (ii) With counts (check "E17 [CX] … with counts"): one 2×2 unit, $a=(100,100)$, 100 yes votes at both elections, and one linear count row coupling the elections, $3y_0-2y_1=140$ ($y_t$ the treated-yes count). $\max D=-0.400$ at $\theta^{PT}(n^+)=1.200$; $\theta^{PT}(n^-)=0.933$; the sharp set under PT is $[-0.933,-0.800]$, not ending at $-0.400$.

*Margins only it is harmless* (check "E17 [V]", 200 designs with R = 2, a named pair and the rest in R = 3, Λ ∈ {0, .03}): under margins only every free count moves $D_t$ and $\theta^{PT}_t$ in opposite directions or leaves $\theta^{PT}$ fixed, so $n^+$ minimises and $n^-$ maximises $\theta^{PT}$, and the missing conditions follow from the stated ones. 2e's check E16 is margins-only, which is why it did not see this.

*Fix.* State (b) as an exact criterion: "the upper end equals $\max_F D_t+M+\Lambda$ **iff** some $n\in\arg\max_FD_t$ has $\theta^{PT}_t(n)\in[M+\Lambda,\,1+M+\Lambda]$ and $B_{s^*ct}(n)\ge\Lambda$; the lower end equals $\min_FD_t-M-\Lambda$ iff some $n\in\arg\min_FD_t$ has $\theta^{PT}_t(n)\in[-(M+\Lambda),\,1-M-\Lambda]$ and $B_{s^*ct}(n)\le1-\Lambda$. Under margins only the outer two conditions are implied." The "only if" is immediate: $\mathrm{ATT}\le D_t(n)+M+\Lambda$, with equality only for $b=-M$, $S=\Lambda$.

### Issue 2 (error, Prop. E6(b)). The emptiness criterion omits $\bar D^+<0$.

E6(b) gives the set $[\max\{0,\bar D^-\},\min\{\bar D^+,\mathrm{cap}\}]$ (correct) and says it is "empty iff $\Delta\tau_{\mathcal H}<\Delta\tau_{\mathcal L}$ … or the cap is below $\bar D^-$". It is also empty when $\bar D^+<0\le\mathrm{cap}$: the mix belief alone is refuted by the margins.

*Counterexample* (check "E6 [CX]"): bin $\mathcal L$ = a pure group-1 unit (20 people) and a pure group-2 unit (80); bin $\mathcal H$ = a pure group-1 unit (80) and a pure group-2 unit (20); yes votes $(10,40,40,4)$ then $(8,40,72,16)$. Then $\bar D=-0.15$ exactly, $\lambda=0.6$, $\Delta\tau_{\mathcal H}-\Delta\tau_{\mathcal L}=0.46>0$, cap $=0.767\ge\bar D^-$, and the LP with both beliefs is infeasible.

*Fix.* "It is empty iff $\max\{0,\bar D^-\}>\min\{\bar D^+,\mathrm{cap}\}$, i.e. iff $\Delta\tau_{\mathcal H}<\Delta\tau_{\mathcal L}$ (the beliefs are jointly refuted), or $\bar D^+<0$ (the mix belief is refuted by the margins), or the cap is below $\bar D^-$ (the rate belief is refuted)." The 2e check E6 compares LP with the interval formula, which is right, so only the sentence changes.

### Issue 3 (gap, Thm E16(b)). The person-level MTR witness needs integer counts; the vector of post-periods and unions in $F$ are not covered as written.

(i) Step 2 takes $z_i=n_{ir^*ct}\theta/B_{r^*ct}(n)$ and gives $Y_{jt}(0)=c$ to "$z_i$ of the $n_{ir^*ct}$ people". For a person-level statement, $z_i$ must be an integer. With an integer realised table and $\theta A_{r^*}\in\mathbb Z$ (Remark 7.1) integers $z_i\in\{0,\dots,n_{ir^*ct}\}$ with $\sum_iz_i=\theta A_{r^*}$ exist (greedy, since $\theta A_{r^*}\le\sum_in_{ir^*ct}$). For continuous tables "person by person" can only mean a statement about masses. The formula also divides by $B_{r^*ct}(n)$, which may be 0 (then $\theta=0$ under MTR and $z\equiv0$).
(ii) The comparison group's counterfactual must also be realisable unit by unit (counts in $[0,a_{is^*}]$ summing to $(B_{s^*ct}-S)A_{s^*}$); true, but say it.
(iii) E16 is stated for one $(\theta_t,S_t)$. E18(a) applies it to "the vector of post-periods". The proof extends verbatim (the constraints for different $t$ involve different $(\theta_t,S_t)$, and the witness is built per $t$), but the statement should say so.
(iv) (c) "a union of at most $K$ intervals": if $F$ is itself a union (CR branches), the bound is $K$ times the number of branches of $F$.

*Fix.* In (b): "With MTR in Γ and an integer realised table, the configuration satisfies MTR person by person; for continuous tables, as an inequality between masses." Replace the formula by "integers $z_i\le n_{ir^*ct}$ with $\sum z_i=\theta A_{r^*}$". State E16 for $(\theta_t,S_t)_{t\ge1}$ and amend (c).

### Issue 4 (gap, Def. 5.1; E14, E20, E21). The labelling rule omits the population; NA is not stated everywhere it is used.

E14 uses (FP): the same people at base and post, with fixed group membership. Without it the ATT of which population is undefined (G1–G4, G9 say so). Def. 5.1 lists groups, trend, spillover, NA and the descriptive base, but not the population.

*Where NA is used.* E14 states it; E15, E16, E19 and E22 state or inherit it. E20 uses $\delta_{-1}=D_{-1}$ (E14(b), hence NA) without saying so. E21's ND is written with the observed base rate $B_{r^*c0}$; it reads as "no decline of the untreated outcome" only under NA ($B_{r^*c0}=B^{(0)}_{r^*c0}$). E21's algebra does not need NA (ATT $=B_{r^*ct}-\theta_t$ is definitional), but its interpretation does.

*Fix.* Add to Def. 5.1: "0. the population for which the ATT is defined: the fixed population (FP), or the stayers with fixed group membership when (FP) fails (§8)". Add "(FP), NA" to the hypotheses of E20 and E21 and one sentence in E21 that ND is ND for the untreated outcome under NA.

### Issue 5 (gap, Prop. E25; E19(d)). The joint region needs the corrected P3.5(b2)(ii) hypotheses and a definition of what is model-implied inside the causal constraints.

E25(a) is the right idea and its union argument is correct. As written it is not yet a theorem:
1. **Positive definiteness.** Since the Batch B correction, P3.5(b2)(ii) needs $\Sigma$ positive definite with $q=\dim\beta$ (ii-a), or an exact reduction (ii-b). E25 does not state either. Stacking "the regressions of all elections, or of all period-pair changes" is fine; stacking levels together with changes, or overlapping changes (e.g. pairs (−1,0), (0,1) and (−1,1)), makes $\Sigma$ singular by an exact identity, which must then be removed as in (ii-b).
2. **Levels under FD.** E25 lists $B_{r^*ct}$ and $\theta^{PT}_t$ among the "model-implied linear functionals of β". Under FD levels are unrestricted, so these are not model-implied functionals of an FD rung; the rate limits $\theta\in[0,1]$ cannot be imposed through FD. E25 must either restrict to rungs that model levels (per-period LCE, E3(d)/E11) or drop the rate limits (an outer set).
3. **Model-implied vs actual inside the constraints.** RM, SD and the dial are assumptions on the **actual** counterfactual gap $b_t$, and $b_s=D_s$ (actual) for $s<0$ (E14(b)). If the branch LPs use model-implied pre-period DiDs, the assumption has silently become one about model-implied trends, and "the model-implied ATT" (which E25(a) says is covered) is not defined in §5. For the actual ATT every actual quantity that enters a **constraint** (each pre-period $D_s$ in $A(n)$, $D_t$, $B_{r^*ct}$), not only the target, needs its allowance; the natural construction adds box variables $E_q\in[-t_q,t_q]$ for each such $q$ (Lemma E13's ρ for each) and takes a union bound over those events and the projection event (6a Remark 4.2). "Add the allowance events of each actual target used" points to this but does not give the construction or its proof.
4. **E25(b)** ("can only be combined by the union bound … they cannot deliver the joint set") is a claim about all possible procedures and is not proved. It is a remark.

*Fix.* Restate E25 with: (ii-a)/(ii-b) of P3.5(b2) for the stacked $\beta$; a rung that models the levels whenever the rate limits are imposed (else drop them and say "outer"); a definition of the causal assumption in the model layer (on actual $b_t$, with actual pre-DiDs entering as model value + $E_q$); the box-variable construction and its union bound; and (b) as a remark. Until then R11b′(d) and R11c's "one joint region" item stay pending. My check "E25 [V]" confirms the projection coverage with correlated estimates (0.998/0.997/1.000 at the lower end, the upper end and the kink $D_{pre}=0$; plug-in 0.917/0.913), but that is the Gaussian skeleton only.

### Issue 6 (gap, Lemma E13, Thm E12, Prop. E23). Under FD the mean condition of P3-S1a is about changes only.

Lemma E13 assumes "P3-S1a for the stacked unit vectors" and applies 6a Lemma 4.1. 6a Lemma 4.1 step 1 needs $E[Z_i\mid V_i]$ to equal the model-implied value. P3-S1a states $E[\beta_i\mid\text{design}]=m(V_i)$; under FD only $E[\Delta\beta_i\mid V_i]=M_{\Delta H}(X_i)$ is assumed and the levels are free. So for an FD rung: (i) the hypothesis should read "$E[\Delta\beta_i\mid\text{design}]=M_{\Delta H}(X_i)$, independence across units, any dependence within a unit"; (ii) $\tilde\varphi$ is defined only for targets whose coefficients at $u$ and $t$ are opposite within each $(r,c)$ (functions of changes), which includes the DiD. The ρ formula, the attained bound, $S_D=(\sum_i(w_{ir}+w_{is})^2)^{1/2}$ and twice the gap allowance are all correct (my check "E13 [V]": ratio 0.342 under a dependent law, extreme law 0.967 ≈ 1 in 4000 replications, tail 0.0005). The Lindeberg condition of 6a Lemma 4.1(a) ($\max\rho_i^2/\sum\rho_i^2\to0$) is inherited and should be listed.

*Fix.* Add the FD form of the mean condition and the restriction on $\varphi$ to E13 and E23(b); list the Lindeberg condition.

### Issue 7 (gap, Lemma G′; Prop. E8). The selector in step 2 is a function of $(w,z)$, not of $w$; the conditioning variable.

Step 2 (⇐) uses "the correspondence $w\mapsto\{\xi\in\mathcal S(w):A\xi=\zeta(w)\}$". In 2d's Lemma G the latent vector is a random vector jointly distributed with $W$, not a function of $W$. The right object is $(w,z)\mapsto\{\xi\in\mathcal S(w):A\xi=z\}$, which has nonempty compact values on $\{z\in A\mathcal S(w)\}$ and a closed graph, hence a Borel selector $s$; put $\xi=s(W,\zeta)$. Then $E[A\xi\mid V]=E[\zeta\mid V]$. Also, $V$ must be a function of $W$ (as in Lemma G); E8(a) conditions on $V$ and E8(d) on $x$ alone: say $V=X$ (or the 2d convention) once. The conclusion is unchanged.

*Fix.* Replace the correspondence; add "$V$ a function of $W$; below $V=X$".

### Issue 8 (gap, Prop. E9(c)). "iff" holds for the parameter set; the LP claim needs finite support.

(i) For $R=2$ the closed form is the set of $\tilde D=I+s\Delta\kappa$, monotone in $\Delta\kappa$, so "sharp for the full law iff $[-q_{\min},q_{\min}]\subseteq\Theta^{FD}_I$" is right. For $R\ge3$ there is no closed form, and "the row set lies in $\Theta^{FD}_I$" is sufficient, not necessary, for the **target's** set to coincide (a projection can coincide while the sets of $\Delta k$ differ). (ii) "Otherwise (b) is strictly narrower and is computed as one LP": $\Theta^{FD}_I$ is an intersection over all $u$ and a.e. $v$; it is an LP (with selection variables per support point) when $(X,T_u,T_t)$ has finite support, e.g. the empirical law. (iii) For $A$ "chosen points", $q_{\min}$ is the minimum over those points, not over two ends.

*Fix.* "(c) … iff, for the set of $\Delta k$; for $R=2$ also for $\tilde D$. With finite support (e.g. the empirical law), (b) is one LP."

### Issue 9 (wording, Prop. E3(c)). The condition on the other elections is unnecessary.

"Whenever $[G^-_0,G^+_0]$ has positive length and at least two other elections have ranges of positive length" — the named corner ($D_1$ at its maximum with $D_{-1}$ at its minimum) needs $g_0=G^-_0$ and $g_0=G^+_0$, impossible when $G^-_0<G^+_0$ whatever the other ranges (check "E3 [W]": gap widths 0, 0.356, 0, corner infeasible). What is needed is $T\ge3$ (two non-base elections). *Fix:* "whenever $T\ge3$ and the base gap range has positive length".

### Issue 10 (wording, Def. 5.2, Prop. E18). Map to Rambachan–Roth; vector form.

Def. 5.2 matches RR's $\Delta^{RM}(\bar M)$: post consecutive differences for $t\ge1$ bounded by $\bar M$ times the largest of the $L$ pre-period consecutive differences, including $\delta_0-\delta_{-1}$, with $\delta_0=0$; 2L branches is right. Two points to state: (i) in RR, $\delta$ is the whole bias of the event-study coefficient, $\beta_t-\tau_t$, which here is $b_t-S_t$ (E14); 2e's $\delta_t:=b_t$ deliberately excludes $S$ (E22). Say so in Def. 5.2, since a reader of RR would otherwise include spillovers in $\delta$. (ii) E18(a) uses Theorem E16 for the vector of post-periods (issue 3(iii)). The union-of-LP form is correct for $L=2$ as well (check "E18 [V] L = 2": four branch LPs equal a grid brute force, with and without rate limits; 2e checked only $L=1$).

### Issue 11 (wording, Prop. E20 "Reading"; §9.3 item 6). "Smoothness is more fragile than relative magnitudes" is not proved.

E20 proves that the base election's width counts twice under SD. The comparison with RM is between different dials ($M$ for SD, $\bar M$ for RM), and E19(c2) shows RM on margins alone is usually vacuous too; in the §9.2 example both are uninformative. *Fix:* replace "more fragile than relative magnitudes" and "Smoothness is worse" by the proved statement: "linear extrapolation of a set-valued pre-trend counts the base election's width twice".

### Issue 12 (wording, §9.2). The row "Relative magnitudes with contextual reinforcement also imposed on the pre-war change" needs CR on both changes.

The printed $[0.135,0.278]$ is reproduced only with CR on the pre-war **and** the post-war change; with CR on the pre-war change alone the set is $[-0.416,0.732]$ (check W). The table's previous RM row says "counts alone", so "also" is ambiguous. *Fix:* "with contextual reinforcement imposed on both changes".

### Issue 13 (wording, N21 narrow variant). "Pre-period totals pinned to ±2 votes" around what?

The script pins each pre-period veterans' total to ±2 around the midpoint of its Duncan–Davis interval. 2e (and NS) do not say. *Fix:* state the centre.

### Issue 14 (wording, §8.1 table). Two counts do not match the script output.

G6: 2e prints 51/300, `verify_phase2e_output.txt` (and my rerun) 53/300. G7: 2e prints 42/200, the output 36/200. *Fix:* update.

### Issue 15 (wording, labels). Results and checks are numbered inconsistently.

Lemma E13 precedes Theorem E12. Checks are tagged with names that differ from the results they test: check "E12" tests Lemma E13, "E13" tests Thm E12, "E3" tests E5/E6/E26, "E23" tests E24, "E24" tests E25, "E25" tests E27. The register rows cite "E1–E11" for checks and results interchangeably. In `verify_phase2e.py`, the comment above the E21 block ("the endpoint-wise formula with $B^L_{s1}$ can be strictly too wide") contradicts its own passing check (the formula equals the LP; 6a Lemma 2.1(b) makes the ends jointly attainable). *Fix:* renumber (E12 ↔ E13) and tag checks by the result tested.

### Issue 16 (wording/scope, Thm E27). A shorter proof; what lies outside.

(i) After the substitution, every column $y_{itc'}$ with $c'\ne c$ appears only in block row $(i,t)$ (a unit column), and with those columns removed each block row is a unit row. Deleting unit rows and columns preserves TU, so TU of the reduced system is TU of the change-row matrix on the column-$c$ variables alone; step 3 is then unnecessary. (ii) The "no cycles" hypothesis is used only to choose the signs $\varepsilon_\pi$ consistently; a cycle whose sign constraints are consistent also works (the triangle $(0,1),(1,2),(0,2)$ with laminar families is TU by exhaustive Ghouila-Houri, check "E27 [V]"). (iii) Three pairs sharing a base (degree 3) gives a matrix that is **not** TU (same check). 2e's conjecture that the class is still integral (no gap in 180 objectives) may be true for these right-hand sides, but it then needs an argument other than TU. No change to the theorem is required.

### Issue 17 (error, minor; last sentence of Prop. E15). "Without a stated comparison group … every trend assumption [is] undefined" is false.

MTR and ND alone use no comparison group and give $\mathrm{ATT}_t\in[0,B_{r^*ct}-B_{r^*c0}]$ (E21 without PS). *Fix:* "Without a stated comparison group, $D_t$, the dial, RM, SD and PS are undefined; MTR and ND alone need none." Def. 5.1 item 1 should then read "the treated group, and the comparison group when the trend assumption uses one".

### Issue 18 (gap, minor; gap register and Prop. G7). Three theorems that "hold without it" need a sentence each.

* **G1** cites "all count-layer results for the period-own populations". That is Lemma E1 with period-specific row margins $a_{irt}$; E1 is stated under (FP). Its proof does not use equal row margins beyond step 1. Say "E1 with $a_{irt}$ in place of $a_{ir}$".
* **G5**: with $N_{it}$ (hence $A_r$) a variable, rates are linear-fractional (2d C1: Charnes–Cooper or a scan), not LP targets. Say so in the "rungs modified" column.
* **G7(c)**: "nondifferential" must mean $E[\Delta\beta\mid X,q]=E[\Delta\beta\mid X]$; identification of $\Delta Q_c$ from the corrected regression needs the linear map $\Delta Q\mapsto\big(\bar q^\top\Delta Q\bar q+\operatorname{tr}(\Delta Q\Sigma_u(q))\big)_q$ to be injective over the design. With $\Sigma_u$ constant it is, under A5 for $\bar q$ (the correction is $c\,\bar q^\top\mathbf 1\mathbf 1^\top\bar q$ with $c=\operatorname{tr}(\Delta Q\Sigma_u)$, and $\Sigma_u\mathbf 1=0$ makes $c$ recoverable); with $\Sigma_u$ varying in $q$ it must be stated.

The naive-fails checks G1–G9 are genuine counterexamples: each computes the naive range or identity on a population where the dropped assumption fails and shows the truth outside it, and each closing device contains the truth in every design (I read the code of G1, G5, G6, G7 and reran all).

### Issue 19 (wording, Prop. E19(a)). Spillover term.

"Without the rate limits, the upper end of $\mathcal A_{RM}$ is $\max_F(D_1+\bar MA)$" holds with $\Lambda=0$; with a spillover allowance add $\Lambda$. One phrase.

---

## 3. Points checked and found correct (no issue)

* **E14** line by line, including the sign of $S$ and $b_0=0$; the estimand is the total effect of the programme on the treated group's rate, relative to the global no-treatment counterfactual; NA is stated; spillover onto the comparison group is the separate term $S_t$.
* **E15**: the set is $[B^L-1,B^U]$, contains 0, is free of $s^*$; with $F$ a union, the union of intervals each containing 0 is still an interval.
* **E16(a), (c), (d)**; **E17(a), (c)**.
* **E18(a)**: "$x\le\bar M\max_ky_k$ iff $x\le\bar My_k$ for some $k$" with one $k$ for all $t$; $A(n)=\max|g_s-g_{s-1}|$ since $\delta_s=g_s-g_0$ for $s<0$.
* **E19**: (a) the inequality and the "iff" ($F$ compact, $D_1$ and $A$ continuous); (b) the vacuity criterion and its sufficient condition (at every $n$ the allowed $\theta$ interval then contains 0 and 1); (c1) monotonicity; (c2) the nesting argument (the ATT range at fixed base and post tables is increasing in $A$) and $\bar A(g_0)\ge w_{-1}/2$; (c3).
* **E20**: coefficients $(1,-2,1)$; width identity.
* **E21(c)**: the ends are jointly attainable for a named row (6a Lemma 2.1(b)) and for the rest; with $\Lambda>B^L_{s^*}$ the formula still holds because $B^L_{r^*c0}\ge0$ absorbs the negative term (check "E21 [V]", Λ = 0.3, 120 cases).
* **E22(a)–(c)**.
* **E4**: (b)'s count forms and slopes; (c)'s two identities derived independently from 2c Def. 1 on brute-forced integer tables (check "E4/E5/E26 [V]").
* **E5(c)–(d)**: the equal-rates witness lies in each unit polytope; the branch argument with $\gamma<1$, $ER\ne0$.
* **E6(a)**: the Kitagawa algebra.
* **E7**; **E8(b)** (Jensen: $L$ convex, $U$ concave; check "E8 [V]"), **(c)**, **(d)** (support function of a random segment with fixed direction).
* **E9(a)**: $g^\Delta\equiv\Delta\kappa$ for $R=2$, $q^\Delta$ affine.
* **E10** algebra; **E11** by [param].
* **E12**: R27's proof uses only R19(a) item 2, an actual-target interval and the union bound.
* **E24**: correct and standard (coverage 0.952 with the truth on the branch boundary, check "E24 [V]").
* **E26**: lattice of step $\kappa_r$ with every point attained and Duncan–Davis ends; integer CR set = lattice points between $NM_\Delta$ and $ER_\Delta$, empty in 3 of 38 brute-forced designs (check "E4/E5/E26 [V]").
* **E27**: the GH signing of step 2 (laminar depth signing per pair, consistent $\varepsilon_\pi$ along paths) and step 3; TU confirmed by exhaustive GH on the theorem's classes.
* **E2**, **E3(a), (b), (d)** (the joint Lemma G condition decomposes into per-election conditions, so the sharp set is the product and depends only on the per-election laws), **Lemma G6**, **G3**'s identity (derived independently), **G2**'s dial identity, **G7(a), (b)**, **G8**'s aggregation.

---

## 4. §11: claimed errors in other documents

All eight items are real or accurately described.

1. PD §7 row "**M-ORD-Δ** (MTR)" exists with that label; the ordering of model-implied changes is descriptive, MTR is about potential outcomes. "M-CEIL-Δ … benchmarked from pre-period variation (relative magnitudes, Rambachan–Roth)" exists; 1b B2 Δ3 calls the calibration a heuristic, as cited.
2. 2d Prop. CL: PS is $\theta\ge B_{n1}$ (observed), (d) "under parallel trends, ATT = DiD" omits NA and no spillover, and "ATT ∈ DiD ⊕ [−M, M] (Rambachan & Roth 2023)" is a bounded-variation dial, not RR's sets, and is outer (rate limits).
3. 2d §11 and 1b B5 use $(\sqrt{m_1(1-m_1)}+\sqrt{m_0(1-m_0)})^2$; 6a Lemma 4.1 restricts Bhatia–Davis to single cells.
4. 2d Δ4(b) says "Never empty, and sign DiD = sign $ER_\Delta$" without $\gamma<1$, $ER_\Delta\ne0$; PLAIN_PHASE2d line 60 says "can never be rejected by the data"; integer CR can be empty (my brute force: 3 of 38 designs).
5. 2d FD4 is outer; status remark, accurately described.
6. NS §3.2 N21 compares with $[-2,2]$, uses $\max_{t\le0}|D_t|$, and calls the ATT set "the descriptive set shifted by the trend dial's set"; all three points stand (N21 numbers reproduced, check N21).
7. NS Part 3A line 879 "L·a·R25 (both periods)" exists; one margin set for a shared composition is right (G7(a)).
8. Status updates, not errors.

---

## 5. What must change before the register rows leave "proof review pending"

| Row | Required | Optional (wording) |
|---|---|---|
| **R11a** | issue 2 (E6(b) emptiness); issue 7 (G′ selector, $V$); issue 8 (E9(c) scope, finite support) | 9, 15 |
| **R11b** | issue 1 (E17(b) criterion); issue 3 (E16(b) integer witness, vector form); issue 4 (population in Def. 5.1; NA in E20/E21); issue 17 (E15 last sentence) | 10 |
| **R11b′** | issue 5 for item (d) (until fixed, (d) stays pending; (a)–(c), E18, E20 can move); issue 12 (worked-example row) | 10, 11, 13, 19 |
| **R11c** | issue 5 (E25 hypotheses and construction); issue 6 (FD mean condition and Lindeberg in E13/E23) | 15 |
| **R11d** | none | 15, 16 |
| **R31** | issue 14 (counts); issue 18 (G1, G5 sentences) | — |
| **R32** | issue 18 (G7(c): nondifferential, rank condition) | — |

Once the required fixes are made, the rows can read "**P; proof reviewed (docs/extension/reviews/PHASE2e_proof_review.md)**", with R11c's sharp-FD-set inference still **O**, CLR versions algorithmic, the degree-3 TU class **C** (and not TU, issue 16), and the gap-register items **O** as 2e already says. The author-flagged open items (inference for the sharp FD set; integer TU beyond paths of pairs; integer causal endpoints; grouped rungs R16; CLR; conditional PT with covariates) are correctly labelled open and none is claimed elsewhere in 2e.

---

## 6. Round 2 confirmation (2e at `claude/did-theory` head `9da2203`)

*Confirmation pass on the author's fixes (2e §12, the issue-by-fix table). `9da2203` was merged into this branch first. Reruns: `verify_phase2e.py` gives **76/76 PASS**, and its output is identical to the committed `verify_phase2e_output.txt`. The new section R covers issues 1, 2, 12 and 14. My `phase2e_review_checks.py` is unchanged and gives **15/15** again. Its [CX] checks test the mathematics behind the old statements, so they now agree with the corrected text: the E17 count instance is exactly 2e's new check E17 [F], and the E6 pure-unit instance falls under the new emptiness criterion.*

### 6.1 Verdict per issue

| # | Fix in 2e | Verdict |
|---|---|---|
| 1 | E17(b) is now an exact criterion | **fixed**. The "only if" is right: ATT ≤ D(n)+M+Λ ≤ max D+M+Λ, with equality only at an argmax with b = −M and S = Λ, and $\mathcal A$ is compact, so the end is attained. The margins-only claim is right as well. Under margins only, every maximiser of $D_t$ has the same aggregates: in R = 2, $B_{s^*}$ is a decreasing function of $B_{r^*}$ in each election; for a named pair or the rest, the maximum of a sum of separately maximised terms forces each term to its extreme (6a Lemma 2.1(b)). Every maximiser therefore minimises $\theta^{PT}$ and every minimiser maximises it, and the chain $\theta^{PT}(n^+)\le\theta^{PT}(n^-)\le1-M-\Lambda$ gives the outer conditions. My check "E17 [V]" (200 margins-only designs) agrees. |
| 2 | E6(b) emptiness: empty iff max{0, D̄⁻} > min{D̄⁺, cap}, with its three cases | **fixed** |
| 3 | E16 is stated in vector form; the integer witness is filled greedily; the case B = 0 is handled; K·K_F | **fixed**. The greedy integer choice is valid because $\theta_tA_{r^*}\le\sum_in_{ir^*ct}$, and both sides are integers. The continuous case is stated as masses. The comparison group's unit masses are now covered. |
| 4 | Def. 5.1 item 0 (population); NA in E20 and E21, with the remark on ND/PS | **fixed** (see new item N2) |
| 5 | E25 restated with (H1)–(H3), allowance boxes and union bound 1 − α_proj − Σα_q; (b) moved to Remark E25.1 | **fixed for the asymptotic statement**. (H1) matches corrected P3.5(b2)(ii), including the singular stackings. Under (H2) the rate limits are dropped under FD and the set is called outer. Under (H3) every actual $q$ is written $\tilde q+E_q$, with a box variable $e_q\in[-t_q,t_q]$. The proof is correct: on $\mathcal E\cap\bigcap_q\{|E_q|\le t_q\}$ the true point is feasible in the branch that holds at the truth at $\beta_0$, and P3.5(b2) steps 1–2 (weak duality on $\mathcal E$) bound every branch on the same event. Treating correlated residuals as separate boxes is conservative and valid. Boxes are counted over all $q$ of all branches and post-periods. The level $1-\alpha_{\rm proj}-\sum_q\alpha_q$ is right. See new items N3 and N4. |
| 6 | Lemma E12 states the FD mean condition, restricts to change targets and lists the Lindeberg condition; E23 updated | **fixed** |
| 7 | Lemma G′: selector of $(w,z)$, $\xi=\sigma(W,\zeta)$; $V$ a function of $W$; $V=X$ in E8 | **fixed** |
| 8 | E9(c): "iff" for the Δk set, and for $\tilde D$ only when R = 2; finite support; $q_{\min}$ over chosen points | **fixed** |
| 9 | E3(c): T ≥ 3 and a base range of positive length | **fixed** |
| 10 | Def. 5.2 maps RR's δ to b − S; E18 cites the vector form of E16 | **fixed** |
| 11 | "Fragile/worse" replaced by "the base width counts twice" | **fixed** |
| 12 | §9.2 row: CR on both changes; the pre-only value [−0.416, 0.732] added | **fixed** |
| 13 | N21 centre stated (DD midpoint) | **fixed** |
| 14 | §8.1 counts are now 53/300 and 36/200; check G [V] compares the note with the run | **fixed** |
| 15 | Renumbered: Lemma E12, Theorem E13; checks retagged; E21 comment corrected | **fixed**, apart from one stale citation (N1) |
| 16 | E27: shorter proof added; consistent cycles noted; degree-3 star downgraded from C to **O** | **fixed**. The consistent-cycle remark holds in general by step 2's argument: each election lies in exactly two pairs with opposite contributions. Downgrading to O is right, because no-gap evidence on 180 objectives cannot support a conjecture whose natural proof route (TU) fails. |
| 17 | E15's last sentence corrected; Def. 5.1 item 1 now conditional | **fixed** |
| 18 | G1 uses E1 with $a_{irt}$; G5 linear-fractional; G7(c) "nondifferential" defined, with injectivity | **fixed**. The injectivity argument is correct. With $\Sigma_u$ constant, the regression identifies $\Delta Q'=\Delta Q+c\mathbf 1\mathbf 1^\top$ under A5 for $\bar q$ (since $\bar q^\top\mathbf 1=1$). Because $u$ sums to zero, $\Sigma_u\mathbf 1=0$, so $\operatorname{tr}(\Delta Q'\Sigma_u)=c$, which recovers $c$ and hence $\Delta Q$. |
| 19 | E19(a): "with Λ = 0; a spillover allowance adds Λ" | **fixed** |

### 6.2 Items introduced or left by the fixes (all wording-level or one sentence)

* **N1 (wording).** The proof of Theorem E13 still cites "[E1; E13; 6a:R27]". It should read "[E1; E12; 6a:R27]" (the lemma).
* **N2 (wording).** The summary bullet in §1.4 ("Rule (Definition 5.1 below)") still says that a causal label requires "a stated treated group and comparison group" and does not mention the population. Align it with the amended Def. 5.1 (items 0 and 1).
* **N3 (wording).** E25 now has items (a) and (c) but no (b) (it became Remark E25.1). Reletter (c) as (b), or say "(b) removed".
* **N4 (gap, one sentence).** E25(a)'s parenthetical says "(or for every $n$ with the Hoeffding radii and a finite-sample estimation event, as in P3.8)". P3.8's event is not a projection event. It relaxes each constraint row and the target, with δ = α_s/(J+2). Two things are therefore needed:
  * The finite-sample version must relax every constraint row of every branch, with $J$ the total number of rows across branches, in place of the projection step.
  * P3.8's Hoeffding radius uses $T_i$ in the simplex; for an FD rung the regressand $\Delta T_i$ lies in a difference of simplices, and the per-unit range doubles.

  *Fix:* state the finite-sample version in those terms, or delete the parenthetical. The asymptotic statement does not depend on it.

None of N1–N4 changes a result.

### 6.3 Register rows: move to "P; proof reviewed (docs/extension/reviews/PHASE2e_proof_review.md)"?

| Row | Move? | Conditions and carried labels |
|---|---|---|
| **R11a** | **yes** | Grouped rungs with model rows stay R16 (O). |
| **R11b** | **yes** | — (apply N2 when convenient) |
| **R11b′** | **yes** | Includes (d) through E25's asymptotic statement. |
| **R11c** | **yes, once N4 is applied** (one sentence, either way; no further review needed) | Inference for the sharp FD set stays **O**; CLR versions are algorithmic (P3.5(b3)). |
| **R11d** | **yes** | The degree-3 class and several columns stay **O**/MILP; causal integer ends stay **O** (Remark 7.1). |
| **R31** | **yes** for the register and Lemma G6 | The items themselves stay **O** (deferred to the uncertain-totals module). |
| **R32** | **yes** for G7(a)–(c) as stated | R25(b)–(c) and the module stay **O**. |

N1 and N3 are labelling fixes and do not block any row.

---

## 7. Round 3: independent-review revisions (2e at `claude/did-theory-r3` head `481c47f`)

*Scope. This round reviews:*
* *the independent econometric review, `PHASE2e_independent_review.md` (findings 1–12 and its Follow-up section), and its round 2, `PHASE2e_independent_review_round2_WPT.md`;*
* *the author's reply;*
* *the revised note: §12 rows IR1–IR12, IR-F/L/A/P/S and R2-1–R2-10, R2-P; the new §§5.8–5.10 (E28–E34); and the revised §9.*

*Branch `481c47f` was merged into this branch first (fast-forward).*

*Reruns.*
* *`verify_phase2e.py`: **102/102 PASS**, output identical to the committed file.*
* *`phase2e_independent_checks.py` and `phase2e_independent_checks_pt.py`: both run to completion, and their printed intervals match the note (four-precinct RM rows, N21, the (99,1)/(1,99) case [.12, .28], the integer RM endpoints, and the six three-precinct LCE fits).*
* *My round-1 script: 15/15 still.*
* *New: `phase2e_review_checks_r3.py` (output `phase2e_review_checks_r3_output.txt`, **4/4 PASS**). It checks only what no passing check already covers:*
  * *E30(c) and ε\* by an independent LP;*
  * *coverage of the E29(b) Hoeffding and E31 empirical-Bernstein radii under precinct-wide shocks, unequal weights and non-identical laws;*
  * *the probability of E32's simultaneous moment event under shocks shared across elections;*
  * *E31(c)'s arithmetic.*

*Maurer–Pontil.* I checked the statement against the paper's text (arXiv 0907.3740, extracted locally):
* Theorem 11 is the non-identically distributed version: independent variables with values in $[0,1]$, bounding $E[P_n(X)]$, the mean of the means.
* $V_n=\frac1{n(n-1)}\sum_{i<j}(X_i-X_j)^2$, the usual divisor $n-1$.
* Theorem 4 is its i.i.d. special case.

E31's use is correct: the common deterministic range is rescaled to $[0,1]$, which multiplies the deviation term by $R$ and $V$ by $R^2$, so $\sqrt{2V(Y)\log(2/\delta)/K}$ and $7R\log(2/\delta)/(3(K-1))$; the reflection $1-X$ has the same sample variance; and $\delta=\alpha/2$ per side gives $\log(4/\alpha)$.

### 7.1 Where the earlier rounds of this review were wrong

The independent review found errors that rounds 1–2 of this review missed or accepted, and in two cases one of my fixes was incomplete. I agree with each of these and say so plainly:
* **E18(b), later horizons** (their finding 1). I accepted the endpoint-wise set for every post-period. Violations accumulate: $|b_h|\le h\bar MA$. Their counterexample has $\mathrm{ATT}_2=-.30$ outside $\pm.15$; it is correct.
* **E23(b), finite sample** (finding 2). I caught the doubled range only for a parenthetical in E25 (my N4) and missed that E23 itself invokes P3.8 with one election's range. The exact coverage .884682 against .997423 is correct.
* **E23(a) and E25 (H1), asymptotic conditions** (finding 3). I accepted "full-rank design" as enough. It does not give the CLT or a nonsingular stacked covariance.
* **E16(b), whole-person witness** (finding 4). My issue 3 fix required integrality for the treated group only. The comparison counterfactual also needs $(B_{s^*}-S)A_{s^*}\in\mathbb Z$ (one treated and one comparison person, $S=1/2$).
* **E16(d), several post-periods** (finding 9). I did not check the joint quantifier; one table must serve every post-period.
* **Remark 7.1** (finding 5). I accepted "beyond one step for the dial" as an open item. Their example (integer set $\{0\}$ against continuous $[-1,1]$, margins only, $M=0$) is correct; I verified the algebra, including feasibility of $y_0=m^2/(2m+1)$.
* **E27, degree-3 star** (finding 6). I showed the matrix is not TU but left integrality conjectural. Their three-row instance ($x_1+x_2\le1$, $x_2+x_3\le1$, $x_1+x_3\le1$: $3/2$ against $1$) is a genuine integrality gap.
* **E7(d)** (finding 7). I accepted the common-composition contrast as the loading. The electorate DiD uses own-group weights; the loading .0075 with all pair contrasts zero is correct.
* **E8, joint stacking** (finding 8). I did not consider several pairs sharing an election. Pairwise sharp sets are only outer for a path; their $x=.4,.5,.6$, $\Delta\kappa=1.25$ example is correct.
* **E19 and §9.3, "only when"** (finding 10). I accepted the applied conclusion that RM helps only when something narrows the pre-trend. The (99,1)/(1,99) case gives $[.12,.28]$ with the counts alone.
* **E3(c), consecutive DiDs** (finding 12). My issue 9 fix ("whatever the other ranges") is right for the base-relative vector, but I left the sentence about consecutive DiDs, which needs an uncertain shared interior election.
* **E21's proof** (finding 9). In round 1 I noticed that $S=\Lambda$ can violate $B_{s^*}-S\ge0$, judged it harmless because the formula survives, and recorded no issue. The proof step was wrong; $S=\min(\Lambda,B_s)$ is the right repair.
* **G5 and G7(a)** (finding 11). My issue 18 said rates over an unknown $N$ are linear-fractional (Charnes–Cooper). A DiD is a *difference* of ratios with different denominators, so a single Charnes–Cooper transformation does not apply; a joint denominator scan is needed. G7(a) also needs the rows to be true.
* **Holes in RM sets** (finding 12). E18(a) always said "union", but I did not flag that endpoint envelopes can include a hole even with one polytope of tables. Their example $[-.90,-.70]\cup[-.30,0]$ is correct.

### 7.2 The fixes to the independent review (IR1–IR12, R2-1–R2-10)

Each fix in §12 is present in the text and is correct. I checked:
* E18(b): the first-post-period form with $\Lambda=0$, the horizon form $\mathcal D_h\oplus[\pm(h\bar M\sup A+\Lambda)]$, and the proof by telescoping.
* E23(a)–(c): design-sequence, leverage and covariance conditions; the doubled range; E12 as a separate event.
* E25 (H1): failures of full rank, and redundant contrasts.
* E16(b), (d).
* Remark 7.1.
* E27: the cycle extension, which holds for every simple cycle by the same signing; the star gap.
* E7(d): the weighted loading; the tree-rank argument under full-dimensional support, which is correct (vanishing edge contrasts force equal rows, and antisymmetry forces them to be zero).
* E8(e): Lemma G′ on levels with an incidence matrix.
* E15: the disconnected-$F$ proof.
* E17(c): the sharp breakdown LP. It is linear, with $|b|\le M$ and $M$ a variable, and it is at least $M^*$ because the sharp set lies inside the Minkowski set.
* E21: $S=\min(\Lambda,B_s)$, with (e) for a shared table.
* E22: scope.
* E19(a): $\bar M>0$; (c1) weak nesting; (c3); (c4); (d).
* E3(c); E6 (both groups in each bin); E10 (density); E12; E24 (fixed true branch).
* G5 and G6 (clean witness); G7(a) ($\varphi_a$, row truth, not a single LP).
* R2-1: the PPT witness $p^0$ is valid. It meets PPT, lies in $[0,1]$ by (b), and the non-veterans' counterfactual equals their factual rate in aggregate.
* R2-2: increments.
* R2-3 to R2-6.

### 7.3 New results, line by line

* **E28.**
  * (a): the identity, the support split and the weights. Correct.
  * (b): rate limits. Correct.
  * (c): the Aumann expectation of the product is the product of the expectations, because the support function of a product is the sum of the support functions. The WPT witness uses constant counterfactual masses; the PPT witness is $p^0$. Both are correct. The set is sharp for continuous masses given a known per-precinct law.
  * (d): correct. It is a special case that needs $a_{iN}=0$ on all of $\mathcal U$; otherwise use the LP.
  * (e): $\Delta m_V-\Delta m_N=q^\Delta(x)+\Delta\kappa$. Correct.
  * (f): one joint LP, with at most $\pi_S\Lambda$ for a fixed fit. Correct.
* **E29.**
  * (a): $E[D_W]-\mathrm{ATT}^{\rm mean}=\sum_{\mathcal U}w_{Vi}\theta_i\in[0,\pi_U]$. Correct.
  * (b): the precinct term has range $4w_{Vi}$ on $\mathcal S$ and $w_{Vi}$ on $\mathcal U$. Independence is across precincts, with any dependence within a precinct. The Hoeffding radius and $\sqrt{8\log(2/\alpha)/n_{\rm eff}}$ (0.992 for 30 equal weights) are correct. The independent-persons version is correct: $a_{ig}$ persons, each with range $2w/a_{ig}$.
  * (c): the event is $\{n^\circ\in F\}\cap\{|D_W-ED_W|\le r_\alpha\}$, with the union bound for belief rows. Correct.
  * (d): RM on increments, within precinct ($(2L)^{|\mathcal S|}$ branches) or after aggregation ($2L$). The zero-box version of E25 is the right one here, since the target and the restriction are both means.
  * (e). Correct.
  * My simulation with precinct-wide shocks: Hoeffding coverage 1.000 at radius 1.07 ($n_{\rm eff}=26$).
* **E30.** (a) is trivial (affine functions form a vector space). In (b), nesting holds and $\epsilon\ge1/2$ makes the rows slack. The breakdown value is one LP and is attained, so strict positivity holds for $\epsilon<\epsilon^*$. (c) is reproduced exactly by my independent LP, including $\epsilon^*=1/600$ and $\epsilon=1/2$ equal to the change-only set.
* **E31.**
  * (a): the event, the independence unit (the precinct), the common deterministic range ($Z_i\in[-2w,2w]$ on $\mathcal S$ and $[0,w]$ on $\mathcal U$, both inside $[-2K\max w,2K\max w]/K$), the rescaling, the two-sided union and the projection over $F$ on $\{n^\circ\in F\}$ are all correct. See issue R3-1 on computation and width.
  * (b): $\sum(Y_i-\bar Y)^2\le\sum Y_i^2$, and $r$ is monotone in $s^2$. Correct.
  * (c): arithmetic 1.41. Correct.
* **E32.** Step 1: $L_{it}\le B_{iVt}\le U_{it}$ a.s., so every element of the identified set has $EL\le m\le EU$ (E8(d)), and $ET=p$. Step 2: the terms $h_iH_{it}\in[0,h_i]$ are independent across precincts, with any dependence across elections; $2\exp(-2c_h^2/\sum h_i^2)=\alpha/J$; union over $J$ rows. Step 3: containment of the identified set, which contains the truth. All correct, and no sharpness is claimed. My simulation with shocks shared across elections gives $P(\text{event})=1.000$ with $J=12$.
* **E33.** (a) and (b) are correct (triangle inequality; opposite violations for the converse). For (c), see issue R3-2.
* **E34.** An affine function has at most one sign change along sorted compositions. Ties and zeros fit either side, and patterns that split a tie are harmless because they force a zero. Exactness holds both ways. The count is $2(|\mathcal S|+1)$. Correct.

### 7.4 Issues (round 3)

* **R3-1 (wording, E31).** E31(a) needs $\inf_F(D_W-r)$ and $\sup_F(D_W+r)$ with $r(n)$ convex in $n$ (a norm of an affine map plus a constant). Both are concave-minimisation problems, not LPs; they need global optimisation (vertex enumeration or a MIQP), or (b) as the practical outer version. Also, $R=4K\max_iw_{Vi}$ grows with the largest weight: with concentrated veterans, E31 is far wider than Hoeffding. In my simulation (unequal weights, $n_{\rm eff}=26$, $K=40$) the mean E31 radius is 4.02, against 1.07 for Hoeffding. *Fix:* one sentence on computation in (a), and one in (c) saying that with unequal weights the additive term scales with $K\max_iw_{Vi}$.
* **R3-2 (gap, minor; E33(c)).** "At $\bar M=0$ … implies aggregated RM at $\bar M=0$" needs the strata to **partition** $\mathcal S$; otherwise the aggregate is not a weighted sum of stratum averages. *Fix:* state that the strata partition $\mathcal S$.
* **R3-3 (wording, §9.3 item 5).** "Without a trend assumption the data say nothing about the effect" overstates E15. The set $[B^L-1,B^U]$ always contains 0, so the sign is never identified, but its width can be well below 2 (in §9.2 it is $[-0.5675,0.7325]$). *Fix:* "…the data cannot sign the effect; they only bound it by the observed post-war rate."

None of R3-1 to R3-3 changes a result. I found no error in the new results.

### 7.5 Verdict per result (this round)

| Result | Verdict |
|---|---|
| E3(c) (consecutive), E6, E7(d), E8(e), E9, E10, E12 | accept |
| E15, E16(b), (d), E17(c), E21, E22, Def. 5.1 | accept |
| E18(a), (b); E19 (a)–(d), (c4) | accept |
| E23(a)–(c), E24 (scope), E25 (H1)–(H3), zero-box clause | accept |
| Remark 7.1, E27 (cycles; star gap) | accept |
| E28 (a)–(f) | accept |
| E29 (a)–(e) | accept |
| E30 (a)–(c) | accept |
| E31 (a)–(c) | accept with fixes (R3-1, wording) |
| E32 (a)–(b) | accept |
| E33 (a)–(b) | accept; (c) accept with fixes (R3-2) |
| E34 | accept |
| G5, G6, G7(a)–(c), register entries G1–G9 | accept |
| §9 (9.1–9.4) | accept with fixes (R3-3); otherwise consistent with the results: the realised/expected labels, the 47.5%/2.5% weights, $n_{\rm eff}$ and the 0.992 allowance, RM on increments and accumulating with the horizon, and the E30 numbers |

### 7.6 Register rows: move to "P; proof reviewed (docs/extension/reviews/PHASE2e_proof_review.md)"?

| Row | Move? | Carried labels |
|---|---|---|
| **R11a** | **yes** | Grouped rungs: R16 (O) |
| **R11b** | **yes**, including E28 | — |
| **R11b′** | **yes**, including E29(d), E33 (after R3-2), E34 | — |
| **R11c** | **yes**, including E29(b)–(c), E31 (R3-1 is wording), E32 | Sharp-FD-set inference **O**; CLR algorithmic; E31/E32 width and coverage simulations remain future work, not a proof gap |
| **R11d** | **yes**, including the cycle extension, the star gap (a negative result) and Remark 7.1 | Causal integer ends: MILP |
| **R31** | **yes** (register and Lemma G6) | Items **O** (deferred) |
| **R32** | **yes** (G7(a)–(c) as revised) | R25(b)–(c), the support of the true compositions and the module **O** |

R3-1 to R3-3 can be applied with the status change; none needs further review.
