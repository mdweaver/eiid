# Round 2: theory plan (Batch B partial identification; Phases 1–3 referee, items 2–5)

*Status (19 Nov 2026): approved by the user in full and implemented on branch `claude/batchB-theory`. Decisions: P3.4 option A and B; projection $q=\dim\beta$; P3.5(c) (c1)/(c2) with the $\bar m$ interval on the rung alone; P3 §8 projection column for $B$ rerun; N2.2′ remark only. The text below is the plan as approved.*

*Plan only (as written). No theory document, script or package file had been changed. Read against `origin/claude/wizardly-volta-lbh4bm` at `9c46fc3`; both referees reviewed `e170be9`. The theory files named below are unchanged between `e170be9` and `9c46fc3`.*

**Abbreviations.** **2b** = `PHASE2b_districts.md`, **2c** = `PHASE2c_monotone.md`, **5a** = `PHASE5a_sign_module.md`, **P3** = `PHASE3_inference.md`. Each row gives:
* **Verdict**: agree / partly agree / dispute, with the reason.
* **Change**: file, location and proposed wording.
* **Check**: a new numerical check, if one is needed.
* **Package impact.**
* **Effort**: S = under a day; M = one to three sessions; L = more.

**Checks run for this plan.**
* The PI referee's three-row equal-rates example reproduces: ER = 1, γ = .24, D = NM = .24, δ_W = .0456.
* In the Phases 1–3 referee's D8 counterexample (2×2, $Q_{1,12}=0$, $A$ = simplex), the vertex rows force $k=0$, as the referee says.

---

## I. Batch B partial-identification referee, round 2 (recommendation: minor revision)

### Corrections the referee requires before acceptance

| # | Point | Verdict | Change | Check | Package | Effort |
|---|---|---|---|---|---|---|
| I.1 | The equal-rates witness shows the CR **union** is feasible, not both branches (2c l.161). | **Agree.** At the witness $\delta_B=0$, but $\delta_W$ has any sign. In 2×2 this is consistent with M3(d): the witness sits at $D_{NM}$, inside the one nonempty branch. The error is only in the R×C sentence. | 2c "Reading", l.161: replace the sentence with the referee's wording: "The equal-rates table ($\beta_{ir}=T_i$ for every row) sets every pairwise between association to zero. It therefore satisfies every weak between sign and, for each specified pairing, the CR union, in the continuous problem without unit restrictions. Which branch contains it depends on its within association; both do only if that is also zero. Example: equal weights on $(.2,.5,.3)$, $(.5,.2,.3)$, $(.8,.1,.1)$ with every rate equal to $X_1$: $ER=1$, $\gamma=.24$, the witness has $D=NM=.24$ and $\delta_W=.0456>0$, so it lies in the positive branch only." Then grep the other copies of this claim and fix any that say "both": PLG §9, the vignette, 5a R3(e)/§1.6, and the EI response l.168. The response letters belong to the coordinator. | **C15** in `verify_phase2c.py`: the three-row table; the witness satisfies $\delta_B=0$ and the CR union; the negative branch is empty at the witness. It fails on the old "both weak CR branches hold". | None, unless a package string says "both branches" (grep `R/` and the vignette). | S |
| I.2 | The scope of R10 contradicts the response: it says "Pooled M0, or M1/M2 in regime R1". | **Agree.** Grouped R3/R10 are R16 (not proved). | 5a Theorem R10, first paragraph: delete "or M1/M2 in regime R1". Add: "The fixed-district baseline intervals of P3 (M1/M2 in R1) are background; the added-row, slab and union statements for district rungs are proof item R16." | None (a scope statement). | None; the print already says "outer … not yet proved" for grouped rungs. | S |
| I.3 | An exact rank reduction is not the same as a numerical diagnostic. The response calls "drop candidates using a numerical null direction" case (ii-b). | **Agree.** (ii-b) needs a fixed $L$ and the exact identities $\beta=L\eta+b$, $\hat\beta=L\hat\eta+b$. A numerical rank cutoff establishes neither, and with $q=\operatorname{rank}_{\rm num}\hat\Sigma$ the $\chi^2_q$ critical value is not calibrated. Dropping candidates can only widen $\hat U$ (a minimum over fewer terms), so it is a safeguard. | P3 §5, after (b2)(ii) (l.179–193): add "*Remark (exact versus numerical reduction).* Case (ii-b) is claimed only for a documented exact reduction, e.g. the adding-up identity that removes column $C$ (Remark P3.4.1), or coordinates fixed by design. When the rank comes from a numerical cutoff on $\hat\Sigma$, the theorem applies only with $q=\dim\beta$ (case (ii-a), assuming the true $\Sigma$ is positive definite). Dropping candidates whose gradients load on a numerically null direction is then a **safeguard**: it can only widen the interval, but it does not calibrate a rank-$q$ critical value." P3 §8 l.379 (rank paragraph): same wording. | Optional: extend the existing "null direction" check in `verify_phase3.py` with a numerically singular $\hat\Sigma$ whose true $\Sigma$ is PD. Using $q=\operatorname{rank}_{\rm num}$ under-covers; using $q=\dim\beta$ does not. | **Yes, small.** `.rxc_projection` (`R/rxc_infer.R:265`) sets `q` to the numerical rank. Recommend `q = dim(beta)` (the $C-1$-column parameterisation is already the exact reduction), keeping the drop rule as a labelled safeguard. WS6 found no rank drop, so no number changes in practice. The response letter's wording needs correcting (coordinator). | S |
| I.4 | N2.3(d) needs the rows to be nondecreasing **on a neighbourhood of x**. N2's monotonicity on the finite support does not give a derivative sign. | **Agree.** Both rows equal to $.5+h$ with $h=(x-\tfrac14)(x-\tfrac12)(x-\tfrac34)$ are constant on the support, yet have derivative $-1/16$ at $1/2$. | 2c N2.3(d): change to "If $m_{1c}$ and $m_{2c}$ are nondecreasing **on a neighbourhood of $x$** (monotonicity on the finite support is not enough), then the local sign is $\ge0$ at $x$." Add the example in one line. | **C16**: the $.5+h$ rows; monotone (constant) on the support, but the local within effect at $1/2$ is $-1/16<0$. It fails on the old premise read as N2 monotonicity. | None. | S |

### Points in the referee's disposition table

| # | Point | Verdict | Change | Check | Package | Effort |
|---|---|---|---|---|---|---|
| I.5 | Lattice: a numerical value of the step needs **both** group totals ($1/N_1+1/N_0$). | **Agree.** Our "about $6\times10^{-6}$ at 370 units of 500–4,000 people" gives no group totals. | 2c Remark M2.1: "The step is $1/N_1+1/N_0=1/(N\mu(1-\mu))$; for example, with $N=10^6$ and $\mu=.3$ it is $4.8\times10^{-6}$." Use a stated example, or the Wisconsin totals if we want a real number. | None. | None. | S |
| I.6 | Binning: the equality uses the **model-implied** bin mean $E_b\tau$; the observed person-weighted bin mean of $T$ estimates it. A prespecified degree still needs correct specification or a bias allowance. | **Agree.** 2c §6 says "$\bar\tau_b$ is then the person-weighted bin mean of $T$", which conflates the two. | 2c Lemma N2.4: define $\bar\tau_b:=E_b[\tau(X)]$ (model-implied); note that $E_b[T]$ equals it in the population and that the sample bin mean estimates it with residual error, which needs inference. 2c §6, degree bullet: "valid if $\tau$ is a polynomial of degree $q$; otherwise add an approximation-bias allowance, or label the bound conditional on the specification." | None. | Step 8 design note (coordinator). | S |
| I.7 | Line 30: narrow "for them every computation is a finite LP" to the regression-only optimisation over $\Theta_O\cap C$. A polyhedral $C$ does not give a finite representation of the full-law $\Theta_I$. | **Agree.** | 5a §1.1, l.30: "…so for them the regression-only computation, optimisation over $\Theta_O\cap C$, is a finite LP. The full-law set $\Theta_I\cap C$ of R3(a) need not have a finite representation even for polyhedral $C$." | None. | None. | S |
| I.8 | M4 example: complementary columns have $\delta_{W,2}=-\delta_{W,1}$. **The same** weak sign on both forces both to zero ($D=ER$); opposite signs duplicate one restriction. | **Agree.** Our "opposite within signs … jointly force $D=ER$" has it backwards. Their own appendix code imposes `weights` and `-weights` on the two complementary columns, i.e. the same sign on both. | 2c §2, M4 paragraph: "…imposing the **same** weak within sign on both columns forces both associations to zero, because the columns are complements ($\delta_{W,2}=-\delta_{W,1}$), and so forces $D=ER$, shrinking $[-.5,0]$ to $\{0\}$. Opposite signs on complementary columns duplicate one restriction." Also check the vignette and PLG: 2×3 joint-sign language near the P-E27 adding-up sentence. | Optional **C3b**: the referee's two-unit example, same sign → $\{0\}$; opposite signs → the single-restriction interval. | None. | S |
| I.9 | M3(d): distinguish a **branch** (a set of latent tables) from its singleton gap **image**. | **Agree.** We did this for $ER=0$ but still write "both branches equal $\{D_{ER}\}$" at $\gamma=1$. | 2c M3(d), $\gamma=1$ bullet: "…both branches are nonempty and have the same gap image $\{D_{ER}\}$". In the second bullet, "both branches give $\{0\}$" becomes "both branches have gap image $\{0\}$". Proof step 6 likewise. | None (C8 already checks the images). | None. | S |
| I.10 | M2(d): distinguish the **set** of slack values that preserves a conclusion from its **frontier** (its boundary). | **Agree.** Our wording calls the set a frontier. | 2c M2(d): "…the set of $\varepsilon$ that preserves a given conclusion (a sublevel set of a piecewise-linear function) is in general not a product of intervals; its boundary is the **breakdown frontier**." Check C7b's label: "the preserving set is not a product". | None (relabel C7b). | None. | S |
| I.11 | N2.2′: a short sufficiency proof exists under the **selection-expectation** reading. For each positive-mass $(x,n)$, take a measurable selection with conditional mean equal to the same $m(x)$, and glue the selections across the finitely many atoms. This does not imply that an $E[T\mid X]$-only procedure achieves it. | **Partly agree / adopt as a remark only** (the user's instruction is to leave N2.2′ open). The argument looks right on the finite support, but it shows sufficiency of the $V$-conditional characterisation, not that any implemented procedure is sharp. | 2c Remark N2.2.1: add "*Remark (sufficiency under the selection-expectation reading; due to the referee, round 2).* If $\mathbb E[S(X,T)\mid V]$ denotes the selection (Aumann) expectation, every $m$ with $m(x)\in\mathbb E[S(X,T)\mid V=(x,n)]$ for all positive-mass $(x,n)$, plus the shape restrictions, is attained. For each atom take a measurable selection of $S(x,T)$ with conditional mean $m(x)$, and glue across the finitely many atoms; accounting holds a.s. and the shape restrictions hold for the common $m$. Conjecture N2.2′ stays open as stated: necessity in general, and in particular whether any procedure that uses only $E[T\mid X]$ attains it (Remark N2.2.1's example shows it does not when size is informative)." Keep the label "Conjecture". | None. | None. | S |
| I.12 | "Each new check fails on the earlier statement" overstates. V10 checks the old containment claim, V11 checks lengths, and C13 gives a smooth witness for a non-implication that was already true. | **Agree.** The plan itself marked C13 and V9 as "nothing fails". | Classify every Batch-B check in three classes: **(F)** fails on an earlier statement (C7, C8, C8b, C9, C10, C10b, C11, C12, C12b, C14, T1c, T4b, S7, V10, [supp], P3 null direction, P3 residual CLT); **(V)** verifies a corrected or new statement (C7a, C7b, C5(c′), S6, V9, V11, C4 rebuilt); **(I)** illustration or strengthening (C13, V12). Edit the 2c header ("Every correction has a check that fails…") and the 5a header to say this. The response letters belong to the coordinator. | Add a tag `[F]`, `[V]` or `[I]` to each check label in the scripts and regenerate the outputs (labels only; numbers unchanged). | None. | S |
| I.13 | R2, R3 (pooled), R10 (restricted scope), district slopes, finite support and derivatives, and the background repairs are accepted. | — | No action beyond I.1–I.4. | — | — | — |

---

## II. Phases 1–3 referee, round 2: theory items 2–5

### Item 2: Theorem P3.4 with changing sample hulls (referee's priority 1)

**Verdict: agree.** The proof is correct for a fixed finite set of affine maps. "Conditional on the design with the sample hull fixed" is not a coherent asymptotic statement once the hull grows with $n$. Near-ties at gap $O(n^{-1/2})$ change the limit. In the referee's example, $f_n(\theta)=\min\{\theta,-\theta+n^{-1/2}\}$ gives $\min\{Z,1-Z\}$, not normal. Only P3.4 and the remarks that rely on it (the "plug-in too narrow" reading) are affected. P3.5(a), (b2), P3.8 and P3.6 do not use P3.4.

**Recommendation: both.** Adopt **A** now as the theorem, and add **B** as a separate proposition. B costs little more than A and covers the realistic case, because it is the form the referee's own counterexample takes.

* **A (Theorem P3.4, fixed constraint system).** Restate the hypothesis as "a fixed finite constraint system": a fixed vertex set $V$ (for example a prespecified polytope $A$ with fixed vertices), or a design sequence along which the sample hull is eventually constant, $\hat A_n=A$ for $n\ge n_0$. The latter holds, for instance, when $A$ is a polytope whose vertices are atoms of the design law. Then $\Lambda$, $\Lambda^*$ and the $g_\lambda$ do not depend on $n$, and the existing proof is unchanged. Delete "with the design and the sample hull $\hat A$ fixed, so that…" and add: "*(Corrected after review, round 2: the theorem does not cover expanding hulls; see P3.4′.)*"
* **B (Proposition P3.4′, triangular arrays).** For each $n$, let $\Lambda_n$ be the dual vertex set of the design-$n$ system, with $\lvert\Lambda_n\rvert\le\bar L$ uniformly, gradients $g_{\lambda,n}$ and localised slacks
  $$c_{\lambda,n}:=\sqrt n\big(\theta_{\lambda,n}(\beta_0)-\varphi_{U,n}(\beta_0)\big)\ \ge0.$$
  **Condition (localisation).** After relabelling, $g_{\lambda,n}\to g_\lambda$ and $c_{\lambda,n}\to c_\lambda\in[0,\infty]$ for each $\lambda$, along the sequence or along subsequences.
  **Conclusion.**
  $$\sqrt n\big(\hat\varphi_{U,n}-\varphi_{U,n}(\beta_0)\big)\Rightarrow\min_{\lambda:\,c_\lambda<\infty}\big(g_\lambda^\top\mathcal Z+c_\lambda\big).$$
  *Proof sketch.* $\sqrt n(\hat\varphi_U-\varphi_U)=\min_\lambda\big(g_{\lambda,n}^\top\sqrt n(\hat\beta-\beta_0)+c_{\lambda,n}\big)$ exactly, by Lemma P3.1. Candidates with $c_{\lambda,n}\to\infty$ drop out with probability $\to1$, because $\sqrt n(\hat\beta-\beta_0)=O_p(1)$ and the gradients are bounded. The rest follows by continuous mapping.
  * The **separation condition** under which A's normal or min-of-normals limit holds is $c_\lambda\in\{0,\infty\}$ for every $\lambda$: each candidate is either exactly active or inactive by more than $n^{-1/2}$, i.e. $\min_{\lambda\notin\Lambda^*_n}\{\theta_{\lambda,n}(\beta_0)-\varphi_{U,n}(\beta_0)\}\cdot\sqrt n\to\infty$, with a stable active set and convergent gradients.
  * The referee's example has $c=1$ for the second piece and gives $\min\{Z,1-Z\}$.
  * State that B is a limit for the sample-hull endpoint along the design sequence. It does not address the population-hull error, which stays as in the P3.2 remark.

**Other P3 changes.** Remark after P3.2 ("Theorem P3.4's limit is stated for the sample-hull endpoint…"): point to A/B. §8 l.436 ("fails where Theorem P3.4 predicts"): add "via P3.4′ when near-ties are of order $n^{-1/2}$". §10: add "uniform (in the design sequence) limit theory for the plug-in endpoint" as open.

* **Check.** `verify_phase3.py` gets **P3.4′**: simulate the referee's $f_n$ at $n=10^4$; the limit is $\min\{Z,1-Z\}$, not normal (mean $\approx-0.17$ rather than 0; KS distance to $N(0,1)$ large). This fails on the old reading, under which a uniquely active first piece gives a normal limit, and matches B's formula with $c=1$.
* **Package impact.** None; the CIs do not use P3.4.
* **Effort.** A is S. B is M (statement, proof and check).

### Item 3: projection interval for $B$ (referee's priority 2)

**Is P3 §8 stale relative to 0.3.0? Yes.**
* P3 §8 l.377 says the projection interval for $B$ takes its Bhatia–Davis (BD) variance bound "from a CLR interval (as in the package)".
* In 0.3.0, `.rxc_fpc_sd` (`R/rxc_infer.R:303–330`) returns two bounds:
  * the bound used by `ci_actual` comes from the CLR interval for $\bar m=\sum_i\pi_im_i$, with $\pi_i=w_i^2/\sum w^2$ (with `ci = "conservative"`, from the conservative interval);
  * the bound used by `proj_actual` comes from the **projection** interval for $\bar m$ (`attr(iv, "projection")`) at level $\alpha_v$.
* So the referee's point does **not** apply to `proj_actual` as implemented. It still applies to `ci_actual` under the CLR default, which is algorithmic anyway, and to §8's text. §8 l.377 should say: "In the package (0.3.0) the projection interval for $B$ takes its BD bound from the projection interval for $\bar m$, so it rests on P3.5(b2)(ii) and the new P3.5(c) below; the §8 runs predate that change and used the CLR-based bound."
* §8 is WS6 text. Changing it needs the coordinator's say-so; I recommend this one-sentence correction plus a pointer.

**Is a BD bound taken from the projection interval covered by the simultaneous event? Yes, with one extra condition.**
* *Coverage.* The projection event $\mathcal E_v=\{n\,\hat e^\top\hat\Sigma^{-1}\hat e\le\chi^2_{q,1-\alpha_v}\}$ covers every linear functional at once, including the model-implied $\bar m$, a linear target in $(\beta,k)$. On $\mathcal E_v$, $\bar m\in\mathrm{CI}_v$. Concavity then gives $s_n^2=\sum_iw_i^2m_i(1-m_i)\le\sum w^2\,\bar m(1-\bar m)\le\hat S_{BD}^2:=\sum w^2\max_{m\in\mathrm{CI}_v}m(1-m)$.
  * Then $\{\lvert E\rvert>z\hat S_{BD}\}\subseteq\{\lvert E\rvert>z\,\bar s\}\cup\mathcal E_v^c$, with $\bar s=(\sum w^2\,\bar m(1-\bar m))^{1/2}$ deterministic.
  * The total level is $\alpha_{\rm est}+\alpha_E+\alpha_v$ by the union bound, which is what the package's `alpha_split` spends.
  * If the same ellipsoid level were used for both the estimation and the $\bar m$ event, $\alpha_v$ could be saved. That is optional.
* *The extra condition.* The residual step must now hold relative to $\bar s$, not $S_1$. When $\bar s/S_1\to0$ (e.g. $\bar m\to0$), the Lindeberg branch of P3.5(a) step 5 is no longer available relative to $\bar s$; Chebyshev gives only $P(\lvert E\rvert>z\bar s)\le s_n^2/(z^2\bar s^2)$, which need not be $\le\alpha_E$. So state a condition: **$\max_iw_i^2/\bar s^2\to0$** (Lindeberg relative to $\bar s$), e.g. $\bar m$ bounded away from 0 and 1. Under it, either $s_n/\bar s\to0$ (Chebyshev; probability $\to0$), or $s_n/\bar s$ is bounded below and Lindeberg applies relative to $s_n$.

**Theory text for the package session to implement against** (new P3.5(c), "Variance-adaptive residual allowance").
* Let $\mathrm{CI}_v$ be a confidence set for $\bar m$ at level $1-\alpha_v$, and $\hat S_{BD}$ as above.
* (c1) *Valid version.* If $\mathrm{CI}_v$ comes from P3.5(a) or (b2)(ii), and $\max_iw_i^2/\bar s^2\to0$, then replacing $z_{\alpha_E}S_1$ by $z_{\alpha_E}\hat S_{BD}$ in (a), or in the projection interval for $B$, gives $\liminf P(B\in\mathrm{CI})\ge1-\alpha_{\rm est}-\alpha_E-\alpha_v$. The Popoviciu version ($S_1$, no $\alpha_v$) needs only $\max w^2/\sum w^2\to0$.
* (c2) *Adaptive version.* If $\mathrm{CI}_v$ is the CLR interval (Conjecture P3.5(b3)), coverage holds conditional on $\mathrm{CI}_v$ covering $\bar m$. The result inherits the algorithmic status.

The package should label `proj_actual` as resting on (b2)(ii) + (c1). It should also fall back to Popoviciu, or print a note, when $\mathrm{CI}_v$ reaches within a small distance of 0 or 1 (where the Lindeberg condition is doubtful). `ci_actual` under the CLR default is (c2): algorithmic, as already labelled. $\tilde B$ is unaffected.

* **Check.** `verify_phase3.py` gets **P3.5(c)**:
  * (i) the referee's illustration: true $m=.5$ while a wrong interval in $[0,.1]$ gives a variance bound $\le.09<.25$, so BD is not conservative if the interval fails;
  * (ii) a small simulation in which BD from a covering interval keeps the allowance valid.
* **Package impact.** Labels and the optional near-boundary fallback only. The numbers are unchanged, since `proj_actual` already uses the projection interval. This goes to the package session.
* **Effort.** M (theory S; §8 sentence S; package S).

### Item 4: Proposition D8 (2b) is sufficient, not necessary (priority 3)

**Verdict: agree.** "Point-identified, whatever the identified set for $K$, iff $\Delta q=0$" is true only when read as "invariant to every $K\in\mathcal K$". As written, "point-identified … iff" is false: with $Q_{1,12}=0$ and LCE on the whole simplex, the vertex rows force $K_{1,12}=0$, so the change is identified even with $\Delta q\ne0$. Step 4's own gloss already hedges this.

**Change** (2b Prop. D8, statement and step 4):
* "(a) **Identification by cancellation.** The change does not depend on $K$ for every $K\in\mathcal K$ iff $\Delta q_{rs}:=\mu^{[\rm post]}_{rs}/\mu^{[\rm post]}_r-\mu^{[\rm pre]}_{rs}/\mu^{[\rm pre]}_r=0$ for all $s\ne r$ (a moment equality).
* (b) **Identification given the feasible set.** The change is point-identified iff $\max_{K\in\Theta}L(K)-\min_{K\in\Theta}L(K)=0$, where $L(K)=\sum_{s\ne r}K_{c,rs}\Delta q_{rs}$ and $\Theta$ is the identified set for $K$. This holds under (a), and also when $\Theta$ fixes $L$ (example: 2×2, $Q_{1,12}=0$, $A$ = simplex, so $\Theta_O=\{0\}$)."

Reading: replace "The DiD of transition rates is point-identified…" with "The change in ecological transition rates between periods is point-identified…". Add: "This is a statement about the difference of two ecological rates; a causal (difference-in-differences) reading needs further structure (2d)." Grep 2d, 1b §B and PD §7 for "iff" or "point-identified exactly when" built on D8.

* **Check.** `verify_phase2b.py` gets **D8b**: the referee's 2×2 example ($\Delta q\ne0$; $\Theta_O=\{0\}$; the LP range of the change has width 0). It fails on the old "iff". It also gets a case with $\Delta q\ne0$ and nondegenerate $\Theta$ (positive width).
* **Package impact.** None in code. Check that the DiD print or vignette text does not say "iff equal composition moments".
* **Effort.** S.

### Item 5: P3 §0.4, independence is substantive (priority 4)

**Verdict: agree.** §0.4 argues from "every unit is observed" toward not needing dependence, then states P3-S1a as "what remains assumed". The referee wants the two arguments separated. Within-district shocks with $\operatorname{Cov}(e_i,e_j\mid\text{design})\ne0$ and no district mean shift violate P3-S1a without violating the pooled mean model, so P3.6(c) does not test district dependence.

**Change** (P3 §0.4; P3.6(c) remark):
* Add a paragraph: "Two arguments must be kept apart. (1) Conditioning on observed district effects and **assuming** that the remaining residuals are independent (P3-S1a) is a coherent model. (2) 'Every unit is observed, so dependence need not be considered' does not follow: observing all units removes sampling uncertainty about which units exist, not dependence among their residuals.
* Within-district shocks can make $\operatorname{Cov}(e_i,e_j\mid\text{design})\ne0$ with no district mean shift. That violates P3-S1a but not the mean model, so the pooling test P3.6(c), which tests mean shifts, is **not** a test of within-district dependence.
* P3-S1a is therefore a substantive assumption. Its plausibility should be argued, and residual correlation within districts and between neighbours can be inspected as a diagnostic.
* We do not reintroduce clustered standard errors, for the reasons above. Under dependence, the variance and residual-allowance statements of P3.5 and P3.8 are not covered."
* Add a sentence after P3.6(c): "(c) tests mean shifts only; it has no power against within-district residual dependence."

* **Check.** None needed. Optional: a two-line simulation showing that P3.6(c) has nominal size under equicorrelated within-district residuals with zero mean shift.
* **Package impact.** Vignette or roxygen wording for the pooling test, if they suggest it detects district dependence (grep "pooling test").
* **Effort.** S.

### Item 7 (no action asked)

The CLR status is accepted. The open question, whether simulated candidate collection supports a uniformly valid critical value, is already Conjecture P3.5(b3) and a P3 §10 item. Add "uniformly" to that §10 line.

---

## III. Items needing a user decision

1. **P3.4: A only, or A and B?** I recommend both: A now (S), B as P3.4′ (M).
2. **Package projection rank (I.3).** Switch `q` from numerical rank to $\dim\beta$, or keep the rank and relabel it as a safeguard? I recommend $\dim\beta$; it is conservative and no observed design changes.
3. **P3.5(c) near-boundary fallback for `proj_actual` (item 3).** Fall back to Popoviciu when $\mathrm{CI}_v$ is within $\epsilon$ of 0 or 1 (choose $\epsilon$), or only print a note? I recommend a note for 0.3.x and the fallback later.
4. **P3 §8 edit (WS6 text).** May the one-sentence 0.3.0 correction at l.377 and the rank-paragraph wording at l.379 be made?
5. **N2.2′.** As instructed, a crediting remark only, with the conjecture kept open. Confirm.

## IV. Order and effort

* **Day 1 (all S):** I.1–I.12, item 4 (D8), item 5 (§0.4), P3.4 option A.
* **Then:**
  * P3.4′ (M);
  * P3.5(c) theory and check (S–M);
  * §8 sentences, if approved.
* **New checks:**
  * C15, C16 (and optional C3b) in `verify_phase2c.py`;
  * D8b in `verify_phase2b.py`;
  * P3.4′ and P3.5(c) (plus the optional rank case) in `verify_phase3.py`;
  * relabelled [F]/[V]/[I] tags across the Batch-B scripts.

  Expected totals: 2c 26 → 28 (29 with C3b); P3 15 → 17 (18); 2b +1.
* **Hand-off to the coordinator:**
  * the response letters (I.1 copies, I.3 wording, I.12 claims);
  * PD/NEXT_BUILD entries for P3.4′ and P3.5(c).
* **Hand-off to the package session:** the `q` choice, the `proj_actual` labels and fallback, and grep-and-fix of "both branches", "iff" and "pooling test" strings.
