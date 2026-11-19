# Response to the referee, second round: Phases 1, 2a, 2b and 3

**Reviewed version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b).
**Revised version:** commit [`1b8b6fb`](https://github.com/mdweaver/eiid/tree/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f) on branch `claude/wizardly-volta-lbh4bm`.
**All changes since the reviewed version:** [compare `e170be9…1b8b6fb`](https://github.com/mdweaver/eiid/compare/e170be9161b9bbc250a17a27ffd39bc523c1c88b...1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f).
**Your report, as received:** [`reviews/A_phases1to3_report_round2.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/reviews/A_phases1to3_report_round2.md).

Abbreviations: **P3** is [`docs/extension/PHASE3_inference.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE3_inference.md), **2b** is [`PHASE2b_districts.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE2b_districts.md), **PD** is `PACKAGE_DESIGN.md`. Each corrected statement carries a note "*(Corrected after review, round 2, Nov 2026 …)*" that says what the earlier text claimed. Tags: **[Theory]** a result or proof, **[Code]** the R package (now 0.3.1), **[Docs]** guide, vignette or design notes, **[Check]** a numerical check.

---

## 1. Summary

Thank you for a second careful report and for the favourable recommendation. We accepted all four priorities without dispute. Each now has a stated result, a check, and, where the package is affected, a code change:

1. **P3.4 with changing sample hulls.** P3.4 is now a fixed-constraint theorem. A new Proposition P3.4′ gives the triangular-array limit under a localisation condition; your $f_n$ is a special case of it.
2. **The projection interval for $B$.** A new Proposition P3.5(c) separates a provably valid variance-adaptive allowance, (c1), from a CLR-dependent one, (c2). The package's `proj_actual` uses (c1); the CLR `ci_actual` uses (c2) and is labelled algorithmic. The Monte Carlo column you questioned was rerun; coverage is unchanged in every design.
3. **D8.** It now states the sufficient "cancellation" condition and the necessary and sufficient condition given the identified set separately. We dropped the causal DiD reading.
4. **Independence.** P3 §0.4 now presents it as a substantive assumption that the pooling test cannot check.

One clarification on item 2. The package had already moved off the CLR-based bound in 0.3.0. The text you read in P3 §8 described the Monte Carlo script, which still used it, and wrongly said "as in the package". Your point was correct for that text and that script, and both are now fixed.

## 2. What changed

| Where | Change | Your item |
|---|---|---|
| P3 §4, Theorem P3.4 | Restated for a fixed constraint system (fixed vertices, or $\hat A_n=A$ eventually); correction note with your $f_n$ | 2 |
| P3 §4, **Proposition P3.4′** (new) | Triangular-array limit $\min_{\lambda:c_\lambda<\infty}(g_\lambda^\top Z+c_\lambda)$ under bounded $\lvert\Lambda_n\rvert$, convergent gradients and localised slacks $c_{\lambda,n}\to c_\lambda\in[0,\infty]$; separation condition $c_\lambda\in\{0,\infty\}$ | 2 |
| P3 §5, **Proposition P3.5(c)** (new) | (c1) valid BD allowance from an interval with established coverage on the rung's own system, with the condition $\max_i w_i^2/\bar s^2\to0$; (c2) the CLR version, conditional | 3 |
| P3 §5, remark after (b2)(ii) | Exact versus numerical rank reduction (raised by the Batch B referee in the same round) | 3 |
| P3 §8, §9 | The stale "as in the package" sentence corrected; projection column for $B$ rerun under (c1); recipe item 3 updated | 3 |
| 2b, Proposition D8 | Split into (a) cancellation and (b) given the identified set; your example; "ecological transition rates", not causal DiD | 4 |
| P3 §0.4 and after P3.6(c) | Independence as a substantive assumption; the pooling test does not test it | 5 |
| P3 §10 | Uniform limit theory for the plug-in listed as open; "uniformly" added to the (b3) question | 2, 7 |
| PD §8 | New register items R17 (P3.4′, proved) and R18 (P3.5(c): (c1) proved, (c2) inherits (b3)) | 2, 3 |
| Package 0.3.1 | `proj_actual` per (c1); critical value $q=\dim\beta$; boundary note; `print()` states which result covers each interval | 3 |
| Checks | `verify_phase3.py` 15 → 17 (P3.4′, P3.5(c)); `verify_phase2b.py` 21 → 22 (D8b). All pass | 2–4 |

## 3. Point-by-point

### Comment 2: Theorem P3.4 and changing sample hulls

> *Clarify whether Theorem P3.4 is intended as (A) a fixed-constraint-system result or (B) a triangular-array result … I would favor explicitly proving the fixed-constraint result first and treating inference under changing active sets as a separate extension.*

Agreed. We did both, in the order you suggest.

* **[Theory] (A).** Theorem P3.4 is now stated for a fixed constraint system: the polytope $A$ is prespecified with fixed vertices, or the sample hull equals $A$ from some $n$ on (for example, when the vertices are atoms of the design law). The proof is unchanged. A correction note records the earlier wording and your counterexample (P3 §4).
* **[Theory] (B), Proposition P3.4′.**
  * **Setup:** the dual vertex sets satisfy $\lvert\Lambda_n\rvert\le\bar L$, the gradients $g_{\lambda,n}$ are bounded, and the localised slacks are $c_{\lambda,n}=\sqrt n\,(\theta_{\lambda,n}(\beta_0)-\varphi_{U,n}(\beta_0))\ge0$.
  * **Conclusion:** if, along subsequences, $g_{\lambda,n}\to g_\lambda$ and $c_{\lambda,n}\to c_\lambda\in[0,\infty]$, then $\sqrt n(\hat\varphi_U-\varphi_U)\Rightarrow\min_{\lambda:\,c_\lambda<\infty}(g_\lambda^\top Z+c_\lambda)$.
  * **Proof:** an exact identity (Lemma P3.1), the $O_p(1)$ estimation error, which removes candidates with $c=\infty$, and continuous mapping.
  * **Separation condition:** your suggested condition is $c_\lambda\in\{0,\infty\}$, meaning every candidate is either exactly active or inactive by more than $n^{-1/2}$. Under it the limit reduces to P3.4's normal or min-of-normals form. Your example has $c=1$.
* **[Theory]** The proposition is about the sample-hull endpoint. The population-hull error is still as described in the P3.2 remark, which now points to P3.4′. Uniform limit theory for the plug-in is listed as open in §10.
* **[Check]** `verify_phase3.py`, "P3.4′ [F]": your $f_n$ at $n=10^4$. The simulated limit has mean about $-0.40$ and Kolmogorov distance 0.31 from $N(0,1)$, matching $\min\{Z,1-Z\}$ and not the normal limit the old reading gives.
* As you note, none of the coverage results in P3 §5 (P3.5(a), (b2), P3.8) or the tests in P3.6 use P3.4. They are unaffected.

### Comment 3: The projection interval for $B$

> *Distinguish two implementations: provably valid … and variance-adaptive … The package should distinguish the two versions accordingly.*

Agreed. This is now Proposition P3.5(c). On the package, the fact we owe you is that 0.3.0 had already taken the Bhatia–Davis bound for `proj_actual` from a projection interval, which is your second acceptable option. P3 §8 said the opposite. That paragraph described the Monte Carlo script, which still used a CLR-based bound, and added "as in the package", which had become false. Both are now corrected.

* **[Theory] Proposition P3.5(c)** (P3 §5).
  * **Setup:** for a cell, set $\pi_i=w_i^2/\sum w^2$ and $\bar m=\sum\pi_i m_i$, which is a linear target. By concavity, $\sum\pi_i m_i(1-m_i)\le\bar m(1-\bar m)\le\max_{m\in \mathrm{CI}_v}m(1-m)$.
  * **(c1):** suppose $\mathrm{CI}_v$ comes from P3.5(a) or (b2)(ii), on the rung's own constraint system, at level $1-\alpha_v$, and $\max_i w_i^2/\bar s^2\to0$. Then the interval with this allowance has coverage at least $1-\alpha_{est}-\alpha_E-\alpha_v$. The proof is a union bound, with the residual step taken relative to the deterministic $\bar s$.
  * **The extra condition:** the condition on $\bar s$ is new; we found it while writing the proof. If $\bar m$ tends to 0 or 1, the Lindeberg step relative to $\bar s$ can fail. The Popoviciu version needs only $\max w^2/\sum w^2\to0$.
  * **(c2):** with a CLR interval, coverage holds conditional on that interval covering $\bar m$. The status is algorithmic, as for (b3). Your numerical illustration ($.09<.25$) is in the text.
* **[Code] Package 0.3.1** (`R/rxc_infer.R`, `R/rxc_api.R`):
  * `proj_actual` takes $\bar m$ from a projection interval computed on the **rung alone**, with no optional restriction rows, as (c1) requires.
  * `ci_actual` (CLR) keeps its CLR-based bound, and `print()` labels it "CLR (algorithmic); its Bhatia–Davis bound also uses a CLR interval".
  * When the $\bar m$ interval lies within $[0,.05]$ or $[.95,1]$, a note says the (c1) condition is doubtful and recommends `fpc = "popoviciu"`.
  * The legend for `proj_actual` cites P3.5(b2)(ii) and the union bound.
  * Tests are in `tests/testthat/test-round2.R`. No reference output changed.
* **[Theory/Check] Monte Carlo rerun** (P3 §8; `phase3_montecarlo.py`, 1000 replications, same seeds).
  * The projection column for $B$ now uses (c1). The old CLR-based column is kept, labelled "not covered".
  * Coverage is unchanged in every design, and median widths move by at most 0.002 (design E: .251 vs .249). So the provably valid version costs essentially nothing here.
* **[Check]** `verify_phase3.py`, "P3.5(c) [V]": (i) your illustration; (ii) a simulation with BD taken from a covering interval, coverage 0.958 at nominal 0.95.
* **[Theory]** In the same round the Batch B referee pointed out a related issue. Reduced-rank case (ii-b) of P3.5(b2) needs an exact reduction; a numerical rank cutoff does not provide one. P3 §5 now says so. The package's critical value uses $q=\dim\beta$, and dropping candidates in numerically null directions is labelled a safeguard.

### Comment 4: Proposition D8

> *Replace the claimed equivalence with a distinction between identification by algebraic cancellation and identification given the feasible set.*

Agreed. Your counterexample is right: with $Q_{1,12}=0$ and the restriction imposed over the whole simplex, probability validity forces $K_{1,12}=0$, and the change is identified even when $\Delta q\ne0$.

* **[Theory]** 2b Proposition D8 now has two parts.
  * (a) **By cancellation:** the change is invariant to every $K\in\mathcal K$ if and only if $\Delta q_{rs}=0$ for all $s\ne r$.
  * (b) **Given the identified set:** the change is point-identified if and only if $\max_{K\in\Theta}L(K)-\min_{K\in\Theta}L(K)=0$, with $L(K)=\sum_{s\ne r}K_{c,rs}\Delta q_{rs}$.
  * Your example follows the proposition. A correction note records the earlier "iff". The proof steps and the cross-reference in 2b §6 are updated.
* **[Theory]** The reading now speaks of a "change in ecological transition rates across periods" and says a causal DiD interpretation needs further structure (developed separately in PHASE2d, which is not yet reviewed).
* **[Check]** `verify_phase2b.py`, "D8b [F]": your example (identified with $\Delta q\ne0$, failing the old "iff") and a case with positive width.

### Comment 5: Independence as a substantive assumption

> *Present independence as a substantive stochastic assumption that enables the proposed inference, not as something implied by observing the entire electorate.*

Agreed. Our earlier text let the second reading through.

* **[Theory]** P3 §0.4 has a new bullet, "Independence is a substantive assumption". It separates two statements:
  1. Conditioning on observed district effects and *assuming* independent residuals (P3-S1a) is coherent.
  2. "Every unit is observed, so dependence need not be considered" does not follow.

  It gives your example: within-district shocks with $\operatorname{Cov}(e_i,e_j\mid\text{design})\ne0$ and no mean shift violate P3-S1a but not the mean model. It also says that P3.5 and P3.8 are not covered under such dependence, and that P3-S1a has to be argued for each application, with residual correlation as a diagnostic. No clustered inference is reintroduced.
* **[Theory]** A sentence after P3.6(c) says the pooling test examines district mean shifts only and is not a test of conditional independence.
* **[Docs]** The vignette's inference paragraph says the same.

### Comment 6: Theorem P3.8

We agree with your characterisation and have kept the emphasis where you suggest. P3.8 stays as the formal fallback guarantee, described as valid but uninformative in the designs simulated. Informative finite-sample inference for many small groups remains an open item (P3 §10), not a claimed solution.

### Comment 7: The CLR implementation

Thank you. CLR remains the package default, labelled algorithmic, with the projection interval printed beside it. From 0.3.1, `print()` also states which result covers each fit's intervals. The open question in P3 §10 now reads as you put it: whether the simulated candidate collection supports a critical value that is valid **uniformly** over the relevant configurations. We make no claim of a generally valid adaptive CLR interval.

## 4. Status of the results you assessed

| Result | Status after this round |
|---|---|
| R×C identification (2a), district extension (2b) | Established; unchanged this round |
| Theorem P3.4 | Established for a fixed constraint system |
| Proposition P3.4′ | New; proved under the localisation condition (register R17) |
| P3.5(a), P3.5(b2) | Established; (b2)(ii-b) only for an exact reduction |
| P3.5(c1) / (c2) | New. (c1) proved; (c2) conditional, algorithmic (register R18) |
| Conjecture P3.5(b3) (CLR default) | Open; now asks for a uniform critical value |
| Theorem P3.8 | Established; formal fallback, not informative in the designs simulated |
| Proposition D8 | Corrected: (a) cancellation, (b) the condition given the identified set |
| P3.6(c) pooling test | A test of mean shifts only; it does not test independence |

## 5. Still open

* Uniform limit theory for the plug-in endpoint under changing active sets, beyond P3.4′'s localisation condition (P3 §10).
* Conjecture P3.5(b3).
* Informative finite-sample inference for many small groups.
* An automatic Popoviciu fallback in the package when the (c1) condition is doubtful; at present a note is printed.
* As before: covariate adjustment and the turnout application.

We are now writing a self-contained technical report collecting the identification and inference results, with every correction from both rounds, along the lines you encourage.
