# Second-Round Referee Report (Phases 1, 2a, 2b, 3)

*Received 18 Nov 2026 from the original referee of Phases 1–3. Reviewed version: commit `e170be9`. Transcribed as received; some display equations were partly lost in transfer and are marked [garbled].*

**Materials:** Revised Phases 1, 2a, 2b and 3
**Recommendation:** Substantial progress; targeted revisions recommended.

## 1. Overall assessment

The revision satisfactorily addresses the large majority of the concerns raised in my previous reports.

I particularly welcome the following changes:

* The corrected probability-validity and compactness assumptions.
* The expanded proof of the conditional-selection characterization.
* The distinction between sharp and outer identified sets.
* The corrected conditioning argument in Proposition D4(b).
* The explicit model-based conditional-independence assumption for inference.
* The separation of the three estimands.
* The distinction between established inference results and conjectural CLR implementations.
* The new finite-sample confidence interval in Theorem P3.8.

The additional distinction between the support of the restricted measure and the intersection of the original support with $A$ is also mathematically appropriate.

I no longer regard the issues raised in my original review as substantial obstacles to the main identification results.

Nevertheless, several claims in the revised documents deserve further attention.

## 2. Major comment: Theorem P3.4 and changing sample hulls

**Location:** Phase 3, Theorem P3.4.

The theorem applies the directional delta method to an LP endpoint represented as a minimum of finitely many affine functions. The proof is correct for a fixed collection of affine functions.

However, the theorem is now described as applying to the sample-hull endpoint, conditional on the observed design. This creates a problem when the sample size increases. As $n$ increases, the sample hull can change. Consequently:

1. The LP constraint matrix can change.
2. The collection of dual vertices can change.
3. The active constraints can change.
4. The differences between nearly binding constraints can shrink at rates comparable to $n^{-1/2}$.

Ordinary directional delta-method theory for a fixed function does not automatically establish the asserted limit for this changing sequence of functions.

**Counterexample.** Consider the sequence of piecewise-linear functions
$$f_n(\theta)=\min\{\theta,\,-\theta+n^{-1/2}\}.$$
Let the true parameter be $\theta_0=0$ and suppose $\sqrt n(\hat\theta-\theta_0)\Rightarrow Z$, $Z\sim N(0,1)$. For every finite $n$, the first affine function is uniquely active at $\theta_0$. Applying the theorem's unique-active-gradient argument would therefore suggest a standard normal limit. But
$$\sqrt n\,\big(f_n(\hat\theta)-f_n(\theta_0)\big)\Rightarrow \min\{Z,\,1-Z\}$$
[display partly garbled in transfer; this is the evident intended statement], which is not normally distributed.

The problem is that the gap between the two candidate functions shrinks at exactly the estimation-error rate. This example has the same finite-minimum-of-affine-functions structure as the LP endpoint. It demonstrates that fixing the design separately for each $n$ is not, by itself, sufficient.

**Recommendation.** Clarify whether Theorem P3.4 is intended as:

* **A. A fixed-constraint-system result.** In this case, the current proof is essentially satisfactory, but the theorem should not be presented as automatically covering a sequence of expanding sample hulls.
* **B. A triangular-array result for sample-hull endpoints.** In this case, additional assumptions are required concerning the stability of the relevant active constraints or the convergence of the corresponding directional derivatives. A possible sufficient condition would be a suitable separation requirement ensuring that asymptotically inactive candidates cannot become relevant under root-$n$ perturbations. More generally, the authors could establish convergence of the sequence of localized endpoint maps.

I would favor explicitly proving the fixed-constraint result first and treating inference under changing active sets as a separate extension.

This issue concerns the claimed asymptotic distribution, not the validity of the conservative confidence-set construction or the new finite-sample result.

## 3. Major comment: The implemented projection interval for $B$ is not necessarily covered by the projection theorem

**Location:** Phase 3, Proposition P3.5(b2) and §8.

The revised projection argument is convincing under the stated positive-definiteness or exact reduced-rank conditions. In particular, the use of a simultaneous ellipsoidal confidence event is a sound way to obtain coverage for arbitrary data-dependent collections of dual-feasible candidates.

However, §8 acknowledges that the implemented projection interval for the realized target $B$ uses a Bhatia–Davis residual allowance derived from CLR intervals. The problem is that the implemented CLR candidate-selection procedure remains conjectural. Consequently, if the CLR-derived bounds on the conditional mean fail to contain the truth, the resulting Bhatia–Davis variance bound need not be conservative. The projection argument does not correct this problem.

For example, suppose the true conditional mean is $m=0.5$, so the Bhatia–Davis variance bound is $m(1-m)=0.25$. If an incorrectly narrowed interval places $m$ between 0 and 0.1, the estimated variance upper bound is at most $0.1(0.9)=0.09$. The resulting residual allowance can be considerably too small.

The authors acknowledge this qualification in §8, but the practical recommendations still describe the projection interval as possessing an asymptotic guarantee.

**Recommendation.** Distinguish two implementations:

* **Provably valid projection interval for $B$:** Use the Popoviciu residual allowance from Theorem P3.5(a), or another variance upper bound derived from a confidence set with formally established coverage.
* **Variance-adaptive projection interval:** Use the Bhatia–Davis allowance computed from estimated conditional-mean bounds, but explicitly characterize its coverage as conditional on the validity of those bounds.

The package should distinguish the two versions accordingly. The projection interval for $\widetilde B$ does not have this particular problem because it does not require the finite-population residual allowance.

## 4. Proposition D8: The condition is sufficient, but not necessary, for point identification

**Location:** Phase 2b, Proposition D8.

The proposition correctly derives [display garbled in transfer]
$$\Delta(\text{transition rate}) = \text{identified component} + \sum_{s\ne r}K_{c,rs}\,\Delta q_{rs},$$
where $q_{rs}$ is the moment ratio $\mu_{rs}/\mu_r$ in each period [garbled; e.g. $q^{[\mathrm{pre}]}_{rs}=\mu^{[\mathrm{pre}]}_{rs}/\mu^{[\mathrm{pre}]}_r$].

The proposed moment-equality condition, $\Delta q_{rs}=0$ for all $s\ne r$, is sufficient to eliminate dependence on the unidentified contextual effects. However, the proposition describes this as an if and only if condition for point identification. That is too strong without an additional qualification. Point identification can also occur when the identified set for $K$ fixes the relevant linear combination, even if some $\Delta q_{rs}\ne0$.

**Counterexample.** Consider the $2\times2$ model with $Q_{1,12}=0$ and the contextual-effects restriction imposed over the whole simplex. Probability validity then requires $K_{1,12}=0$. The relevant contextual-effect parameter is point-identified. Now allow the distributions of $X$ to differ between the pre- and post-periods so that
$$\frac{\mu^{[\mathrm{post}]}_{12}}{\mu^{[\mathrm{post}]}_1}\ne\frac{\mu^{[\mathrm{pre}]}_{12}}{\mu^{[\mathrm{pre}]}_1}.$$
The moment-equality condition fails, but the difference in transition rates is nevertheless point-identified because $K_{1,12}=0$. The proof's final step effectively acknowledges this possibility.

**Recommendation.** Replace the claimed equivalence with a distinction between:

1. **Identification by algebraic cancellation:** The moment-equality condition makes the change invariant to every possible value of $K$.
2. **Identification given the feasible set:** The change is point-identified whenever the linear functional of $K$ is constant over the identified set.

Formally, the latter condition is $\max_{K\in\Theta}L(K)-\min_{K\in\Theta}L(K)=0$ [garbled in transfer; this is the evident intended statement], where $L(K)=\sum_{s\ne r}K_{c,rs}\Delta q_{rs}$. This gives an exact necessary and sufficient condition for point identification, conditional on the relevant identified set.

The moment-equality restriction remains an attractive, empirically interpretable sufficient condition.

I would also avoid describing this as a causal difference-in-differences identification result without additional assumptions. The established result concerns differences in ecological transition rates across periods; a causal interpretation would require additional structure.

## 5. The independence assumption is now explicit, but its substantive justification should be qualified

**Location:** Phase 3, §0.4.

The revision correctly recognizes that the previous combination of cluster-robust regression inference and independent unit-level residuals was inconsistent. The new framework resolves the mathematical inconsistency by imposing conditional independence. That is an acceptable modeling choice.

However, I would distinguish more carefully between two claims:

1. The model conditions on observed district effects and assumes that the remaining unit-level residuals are independent.
2. Because every electoral unit is observed, residual dependence need not be considered.

The first is a coherent assumption. The second does not follow automatically. For example, unobserved within-district shocks could induce $\operatorname{Cov}(e_i,e_j\mid\text{design})\ne0$ without necessarily producing nonzero district-specific mean shifts. Such dependence would violate P3-S1a without necessarily violating the pooled conditional-mean restriction. Consequently, Proposition P3.6(c)'s pooling test cannot be interpreted as a general test for all forms of district-level dependence. The test examines a specific class of mean restrictions. It does not validate conditional independence.

I recommend making this distinction explicit in §0.4. I am not recommending that clustered inference be reintroduced into the current theory. Rather, the authors should present independence as a substantive stochastic assumption that enables the proposed inference, not as something implied by observing the entire electorate.

## 6. Theorem P3.8: A successful addition

**Location:** Phase 3, §5.1.

I regard the new finite-sample confidence interval as a meaningful addition. The proof appropriately exploits the boundedness of the observed outcome shares and the conditional independence assumption. The use of Hoeffding's inequality to construct simultaneous bounds for linear functionals of the regression coefficients is legitimate.

The resulting interval does not require a Gaussian approximation, Slater's condition, a fixed regression dimension, or a fixed number of districts. The simultaneous confidence event contains the true contextual-effect parameter with the advertised probability, and the separate Hoeffding allowance controls the realized residual contribution. Thus, the logic of the finite-sample coverage theorem is convincing under its stated assumptions.

The limitation is practical rather than mathematical. As the authors acknowledge, the interval covers the entire logical range in the simulations conducted so far. I would therefore retain Theorem P3.8 as the formal fallback guarantee for difficult asymptotic regimes, but would not emphasize it as a practically informative solution for many-small-district inference. This is an appropriate distinction between mathematical validity and statistical informativeness.

## 7. The CLR implementation: The revised treatment is appropriate

I accept the separation of Proposition P3.5 into the established theoretical results and Conjecture P3.5(b3). In particular:

* The weak-duality argument for arbitrary candidate sets is sound.
* The projection-based simultaneous guarantee is a useful contribution.
* The candidate-generation algorithm remains a conjectural approximation to the complete CLR procedure.
* The Monte Carlo evidence provides useful evidence about its practical behavior.

The 1,000-replication simulation study is a substantial improvement. I would retain the experimental CLR procedure in the package provided that its inferential status is clearly communicated.

The remaining research question is whether the simulation-based candidate collection can support a valid critical value uniformly over relevant parameter configurations. This does not need to be resolved before publishing the core identification results, but it would need to be resolved before claiming a generally valid adaptive CLR confidence interval.

## 8. Final assessment

The revised theory is much stronger than the original submission. I regard the principal $R\times C$ identification results and the district extension as substantially established. The author has also satisfactorily addressed most of the concerns about the statistical inference framework.

The most important outstanding questions concern the transition from fixed LP constraints to changing sample-hull constraints, and the relationship between formally established projection intervals and variance-adaptive implementations.

I would recommend the following priorities:

1. Clarify or restrict Theorem P3.4's asymptotic claim for changing sample hulls.
2. Make the formally valid projection interval for $B$ operationally distinct from its CLR-dependent variance-adaptive version.
3. Correct Proposition D8's necessity claim.
4. Clarify that conditional independence is a substantive assumption not established by the pooling test.

I do not believe these points require a fundamental reworking of the identification framework.

My revised recommendation is favorable, conditional on these targeted corrections. I regard the principal contribution as mathematically credible, and I would encourage the authors to proceed toward a methodological manuscript while keeping the unproved inference extensions explicitly separate from the established results.
