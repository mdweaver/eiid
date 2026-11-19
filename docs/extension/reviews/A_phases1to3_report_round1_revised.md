# Referee report, round 1 (revised): Phases 1, 2a, 2b and 3

*Received Nov 2026 from the original referee of Phases 1–3 (round 1). Transcribed as received (pasted text); formatting from the original may be lost. The response plan is `../REVIEW_RESPONSE_PLAN.md`; the response is `../REVIEW_RESPONSE.md`.*

Revised Reviewer Report
Manuscript: Extension of linear contextual-effects partial identification bounds in ecological inference to the \(R\times C\) setting, with district heterogeneity and statistical inference.

Materials reviewed:
- PHASE1_foundations.md
- PHASE2a_rxc.md
- PHASE2b_districts.md
- PHASE3_inference.md
Recommendation: Major revision, with a positive assessment of the overall theoretical contribution.
1. Revised overall assessment
Having now read Phases 2b and 3, I am more positive about the scope and potential contribution of the project than I was after reviewing Phase 2a alone.

The additional work does three important things.

First, Phase 2b gives the \(R\times C\) extension a broader application. The district framework is not simply an application of the pooled model with district-specific parameters. It develops a hierarchy of restrictions that connects fully district-specific contextual effects, common contextual effects with district-specific baselines, bounded contextual-effect heterogeneity, and the Goodman constancy assumption.

Second, Phase 2b establishes an interesting result about the relationship between district-level information and aggregate identification: district-level restrictions can sharpen aggregate bounds even when no additional unidentified parameters are eliminated.

Third, Phase 3 addresses a particularly important weakness of ecological-inference applications: the difference between identification with known population quantities and inference when regression coefficients, constraints, and the composition hull must all be estimated.

I regard these as genuine methodological contributions.

However, reading the three phases together also reveals some problems that were not apparent from Phase 2a alone. In particular, I have concerns about the asymptotic regimes covered by Phase 3, the validity of the proposed confidence intervals under the full set of sampling assumptions, and the claimed validity of the computational implementation of the CLR procedure.

My revised assessment is:
Component	Assessment
Phase 1: identification framework	Strong, with a few remaining assumption issues
Phase 2a: \(R\times C\) identification	Strong; principal mathematical issues remain as previously identified
Phase 2b: district extension	Very promising; main algebra and nesting results appear correct
Phase 3: estimation and inference	Important contribution, but several inferential claims need qualification or additional proof
Overall	Worth developing into a publishable methodological contribution

My most significant new concern is that Phase 3 currently makes broader claims about the validity of inference than its stated sampling assumptions justify. I would address this before treating the project as theoretically complete.

2. Revisions to my previous Phase 1 / Phase 2a comments
Before discussing the new material, I want to update three points from my previous report.
2.1. The population support versus sample hull issue is substantially better addressed
My previous review identified the relationship between the population support and the observed sample hull as an important unresolved issue.

Phase 2b now explicitly defines
\[
A_d=\operatorname{conv}(\operatorname{supp}P_{X\mid D=d}),
\]

with the conditional support region defined separately for each district.

This is an appropriate extension of the pooled model. The condition
\[
P(X\in A_D)=1
\]

also ensures that the district-specific affine representation can be used in deriving the aggregate estimand.

Phase 3 then explicitly distinguishes the population hull \(A_d\) from its sample counterpart \(\hat A_d\), establishing that
\[
\hat A_d\subseteq A_d.
\]

It shows that using the sample hull produces a weakly larger outer identified set and that the corresponding bounds converge when the sample hull converges to the population hull in Hausdorff distance.

This is exactly the distinction I wanted to see. It substantially resolves my earlier objection at the level of identification theory. 
GitHub+1


There is still an inferential qualification, which I discuss in §4.3: convergence of the hull is not the same as root-\(n\) convergence, and Phase 3 does not establish that the sample-hull effect is negligible at the root-\(n\) scale.
2.2. The S3 qualification identified in my previous review remains unresolved
Phase 2a still states that S3 is implied when
\[
A\subseteq\operatorname{conv}(\operatorname{supp}P_X).
\]

As I noted previously, the proof establishes the result when
\[
A=\operatorname{conv}(\operatorname{supp}P_X),
\]

not for an arbitrary subset of that convex hull.

The distinction remains important.

Phase 2b's default specification,
\[
A_d=\operatorname{conv}(\operatorname{supp}P_{X\mid D=d}),
\]

avoids the problem. But the more general Phase 2a statement should still be corrected.

I recommend replacing the claim with:
When the contextual-effects restriction is imposed on the entire convex hull of the population support, probability validity follows from the accounting identity and the conditional-mean restriction. For a more general region \(A\), probability validity must be imposed separately unless \(A\) is contained in the convex hull of the support on which the conditional-mean restriction is imposed.

This is a relatively easy correction. 
GitHub

2.3. The LP value-function regularity claim also remains incorrect
Phase 2a Proposition 7 still states that the LP value function is concave and piecewise linear jointly in \((Q,\mu)\).

This claim needs to be corrected.

Phase 3 actually makes the important distinction that was missing in Phase 2a: conditional on the design, the objective coefficients and constraint matrix are fixed, and the endpoint is a minimum of finitely many affine functions of the regression coefficients.

That is a valid and useful result.

However, it does not retroactively establish joint concavity or piecewise linearity in \((Q,\mu)\).

I recommend retaining Phase 3 Lemma P3.1 but changing the Phase 2a Proposition 7 concluding paragraph to distinguish:
- Variation in the constraint right-hand side;
- Variation in the objective coefficients;
- Joint variation in both.
The joint claim is generally false. 
GitHub+1


3. Review of Phase 2b: District heterogeneity and the assumption ladder
3.1. Overall assessment
I consider Phase 2b a strong extension.

The central conceptual contribution is the distinction between three situations:
- M2: completely district-specific contextual effects;
- M1: common contextual effects but district-specific baselines;
- M1(\(\delta\)): bounded deviations from common contextual effects.
This distinction is important because it separates the identifying contribution of cross-district restrictions from the identifying contribution of the individual district distributions.

The resulting assumption ladder is also a useful way to organize sensitivity analysis.

I found no fundamental algebraic error in the main derivations. There are, however, several qualifications that should be incorporated.
3.2. Proposition D0: Districts add nothing without cross-district restrictions
Phase 2b, §1, Proposition D0.

This result is correct and important.

The proposed bounds are
\[
DD^*_{rc}
=
\left[
\frac{E[N\max\{0,X_r+T_c-1\}]}{\mu_r},
\frac{E[N\min\{X_r,T_c\}]}{\mu_r}
\right].
\]

The claim that these sharp bounds are unaffected by conditioning on districts is correct.

The proof through the transportation-polytope representation is convincing.

In particular, the construction fixing cell \((r,c)\) at either its Fréchet upper or lower limit and then solving the residual transportation problem is appropriate. The application of Gale's supply-demand theorem correctly establishes feasibility of the residual table.

The measurable-selection step then allows the pointwise extremes to be aggregated and attained.

I also agree with the extension to the joint set for the vector \(B^*\). The district variable does not enter the unit-level feasibility condition, so it cannot independently tighten the bounds without an additional restriction linking units across districts.

Recommendation: Retain this result essentially as written.

One refinement: distinguish the sharp marginal interval for each component from the joint identified set for the entire transition matrix. The joint set is not generally the Cartesian product of the individual Duncan–Davis intervals. The document correctly notes the joint selection characterization, but this distinction deserves emphasis. 
GitHub

3.3. Lemma D1: Identification of district baselines
Phase 2b, §2, Lemma D1.

I find this result correct.

Under M1,
\[
H_c^{(d)}
=
H_c+\zeta_{c,d}\mathbf 1^\top.
\]

Taking symmetric parts gives
\[
Q_c^{(d)}
=
Q_c+\frac12
\left(
\zeta_{c,d}\mathbf 1^\top+
\mathbf 1\zeta_{c,d}^\top
\right).
\]

Consequently,
\[
\boxed{
\zeta_{c,d}
=
\operatorname{diag}
\left(Q_c^{(d)}-Q_c^{(1)}\right).
}
\]

The baseline shifts are identified.

The off-diagonal restrictions,
\[
(Q_c^{(d)}-Q_c^{(1)})_{rs}
=
\frac12(\zeta_{c,d,r}+\zeta_{c,d,s}),
\]

are also correct.

The important insight here is that the district-specific symmetric components do not introduce additional unidentified parameters under M1. The unidentified dimension remains
\[
d_0=(C-1)\frac{R(R-1)}2.
\]

This is a useful result because it means the district extension does not compound the dimensionality of the unidentified parameter space.
One important qualification: distinguish identification from testability
Lemma D1 establishes population identification of the district shifts and overidentifying restrictions.

It does not establish that these quantities can be estimated precisely in small districts.

The distinction is particularly consequential for the intended riding-level application.

Remark D1.1 correctly notes that, under M1, the district-specific shift \(\zeta_{c,d}\) requires less within-district variation than the saturated M2 model requires to identify \(Q_c^{(d)}\).

I would elevate this point in the exposition.

Under M1, the district-specific parameters enter linearly:
\[
E[T_c\mid X,D=d]
=
X^\top Q_cX+\zeta_{c,d}^\top X.
\]

Once \(Q_c\) is identified from the pooled variation, identifying \(\zeta_{c,d}\) requires only enough within-district composition variation to identify a linear function on the simplex.

This is a substantial practical advantage.

However, the precise rank condition for the joint M1 regression should be stated. It is not sufficient simply to say that each district has \(R\) distinct compositions; those compositions must be affinely independent, and the pooled quadratic and district-shift design matrix must have full column rank.

I would replace the informal sufficiency claim with an explicit regression-design rank condition. 
GitHub+1

3.4. Lemma D2: The aggregate estimand under district heterogeneity
Phase 2b, §2, Lemma D2.

This derivation is correct.

The aggregate estimand is
\[
B^*_{rc}
=
\frac1{\mu_r}
\sum_d
\left[
\sum_s Q^{(d)}_{c,rs}\mu^{(d)}_{rs}
+
\sum_{s\ne r}K^{(d)}_{c,rs}\mu^{(d)}_{rs}
\right].
\]

Under M1, the unidentified part of each district's contextual-effects matrix differs from the common \(K\) only by the identified baseline shift:
\[
K_c^{(d)}
=
K_c+\operatorname{skew}
(\zeta_{c,d}\mathbf 1^\top).
\]

The coefficient on the common unidentified parameter therefore aggregates to
\[
\frac{\mu_{rs}}{\mu_r},
\]

exactly as in the pooled case.

This is an elegant result: district-specific baselines change the identified component of the estimand but do not change its sensitivity to the remaining common contextual effects.

I would emphasize this more strongly. It demonstrates that M1 preserves the interpretability of the original \(R\times C\) unidentified parameters, even while allowing meaningful between-district heterogeneity.

No substantive correction is needed here. 
GitHub

3.5. Theorem 1′ and Theorem 2: Sharp identified sets for the different rungs
Phase 2b, §3.

I find the generalization of the conditional-selection theorem appropriate.

Theorem 1′ is a useful abstraction: the proof of Phase 2a Theorem 1 does not depend specifically on the mean being \(M_H(X)\). It requires an appropriate measurable, bounded, finite-dimensional affine conditional-mean specification.

That allows the district restrictions to be imposed on the parameter space separately from the underlying conditional-selection argument.

The resulting characterization,
\[
\Theta_I^\rho
=
\left\{
\theta\in\mathcal R_\rho:
\text{each district satisfies its conditional support inequalities}
\right\},
\]

is correct under the stated assumptions.

I particularly like the distinction between the sharp set and outer set for each rung:
\[
\Theta_I^\rho\subseteq\Theta_O^\rho.
\]

The outer set is polyhedral, while the sharp set generally requires conditional distributional information.

This carries forward the strongest feature of Phase 2a.
A qualification concerning Theorem 1′
The proof says that replacing \(M_H(X)\) with \(m_\theta(V)\) carries over directly.

I agree with the conclusion, but the authors should be more explicit about the measurable-selection conditions needed for the generalized theorem.

In particular, state that:
- The conditional-mean correspondence is measurable in \(V\);
- The parameterized mean is uniformly bounded on the relevant domain;
- The parameter-to-mean map is affine and continuous;
- The underlying selection set remains nonempty and compact.
The current proof references the earlier functional-analytic argument but does not restate enough of the assumptions for this to be fully self-contained.

This is not a fundamental objection. 
GitHub

3.6. The assumption ladder is mathematically sound
Phase 2b, Theorem 2(c).

The nesting
\[
\begin{aligned}
\mathcal B^{L6}
&\subseteq\mathcal B^{L5}
\subseteq\mathcal B^{L4(r)}
\subseteq\mathcal B^{M1}\\
&\subseteq\mathcal B^{M1(0)}
\subseteq\mathcal B^{M1(\delta)}
\subseteq\mathcal B^{M2}
\subseteq DD^*
\end{aligned}
\]

is correct, subject to the stated nonemptiness qualifications.

The proof appropriately separates:
- Nesting of model classes;
- Nesting of identified parameter sets;
- Nesting of the resulting estimand sets.
I particularly appreciate the clarification that M1(0) is not quite M1: it preserves the common-skew-parameter relationship but does not impose the testable symmetric-part restrictions.

This is conceptually useful because it allows researchers to distinguish a failure of the testable implications of M1 from uncertainty about the unidentified contextual effects.

One recommendation: Make this distinction prominent in the presentation of the ladder. M1(0) is a sensitivity benchmark rather than simply another substantive model in the same sense as M1.

GitHub+1
3.7. Proposition D3: Closed form under M1 with two rows
This result is correct.

The district-specific inequalities generate
\[
\ell_c=
\max_d\max_{t\in\{l_d,u_d\}}
\left[
\lambda_c(t)-\frac{\zeta_{c,d,1}}{1-t}
\right],
\]

and
\[
h_c=
\min_d\min_{t\in\{l_d,u_d\}}
\left[
\rho_c(t)+\frac{\zeta_{c,d,2}}{t}
\right].
\]

The resulting feasible set is
\[
\Theta_O^{M1}
=
\left\{
\kappa:
\ell_c\leq\kappa_c\leq h_c,\quad
\sum_c\kappa_c=0
\right\}.
\]

I checked the sign logic in the two row-specific inequalities. The district shifts enter with the correct signs.

The interpretation is also compelling: a district with an extreme baseline and a composition range near a simplex boundary can strongly constrain the shared contextual effect.

This provides a particularly accessible substantive explanation of why district fixed effects can sharpen ecological-inference bounds.

The condition \(0<l_d\leq u_d<1\) should be emphasized as necessary for the displayed closed form, since the expression divides by \(t\) and \(1-t\). For districts containing pure-composition observations, the limiting constraints should be derived separately rather than applying the formula mechanically. 
GitHub

3.8. Propositions D4–D7: Additional district results
These propositions are, in general, convincing.

D4 (when dummies matter) is particularly interesting.

The argument that M0 holds under M1 if and only if
\[
x\longmapsto E[Z_D\mid X=x]
\]

is affine is correct, subject to the simplex representation and the coherence restriction that the district shifts sum to zero across columns.

It also provides a useful interpretation: a pooled regression can fail because district baselines vary nonlinearly with composition.

I would retain this result, while slightly qualifying the phrase that dummies are needed for validity "exactly when" this conditional expectation is non-affine. This is correct for the conditional-mean representation, but the resulting outer bounds also require probability validity over the specified pooled support region.

D5 (M2 adds district intervals) is correct. Under M2, the district-specific feasible sets are a Cartesian product, and the linear target is separable. The Minkowski-sum result follows.

D6 (bounded deviation interpolates) is also correct, conditional on the boundedness assumptions established for the district-specific feasible sets. Its piecewise-linearity argument follows from ordinary parametric LP theory.

D7 (the most restrictive rungs) is a useful result.

The L5 point-identification result and the count
\[
(C-1)(R-1)
\]

of unidentified \(K\)-coordinates pinned down by L4(r) are correct.

One point of exposition deserves emphasis: L5 is not an innocuous alternative estimator. It imposes zero contextual effects. The document acknowledges this in Remark D7.1, and I strongly agree with showing the L5 point estimate alongside the less restrictive bounds as a diagnostic, provided it is clearly identified as assumption-dependent. 
GitHub+1

3.9. The DiD connection is a worthwhile contribution
Phase 2b, Remark 6.2.

I particularly like this connection.

Treating election periods as districts produces a direct interpretation of M1 as an unchanged-contextual-effects assumption with period-specific baselines.

When the composition and size moments are held fixed between periods, the common unidentified \(K\) component cancels in the difference:
\[
B^*_{\text{post}}-B^*_{\text{pre}}.
\]

This is a genuine point-identification result under the specified restrictions.

I would elevate this application beyond a remark. It provides a bridge between ecological inference and difference-in-differences methodology.

One caution: the statement requires equality of the relevant \(N X_r X_s\) and \(N X_r\) moments, not merely superficially similar distributions of vote shares. The existing notation \(\mu^{(\mathrm{pre})}=\mu^{(\mathrm{post})}\) captures the necessary condition; I would state it explicitly in the verbal interpretation.

GitHub

Correction to §3.8: Proposition D4(b) contains an important conditioning error
On closer examination, I need to revise my positive assessment of D4(b). Part (a) is correct, but part (b) is not correct as stated. This is one of the more important issues I found in Phase 2b.

The problem is that D4 uses
\[
E[Z_D\mid X],
\]

whereas the pooled M0 model inherited from Phase 2a conditions on both \(X\) and \(N\).

Under M1, the district-specific conditional mean is
\[
E[\beta\mid X,N,D]
=
M_H(X)+Z_D.
\]

Taking conditional expectations over districts, the pooled conditional mean is actually
\[
\boxed{
E[\beta\mid X,N]
=
M_H(X)+E[Z_D\mid X,N].
}
\]

For M0 to hold, this expression must be representable as \(M_{H'}(X)\). Therefore, the correct condition is that
\[
\boxed{
E[Z_D\mid X,N]
}
\]

must be an affine function of \(X\) that does not depend on \(N\).

The current proof establishes a result about \(E[\beta\mid X]\), which is not the same thing.

There is also a concrete problem with the stated sufficient condition:
In particular, \(D\perp X\) suffices.

This is insufficient. District membership could be independent of composition but correlated with precinct size. In that case, \(E[Z_D\mid X,N]\) could still depend on \(N\).

A sufficient independence condition would instead be
\[
D\perp (X,N).
\]

More generally, it is sufficient that the conditional expectation of the district shift given \((X,N)\) be affine in \(X\) and invariant to \(N\).

Recommended correction: Rewrite Proposition D4(b) using the conditional expectation given \((X,N)\). If the authors intended a pooled model conditioning only on \(X\), they should explicitly distinguish that model from the M0 specification in Phase 2a.

This matters because the current version gives an incomplete characterization of when district fixed effects can be omitted without misspecifying the pooled conditional-mean model. 
GitHub+1


4. Review of Phase 3: Estimation and inference
4.1. Overall assessment
Phase 3 is ambitious and addresses the right problems.

I particularly appreciate the decision to distinguish the two targets:
\[
B_{rc}
=
\frac{\sum_i N_iX_{ir}\beta_{irc}}
{\sum_iN_iX_{ir}},
\]

and
\[
\widetilde B_{rc}
=
\sum_iw_{ir}m_{rc}(V_i),
\qquad
w_{ir}=
\frac{N_iX_{ir}}{\sum_jN_jX_{jr}}.
\]

This distinction matters for inference.

The first is the realized finite-population transition rate. The second is the conditional-mean transition rate implied by the model at the observed compositions.

The decomposition
\[
B_{rc}=\widetilde B_{rc}+E_{rc}
\]

is useful, with
\[
E_{rc}=\sum_iw_{ir}e_{irc},
\qquad
e_{irc}=\beta_{irc}-m_{rc}(V_i).
\]

This makes clear that statistical inference must accommodate two sources of uncertainty:
1. Estimation of the coefficients defining the identified set;
2. Residual variation between the realized finite-population rates and the conditional means.
That is the right conceptual framework.

However, I identify four significant problems with the inferential claims:
- The sampling assumptions do not support all the stated asymptotic regimes.
- The confidence interval for the finite-population target is not valid under every dependence structure permitted by P3-S1.
- The CLR implementation is not yet formally justified.
- The asymptotic distribution is stated for a fixed design and hull, while the actual estimator uses an estimated, expanding hull.
I discuss each below.
4.2. Lemma P3.1: The LP endpoint is a minimum of finitely many linear functions
Phase 3, §2, Lemma P3.1.

This is a correct and useful result.

For a fixed design, the upper endpoint can be written as
\[
\varphi_U(\beta)
=
c(\beta)+\min_{\lambda\in\Lambda}
\lambda^\top s(\beta),
\]

where
\[
\Lambda=
\operatorname{vert}
\{\lambda\geq0:A^\top\lambda=-\gamma\}.
\]

The dual feasible region does not depend on the regression coefficients. Consequently, the upper endpoint is a minimum of finitely many affine functions of those coefficients.

Thus,
\[
\varphi_U(\beta)
\]

is concave and piecewise linear in \(\beta\).

The lower endpoint has the corresponding convex representation.

This is precisely the right way to frame the nonsmoothness of the general \(R\times C\) LP.

It also resolves an important issue from my previous review: the authors are now correctly conditioning on the design, thereby fixing the objective and constraint matrices.

I would retain this lemma.

One qualification: it is a result for a fixed constraint system and fixed objective coefficients. In particular, it is not a general statement about joint variation in the identified quadratic coefficients, population moments, and LP constraints.

That distinction should remain explicit. 
GitHub

4.3. Proposition P3.2: Consistency and the estimated composition hull
This is an important result, and I am generally convinced by its logic.

The authors establish that
\[
\hat\varphi_U-
\varphi_U(\beta_0;\hat A)
=
O_p(n^{-1/2}),
\]

and that
\[
\varphi_U(\beta_0;\hat A)
\longrightarrow
\varphi_U(\beta_0;A)
\]

when the sample hull converges to the population hull in Hausdorff distance.

The one-sided effect of the sample hull is also correct:
\[
\hat A_d\subseteq A_d
\]

implies that the sample-hull constraints are weaker than the population-hull constraints.

Consequently,
\[
\varphi_U(\beta_0;\hat A)
\geq
\varphi_U(\beta_0;A),
\]

and
\[
\varphi_L(\beta_0;\hat A)
\leq
\varphi_L(\beta_0;A).
\]

This is a valuable result. It provides a clear explanation of why estimation using the empirical hull tends to widen, rather than artificially narrow, the outer identified set.
Three qualifications
First: The convergence result needs a more explicit sampling argument.

The proof uses a decreasing sequence of feasible sets as the sample hull expands.

That is straightforward when the samples are nested, as in the usual construction where each additional observation enlarges the existing sample. However, Hausdorff convergence alone does not imply that an arbitrary sequence of sample hulls is nested.

I would explicitly state that the result is being established along a nested sampling sequence, or supply a separate argument for non-nested Hausdorff-convergent sets.

Second: The Slater condition is restrictive.

The proof requires a strictly feasible point satisfying
\[
s_j(\beta_0)+a_j^\top\widetilde k
\geq\varepsilon>0
\]

for every constraint over the population hull.

This is a genuine interiority assumption.

It can fail when one or more implied transition probabilities are exactly zero at relevant compositions. Such boundary cases are not particularly exotic in ecological inference.

The authors should acknowledge that the root-\(n\) Lipschitz result is established for an interior case. Boundary cases may require different regularity arguments, and the confidence-interval results should not be represented as depending on Slater unless that is actually necessary.

Third: The hull error is not shown to be root-\(n\).

The authors correctly recognize that the estimation component can be root-\(n\) even while the hull component converges more slowly.

This distinction is crucial.

The result establishes
\[
\hat\varphi_U
-
\varphi_U(\beta_0;A)
=
O_p(n^{-1/2})
+
\underbrace{
\left[
\varphi_U(\beta_0;\hat A)
-
\varphi_U(\beta_0;A)
\right]
}_{\text{hull approximation error}}.
\]

The second term is nonnegative but is not established to be \(O_p(n^{-1/2})\).

This means the population-hull endpoint may not have the same root-\(n\) distribution as the plug-in endpoint.

I recommend retaining the convergence result, but being very careful about which endpoint is the target in the subsequent asymptotic distribution theorem. 
GitHub+1

4.4. Theorem P3.4: Directional asymptotics and non-unique LP solutions
I find the central argument correct.

The directional derivative of the upper endpoint is
\[
\varphi'_U(\beta_0;h)
=
\min_{\lambda\in\Lambda^*}
g_\lambda^\top h,
\]

where
\[
g_\lambda=c+S^\top\lambda
\]

and \(\Lambda^*\) denotes the optimal dual vertices at the true parameter.

Applying the directional delta method gives the stated minimum-of-correlated-normals limit.

This is precisely the relevant asymptotic structure for the general \(R\times C\) case.

I agree with the substantive interpretation:
- Unique active dual solution: ordinary differentiability and a normal limit.
- Multiple distinct active gradients: a non-Gaussian directional limit.
- Standard bootstrap: generally inconsistent when the derivative is genuinely nonlinear.
This is a valuable theoretical contribution because the unique-optimal-basis condition becomes increasingly demanding as the number of contextual-effect parameters and district constraints increases.
One mathematical qualification
The statement that a non-unique \(\Lambda^*\) necessarily produces a negative-mean minimum of correlated normals needs a qualification.

What matters is whether at least two active gradients produce distinct random variables under the limiting covariance matrix.

If
\[
(g_\lambda-g_{\lambda'})^\top
\Sigma
(g_\lambda-g_{\lambda'})
=0,
\]

then those gradients are indistinguishable in the directions in which the estimator varies.

Consequently, it is possible to have multiple optimal dual vertices without a genuinely nonlinear directional derivative.

The strict negative-mean result should therefore be stated in terms of distinct projected gradients, rather than simply a non-singleton optimal dual set.

This is especially relevant because the regression coefficients have adding-up restrictions, so the joint covariance matrix can be singular when all outcome columns are estimated together.

A useful solution would be to parameterize only \(C-1\) outcome columns and recover the final column through the adding-up restriction, or explicitly formulate the result using the support of the limiting Gaussian vector.

The rest of the argument is sound. 
GitHub+1

4.5. Theorem P3.5(a): The conservative confidence interval
This is where my most important objection arises.

The basic construction is sound under conditional independence of the relevant unit-level residuals.

The authors construct a simultaneous confidence event for the regression coefficients and the constraints defining the feasible contextual-effect parameters.

On that event, the true parameter \(k_0\) belongs to the expanded feasible set
\[
\hat\Theta^+
=
\left\{
k:
s(\hat\beta)+\hat k\hat\sigma_s+Ak\geq0
\right\}.
\]

The resulting interval therefore covers the conditional-mean target \(\widetilde B\), subject to the stated asymptotic Gaussian approximation.

The additional residual allowance for \(B\) is also reasonable.

Under conditional independence,
\[
\operatorname{Var}(E_{rc}\mid\text{design})
=
\operatorname{Var}
\left(
\sum_iw_{ir}e_{irc}
\mid\text{design}
\right)
\leq
\frac14\sum_iw_{ir}^2.
\]

Thus
\[
S_1=
\frac12
\left(\sum_iw_{ir}^2\right)^{1/2}
\]

is an appropriate variance bound.

The Lyapunov condition is a suitable way to obtain the stated asymptotic normal approximation for this weighted residual sum.

However, the theorem is not valid under all the sampling assumptions currently permitted by P3-S1.
The problem: dependence within districts
P3-S1 explicitly allows arbitrary dependence within districts for the pooled M0 model, with cluster-robust variance estimation.

But the proof of Theorem P3.5(a), steps 4–5, assumes that the \(e_i\) are independent conditional on the design.

Those are different assumptions.

Cluster-robust standard errors for the OLS regression do not solve the problem with the finite-population residual term.

For example, suppose all precincts within district \(d\) share a common random residual shock. Then
\[
e_{irc}=U_d
\]

for every precinct in that district, with mean zero and variance \(1/4\).

The residual contribution from district \(d\) is
\[
\sum_{i:D_i=d}w_{ir}U_d.
\]

Its variance is
\[
\frac14
\left(
\sum_{i:D_i=d}w_{ir}
\right)^2.
\]

This does not generally shrink at the rate implied by
\[
\frac14\sum_{i:D_i=d}w_{ir}^2.
\]

For equal weights among \(n_d\) precincts, the difference is a factor of \(n_d\).

Consequently, the stated \(S_1\) allowance can substantially understate uncertainty when within-district residuals are correlated.
Recommended correction
There are three possible approaches.

Option A: Restrict Theorem P3.5(a) for \(B\) to independent unit-level residuals.

This is the simplest solution. State explicitly that the finite-population confidence interval requires conditional independence across units, not merely heteroskedasticity-robust OLS inference.

Option B: Develop a cluster-level residual allowance.

If residuals are arbitrary within districts but districts are independent, define
\[
W_{dr}=\sum_{i:D_i=d}w_{ir}.
\]

A valid variance upper bound is
\[
\operatorname{Var}(E_{rc}\mid\text{design})
\leq
\frac14\sum_d W_{dr}^2.
\]

This follows by bounding the weighted residual sum within each district using the fact that each individual residual has a conditional support interval of length at most one.

An appropriate CLT would then require conditions on the district-level contributions. In particular, a growing number of independent districts and a no-dominant-district condition would be needed.

This would be a more coherent approach for the clustered M0 specification.

Option C: Restrict the confidence interval to \(\widetilde B\).

The conditional-mean target does not require the additional residual allowance. The existing argument for that target is unaffected by this particular issue.

My preference would be Option A as the initial formal result, with Option B as a subsequent extension.

GitHub+1
4.6. Theorem P3.5(b): The CLR procedure
I believe the underlying use of Chernozhukov, Lee and Rosen is appropriate.

The upper endpoint is an intersection bound:
\[
\varphi_U(\beta)
=
\min_{\lambda\in\Lambda}
\theta_\lambda(\beta).
\]

That is exactly the kind of object for which intersection-bound inference is useful.

The candidate-specific standard errors and selection step are sensible ways to avoid the severe conservatism of taking a simultaneous maximum over every possible dual vertex.

However, the manuscript does not yet prove the validity of the implemented CLR procedure.

The reason is explicitly acknowledged in the document.

The complete candidate set \(\Lambda\) is obtained by collecting optimal dual vertices at simulated draws
\[
\beta^\star\sim N(\hat\beta,\hat\Sigma).
\]

The authors then state that this collection is consistent for \(\Lambda^*\) as the number of draws increases, but acknowledge that its formal justification remains open.

I agree that this is a plausible computational strategy. But it is not yet a theorem.

There are two issues to resolve.

First: Candidate recovery.

The authors need to establish that the simulated optimization procedure recovers every candidate that is relevant for the directional derivative.

It is not enough to show that the procedure frequently returns the optimal vertex in finite samples.

The result must address the possibility of:
- Optimal vertices that are rarely selected;
- Degenerate optimal faces;
- Singular covariance matrices;
- Candidates that are relevant only in particular directions of the Gaussian perturbation.
The relevant target is the collection of gradients that actually determine the directional derivative on the support of the limiting Gaussian distribution.

Second: Increasing numbers of candidates.

For fixed \(R,C,D\) and a fixed sample hull, the set of dual vertices is finite.

But as the sample size and number of district constraints grow, the number of candidate vertices can also grow.

The CLR theorem invoked in the document is stated for a finite index set. The authors need to establish the conditions under which their increasing-index implementation inherits the required coverage guarantee.
Recommended revision
Separate the theoretical result from the proposed computational algorithm.

State a theorem establishing CLR validity when the full candidate set is known and the required Gaussian approximation holds.

Then state the simulation-based candidate-generation procedure as an approximation, with a separate proposition or conjecture concerning its consistency.

Until that latter result is proved, I would not describe the implemented CLR interval as formally established.

The conservative version in Theorem P3.5(a) is a more defensible baseline for the paper's formal contribution.

GitHub+1
4.7. Proposition P3.6: Model tests
I agree with the basic tests.

The emptiness test is valid as an asymptotic test: under the model, the true contextual-effect parameter belongs to the expanded feasible set on the simultaneous confidence event.

The Wald/F test of M1 against M2 is also appropriate.

There is, however, an important cross-phase discrepancy in the number of overidentifying restrictions.

Phase 3 gives
\[
(D^\star-1)(C-1)\frac{R(R-1)}2.
\]

This is the correct number of independent restrictions, accounting for the cross-column adding-up restriction.

Phase 1, §1.6, instead states \(R(R-1)/2\) restrictions for each column and district difference. This overcounts the independent restrictions by ignoring their dependence across columns.

I recommend correcting Phase 1 to match Phase 3.

The degrees-of-freedom calculation can be derived directly by comparing the dimension of the district-specific symmetric components with the dimension of the common symmetric component plus district shifts:
\[
(C-1)\frac{R(R+1)}2
-
(C-1)R
=
(C-1)\frac{R(R-1)}2.
\]

This confirms the Phase 3 count.

A separate qualification concerns the Wald/F approximation: it is straightforward under fixed-dimensional asymptotics, but needs additional conditions when the number of restrictions grows.

GitHub+1
4.8. Proposition P3.7: Many districts with few units
This section correctly identifies the incidental-parameter problem.

The authors are right that, when the number of districts increases but the number of units per district remains bounded, estimates of district-specific baseline shifts do not become consistent.

I also agree with the deterministic claim that dropping district constraints weakly enlarges the identified set.

However, I do not think the second part of Proposition P3.7 is established as stated.

It says that if a subset of districts is selected using the design alone, every result in §§3–5 holds conditionally on that subset.

Design-based selection is indeed helpful because it avoids selecting constraints based on their estimated outcomes.

But design-based selection does not, by itself, establish the asymptotic approximation or root-\(n\) consistency assumptions used in the preceding theorems.

In regime R2, with \(n_d\) bounded:
- District-specific shift estimates remain noisy.
- The number of district constraints grows.
- The dimension of the candidate set can grow.
- The Gaussian approximation for the maximum over candidates requires additional justification.
Dropping some constraints does not solve the underlying estimation problem if the retained districts still have bounded sample sizes.
Recommended correction
Retain the first part of P3.7 as a valid identification result.

Recast the second part as a conditional validity statement:
If the inferential procedure used on the retained districts is valid under the resulting sampling regime, selecting those districts as a function of the design alone preserves conditional coverage.

Then explain that standard fixed-dimensional OLS and CLR theory do not automatically deliver that validity in R2.

The alternatives in §6 are sensible research directions:
- Restricting attention to districts with adequate sample information;
- Using M1(\(\delta\));
- Pooling districts;
- Targeting \(\widetilde B\) instead of \(B\).
However, each should be accompanied by the appropriate asymptotic or finite-sample justification.

The current statement that the CLR procedure "restores coverage" should be understood as a simulation finding, not yet a general theorem. 
GitHub+1

4.9. The asymptotic regimes R1 and R3 need to be distinguished
There is a broader issue underlying the preceding concerns.

The manuscript identifies three regimes:
- R1: fixed number of districts, increasing observations per district;
- R2: increasing number of districts, bounded observations per district;
- R3: both the number of districts and observations per district increase.
It says that the results in Phase 3 apply to R1 and R3.

I do not think this is established.

Under R1, with fixed \(R,C,D^\star\), the regression parameter vector is finite-dimensional. Standard OLS asymptotics are appropriate, subject to the independence and rank conditions already discussed.

Under R3, however, the number of district-specific coefficients increases with the number of districts.

For M1, this includes the district-specific shifts \(\zeta_{c,d}\). For M2, it includes the entire district-specific collection of symmetric regression coefficients.

Consequently:
\[
\dim(\widehat\beta)\longrightarrow\infty.
\]

The standard finite-dimensional OLS CLT and directional delta method do not automatically apply.

The manuscript needs additional results concerning:
- The growth rate of \(D^\star\) relative to \(n_d\);
- Uniform consistency of the estimated regression coefficients;
- High-dimensional Gaussian approximation;
- The behavior of the number of relevant LP dual vertices;
- Uniform control of the variance estimators.
I would therefore recommend the following division:
Model / regime	What I think is currently supportable
M1/M2, R1	Standard asymptotic results are plausible and largely established under the stated rank and independence assumptions.
M0, R1 with independent units	Standard heteroskedasticity-robust inference is plausible.
M0, R1 with arbitrary within-district dependence	Cluster-robust inference is not justified with a fixed number of clusters.
M0, R3 with growing independent districts	Cluster-robust inference is possible under additional cluster asymptotics.
M1/M2, R3	Requires additional high-dimensional theory.
M1, R2	Requires an explicit many-moments/increasing-index inference argument.

This is an important qualification because the manuscript currently treats R1 and R3 together in its introductory statement of the results.

I would make R1 the principal formal asymptotic result and present R2/R3 as extensions requiring their own conditions.

GitHub+1
4.10. The two-target discussion should be made more precise
I like the introduction of \(\widetilde B\), but there are actually three related objects in the three-phase framework:

1. The superpopulation target from Phase 1:
\[
B^*_{rc}
=
\frac{E[NX_r\beta_{rc}]}{E[NX_r]}.
\]

2. The realized finite-population target:
\[
B_{rc}
=
\frac{\sum_iN_iX_{ir}\beta_{irc}}
{\sum_iN_iX_{ir}}.
\]

3. The realized-design conditional-mean target:
\[
\widetilde B_{rc}
=
\sum_iw_{ir}E[\beta_{irc}\mid V_i].
\]

The third is the conditional expectation of the second, given the observed design, under the stipulated sampling model.

It is a sample-weighted conditional-mean analogue of \(B^*\), but it is not generally identical to \(B^*\).

If the goal is inference conditional on the observed precinct composition and size information, the present approach is appropriate.

If the goal is inference for the superpopulation quantity \(B^*\), however, uncertainty associated with sampling the \(X_i,N_i,D_i\) must also be considered.

I would therefore change the description of \(\widetilde B\) to:
The finite-design conditional-mean transition rate, which is the conditional expectation of the realized finite-population rate given the observed design. Under suitable sampling conditions it estimates the superpopulation target.

There is also a small mathematical wording issue at the end of §0.2. The document says that both \(B\) and \(\widetilde B\) are linear in \(K\).

Strictly, only \(\widetilde B\) is affine in \(K\). The realized \(B\) also contains the residual component \(E\), which is treated separately.

This is a minor correction but worth making. 
GitHub

4.11. Sensitivity Proposition S1: Bounded nonlinearity
This is a good result.

The proposed relaxation
\[
E[\beta\mid V]
=
M_H(X)+\Delta(X)
\]

with
\[
X^\top\Delta(X)=0,
\qquad
\Delta(X)\mathbf 1=0,
\qquad
\|\Delta(X)\|_\infty\leq\eta
\]

is a natural way to bound the consequences of observationally invisible nonlinearities.

The logic is correct:
\[
E[T\mid X]
=
X^\top QX+X^\top\Delta(X)
=
X^\top QX.
\]

Thus the identified symmetric component remains unchanged.

Probability validity implies
\[
M_H(x)\geq-\eta
\]

elementwise, which produces the relaxed outer constraints.

The estimand itself differs from its linear-model counterpart by at most \(\eta\) for a single transition rate, because the weights sum to one.

For a general linear functional, the proposed widening
\[
\eta\sum_{rc}|a_{rc}|
\]

is correct.

I would retain this proposition.

One qualification: the result is an outer bound under bounded deviations, not a characterization of the sharp identified set under the relaxed model.

That distinction should be maintained in the presentation. 
GitHub


5. Cross-phase consistency issues
Reading all four documents together revealed several inconsistencies that I would resolve before submission.
Issue	Location	Required correction
General covariates versus district conditioning	P1 §§1.3, 1.5; P2b §0	Reconcile conditioning on \(Z\) with the specified mean-independence restriction.
Overidentifying restriction count	P1 §1.6; P3 §5	Use the corrected count \((D^\star-1)(C-1)R(R-1)/2\).
Pooled M0 conditional expectation	P2b Proposition D4(b)	Condition on both \(X\) and \(N\), not just \(X\).
Probability-validity assumption	P2a §0.2, S3	Correct the subset-versus-equality statement.
Turnover sensitivity	P1 A1; P3 §7	P1 says turnover is treated in Phase 3, but P3 explicitly leaves it open.
Cluster dependence	P1 A6; P3 P3-S1 and P3.5	Distinguish regression-error inference from the finite-population residual allowance.
Asymptotic regime	P1 A10; P3 §§0, 6	Separate fixed-dimensional results from growing-district results.

I want to elaborate on two of these.
5.1. The covariate specification in Phase 1 is still not fully reconciled
Phase 1 initially defines
\[
\mathcal G=\sigma(X,N,Z)
\]

and states
\[
E[\beta\mid X,N,Z]
=
E[\beta\mid X,Z].
\]

But the subsequent contextual-effects model specifies
\[
E[\beta\mid X,N,Z]=M_H(X).
\]

The former allows the conditional mean to depend on \(Z\); the latter does not.

Phase 2a resolves this by explicitly working without \(Z\), and Phase 2b develops the discrete district extension.

That is fine as a development strategy, but the paper should not imply that the general covariate extension has already been proved.

If the intended scope includes continuous covariates such as \(\log N\) or census characteristics, that remains a separate extension.

GitHub+2
5.2. The turnover result is not actually established in Phase 3
Phase 1 A1 says that Phase 3 treats bounded turnover,
\[
\|\varepsilon_i\|_1\leq\epsilon,
\]

and widens the bounds accordingly.

Phase 3 §7, however, says that turnover requires relaxing both the Duncan–Davis sets and the regression identity and is left open.

I agree with the more cautious statement in Phase 3.

The two changes are linked: if the accounting identity fails, both the unit-level feasible transportation polytope and the observable quadratic regression identity change.

It is not sufficient simply to widen the final estimand interval while retaining the original identified \(Q\).

I recommend describing turnover as an open extension, rather than a completed Phase 3 result. 
GitHub+1


6. Consolidated recommendations to the author
I would prioritize the revisions as follows.
Priority	File / location	Recommendation
Critical	P2b D4(b)	Correct the conditioning error: M0 requires the conditional district shift given \(X,N\) to be affine in \(X\) and independent of \(N\).
Critical	P3 P3.5(a)	Restrict finite-population CI coverage to independent unit residuals, or derive a cluster-level residual allowance.
Critical	P3 P3.5(b)	Do not claim formal CLR validity for the simulated candidate set until candidate recovery and coverage are proved.
Critical	P3 §§0, 6	Restrict established asymptotics to the regimes supported by the sampling assumptions; develop separate R2/R3 theory.
Major	P3 P3.2–P3.4	Distinguish sample-hull-specific asymptotics from inference for the population-hull endpoint.
Major	P3 §0.2	Distinguish the superpopulation, finite-population and finite-design conditional-mean targets.
Major	P3 P3.2	State Slater as a strict-interiority condition and acknowledge boundary cases.
Major	P3 P3.4	Qualify the negative-mean result in terms of distinct projected gradients.
Major	P2b D1.1 / P3 §1	State the precise rank conditions for M1 and M2 estimation.
Medium	P1 §1.6	Correct the count of independent M1 overidentifying restrictions.
Medium	P1 A1 / P3 §7	Reconcile the treatment of bounded turnover.
Medium	P1 A3 / P2b	Reconcile the continuous-covariate and district-effects specifications.
Medium	P2a S3	Correct the probability-validity qualification.
Medium	P2a Proposition 7	Correct the joint value-function regularity claim.
Low	P3 simulations	Report Monte Carlo standard errors and increase the number of replications for coverage evaluation.


7. What I now regard as established
Given the additional material, I would update the original assessment as follows.
Result	Revised assessment
R×C identification decomposition	Established
Dimension of unidentified contextual effects	Established
Conditional-selection characterization of the sharp set	Established, subject to stated functional-analytic conditions
District-specific identification decomposition	Established
Identification of district baseline shifts under M1	Established
M1 overidentifying restrictions	Established
Aggregate estimand under M1	Established
M1/M1(\(\delta\))/M2 nesting	Established
Two-row district closed form	Established
M2 Minkowski-sum result	Established
L4 and L5 identification results	Established
D4(a), equivalence when district shifts vanish	Established
D4(b), general characterization of pooled-model validity	Incorrect conditioning
Sample-hull one-sided conservatism	Established
Root-n plug-in estimation under Slater	Established under fixed-dimensional conditions
Directional delta-method result	Established for fixed design/hull
Conservative CI for \(\widetilde B\)	Established under the stated Gaussian approximation
Conservative CI for \(B\) under independent residuals	Established under additional sampling conditions
Conservative CI for \(B\) with arbitrary within-district residual dependence	Not established
CLR inference with the complete dual candidate set	Theoretical route is appropriate, subject to conditions
Simulated CLR candidate-generation algorithm	Formal justification outstanding
General inference in R2/R3	Requires additional asymptotic theory
Bounded nonlinear sensitivity result	Established as an outer bound
Bounded-turnover sensitivity	Open


8. Final recommendation
I remain positive about this project, and the additional phases make me more convinced that there is a substantial methodological paper here.

In fact, I think the overall contribution is becoming clearer.

The paper is developing three complementary methodological results:

First: A general theory of contextual-effects identification.

Phase 2a shows that the unidentified component of the linear contextual-effects model has a precise algebraic structure in the \(R\times C\) setting.

Second: A theory of partial pooling across contextual environments.

Phase 2b shows how restrictions across districts can sharpen identification without necessarily requiring common district-specific transition probabilities. The assumption ladder is an especially useful contribution because it gives researchers a way to study how much identification depends on the strength of cross-district homogeneity assumptions.

Third: A framework for inference on partially identified ecological transition rates.

Phase 3 addresses the practical difficulties of estimated constraints, nonsmooth LP endpoints, finite-population variation and many district-specific inequalities.

That combination is potentially quite valuable.

I would not recommend abandoning any of these directions. However, I would distinguish between results that are now ready to be presented as theorems and results that remain promising but provisional.
My suggested strategy for the next revision
I would make five decisions before moving to implementation or a full manuscript draft.

1. Correct D4(b) first. This is a substantive issue in the district model, not just an inference qualification.

2. Narrow the formal inference claims to R1 initially. Establish the fixed-dimensional results cleanly for M1 and M2, with an explicit independence assumption. Treat clustered M0 separately.

3. Make Theorem P3.5(a) the principal confidence-interval result. Its construction is relatively transparent and can support a strong contribution even without the simulated CLR procedure.

4. Retain CLR as an important computational extension, but distinguish the exact theoretical procedure from the candidate-generation approximation. I think the latter is worth pursuing, but it should not be represented as fully proved yet.

5. Make the sample-hull argument a central part of the inferential framework. It is a strength of the project that using the empirical hull is conservative. But the distinction between sample-hull inference and population-hull inference should be explicit throughout.

Bottom line: I would characterize Phase 2a and most of Phase 2b as mathematically successful, with one substantive error in the characterization of pooled-model validity. Phase 3 contains a strong and appropriate inferential strategy, but its formal validity is narrower than currently claimed, particularly for clustered observations and growing numbers of districts.

I think the project is worth pursuing as a methodological journal article. The revisions I am recommending are intended to make the claims more defensible, not to diminish the underlying contribution.

The most important thing is to avoid presenting the promising R2/R3 CLR implementation as a general theorem before its asymptotic foundations are established. The fixed-dimensional results, the district identification theory, and the conservative inference construction already provide a substantial contribution.
