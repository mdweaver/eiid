# Referee report, round 1: Phases 1 and 2a

*Received Nov 2026 from the original referee of Phases 1–3 (round 1). Transcribed as received (pasted text); formatting from the original may be lost. The response plan is `../REVIEW_RESPONSE_PLAN.md`; the response is `../REVIEW_RESPONSE.md`.*

Reviewer Report
Manuscript: Extension of linear contextual-effects partial identification bounds in ecological inference to the \(R\times C\) setting

Materials reviewed: Phase 1: Structural Foundations and Assumptions; Phase 2a: \(R\times C\) Bounds, Derivations and Proofs.

Recommendation: Major revision, with a positive assessment of the underlying mathematical contribution.
1. Overall assessment
The proposed extension is mathematically promising. In particular, I think the move from a single unidentified contextual-effect parameter in the \(2\times2\) case to a structured, higher-dimensional unidentified parameter space in the \(R\times C\) case is well motivated and, for the most part, correctly developed.

My overall assessment is that the central algebraic generalization and the main sharp-identification theorem are substantially correct under the stated model. I have not found a fundamental error in the construction of the unidentified parameter space, the expression for the estimand, or the principal conditional-selection argument.

There are, however, several issues that I would want resolved before treating these results as publication-ready.

My principal concerns are:
1. A discrepancy between the assumptions in Phase 1 and those actually used in Phase 2a, particularly concerning covariates and conditioning.
2. An incorrect or insufficiently qualified claim about when the probability-validity assumption (S3) follows automatically.
3. An important support condition that needs to be carried through to the propositions providing bounds on the estimand.
4. An incorrect characterization of the regularity of the parametric LP value function.
5. Insufficient precision concerning the dimensionality of \(A\) and the compactness of the identified sets.
6. A distinction between population identification and what can actually be learned from a finite sample that needs more explicit treatment.
I discuss these below, distinguishing errors in mathematical claims from issues of exposition or scope.

2. Major comments
2.1. The central identification decomposition is correct
Phase 2a, §2, Lemma 1.

This is the most important algebraic result, and I find it convincing.

The decomposition
\[
H_c=Q_c+K_c,
\]

where
\[
Q_c=\frac12(H_c+H_c^\top),\qquad
K_c=\frac12(H_c-H_c^\top),
\]

is exactly the right parameterization.

The observable regression identifies
\[
E[T_c\mid X=x]=x^\top H_cx=x^\top Q_cx,
\]

while the skew-symmetric component disappears from the quadratic form.

The argument establishing
\[
\sum_c H_c=J
\]

is also correct, provided the stated richness assumption A5 holds on the relevant support. This gives
\[
\sum_cK_c=0.
\]

Consequently, the dimension
\[
\boxed{
\dim(\mathcal K)=(C-1)\frac{R(R-1)}2
}
\]

is correct.

The interpretation is particularly useful: for every pair of row categories, the symmetric component identifies the average of two contextual effects, while their difference remains unidentified.

This gives the following dimensions:
Configuration	Unidentified dimensions
\(2\times2\)	1
\(2\times3\)	2
\(3\times2\)	3
\(3\times3\)	6
\(4\times3\)	15

I would retain this result essentially as written.

One qualification: Phase 1, §1.3, initially argues that \(\sum_c H_c=J\) follows because the row sums equal one on an open subset of the simplex. That is not an assumption explicitly established there. Phase 2a correctly replaces this reasoning with the A5 argument. I would revise Phase 1 to use the latter, or explicitly state the required interior condition.

In the published treatment, A5 should also be stated in its matrix form:
\[
x^\top D x=0
\quad\forall x\in\operatorname{supp}(P_X)\cap A
\quad\Longrightarrow\quad D=0
\]

for symmetric \(D\).

This is important because the uniqueness of \(Q\) is not merely a consequence of observing variation in \(X\). It requires a specific quadratic rank condition. In particular, for \(R=3\), having three distinct compositions is not sufficient. The support must be rich enough to identify all six coefficients of the quadratic regression.

Assessment: Correct, subject to the explicit rank assumption.

GitHub

2.2. The estimand is correctly expressed as an affine function of the unidentified parameters
Phase 2a, §3, Lemma 2.

The expression
\[
B^*_{rc}
=
\frac{1}{\mu_r}
\left[
\sum_s Q_{c,rs}\mu_{rs}
+
\sum_{s\neq r}K_{c,rs}\mu_{rs}
\right]
\]

is correct under the stated condition \(P(X\in A)=1\).

The derivation is straightforward:
\[
E[NX_r\beta_{rc}]
=
E\left[NX_r E[\beta_{rc}\mid X,N]\right],
\]

followed by substitution of the contextual-effects restriction.

The weighting by \(NX_r\) is handled correctly, as are the signs associated with the skew-symmetric entries.

This is an important result because it means that any linear functional of the transition rates can be optimized over a convex identified set.

There is one minor mathematical wording issue.

The document writes
\[
\frac{\partial B^*_{rc}}{\partial K_{c,rs}}
=
\frac{\mu_{rs}}{\mu_r}.
\]

Strictly speaking, the entries of \(K\) are not independent coordinates: skew symmetry and the cross-column restriction impose constraints. I would describe this instead as the coefficient on the corresponding entry in the affine representation, or specify that the derivative is taken in an ambient parameterization before imposing the restrictions.

This becomes especially important in the \(R\times C\) case, where the natural free parameters are the entries for \(c<C\) and \(r<s\).

Assessment: Correct. A minor refinement of the derivative notation is advisable.

GitHub

2.3. Theorem 1: The sharp-identification argument is convincing
Phase 2a, §4, Theorem 1.

This is the most mathematically demanding part of the proposal, and it is where I spent the most time.

My assessment is positive: the conditional-selection argument appears valid. I do not see a fatal gap in the proof of the equivalence between feasibility of a contextual-effects parameter and the conditional support-function inequalities.

The underlying argument is as follows.

For each observation \((X,T)\), define the feasible Duncan–Davis set
\[
S(X,T)=\{b\in\mathfrak B:b^\top X=T\}.
\]

A candidate contextual-effects parameter \(H=Q+K\) is feasible if there exists a measurable selection \(\beta(W)\in S(X,T)\) such that
\[
E[\beta\mid V]=M_H(X),
\qquad V=(X,N).
\]

The theorem establishes that this is equivalent to
\[
\langle u,M_H(x)\rangle
\leq
E[h_S(u;X,T)\mid V=v]
\]

for every direction \(u\), almost everywhere in \(v\) with \(x\in A\), together with S3.

This is a natural application of conditional selection expectations and support-function characterizations of closed convex sets.

I checked the main components of the proof.

First, the necessary condition is correct.

Since
\[
\beta\in S(X,T),
\]

we have, for every \(u\),
\[
\langle u,\beta\rangle
\leq h_S(u;X,T).
\]

Taking conditional expectations gives the desired support-function inequality.

Second, the common-null-set argument is correct.

The proof initially establishes the inequality for each fixed \(u\), which produces a potentially different null set for each direction.

The use of a countable dense subset of \(\mathbb R^{R\times C}\), followed by continuity and the uniform Lipschitz bound on \(h_S\), legitimately extends the result to all directions outside one common null set.

This is an important technical detail, and it is handled correctly.

Third, the sufficiency argument is appropriate.

The proof defines the set of measurable selections
\[
\mathcal S
=
\{\xi\in L^1:\xi(W)\in S(X,T)\text{ a.s.}\}
\]

and the conditional-expectation image
\[
\mathcal C
=
\{E[\xi\mid V]:\xi\in\mathcal S\}.
\]

The argument that \(\mathcal S\) is convex, norm-closed and uniformly integrable, and therefore weakly compact, is sound.

Conditional expectation is a bounded linear operator, so its image \(\mathcal C\) is also weakly compact and convex.

Separation then provides the appropriate converse: if the candidate mean is outside \(\mathcal C\), some bounded measurable direction \(g(V)\) separates it from every feasible conditional expectation.

The measurable maximum argument correctly identifies the corresponding support function.

I particularly appreciate the construction of
\[
\eta
=
M_H(X)1\{X\in A\}
+
E[\xi^0\mid V]1\{X\notin A\},
\]

where \(\xi^0\) is an arbitrary measurable feasible selection outside \(A\).

This handles the fact that the contextual-effects restriction is imposed only on \(A\), rather than unnecessarily requiring it outside \(A\).

Fourth, the construction of a witness distribution is correct.

Once the proof obtains a selection \(\xi^\star\) satisfying
\[
E[\xi^\star\mid V]=\eta,
\]

setting
\[
\beta=\xi^\star(W)
\]

constructs a joint law with the required observable marginal and conditional mean restriction.

I agree with Remark 4.1: auxiliary randomization is unnecessary. Because the Duncan–Davis set is convex, the relevant conditional expectation can be realized by a measurable selector.

This is a useful correction to Phase 1.

GitHub+4
What I would improve in this proof
Although I believe the argument is correct, it should be made more accessible and self-contained.

In particular, I would explicitly state the following facts rather than relying on the sequence of justification tags:
1. The measurable-selection set \(\mathcal S\) is a subset of vector-valued \(L^1(P)\).
2. Its boundedness in \(L^\infty\) implies uniform integrability.
3. Norm-closedness and convexity imply weak closedness.
4. Dunford–Pettis therefore gives weak compactness.
5. The conditional-expectation operator is weak-to-weak continuous.
6. The support function of \(\mathcal C\) is obtained using a measurable maximizing selector.
The references to measurable maximum and separation theorems should be checked against the exact hypotheses and theorem numbering before publication. The document itself appropriately identifies some references as requiring verification.

I would also explicitly explain why the support-function inequality must hold for all directions simultaneously, rather than simply saying that the inequality holds for each \(u\).

These are improvements to presentation, rather than objections to the underlying proof.

Assessment: I regard Theorem 1 as a valid and potentially important result.

2.4. An important error concerning assumption S3
Phase 2a, §0.2; Lemma 3(c).

This is the first place where I believe the text contains an actual mathematical error.

The document states that
Lemma 3(c) shows this is implied when \(A\subseteq\operatorname{conv}\operatorname{supp}P_X\).

However, Lemma 3(c) proves the result under the condition
\[
\boxed{A=\operatorname{conv}(\operatorname{supp}P_X).}
\]

These are not equivalent conditions.

The distinction matters because S2 is imposed only where \(X\in A\).

The proof of Lemma 3(c) proceeds by establishing
\[
E[\beta\mid V]\in\mathfrak B
\]

almost surely, then extending the validity of \(M_H(x)\) by continuity to the support of \(X\), and finally using linearity and convexity to extend it to the convex hull of the support.

This establishes S3 over the convex hull of the support. It does not establish S3 over an arbitrary convex subset of that hull when the contextual-effects restriction is only imposed on that subset.

For example, suppose the support of \(X\) is sufficiently rich, but \(A\) is a singleton located strictly inside its convex hull and not in the support. Then \(A\subseteq\operatorname{conv}(\operatorname{supp}P_X)\), but S2 is imposed on a probability-zero event. There is no reason for the extrapolated value \(M_H(x)\) at that singleton to be a valid transition-probability matrix.

The proof does not rule this out.
Recommended correction
Change the statement in §0.2 to something like:
If \(A=\operatorname{conv}(\operatorname{supp}P_X)\), S3 follows automatically from S1–S2. More generally, the result holds when \(A\subseteq\operatorname{conv}(\operatorname{supp}P_X\cap A)\), provided the contextual-effects restriction is imposed on that support.

Alternatively, retain S3 as a separate assumption whenever the contextual-effects model is imposed on a region extending beyond the convex hull of the support on which it is actually identified.

This also requires revising the sentence that describes S3 as a genuine extra assumption only when \(A\) extrapolates beyond the data. The more precise distinction is whether the probability-validity region is covered by the convex hull of the observations to which S2 applies.

GitHub+1

2.5. The support condition for the estimand must be carried through the propositions
This is an important issue of scope that appears in several places.

Lemma 2 explicitly assumes
\[
P(X\in A)=1.
\]

That assumption is necessary for the derivation of the estimand as an affine function of \(K\).

If observations occur outside \(A\), then the contextual-effects model does not determine their conditional transition means.

The estimand becomes
\[
B^*_{rc}
=
\frac{1}{\mu_r}
\left\{
E[NX_rM_{H,rc}(X)1\{X\in A\}]
+
E[NX_r\beta_{rc}1\{X\notin A\}]
\right\}.
\]

The second term is not generally a function of \(K\).

The authors recognize this in Remark 2.1, which is good. However, the qualification is not consistently carried into the subsequent results.

In particular:

Proposition 7 (§6.5).

The proposition provides an LP for bounding a linear functional of \(B^*\), using the affine representation from Lemma 2. Its statement should explicitly impose \(P(X\in A)=1\), or refer to the more general bounds incorporating the observations outside \(A\).

The same qualification should be made wherever Proposition 6 uses the affine expression for \(B^*\) with a restricted interval \(A=[l,u]\).

Proposition 5 avoids this issue because \(A=\Delta^{R-1}\) necessarily contains the support.
Recommended correction
Either:
1. State throughout §6 that \(A\) contains the support of \(X\); or
2. Develop the aggregate bound that adds the model-free Duncan–Davis contribution from observations outside \(A\).
I would favor the first option for Phase 2a, while retaining Remark 2.1 as the explicit transition to Phase 2b.

This is more than a technical qualification. Without it, the document risks claiming that the LP gives sharp bounds for a superpopulation estimand even when some of the population lies outside the region where the identifying restriction is imposed.

GitHub+2

2.6. The claimed regularity of the LP value function is incorrect
Phase 2a, Proposition 7, concluding paragraph.

I believe there is a genuine mathematical error in the statement:
The value function \(Q\mapsto\max\{\gamma^\top k:Gk\le g(Q)\}\) is concave and piecewise linear in \((Q,\mu)\).

This is not generally true when both the objective coefficients and the constraint right-hand sides vary.

Recall that
\[
\gamma_{c,rs}
=
\mu_{rs}
\left[
\frac{a_{rc}-a_{rC}}{\mu_r}
-
\frac{a_{sc}-a_{sC}}{\mu_s}
\right].
\]

Thus, the objective coefficients themselves depend on the moments \(\mu\), generally through ratios.

There are two distinct parametric-LP results that should not be conflated:
- Holding the objective coefficients fixed, the optimal value of a maximization LP is concave and piecewise linear in its right-hand-side parameters, when the relevant feasibility and boundedness conditions hold.
- Holding the feasible set fixed, the optimal value is a support function of the objective coefficients. It is convex and piecewise linear in those coefficients.
Neither result establishes joint concavity or joint piecewise linearity when both the objective and constraints vary.
Simple illustration
Consider
\[
\max_k \ \mu k
\]

subject to
\[
-Q\leq k\leq Q,
\]

with \(\mu,Q>0\).

The optimal value is
\[
V(Q,\mu)=\mu Q.
\]

This is bilinear, not jointly piecewise linear or jointly concave.

The same issue arises in the proposed ecological-inference LP.
Recommended correction
I would replace the current assertion with:
For fixed \(\mu\), the optimal value is concave and piecewise linear in \(Q\). For fixed \(Q\), it is convex and piecewise linear in the objective coefficients \(\gamma\). Jointly in \((Q,\mu)\), the value function is generally neither concave nor piecewise linear.

You can still establish local smoothness when the optimal basis remains locally unchanged and the relevant nondegeneracy and regularity conditions hold.

However, this needs to be established for the joint parameterization, taking account of the dependence of \(\gamma\) on \(\mu\).

This correction is particularly important because this paragraph motivates the subsequent inference theory in Phase 3.

GitHub+1

2.7. The dimensionality condition for compactness needs to be stated more precisely
Phase 2a, Theorem 1 and Lemma 3(b).

The document states that \(\Theta_O\) is bounded when \(A\) has nonempty relative interior.

The proof, however, uses the existence of \(R\) linearly independent points in \(A\), obtained from a relatively open subset of the simplex.

This requires \(A\) to be full-dimensional relative to the simplex.

Under the usual definition, the relative interior of a nonempty convex set is its interior relative to its own affine hull. A singleton, for example, has nonempty relative interior in this sense.

That condition alone does not imply boundedness.

To see the problem, consider \(R=3\), \(C=2\), and
\[
A=\{x_0\},
\qquad
x_0=(1/3,1/3,1/3).
\]

Take
\[
Q_1=Q_2=\frac12J
\]

and let
\[
K_1=tK,\qquad K_2=-tK,
\]

where
\[
K=
\begin{pmatrix}
0&1&-1\\
-1&0&1\\
1&-1&0
\end{pmatrix}.
\]

Notice that
\[
Kx_0=0.
\]

Therefore,
\[
(Q_c+K_c)x_0
=
\frac12\mathbf 1
\]

for every value of \(t\).

All these parameter values satisfy the outer-set constraints, so \(\Theta_O\) is unbounded.

Yet \(A\) has nonempty relative interior in its own affine hull.
Recommended correction
State the condition as
\[
\operatorname{affdim}(A)=R-1,
\]

or equivalently that \(A\) has nonempty interior relative to the simplex \(\Delta^{R-1}\).

The proof in steps 5–6 of Lemma 3(b) then works.

There is an additional point worth making: under A5 and the default choice \(A=\operatorname{conv}(\operatorname{supp}P_X)\), the required full dimensionality follows from the richness condition. But this relationship should be stated explicitly rather than leaving the compactness condition ambiguous.

This also affects the use of the word polytope. A finite intersection of half-spaces is a polyhedron, but it is not necessarily bounded. The authors should establish boundedness under the stated conditions before describing \(\Theta_O\) as a polytope.

GitHub

3. Important issues concerning the identification framework
These are not necessarily errors in Phase 2a's main proof, but they matter for the coherence of the proposed extension.
3.1. The treatment of covariates in Phase 1 is inconsistent
Phase 1 introduces covariates \(Z\) and defines the conditioning sigma-algebra as
\[
\mathcal G=\sigma(X,N,Z).
\]

It then writes the structural restriction as
\[
E[\beta\mid\mathcal G]=M_H(X).
\]

However, the corresponding A3 assumption is
\[
E[\beta\mid X,N,Z]
=
E[\beta\mid X,Z].
\]

This eliminates dependence on \(N\), conditional on \(X,Z\), but it does not eliminate dependence on \(Z\).

Similarly, the pooled model M0 is initially formulated as
\[
E[\beta_{irc}\mid X,N]
=
e_r^\top H_cX,
\]

not as a conditional expectation given \(X,N,Z\).

Consequently, the claim that S2 follows from A2 and A3 is not quite correct as written.

There are two possible formulations:

Option A: No covariate effects in the conditional mean.

Explicitly assume
\[
E[\beta\mid X,N,Z]=M_H(X).
\]

This imposes mean independence from both \(N\) and \(Z\).

Option B: Covariate-dependent conditional means.

Define an expanded contextual-effects model such as
\[
E[\beta\mid X,N,Z]
=
M_{H}(X,Z).
\]

The identified quadratic structure and unidentified parameter space would then need to be developed for that expanded model.

Phase 2a is internally consistent because it removes \(Z\) and works with \(V=(X,N)\). I would keep that distinction clear. In particular, I would not describe the Phase 2a theorem as already establishing the covariate extension.

GitHub+3
3.2. The distinction between \(\Theta_I\), \(\Theta_X\), and \(\Theta_O\) is essential
I think this is one of the strongest conceptual features of the proposal.

The nesting
\[
\Theta_I\subseteq\Theta_X\subseteq\Theta_O
\]

is correct.

The sets correspond to three different informational environments:
Set	Conditioning information	Interpretation
\(\Theta_I\)	\(X,N\)	Sharp set under the full conditional-mean restriction
\(\Theta_X\)	\(X\)	Sharp set when only composition-conditional information is used
\(\Theta_O\)	Implied regression and probability validity	Outer set obtained by the Jensen/probability-validity argument

This distinction is not merely terminological.

In particular, the Phase 2a LP in Proposition 7 computes bounds over \(\Theta_O\), not generally over the sharp identified set \(\Theta_I\).

The sharp set requires the conditional support-function restrictions in Theorem 1. Those restrictions incorporate the distribution of \(T\) conditional on \(X,N\), and need not reduce to finitely many linear constraints on \(K\).

The authors understand this distinction, and I strongly recommend preserving it throughout the paper.

One refinement: Proposition 2 establishes that \(\Theta_O\) is sharp when only the regression function is retained, rather than when the complete observable law is retained. This is an important qualification that should accompany any claim that \(\Theta_O\) is sharp.

Similarly, Proposition 3 provides a sufficient condition for
\[
\Theta_I=\Theta_O.
\]

It does not provide a necessary and sufficient characterization. The document correctly calls this a sufficient condition, and I would retain that language.

GitHub+2
3.3. The meaning of the population support and the sample hull needs attention
This is perhaps the most important issue for the eventual empirical application.

The theoretical model is formulated in terms of a probability law \(P\), and the target is
\[
B^*_{rc}
=
\frac{E[NX_r\beta_{rc}]}{E[NX_r]}.
\]

The natural default for the contextual-effects region is
\[
A=\operatorname{conv}(\operatorname{supp}P_X).
\]

The computational discussion, however, frequently uses
\[
A_n=\operatorname{conv}\{X_1,\ldots,X_p\},
\]

the observed sample hull.

These are not generally identical.

For a finite sample from a continuous population distribution, the sample hull ordinarily does not contain the full population support.

Consequently, a constraint imposed over the sample hull is not automatically a population assumption valid over the entire support needed to derive the formula for \(B^*\).

There are two distinct uses of the sample hull:
1. As the support hull of the empirical distribution, in a finite-sample identification analysis.
2. As a plug-in estimate of the population support hull, in estimating population bounds.
These should not be conflated.

The first is mathematically straightforward. The second requires assumptions and asymptotic arguments concerning support estimation and the behavior of the hull vertices.

This is particularly relevant because the outline already recognizes that the vertices are estimated and that outlying precincts can have disproportionate effects on the constraints.

I would recommend making the distinction explicit in Phase 2a, even if the full treatment is deferred to Phase 3.

GitHub+2
3.4. The sharpness claim should be explicitly defined relative to the continuous probability model
The latent transition probabilities are introduced as
\[
\beta_{irc}=\frac{N_{irc}}{N_{ir}}.
\]

With integer cell counts, these probabilities lie on discrete grids determined by the row denominators \(N_{ir}\).

The mathematical model in Phase 2a instead takes
\[
\beta\in(\Delta^{C-1})^R,
\]

which allows all real-valued probability vectors.

This is standard and reasonable for a superpopulation ecological-inference model, but it is not exactly the same as maintaining integer restrictions on every unobserved cell count in a finite population.

Thus, the sharpness claim should be framed as:
Sharpness relative to the continuous probability model of latent transition rates satisfying the accounting identity and the stated conditional-mean restrictions.

It should not be interpreted as demonstrating sharpness of the exact finite-integer feasible set.

I do not regard this as a defect in the proposed theory. I regard it as a necessary clarification of the estimand and model.

4. Assessment of the closed-form results
I reviewed the algebra in §6 separately.

Overall, I find the derivations convincing.
4.1. Proposition 4: Equivalence with JKST
The \(2\times2\) mapping is correct.

In particular, the derivation
\[
K=\frac{c_1}{2}-w_1
\]

and the reconstruction of the JKST restrictions from the vertex constraints are algebraically consistent.

The expression for \(B^*(K)\) also matches the population analogue of JKST equation (16).

I found no substantive error here.

One important qualification: the bounds are bounds over the outer set, not necessarily the sharp set under the richer conditional information in Theorem 1. The distinction should remain explicit when describing the equivalence with JKST.

GitHub
4.2. Proposition 5: Whole-simplex case
This is a particularly useful result.

When \(A=\Delta^{R-1}\), the constraints reduce to
\[
|K_{c,rs}|\leq Q_{c,rs}.
\]

Together with
\[
\sum_cK_{c,rs}=0,
\]

the coordinate projection becomes
\[
K_{c,rs}\in[-m_{c,rs},m_{c,rs}],
\]

where
\[
m_{c,rs}=\min\{Q_{c,rs},1-Q_{c,rs}\}.
\]

The resulting interval for \(B^*_{rc}\) is correct.

I checked the logic of the independence argument: restrictions associated with different row pairs are separable, and for a given row \(r\), the unidentified parameters entering \(B^*_{rc}\) belong to different row pairs.

Consequently, the endpoints can be attained simultaneously.

This is a sound result.

I would make one minor point clearer: the expression involving \(m_{c,rs}\) presupposes a nonempty feasible set. If the identified quadratic coefficients violate the necessary probability restrictions, the relevant interval and feasible set are empty.

GitHub
4.3. Proposition 6 and Lemma 4
The box-intersected-with-hyperplane argument is correct.

In particular, the coordinate projection
\[
\left[
\max\left\{\ell_c,-\sum_{c'\ne c}h_{c'}\right\},
\;
\min\left\{h_c,-\sum_{c'\ne c}\ell_{c'}\right\}
\right]
\]

is correct.

The greedy optimization rule is also correct: beginning at the lower bounds and allocating the required remaining mass to coordinates in descending order of their objective coefficients is the appropriate solution.

I see no mathematical issue here.

GitHub
4.4. Corollary 6.1: Joint versus collapsed analysis
I consider this a useful contribution, particularly for the \(2\times3\) application.

The result that the joint interval is contained in the collapsed interval is supported by the algebra.

The sufficient condition for equality is also correctly motivated.

I would, however, sharpen the language in part (c).

Rather than saying that strict tightening requires "strong contextual effects," state the precise condition:
Strict tightening requires the failure of the common-sign condition in part (b); given the sum restrictions on the identified coefficients, this implies that at least one of the relevant pure-poll extrapolated rates is negative.

This is more precise and avoids an undefined claim about the magnitude of contextual effects.

The result that the joint analysis can produce strict tightening even when all extrapolated rates are valid for \(R\geq3\) is interesting. The authors appropriately leave a general characterization open.

GitHub
4.5. Proposition 7: The general \(R\times C\) LP
The parameterization of \(\mathcal K\) is correct.

The number of free parameters is
\[
d=(C-1)\frac{R(R-1)}2.
\]

The construction of the constraint matrix \(G\) is correct, and I also checked the signs in the objective coefficient
\[
\gamma_{c,rs}
=
\mu_{rs}
\left[
\frac{a_{rc}-a_{rC}}{\mu_r}
-
\frac{a_{sc}-a_{sC}}{\mu_s}
\right].
\]

The signs correctly reflect both the skew symmetry of \(K_c\) and the requirement \(K_C=-\sum_{c<C}K_c\).

This is a strong general computational result.

My principal objection is the value-function regularity claim discussed in §2.6, rather than the LP construction itself.

I would also emphasize in the proposition that this is a finite-dimensional LP for the outer set. It should not be presented as a finite-dimensional LP representation of the general sharp identified set.

GitHub

5. The untestability result is important and convincing
Phase 2a, Proposition 9.

I think this result deserves greater prominence in the eventual paper.

The proposition establishes that a nonlinear deviation
\[
\Delta(X)
\]

can satisfy
\[
X^\top\Delta(X)=0,
\qquad
\Delta(X)\mathbf 1=0,
\]

while preserving valid transition probabilities.

The resulting transition structure
\[
\beta=M_H(X)+\Delta(X)
\]

generates exactly the same observable distribution as
\[
\beta'=M_H(X).
\]

Nevertheless, the aggregate transition rates differ by
\[
B^*_{rc}(\beta)-B^*_{rc}(\beta')
=
\frac{E[NX_r\Delta_{rc}(X)]}{\mu_r}.
\]

This is a convincing demonstration that the linear contextual-effects restriction is not generally testable from the observed margins.

The argument also provides an important interpretation of the unidentified skew-symmetric component: it represents only the linear part of a much larger space of observationally invisible deviations.

I would retain this result and consider using it earlier in the paper to motivate why the contextual-effects assumption must be treated as a substantive structural restriction rather than something fully validated by goodness-of-fit tests.

One wording qualification: the characterization of linear \(\Delta\) as exactly \(M_{K'}\) with \(K'\in\mathcal K\) is correct when the two restrictions on \(\Delta\) hold throughout the simplex. The domain on which those restrictions are imposed should be explicit.

GitHub

6. Additional technical and expositional recommendations
I would make the following changes before submission.
Priority	Location	Recommended change
High	P2 §0.2, S3	Correct the \(A\subseteq\operatorname{conv}(\operatorname{supp}P_X)\) statement.
High	P2 Proposition 7	Correct the joint LP value-function regularity claim.
High	P2 Lemma 2, Propositions 6–7	Explicitly carry the \(P(X\in A)=1\) condition wherever \(B^*\) is expressed as an affine function of \(K\).
High	P1 §1.5 and A3	Reconcile conditioning on \(Z\) with the stated mean-independence assumption.
Medium	P2 Lemma 3(b)	Specify full dimensionality relative to the simplex.
Medium	P1 §1.3	Use A5 rather than an unstated open-subset argument to establish \(\sum H_c=J\).
Medium	P2 Theorem 1	Expand the functional-analytic proof and verify the precise cited theorem conditions.
Medium	P2 Proposition 7	Distinguish the outer-set LP from the sharp identified set.
Medium	P1/P2 sampling framework	Distinguish the population support hull from the empirical sample hull.
Low	P2 Lemma 2	Replace unconstrained partial-derivative language with coefficients in a free parameterization.
Low	P2 Corollary 6.1	Replace qualitative language about strong contextual effects with the precise sign condition.

I would also make the following conceptual point more explicit in the paper.

Theorem 1 and Proposition 7 are doing fundamentally different things.

Theorem 1 establishes sharp identification, but requires conditional distributional information about \(T\) given \(X,N\).

Proposition 7 establishes an accessible finite-dimensional LP for the outer identified set, using the quadratic regression and probability-validity restrictions.

That is an important methodological trade-off. The sharp set is theoretically valuable, while the outer set is computationally accessible and extends the existing JKST implementation.

The paper should be explicit about what information is required to implement each.

In particular, if \(X,N\) are effectively unique for every precinct, the empirical distribution does not provide repeated observations from which to estimate \(P(T\mid X,N)\) nonparametrically. Estimating the sharp set in practice will therefore require additional structure, such as smoothing, conditional distribution estimation, or pooling over relevant composition/size neighborhoods.

I would regard this as an important problem for the empirical implementation, although it does not invalidate the population identification theorem.

7. What I would accept as established
It is useful to separate the results that I think have been proved from those that need additional work.
Result	Assessment
Symmetric/skew decomposition of \(H_c\)	Established
Dimension of unidentified parameter space	Established
Identification of \(Q_c\) under A5	Established
Affine representation of \(B^*\) in \(K\)	Established
, with support condition
Conditional-selection characterization of \(\Theta_I\)	Appears established
Convexity and closedness of \(\Theta_I\)	Established
Compactness	Requires precise dimensionality condition
Nesting \(\Theta_I\subseteq\Theta_X\subseteq\Theta_O\)	Established
Sharpness of \(\Theta_O\) relative to the regression-function information	Established
Sufficient condition for \(\Theta_I=\Theta_O\)	Established
\(2\times2\) equivalence with JKST	Established
Whole-simplex closed form	Established
\(R=2\), arbitrary \(C\), closed form	Established
Joint-versus-collapsed comparison for \(R=2\)	Established
General \(R\times C\) outer-set LP	Established
Joint LP value-function concavity and piecewise linearity in \((Q,\mu)\)	Incorrect as stated
Nonlinear observationally invisible deviations	Established

The numerical verification is useful, especially for the closed forms and parameterization. But it should be treated as supplementary evidence. The 22 passing tests do not independently validate the measure-theoretic sufficiency argument, and the document correctly acknowledges that limitation.

GitHub

8. Final recommendation
I am positive about the underlying contribution.

The strongest element is not simply that the JKST bounds can be generalized to more than two rows or columns. It is that the proposed framework identifies the structure of the unidentified contextual effects in the general \(R\times C\) case and distinguishes three different levels of identification:
\[
\boxed{
\Theta_I\subseteq\Theta_X\subseteq\Theta_O.
}
\]

This distinction is theoretically meaningful. It also produces a compelling research agenda: characterize the sharp set, derive computationally convenient outer bounds, and determine when the latter coincide with the former.

I would particularly emphasize three contributions in an eventual article:

1. The general identification decomposition. The unidentified space has dimension \((C-1)R(R-1)/2\), with a clear interpretation in terms of unidentified directional contextual effects.

2. The conditional-selection characterization. This is the key theoretical result. It establishes what the full observable distribution can identify under linear contextual effects.

3. The general LP representation. This provides a practical extension of the JKST probability-validity bounds to arbitrary \(R\times C\) tables.

My main reservation is that the paper currently moves somewhat too quickly between the sharp population identified set, the outer set, and the bounds actually obtainable from a finite sample.

Those are three distinct objects, and maintaining the distinctions is essential to making the contribution credible.

My recommendation would be to retain the main theoretical results, correct the issues identified above, and proceed with the extension. I would not abandon the approach or substantially restructure the mathematics. The principal results are sufficiently promising that these revisions seem worthwhile.

The one issue I would resolve before developing Phase 3 is the precise relationship between the support of the population distribution, the contextual-effects region \(A\), and the finite-sample hull used by the estimator. That relationship determines whether the proposed computational bounds are estimating the population object that the identification theorem actually characterizes.


---
