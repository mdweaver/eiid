# Follow-up review of within-precinct parallel trends

**Verdict.** Keep within-precinct parallel trends in expectation (WPT) as an explicit choice alongside pooled parallel trends (PPT), with the same full veteran ATT as the target. E28's main identity and the local counterfactual rate limits implement the requested comparison. E29(a)–(c) give a valid conservative bridge from realised counts to that expected effect under the stated independence and count-row conditions. E30's level/change equivalence and numerical example stand. I would revise the pooled witness in E28 and the multi-post sensitivity formula in E29(d) before clearing these new sections. This is a focused review of §§5.8–5.9 and the three open questions.

[Opus reply](https://github.com/mdweaver/eiid/blob/5d6c10658d40423f4c6c01be7c764d21eb4764ce/docs/extension/reviews/PHASE2e_reply_to_independent_review.md) · [Theory reviewed at 44cd35f](https://github.com/mdweaver/eiid/blob/44cd35f5fcbf55f944aff700dbb44d9aa0e8fee4/docs/extension/PHASE2e_did_general.md#58-parallel-trends-within-each-precinct-in-expectation-e28e29)

## 1. Confidence for the expected ATT from one realised election pair

The expenditure of precision in E29 has a cause. If all voters in a precinct can experience a common shock, more ballots in that precinct need not give independent observations of its trend. The same problem arises under pooled PT, although WPT can give substantial influence to a small local comparison group. With all precincts supported, define the effective number of veteran-weighted precincts as $n_{eff}=1/\sum_i w_{Vi}^2$. The two-sided E29 radius is $\sqrt{8\log(2/\alpha)/n_{eff}}$. At 95% it is 0.992 for 30 equal veteran weights. The reported 1.09 comes from a simulation with unequal veteran weights. Both illustrate the substantive problem. The simulation uses independent binomial draws, so it illustrates the two allowances under one data-generating process. It does not measure their performance under precinct-wide shocks.

**Answer on empirical Bernstein.** There is a finite-sample way to use the count rows. Let $K\ge2$ be the number of independent precincts and $Z_i$ the realised contribution of precinct $i$ to $D_W$, including its weight. Write $Y_i=KZ_i$, so the sample mean of the $Y_i$ is $D_W$. Choose a common deterministic range of width $R$ for all $Y_i$ before looking at the ballots. For example $R=4K\max_i w_{Vi}$ works. Let $s^2(y)$ be the usual sample variance with divisor $K-1$. For each feasible realised table $n$, define
$$
r(n)=\sqrt{2s^2(Y(n))\log(4/\alpha)/K}
       +\frac{7R\log(4/\alpha)}{3(K-1)}.
$$
Applying Maurer and Pontil's Theorem 11 to the true latent vector, after rescaling to the unit interval, and to its reflection gives two-sided coverage by a union bound. The theorem permits independent, nonidentically distributed precinct terms. Project over feasible realised tables $F$ and report
$$
[\inf_{n\in F}(D_W(n)-r(n))-\pi_U,\quad
  \sup_{n\in F}(D_W(n)+r(n))]\ \cap\ [-1,1].
$$
On the concentration event and the event that the true table lies in $F$, this covers the expected ATT by E29(a). A computationally simpler conservative version uses $[L,U]$ and an upper bound on $\sup_{n\in F}s^2(Y(n))$. For example, bound each $Y_i$ by $[a_i,b_i]$ using the observed counts. Since $\sum_i(Y_i-\bar Y)^2\le\sum_iY_i^2$, $\sum_i\max(a_i^2,b_i^2)/(K-1)$ is safe. These are bounds on the realised sample variance entering a concentration theorem. Maximising variance over allocations of just the observed returns would not by itself bound the sampling variance over unseen elections. With belief rows use the coverage penalty in E29(c). To choose the tighter of several intervals after seeing the returns, construct them with joint coverage or use a prespecified choice.

**Practical limitation.** For 30 equal weights the additive term alone is 1.41 at 95%, even with zero empirical variance, for the displayed global range $R=4$. This particular empirical-Bernstein theorem cannot improve the 0.992 Hoeffding allowance in that example. More precincts and informative count rows can change the comparison. Report $n_{eff}$, $\pi_U$, local comparison counts, and the independence unit. A proposed intraclass-correlation cap should bound the variance of the two-group, two-election contrast, including cross-group and cross-time covariances or conservative bounds on them. A single-election within-group correlation cap alone does not specify that variance. Treat such a cap as sensitivity unless supported by data. A randomisation analysis needs a credible assignment or sampling mechanism for this electorate. The realised count bound is a useful descriptive output in its own right.

[Maurer and Pontil, Theorem 11](https://arxiv.org/pdf/0907.3740) · [Hoeffding (1963)](https://doi.org/10.1080/01621459.1963.10500830)

## 2. A confidence projection using observable returns

E28(c) correctly distinguishes population identification from sampling inference. One return per precinct-election cannot recover a separate unrestricted ballot law for that precinct-election. There is a middle option that preserves the fixed electorate and WPT in every supported precinct. It does not require a common ballot law for precincts in the same composition band.

For a concrete finite-sample construction, condition on the fixed design, assume independent precincts, and preselect $J$ scalar moment rows and their design weights. For precinct $i$, election $t$, let $T_{it}$ be the observed yes rate and $L_{it},U_{it}$ the observed Duncan–Davis bounds on the veteran rate. They are functions of the published returns. Under the true law, $E L_{it}\le m_{iVt}\le E U_{it}$. The aggregate expected rate is $p_{it}=x_i m_{iVt}+(1-x_i)m_{iNt}$. For each of $H_{it}=L_{it}, U_{it}, T_{it}$ choose nonnegative weights $h_i$ with sum 1, fixed from the design. The observable weighted sum estimates $\sum_i h_i E H_{it}$ within
$$
c_h=\sqrt{\tfrac12\log(2J/\alpha)\sum_i h_i^2}
$$
simultaneously over the $J$ rows with probability at least $1-\alpha$, by Hoeffding and a union bound. Dependence between elections in one precinct is allowed. In each chosen row impose the corresponding linear condition:
$$
\sum_i h_i m_{iVt}\ge\sum_i h_i L_{it}-c_h,\qquad
\sum_i h_i m_{iVt}\le\sum_i h_i U_{it}+c_h,\qquad
\left|\sum_i h_i p_{it}-\sum_i h_i T_{it}\right|\le c_h .
$$
Each display term uses the weights and allowance of its own row. Omit nonexistent group rates at pure compositions, use the observed pure-group mixture there, and carry unsupported veteran counterfactuals as free rates. Retain E28's local WPT rate constraints, factual FD-LCE if maintained, and the full-ATT objective. On the simultaneous event the true mean vector is feasible, so its projected interval covers the true expected ATT. Indeed the event contains the known-law identified set for those same maintained restrictions. No severity-of-selection or endpoint differentiability argument is needed for this conservative projection.

For a first comparison, use the entire electorate and a few predeclared veteran-share bands. The bands define which *observable moments* to average. The causal assumption still compares veterans to non-veterans within each precinct. Broad bands can discard identifying information. Very narrow bands incur large allowances. This construction establishes coverage, not convergence to the E28 sharp set or useful width with 30 precincts. Approaching sharpness would require a sufficiently rich sequence of well-estimated moments, enough comparable precincts and explicit conditions on the design or conditional means. Compare coverage and width against E29 in simulations with common precinct shocks, concentrated veteran weights, and thin comparison cells before recommending it.

[Andrews and Shi on instrumenting conditional moment inequalities](https://doi.org/10.3982/ECTA9370)

## 3. Sensitivity to violations of WPT

For the owner's question, local WPT is the natural starting assumption. A positive untreated trend gap in one precinct cannot compensate for a negative one in another under exact WPT. Aggregate relative magnitudes (RM) permits that cancellation. Both specifications can be informative, but they represent distinct historical claims.

**Correction to E29(d).** The paper defines $b_{it}$ relative to baseline, so the Rambachan–Roth-style restriction from E18 applies to $b_{it}-b_{i,t-1}$ at each post-election. The display in E29(d) uses $b_{it}$ itself. It agrees at the first post-election. At the second, a pre-period maximum of 0.10 and post gaps of 0.10 and 0.20 satisfy the consecutive-increment restriction with $\bar M=1$, but violate the displayed cap on the second level. Make the same repair after aggregation. A cap on levels is a legitimate alternative if explicitly labelled. RM with a positive allowance relaxes WPT. Limit the sums to supported precincts and retain the unrestricted contribution of the others.

**An intermediate assumption.** Prespecify a few strata by pre-treatment composition or region. Within each stratum, average the supported precincts' gap increments using fixed veteran weights, and benchmark each stratum's post increment against its own largest absolute pre increment. Keep the local baseline and counterfactual rate limits, writing the counterfactual in terms of the local violation. This permits cancellation within a stratum. At zero allowance it enforces parallel trends on stratum averages, which is weaker than exact precinct WPT. If the aim is to relax exact WPT continuously, use a local absolute dial $|b_{it}-b_{i,t-1}|\le M$ in every supported precinct. This is one LP for a fixed $M$; the strict WPT model is recovered at zero. A stratum-specific **fixed** tolerance for all its precincts works the same way. Rates inferred from pre-periods require a joint uncertainty calculation.

**A tractable exact RM case.** With one pre-election transition and two-group FD-LCE on that transition, the expected factual group-gap increment is affine in the same scalar composition. Sorting supported precincts by composition, the sign of that affine function can switch at most once. Enumerate the split and orientation, enforcing the corresponding weak signs. On each branch the absolute pre increment is linear, so the local RM limits on all post increments are linear. This needs at most $2(|\mathcal S|+1)$ LPs per underlying convex branch, including constant-sign and zero cases, versus $2^{|\mathcal S|}$ arbitrary sign choices. Keep the union of branch projections if there are holes. This reduction relies on a common affine pre-change profile. It does not extend unchanged to arbitrary local pre-trends. With known benchmarks no sign search is needed.

[Rambachan and Roth (2023)](https://doi.org/10.1093/restud/rdad018)

## 4. Small repairs that affect interpretation

* **Pooled witness in E28(c).** Step 3 constructs the untreated rate $m_{iV0}+d_{iN}$ locally. That is the WPT witness. Under PPT only the weighted mean $p^0=\sum_i w_{Vi}m_{iV0}+\sum_i w_{Ni}d_{iN}$ is constrained. Assign the untreated post veteran rate $p^0$ to every veteran-containing precinct. This meets PPT and the probability limits for every feasible $p^0$. For a check, take equal weights for both groups across two precincts, veteran baselines $(.1,.9)$ and comparison changes $(-.3,+.3)$. The local imputation gives $(-.2,1.2)$, while the pooled imputation $(.5,.5)$ is valid. Step 1 should invoke mean zero spillover only at the aggregation level assumed. The stated pooled interval survives the proof repair.
* **Support.** E28(a)'s displayed full-ATT identity is correct. Its following shorthand $\sum_i w_{Vi}d_{iN}$ for the WPT-imputed trend requires full support. Otherwise use the supported sum plus $\sum_{i\in\mathcal U}w_{Vi}(\theta_i-m_{iV0})$. Omit zero-veteran rows from objectives and avoid evaluating nonexistent group rates.
* **Spillover.** For sensitivity let $S_i=m_{iNt}-m^{(0)}_{iNt}$ on supported precincts and impose $|S_i|\le\Lambda$, $m_{iNt}-S_i\in[0,1]$, and $m_{iV0}+d_{iN}-S_i\in[0,1]$. Optimise jointly. For a fixed factual fit the largest absolute additive contribution is $\pi_S\Lambda$ before the rate limits. Relaxing spillover can also admit new factual fits. Simple dilation of the sharp zero-spillover interval has no general sharpness guarantee.
* **E30.** The sign holds for $\epsilon<1/600$. At $\epsilon=1/600$, zero first becomes feasible, so strict positivity is lost. The identified sets still nest. In general $\epsilon=1/2$ already removes the baseline constraint because the constant line $1/2$ fits every rate to within $1/2$.

## Verification scope

I reran the verifier at the pinned theory commit. All 98 checks passed, including both five-precinct rows, all six three-precinct LCE fits and the curvature breakdown. In ten additional six-precinct LP comparisons with one affine pre-increment, the split-and-orientation enumeration and exhaustive sign enumeration agreed on both extrema. The finite-sample moment projection and empirical-Bernstein construction above are proposed derivations. Their application-specific coverage and width have not been simulation-tested here.
