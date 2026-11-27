# Reply to the independent review of PHASE2e: new sections and open problems

**For:** the author of `PHASE2e_independent_review.md`.
**Revised theory:** [`PHASE2e_did_general.md` at commit `44cd35f`](https://github.com/mdweaver/eiid/blob/44cd35f5fcbf55f944aff700dbb44d9aa0e8fee4/docs/extension/PHASE2e_did_general.md) (branch `claude/did-theory-r3`). Line numbers below refer to that file.
**Checks:** [`verify_phase2e.py`](https://github.com/mdweaver/eiid/blob/44cd35f5fcbf55f944aff700dbb44d9aa0e8fee4/docs/extension/verify_phase2e.py) passes 98/98. Your two scripts are ported into it (sections IR and PT) and both reproduce unchanged: [`phase2e_independent_checks.py`](https://github.com/mdweaver/eiid/blob/44cd35f5fcbf55f944aff700dbb44d9aa0e8fee4/docs/extension/reviews/phase2e_independent_checks.py) and [`phase2e_independent_checks_pt.py`](https://github.com/mdweaver/eiid/blob/44cd35f5fcbf55f944aff700dbb44d9aa0e8fee4/docs/extension/reviews/phase2e_independent_checks_pt.py).

Thank you for the review. This note does three things:
1. It points you to the sections that implement your follow-up requests: §5.8, §5.9 and the revised §9.
2. It lists how each numbered finding was handled.
3. It sets out the problems we could not solve, with specific questions.

We would value your view on whether the new sections implement what you proposed, and any route you see for the open problems. The revisions have not yet had a further proof review on our side.

## 1. Owner's decisions that shaped the revision

* **Comparison choice.** Within-precinct parallel trends in expectation (WPT) is a separate option beside electorate-wide parallel trends (PPT). **There is no default**: the user must choose.
* **Unsupported precincts.** The target stays the **full** veteran ATT. Precincts with fewer than $a_{\min}$ non-veterans are bounded with no assumption. The veteran share living there, $\pi_U$, is always reported.
* **Joint level/change LCE.** Adopted as §5.9 (E30).
* **Age adjustments.** Deferred to a later build step on covariates. That step now covers both ordinary EI and DiD, including both of your age options. Gap item G8 points to it.

## 2. The new and revised sections

### §5.8 Parallel trends within each precinct, in expectation (E28–E29), lines 444–504

* **Expected rates (l. 450–452).** $m_{igt}=E[B_{igt}\mid\mathcal X]$ under a probability law for ballots in each regime, with the design held fixed. A random draw of persons from a fixed schedule of potential ballots is not used; that is E16's realised version.
* **Estimand (l. 454).** $\mathrm{ATT}^{\rm mean}_t=\sum_iw_{Vi}(m_{iVt}-m^{(0)}_{iVt})=E[\mathrm{ATT}_t\mid\mathcal X]$, the expected effect on the same veterans.
* **Assumptions (l. 456–460).** WPT or PPT; mean spillover zero, or bounded by $\Lambda$; support $\mathcal S=\{i:a_{iN}\ge a_{\min}\}$.
* **Proposition E28 (l. 462–476):**
  * (a) the identity under WPT, including the unsupported terms $\theta_i$, and under PPT;
  * (b) the rate limits $m_{iV0}+d_{iN}\in[0,1]$ in every supported precinct;
  * (c) the sharp set for a **known** law of aggregate totals, as an Aumann-expectation LP;
  * (d) the case $a_{iN}=0$: the interval is shifted and extended downward by $\pi_U$;
  * (e) the FD-LCE projection: $\sum_{\mathcal S}w_{Vi}q^\Delta(X_i)+\pi_S\Delta\kappa$ under WPT, against $I+s\Delta\kappa$ under PPT.
* **Your averaging construction (l. 478).** Your script averages per-state Duncan–Davis ends over equally likely totals. We state this as **one** admissible construction of E28(c) for a known law, not as the general definition.
* **Worked comparison (l. 480–487).** Your five-precinct table reproduces:

  | Margins | WPT | PPT |
  |---|---|---|
  | Rare | $[+10.39,+12.25]$ | $[-3.96,-2.22]$ |
  | Ordinary | $[-67.33,+43.78]$ | $[-80.56,+52.13]$ |

  We cite, but did not recompute, your true-ATT column and generating values (+11.45, +6.00, 19.09, 30.55, 16.00). The four-precinct data give $D_W\in[-0.9387,1.4472]$, against the pooled DiD $[-0.4075,0.7925]$. Precinct 4 carries 47.5% of the comparison trend under WPT and 2.5% under PPT.
* **Proposition E29 (l. 489–504):**
  * the M-Count target $D_W$, including the unsupported term;
  * (a) $E[D_W\mid\mathcal X]-\mathrm{ATT}^{\rm mean}\in[0,\pi_U]$;
  * (b) the allowance $r_\alpha$: Hoeffding with $\rho_i=4w_{Vi}$ under arbitrary within-precinct dependence, or a person-level version with $1/a_{iN}$ terms;
  * (c) coverage of $[L-\pi_U-r_\alpha,U+r_\alpha]$, with identity rows and belief rows separated;
  * (d) the dial and RM on expected trends, with the two RM versions (within precinct, or after aggregation);
  * (e) E25 does not transfer to this target.

### §5.9 Baseline, post-election and change LCE in one calculation (E30), lines 506–515

* **(a) Any two imply the third.** Any two of baseline, post and change LCE imply the third when the compositions are the same at both elections. This can fail when compositions change.
* **(b) Baseline-curvature allowance $\epsilon$.** The set is monotone in $\epsilon$; $\epsilon=0$ gives the two-election LCE set and $\epsilon\ge1$ the change-only set. The breakdown value $\epsilon^*$ is one LP.
* **(c) Your three-precinct check reproduces exactly:**
  * no LCE: $[-58.00,27.33]$;
  * change LCE only: $[-49.33,7.33]$;
  * any two of the three, or all three: $[2.67,7.33]$, exactly $[2/75,11/150]$.

  The sign survives baseline curvature up to $\epsilon^*=1/600$.

### §9, plain language, lines 629–692

* §9.1 adds the fifth choice: whose trend stands in for the veterans. No anticipation is explained as "treatment did not affect votes before it began".
* §9.2:
  * the ladder gains a $D_W$ row;
  * the RM row now says the outcome depends on these data;
  * the RM-with-CR row states that CR must be on both changes, and gives $[-0.4155,0.7325]$ with CR on the pre-war change only.
* §9.3:
  * item 6 is rewritten around "admissible pre-trend against post signal", with your (99,1)/(1,99) example $[.12,.28]$ and the growing allowance at later horizons;
  * item 9 covers the comparison choice;
  * item 10 notes that there is no one-step integer guarantee.
* §9.4 (new) compares the two comparison choices in a table, with your five-precinct example and the E30 example.

## 3. How each finding was handled

Every fix is logged in §12 under "Independent review" (rows IR1–IR12 and IR-F, IR-L, IR-A, IR-P, IR-S).

| Your finding | Change |
|---|---|
| 1. RM horizons accumulate | E18(b) is restricted to the first post-period with $\Lambda=0$. The marginal outer set at horizon $h$ is $D_h\oplus[-h\bar M\sup A-\Lambda,h\bar M\sup A+\Lambda]$. Joint statements keep path constraints. Branch counts are stated: $2LK_F$ for RM and $K_F$ for SD. |
| 2. Finite-sample FD | E23(b) doubles each summand's range in P3.8. E12's allowance is a separate event. Your .884682/.997423 check is included. |
| 3. Regression CLT and covariance | E23(a) and E25(H1) now state the design-sequence (Gram limit), Lindeberg/leverage and clustered-covariance conditions. "Any regime" is removed. E23(c) requires a correctly specified level mean, and grid rows are outer. The E25 zero-box clause applies only with the restriction at the mean level; your ±.10 example is checked. |
| 4. Integer witness | E16(b) adds $(B_s-S)A_s\in\mathbb Z$, greedy allocation per group, one designated outcome, and "net" spillover. |
| 5. No one-step integer guarantee | Remark 7.1 now has your $\{0\}$ against $[-1,1]$ example. |
| 6. Star gap | E27 records the strict gap ($3/2$ against $1$). Cycles are TU, with a short proof. |
| 7. E7(d) loading | The weighted loading replaces the common-composition claim (your .0075 example). The tree-rank proof is added. |
| 8. E8 across shared elections | New E8(e): pairwise stacking is outer, and the joint sharp set comes from Lemma G′ on levels. "Strict is the rule" is labelled numerical experience. |
| 9. Causal clarifications | E16(d) requires one shared table across post-periods. E15's proof is rewritten for disconnected $F$. E17(c) gets an optional sharp breakdown LP. E21 uses $S=\min(\Lambda,B_s)$, with one hypograph LP per branch and a shared table for several post-periods. E22's scope is stated. Def. 5.1 is labelled a reporting convention. |
| 10. When RM helps | New E19(c4) with your example. (a)'s "iff" now requires $\bar M>0$; (c1) is weak nesting; (c3) is fixed; (d) now says plug-ins "can" under-cover. §9.3 item 6 is rewritten. |
| 11. Uncertain margins, gap register | G7(a)–(b) include the truth of the rows; the target is written $\varphi_a$; the union is not one LP; G7(c) states its domain. G1–G9 are revised, including the G6 witness (90,10)/(10,90) and a joint denominator scan for G5. |
| 12. Small repairs | E3(c) for consecutive differences; E4, E6, E10, E12; E24 holds for a fixed true branch. The RM set with a hole from one polytope, $[-.90,-.70]\cup[-.30,0]$, is a check, and branches are merged, never enveloped. |

## 4. Open problems: where we would most value your view

### Problem 1. The coverage allowance for the expected ATT can be uninformative (E29(b))

**The issue.**
* The count layer gives $[L,U]$ for the realised target $D_W$ with no sampling error.
* Moving to the expected ATT adds an allowance $r_\alpha$ for $\lvert D_W-E[D_W\mid\mathcal X]\rvert$.
* With independence across precincts and arbitrary dependence within a precinct, across groups and elections, each precinct is one bounded draw. Hoeffding with ranges $\rho_i=4w_{Vi}$ gives $r_\alpha\approx1.09$ with 30 equal precincts (check E29 [F]). That is wider than the logical range.
* The normal version needs $\max_i\rho_i^2/\sum_j\rho_j^2\to0$, which fails exactly when veterans are concentrated in a few precincts.
* Independence of persons' ballot pairs gives $r_\alpha\approx0.12$ in the same check, but adds $1/a_{iN}$ terms that are large where comparison voters are few.

**What we are considering.**
* (a) A variance-based allowance (Bernstein or empirical Bernstein) in place of the range-based one. Within-precinct variances of rate changes are probably far below the worst-case range.
* (b) A bounded within-precinct dependence assumption, for example an intraclass-correlation bound or a design-effect cap, interpolating between person-level independence and arbitrary dependence.
* (c) Reporting the realised version ($D_W$, or E16) as the primary result, with the expected effect secondary.

**Questions.**
* Is (a) sound here, given that the variance of $D_W$'s precinct terms is not identified from aggregate returns? Would a worst-case variance under the count rows work?
* Is there a standard dependence assumption for (b) that is defensible for voting data?
* Is there a better route we have missed, for example randomisation-style or design-based arguments?

### Problem 2. No coverage result for the sharp set when the ballot law is estimated (E28(c))

**The issue.** E28(c) is sharp when the law of each precinct's aggregate totals is known, as in your examples with two equally likely totals. In an application we see one realised total per precinct per election, so the law, and hence the per-state Duncan–Davis averaging, must be estimated. We have no result that an estimated-law version of E28(c) covers $\mathrm{ATT}^{\rm mean}$. The only route with a coverage proof is E29 (count bounds plus $r_\alpha$), which runs into Problem 1.

**Questions.**
* What would you take as the estimable object here? For example, $E[L(T;x)\mid X=x]$ and $E[U(T;x)\mid X=x]$ modelled across precincts with similar compositions, as in our open inference problem for the sharp FD set (E23(c)).
* Is there a construction that keeps the Aumann-expectation sharpness and has a coverage theorem? Or should the sharp set be presented only as an identification result, with E29 the inferential statement?

### Problem 3. Relative magnitudes for expected trends: two versions, no relation, no specific inference (E29(d))

**The issue.** Under WPT, the RM violation can be measured in two ways:
* **within each precinct:** $\lvert E[b_{i,t}]\rvert\le\bar M\max_s\lvert E[D_{i,s}]-E[D_{i,s-1}]\rvert$ for each supported $i$, with $(2L)^{\lvert\mathcal S\rvert}$ branches (a MILP with one binary per precinct and branch);
* **after aggregation:** $\lvert\sum_iw_{Vi}E[b_{i,t}]\rvert\le\bar M\max_s\lvert\sum_iw_{Vi}(E[D_{i,s}]-E[D_{i,s-1}])\rvert$, with $2L$ branches.

Neither implies the other:
* pre-period gaps of $+.1$ and $-.1$ cancel in the aggregate, so the aggregated version forces zero bias while the precinct version allows $.1$ each;
* gaps of $0$ and $.2$ make the precinct version force zero bias in the first precinct, while the aggregated version allows $.1$ anywhere.

The expected pre-trends are identified only through a model of expected rates, so inference needs a joint region for them. We point only to E25's mean-level version under (H1)–(H2), with no specific theorem.

**Questions.**
* Which version matches the substantive within-precinct comparison you had in mind? Is there a natural intermediate version?
* With $(2L)^{\lvert\mathcal S\rvert}$ branches the within-precinct version may be computationally heavy and statistically demanding. Do you see a tractable formulation, or a reason to prefer the aggregated version for inference?

### Other items still open (for completeness)

* Inference for the sharp FD set (E8, E23(c)).
* CLR versions, which remain algorithmic.
* Integer theory for several columns and mixtures beyond E27.
* Grouped rungs with model rows on changes (R16).
* Uncertain margins, R25(b)–(c): deferred to a later module.

## 5. What would help most

1. Confirmation that §5.8 and §5.9 implement your follow-up as intended, or corrections.
2. Your view on Problems 1–3, with any references.
3. Anything in the §12 "Independent review" fixes that does not match what you asked.
