# Review: `eco`'s Dirichlet-process Gibbs sampler versus the `eiDPM` blocked Gibbs sampler


*The source code reviewed is `eco` 4.0-6:*
* *`src/gibbsDP.c`, the CAR DP sampler behind `ecoNP`;*
* *the shared routines `rand.c` (`rWish`, `dMVT`), `sample.c` (`rMH`, `rGrid`) and `bayes.c` (`NIWupdate`).*

*It is compared with `package/eiDPM/src/gibbs.cpp`.*

*Status labels follow `PHASE1_proposal.md`:*
* ***verified*** = checked in the source and by a run;
* ***proved*** = short argument given here;
* ***conjecture***.

*Theorem numbers marked † are cited from memory and must be checked.*

## 1. Step by step

| Step | `eco` (`gibbsDP.c`) | `eiDPM` (`gibbs.cpp`) | Same target? |
|---|---|---|---|
| Data and boundary rule | Units with $Y\in\{0,1\}$ fixed at $10^{-4}$ or $1-10^{-4}$; $X\in\{0,1\}$ units partly observed | Same rule (`dpmData`) | yes |
| Stated prior | $\Sigma\sim IW(\nu_0,S_0)$, $\mu\mid\Sigma\sim N(\mu_0,\Sigma/\tau_0)$, $\alpha\sim\Gamma(a_0,b_0)$ | Same | yes, as stated (but see B1) |
| DP representation | Pólya urn: Escobar–West / Neal (2000) Algorithm 1, over **units**, so $O(n^2)$ per sweep, plus Bush–MacEachern "remixing" | Truncated stick-breaking, blocked Gibbs (Ishwaran & James 2001), $O(nK)$ per sweep | same in the limit $K\to\infty$ (§3a) |
| Latent rate $W_i$ on its line | One independence Metropolis–Hastings step, uniform proposal on $[L_i,U_i]$; target $\phi(\operatorname{logit}w)/\prod w_j(1-w_j)$ | One slice-sampling step on $(L_i,U_i)$, same target | yes (both leave the conditional invariant) |
| Assignment weights | Existing value: $\propto N(W_i^*\mid\theta_j)$ for each other unit $j$. New value: $\propto\alpha\,t_{\nu_0-p+1}(W_i^*\mid\mu_0,S_0(1+\tau_0)/(\tau_0(\nu_0-p+1)))$, the NIW($\nu_0$) prior predictive | $\propto\pi_k N(W_i^*\mid\theta_k)$ | see B1 |
| Cluster parameters | `NIWupdate` (conjugate) | Conjugate NIW update | formulas identical; the Wishart draw differs (B1) |
| $\alpha$ | Escobar–West (1995) auxiliary variable $\eta$ | Conjugate Gamma given the sticks | see B2 |
| Initial state | Every unit in its own cluster | All units in one component, $\alpha=1$ | irrelevant at stationarity |

## 2. Two departures in `eco` from its stated model

### B1. `rWish` draws Wishart($\mathrm{df}-1$), not Wishart($\mathrm{df}$) (proved; verified)

`rand.c`, `rWish(Sample, S, df, size)`:
* $V_i\sim\chi^2(\mathrm{df}-i-1)$ for $i=0,\dots,p-1$, and $N_{ij}\sim N(0,1)$ for $i<j$.
* It returns $C\,T^\top T\,C^\top$, with $T$ upper triangular, $T_{ii}=\sqrt{V_i}$, $T_{ij}=N_{ij}$, and $C=\operatorname{chol}(S)$.

*Proof that this is Wishart($\mathrm{df}-1$, $S$).*
1. Bartlett: if $T$ is upper triangular with independent $T_{ii}^2\sim\chi^2(m-i)$ ($i=0,\dots,p-1$) and $T_{ij}\sim N(0,1)$, then $T^\top T\sim W_p(m,I)$. [Bartlett decomposition; Anderson 2003 Thm 7.2.1†; Muirhead 1982 Thm 3.2.14†]
2. `eco` uses $\chi^2(\mathrm{df}-i-1)=\chi^2((\mathrm{df}-1)-i)$, so $T^\top T\sim W_p(\mathrm{df}-1,I)$. [1 with $m=\mathrm{df}-1$]
3. $C\,W_p(m,I)\,C^\top\sim W_p(m,CC^\top)=W_p(m,S)$. [linear transformation of a Wishart]
4. Check: $E[(T^\top T)_{jj}]=(\mathrm{df}-j-1)+j=\mathrm{df}-1$. A transcription of `rWish` in R gives mean$/S=\mathrm{df}-1$ (`diagnostics/nu_test_output.txt`). ∎

**Consequences.**
* `NIWupdate` calls `rWish(..., nu0 + n)`, so every conjugate draw of a covariance comes from $IW(\nu_0+n-1,S_n)$ instead of $IW(\nu_0+n,S_n)$. The initial draws in `gibbsDP.c` use $\nu_0-1$ too.
* **Parametric `eco` (proved).**
  * Prior and updates are off by one consistently, so the sampler is a valid Gibbs sampler for the stated model with $\nu_0$ replaced by $\nu_0-1$ (3 instead of 4). [conjugacy: prior df + $n$ = posterior df]
  * With $n$ in the hundreds the effect is negligible. This is consistent with the benchmark, where parametric `eco` and `eiDPM` agree at unit level (correlation 0.997–1.000).
* **`ecoNP` (conjecture, with a mechanism).**
  * The new-cluster weight uses the predictive $t$ of the NIW($\nu_0$) prior, while the cluster parameters are drawn as if the prior were NIW($\nu_0-1$).
  * The conditionals are therefore not those of any single joint law, and the chain's stationary distribution (if it exists) is not the posterior of either prior.
  * Direction: the drawn covariances are inflated most for small clusters. For $p=2$, the mean of $IW(\nu,S)$ is $S/(\nu-p-1)$; a singleton cluster has $\nu_0+1=5$ intended against 4 drawn, so its mean covariance doubles.
  * Wider clusters absorb more units, so we expect fewer clusters and a smaller $\alpha$ than the stated model implies.
  * That is the pattern on `reg`: `eco` gives 4 clusters and $\alpha\approx0.7$; `eiDPM` gives 5–6 clusters and $\alpha\approx1$.
  * Test: run `eiDPM` with $\nu_0=3$ (§4).
* For NCAR ($p=3$, $\nu_0=4$), the effective $\nu_0-1=3$ is still proper ($>p-1$), but the inverse-Wishart then has no mean. The stated prior's mean is already infinite for $p=3$.

**Bug or convention? (reasoned; see also the paper.)**
* The current upstream source (GitHub `kosukeimai/eco`, `src/rand.c:129`) is unchanged.
* For $p=1$, `rWish` returns $S\chi^2(\mathrm{df}-1)$, whereas $W_1(\mathrm{df},S)=S\chi^2_{\mathrm{df}}$. No standard parametrisation shifts df by a constant $-1$ in every dimension. Dawid's shifts by $p-1$, and in the other direction.
* The `eco` documentation calls `nu0` "the prior degrees of freedom".
* ILS App. B states the posterior $IW(\nu_0+n,\cdot)$ and the new-cluster weight as a bivariate $t$ with $\nu_0-1$ df and scale $S_0(1+\tau_0^2)/(\tau_0^2(\nu_0-1))$. That is the standard predictive for an $IW(\nu_0)$ prior when $p=2$, and the code implements it (`dMVT(..., nu0 - n_dim + 1)`).
* One sentence in ILS §3.2 ("for the inverse-Wishart distribution to be proper, $\nu_0$ needs to be greater than 3" for $p=3$) fits a df$-1$ convention, since the standard condition is $\nu_0>2$. It contradicts App. B.
* **Decisive point:** in `ecoNP` the new-cluster weight assumes standard df $\nu_0$, while the covariance draws use $\nu_0-1$.
  * If $\nu_0$ is the standard df, the draws are wrong.
  * If $\nu_0$ is df$+1$, the predictive should have $\nu_0-2$ df.
  * So no convention makes the DP sampler's conditionals come from one model.
* Separately, ILS write $\mu\mid\Sigma\sim N(\mu_0,\Sigma/\tau_0^2)$, while the code uses $\Sigma/\tau_0$ with a matching predictive scale $(1+\tau_0)/\tau_0$. `eiDPM` follows the code.
* ILS eq. (7) gives yet another predictive ($\nu_0$ df, scale over $1+\nu_0$), inconsistent with App. B.

**The package paper (Imai, Lu & Strauss 2011, J. Stat. Softw. 42(5), §§2.2–2.3; checked):**
* It states the same prior for the parametric CAR and CCAR models and for the DP base measure: $\mu\mid\Sigma\sim N(\mu_0,\Sigma/\tau_0^2)$, $\Sigma\sim\text{InvWish}(\nu_0,S_0^{-1})$, with $\nu_0$ "the prior degrees of freedom parameter".
* It does not describe the sampler (it refers to ILS 2008) and gives no shifted-df convention.
* Its version history (to 3.1-4, 2009) lists no change to the prior or the Wishart routine.
* It says the defaults ($\tau_0=2$, $\nu_0=4$, $S_0=10I$) give "an approximately uniform prior distribution of $W_1$ and $W_2$". Prior-predictive deciles of $W_1$ (uniform = 0.10):
  * stated $IW(4,10I)$ with $\Sigma/\tau_0$: 0.20 in each outer decile, about 0.07 in the middle;
  * as drawn, $IW(3,10I)$: 0.25 outer, about 0.05 middle;
  * paper notation $\Sigma/\tau_0^2$: 0.18 outer.

  None is uniform. The stated prior is closer than what `eco` draws, so the remark gives no support to an intended df$-1$.

**Patched `eco` (verified; `diagnostics/` scripts to be committed).** Changing `df-i-1` to `df-i` in a local copy, on `reg`, CAR, DP, two runs each (`alpha` = posterior median; RMSD is of unit posterior means):

| | $\alpha$ | clusters | unit cor. with Gibbs | unit RMSD vs Gibbs |
|---|---|---|---|---|
| `eco` original | 0.68, 0.73 | 4.2, 4.4 | 0.905 | 0.068 |
| `eco` patched | 0.72, 0.89 | 4.4, 5.2 | 0.932 | 0.058 |
| Gibbs | about 0.9–1.0 | 5.2–5.5 | — | — |

The fix moves `eco` toward `eiDPM`, and one patched run matches it, but it does not close the gap. The rest is consistent with slow, seed-dependent mixing over the number of clusters in both samplers (§3b).

The $\nu_0=3$ mimic in `eiDPM` (not an exact mimic, since a blocked sampler has no new-cluster weight) did not move it toward `eco`: unit correlation 0.90 either way.

### B2. The $\alpha$ update uses the odds as a probability (verified; negligible)

Escobar & West (1995, eq. 13†):
* draw $\eta\sim\text{Beta}(\alpha+1,n)$;
* then draw $\alpha$ from the mixture $\pi_\eta\,\Gamma(a+k,\ b-\log\eta)+(1-\pi_\eta)\,\Gamma(a+k-1,\ b-\log\eta)$, with
$$\frac{\pi_\eta}{1-\pi_\eta}=\frac{a+k-1}{n(b-\log\eta)}.$$

`gibbsDP.c` computes this odds ratio $r$ and selects the first component when `unif_rand() < r`. So it uses $\min(r,1)$ as the probability instead of $r/(1+r)$.

On `reg` ($n=275$, $k\approx4$, $\alpha\approx0.7$): $-\log\eta\approx5$, so $r\approx0.003$, and the error in the mixing probability is about $r^2\approx10^{-5}$. **Negligible.** It matters only when $n(b-\log\eta)$ is small relative to $k$.

## 3. What `eiDPM` changed, and what that implies

**(a) Truncation at $K$ (a change of model).**
* The truncated stick-breaking prior approximates the DP. Ishwaran & James (2001, Thm 2†) bound the total-variation error in the marginal law of the data by about $4n\exp(-(K-1)/\alpha)$.
* At $\alpha\approx1$ this is negligible. On `reg`, $K=30$ and $K=100$ gave the same answer (`diagnostics/reg_dp_modes_output.txt`).
* But the prior $\alpha\sim\Gamma(1,0.1)$ has mean 10, so when data favour a large $\alpha$, $K=30$–$50$ is not enough. SBC design B shows this: the number of clusters is badly miscalibrated ($p\approx6\times10^{-32}$), while unit rates and means are calibrated.

**(b) Blocked Gibbs over ordered stick-breaking weights (a mixing issue, not a change of target).**
* The weights are stochastically ordered by label. The sampler has no move that permutes labels, so it can move slowly between orderings.
* Through the dependence of $\alpha$ on $\sum_k\log(1-v_k)$, this slows the mixing of $\alpha$ and of the number of clusters.
* This is documented: Papaspiliopoulos & Roberts (2008, Biometrika) and Hastie, Liverani & Richardson (2015, Statistics and Computing) propose label-switching moves for exactly this.
* Evidence here:
  * with the true rates treated as observed, one seed stayed at $\alpha\approx1.9$ while another found `eco`'s $\alpha\approx0.9$;
  * SBC design A, which simulates from the sampler's own target, shows a mild excess of true $\alpha$ above all posterior draws ($p=0.017$).
* **Remedy:** add the label-switching moves, or use a collapsed sampler (Neal 2000, Algorithm 3: conjugate, label-free, exact DP, $O(nk)$ per sweep).

**(c) Slice sampling versus independence MH for $W_i$.** Same target. Slice sampling adapts to peaked conditionals; independence MH with a uniform proposal can have low acceptance there. This affects mixing only.

**(d) Bounds shrunk by $10^{-8}$ of their width, and the initial state.** No practical effect.

## 4. Tests

* **SBC** (`benchmarks/eco_2x2/sbc.R`; results in `results/sbc_summary.csv`):

  | Design | Unit rates and means | $\alpha$ / cluster count |
  |---|---|---|
  | A: exact target, $K=30$ | calibrated ($p=0.36$–$0.90$) | mildly off ($p=0.017$ / $0.08$; consistent with 3b) |
  | B: DP data, $K=50$ | calibrated ($p=0.08$–$0.81$) | clusters badly miscalibrated (3a) |
  | C: `ecoNP`, 400 of the B simulations | calibrated ($p=0.12$–$0.95$) | **both miscalibrated, systematically** ($p=1.2\times10^{-4}$ / $4\times10^{-5}$; mean rank 0.574 / 0.578: too few clusters and too small an $\alpha$, the direction B1 predicts) |

* **The $\nu_0$ test and the patched-`eco` test:** see §2, B1.

**Same 400 simulations, Gibbs (design B, $K=50$) vs `ecoNP` (design C)** (`alpha` and `clusters` rows; mean rank 0.5 = unbiased):

| | Gibbs | `ecoNP` |
|---|---|---|
| unit $W_1$, mean $W_1$ | calibrated ($p=0.94$, $0.38$) | calibrated ($p=0.12$, $0.95$) |
| $\alpha$ | $p=0.075$, mean rank 0.487 | $p=1.2\times10^{-4}$, mean rank 0.574 (underestimates) |
| clusters | $p=1.8\times10^{-10}$: a tail failure, with 15% of ranks extreme (true count above all draws); mean rank 0.524 | $p=4\times10^{-5}$: a systematic shift, mean rank 0.578 (too few clusters); 1% extreme |

## 5. Conclusions

1. **Both samplers are calibrated for the quantities reported in EI:** unit rates and in-sample aggregates (SBC A, B, C).
2. **`ecoNP` systematically underestimates $\alpha$ and the number of clusters.** That is the direction predicted by the Wishart off-by-one (B1), and patching it moves `ecoNP` toward `eiDPM` on `reg`.
3. **`eiDPM` is calibrated for $\alpha$ under its own truncated model (A, $p=0.017$, mild)** and roughly so against DP data (B). It has a tail failure for the cluster count when the true $\alpha$ is large: truncation at $K$ (3a) plus slow label mixing (3b).
   * Remedies: a larger or adaptive $K$; label-switching moves (Papaspiliopoulos & Roberts 2008; Hastie, Liverani & Richardson 2015); or an exact collapsed sampler (Neal 2000, Alg. 3).
4. **Effect on the `reg` discrepancy.** On data like `reg`, where the DP clustering drives unit-level estimates, the two samplers' different $\alpha$ and cluster counts translate into the unit-level differences seen there (RMSD about 0.06). On the SBC evidence, `eco`'s $\alpha$ and cluster count are the biased ones. SBC averages over the prior, though, so it does not *prove* which answer is right on `reg` specifically.
   * **Revised (2026-11-02, `benchmarks/eco_2x2/eco_fixed/`):** a corrected re-implementation of `eco` (B1 and B2 fixed) reproduces most of the `reg` discrepancy at `eco`'s run length when it keeps `eco`'s uniform-proposal MH step for the rates, and removes it when that step is replaced by slice sampling. So the `reg` difference is mainly slow mixing, with B1 the smaller part. The SBC bias of `ecoNP` in item 2 is likewise consistent with B1 but may be partly mixing.

## 6. The exact sampler checked against Neal (2000)

Source: Neal, "Markov Chain Sampling Methods for Dirichlet Process Mixture Models", Technical Report 9815, University of Toronto, 1998 (published in *JCGS* 9, 2000), §3. Neal's Algorithm 1 is the Escobar (1994) / Escobar–West (1995) per-unit sampler. Algorithm 2 (cluster labels plus cluster parameters) is attributed to Bush & MacEachern (1996) and West, Müller & Escobar (1994); `eco` is Algorithm 1 followed by an Algorithm-2-style redraw of each cluster's parameters (the "remixing"). Algorithm 3, attributed to MacEachern (1994) and Neal (1992), keeps only the labels $c_i$ and updates each by Gibbs sampling from his eq. (12):

$$P(c_i=c\mid c_{-i},y_i)\propto \frac{n_{-i,c}}{n-1+\alpha}\int F(y_i,\phi)\,dH_{-i,c}(\phi),\qquad P(c_i\text{ new}\mid c_{-i},y_i)\propto\frac{\alpha}{n-1+\alpha}\int F(y_i,\phi)\,dG_0(\phi).$$

**Code check (`src/collapsed.cpp`), item by item:**

| Item | Neal / standard result | Code | Status |
|---|---|---|---|
| Assignment weights | $n_{-i,c}\times$ predictive of $c$ without $i$; $\alpha\times$ prior predictive (eq. 12; the common denominator cancels) | `log(m_c) + lt(x, cl[c])`; `log(alpha) + lt(x, prior_c)`, after removing $i$ | matches |
| Predictive under the NIW base | $t_{\nu_m-p+1}\bigl(\mu_m,\;S_m(\tau_m+1)/(\tau_m(\nu_m-p+1))\bigr)$, $S_m=S_0+\sum(x-\bar x)(x-\bar x)^\top+\frac{\tau_0 m}{\tau_m}(\bar x-\mu_0)(\bar x-\mu_0)^\top$ | `cache()` | matches |
| Empty cluster | removed from the state (Neal: "remove $\phi_{c_i}$") | dropped and relabelled | matches |
| Latent rates (not in Neal: data augmentation) | $W_i\mid c,W_{-i}\propto t_{c_i,-i}(W_i^*)\,/\prod w(1-w)$ on the line | slice sampler on $W_1$ with the cluster's predictive without $i$; the prior predictive if $i$ was alone | correct Gibbs step |
| Partially observed units | missing coordinate $\mid$ observed coordinate under the cluster's predictive $t$: a $t$ with df $+1$ and the scale inflated by $(\nu+d^2/s_{oo})/(\nu+1)$ | `rcond_t()` | matches |
| $\alpha$ | Escobar–West (1995): $\eta\sim\mathrm{Beta}(\alpha+1,n)$, mixture weight $\pi/(1-\pi)=(a_0+k-1)/\{n(b_0-\log\eta)\}$ | as stated (odds converted to a probability) | matches |

So the code implements Algorithm 3 exactly, with two additions Neal's setting does not need: a Gibbs step for the latent rates on each tomography line, and the Escobar–West update for $\alpha$.

**Brute-force test (`verify_collapsed_exact.R`, output in `verify_collapsed_exact_output.txt`).** With five units, all 52 partitions can be enumerated. The exact posterior of the number of clusters, of $E[\alpha]$, and of a latent unit's $W_1$ is computed from the closed-form NIW marginal likelihood (independent of the sampler's predictive-$t$ code), with $\alpha$ integrated numerically and a partially observed unit integrated over its missing coordinate.

* CAR ($p=2$; one latent, three observed, one partially observed unit): all $|z|\le1.06$; maximum CDF difference for $W_1$ 0.007 (grid resolution 0.00125). PASS.
* NCAR ($p=3$; one latent, four observed units): one chain gave $|z|$ up to 2.8 with batch-means SEs. Sixteen independent chains of $4\times10^5$ sweeps (`verify_collapsed_exact_rep.R`) give all $|z|\le1.7$, e.g. $P(k=1)$ 0.23020 vs 0.23025 exact. PASS.

This is a stronger check than SBC for the update rules themselves: it compares the sampler to the exact posterior, not to uniformity of ranks.

## 7. Algorithm 1 with efficient code (corrected `eco`)

See `benchmarks/eco_2x2/eco_fixed/README.md`. Per sweep, an efficient corrected Algorithm 1 with remixing is faster than the collapsed sampler (0.13–0.46 ms against 0.30–0.99 ms), so the speed advantage over `eco` comes entirely from implementation. Algorithm 3's advantage is mixing: on `census`, 3 of 8 (CAR) and 4 of 8 (NCAR) Algorithm 1 chains stayed for 22,000 sweeps in a state with about 2.5 clusters that no collapsed chain visited.
