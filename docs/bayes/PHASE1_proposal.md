# Phase 1: Bayesian nonparametric EI after Imai, Lu & Strauss (2008): proposal and assumptions

*This proposal covers Imai, Lu & Strauss (2008, Political Analysis 16:41–69; ILS below), the CRAN package `eco` (4.0-6), and a new stand-alone package, working name `eiDPM`, that re-implements the ILS models in Stan and in a C++ blocked Gibbs sampler.*

*Status labels:*
* ***verified*** = checked against the source code or by a run;
* ***to prove*** = a precise statement whose proof belongs in Phase 2;
* ***conjecture***;
* ***heuristic***.

*Nothing here is a derivation.*

*Decisions taken with the owner (2026-11-01):*
* the package is separate from `eiPartialID`, with no dependency in either direction for now;
* Stan is used through `rstan`;
* both samplers are prototyped in 2×2 first and benchmarked against `eco`;
* boundary data follow `eco`'s rule;
* R×C is designed but not yet built.

---

## 1. The ILS models (2×2)

Unit $i$ has composition $X_i$, outcome share $Y_i$ and unknown rates $W_i=(W_{i1},W_{i2})$ with $Y_i=X_iW_{i1}+(1-X_i)W_{i2}$. $W^*=\operatorname{logit}W$ elementwise.

| ILS label | Assumption | Parametric model | Nonparametric model |
|---|---|---|---|
| CAR | $W\perp X$ | $W^*\sim N_2(\mu,\Sigma)$ | $W^*\sim$ DP mixture of $N_2$ |
| CCAR | $W\perp X\mid Z$ | $W^*\sim N_2(Z^\top\beta,\Sigma)$ | — |
| NCAR | none on $(W,X)$; model the joint law | $(W^*,X^*)\sim N_3(\mu,\Sigma)$ | $(W^*,X^*)\sim$ DP mixture of $N_3$ |

Priors (ILS §3; `eco` defaults):
* $\mu\mid\Sigma\sim N(\mu_0,\Sigma/\tau_0)$ and $\Sigma\sim\text{InvWishart}(\nu_0,S_0)$, with $\mu_0=0$, $\tau_0=2$, $\nu_0=4$, $S_0=10I$. The parametrisation is read off ILS App. B: the posterior scale is $S_0+\text{scatter}$, so $S_0$ is the inverse-Wishart *scale*.
* The DP base measure is the same normal–inverse-Wishart law.
* Concentration $\alpha\sim\text{Gamma}(1,0.1)$ (shape, rate), with mean 10.

## 2. What `eco` actually provides (verified against the 4.0-6 source)

* **Exported:**
  * `eco` (parametric CAR, and NCAR with `context = TRUE`);
  * `ecoNP` (DP CAR / NCAR);
  * `ecoML` (EM + SEM, CAR / NCAR);
  * `ecoBD`, the deterministic bounds. `ecoBD` is the **only exported R×C function**.
* **Unexported:**
  * `eco:::ecoRC`: an R×C parametric logit-normal CAR Gibbs sampler (C routines `cBase2C`, `cBaseRC`), undocumented;
  * `eco:::ecoX` in `R/ecoCV.R`: the covariate (CCAR) sampler.
* **No nonparametric R×C exists in `eco`,** and none was found elsewhere on CRAN:
  * `eiPack` (Rosen–Jiang–King–Tanner multinomial–Dirichlet with covariate-dependent parameters) is parametric;
  * `RxCEcolInf` (Greiner–Quinn logistic-normal) is parametric and was archived from CRAN in 2022;
  * `ei` uses `eiPack` for R×C.
* **The latent-rate draw.**
  * The Gibbs samplers (`gibbsBase.c`, `gibbsDP.c`) draw $W_i$ on its tomography line from $\phi(\operatorname{logit}w)/\prod_j w_j(1-w_j)$, uniformly in $w_1$, using independence Metropolis–Hastings with a uniform proposal by default, or a grid. This is the correct conditional density (§4.1).
  * `ecoML` instead integrates the normal kernel of $W^*$ against **Euclidean arc length in logit space**: `pfact = sqrt(W1*'^2 + W2*'^2)` in `fintegrate.c`. That factor is not proportional to the correct Jacobian $\lvert W_1^{*\prime}\rvert\lvert W_2^{*\prime}\rvert$ along the line. *(The source reading is verified; whether it biases the ML estimates is not checked. Not needed for the nonparametric path; recorded as an open item.)*
* **Boundary rule** (`checkdata.R`, `gibbsDP.c`):
  * units with $X\in\{0,1\}$ have one rate observed ($W_1=Y$ if $X=1$; $W_2=Y$ if $X=0$);
  * units with $Y\in\{0,1\}$ have both rates fixed;
  * every exact 0 or 1 is replaced by $10^{-4}$ or $1-10^{-4}$ before the logit.

  In the example data: `reg` has 7 units with $Y=1$; `census` has none.
* **Speed** (verified; one core; 1,000 draws; `census`, $n=1040$):

  | Model | Time |
  |---|---|
  | parametric CAR | 0.8 s |
  | parametric NCAR | 0.7 s |
  | `ecoNP` CAR | 195 s |
  | `ecoNP` NCAR | 257 s |

  The slow case is the DP sampler: a collapsed Pólya-urn update for every unit against every cluster, plus a "remixing" step.

## 3. Identification: where the ILS models stand

The rungs are those of `docs/extension` (DD, M0, L5, L6, …). Here they serve only to locate the ILS assumptions; the new package does not use them.

**C1 (to prove).** If $W\perp(X,N)$, the CAR model class is contained in L6 (Goodman constancy), and the aggregate $B^*_1=E[NXW_1]/E[NX]$ equals $E[W_1]$. CCAR with $Z$ = district dummies is contained in L5. *(The condition on $N$ matters: ILS state independence of $X$ only.)*

**C2 (to prove).** Under nonparametric CAR with $W\in[0,1]^2$, if $\operatorname{supp}P_X$ is infinite, the whole law $P_W$ is point-identified.
* *Sketch:* $E[Y^k\mid X=x]=\sum_j\binom kj x^j(1-x)^{k-j}E[W_1^jW_2^{k-j}]$ is a polynomial of degree $k$ in $x$. Its values at $k+1$ distinct $x$ fix all mixed moments of order $k$, and a law on a compact set is determined by its moments.
* *Consequences:*
  * population quantities are identified under CAR without any parametric form;
  * the conditional law of $W_i$ given $(X_i,Y_i)$ is identified, so sample (in-sample) EI has an identified predictive distribution;
  * the nonparametric CAR model is therefore not prior-driven asymptotically. Its weak point is CAR itself.

**C3 (to prove).** Under nonparametric NCAR (the joint law of $(W,X)$ unrestricted), the identified set of $E[W_1]$, or of any aggregate, is the averaged Duncan–Davis interval.
* *Sketch:* any rule that places $W_i$ anywhere on its line, e.g. always at the lower bound, defines a valid joint law of $(W,X,Y)$ with the observed law of $(X,Y)$. This is the 2×2 analogue of `docs/extension` Prop. D0.

**C4 (to prove).** Parametric NCAR implies that $x\mapsto E[W_r\mid X=x]$ is monotone, with the sign of $\operatorname{Cov}(W_r^*,X^*)$.
* *Reason:* $W_r^*\mid X^*$ is normal with linear mean and constant variance, and $\operatorname{expit}$ is increasing.
* So parametric NCAR is **not** nested in linear contextual effects (M0): its mean function is sigmoidal.

**C5 (conjecture).** Under DP-NCAR, the posterior of an aggregate does not concentrate at a point as $n\to\infty$, and its limit depends on the prior through $\alpha$ and the base measure.

**H1 (heuristic): why DP-NCAR still gives narrow, often good, answers.**
* Within one mixture component, $W^*\mid X^*$ is normal with linear mean: rates vary smoothly with composition.
* Units with similar $X$ tend to share components, so the model behaves like *local CAR*: within a window of compositions, $W$ is treated as nearly independent of $X$. Different $x$ inside the window give different line slopes, and C2's mechanism then applies locally.
* The prior sets the window width. The data cannot check it: by C3 they are equally consistent with every placement on the lines.
* ILS's good literacy results (their Table 5) fit with smooth true contextual effects in that application.

**R×C.**
* C2 should extend: the moments of $\beta^\top x$ are polynomials in $x$ *(conjecture)*.
* C3 extends directly: the per-unit fibre replaces the line *(to prove)*.
* The case for a nonparametric NCAR in R×C is weaker. Data constrain a $(C-1)$-dimensional projection of an $R(C-1)$-dimensional object per unit, so the prior fills in more. ILS's own conclusion says as much.

## 4. The 2×2 prototypes

### 4.1 The density on the tomography line (to prove; check numerically in Phase 4)

With $L_i,U_i$ the Duncan–Davis bounds of $W_{i1}$ and $w_2(w_1)=(Y_i-X_iw_1)/(1-X_i)$:
$$p(Y_i\mid X_i,\theta)=\frac1{1-X_i}\int_{L_i}^{U_i}f_W\big(w_1,w_2(w_1)\mid\theta\big)\,dw_1,\qquad f_W(w)=\frac{\phi_\Sigma(\operatorname{logit}w-\mu)}{\prod_jw_j(1-w_j)}.$$

So the augmented posterior over $(\theta,W_{11},\dots,W_{n1})$ is the prior times $\prod_if_W(W_{i1},w_2(W_{i1}))$. The factor $1/(1-X_i)$ is a data constant.
* **Stan:** declare $W_{i1}$ with vector bounds $(L_i,U_i)$; Stan adds the bounding transform's Jacobian itself.
* **NCAR:** replace $f_W$ by the conditional density of $W$ given $X^*_i$ under the component. The marginal of $X^*$ enters the likelihood for the $X$-block of $(\mu,\Sigma)$.

### 4.2 Samplers

Both use the truncated stick-breaking form of the DP (Ishwaran–James): $v_k\sim\text{Beta}(1,\alpha)$ for $k<K$, $v_K=1$, and $\pi_k=v_k\prod_{l<k}(1-v_l)$.
* With $\alpha$ near 10, the expected leftover mass $(\alpha/(1+\alpha))^{K-1}$ is about 0.01 at $K=50$. The default is $K=50$, and the posterior mass of the last component is reported as a truncation check.

**(A) Stan, latent rates.**
* Parameters: $W_{i1}$ per unit, plus $v$, $\alpha$, and $(\mu_k,\Sigma_k)$ per component.
* The component indicators are summed out: $\log f(W_i)=\operatorname{logsumexp}_k[\log\pi_k+\log\phi_{\Sigma_k}(W_i^*-\mu_k)]$.
* Cost per gradient is about $nK$ evaluations of a 2-D (CAR) or 3-D (NCAR) normal.
* Variants: parametric CAR and NCAR ($K=1$).
* *Risk (heuristic):* HMC cannot jump between cluster configurations. Label switching is harmless for the reported functionals; genuine multimodality is not. Diagnostics are reported on functionals of $W$ only.

**(B) C++ blocked Gibbs (Rcpp).** Each sweep:
1. Update $W_{i1}$ on its line given its component (slice sampling, tuning-free).
2. Draw $z_i$ with probability proportional to $\pi_k\phi(W_i^*\mid\mu_k,\Sigma_k)$.
3. Draw $v_k\mid z\sim\text{Beta}(1+n_k,\ \alpha+\sum_{l>k}n_l)$.
4. Draw $(\mu_k,\Sigma_k)$ from its conjugate NIW posterior.
5. Draw $\alpha\sim\text{Gamma}\big(1+K-1,\ 0.1-\sum_{k<K}\log(1-v_k)\big)$.

Cost is about $nK$ per sweep, with no remixing step. NCAR uses the 3-D NIW with $X^*$ observed, and step 1 uses the conditional normal of $W^*$ given $X^*$ within the component.

### 4.3 Data rules

These follow `eco` (§2):
* exact 0/1 are replaced by $10^{-4}$ or $1-10^{-4}$;
* $Y\in\{0,1\}$: both rates are fixed;
* $X\in\{0,1\}$: one rate is observed and enters only through its marginal density.

To keep HMC away from $w_2\in\{0,1\}$ at the segment ends, the bounds are shrunk by $10^{-8}(U_i-L_i)$.

## 5. Benchmark (Phase 4)

**Data.** `census` (literacy, $n=1040$) and `reg` (registration, $n=275$). Both ship the true unit rates $W_1,W_2$.

**Models.** CAR and NCAR, each parametric and DP. For each: `eco`/`ecoNP`, the Stan version (A) and the Gibbs version (B). The parametric models use (A) only, because `eco`'s parametric sampler is already fast.

**Reported** (all with posterior means, 90% intervals and Monte Carlo SE):
* **Aggregates:** the in-sample unweighted means of $W_1$ and $W_2$, and the $N$-weighted aggregates $\sum_iN_iX_iW_{i1}/\sum_iN_iX_i$ and the analogue for $W_2$, against the truth and the averaged Duncan–Davis interval.
* **Unit level:** RMSE of the posterior-mean $W_i$ against the truth; coverage of the unit 90% intervals; the correlation of unit posterior means with `eco`'s.
* **Agreement with `eco`:** the difference in aggregate posterior means, in units of the combined Monte Carlo SE.
* **Speed:** wall time and effective sample size per second for the aggregates (`posterior::ess_bulk`). Stan's compile time is reported separately.
* **Diagnostics:** $\hat R$ (4 chains for Stan; 2 chains for Gibbs and `eco`), divergences, and the truncation mass.

**Pass criteria for agreeing with `eco`** (stated in advance):
* aggregate posterior means within 3 combined MC SEs, or within 0.005 if that is larger;
* unit posterior means correlated above 0.98 with `eco`'s.

A DP sampler that fails these is reported as failing, not tuned until it passes.

## 6. R×C (design only; not in this round)

**Latent table.** In unit $i$, the free cells $n_{rc}$, $r<R$, $c<C$, are filled in sequence. Each cell is drawn as $u_{rc}\in(0,1)$ mapped linearly into its Fréchet bounds given the cells already filled.
* The map is triangular, so $\log\lvert J\rvert=\sum\log(\text{ub}_{rc}-\text{lb}_{rc})$.
* Parametrising the fibre by its free cells is an affine chart, so the fibre's flat measure differs only by a data constant *(to prove)*. In 2×2 this reduces to §4.1.
* The bounds are continuous with kinks, so divergent transitions are a risk to test.

**Rates.** Additive log-ratios of each row's rate vector, plus those of $X$ for NCAR: $R(C-1)+(R-1)$ dimensions (8 for 3×3).

**Leaner NCAR variant.** A mixture of regressions for $\operatorname{alr}\beta\mid\operatorname{alr}X$, with diagonal or low-rank component covariances, instead of the full joint. It keeps the NCAR dependence while cutting parameters.

**Benchmark.** The `ei.Datasets` NZ tables (electorate-level truth) once the 2×2 prototypes pass.

## 7. Verification and code standards

**Tests.** `testthat` tests:
* the line density against direct numerical integration;
* recovery of known parameters on simulated data;
* the boundary rule.

**Benchmark scripts.** Committed under `benchmarks/eco_2x2/`, with their output.

**Packaging.** `R CMD check` must be clean, with `rstan` in Suggests. The models are compiled at first use, and the Gibbs sampler needs only Rcpp.

## 8. Open items

* Whether `ecoML`'s arc-length measure biases its estimates (§2).
* Proofs of C1–C4; a counterexample or proof for C5.
* DP versus an overfitted finite mixture with a sparse Dirichlet prior (often better behaved under HMC). Decide after the 2×2 benchmark.
* R×C (§6).
