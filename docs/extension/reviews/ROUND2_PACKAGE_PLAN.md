# Round 2: package-side plan (EI report, round 2; Phases 1–3 report, item 3)

*Approved by the user (19 Nov 2026) and implemented in eiPartialID 0.3.1. Decisions: D1 = A; D2 = rung-only variance bound for the projection column only. Added at approval: the projection critical value uses q = dim(beta), with the V-hat rank printed only as a diagnostic (theory decision), and a near-boundary note for the Bhatia-Davis allowance (P3.5(c1)). D3 is the theory session's.*

**Basis.**
- Worked from `claude/wizardly-volta-lbh4bm` at 9c46fc3. The referees reviewed e170be9.
- Between e170be9 and 9c46fc3 these files are unchanged: `PLAIN_LANGUAGE_GUIDE.md`, the vignette, `PHASE5a_sign_module.md`, `PHASE2c_monotone.md`, `verify_R10.py` and the whole of `package/`. So the referees' line numbers hold.
- **Ownership:**
  - **mine:** guide §9 (PLG), the vignette (Vig) and package code;
  - **theory session:** 5a, 2c, P3 and `phase3_montecarlo.py`. For these files I give proposed text only.
- **Alignment with the theory session:** R10 becomes pooled M0 only. Every status label below follows that. Grouped M-SIGN confidence intervals, projection included, are not covered by any theorem (R16 is open).

**Checks run for this plan.**

| Check | Result |
|---|---|
| The referee's reproduction script against `verify_R10.py` at e170be9 (NumPy/SciPy, 5 s) | `Literal intersection: 304/600; reference CR-naive: 550/600`. **Confirmed exactly.** |
| W1 example: two equal precincts, A's share .25/.75, A's rate .40→.60, the rest .55→.45 | W1 $\delta_W$ = **+0.0046875**: A's part +0.009375, the rest's part −0.0046875. The package's `.rxc_assoc_W` gives the same value. **Confirmed.** |
| How 0.3.0 computes the Bhatia–Davis bound in `proj_actual` (`R/rxc_infer.R` l.308–330, l.434–437, l.472–473 at 9c46fc3) | m̄ is bounded by the **projection** interval (`attr(iv, "projection")`, at level α_v) from the same candidates. CLR is used only for `ci_actual`. In conservative mode, both columns use the conservative interval (since 6e24c5e). |
| The Monte Carlo script `phase3_montecarlo.py` l.267–282 | `proj_B` uses `S1bd` built from **CLR** bounds (`clr_upper`). P3 §8 l.377 says "as in the package", which is stale since 0.3.0. |

---

## I. EI round-2 report

### EI-1. The W1 example, and pairwise CR beside it

**(1) Agree.**
- The guide (l.215) and the vignette (l.173–176) still say a positive W1 sign means both A's rate and everyone else's rate rise with A's share. That is false; the counterexample is confirmed above.
- The W2 text also never says that `pair` changes only the *between* side: the contextual term stays the pooled W2 sum over all original groups, C included.
- I wrote the current wording myself in 0.3.0; it fixed the predictor but kept the false "both rise".

**(2) Changes.**
- **PLG §9 l.215.** Replace the bullet with:
  > "Using A's share as the predictor for both groups, the within association of A's rate and that of everyone else's rate have a nonnegative weighted sum": `type = "dW"`. Either part may be negative. In two equal precincts with A's share 25% and 75%, A's rate goes .40 → .60 while everyone else's goes .55 → .45, and the sign still holds (+0.0047). A sign on one group's own slope is not available.
- **PLG §9 l.217.** After "A against one other group B (`version = "pair"`, `other = "B"`)", insert:
  > For reinforcement with `"pair"`, only the *between* term compares A with B; the *contextual* term is still the pooled W2 sum over all the original groups, including any third group C. Choosing `"pair"` changes one side of the pairing, not both.
- **Vig l.173–176.** Replace with:
  > the **within** (contextual) association $\delta_W$: with the default one-vs-rest version (`"W1"`), using group $r$'s share as the predictor for both groups, the within-group associations of $r$'s rate and of the rest's rate have a nonnegative (or nonpositive) weighted sum. Either part may have the other sign. In two equal units with $r$'s share .25 and .75, $r$'s rate .40 → .60 and the rest's .55 → .45, W1 $\delta_W = +0.0047$.
- **Vig l.182–183.** After "**contextual reinforcement (CR)** … without saying which", add:
  > With `version = "pair"` the between term compares $r$ with the named group only, while the contextual term is the pooled W2 sum over all original groups.
- **5a §1.4 l.101** (theory session; proposed text). After "(S5).", add:
  > For the pairwise version only $\delta^*_B(r,s)$ is pairwise; $\Delta^*_r$ pools every original group, including those outside $\{r,s\}$.

**(3) Code and tests.**
- The CR label (`.rxc_cr_label`, `R/rxc_sign.R`) for `pair` reads "same sign of delta_B(A, B), Delta_W(A) [W2]". Change it to "same sign of delta_B(A, B) [pair] and Delta_W(A) [W2, pooled over all groups]".
- New test, `test-batchB.R`: the two-precinct W1 example. Check that `.rxc_assoc_W` gives +0.0046875 with a negative "rest" part, so that the documented example is pinned down.

**(4) Effort:** S.

### EI-2. Name the procedure behind each coverage figure

**(1) Agree; reproduced.**
- The "8–10%" is the reference `CR-naive` procedure, which expands for target uncertainty after intersecting.
- The *literal* intersection of the estimated bracket with a completed JKST CI covered 304/600 (miss 49.3%). That fits 5a §3.5 item 4: asymptotic coverage at most 1/2 at this edge.
- The guide and vignette attach 8–10% to the literal procedure.

**(2) Changes.**
- **PLG §9 l.225.** Replace the parenthesis "(in the PHASE5a design, 8–10% of the time instead of 5%; the exact rate depends on the design)" with:
  > How badly depends on how it is done. In the PHASE5a design with 150 precincts and the truth at the edge, intersecting the estimated bracket with a finished confidence interval missed the truth 49% of the time (296 of 600). The theory says such an interval covers at most half the time in large samples. Intersecting first and then widening for the target's uncertainty missed 8% of the time (50 of 600) instead of 5%. These rates are for that design.
- **PLG §9 l.235.** In the last sentence of "Reading the printout", replace "The 8–10% miss rates quoted for the hand-made interval are for one design; the failure itself is proved to persist in large samples." with:
  > The miss rates quoted for the hand-made interval (49% for the literal intersection, 8% with widening afterwards) are for one design; the failure itself is proved to persist in large samples.
- **Vig l.257–264.** Replace "(the 8–10% miss rates in PHASE5a are for its design)" with:
  > (in the PHASE5a edge design with 150 units: the literal intersection covered 304 of 600 times, at most one half asymptotically; intersecting first and widening for target uncertainty afterwards covered 550 of 600)
- **5a §3.4/§3.5 l.360–365** (theory session; proposed text). In item 4, after "at most one half", add:
  > In the 150-precinct positive-edge design the literal intersection covered 304/600, while `CR-naive` in `verify_R10.py`, which expands the intersected interval for target uncertainty, covered 550/600. The 8–10% miss rates in §3.4 refer to `CR-naive` only.

  Also suggested: add the literal count to `verify_R10.py` as V13, using the referee's script.

**(3) Code:** none. The package never computes this interval.

**(4) Effort:** S.

### EI-3. Give both false-belief results

**(1) Agree.**
- "Never contained the truth under a moderate departure" is wrong. At 600 precincts, *s* = .08 covered 110/400; only *s* = .25 covered 0/400.
- At 150 precincts: 399, 336 and 2 of 400.
- No interval was empty, and JKST covered 400/400 throughout (5a §3.4, V12).

**(2) Changes.**
- **PLG §9 l.223.** Replace the last sentence, from "If reinforcement is false…", with:
  > If reinforcement is false the interval stays narrow and never comes out empty, but it can miss. In one referee's experiment with 600 precincts, in which only the contextual slopes changed (by *s*), the conservative interval contained the truth 398 times out of 400 with reinforcement just holding (*s* = 0), 110 times at *s* = .08 and never at *s* = .25. The straight-line rung alone contained it every time. Falsity does not mean every interval misses; it means the coverage guarantee is gone.
- **Vig l.254–256.** Replace "If CR is false, … covered every time." with:
  > If CR is false the interval stays narrow and non-empty but can be wrong. In the referee's experiment (5a §3.4, V12: 600 units, contextual slopes $-s$) the conservative slab covered 398, 110 and 0 of 400 at $s$ = 0, .08 and .25, while the linear rung alone covered 400 each time.

**(3) Code:** none.

**(4) Effort:** S.

### EI-4. Scope the projection guarantee to the fitted model

**(1) Agree.**
- The guide (l.223, "the projection interval, which carries a proof") and the vignette (l.253) do not say which fit the proof covers.
- After the R10 correction, the scope is:
  - **pooled M0 without sign/CR rows:** P3.5(b2)(ii), with $\hat\Sigma$ nonsingular, or the reduced rank as implemented by the guard;
  - **pooled M0 with sign/CR rows:** R10 (b2)(ii);
  - **grouped rungs without M-SIGN:** P3.5 in regime R1 only (few, large groups); with many small groups (R2) it is a simulation finding;
  - **grouped rungs with M-SIGN:** no theorem (R16 is open).
- The package prints one projection label for every fit. It also gives grouped M-SIGN fits a sharpness warning but no statement about their CIs.

**(2) Changes.**
- **PLG §9 l.211** (status paragraph). Append:
  > Confidence intervals with sign or reinforcement beliefs carry a proof only for the pooled model M0 (PHASE5a R10, for the projection and conservative intervals). With groups they have no proof yet (R16); their coverage is untested.
- **PLG §9 l.223.** "…the package prints the projection interval, which carries a proof, beside it" becomes:
  > …beside it the package prints the projection interval, which carries a proof for the pooled model M0 when the regression's estimated covariance matrix is nonsingular (the package checks this). It carries no proof for grouped fits with these beliefs.
- **Vig l.194–197.** Replace "but sharpness is proved only for the pooled model M0, and a note says so at fit time" with:
  > but sharpness (R3(b)) and inference (R10) are proved only for the pooled model M0. For grouped fits with signs, no confidence interval, the projection included, is covered by a theorem (R16 open), and the fit prints this.
- **Vig l.252–254.** "`print()` shows the projection interval beside it" becomes:
  > `print()` shows the projection interval beside it, proved for pooled M0 with a nonsingular estimated covariance (R10), not for grouped fits
- **Vig "What is reported" (l.~45).** After "…usually only slightly wider." add one sentence on scope: for pooled M0 in any full-rank design; for grouped models with few, large groups (regime R1); a simulation finding with many small groups; and no theorem when sign or CR restrictions are combined with groups.
- **5a R10 l.241 and 2c M5 l.184** (theory session). This follows from their scope correction. Proposed status row for 5a §5:
  > Package labels (0.3.x): `print()` states the CI status per fit — M0: P3.5(b2)(ii)/R10 with rank q; grouped without M-SIGN: R1 only; grouped with M-SIGN: no theorem (R16).

**(3) Code** (`R/rxc_api.R`, `print.eiRxC` and the notes):
- Store the projection rank: `est$proj_rank`, the numerical rank of $\hat V$, beside `dim(beta)`.
- Replace the single `proj_*` legend with one line chosen by the fit:
  - pooled M0, no sign/CR rows: `proj_*: projection intervals; asymptotic guarantee (P3.5(b2)(ii)): pooled model, V-hat rank q of p`, with `(rank guard dropped k candidates)` added when k > 0;
  - pooled M0 with sign/CR rows: the same line, citing `PHASE5a R10 (pooled M0)`;
  - grouped, no M-SIGN: `asymptotic guarantee only with few, large groups (P3 regime R1); with many small groups a simulation finding`. This merges with the existing few-units note.
  - grouped with M-SIGN: `proj_* and ci_*: no coverage theorem for sign/CR restrictions under grouped rungs (R16 open)`.
- The grouped M-SIGN fit-time note gets, after "sharpness … proved only for the pooled model M0":
  > ; nor are its confidence intervals, the projection and conservative ones included (R10 covers pooled M0; R16 open)
- The roxygen of `eiBounds` (Details, status paragraph): the same scoping.
- **Tests:**
  - a pooled M0 fit prints "rank" and "P3.5(b2)(ii)";
  - an M0 + CR fit prints "R10 (pooled M0)";
  - an M1 + CR fit prints "no coverage theorem" and its note contains "nor are its confidence intervals";
  - an M1 fit without restrictions prints "regime R1".

**(4) Effort:** S–M.

### EI-also. Open items visible; rung-only comparison as a joint effect

**(1) Agree.**
- The `base_*` block removes *all* optional restrictions at once. With several restrictions, one belief's own increment is the added-last gain in `eiRedundancy()`.
- The open items are listed in NEWS but not in one status key in the guide or vignette.

**(2) Changes.**
- **PLG §9 l.222** ("Compare with the rung alone"). Append:
  > With several beliefs, this comparison shows the joint effect of removing them all. For what one belief adds on top of the others, use the "added last" gain from `eiRedundancy(fit)`.
- **Vig l.232–234.** Add the same sentence.
- **PLG §9, new last bullet "Not yet available":**
  > a held-out benchmark with known truth and signs chosen in advance; the drop-one display (each belief removed in turn, intervals shown); whole-number (integer) tables; inference on the width change a belief buys; a direction for each group under W2; and proofs for grouped rungs (R16).
- **Vig:** the same list, at the end of the "How to read" bullets.

**(3) Code:**
- In `print.eiRxC`, when more than one optional restriction is given, the base block header becomes:
  > The rung alone (M0, all N optional restrictions removed; joint effect: see eiRedundancy()'s added-last gains for each one):
- **Test:** with two restrictions, the header contains "all 2 optional restrictions removed".

**(4) Effort:** S.

---

## II. Phases 1–3 report, item 3: the projection interval for B and its Bhatia–Davis allowance (code and API)

**What 0.3.0 does** (verified at 9c46fc3):
- `proj_actual` = projection interval for the kept part (level α_est) ± z_{α_E}·S, where S is the Bhatia–Davis sd bound.
- In S, m̄ (the π-weighted conditional-mean rate) is bounded by the **projection** interval at level α_v, built from the same candidates (`.rxc_fpc_sd`, attribute `"projection"`). Popoviciu is used if that interval is missing or empty.
- `ci_actual` (CLR) uses the CLR bound for m̄. In conservative mode, both columns use the conservative interval.

**(1) Partly agree.**
- **The referee is right about P3 §8 and the Monte Carlo.** `phase3_montecarlo.py` l.267–282 builds `proj_B` with a CLR-based BD bound, and P3 §8 l.377 says "as in the package", which has been stale since 0.3.0. The practical recommendations (§9) then call that projection interval "guaranteed".
- **The concern does not apply to the package's `proj_actual`.** Its variance bound comes from "a set with formally established coverage", the referee's own second acceptable option. P3.5(b2)(ii) covers the m̄ projection interval at level α_v, under its conditions: fixed dimension, $\hat\Sigma$ nonsingular or the reduced parameterisation. On that event, $\sum\pi_i m_i(1-m_i) \le \bar m(1-\bar m) \le \max_{m\in \mathrm{CI}} m(1-m)$. A union bound over (α_est, α_E, α_v) then gives the stated level. So `proj_actual` is provably valid wherever the projection interval for B̃ is, and the residual step needs the same CLT repair either way (P-B41, Lindeberg/Chebyshev).
- **Two residual gaps:**
  1. **m̄ is bounded on the "rung + linear rows" system** (`base` in `.rxc_estimate`), which includes `assoc` sign rows.
     - Pooled M0 with sign rows: its validity is R10 (b2)(ii), which is fine.
     - Grouped rungs with sign rows: R16 is open, so even the variance bound inherits an unproved step.
     - Computing m̄'s interval on the **rung alone** (no optional rows) removes this. That interval is a superset under the null, so validity needs only P3.5(b2) for the rung. It is slightly wider, and only in sign fits.
  2. The projection for m̄ and the projection for B are two separate events with separate budgets (α_v and α_est). That is valid; a single joint event could save α_v but is not needed now.
- **The B̃ projection interval (`proj_model`) is unaffected**, as the referee says.

**(2) Proposed API, with a recommendation.**
- **Option A (recommended): keep one `proj_actual`, BD with m̄ from a projection interval; make it independent of optional rows; label its provenance.**
  - Default: `fpc = "bhatia-davis"`. `proj_actual` uses BD with m̄ from a rung-only projection interval.
  - `fpc = "popoviciu"` (already available) gives Popoviciu in every column.
  - Print legend:
    > `proj_actual: projection interval for B + unit-level allowance (Bhatia–Davis, variance bound from a projection interval at level 0.005; P3.5(b2)(ii) and a union bound). fpc = "popoviciu" uses the bound 1/4.`
  - `ci_actual` (CLR) keeps the CLR-based bound, labelled:
    > `ci_actual: CLR (algorithmic); its Bhatia–Davis bound also uses a CLR interval`
  - Why recommended: it gives a single provably valid column for B under the same conditions as for B̃, and no extra column. It is narrower than Popoviciu only where it matters (boundary designs; P3 §10 item 9: 0.057 against 0.086 in design C).
- **Option B (the coordinator's suggestion): `proj_actual` = Popoviciu (always valid), plus `proj_actual_adaptive` = BD from a projection interval, labelled conditional.**
  - Under 0.3.0's construction, the "adaptive" column is not conditional: it has the same formal status as Option A's column. Labelling it "conditional" would understate it.
  - The extra column costs print width for no gain in rigour.
  - If the user prefers maximum simplicity of the guarantee (no α_v split, no second interval), choose B but call the column `proj_actual_bd`, with the legend "Bhatia–Davis; valid under P3.5(b2)(ii) for the variance bound as well".
- **Option C:** Popoviciu-only projection column, dropping BD there. It is simplest, and wider by up to 0.03 in boundary designs.

**(3) Code and tests (Option A).**
- `.rxc_estimate`: compute `fpc_cache` from the rung-only system (`st0 <- st; st0$extra <- NULL; .rxc_constraints(st0)`), not `base`. This affects only fits with linear optional rows (`assoc`, `order`, `sign`, `slope`, …). Their `proj_actual` and `ci_actual` get slightly wider; M-SIGN reference values for estimates and `ci_model` do not change.
  - Alternative: rung-only for the projection column only, keeping `base` for CLR `ci_actual` so that its numbers do not change. **User decision.**
- Print legend as above. Roxygen `@param fpc` and the Details status paragraph get the same text.
- **Tests:**
  - with an `assoc` row, the BD sd used for `proj_actual` equals the one computed on the rung alone;
  - `proj_actual` with `fpc = "popoviciu"` contains `proj_actual` with BD in a design where m̄'s interval excludes .5 (design C-like), and equals it when the interval contains .5;
  - the print legend names "variance bound from a projection interval".
- **Theory session (proposed text, P3 §8 l.377):** replace the bullet with:
  > In this Monte Carlo the projection interval for $B$ takes its Bhatia–Davis bound from a CLR interval; only the Popoviciu version, or a BD bound taken from a projection interval (as the package does since 0.3.0), is covered by P3.5(b2). The package's `proj_actual` is the latter.

  Also: either rerun the `proj_B` column in `phase3_montecarlo.py` with the bound from `proj_upper` at level α_v (one-line change at l.272–273), or relabel it "projection + CLR-based BD (not covered)". In P3 §9 (recipe), the projection interval for B is "guaranteed" only for the projection-based or Popoviciu allowance.

**(4) Effort:**
- Option A code, tests and labels: S.
- MC rerun (theory session): M, about 40 minutes of compute.

---

## III. Does report I need code or test changes beyond wording?

Yes, five small ones, listed above:
1. **EI-1:** the pair-CR label names both sides; a regression test of the W1 two-precinct example.
2. **EI-4:**
   - per-fit CI-status legend, with the V-hat rank for M0 and R10 for M0 with sign rows;
   - grouped M-SIGN fits: "no coverage theorem" in print and notes;
   - tests.
3. **EI-also:** the base-block header says "all N restrictions removed; joint effect", and points to the added-last gains; a test.
4. **II (linked):** the BD m̄ interval computed on the rung alone. This applies to both `ci_actual` and `proj_actual`, or to `proj_actual` only; see decision D2.
5. **EI-2 and EI-3:** no code, since the package never computes the naive interval. The 304/600 count belongs in `verify_R10.py` (V13, theory session).

No numeric output changes except under item 4: slightly wider actual-rate intervals, only in fits with linear optional rows.

**Total effort:** S–M, about one session, plus the theory session's MC relabel or rerun.

## Decisions for the user

- **D1. API for `proj_actual`** (item II).
  - **A (recommended):** one column; BD with the variance bound from a rung-only projection interval; provably valid under P3.5(b2)(ii); `fpc = "popoviciu"` for the 1/4 bound.
  - **B:** Popoviciu `proj_actual` plus a separate BD column.
  - **C:** Popoviciu only.
- **D2. Scope of the rung-only variance bound.** Apply it to the projection column only, so CLR `ci_actual` is unchanged (**recommended**: numbers move only where the proof needs it). Or apply it to both columns.
- **D3. P3 §8 Monte Carlo.** Rerun `proj_B` with a projection-based BD bound, or relabel the existing column as not covered (theory session's call; relabelling is enough for the response).
