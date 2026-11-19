# Response to the referee, second round: Batch B, ecological inference

**Reviewed version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b).
**Revised version:** commit [`1b8b6fb`](https://github.com/mdweaver/eiid/tree/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f) on branch `claude/wizardly-volta-lbh4bm`.
**All changes since the reviewed version:** [compare `e170be9…1b8b6fb`](https://github.com/mdweaver/eiid/compare/e170be9161b9bbc250a17a27ffd39bc523c1c88b...1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f).
**Your report, as received:** [`reviews/B_ecological_inference_report_round2.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/reviews/B_ecological_inference_report_round2.md).

Abbreviations:
* **Guide** is [`PLAIN_LANGUAGE_GUIDE.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PLAIN_LANGUAGE_GUIDE.md) §9.
* **Vignette** is [`eiPartialID_rxc.Rmd`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/package/eiPartialID/vignettes/eiPartialID_rxc.Rmd).
* **5a** is [`PHASE5a_sign_module.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE5a_sign_module.md).
* **2c** is [`PHASE2c_monotone.md`](https://github.com/mdweaver/eiid/blob/1b8b6fb6677b30e50cc8f9fba9f548e31cc8718f/docs/extension/PHASE2c_monotone.md).

Tags: **[Theory]**, **[Code]** (the package is now 0.3.1), **[Docs]**, **[Check]**. Line numbers below refer to the revised version.

---

## 1. Summary

Thank you for the second report, for rerunning the scripts, and for the recommendation of minor revision. We made all four corrections, with no disputes, and the guide, the vignette and the package's printed labels now agree with the revised mathematics. We reproduced both of your new numbers exactly. Your reproduction script gives 304/600 against 550/600, and it is now check V13 in `verify_R10.py`. The W1 example gives +0.0046875, and the package's own association function returns the same value. Thank you also for accepting the correction about sign slack.

## 2. What changed

| Where | Change | Your item |
|---|---|---|
| Guide l.215, l.217; vignette l.189, l.197 | W1 explained as a nonnegative weighted sum, "either part may be negative", with your two-precinct example; `pair` changes only the between side, the contextual term stays pooled over all groups | 1 |
| 5a §1.4 (pairings) | "Only $\delta^*_B(r,s)$ is pairwise; $\Delta^*_r$ pools every original group, including those outside $\{r,s\}$" | 1 |
| Package label | Pairwise CR label now names both sides: "delta_B(A, B) [pair] and Delta_W(A) [W2, pooled over all groups]" | 1 |
| Guide l.225, l.235; vignette l.282; 5a §3.5 item 4 | Each coverage figure attached to its own procedure: literal intersection 304/600, intersect-then-widen 550/600 | 2 |
| Guide l.223; vignette l.276 | All three false-belief results with parameter values and counts | 3 |
| Guide l.211, l.223; vignette l.49, l.213, l.273; 5a R10, §5; 2c M5 | Projection guarantee scoped to the fitted model; CI status per kind of fit | 4 |
| Package `print()` | A "status:" line names the result covering each fit's intervals; grouped sign fits print that no theorem covers them | 4 |
| Guide l.222, l.237; vignette l.251, l.298; package base-block header | Rung-only comparison described as the joint effect of removing all beliefs; added-last gain for one belief; a "Not yet available" list | final paragraph |
| Checks and tests | V13 in `verify_R10.py`; new tests in `tests/testthat/test-round2.R`; all pass | 1, 2, 4 |

## 3. Point-by-point

### Correction 1: the W1 example and pairwise CR

Agreed. In 0.3.0 we fixed the predictor but kept the false "both rise".

* **[Docs]** The guide's W1 bullet (l.215) now reads, following your suggestion:
  > "Using A's share as the predictor for both groups, the within association of A's rate and that of everyone else's rate have a nonnegative weighted sum. Either part may be negative."

  It continues with your example: two equal precincts at 25% and 75% A, A's rate .40 → .60, everyone else's .55 → .45, and the sign still holds (+0.0047). The vignette (l.189) says the same.
* **[Docs]** The guide (l.217) and the vignette (l.197) now say that with `version = "pair"` only the *between* term compares A with B. The contextual term is still the pooled W2 sum over all the original groups, including any third group C. So choosing `pair` changes one side of the pairing, not both.
* **[Theory]** 5a §1.4 adds the same statement to the definition of the pairings.
* **[Code]** The printed label for pairwise CR (`R/rxc_sign.R`) now names both sides.
* **[Check]** Two new tests:
  * "a positive W1 sign does not mean both groups' rates rise", which checks your example: total +0.0046875, A's part +0.009375, the rest's part −0.0046875;
  * "the pairwise CR label names both sides of the pairing".

### Correction 2: the coverage figure for the procedure described

Agreed, and thank you for the comparison. We reran your script against the pinned `verify_R10.py` and got 304/600 and 550/600 exactly.

* **[Docs]** The guide (l.225) now separates the two procedures.
  * Intersecting the estimated bracket with a finished confidence interval missed the truth 49% of the time (296 of 600), and the theory says such an interval covers at most half the time in large samples.
  * Intersecting first and then widening for the target's uncertainty missed 8% of the time (50 of 600).
  * Both figures are for that design.

  "Reading the printout" (l.235) and the vignette (l.282) give the same figures.
* **[Theory]** 5a §3.5 item 4 gives both counts and says the 8–10% miss rates in §3.4 refer to `CR-naive` only. The warning itself is unchanged.
* **[Check]** `verify_R10.py`, V13 [V]: your literal-intersection count beside `CR-naive`, using the same random stream.
* The package never builds this hand-made interval, so no code changed.

### Correction 3: both false-belief results

Agreed. "Never contained the truth under a moderate departure" was wrong.

* **[Docs]** The guide (l.223) now reports your experiment with 600 precincts, in which only the contextual slopes change, by $s$. The conservative interval contained the truth:
  * 398 of 400 times with reinforcement just holding ($s=0$);
  * 110 times at $s=.08$;
  * never at $s=.25$.

  It was never empty, and the linear rung alone covered every time. The guide then says: "Falsity does not mean every interval misses; it means the coverage guarantee is gone." The vignette (l.276) gives the same figures. 5a §3.4 already contains the table.

### Correction 4: scope the projection guarantee to the fitted model

Agreed. The partial-identification referee raised the same point from the theory side, so we fixed both together.

* **[Theory]** 5a Theorem R10 now covers **pooled M0 only**. The fixed-district intervals of P3 cover the rung alone and are background. Sign and CR rows under district rungs are open item R16. 2c M5 and a new 5a §5 row set out the CI status for each kind of fit.
* **[Docs]** The guide (l.211, l.223) and the vignette (l.49, l.213, l.273) now say:
  * the projection interval carries a proof for the pooled model M0 when the estimated covariance is nonsingular (P3.5(b2)(ii), and R10 with sign or CR rows);
  * with groups and no sign rows, a proof only with few, large groups;
  * with groups and sign or CR rows, no theorem yet (R16), and coverage there is untested.
* **[Code]** `print()` now prints a "status:" line chosen by the fit:
  * pooled M0: "asymptotic guarantee (P3.5(b2)(ii)), pooled model; V-hat rank r of p (diagnostic)", with R10 named when sign or CR rows are present;
  * grouped without signs: "only with few, large groups (P3 regime R1); with many small groups a simulation finding";
  * grouped with signs: "no coverage theorem for sign/CR restrictions under grouped rungs (R16 open)".

  A grouped M-SIGN fit also prints a note that its confidence intervals, projection and conservative included, are not covered. The help page for `eiBounds()` says the same. Tests cover each case.

### Final paragraph: open items and the rung-only comparison

Agreed on both.

* **[Docs] [Code]** When several optional restrictions are entered, the base block is labelled as "the rung alone (M0, all N optional restrictions removed: their joint effect; for one restriction's own increment see the added-last gains in eiRedundancy())". The guide (l.222) and the vignette (l.251) explain the difference. A test checks the label.
* **[Docs]** The guide (l.237) and the vignette (l.298) have a new "Not yet available" list, so their status stays visible:
  * a held-out benchmark with known truth and signs chosen in advance;
  * the drop-one display;
  * integer tables;
  * inference on the width change a belief buys;
  * a per-group W2 direction;
  * proofs for grouped rungs.

## 4. Checks

* `verify_phase2c.py` 29/29 PASS.
* `verify_R10.py`: V9, V10, V11 and V13 PASS; main table unchanged.
* The package (0.3.1) passes `R CMD check` with 0 test failures. No reference output changed. The only output differences are the new labels and notes.

## 5. Still open

These are visible in the guide's new list and in the package's printed notes:

* the held-out benchmark with pre-specified signs;
* the drop-one display;
* integer tables (R15);
* inference on width changes;
* per-group W2;
* grouped-rung proofs (R16).

M-SIGN stays labelled experimental.
