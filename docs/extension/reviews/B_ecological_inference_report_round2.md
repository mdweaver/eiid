# Second-round referee report on M-SIGN

**Recommendation: minor revision for the experimental release.** I would close the major revision once the guide, vignette and status labels agree with the revised mathematics. The authors corrected the two errors that affected the first report's substantive conclusions about attainability and interval width. The remaining issues could still lead an applied researcher to impose a stronger belief than intended or quote coverage for the wrong procedure. They can be addressed without changing the pooled estimator.

I compared the response with commit `e170be9161b9bbc250a17a27ffd39bc523c1c88b`, reran both revised Python verification scripts, and compared two plug-in procedures on the same simulated elections. I inspected the cited R tests. R was unavailable in this environment, so the print and entry-point tests were not independently executed here. [Revised commit](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b)

## Findings I would close

* **Continuous and counted tables.** M2 now states what its LP attains. M3 distinguishes a weak between sign, which still has an integer witness on the appropriate side in the two-row margin-only problem, from reinforcement, which can exclude every integer table. The revised text also handles pure units and distinguishes a population contradiction, an empty plug-in and a calibrated rejection. The W2 definition and sign-reversal example are correct. The revised 2c script passed all 26 checks. [2c, §§1–3](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2c_monotone.md#L49-L56)
* **Union length and false reinforcement.** I reproduced median occupied lengths of 0.357, 0.242 and 0.480. The edge comparison of 0.359 versus 1.106 still stands for the conservative slab and the unrestricted linear baseline. The false-belief table also reproduces. [5a, §3.4](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md#L329-L357)
* **Reporting and slack.** Labelling the default CLR branch union algorithmic, displaying the projection result, preserving disconnected pieces and printing the base result address the main presentation requests. I accept the authors' correction about sign slack. My earlier report relied on a stale implementation-status line. The revised status table and dedicated tests document per-restriction slack and the sensitivity dial. [Reporting tests](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/tests/testthat/test-batchB.R#L205-L258) · [Slack tests](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/tests/testthat/test-msign-dial.R#L30-L91)
* **N2.** Row-specific person weights repair the binned accounting identity. The revised text limits derivative arguments to settings where the derivative has empirical content and explicitly leaves sharpness under size conditioning open. This resolves my concern for the scope of this release. [2c, §4](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2c_monotone.md#L237-L273)

## Corrections to request before closing

**1. Rewrite the W1 example; spell out pairwise CR beside it.** Guide line 215 and vignette lines 173–176 still describe a positive W1 sign as both A's and everyone else's rates rising with A's share. W1 signs the sum of two person-weighted covariances. Imagine two equally large precincts, 25% and 75% A. A's P rate goes from 40% to 60%, while everyone else's goes from 55% to 45%. The weighted within association is positive (+0.0046875) even though the second group's rate falls. Suggested wording: “Using A’s precinct share as the predictor for both groups, the two within-group associations have a nonnegative weighted sum. Either group’s association may be negative.” For pairwise CR, add that the between term compares A with the named B, while the contextual term pools the original groups, including C. The restored W2 explanation does not tell readers that choosing `pair` changes only one side of that pairing. [Guide §9](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PLAIN_LANGUAGE_GUIDE.md#L213-L218) · [Vignette](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/vignettes/eiPartialID_rxc.Rmd#L168-L192) · [Pairing definition](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md#L99-L103)

**2. Give the coverage figure for the procedure being described.** The guide and vignette attach an 8–10% miss rate to intersecting the estimated reinforcement bracket with a completed JKST confidence interval. The reference code's `CR-naive` procedure allows a subsequent expansion for target uncertainty. That changes the event being counted. I compared them on the same 600 replications of the 150-precinct, positive-edge design.

| Procedure | Covered the true model-implied gap | Miss rate |
|---|---:|---:|
| Literal intersection of the fitted bracket with the completed JKST interval | 304/600 (50.7%) | 49.3% |
| Reference script's `CR-naive` procedure with target expansion | 550/600 (91.7%) | 8.3% |

At this edge the true gap equals the population regression slope. Whenever its estimate falls below truth and remains positive, the entire literal intersection stops below truth, however wide the JKST interval is. This explains why 5a §3.5 now gives asymptotic coverage of at most one half in its boundary example. Keep the warning and attach each numerical result to its own procedure. [Guide §9](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PLAIN_LANGUAGE_GUIDE.md#L223-L235) · [Corrected theoretical argument](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md#L360-L365)

**3. Give both false-belief results.** Guide line 223 and vignette lines 249–256 describe zero coverage of the true gap under a moderate departure. At 600 precincts the conservative slab covered 398/400 at the boundary, 110/400 for the smaller departure ($s=.08$), and 0/400 for the larger departure ($s=.25$). All were nonempty. Report the parameter values and counts. A false belief *can* produce narrow, nonempty intervals that miss the target. Its falsity alone does not imply that every projection will miss. [Simulation table](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md#L344-L357)

**4. Scope the projection guarantee to the fitted model.** The guide describes the projection interval as carrying a proof. 2c explicitly lists inference for signed district models as open under R16. In 5a, §3.1 defaults to pooled M0, whereas the R10 statement also names fixed-district M1/M2. Clarify which result is an abstract coverage theorem conditional on augmented-system assumptions and which district implementation has been shown to meet them. Until that mapping is explicit, label pooled M0 with its rank and inferential conditions, and identify the status of each grouped fit separately. The warning about sharpness for grouped models addresses identification; readers also need the status of their CI. [2c, M5](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2c_monotone.md#L180-L186) · [5a, R10](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md#L229-L251)

A held-out benchmark with pre-specified signs would bear on routine use in elections. It can remain a later validation step for a release labelled experimental. Likewise the drop-one display, integer-table extension, and inference for width changes can remain open if their status stays visible. When several optional restrictions are entered, describe the rung-only comparison as the joint effect of removing them all. A reader seeking one belief's incremental contribution needs the added-last comparison.

## Reproduction of the literal-intersection comparison

Run the following script beside the pinned `verify_R10.py` using NumPy and SciPy. It records the two beta draws already made on each replication, preserving the reference simulator's random stream. It uses the simulator's unit-weighted quadratic fit and the person-weighted slope of its fitted returns for the descriptive bracket, so both procedures use the same slope estimate.

```python
import runpy
import numpy as np

ref = runpy.run_path("verify_R10.py")
x, p = ref["design"](150)
theta, ab = ref["make_truth"](x, p, "boundary, ER far from 0")
parts = ref["linear_parts"](x, p)
truth = parts[3](theta)["D"]
gamma = parts[2]
F = np.column_stack(((1-x)**2, 2*x*(1-x), x**2))
fit = np.linalg.solve(F.T @ F, F.T)
er_weights = (p*(x-p@x)) @ F @ fit / (p @ ((x-p@x)**2))

class CaptureRates:
    def __init__(self, generator):
        self.generator = generator
        self.rates = []
    def beta(self, *args):
        result = self.generator.beta(*args)
        self.rates.append(result)
        return result
    def __getattr__(self, name):
        return getattr(self.generator, name)

rng = CaptureRates(ref["RNG"])
ref["one_rep"].__globals__["RNG"] = rng
literal = expanded = 0
for _ in range(600):
    rng.rates.clear()
    out = ref["one_rep"](x, p, theta, ab, 15, parts)
    b0, b1 = rng.rates
    er = er_weights @ ((1-x)*b0 + x*b1)
    lo, hi = sorted([gamma*er, er])
    base = out["JKST"]
    literal += base is not None and max(lo, base[0]) <= truth <= min(hi, base[1])
    ci = out["CR-naive"]
    expanded += ci is not None and ci[0] <= truth <= ci[1]
print(f"Literal intersection: {literal}/600; reference CR-naive: {expanded}/600")
```
