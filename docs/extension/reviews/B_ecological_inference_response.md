# Response to the referee: Batch B, ecological inference

**Reviewed version:** commit [`8b98589`](https://github.com/mdweaver/eiid/tree/8b98589653db71849fc7cdca500312754c26e95d).
**Revised version:** commit [`e170be9`](https://github.com/mdweaver/eiid/tree/e170be9161b9bbc250a17a27ffd39bc523c1c88b) on branch `claude/wizardly-volta-lbh4bm`.
**All changes since the reviewed version:** [compare `8b98589…e170be9`](https://github.com/mdweaver/eiid/compare/8b98589653db71849fc7cdca500312754c26e95d...e170be9161b9bbc250a17a27ffd39bc523c1c88b).

**2c** is `docs/extension/PHASE2c_monotone.md`; **5a** is `PHASE5a_sign_module.md`; **Guide** is `PLAIN_LANGUAGE_GUIDE.md`; **Vig** is `package/eiPartialID/vignettes/eiPartialID_rxc.Rmd`. Line numbers in quoted comments refer to the reviewed version.

---

## 1. Summary

Thank you for an applied reading that changed what the package reports and how it explains itself.

* **Reproduction.** Your appendix reproduces exactly: the union-length table, the W1/W2 sign reversal and the false-CR table (399/336/2; 398/110/0). It is saved as `reviews/referee_checks_B_ecological_inference.py`.
* **Recommendation.** We agree with major revision. The module stays in the package, **marked experimental**, with revised claims, examples and default reporting.
* **Overlap.** A partial-identification referee reviewed the same material at the same time. Where your points overlap (integer tables, M3(d), refutation, binning, size conditioning, the union lengths), the corrections are shared and described here as well.

**One partial dispute.** You note (5a line 385) that the slack on beliefs was outside the implementation. It is implemented (`sensitivity(dial = "slack")`, with a `slack` element in each restriction). Our implementation-status section was stale and said otherwise; it is corrected (5a §5).

**What changed for users** (release 0.3.0 of `eiPartialID`):
* **W2 is described correctly:** one pooled sum over the original groups, not a direction for each group. Your two-unit example is now a check.
* **Refutation is stated with its base model,** and the output says whether the linear rung alone is rejected next to any joint rejection.
* **The default reinforcement interval is labelled algorithmic** on its own line. A projection interval, which has a proof, is printed beside it. Quoted simulation figures name their method.
* **Every union piece is printed.**
* **A rung-only (base) column sits beside every restricted result.**
* **Weight concentration** is reported for each signed association.
* **The guide and the vignette now cover:**
  * writing each belief down before looking at the result;
  * what a sign does not mean;
  * columns and denominators, including turnout;
  * how to read sets, intervals, gaps and hulls;
  * the relation to King's EI and hierarchical R×C methods;
  * your false-belief result.

---

## 2. What changed

### 2.1 Theory

| Document | Changes relevant to your report |
|---|---|
| [2c](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE2c_monotone.md) | <ul><li>**Definition 1 remarks** (new): whose behaviour is averaged and with which weights; W2 as a pooled sum; δ_B weights $p_ix_i(1-x_i)$; W2, pairwise and their CR marked as extensions of EGGHM; $g-B=M^{-1}\Delta$ (signing $\Delta_A$ does not sign the Goodman bias).</li><li>**M2** sharp for continuous unit tables, with **Remark M2.1** giving your two-precinct example and the partial-identification referee's.</li><li>**M3** hypotheses; (d) with $\gamma=1$ and $ER=0$ (the singleton is the gap's image).</li><li>**Reading** gives your equal-rates witness for R×C and the three distinct emptiness events.</li><li>**M5** per rung: under district regressions a sign is one pooled restriction.</li><li>**N2**: support convention, size conditioning (Conjecture N2.2′), the derivative not identified on a finite support (N2.3), person-weighted bins (Lemma N2.4).</li><li>§5 ladder: "both signs" means both weakly positive or both weakly negative.</li></ul> |
| [5a](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE5a_sign_module.md) | <ul><li>R3 scoped to the pooled model; R3(e) "rejects the conjunction; does not attribute".</li><li>R2(e) "zero incremental gain for this target"; §2.5 width gains need joint endpoint inference.</li><li>R10(d): conditioning on non-emptiness gives no protection against choosing a belief for its result.</li><li>§3.4: union column = occupied length, hull span, unmerged sum and overlap counts (your table A); new false-CR table (your table B).</li><li>§3.5: quoted numbers name their method.</li><li>§5: release 0.3.0 table, slack line corrected.</li></ul> |
| [1b](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PHASE1b_monotone_did.md) | Your four outside-scope items (see §4). |

### 2.2 Code (`eiPartialID` 0.3.0)

* **Status.** `assoc` and `cr` restrictions are **experimental**, and `print()` shows the flag.
* **Grouped rungs.** A fit-time note says each sign is one pooled restriction across groups.
* **W2** is labelled "pooled over groups" in the metadata, notes, print and help.
* **Refutation messages.**
  * They name the base model.
  * They say whether the rung alone is rejected.
  * An empty plug-in is called a sample diagnostic, not a test.
* **Reinforcement intervals.** The CLR branch union stays the default, labelled **algorithmic** on its line, with the projection branch union printed beside it.
* **Union pieces.** All of them are printed, including three or more.
* **Base columns.** New `base_*` columns hold the rung-only estimate and intervals, printed beside restricted results.
* **Weight concentration.** `eiRedundancy()$concentration` gives, per signed association, the largest unit weight, its unit, and the effective number of units $1/\sum u_i^2$.
* **`eiRedundancy()` verdicts.** These read:
  * "empty at the point estimate (see the emptiness test)";
  * "cuts the set, not this target", which means zero gain for this target;
  * a header naming the rung and explaining the first and last orderings.
* **Tests:** [`tests/testthat/test-batchB.R`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/tests/testthat/test-batchB.R). **Changes and known limits:** [`NEWS.md`](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/NEWS.md).

### 2.3 Documentation

* **[Guide](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/docs/extension/PLAIN_LANGUAGE_GUIDE.md)**
  * **§2:** which requirements each rung actually uses.
  * **§9, rewritten.**
    * Status: experimental.
    * The δ_B weighting; the W1 predictor (A's share in both parts); W2 pooled.
    * Refutation with its base model.
    * Comparison with the rung alone.
    * Gains and the cost of a wrong belief.
    * Separate ranges; "do not do this by hand".
    * Redundancy verdicts.
    * **New:** write each belief down first; what a sign does not mean; columns and denominators; reading the printout; other methods.
* **[Vignette](https://github.com/mdweaver/eiid/blob/e170be9161b9bbc250a17a27ffd39bc523c1c88b/package/eiPartialID/vignettes/eiPartialID_rxc.Rmd)**, section "Signs on associations and contextual reinforcement".
  * **Corrections:**
    * the W1 predictor;
    * the δ_B weights;
    * the slab written as $\operatorname{hull}\{\gamma ER, ER\}$ (your line-163 point);
    * "exact" removed;
    * the 0.36/99% figures attributed to the conservative method;
    * the false-CR price added.
  * **New subsection "Reading the output":** sets versus intervals; gaps and hulls; writing beliefs down first; what a sign does not mean; columns, orientation and turnout; other methods and benchmarks.
  * **Descriptive values:** named as the weighted aggregate regression slope and the neighbourhood-model gap.

### 2.4 Verification scripts

* `verify_phase2c.py`: new checks, among them:
  * **C14**, your W1/W2 sign reversal;
  * **C7–C7b**, the integer examples;
  * **C10–C10b**, the grid and size examples;
  * **C12**, binning;
  * **C13**, the smooth example.
* `verify_R10.py`:
  * union lengths merged before measuring, reproducing your table A (0.357 / 0.242 / 0.480, with 134 / 591 / 0 overlapping runs);
  * **V12** reproduces your table B exactly.

---

## 3. Assessment table: your rows, with what changed

| Belief or feature | Your concern | What changed |
|---|---|---|
| Between sign, δ_B | A weighted average, not an ordering everywhere or against each rival | **[Theory]** 2c Definition 1 remarks (weights; pure units get zero weight). **[Docs]** Guide §9, Vig: weighting stated; "what a sign does not mean". **[Code]** weight concentration reported. |
| Within sign, δ_W (W1) | Predictor is A's share in both parts; either covariance can be negative; causal reading needs a design | **[Docs]** Guide §9 and Vig corrected ("rates of A, and of everyone else, rise on average with A's share"); causal and sorting caveat added. |
| W2 | A pooled sum, not per group; W1 and W2 can disagree | **[Theory]** 2c Definition 1. **[Code]** "pooled over groups" labels. **[Docs]** Guide §9, Vig, 1b §D.1. **[Check]** C14. A per-group version is planned. |
| Pairwise CR | Print the two terms separately | **Planned (0.3.x):** the CR table shows branch intervals but not yet each term. **[Docs]** The guide and vignette explain the pairing. |
| CR | Belief selection; a wrong belief gives narrow, nonempty, wrong intervals | **[Theory]** 5a §3.4 false-CR table; R10(d) belief-selection sentence. **[Docs]** Guide §9 and Vig state the cost and recommend the base column and pre-specification. **[Code]** base columns. |
| N1 "exact" bounds | Fails with counted ballots | **[Theory]** 2c M2 continuous label, Remark M2.1 (your example). |
| N2 | Support, size conditioning, derivatives, binning | **[Theory]** 2c §4: support convention, Conjecture N2.2′, N2.3 not an observable restriction, Lemma N2.4. Local rows dropped from the build. |
| Redundancy diagnostic | Label base, target and ordering; zero last-added gain is not logical redundancy | **[Theory]** 5a R2(e), §2.5. **[Code]** verdict wording and header. |
| Two-range output | Keep every component; separate set and CI | **[Code]** every piece printed. **[Docs]** "Reading the output". |
| Descriptive values | Name them; arithmetic check only | **[Docs]** Vig and Guide name them; the `voter_vax` agreement is called an arithmetic check. Renaming the table columns in the printout is planned. |

---

## 4. Point-by-point

Tags: **[Theory]**, **[Code]**, **[Docs]**, **[Check]**.

### Major comment 1: state whose behaviour is averaged, and with which weights

> *δ_B weights $p_ix_i(1-x_i)$ and ignores pure precincts; W1's predictor is A's share in both covariances; W2 signs one pooled sum (the two-unit example gives W1 −0.02633, W2 +0.00110); signing $\Delta_A$ does not sign the Goodman bias; print the two pairwise-CR terms.*

**Response.** Accepted.
* **[Theory]** New remarks under 2c Definition 1 state the weights, the W2 pooled sum and $g-B=M^{-1}\Delta$. 1b §D.1 now describes W2 as pooled.
* **[Docs]**
  * Guide §9 and the vignette: "rates of A, and of everyone else, rise on average with A's share"; W2 is "one pooled sum"; a direction for each named group needs separate restrictions.
* **[Code]** W2 is labelled "pooled over groups" in metadata, notes, print and help.
* **[Check]** C14 reproduces your example.
* **Not yet done:**
  * printing the two pairwise-CR terms separately (planned for 0.3.x);
  * per-group W2 restrictions (planned);
  * placing the two-unit example beside the selector in the vignette (it is cited from the guide and checked in C14).

### Major comment 2: an argument for a between sign rarely supplies one for reinforcement

> *Pre-specify each restriction (row, column, population, geography, pairing, reason); show results without the belief and under the alternative; under district regressions a sign is one pooled restriction.*

**Response.** Accepted.
* **[Docs]** Guide §9 and the vignette now open the advice with "**Write each belief down first**": row, column, population and denominator, geography, comparison, the independent reason, and how well that evidence transfers. Report the result without the belief and, if contested, with the opposite sign.
* **[Code]**
  * The base (rung-only) columns are printed beside every restricted result, so the price of a belief is always visible.
  * Under grouped rungs, a fit-time note says each sign is one pooled restriction across groups.
* **[Theory]** 2c M5 per rung.

### Major comment 3: distinguish three kinds of sharpness and two targets

> *EGGHM mappings are faithful for W1; W2, pairwise and their CR are extensions. JKST plus a sign is a coherent intersection; an implied sign adds nothing. R3 gives an estimated bound for the model-implied gap using the regression and design, which can be wider than the full-law sharp set. M2 fails with counted ballots.*

**Response.** Accepted.
* **[Theory]**
  * W2, pairwise and their CR are marked as extensions of EGGHM (2c Definition 1 remarks).
  * R3 is described as sharp for the information in the regression and the design. Remark R3.1 gives the full-law-versus-regression example, which is the one you describe.
  * M2 is labelled continuous, and Remark M2.1 gives your two-precinct example.
* **[Docs]** The guide says what the bound is: an estimate, from the fitted regression and the observed compositions, of the values the model allows, sharp only relative to that information. It also says what a sign does not mean, including that an implied sign adds nothing (see the "already implied" table).

### Major comment 4: say exactly what a refutation establishes

> *With margins and continuous rates the non-refutability is sound, but known cells, zeros or integer counts can also rule out the witness; distinguish population emptiness, an empty plug-in and a test rejection; show the base result; M3(d) needs a mixed unit.*

**Response.** Accepted.
* **[Theory]** 2c "Reading" gives your equal-rates witness for R×C and lists the three emptiness events. 5a R3(e) says that population emptiness "rejects the conjunction; does not attribute". M3(d) covers $\gamma=1$ and $ER=0$.
* **[Code] [Docs]**
  * The non-refutability wording now names the base model: "always feasible only in the model-free problem with continuous tables; with the rung, known cells, integer counts or other restrictions it can be refuted".
  * The output says whether the rung alone is rejected next to a joint rejection, and calls an empty plug-in a sample diagnostic.
* **[Check]** C8, C8b.

### Major comment 5: report disconnected sets and uncertainty in a usable form

> *Fix the slab ordering in Vig line 163; several CRs can give more than two pieces; print set and CI separately, keep every component; the default CLR CR interval is conjectural, so say so on the line and show the conservative result; the 0.359/0.992 figures are for the conservative method; name the descriptive values; label redundancy verdicts with base, target and ordering; report occupied length and hull span; show DD, the base LCE result and the restricted result together.*

**Response.** Accepted, except two printout items that are planned.
* **[Code]**
  * All union pieces are printed.
  * The CLR line says "algorithmic", with the projection union beside it.
  * Base columns are shown.
  * Weight concentration is reported.
  * Redundancy verdicts carry a header with the rung and orderings.
* **[Docs]**
  * The vignette writes the slab as $\operatorname{hull}\{\gamma ER, ER\}$.
  * The guide and vignette attribute 0.36/99% to the conservative method.
  * "Reading the output" covers set versus CI, gaps, probability weights and hulls.
  * The descriptive values are named as the weighted aggregate regression slope and the neighbourhood-model gap, and the `voter_vax` agreement is called an arithmetic check.
* **[Theory]** 5a §3.4 reports occupied length and hull span; §3.5 numbers name their method.
* **Planned:**
  * renaming the descriptive table's columns in the printout;
  * the L5/L6 sign check at the fitted point;
  * a confidence interval for width gains (5a §2.5; this needs joint endpoint inference).

### Major comment 6: correct the simulation table and show the price of the narrower result

> *The union column double-counts overlapping branches (0.242, not 0.430, at the near-zero design); the design is narrow; a matched false-CR experiment shows coverage 398, 110, 0 of 400 while always nonempty; show a misspecification comparison and a drop-one display.*

**Response.** Accepted.
* **[Theory]** 5a §3.4 replaces the union column with your three statistics and adds the false-CR table, with a "does not test" list (fixed designs, continuous rates, no segregation or spatial dependence).
* **[Check]** `verify_R10.py` merges intervals before measuring, and V12 reproduces your table B.
* **[Docs]** The guide and vignette state the cost of a wrong belief and recommend showing the base result.
* **Planned:** the drop-one display (each belief removed in turn).

### Major comment 7: N2 requires an explicit support and conditioning argument

> *Finite support versus a grid; conditioning on $(X,N)$; the derivative is not identified on isolated support; binned means need the row's own population weights.*

**Response.** Accepted.
* **[Theory]**
  * A support convention in 2c §4.
  * The sharp version is recorded as **Conjecture N2.2′**, with the size example as Remark N2.2.1. We leave it open for now, because a separate review of the rungs bears on it.
  * N2.3 is not an observable restriction on a finite support, and the local rows are dropped from the build.
  * **Lemma N2.4** proves the person-weighted bin repair.
* **[Check]** C10, C10b, C11, C12, C12b.

### Minor comments

1. **Turnout example.** **[Docs]** Guide §9 and the vignette explain that with abstention as a column the eligible population is the denominator, and that loyalty among voters can coexist with a negative eligible-population gap. A fully worked numerical example is not yet included.
2. **Weight concentration and the residual claim.**
   * **[Code]** The largest weight and the effective number of units are reported per signed association.
   * **[Theory]** The residual claim in 2c (Remark M5.1) now states its conditions.
   * **Planned:** the leave-largest-unit-out display, with the drop-one display.
3. **King's EI and hierarchical R×C; a held-out benchmark.**
   * **[Docs]** A paragraph in the guide and vignette: a posterior mean inside a sign bound is not independent support.
   * **Planned:** the held-out benchmark with locked sign choices, as part of the package's real-data step (`PACKAGE_DESIGN.md` S11).
4. **2×3 versus 3×2; adding-up; complementary outcomes; "both signs".**
   * **[Docs]** Added to the guide and vignette.
   * **[Theory]** 2c §5 defines "both signs" as both weakly positive or both weakly negative.
   * An automatic check for contradictory signs on complementary columns is not implemented.
5. **One status key; the planned API versus the shipped syntax.**
   * **[Theory]** 5a §5 has one status table for 0.3.0, and 2c's header describes the shipped `restrictions` syntax.
   * **[Docs]** The guide says what is experimental.

### Outside scope (background documents)

All four are fixed in 1b. **[Theory]**
* **1b §A2.** "CR needs two LPs" is marked superseded (one slab for W1 and 2×2). "Different moments" is no longer offered as proof of non-nesting.
* **1b §D.1 and §A1.**
  * W2 is described as a pooled sum.
  * Table 1 is attributed to Theorem 1, and Tables 2–3 to Theorems 2–3.
  * The vaccination evidence now carries the transport and selection caveats.
* **Guide §2.** The table says which requirements each rung uses: Duncan–Davis needs only the accounting identity, and N1 needs no regression.
* **1b §D.2.** A reversal in the aggregate yes-share is not evidence that a group's rate reverses, and a monotone aggregate does not show that each group's rate is monotone. The degree of the fitted curve must be prespecified.

### Comparison with the earlier internal review

Agreed. 5a §4 now records what the internal review missed (the union length, counted ballots, the W1/W2 reversal, false-CR coverage). R10(d) now says that conditioning on a non-empty interval gives no protection against choosing a belief because it gives an attractive result.

---

## 5. Not yet done

**Planned for the package (0.3.x):**
* the two pairwise-CR terms printed separately;
* renaming the descriptive table's columns;
* the drop-one and leave-largest-unit-out displays;
* the L5/L6 sign check at the fitted point;
* per-group W2 restrictions;
* inference on width gains.

**Documentation:** a fully worked turnout example.

**Research:**
* Conjecture N2.2′;
* integer unit tables beyond 2×2 margins (R15);
* sharpness of the sign restrictions under district rungs (R16);
* the held-out benchmark with known truth.
