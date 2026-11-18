# Review brief, Batch B: reviewer 2 (ecological-inference specialist)

## Your role

You are a methodologist who specialises in ecological inference (EI) and its applications to voting:
* the Duncan–Davis method of bounds and Goodman regression;
* King's EI and its R×C extensions;
* the Jiang–King–Schmaltz–Tanner (JKST) partial-identification bounds;
* the recent Monotone EI work of Elzayn et al.

You have also seen many applied EI analyses go wrong.

You have been asked to referee a new module of an R package for EI (the "M-SIGN" module). It lets users impose beliefs about the **sign** of associations, and **contextual reinforcement**, on top of linear-contextual-effects bounds. Read it as a **critical, skeptical referee**.

Your focus is substantive and applied rather than proof by proof; a second reviewer checks the proofs. Your questions are:
* Do these assumptions mean what the documents say they mean, in real EI applications?
* Are they credible, and when?
* Are the mappings to the EI literature faithful?
* Would the package's output be read correctly, or does it invite misuse?

When something seems mathematically wrong, say so, and test it numerically if you can.

You report only. Do not edit the documents or the code.

## What you are reviewing

**Repository:** `https://github.com/mdweaver/eiid`, at **commit `8b98589`**. Check out exactly that commit (`git checkout 8b98589`). If you have no repository access, ask the person who gave you this brief for the files listed below.

**Primary documents (review these in full):**
1. `docs/extension/PHASE2c_monotone.md`: monotone EI in R×C. Associations within and between groups (W1 one-vs-rest, W2, pairwise); exact bounds under sign beliefs (N1); sign beliefs under linear contextual effects (K-rows, Prop. M5); monotone nonparametric contextual effects (N2); the ladder (§5).
2. `docs/extension/PHASE5a_sign_module.md`, **§§1–3 and §5**:
   * Theorem R3, read for what "sharp for the regression" means in practice, and for contextual reinforcement giving one range or two;
   * Prop. R2, read for what "this belief buys nothing / narrows the set but not this quantity / narrows this quantity" means;
   * Theorem R10, read for which intervals users get;
   * §5, the implementation status (what is in the package).

   Read §4 (an earlier internal review) **only after** you have written your own verdicts.
3. **User-facing descriptions of the same module:**
   * `docs/extension/PLAIN_LANGUAGE_GUIDE.md`, **§9 only**;
   * `package/eiPartialID/vignettes/eiPartialID_rxc.Rmd`, the section **"Signs on associations and contextual reinforcement"** only.

**Background (already reviewed; read only what you need to follow the primary documents):**
* `docs/extension/PHASE1b_monotone_did.md`, **§0, §A1, §A2, §A3 and §D only**: the summary of the Elzayn et al. paper, the proposed rungs, and the owner's decisions that PHASE2c implements.
* `docs/extension/PHASE2a_rxc.md`, **§0, §2 (Lemma 1), §3 (Lemma 2) and §6.2 only**: the linear-contextual-effects model, its unidentified part $K$, and the estimand.
* `docs/extension/PLAIN_LANGUAGE_GUIDE.md` **§§1–3**: plain-language versions of the base model and rungs.

If you find a problem in a background document, list it separately as "outside scope".

**Numerical checks** (optional; useful for testing a doubt):
* `docs/extension/verify_phase2c.py`, `docs/extension/verify_phase1_mei_did.py`, `docs/extension/verify_R10.py`, and their `_output.txt` files.
* Run them with `python3 -I <script>` from `docs/extension/` (numpy and scipy are needed).

**External literature you may consult:**
* Elzayn, Goldin, Guage, Ho & Morton, "Monotone Ecological Inference", arXiv:2504.14752 (**v1**; cited as [EGGHM:x]);
* Jiang, King, Schmaltz & Tanner (2019), *Political Analysis*;
* Duncan & Davis (1953); Goodman (1953, 1959);
* King (1997), *A Solution to the Ecological Inference Problem*;
* Cho & Manski (2008); Cross & Manski (2002);
* Wakefield (2004), "Ecological inference for 2×2 tables", *JRSS A*;
* any other published EI literature you judge relevant. Cite what you use.

## Read nothing else

Do **not** read any other file in the repository. In particular, do not read:
* the review logs or plans: `MSIGN_REVIEW_LOG.md`, `REVIEW_RESPONSE_PLAN.md`, `EXTERNAL_REVIEW_PLAN.md`, `NEXT_BUILD.md`;
* `PACKAGE_DESIGN.md`, `NEW_AVENUES.md`, `PHASE2d_did.md`, `PHASE3_inference.md`;
* the other plain-language files;
* anything under `docs/extension/reviews/` other than this brief;
* the R code;
* commit messages or other branches.

They contain other people's conclusions about these documents, and your review must be independent of them.

## Questions to answer

1. **Meaning of the beliefs.**
   * For a typical voter-transition table (rows = previous vote or group; columns = current vote), what does each belief say about voters?
     * $\delta_B$ (between): "group A is more likely than others *in the same unit* to choose P";
     * $\delta_W$ (within, contextual): "rates for P rise with A's share of the unit";
     * W1 vs W2 vs pairwise;
     * contextual reinforcement (CR): $\delta_B$ and $\delta_W$ share a sign, which is not assumed known.
   * Are the documents' interpretations right?
   * When would a careful analyst believe each one, and when not?
     * Racial bloc voting.
     * Partisan transitions with abstention as a category.
     * Referendum support by veteran status.
2. **Faithfulness to Elzayn et al.**
   * Are their definitions, their results ([EGGHM:x], v1 numbering) and their contextual-reinforcement framing represented faithfully?
   * Where PHASE2c generalises them to R×C (W1/W2), is the generalisation natural, or does it change the meaning?
3. **Relation to JKST and the linear-contextual-effects model.**
   * Adding sign rows to the linear model's bounds (Prop. M5, Theorem R3): does this combine two assumptions sensibly?
   * Does the combination smuggle in anything stronger than either alone?
4. **Refutability.**
   * The documents say that with two rows a sign on $\delta_B$, and CR, can never be refuted by aggregate data on their own; only their combination with linear contextual effects can be. Is this correct?
   * Is it communicated so that users do not over-read a non-rejection?
5. **What the user sees.**
   * Two disjoint ranges (`[a, b] U [c, d]`).
   * The "already implied / narrows the set but not this quantity / narrows this quantity" diagnostic.
   * The "descriptive only" Goodman and neighbourhood values.
   * The warning against intersecting the estimated bracket with a confidence interval.

   Will applied researchers read these correctly? What will they get wrong? What should be reported by default?
6. **Gains and their price.** The documents report large narrowing (interval width about 0.36 against 1.1 for linear contextual effects alone) in a simulation where the beliefs hold.
   * Is that simulation design representative of EI data?
   * What happens to coverage when the belief is wrong, and does the module make that risk visible enough? Explore numerically if you can.
7. **Missing pieces.** What would an EI specialist expect here that is absent? For example:
   * ties to King's EI or to Bayesian R×C methods;
   * aggregation over units of very different sizes;
   * turnout as a column;
   * a known-truth benchmark.

## Order of work

1. Read the primary documents and form your own verdicts. Write them down.
2. Only then read PHASE5a §4 (the earlier internal review). Add a short section saying where you agree, disagree, or found something it missed.

## Deliverable

One Markdown file, `docs/extension/reviews/B_ecological_inference_report.md`, with:

1. **Recommendation:** accept / minor revision / major revision, with a one-paragraph overall assessment.
2. **Assessment table:** one row per belief or feature (δ_B sign, δ_W sign, W1/W2/pair, CR, N1 exact bounds, N2, the redundancy diagnostic, two-range output, the descriptive values). Columns:
   * meaning correct? (yes / qualified / no);
   * credible in which applications;
   * risk of misuse;
   * location (document § and line).
3. **Major comments:** substantive errors, misinterpretations, or misleading presentation, each with a suggested fix.
4. **Minor comments.**
5. **Numerical explorations** (if any): what you ran and found. Include the code in an appendix or as files next to the report.
6. **Outside scope:** problems noticed in the background documents.
7. **Comparison with the earlier internal review** (PHASE5a §4).

If you can push to the repository, commit only your report (and any scripts) on a new branch named `review/B-ecological-inference`. Otherwise return the file to the person who gave you this brief.
