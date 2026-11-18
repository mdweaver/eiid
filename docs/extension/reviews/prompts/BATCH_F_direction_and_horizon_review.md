# Review brief, batch F: direction and new avenues

*For an external review agent. Read-only. Review the materials at commit [`27557df`](https://github.com/mdweaver/eiid/tree/27557df) of `mdweaver/eiid` (all paths below are relative to the repository root). This batch is **not** a proof check. Batches A–E (`docs/extension/EXTERNAL_REVIEW_PLAN.md`) already cover correctness, result by result. Here we want a fresh view of **where the project is heading** and **what it has not yet thought of**.*

---

## 1. Your role

Act as three reviewers in one:

1. An **econometrician working on partial identification**: sharp identified sets, random-set theory, moment-inequality inference, sensitivity and breakdown analysis, data combination, and shape restrictions.
2. A specialist in **ecological inference (EI)**: Goodman, Duncan–Davis, King, Wakefield, Greiner–Quinn, Imai–Lu–Strauss, Cho–Manski, Jiang–King–Schmaltz–Tanner (JKST), Elzayn et al. (Monotone EI), and the recent ML / sensitivity work (e.g. McCartan–Kuriwaki).
3. A **mathematician** with broad reach: convex and polyhedral geometry, optimal transport and multi-marginal problems, Fréchet / copula bounds, contingency-table combinatorics, moment problems and sum-of-squares, total positivity and order restrictions, inverse problems and tomography, and stochastic orders.

Say which hat you are wearing when it matters.

## 2. What the project is (a short orientation; the documents are authoritative)

* **Problem.** For units $i$ (polling places, districts), we see row shares $X_i\in\Delta^{R-1}$ (group composition), column totals $T_i\in\Delta^{C-1}$ (outcomes), and sizes $N_i$. The targets are the population cell rates $B_{rc}$ (for example, turnout by group) and contrasts such as gaps between groups.
* **Core framework (Phases 1–3).**
  * The conditional mean of unit rates is $M_H(x)$ with $H=Q+K$. Here $Q$ is identified from the regression $E[T\mid X]$. $K$ is skew-symmetric with $(C-1)\binom R2$ free coordinates $k$, which the data never pin down.
  * Each assumption becomes LP rows in $k$, and bounds are LP optima. "Linear contextual effects" (LCE) is one rung of a ladder that runs from Duncan–Davis (DD) to Goodman. There are district variants (M0, M2, M1δ, M1).
  * The sharp set $\Theta_I$ uses the whole law of $T\mid X$. The LP set $\Theta_O$ uses only the regression.
  * Targets come in two kinds. Model-implied rates $\tilde B$ are affine in $k$. Actual-population rates $B$ need a finite-population allowance.
  * Inference: CLR-type intervals (algorithmic), projection intervals (proved), conservative intervals, an emptiness test, and sensitivity dials with breakdown values. It is conditional on the design (decision D3: no clustering).
* **Monotone EI integration (Phase 2c, Phase 5a).** EGGHM's sign beliefs on between- and within-group associations ($\delta_B$, $\delta_W$, and contextual reinforcement, CR) become rows in $k$, or unions of rows. Under LCE this tightens the JKST bounds.
* **Package.** `package/eiPartialID` (R) ships the core, the M-SIGN module (sign rows, CR, a redundancy diagnostic) and a module registry.
* **Agenda.**
  * `docs/extension/PACKAGE_DESIGN.md`: architecture, module catalogue, layering matrix, two ladders, DiD (§7), proof register R1–R14.
  * `docs/extension/NEW_AVENUES.md`: avenues A–J; §§12–14 respond to outside suggestions.
  * `docs/extension/NEXT_BUILD.md`: the next build choice, M-SHARP (a series-estimated sharp set) or M-LOC (local contextual monotonicity).
* **Structural facts the project treats as settled** (please challenge them if you think they are mis-stated or too narrow):
  * joint ballot-level margins across contests add nothing without a cross-contest assumption;
  * a known aggregate table pins only $(R-1)(C-1)$ of the unknown directions;
  * independent unit-level shocks, however well deconvolved, never reveal $K$;
  * in two groups, the $\delta_B$ sign and CR are never refutable by aggregates.
* **Set aside so far:** label symmetry (G), deconvolution (H), entropy / OT / Choquet as identifying tools (I–J). The reasons are in `NEW_AVENUES.md` §7. We are open to being told this was premature.
* **Recent referee findings that bear on direction** (`docs/extension/reviews/B_partial_identification_report.md`, major comments 1–5). These are not yet propagated to the design documents.
  1. A continuous LP and an integer electorate can disagree about refuting CR.
  2. District intercepts can reverse the pooled within-association sign that M5(c) predicts. The rows stay affine, but sign shortcuts are M0-only.
  3. A local derivative is not identified at isolated support points. This affects N2.3 and any nonparametric local monotonicity.
  4. Full-law sharpness needs conditioning on $(X,N)$, not just $X$. This bears directly on M-SHARP, the next planned build.
  5. Inference guarantees depend on how candidates are constructed and on the asymptotic regime.

## 3. Reading list (in priority order; skim where marked)

1. `docs/extension/PACKAGE_DESIGN.md`: all of it. It is the main object of this review.
2. `docs/extension/NEW_AVENUES.md`: §0 (summary and ranking), §§7–9, §§12–14. Skim §§1–6.
3. `docs/extension/NEXT_BUILD.md`, and `docs/extension/EXTERNAL_REVIEW_PLAN.md` §1 (what is reviewed).
4. `docs/extension/ASSUMPTION_MENU.md` and `docs/extension/PLAIN_LANGUAGE_GUIDE.md`: how assumptions are pitched to applied users.
5. Skim for the model and results:
   * `docs/extension/PHASE2a_rxc.md` §0 and Props. 2, 3 and 7;
   * `docs/extension/PHASE2c_monotone.md` §0 and M1–M5;
   * `docs/extension/PHASE5a_sign_module.md` §1.2, §2.2 and §3.2 (statements only);
   * `docs/extension/PHASE3_inference.md` §0;
   * `docs/extension/PHASE2d_did.md` §0 and §10.
6. `docs/extension/reviews/B_partial_identification_report.md` §§1–2.

You do not need to check proofs or run code. Do quote section numbers when you rely on a stated result.

## 4. Part I: problems with, and disagreements about, the proposed direction

We want judgements about direction, not line edits. For each point, say what you would do instead. Questions to consider (you need not answer all, and you may raise others):

* **Target and estimand.**
  * Is bounding the model-implied rate $\tilde B$ (with an allowance for $B$) the right primary estimand for users?
  * Should the package lead with the actual-population rate, or with a superpopulation parameter? (Note decision D3: we condition on the design and on district effects.)
* **The LCE-centred architecture.** Everything is an LP in the skew part $k$ of one parametric family.
  * Is anchoring on LCE a strength (tractability, interpretable dials) or a trap (every module inherits LCE)?
  * What would a genuinely nonparametric core look like? Would it be estimable?
* **The two ladders.**
  * Are the plausibility tiers and the nesting and non-nesting claims ordered sensibly from a user's point of view?
  * Are some rungs ones that no applied researcher would believe, or are credible rungs missing?
* **Monotone EI integration.**
  * Is turning EGGHM's sign beliefs into restrictions on a model-implied mean the right reading of what their users believe?
  * Is CR's unrefutability in 2×2 a reason to make it less prominent, or more?
* **Information use.**
  * Are we leaving identifying information in $P_{T\mid X,N}$ unused? The review found that conditioning on $N$ matters.
  * Is the M-SHARP plan (series regression of unit-level DD limits) the right way to use it?
* **Inference choices.**
  * CLR is the algorithmic default, with projection alongside. We condition on the design. We give finite-sample bounds that are valid but uninformative for many small units.
  * Is this a defensible stance? What would you change first?
* **Discreteness.** Electorates are integers, and the continuous relaxation can disagree with integer feasibility (report §2.1). Should integrality be a module, a diagnostic, or ignored at realistic $N$?
* **Build order and priorities** (`PACKAGE_DESIGN.md` §6, `NEXT_BUILD.md`, `NEW_AVENUES.md` §13.8).
  * What would you cut, promote or demote?
  * M-SHARP or M-LOC next?
* **DiD extension** (`PACKAGE_DESIGN.md` §7). Is mapping each module onto changes in unit rates coherent with modern DiD practice (parallel trends, Rambachan–Roth-style relaxations)?
* **Communication.** Do the plain-language framings in `PLAIN_LANGUAGE_GUIDE.md` overclaim anywhere?

Rate each disagreement **fundamental** (changes the direction), **significant** (changes priorities or framing) or **minor**.

## 5. Part II: avenues not yet explored

The project's method so far has been to find **observable information the current bounds discard**, or **weak, plausible, checkable assumptions**, and to turn either into restrictions on the unidentified part. Please bring fresh eyes to this. We want ideas that are new relative to `NEW_AVENUES.md` and `PACKAGE_DESIGN.md`, or old ideas that we dismissed too quickly.

For each idea give:

1. **Idea:** one paragraph, with the mathematical object involved.
2. **Information or assumption used:** which observable it exploits (e.g. $P_{T\mid X,N}$, unit sizes, spatial or hierarchical structure, multiple contests or periods, auxiliary surveys or tables, turnout as a second outcome, integrality), or what it assumes. Say how plausible the assumption is and whether it can be checked.
3. **Expected identifying power:** none / small / medium / large / point. Give your reasoning, and a toy example if you can.
4. **Fit with the framework:** do we get LP rows in $k$, a union, a convex non-polyhedral set, or an integer program? Does it compose with existing modules?
5. **Estimability and inference:** sample needs and known inference tools.
6. **Key references:** verified, with author, year and venue. Mark any reference you are not certain of as *unverified*.
7. **Priority:** high / medium / low, with one line of justification.

Directions we have **not** worked through, as prompts only. Do not let this list limit you.

* **Partial identification:**
  * monotone instrumental variables and monotone treatment response analogues for composition (Manski–Pepper);
  * a rigorous random-set / support-function computation of the sharp set (Beresteanu–Molchanov–Molinari), including core-determining classes;
  * robust-Bayes or prior-free posteriors for set-identified models (Giacomini–Kitagawa; Kline–Tamer);
  * breakdown frontiers in place of scalar breakdown values (Masten–Poirier);
  * data combination (e.g. Fan–Sherman–Shum; D'Haultfœuille et al.);
  * testable shape restrictions (Chetverikov–Santos–Shaikh);
  * many-moment-inequality inference.
* **Ecological inference:**
  * hybrid designs that combine aggregates with individual surveys (Glynn–Wakefield and later work), and survey-calibrated EI;
  * hierarchical or spatial borrowing as an *identifying* assumption, not just a smoothing device;
  * the Cho–Manski tomography view;
  * turnout and vote choice jointly (nested tables);
  * using several elections or offices in the same units as replicates.
* **Mathematics:**
  * multi-marginal and martingale optimal transport, and Fréchet–Hoeffding / copula bounds with extra structure;
  * transportation polytopes, Markov bases and integer-table counting (Diaconis–Sturmfels; Barvinok) for the integer electorate;
  * moment problems and Lasserre / sum-of-squares hierarchies for polynomial rungs (cf. the Bernstein ladder, avenue D);
  * total positivity and variation-diminishing properties for monotone or likelihood-ratio orders across groups;
  * stochastic orders on unit-rate distributions;
  * isotonic and convex regression;
  * inverse problems and tomography (Radon-type views of the $X\to T$ map), sparsity or low-rank priors on $K$, and their identification analogues;
  * information geometry beyond the max-entropy selection we set aside.

Please also say plainly which **set-aside avenues** (G–J) or **low-ranked avenues** you would revive, and why.

## 6. Output format

Write a single Markdown report, `F_direction_review_<date>.md`, with:

1. **Summary** (≤ 1 page). Your overall view of the direction, the three most important disagreements and the three most promising new avenues.
2. **Part I table:** disagreement, severity, where it bites (doc and section), and the recommended change.
3. **Part I discussion:** one subsection per fundamental or significant item.
4. **Part II:** one entry per idea in the 7-field format above, ordered by priority.
5. **Recommended build order:** your version of `PACKAGE_DESIGN.md` §6 / `NEXT_BUILD.md`, with the reasons for each change.
6. **References:** verified ones separated from *unverified* ones.

**Norms.**
* Be direct. If you think a direction is not worth pursuing, say so and explain why.
* Separate what you know from what you conjecture.
* Proof sketches and small counterexamples are welcome but not required.
* Do not repeat the batch B referee's points. Build on them where they matter for direction.
* Do not edit the repository.
