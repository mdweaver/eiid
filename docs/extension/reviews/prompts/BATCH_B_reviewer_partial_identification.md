# Review brief, Batch B: reviewer 1 (econometrician, partial identification)

## Your role

You are an econometrician who works on partial identification and on inference for moment inequalities and intersection bounds (Manski; Chernozhukov–Lee–Rosen; Beresteanu–Molchanov–Molinari; Fang–Santos).

You have been asked to referee the theory behind a new module of an R package for ecological inference (EI). Read it as a **critical, skeptical referee**:
* Assume nothing is right until you have checked it.
* "Checked numerically" is evidence, not proof.
* A justification tag is a claim to be verified, not a reason to believe the step.
* Where you suspect a gap, try to build a counterexample, analytically or numerically.

You report only. Do not edit the documents or the code.

## What you are reviewing

**Repository:** https://github.com/mdweaver/eiid, at **commit `8b98589`** (full hash `8b98589653db71849fc7cdca500312754c26e95d`). Everything you need is on GitHub; you do not need any other environment. If the repository is not visible to you, ask the person who gave you this brief for access or for copies of the files below.

**Direct links (GitHub, pinned to the reviewed commit; open these in a browser or fetch the raw file):**

| File | What to read |
|---|---|
| [`docs/extension/PHASE2c_monotone.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE2c_monotone.md) | primary: all |
| [`docs/extension/PHASE5a_sign_module.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE5a_sign_module.md) | primary: §§1–3; §4 only after your verdicts; §5 optional |
| [`docs/extension/verify_phase2c.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_phase2c.py) | check script |
| [`docs/extension/verify_phase2c_output.txt`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_phase2c_output.txt) | its output |
| [`docs/extension/verify_R2.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R2.py) | check script |
| [`docs/extension/verify_R2_output.txt`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R2_output.txt) | its output |
| [`docs/extension/verify_R3.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R3.py) | check script |
| [`docs/extension/verify_R3_output.txt`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R3_output.txt) | its output |
| [`docs/extension/verify_R10.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R10.py) | check script |
| [`docs/extension/verify_R10_output.txt`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_R10_output.txt) | its output |
| [`docs/extension/verify_package_design.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_package_design.py) | helper imported by verify_R2/R3 (as needed) |
| [`docs/extension/verify_phase1_mei_did.py`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_phase1_mei_did.py) | the '1b V1–V8' checks |
| [`docs/extension/verify_phase1_mei_did_output.txt`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/verify_phase1_mei_did_output.txt) | its output |
| [`docs/extension/PHASE2a_rxc.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE2a_rxc.md) | background: sections named below only |
| [`docs/extension/PHASE2b_districts.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE2b_districts.md) | background: §0 and Theorem 1′ only |
| [`docs/extension/PHASE3_inference.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE3_inference.md) | background: §0, §2, §5 only |
| [`docs/extension/PHASE1b_monotone_did.md`](https://github.com/mdweaver/eiid/blob/8b98589653db71849fc7cdca500312754c26e95d/docs/extension/PHASE1b_monotone_did.md) | background: §0, §A1, §A2, §D only |

To work locally instead: `git clone https://github.com/mdweaver/eiid.git && cd eiid && git checkout 8b98589653db71849fc7cdca500312754c26e95d`. Browse the whole commit at https://github.com/mdweaver/eiid/tree/8b98589653db71849fc7cdca500312754c26e95d (but read only the files listed here).


**Primary documents (review these in full):**
1. `docs/extension/PHASE2c_monotone.md`: monotone EI in R×C. Lemma M1; Propositions M2, M3, M4, M5; Lemma N2.1; Propositions N2.2, N2.3; §5 ladder; §6.
2. `docs/extension/PHASE5a_sign_module.md`, **§§1–3 only**: Theorem R3 (sharpness of the linear-contextual-effects model plus restrictions on $k$; contextual reinforcement as a union or a slab), Proposition R2 (when an added restriction tightens the bound: the face test), and Theorem R10 (inference with added rows, the parallel disjunction, the union of branch confidence intervals, emptiness, and the counterexample in (e)). This includes the status notes inside these sections.
   * §4 (an earlier internal review) **only after** you have written your own verdicts (see "Order of work").
   * §5 (implementation status) is optional context.

**Numerical checks that accompany them** (read them, run them, and extend them if useful):
* `docs/extension/verify_phase2c.py` and `verify_phase2c_output.txt`
* `docs/extension/verify_R2.py`, `verify_R3.py`, `verify_R10.py` and their `_output.txt` files
* `docs/extension/verify_package_design.py`: a helper module that `verify_R2.py` and `verify_R3.py` import. Read it only as needed.
* `docs/extension/verify_phase1_mei_did.py` and its output: the "1b V1–V8" checks cited in PHASE2c.

Run them with `python3 -I <script>` from `docs/extension/` (numpy and scipy are needed). `verify_R3.py` takes about 12 minutes.

**Background (already externally reviewed; read only the parts the primary documents cite):**
* `docs/extension/PHASE2a_rxc.md`:
  * §0 (notation, tags, citation table §0.4);
  * Lemmas 0–3;
  * Theorem 1;
  * Propositions 1, 2 and 7, with Prop. 7′;
  * Remarks 5.1–5.2.
* `docs/extension/PHASE2b_districts.md`: §0 (the rungs) and Theorem 1′.
* `docs/extension/PHASE3_inference.md`:
  * §0 (sampling, targets, tags, §0.4 source of randomness, §0.5 regime table);
  * §2 (Lemma P3.1);
  * §5: P3.5(a), (b1), (b2), Conjecture (b3), P3.6.

  PHASE5a's Theorem R10 is built on these.
* `docs/extension/PHASE1b_monotone_did.md`, **§0, §A1, §A2 and §D only**. These give PHASE2c's summary of the paper it follows, its theorem numbering, and the owner's decisions.

If you find a problem in a background document, list it separately as "outside scope". Do not review those documents otherwise.

**External literature you may consult:**
* Elzayn, Goldin, Guage, Ho & Morton, "Monotone Ecological Inference", arXiv:2504.14752 (**v1**; PHASE2c cites its numbering as [EGGHM:x]);
* Jiang, King, Schmaltz & Tanner (2019), "Ecological Regression with Partial Identification", *Political Analysis*;
* Chernozhukov, Lee & Rosen (2013), *Econometrica* 81:667–737;
* Fang & Santos (2019), *Review of Economic Studies* 86:377–412;
* Cross & Manski (2002), *Econometrica* 70:357–368;
* standard convex-analysis and LP references.

## Read nothing else

Do **not** read any other file in the repository. In particular, do not read:
* the review logs or plans: `MSIGN_REVIEW_LOG.md`, `REVIEW_RESPONSE_PLAN.md`, `EXTERNAL_REVIEW_PLAN.md`, `NEXT_BUILD.md`;
* `PACKAGE_DESIGN.md`, `NEW_AVENUES.md`, `PHASE2d_did.md`;
* the plain-language guides;
* anything under `docs/extension/reviews/` other than this brief;
* commit messages or other branches.

They contain other people's conclusions about these documents, and your review must be independent of them. Do not read the R package code either; a separate code review covers it.

## Questions to answer (in addition to checking every result)

1. **Theorem R3.**
   * Is the sharpness claim (a) correct as stated, and is the claim "sharp for the regression" precisely defined?
   * In (c), is the contextual-reinforcement set exactly the union of the two branch images? Does the stated condition for the hull to be sharp hold?
   * Are the hypotheses sufficient, in particular the dimension condition on $A$, admissibility of $C$, and $P(X\in A)=1$?
2. **Proposition R2.**
   * Is the equivalence "an added restriction moves an endpoint iff it removes the whole optimal face" correct, including for non-singleton optimal faces?
   * Is the claim that a solver's dual alone does not certify correct?
3. **Theorem R10.**
   * (a) Added rows join the critical-value maximum; does coverage really carry over?
   * (b) Is the "parallel disjunction" (inflate each piece, then take the maximum) valid for the conservative interval?
   * (c) Is the union of branch intervals, each at the full level $\alpha$, valid without a Bonferroni split? Is the true branch $j^*$ non-random, given how the branch rows are built?
   * (d) Is the emptiness bound right for both versions?
   * (e) Does the counterexample show persistent under-coverage of the naive bracket?
   * Which statements remain valid given P3's split of the CLR result into (b1) theorem / (b2) proved / (b3) conjecture? Do the status notes in §3.2 get this right?
4. **PHASE2c.**
   * Are Lemma M1's identities and the W1/W2 associations correctly derived from the accounting identity?
   * Is M2's exactness claim right?
   * Are M3's closed forms and refutability statements correct?
   * Does M5 correctly turn sign beliefs into linear rows in $k$ under every linear-contextual-effects rung?
   * Are N2.1–N2.3 correct, and is N2's status (outer or sharp) stated accurately?
5. **Mapping to EGGHM.** Are the cited EGGHM results ([EGGHM:x], v1 numbering) used correctly, with their conditions met?
6. **Numerical evidence.** Do the scripts test what the documents say they test? Are there adversarial cases they miss? Run at least one check of your own design for each major result.

## Order of work

1. Read the primary documents and form your own verdicts. Write them down.
2. Only then read PHASE5a §4 (the earlier internal review). Add a short section saying where you agree, disagree, or found something it missed.

## Deliverable

One Markdown file, `docs/extension/reviews/B_partial_identification_report.md`, with:

1. **Recommendation:** accept / minor revision / major revision, with a one-paragraph overall assessment.
2. **Verdict table:** one row per result (M1, M2, M3, M4, M5, N2.1, N2.2, N2.3, R3(a)–(e), R2 parts, R10(a)–(e)). Columns:
   * the result;
   * your verdict: *established*, *established under additional or stated conditions* (say which), *incorrect*, or *not established*;
   * severity, if not established;
   * the location (document § and line).
3. **Major comments:** errors in mathematical claims, each with a counterexample or the failing step, and a suggested correction.
4. **Minor comments:** exposition, notation, precision of statements.
5. **Numerical checks run:** what each check tests and its result. Include the code in an appendix, or as files you add next to the report.
6. **Outside scope:** problems noticed in the background documents.
7. **Comparison with the earlier internal review** (PHASE5a §4).

If you can push to the GitHub repository, commit only your report (and any check scripts) on a new branch named `review/B-partial-identification`. Otherwise return the file to the person who gave you this brief.
