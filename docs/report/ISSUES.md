# Issues found in the source documents while compiling the report

Found while compiling the report against `b270285`, recorded here (the report's brief forbids changing theory documents, scripts or package files). **All ten were fixed in the sources at `ec3644f`** (commit "Source consistency fixes found while writing the technical report"), to which the report is now pinned; the appendices are regenerated from `ec3644f`. None affected a proved result; all were labels, counts, stale cross-references or wording.

| # | Where | Issue | Severity | Status |
|---|---|---|---|---|
| 1 | `PHASE2a_rxc.md` header vs §9 | Header says `verify_phase2a.py` is 32/32 PASS; §9 says "31 PASS, 0 FAIL". The committed output has 32 PASS lines. | count (stale) | fixed at `ec3644f` |
| 2 | `PHASE2b_districts.md` header vs §7 | Header says 22/22 PASS; §7 says "21 PASS, 0 FAIL". The committed output has 22 (D8b was added in round 2). | count (stale) | fixed at `ec3644f` |
| 3 | `notation.tex` l.21; `PHASE1_foundations.md` §3 dictionary (l.279) | Both still give $\mathcal G=\sigma(X,N,Z)$. After review R-15, $Z$ was removed from the conditioning σ-algebra (P1 §1.5 and A3 now use $\sigma(X,N)$, with $D$ added in 2b). The report's notation table uses $\sigma(V)$, $V=(X,N,D)$. | notation (stale) | fixed at `ec3644f` |
| 4 | `PACKAGE_DESIGN.md` §0 item 5 | "Proof register (§8): fourteen claims"; §8 now has R1–R18 (R15–R16 Batch B, R17–R18 round 2). `EXTERNAL_REVIEW_PLAN.md` §1 and batch D also say "R1–R14". | count (stale) | fixed at `ec3644f` |
| 5 | `PACKAGE_DESIGN.md` §4.1 Part A, row M-SIGN, column "Sharp" | Still reads "sharp (M2)". Ledger row P-B1 (`B_PI_RESPONSE_PLAN.md`) listed this cell for relabelling to "sharp for continuous unit tables"; the T1.6 and M-PER cells were relabelled, this one was not. | wording (missed by P-B1) | fixed at `ec3644f` |
| 6 | `PLAIN_LANGUAGE_GUIDE.md` §8 (l.195) | "With two groups, only the 'within' sign can be contradicted by aggregate data. The 'between' sign and contextual reinforcement never produce an empty bound." This is the unqualified claim that P-B5/P-B16 had removed elsewhere: it holds only in the model-free continuous problem with no other restrictions (2c M3 Reading; 5a Remark R3.2). The same guide's §9 has the corrected wording. | wording (missed by P-B5) | fixed at `ec3644f` |
| 7 | `PHASE4_implementation.md` l.9 and l.78 | Theory column cites "P3 Thm 5" (now Theorem P3.5) and the `senc` bullet "P3 Thm 4" (now Theorem P3.4); MSIGN_REVIEW_LOG item M-6 renumbered these elsewhere. The top-level "Projection intervals" bullet still gives $q=\operatorname{rank}\hat\Sigma$; its own 0.3.1 sub-bullet supersedes this with $q=\dim\beta$. | cross-reference (stale) | fixed at `ec3644f` |
| 8 | `EXTERNAL_REVIEW_PLAN.md` §1 table | Lists `PHASE2c_monotone.md` as "none" reviewed and `PHASE5a` §§1–3 as "fresh-session adversarial review only". Both were reviewed externally in Batch B, two rounds (`reviews/B_*`). | status table (stale) | fixed at `ec3644f` |
| 9 | `PHASE3_inference.md` §8.3/§8.4 item 7 vs §5 calibration paragraph (and `PHASE4`) | Two sizes for the χ² pooling test in design S0: 0.090 (SE 0.009; F 0.079), from the 1000-replication Monte Carlo, and 0.078 (SE 0.006; F 0.072), from the 2000-replication `benchmarks/pooling_test_size.R`. Both are correct for their runs, but the text does not say they come from different runs; they differ by about 1.1 combined standard errors. | presentation | fixed at `ec3644f` |
| 10 | `benchmarks/results/ei_datasets_nz_summary.md` header | Says it was generated with "eiPartialID 0.2.0" on 2026-11-18. `REVIEW_RESPONSE_PLAN.md` WS6 and P3 §8.4 item 10 describe it as rerun with the current package (projection columns). If the run predates 0.3.0, its `proj` column uses the earlier Bhatia–Davis source (from the CLR interval) and $q=\operatorname{rank}\hat\Sigma$. Coverage figures are unlikely to change, but the version label and the text should agree. | provenance | fixed at `ec3644f` |

## After `ec3644f`

* The completeness review noted a stale code comment in `package/eiPartialID/R/rxc_infer.R` ("q = rank of V" above code that uses `ncol(V)`). It was fixed upstream at `903c215` (comment only); the report's pin, `ec3644f`, is unaffected, since `903c215` changes no document, script output or behaviour.
* No new inconsistency in the source documents was found while making the second-round report changes.

These items were noticed while compiling the report; this was not an independent re-derivation of every proof step.
