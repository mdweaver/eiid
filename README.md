### Code Version 0.3.1

This repo includes the R package for replicating the results in:

Wenxin Jiang, Gary King, Allen Schmaltz, and Martin A. Tanner. Forthcoming. “Ecological Regression with Partial Identification.” Political Analysis.

A step-by-step guide for reproducing the results in the paper is available in the associated Dataverse archive at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/W7KZVL

The Dataverse archive also includes all of the data analyzed in the paper.

For convenience, this repo only includes the standalone R package (eiPartialID) in the folder 'package' (which is similar to that in the Dataverse archive, with some minor differences to conform to CRAN policies). The package should be installed locally from source. An overview of usage is provided in the following HTML vignette, which can be viewed in your browser:

[package/eiPartialID/vignettes/rendered_html_archive/2019_02_04/eiPartialID_usage_2019_02_04.html](package/eiPartialID/vignettes/rendered_html_archive/2019_02_04/eiPartialID_usage_2019_02_04.html)

### New in 0.3.1: inference status and corrections after the second review round

`print()` now states which result covers each fit's confidence intervals: pooled M0 (proved), grouped fits without signs (proved only with few, large groups) and grouped fits with signs or reinforcement (no theorem yet). The projection critical value uses the full parameter dimension, and the rank guard is a safeguard. `proj_actual` takes its variance bound from a projection interval for the rung alone, with a note when that interval lies near 0 or 1. The guide and vignette correct the W1 explanation and the coverage figures. M-SIGN stays experimental. A technical report collecting the methods and proofs is in `docs/report/` (`eiPartialID_report.pdf`). See `package/eiPartialID/NEWS.md`.

### New in 0.3.0: signs on associations and contextual reinforcement (experimental)

`eiBounds(restrictions = list(assoc = ..., cr = ...))` adds signs on mean associations and contextual reinforcement (the M-SIGN module, `docs/extension/PHASE5a_sign_module.md`). The module is **experimental**: its sharpness and inference results are proved for the pooled model only (`docs/extension/reviews/` holds the two referees' reports and our responses); `print()` flags it, and under grouped rungs each sign is one pooled restriction across groups. This release also adds a rank guard for the projection intervals, a projection-based variance bound in `proj_actual`, a wild-bootstrap pooling test, every piece of a union in the printout, and the rung-only result beside restricted results. See `package/eiPartialID/NEWS.md`.

### New in 0.2.0: R x C tables, grouped units and an assumption ladder

`eiBounds()` extends the method to R x C tables (e.g. voter transitions), to units nested in groups (e.g. polls within ridings, with group baselines or bounded group heterogeneity), and to a ladder of assumptions from Duncan-Davis bounds to Goodman constancy. It reports intersection-bound (CLR-type) confidence intervals, with projection intervals alongside that carry an asymptotic guarantee, a pooling test for grouped data, optional additional restrictions, and sensitivity and breakdown analysis. The original 2x2 functions are unchanged, and `eiBounds()` reproduces them exactly in the 2x2 case. See the vignette `eiPartialID_rxc.Rmd`, and in `docs/extension/` the derivations (`PHASE1_foundations.md`, `PHASE2a_rxc.md`, `PHASE2b_districts.md`, `PHASE3_inference.md`), the implementation notes (`PHASE4_implementation.md`) and a plain-language guide (`PLAIN_LANGUAGE_GUIDE.md`).

