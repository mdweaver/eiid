# Technical report: partial identification for R×C ecological inference

`eiPartialID_report.tex` (with `preamble.tex`, `sections/`, `generated/`, `references.bib`) compiles to
`eiPartialID_report.pdf`:

    pdflatex eiPartialID_report && bibtex eiPartialID_report && pdflatex eiPartialID_report && pdflatex eiPartialID_report

* `sections/` — the main text (written for this report) and the consolidated revision table (Appendix A).
* `generated/` — Appendices B–H: the repository's theory documents reproduced verbatim, produced by
  `python3 -I tools/make_appendices.py` (pandoc). Do not edit by hand; regenerate after changing a source.
* `tools/coverage_check.py` — writes `COVERAGE.md` (every numbered result and R1–R18 found in the report).
* `ISSUES.md` — inconsistencies noticed in the source documents (recorded, not fixed).
* `OUTLINE.md` — the stage-1 outline.

Pinned to commit `b270285` of `claude/wizardly-volta-lbh4bm`. Requires TeX Live (latex-extra) and pandoc ≥ 3.
