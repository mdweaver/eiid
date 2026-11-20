#!/bin/sh
# Build the paper and its supplement (cross-references between them need several passes).
set -e
cd "$(dirname "$0")"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error eiPartialID_paper.tex >/dev/null
  pdflatex -interaction=nonstopmode -halt-on-error eiPartialID_supplement.tex >/dev/null
done
bibtex eiPartialID_paper >/dev/null
bibtex eiPartialID_supplement >/dev/null
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error eiPartialID_paper.tex >/dev/null
  pdflatex -interaction=nonstopmode -halt-on-error eiPartialID_supplement.tex >/dev/null
done
echo "build ok"
