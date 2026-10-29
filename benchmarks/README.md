# Benchmarks

Benchmarks of `eiPartialID` (0.2.0 and later) on real data where the true cross-tabulation is known.

## ei.Datasets: New Zealand general elections 2002–2020

Source: the CRAN package [`ei.Datasets`](https://CRAN.R-project.org/package=ei.Datasets) (Pavia 2023, *Social Science Computer Review*).

In New Zealand's mixed-member system each voter casts a party vote and an electorate (candidate) vote. For every electorate, polling-place counts of both votes are published, together with the official electorate-level cross-tabulation of party vote by candidate vote. That cross-tabulation aggregates exactly the same voters as the polling places, so it is an exact truth for an ecological-inference problem:

* **units:** polling places in the electorate (median 64 per electorate);
* **rows:** party vote (Labour, National, all other parties and informal votes);
* **columns:** electorate vote (Labour candidate, National candidate, all others and informal);
* **truth:** the official cross-tabulation, collapsed to the same categories;
* **targets:** all nine row-by-column rates, i.e. the proportions of each party's voters who voted for each candidate (split-ticket voting).

There are 492 tables (7 elections × 69–72 electorates). In 42 electorates National fielded no candidate, so those tables are 3×2.

`eiDatasetsTables()` builds the tables. `benchmarks/run_ei_datasets_benchmark.R` runs every table through the rungs and writes:

* `results/ei_datasets_nz.csv.gz`: one row per table × rung × cell, with the truth, the plug-in estimate, the CIs, the Duncan–Davis bounds and coverage indicators;
* `results/ei_datasets_nz_summary.md`: coverage and width summaries, overall, for the diagonal cells, and by election.

```sh
# from the repository root, with eiPartialID and ei.Datasets installed
Rscript benchmarks/run_ei_datasets_benchmark.R 4 100   # 4 cores, 100 draws for candidate endpoints
```

The package test `tests/testthat/test-benchmark-ei-datasets.R` checks the table builder (exact truth margins) and runs the machinery on a few electorates. It is skipped when `ei.Datasets` is not installed and on CRAN.

### Adding benchmarks

Any dataset with known cell counts can be added as an entry `list(data, rows, cols, truth, group = NULL, meta = list(...))` and passed to `eiBenchmark()`. `truth` is the R×C matrix of true aggregate counts, with rows and columns in the order of `rows` and `cols`.

Not included: `ei_SCO_2007`, also in `ei.Datasets`, whose candidate columns carry names but not party affiliations.
