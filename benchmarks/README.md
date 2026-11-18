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

`eiDatasetsTables()` builds the tables. `benchmarks/run_ei_datasets_benchmark.R` has two parts:

1. every electorate table through the rungs DD, M0, M0+order, M0 trim and L6;
2. one table per election, pooling the polling places of the electorates where National stood, fitted by M0 and by M1 with electorates as groups. The truth is the sum of the electorates' cross-tabulations. This part also reports the pooling test (M0 against M1), the M1 test (M1 against M2) and the emptiness test.

Both parts report the CLR interval (the default, algorithmic) and the projection interval built from the same candidates (`proj_*`, P3 Prop. P3.5(b2)(ii)). The script writes:

* `results/ei_datasets_nz.csv.gz`: one row per table × rung × cell (`part` is `electorate` or `election`), with the truth, the plug-in estimate, the CLR and projection CIs, the Duncan–Davis bounds, coverage indicators and, for part 2, the test results;
* `results/ei_datasets_nz_summary.md`: coverage and width summaries, overall, for the diagonal cells and by election (part 1), and the tests and every cell (part 2).

Set `BENCH_SUMMARY_ONLY=1` to rebuild the summary from the CSV without refitting.

Headline (Nov 2026, eiPartialID 0.2.0 after WS5): per electorate, M0's CLR and projection intervals cover all 4302 cells, at about Duncan–Davis width. Pooled by election with electorates as groups, the data reject both M0 and M1 (empty estimates in every election; pooling and M1 tests with $p\approx0$). There the CLR interval is rejected in 21% of cells and misses the truth in 8%, while the projection interval never misses but is as wide as Duncan–Davis.

```sh
# from the repository root, with eiPartialID and ei.Datasets installed
Rscript benchmarks/run_ei_datasets_benchmark.R 4 100   # 4 cores, 100 draws for candidate endpoints (about 40 minutes on 2 cores)
```

The package test `tests/testthat/test-benchmark-ei-datasets.R` checks the table builder (exact truth margins) and runs the machinery on a few electorates. It is skipped when `ei.Datasets` is not installed and on CRAN.

### Adding benchmarks

Any dataset with known cell counts can be added as an entry `list(data, rows, cols, truth, group = NULL, meta = list(...))` and passed to `eiBenchmark()`. `truth` is the R×C matrix of true aggregate counts, with rows and columns in the order of `rows` and `cols`.

Not included: `ei_SCO_2007`, also in `ei.Datasets`, whose candidate columns carry names but not party affiliations.

## Size of the pooling test

`pooling_test_size.R` simulates design S0 of `docs/extension/PHASE3_inference.md` §8.3 (2×2, 20 districts × 50 polls, no shifts; 38 restrictions) and reports the rejection rate of the pooling test at nominal 5% for the HC3 Wald statistic against $\chi^2$, the F form and the restricted wild bootstrap (the default from 0.3.0). With 2,000 replications and $B=399$ (`results/pooling_test_size_output.txt`): $\chi^2$ 0.078, F 0.072, wild bootstrap 0.043 (Monte Carlo SE 0.005–0.006).

```sh
Rscript benchmarks/pooling_test_size.R 2000 399 4   # replications, bootstrap draws, cores (about 5 minutes on 4 cores)
```

