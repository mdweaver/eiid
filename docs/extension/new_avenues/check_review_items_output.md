# Checks for the review items

## (i) multiple validation tables, shared K
| R | C | dim K | rank with J = 1, 2, 3, 4 tables | tables needed | overidentifying restrictions at that J |
|---|---|---|---|---|---|
| 3 | 2 | 3 | 2, 3, 3, 3 | 2 | 1 |
| 3 | 3 | 6 | 4, 6, 6, 6 | 2 | 2 |
| 4 | 2 | 6 | 3, 6, 6, 6 | 2 | 0 |
| 4 | 3 | 12 | 6, 12, 12, 12 | 2 | 0 |
| 5 | 3 | 20 | 8, 16, 20, 20 | 3 | 4 |

## (ii) crossing partitions of the same individuals (2 groups x yes/no; 60 atoms; 12 units each way)
| design | DD, partition 1 | DD, partition 2 | intersection | joint atom LP | truth |
|---|---|---|---|---|---|
| homogeneous atoms within units | [0.050, 0.964] | [0.098, 0.939] | [0.098, 0.939] | [0.118, 0.930] | 0.722 |
| heterogeneous atoms | [0.044, 0.961] | [0.053, 0.967] | [0.053, 0.961] | [0.072, 0.954] | 0.677 |
Atoms (e.g. census blocks) have KNOWN group counts; outcomes are observed only as unit totals, in both partitions. DD per partition uses each unit's group and outcome totals.

## (iii) age-stratified margins inside units (2 groups x yes/no x 4 age bands), 200 units
DD width for group A's aggregate yes rate: one-way margins 0.748; with group x age and outcome x age margins 0.511 (ratio 0.68).

## (iv) Wisconsin: covariance of unit-level yes rates from Var(T|x) = x' Sigma x (units >= 100 men)
| year | sd veterans | sd non-veterans | correlation | PSD | cluster-bootstrap 90% interval for the correlation |
|---|---|---|---|---|---|
| 1857 | 0.084 | 0.129 | 3.35 | False | [1.01, 11.91] |
| 1865 | 0.218 | 0.106 | -0.16 | True | [-0.49, 0.70] |
