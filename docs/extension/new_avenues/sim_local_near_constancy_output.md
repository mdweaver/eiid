# Local near-constancy with a data-implied floor (2x2)

## low heterogeneity (sd 0.03), large units: 200 units, 5 nearest-composition neighbours; truth 0.397; DD [0.203, 0.719]
| version | eps_min | eps | bound on group W's aggregate rate | width / DD | covers |
|---|---|---|---|---|---|
| sup | 0.113 | 1.0x = 0.113 | [0.225, 0.718] | 0.96 | True |
| sup | 0.113 | 1.5x = 0.170 | [0.203, 0.719] | 1.00 | True |
| sup | 0.113 | 2.0x = 0.227 | [0.203, 0.719] | 1.00 | True |
| sup | 0.113 | 3.0x = 0.340 | [0.203, 0.719] | 1.00 | True |
| mean | 0.019 | 1.0x = 0.019 | [0.448, 0.448] | 0.00 | False |
| mean | 0.019 | 1.5x = 0.028 | [0.203, 0.719] | 1.00 | True |
| mean | 0.019 | 2.0x = 0.038 | [0.203, 0.719] | 1.00 | True |
| mean | 0.019 | 3.0x = 0.057 | [0.203, 0.719] | 1.00 | True |

## high heterogeneity (sd 0.08), large units: 200 units, 5 nearest-composition neighbours; truth 0.391; DD [0.198, 0.733]
| version | eps_min | eps | bound on group W's aggregate rate | width / DD | covers |
|---|---|---|---|---|---|
| sup | 0.317 | 1.0x = 0.317 | [0.214, 0.732] | 0.97 | True |
| sup | 0.317 | 1.5x = 0.476 | [0.198, 0.733] | 1.00 | True |
| sup | 0.317 | 2.0x = 0.634 | [0.198, 0.733] | 1.00 | True |
| sup | 0.317 | 3.0x = 0.951 | [0.198, 0.733] | 1.00 | True |
| mean | 0.051 | 1.0x = 0.051 | [0.569, 0.569] | -0.00 | False |
| mean | 0.051 | 1.5x = 0.077 | [0.198, 0.733] | 1.00 | True |
| mean | 0.051 | 2.0x = 0.103 | [0.198, 0.733] | 1.00 | True |
| mean | 0.051 | 3.0x = 0.154 | [0.198, 0.733] | 1.00 | True |

## low heterogeneity (sd 0.03), units of ~200 (binomial noise): 200 units, 5 nearest-composition neighbours; truth 0.395; DD [0.196, 0.734]
| version | eps_min | eps | bound on group W's aggregate rate | width / DD | covers |
|---|---|---|---|---|---|
| sup | 0.676 | 1.0x = 0.676 | [0.198, 0.733] | 1.00 | True |
| sup | 0.676 | 1.5x = 1.015 | [0.196, 0.734] | 1.00 | True |
| sup | 0.676 | 2.0x = 1.353 | [0.196, 0.734] | 1.00 | True |
| sup | 0.676 | 3.0x = 2.029 | [0.196, 0.734] | 1.00 | True |
| mean | 0.036 | 1.0x = 0.036 | [0.484, 0.484] | -0.00 | False |
| mean | 0.036 | 1.5x = 0.054 | [0.196, 0.734] | 1.00 | True |
| mean | 0.036 | 2.0x = 0.072 | [0.196, 0.734] | 1.00 | True |
| mean | 0.036 | 3.0x = 0.108 | [0.196, 0.734] | 1.00 | True |


## All-pairs version: |b_i - b_j| <= e + L|x_i - x_j| for every pair (100 units)
Truth: group slopes +0.30 and -0.20 in x, unit heterogeneity sd s.  e_min(L) = smallest feasible e.
| design | L | e_min(L) | e used | bound on W's aggregate rate | width / DD | covers |
|---|---|---|---|---|---|---|
| sd 0.03, large units: e the truth needs at L=0.0 = 0.344 (rows with e below this assume something false) | | | | | | |
| sd 0.03, large units (truth 0.394) | 0.0 | 0.175 | 1.0x = 0.175 | [0.381, 0.454] | 0.13 | True |
| sd 0.03, large units (truth 0.394) | 0.0 | 0.175 | 1.5x = 0.262 | [0.312, 0.522] | 0.38 | True |
| sd 0.03, large units (truth 0.394) | 0.0 | 0.175 | 2.0x = 0.349 | [0.255, 0.581] | 0.59 | True |
| sd 0.03, large units: e the truth needs at L=0.3 = 0.119 (rows with e below this assume something false) | | | | | | |
| sd 0.03, large units | 0.3 | 0.106 | 1.0x | infeasible | | |
| sd 0.03, large units (truth 0.394) | 0.3 | 0.106 | 1.5x = 0.160 | [0.285, 0.552] | 0.49 | True |
| sd 0.03, large units (truth 0.394) | 0.3 | 0.106 | 2.0x = 0.213 | [0.245, 0.593] | 0.63 | True |
| sd 0.03, large units: e the truth needs at L=0.6 = 0.109 (rows with e below this assume something false) | | | | | | |
| sd 0.03, large units (truth 0.394) | 0.6 | 0.104 | 1.0x = 0.104 | [0.252, 0.386] | 0.24 | False |
| sd 0.03, large units (truth 0.394) | 0.6 | 0.104 | 1.5x = 0.156 | [0.213, 0.656] | 0.81 | True |
| sd 0.03, large units (truth 0.394) | 0.6 | 0.104 | 2.0x = 0.208 | [0.199, 0.685] | 0.89 | True |
| sd 0.08, large units: e the truth needs at L=0.0 = 0.506 (rows with e below this assume something false) | | | | | | |
| sd 0.08, large units (truth 0.393) | 0.0 | 0.363 | 1.0x = 0.363 | [0.327, 0.471] | 0.30 | True |
| sd 0.08, large units (truth 0.393) | 0.0 | 0.363 | 1.5x = 0.544 | [0.229, 0.598] | 0.76 | True |
| sd 0.08, large units (truth 0.393) | 0.0 | 0.363 | 2.0x = 0.726 | [0.223, 0.692] | 0.96 | True |
| sd 0.08, large units: e the truth needs at L=0.3 = 0.397 (rows with e below this assume something false) | | | | | | |
| sd 0.08, large units (truth 0.393) | 0.3 | 0.296 | 1.0x = 0.296 | [0.305, 0.435] | 0.27 | True |
| sd 0.08, large units (truth 0.393) | 0.3 | 0.296 | 1.5x = 0.445 | [0.224, 0.632] | 0.84 | True |
| sd 0.08, large units (truth 0.393) | 0.3 | 0.296 | 2.0x = 0.593 | [0.223, 0.699] | 0.98 | True |
| sd 0.08, large units: e the truth needs at L=0.6 = 0.358 (rows with e below this assume something false) | | | | | | |
| sd 0.08, large units (truth 0.393) | 0.6 | 0.294 | 1.0x = 0.294 | [0.236, 0.495] | 0.53 | True |
| sd 0.08, large units (truth 0.393) | 0.6 | 0.294 | 1.5x = 0.440 | [0.223, 0.691] | 0.96 | True |
| sd 0.08, large units (truth 0.393) | 0.6 | 0.294 | 2.0x = 0.587 | [0.223, 0.710] | 1.00 | True |
| sd 0.03, units ~200: e the truth needs at L=0.0 = 0.571 (rows with e below this assume something false) | | | | | | |
| sd 0.03, units ~200 (truth 0.390) | 0.0 | 0.238 | 1.0x = 0.238 | [0.357, 0.458] | 0.19 | True |
| sd 0.03, units ~200 (truth 0.390) | 0.0 | 0.238 | 1.5x = 0.357 | [0.273, 0.549] | 0.51 | True |
| sd 0.03, units ~200 (truth 0.390) | 0.0 | 0.238 | 2.0x = 0.476 | [0.206, 0.627] | 0.78 | True |
| sd 0.03, units ~200: e the truth needs at L=0.3 = 0.566 (rows with e below this assume something false) | | | | | | |
| sd 0.03, units ~200 (truth 0.390) | 0.3 | 0.182 | 1.0x = 0.182 | [0.321, 0.407] | 0.16 | True |
| sd 0.03, units ~200 (truth 0.390) | 0.3 | 0.182 | 1.5x = 0.273 | [0.229, 0.594] | 0.68 | True |
| sd 0.03, units ~200 (truth 0.390) | 0.3 | 0.182 | 2.0x = 0.364 | [0.196, 0.660] | 0.86 | True |
| sd 0.03, units ~200: e the truth needs at L=0.6 = 0.561 (rows with e below this assume something false) | | | | | | |
| sd 0.03, units ~200 | 0.6 | 0.158 | 1.0x | infeasible | | |
| sd 0.03, units ~200 (truth 0.390) | 0.6 | 0.158 | 1.5x = 0.237 | [0.198, 0.664] | 0.87 | True |
| sd 0.03, units ~200 (truth 0.390) | 0.6 | 0.158 | 2.0x = 0.316 | [0.191, 0.703] | 0.95 | True |
