# NZ exploratory diagnostics

tables: 450, polling places: 31359, median size 262, 10%/90% sizes 24/1195, share with n<100: 0.29
share of units whose party-vote and candidate-vote totals differ by >2%: 0.00

## D1: is the 1/n (binomial) term visible, and is it binomial-sized?
| column | slope on b/n | cluster SE | units | reading |
|---|---|---|---|---|
| LAB | 0.693 | 0.104 | 31359 | binomial-sized |
| NAT | 0.808 | 0.115 | 31359 | binomial-sized |
| OTH | 0.690 | 0.098 | 31359 | binomial-sized |

## D2: size and persistence of the part of the truth that M0 attributes to K
| cell (row->col) | median abs v | corr with same electorate, previous election | pairs |
|---|---|---|---|
| LAB->LAB | 0.150 | 0.35 | 357 |
| LAB->NAT | 0.127 | 0.26 | 357 |
| LAB->OTH | 0.036 | 0.26 | 357 |
| NAT->LAB | 0.113 | 0.48 | 357 |
| NAT->NAT | 0.136 | 0.51 | 357 |
| NAT->OTH | 0.028 | 0.03 | 357 |
| OTH->LAB | 0.052 | 0.28 | 357 |
| OTH->NAT | 0.058 | 0.28 | 357 |
| OTH->OTH | 0.071 | 0.22 | 357 |
