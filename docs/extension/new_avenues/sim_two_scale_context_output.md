# Avenue B verification: context measured at a second scale (2x2)

a1=0.6, a0=0.03, g1=0.15, g0=0.25; own-poll effects (b1,b0) = lambda (g1,g0).

## lambda = 0.4 (truth B* = 0.727; DD at the mean shares [0.309, 0.880])
Identification (population coefficients):
| restriction | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.495, 0.805] | 0.310 | True |
| RM(M=2) | [0.637, 0.805] | 0.168 | True |
| RM(M=1) | [0.693, 0.805] | 0.112 | True |
| RM(M=0.5) | [0.722, 0.778] | 0.057 | True |
| RM(M=0) = no own effect | empty | — | False (refuted) |
| PROP | [0.727, 0.727] | 0.000 | True |

Estimation (93 neighbourhoods x 10 polls, 300 reps; plug-in):
| restriction | median width | covers truth | empty |
|---|---|---|---|
| O | 0.332 | 1.00 | 0.00 |
| RM(M=2) | 0.176 | 0.97 | 0.00 |
| RM(M=1) | 0.104 | 0.84 | 0.00 |
| RM(M=0.5) | 0.045 | 0.30 | 0.00 |
| RM(M=0) = no own effect | nan | 0.00 | 1.00 |
| PROP | 0.000 | 0.00 | 0.17 |
PROP point estimate (ignoring the O rows): median abs error 0.041, 5-95% error range [-0.542, 0.392]

## lambda = 1.5 (truth B* = 0.828; DD at the mean shares [0.442, 0.996])
Identification (population coefficients):
| restriction | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.689, 0.944] | 0.255 | True |
| RM(M=2) | [0.800, 0.944] | 0.144 | True |
| RM(M=1) | [0.876, 0.944] | 0.068 | False |
| RM(M=0.5) | [0.924, 0.943] | 0.019 | False |
| RM(M=0) = no own effect | empty | — | False (refuted) |
| PROP | [0.828, 0.828] | 0.000 | True |

Estimation (93 neighbourhoods x 10 polls, 300 reps; plug-in):
| restriction | median width | covers truth | empty |
|---|---|---|---|
| O | 0.267 | 1.00 | 0.00 |
| RM(M=2) | 0.137 | 0.67 | 0.00 |
| RM(M=1) | 0.072 | 0.00 | 0.01 |
| RM(M=0.5) | 0.023 | 0.00 | 0.21 |
| RM(M=0) = no own effect | nan | 0.00 | 1.00 |
| PROP | 0.000 | 0.00 | 0.29 |
PROP point estimate (ignoring the O rows): median abs error 0.074, 5-95% error range [-0.606, 0.704]

