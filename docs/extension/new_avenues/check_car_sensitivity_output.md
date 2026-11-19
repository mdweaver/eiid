# CAR sensitivity checks (Section 14)

## Part 1: invisible CAR violation

max |E[Y|X] - eta'X| = 1.11e-16  (exactly linear: invisible)
true beta_1 = 0.6800;  CAR/Goodman limit = 0.6000;  bias = +0.0800
formula k*(1 - E[x1^2]/E[x1]) = +0.0800

## Part 2: omitted-confounder bound when the confounder tracks composition

DGP: x1 ~ Beta(2,2); A = x1 + tau*nu, nu ~ N(0,1); B = eta + k*(1-A, -A); Y = B'X.
CAR holds given A; the short (no-A) CEF is exactly eta'X, so the short model is correct
in form and the CAR estimate converges to eta_1. Target beta_1 (equal unit sizes).

| tau | R2(x1~A) | C_alpha | rho | S | true bias | bound (rho=1, true C) | bound with C_alpha=0.65 | DD beta_1 interval | delta needed | kappa needed |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 0.992 | 8.42 | +0.97 | 0.0098 | +0.0800 | 0.0827 | 0.0064 | [0.197, 0.813] | 0.20 | 0.44 |
| 0.05 | 0.952 | 3.36 | +0.97 | 0.0245 | +0.0800 | 0.0822 | 0.0159 | [0.197, 0.813] | 0.20 | 0.43 |
| 0.1 | 0.833 | 1.67 | +0.98 | 0.0490 | +0.0801 | 0.0816 | 0.0318 | [0.197, 0.813] | 0.20 | 0.41 |
| 0.2 | 0.556 | 0.84 | +0.98 | 0.0979 | +0.0800 | 0.0818 | 0.0637 | [0.198, 0.810] | 0.20 | 0.33 |
| 0.3 | 0.357 | 0.56 | +0.98 | 0.1470 | +0.0799 | 0.0816 | 0.0955 | [0.201, 0.806] | 0.20 | 0.25 |

Reading: the true bias (0.4k = 0.08) does not depend on tau. As tau -> 0 the confounder
becomes composition itself (a contextual effect): S -> 0 but C_alpha -> infinity. A
C_alpha calibrated on an ordinary covariate (0.65 here, i.e. R2 about 0.3) gives a bound
far below the true bias when tau is small. The Duncan-Davis interval is always valid.

## Part 3: CAR-delta (sup-norm) width and overlap, 2x2, no covariates

Composition confined to [a, b]; CAR holds in truth (g(x) = eta'x); |Delta| <= delta.
Half-widths from an LP on a 41-point grid with equal weights.

| a | b | delta | eta_1 half-width (LP) | closed form delta(2-a-b)/(b-a) | beta_1 half-width (LP) |
|---|---|---|---|---|---|
| 0.0 | 1.0 | 0.02 | 0.0200 | 0.0200 | 0.0195 |
| 0.0 | 1.0 | 0.05 | 0.0500 | 0.0500 | 0.0488 |
| 0.1 | 0.9 | 0.02 | 0.0250 | 0.0250 | 0.0266 |
| 0.1 | 0.9 | 0.05 | 0.0625 | 0.0625 | 0.0664 |
| 0.2 | 0.8 | 0.02 | 0.0333 | 0.0333 | 0.0369 |
| 0.2 | 0.8 | 0.05 | 0.0833 | 0.0833 | 0.0924 |
| 0.3 | 0.7 | 0.02 | 0.0500 | 0.0500 | 0.0557 |
| 0.3 | 0.7 | 0.05 | 0.1250 | 0.1250 | 0.1391 |
| 0.4 | 0.6 | 0.02 | 0.1000 | 0.1000 | 0.1077 |
| 0.4 | 0.6 | 0.05 | 0.2500 | 0.2500 | 0.2693 |

Reading: under CAR-delta, poor overlap multiplies delta by (2-a-b)/(b-a) for the
vertex rate. This is the partial-identification counterpart of a large Riesz-representer
norm (McCartan & Kuriwaki's positivity diagnostic). The aggregate target beta_1 is
affected about as much: it too needs the two groups' rates separated at each composition.
