# Avenue A verification: second moments and unit size (2x2)

Truth: loyalty m1(x)=0.75+0.15x, crossover m0(x)=0.05+0.20x; kappa0=0.250.
In the estimation rows Q is treated as known so that only the variance step is noisy.

### Design low heterogeneity, NZ-like sizes: phi=60 (Beta dispersion), rho=0.5, unit size median 300, sdlog 0.9
truth B* = 0.8399;  mean Var(pi1-pi0|x) over grid = 0.0023
Identification (population moments):
| constraint set | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.521, 0.941] | 0.421 (1.00 of O) | True |
| O+V0 | [0.521, 0.941] | 0.421 (1.00 of O) | True |
| O+V1 | [0.521, 0.841] | 0.320 (0.76 of O) | True |
| O+ord | [0.578, 0.941] | 0.363 (0.86 of O) | True |
| O+V1+ord | [0.578, 0.841] | 0.263 (0.62 of O) | True |
| O+V1+R+ord | [0.839, 0.841] | 0.002 (0.00 of O) | True |
| O+V1+ord, integrated (one constraint) | [0.578, 0.840] | 0.263 | True |
| O+V1+R+ord, integrated | [0.835, 0.840] | 0.006 | True |

Estimation, reps=200.  'pointwise plug-in': grid constraints with estimated A,B (min over grid, no correction).
'integrated, 95%': one averaged constraint; its right-hand side moved by 1.645 HC0 standard errors in the conservative direction.
| p | pointwise plug-in O+V1+ord: width / cover / empty | integrated 95% O+V1+ord: width / cover / empty | integrated 95% O+V1+R+ord: width / cover / empty |
|---|---|---|---|
| 64 | 0.271 / 0.09 / 0.83 | 0.363 / 0.93 / 0.01 | 0.239 / 0.67 / 0.08 |
| 500 | 0.228 / 0.17 / 0.56 | 0.363 / 0.96 / 0.00 | 0.184 / 0.74 / 0.00 |
| 5000 | 0.203 / 0.06 / 0.32 | 0.305 / 0.97 / 0.00 | 0.086 / 0.88 / 0.00 |

### Design high heterogeneity, NZ-like sizes: phi=8 (Beta dispersion), rho=0.5, unit size median 300, sdlog 0.9
truth B* = 0.8400;  mean Var(pi1-pi0|x) over grid = 0.0173
Identification (population moments):
| constraint set | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.521, 0.941] | 0.420 (1.00 of O) | True |
| O+V0 | [0.521, 0.941] | 0.420 (1.00 of O) | True |
| O+V1 | [0.521, 0.845] | 0.324 (0.77 of O) | True |
| O+ord | [0.578, 0.941] | 0.363 (0.86 of O) | True |
| O+V1+ord | [0.578, 0.845] | 0.266 (0.63 of O) | True |
| O+V1+R+ord | [0.832, 0.845] | 0.013 (0.03 of O) | True |
| O+V1+ord, integrated (one constraint) | [0.578, 0.845] | 0.266 | True |
| O+V1+R+ord, integrated | [0.805, 0.845] | 0.040 | True |

Estimation, reps=200.  'pointwise plug-in': grid constraints with estimated A,B (min over grid, no correction).
'integrated, 95%': one averaged constraint; its right-hand side moved by 1.645 HC0 standard errors in the conservative direction.
| p | pointwise plug-in O+V1+ord: width / cover / empty | integrated 95% O+V1+ord: width / cover / empty | integrated 95% O+V1+R+ord: width / cover / empty |
|---|---|---|---|
| 64 | 0.362 / 0.07 / 0.91 | 0.363 / 0.93 / 0.05 | 0.363 / 0.81 / 0.12 |
| 500 | 0.264 / 0.08 / 0.84 | 0.362 / 0.96 / 0.03 | 0.362 / 0.87 / 0.05 |
| 5000 | 0.165 / 0.04 / 0.67 | 0.362 / 0.97 / 0.00 | 0.310 / 0.93 / 0.00 |

### Design low heterogeneity, large units: phi=60 (Beta dispersion), rho=0.5, unit size median 3000, sdlog 0.5
truth B* = 0.8398;  mean Var(pi1-pi0|x) over grid = 0.0023
Identification (population moments):
| constraint set | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.520, 0.941] | 0.422 (1.00 of O) | True |
| O+V0 | [0.520, 0.941] | 0.422 (1.00 of O) | True |
| O+V1 | [0.520, 0.840] | 0.321 (0.76 of O) | True |
| O+ord | [0.577, 0.941] | 0.364 (0.86 of O) | True |
| O+V1+ord | [0.577, 0.840] | 0.263 (0.62 of O) | True |
| O+V1+R+ord | [0.838, 0.840] | 0.002 (0.00 of O) | True |
| O+V1+ord, integrated (one constraint) | [0.577, 0.840] | 0.263 | True |
| O+V1+R+ord, integrated | [0.834, 0.840] | 0.006 | True |

Estimation, reps=200.  'pointwise plug-in': grid constraints with estimated A,B (min over grid, no correction).
'integrated, 95%': one averaged constraint; its right-hand side moved by 1.645 HC0 standard errors in the conservative direction.
| p | pointwise plug-in O+V1+ord: width / cover / empty | integrated 95% O+V1+ord: width / cover / empty | integrated 95% O+V1+R+ord: width / cover / empty |
|---|---|---|---|
| 64 | 0.365 / 0.07 / 0.92 | 0.364 / 0.94 / 0.05 | 0.364 / 0.79 / 0.17 |
| 500 | 0.364 / 0.06 / 0.91 | 0.364 / 0.97 / 0.02 | 0.364 / 0.89 / 0.09 |
| 5000 | 0.234 / 0.04 / 0.88 | 0.364 / 0.95 / 0.01 | 0.364 / 0.89 / 0.04 |

### Design low heterogeneity, independent shocks: phi=60 (Beta dispersion), rho=0.0, unit size median 300, sdlog 0.9
truth B* = 0.8401;  mean Var(pi1-pi0|x) over grid = 0.0044
Identification (population moments):
| constraint set | bound on loyalty B* | width | contains truth |
|---|---|---|---|
| O | [0.522, 0.941] | 0.419 (1.00 of O) | True |
| O+V0 | [0.522, 0.941] | 0.419 (1.00 of O) | True |
| O+V1 | [0.522, 0.841] | 0.319 (0.76 of O) | True |
| O+ord | [0.579, 0.941] | 0.362 (0.86 of O) | True |
| O+V1+ord | [0.579, 0.841] | 0.262 (0.63 of O) | True |
| O+V1+R+ord | [0.840, 0.841] | 0.001 (0.00 of O) | True |
| O+V1+ord, integrated (one constraint) | [0.579, 0.841] | 0.262 | True |
| O+V1+R+ord, integrated | [0.837, 0.841] | 0.004 | True |

Estimation, reps=200.  'pointwise plug-in': grid constraints with estimated A,B (min over grid, no correction).
'integrated, 95%': one averaged constraint; its right-hand side moved by 1.645 HC0 standard errors in the conservative direction.
| p | pointwise plug-in O+V1+ord: width / cover / empty | integrated 95% O+V1+ord: width / cover / empty | integrated 95% O+V1+R+ord: width / cover / empty |
|---|---|---|---|
| 64 | 0.238 / 0.13 / 0.70 | 0.362 / 0.94 / 0.01 | 0.198 / 0.61 / 0.09 |
| 500 | 0.197 / 0.10 / 0.56 | 0.347 / 0.95 / 0.00 | 0.182 / 0.80 / 0.01 |
| 5000 | 0.210 / 0.06 / 0.36 | 0.300 / 0.98 / 0.00 | 0.078 / 0.91 / 0.00 |
