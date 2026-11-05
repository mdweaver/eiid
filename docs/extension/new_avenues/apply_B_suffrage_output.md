# Avenue B applied to the veterans / suffrage data (quick look, plug-in, no CIs)

### Wisconsin 1857: 364 units, 52 counties; veteran share mean 0.353
separation of own and county share: SD(x) 0.196, SD(x - z) 0.192, corr(x, z) 0.27; units alone in their county: 0.04
identified county-share effects on the YES rate: veterans g_v = -0.691, non-veterans g_n = +0.128; gap in own-share slopes d = +0.025
**all units as constraint points**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.015, 0.463] | [0.046, 0.290] | [-0.275, 0.417] |
| M0 (pooled LCE) | [0.122, 0.365] | [0.097, 0.230] | [-0.108, 0.268] |
| B0 (+ county share) | [0.209, 0.351] | [0.099, 0.176] | [0.033, 0.252] |
| B1(M=2), yes column only | [0.209, 0.351] | [0.099, 0.176] | [0.033, 0.252] |
| B1(M=1), yes column only | [0.209, 0.347] | [0.101, 0.176] | [0.033, 0.246] |
| B1(M=0.5), yes column only | [0.243, 0.312] | [0.120, 0.158] | [0.085, 0.192] |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): 0.05
M_min (all columns, 0.05 grid up to 8): 3.85

**trimmed constraint points (N >= 100, x and z inside 2.5-97.5%; 295 units)**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.015, 0.463] | [0.046, 0.290] | [-0.275, 0.417] |
| M0 (pooled LCE) | [0.119, 0.397] | [0.080, 0.231] | [-0.112, 0.317] |
| B0 (+ county share) | [0.168, 0.385] | [0.080, 0.198] | [-0.030, 0.305] |
| B1(M=2), yes column only | [0.168, 0.385] | [0.080, 0.198] | [-0.030, 0.305] |
| B1(M=1), yes column only | [0.208, 0.347] | [0.101, 0.177] | [0.031, 0.246] |
| B1(M=0.5), yes column only | [0.243, 0.312] | [0.120, 0.158] | [0.085, 0.192] |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): 0.05
M_min (all columns, 0.05 grid up to 8): 3.85


### Wisconsin 1865: 364 units, 52 counties; veteran share mean 0.268
separation of own and county share: SD(x) 0.150, SD(x - z) 0.147, corr(x, z) 0.27; units alone in their county: 0.04
identified county-share effects on the YES rate: veterans g_v = -0.704, non-veterans g_n = +0.139; gap in own-share slopes d = +0.245
**all units as constraint points**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.003, 0.637] | [0.044, 0.276] | [-0.273, 0.594] |
| M0 (pooled LCE) | [0.161, 0.406] | [0.131, 0.221] | [-0.061, 0.275] |
| B0 (+ county share) | empty | empty | empty |
| B1(M=2), yes column only | empty | empty | empty |
| B1(M=1), yes column only | empty | empty | empty |
| B1(M=0.5), yes column only | empty | empty | empty |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): None
M_min (all columns, 0.05 grid up to 8): None

**trimmed constraint points (N >= 100, x and z inside 2.5-97.5%; 309 units)**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.003, 0.637] | [0.044, 0.276] | [-0.273, 0.594] |
| M0 (pooled LCE) | [0.156, 0.536] | [0.084, 0.223] | [-0.066, 0.452] |
| B0 (+ county share) | [0.186, 0.525] | [0.085, 0.209] | [-0.023, 0.440] |
| B1(M=2), yes column only | [0.186, 0.510] | [0.090, 0.209] | [-0.023, 0.420] |
| B1(M=1), yes column only | [0.238, 0.420] | [0.123, 0.190] | [0.049, 0.297] |
| B1(M=0.5), yes column only | [0.284, 0.374] | [0.140, 0.173] | [0.111, 0.235] |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): 0.3
M_min (all columns, 0.05 grid up to 8): 6.65


### Iowa 1857: 205 units, 80 counties; veteran share mean 0.488
separation of own and county share: SD(x) 0.255, SD(x - z) 0.228, corr(x, z) 0.46; units alone in their county: 0.31
identified county-share effects on the YES rate: veterans g_v = +0.076, non-veterans g_n = +0.025; gap in own-share slopes d = +0.047
**all units as constraint points**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.006, 0.147] | [0.005, 0.138] | [-0.132, 0.142] |
| M0 (pooled LCE) | [0.069, 0.129] | [0.043, 0.100] | [-0.031, 0.086] |
| B0 (+ county share) | [0.074, 0.128] | [0.047, 0.098] | [-0.024, 0.081] |
| B1(M=2), yes column only | [0.075, 0.118] | [0.056, 0.097] | [-0.022, 0.062] |
| B1(M=1), yes column only | [0.086, 0.107] | [0.066, 0.087] | [-0.001, 0.041] |
| B1(M=0.5), yes column only | [0.100, 0.102] | [0.071, 0.073] | [0.027, 0.031] |
| B1(M=2), all columns | [0.075, 0.118] | [0.056, 0.097] | [-0.022, 0.062] |
| B1(M=1), all columns | [0.086, 0.107] | [0.066, 0.087] | [-0.001, 0.041] |
M_min (yes column only, 0.05 grid up to 8): 0.5
M_min (all columns, 0.05 grid up to 8): 1.0

**trimmed constraint points (N >= 100, x and z inside 2.5-97.5%; 155 units)**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.006, 0.147] | [0.005, 0.138] | [-0.132, 0.142] |
| M0 (pooled LCE) | [0.067, 0.129] | [0.043, 0.102] | [-0.035, 0.086] |
| B0 (+ county share) | [0.072, 0.128] | [0.047, 0.100] | [-0.028, 0.081] |
| B1(M=2), yes column only | [0.075, 0.118] | [0.056, 0.097] | [-0.022, 0.062] |
| B1(M=1), yes column only | [0.086, 0.107] | [0.066, 0.087] | [-0.001, 0.041] |
| B1(M=0.5), yes column only | [0.100, 0.102] | [0.071, 0.073] | [0.027, 0.031] |
| B1(M=2), all columns | [0.075, 0.118] | [0.056, 0.097] | [-0.022, 0.062] |
| B1(M=1), all columns | [0.086, 0.107] | [0.066, 0.087] | [-0.001, 0.041] |
M_min (yes column only, 0.05 grid up to 8): 0.5
M_min (all columns, 0.05 grid up to 8): 1.0


### Iowa 1868: 223 units, 97 counties; veteran share mean 0.272
separation of own and county share: SD(x) 0.202, SD(x - z) 0.189, corr(x, z) 0.36; units alone in their county: 0.36
identified county-share effects on the YES rate: veterans g_v = -0.210, non-veterans g_n = -0.271; gap in own-share slopes d = -0.632
**all units as constraint points**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.016, 0.929] | [0.341, 0.681] | [-0.664, 0.588] |
| M0 (pooled LCE) | empty | empty | empty |
| B0 (+ county share) | empty | empty | empty |
| B1(M=2), yes column only | empty | empty | empty |
| B1(M=1), yes column only | empty | empty | empty |
| B1(M=0.5), yes column only | empty | empty | empty |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): None
M_min (all columns, 0.05 grid up to 8): None

**trimmed constraint points (N >= 100, x and z inside 2.5-97.5%; 196 units)**
| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |
|---|---|---|---|
| Duncan-Davis | [0.016, 0.929] | [0.341, 0.681] | [-0.664, 0.588] |
| M0 (pooled LCE) | empty | empty | empty |
| B0 (+ county share) | empty | empty | empty |
| B1(M=2), yes column only | empty | empty | empty |
| B1(M=1), yes column only | empty | empty | empty |
| B1(M=0.5), yes column only | empty | empty | empty |
| B1(M=2), all columns | empty | empty | empty |
| B1(M=1), all columns | empty | empty | empty |
M_min (yes column only, 0.05 grid up to 8): None
M_min (all columns, 0.05 grid up to 8): None


### DiD of the veteran - non-veteran yes gap, no cross-period link (interval difference; trimmed points; B1 on the yes column)
| state | M0 | B0 | B1(M=2) | B1(M=1) | B1(M=0.5) |
|---|---|---|---|---|---|
| Wisconsin | [-0.383, 0.564] | [-0.328, 0.470] | [-0.328, 0.451] | [-0.198, 0.266] | [-0.082, 0.150] |
| Iowa | empty | empty | empty | empty | empty |
