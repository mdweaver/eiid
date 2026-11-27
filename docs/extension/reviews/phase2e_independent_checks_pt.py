"""Reproduce the two examples in the PHASE2e independent review.
Run with Python 3, NumPy, and SciPy. Results are percentage points.
"""

import numpy as np
from scipy.optimize import linprog


def att_bounds(size, veterans, pre, post, comparison, lce):
    """Expected veteran ATT under the specified PT and LCE assumptions.

    Each row of post contains equally probable aggregate yes-vote totals
    for one precinct. The baseline aggregate total is fixed in each precinct.
    No group-specific votes or audit restrictions are supplied.
    """
    veterans = np.asarray(veterans, dtype=float)
    controls = size - veterans
    x = veterans / size
    pre = np.asarray(pre, dtype=float)
    post = np.asarray(post, dtype=float)
    n = len(x)

    # Decision vector: veteran pre, non-veteran pre, veteran post,
    # non-veteran post. Each block contains n expected voting rates.
    v0, c0, v1, c1 = np.eye(4 * n).reshape(4, n, 4 * n)
    dv, dc = v1 - v0, c1 - c0
    wv, wc = veterans / veterans.sum(), controls / controls.sum()

    eq, rhs = [], []
    var_bounds = [(0.0, 1.0)] * (4 * n)

    for t, (totals, v, c) in enumerate(((pre[:, None], v0, c0),
                                          (post, v1, c1))):
        # Bound veteran votes separately in each observable state, THEN
        # average the endpoints. Bounding only the average total loses
        # information when the state totals straddle a capacity threshold.
        low = np.maximum(0, totals - controls[:, None]).mean(axis=1) / veterans
        high = np.minimum(veterans[:, None], totals).mean(axis=1) / veterans
        start = 2 * t * n
        var_bounds[start:start + n] = list(zip(low, high))

        eq.extend(x[:, None] * v + (1 - x[:, None]) * c)
        rhs.extend(totals.mean(axis=1) / size)

    # A group profile is affine in x iff its adjacent slopes are equal.
    # Both examples have distinct compositions ordered from low to high.
    profiles = {"pre": (v0, c0), "post": (v1, c1), "change": (dv, dc)}
    for period in lce:
        for profile in profiles[period]:
            slopes = np.diff(profile, axis=0) / np.diff(x)[:, None]
            eq.extend(np.diff(slopes, axis=0))
            rhs.extend(np.zeros(n - 2))

    # Under no anticipation and zero comparison-group spillover,
    # baseline veteran rate + comparison trend is the untreated post rate.
    if comparison == "within":
        untreated_post = v0 + dc                   # one rate per precinct
        objective = wv @ (dv - dc)
    elif comparison == "pooled":
        untreated_post = (wv @ v0 + wc @ dc)[None, :]  # one electorate rate
        objective = wv @ dv - wc @ dc
    else:
        raise ValueError(comparison)

    # Every imputed untreated post rate must lie in [0, 1].
    ub = np.vstack((untreated_post, -untreated_post))
    cap = np.r_[np.ones(len(untreated_post)), np.zeros(len(untreated_post))]

    endpoints = []
    for sign in (1, -1):
        fit = linprog(sign * objective, A_eq=eq, b_eq=rhs,
                      A_ub=ub, b_ub=cap, bounds=var_bounds, method="highs")
        if not fit.success:
            raise RuntimeError(fit.message)
        endpoints.append(100 * sign * fit.fun)
    return np.array(endpoints)


def check(label, actual, expected):
    # Published endpoints are rounded to two decimal places.
    np.testing.assert_allclose(actual, expected, rtol=0, atol=0.0051)
    print(f"{label:32s} [{actual[0]:+.2f}, {actual[1]:+.2f}]")


# 1. Five-precinct PT table. Each post column has probability 1/2.
# These are the aggregate-only FD rows in the review.
rare = (
    20000, [5000, 8000, 11000, 14000, 17000],
    [40, 12000, 12000, 12000, 12000],
    [[50, 30], [9964, 9956], [8188, 8192], [6682, 6698], [5446, 5474]],
)
ordinary = (
    1000, [100, 300, 500, 700, 900],
    [500, 500, 500, 500, 500],
    [[682, 650], [606, 590], [530, 530], [454, 470], [378, 410]],
)

for name, data, within, pooled in (
    ("Rare", rare, [10.39, 12.25], [-3.96, -2.22]),
    ("Ordinary", ordinary, [-67.33, 43.78], [-80.56, 52.13]),
):
    for pt, expected in (("within", within), ("pooled", pooled)):
        check(f"{name}, {pt} PT",
              att_bounds(*data, comparison=pt, lce=("change",)), expected)


# 2. Three-precinct baseline check. One fixed aggregate total per election.
# All specifications maintain the SAME within-precinct mean PT assumption.
three = (
    800, [200, 400, 600],
    [4, 164, 324],             # 0.5%, 20.5%, 40.5%
    [[204], [124], [4]],       # 25.5%, 15.5%, 0.5%
)

for restrictions, expected in (
    ((),                         [-58.00, 27.33]),
    (("change",),                [-49.33,  7.33]),
    (("pre", "change"),          [  2.67,  7.33]),
    (("pre", "post"),            [  2.67,  7.33]),
    (("post", "change"),         [  2.67,  7.33]),
    (("pre", "post", "change"),  [  2.67,  7.33]),
):
    label = " + ".join(restrictions) or "No LCE"
    check(label, att_bounds(*three, comparison="within", lce=restrictions),
          expected)
