"""Reproduce the examples in the independent PHASE2e review. Run with Python, NumPy and SciPy."""
from fractions import Fraction
from itertools import product
from math import comb

import numpy as np
from scipy.optimize import linprog


def pair_loading():
    # Vertices and edge midpoints identify a homogeneous quadratic on the simplex.
    unit = np.eye(3)
    points = [*unit] + [(unit[i] + unit[j]) / 2 for i, j in [(0, 1), (0, 2), (1, 2), (0, 2)]]
    X = np.array(points)
    K = np.array([[0, 0, 1], [0, 0, 1], [-1, -1, 0]]) / 10
    mean_change = X @ K.T
    assert np.allclose(mean_change[:, 0], mean_change[:, 1])
    own_means = (X * mean_change).sum(axis=0) / X.sum(axis=0)
    assert np.isclose(own_means[0] - own_means[1], .0075)
    print("E7(d): zero same-composition contrasts; electorate DiD loading = 0.0075")


def shared_baseline():
    def feasible(x, increments):
        # One treated and one control rate at each time; observed aggregate is always 1/2.
        T = len(increments) + 1
        rows = np.kron(np.eye(T), np.array([x, 1 - x]))
        changes = np.kron(np.diff(np.eye(T), axis=0), np.eye(2))
        rhs = np.r_[np.full(T, .5), np.tile([1 - x, -x], len(increments)) * np.repeat(increments, 2)]
        return linprog(np.zeros(2*T), A_eq=np.vstack([rows, changes]),
                       b_eq=rhs, bounds=(0, 1)).success

    for x in [.4, .5, .6]:
        assert feasible(x, [1.25])
        assert not feasible(x, [1.25, 1.25])
    print("E8: each of two adjacent changes feasible at all three compositions; joint path infeasible")


def causal_paths():
    # Pure groups identify these realized rates. Regime 0 has no treatment anywhere.
    r = np.array([.35, .5, .5, .5])
    s = np.full(4, .5)
    r0 = np.array([.35, .5, .65, .80])
    b = (r0 - s) - (r0 - s)[1]
    A = abs(b[1] - b[0])
    assert np.all(np.abs(np.diff(b[1:])) <= A + 1e-12)
    assert np.isclose(r[-1] - r0[-1], -.30) and np.isclose(A, .15)
    print("E18(b): at horizon 2, D = 0, A = .15, admissible ATT = -.30")

    # E16(d): g_base=.5, control post rates s1+s2=1.2, individually in [.2,1].
    for restrict in [[0], [1], [0, 1]]:
        fit = linprog(np.zeros(2), A_eq=[[1, 1]], b_eq=[1.2],
                      A_ub=np.eye(2)[restrict], b_ub=[.5]*len(restrict), bounds=(.2, 1))
        assert fit.success == (len(restrict) == 1)
    print("E16(d): PT feasible at each post separately, jointly refuted")
    feasible_integer = any(Fraction(j, 1) == Fraction(1, 2) for j in [0, 1])
    assert not feasible_integer
    print("E16(b): one comparison person cannot have untreated rate 1/2")


def disconnected_rm():
    # x is the pre-period gap change; theta_PT = 1.1. Optimize .5 + 2*x - theta.
    branches = []
    obj = np.array([2., -1.])
    for sign in [-1, 1]:
        rows = [[-sign, 1], [-sign, -1]]
        fits = [linprog(k*obj, A_ub=rows, b_ub=[1.1, -1.1],
                        bounds=[(-.2,.2), (0,1)]) for k in [1,-1]]
        assert all(f.success for f in fits)
        branches.append((.5+fits[0].fun, .5-fits[1].fun))
    assert np.allclose(branches, [[-.9,-.7],[-.3,0]])
    print("E18: sharp RM projection has a hole:", np.round(branches, 6))


def finite_electorate():
    for m in [2, 3, 5, 10, 25, 100]:
        # y0,y1 are treated yes counts; m total yes at each election.
        possible = set()
        for y0, y1 in product(range(m + 1), repeat=2):
            theta = Fraction(y0, m) - Fraction(y1 - y0, m + 1)
            if 0 <= theta <= 1 and (m * theta).denominator == 1:
                possible.add(Fraction(y1, m) - theta)
        assert possible == {0}
        # Variables y0, y1, counterfactual treated rate.
        row = [-(1/m + 1/(m+1)), 1/(m+1), 1]
        obj = np.array([0, 1/m, -1])
        for sign in [-1, 1]:
            fit = linprog(sign*obj, A_eq=[row], b_eq=[0], bounds=[(0, m), (0, m), (0, 1)])
            assert fit.success and np.isclose(fit.fun, -1)
    print("Remark 7.1: continuous [-1,1], integer {0}, for m = 2,3,5,10,25,100")


def degree_three():
    pairs = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])
    opt = linprog(-np.ones(3), A_ub=pairs, b_ub=np.ones(3), bounds=(0, 1))
    integer_max = max(sum(v) for v in product(range(2), repeat=3) if np.all(pairs @ v <= 1))
    assert opt.success and np.isclose(-opt.fun, 1.5) and integer_max == 1
    print("E27 extension: three-pair star, maximum 1.5 continuous, 1 integer")


def finite_sample():
    n = 20
    bound = np.sqrt(np.log(1/.025)/(2*n))
    values = 2*np.arange(n+1)/n - 1
    weights = np.array([comb(n, i)/2**n for i in range(n+1)])
    exact = [weights[np.abs(values) <= k*bound].sum() for k in [1, 2]]
    rng = np.random.default_rng(20261123)
    sample = 2*rng.binomial(n, .5, 100000)/n - 1
    sim = [np.mean(np.abs(sample) <= k*bound) for k in [1, 2]]
    assert exact[0] < .95 <= exact[1] and np.allclose(sim, exact, atol=.004)
    print("E23: coverage, level/change radii:", np.round(exact, 6), "MC:", np.round(sim, 6))


def rm_branches(groups, votes, reinforce=()):
    # Independent reduction for two groups, one pre contrast, one post. All functions
    # below are affine in (g_early, g_base, g_post, untreated_treated_post_rate).
    counts = np.asarray(groups)
    votes = np.asarray(votes)
    N = counts.sum(axis=1)
    r, s = counts.sum(axis=0)
    total = r + s
    def gap(t, y):
        return y/r - (votes[t].sum()-y)/s
    intervals = [(gap(t, np.maximum(0, counts[:,0]+Y-N).sum()),
                  gap(t, np.minimum(counts[:,0], Y).sum())) for t, Y in enumerate(votes)]
    bounds = intervals + [(0, 1)]
    x = counts[:,0]/N
    w = N / total
    xbar = r/total
    var = np.dot(w, (x-xbar)**2)
    extra = []
    for t in reinforce:
        cov = np.dot(w*(x-xbar), (votes[t]-votes[t-1])/N)
        lo, hi = sorted([cov/var, cov/(xbar*(1-xbar))])
        vec = np.eye(4)[t]-np.eye(4)[t-1]
        extra.extend([(vec, hi), (-vec, -lo)])
    post = votes[2].sum()/total
    target = np.array([0, 0, s/total, -1])
    # b = theta - g_base - B_s,post ; B_s,post = observed_post - (r/total)*g_post
    bcoef = np.array([0, -1, r/total, 1])
    pre = np.array([-1, 1, 0, 0])
    branches = []
    for sign in [1, -1]:
        constraints = extra + [(bcoef-sign*pre, post), (-bcoef-sign*pre, -post)]
        A, rhs = zip(*constraints)
        fits = [linprog(k*target, A_ub=A, b_ub=rhs, bounds=bounds) for k in [1, -1]]
        if fits[0].success:
            assert fits[1].success
            branches.append((post+fits[0].fun, post-fits[1].fun))
        else:
            assert fits[0].status == fits[1].status == 2
    merged = []
    for lo, hi in sorted(branches):
        if merged and lo <= merged[-1][1] + 1e-9:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def worked_examples():
    n21 = rm_branches([[20, 80], [60, 40]], [[30, 45], [32, 48], [30, 60]])
    assert np.allclose(n21, [[-.75, 1.]])
    groups = [[10,190], [50,150], [150,50], [190,10]]
    votes = [[60,62,61,56], [62,66,63,61], [74,86,113,120]]
    for rungs, expected in [((), [-.568,.732]), ((1,), [-.416,.732]), ((1,2), [.135,.278])]:
        intervals = rm_branches(groups, votes, rungs)
        assert len(intervals) == 1 and np.allclose(intervals[0], expected, atol=.00051)
        print("Four-precinct RM, CR periods", rungs, ":", np.round(intervals, 6))
    mixed = rm_branches([[99,1], [1,99]], [[50,50], [50,50], [70,50]])
    assert len(mixed) == 1 and mixed[0][0] > 0
    print("N21 RM:", n21, "; margins-only, two mixed precincts:", np.round(mixed, 6))


def integer_worked_example():
    # With two groups of 400, gap in units of 1/400 is 2*treated_yes - total_yes.
    a = np.array([10, 50, 150, 190])
    votes = np.array([[60,62,61,56], [62,66,63,61], [74,86,113,120]])
    yy = [np.arange(np.maximum(0,a+v-200).sum(), np.minimum(a,v).sum()+1) for v in votes]
    gaps = [2*y-v.sum() for y,v in zip(yy, votes)]
    x = a/200
    var = np.var(x)
    rows = []
    for change in np.diff(votes, axis=0):
        cov = np.mean((x-.5)*change/200)
        rows.append(400*np.sort([cov/var,cov/.25]))
    for cr, expected in [((), (-227,293)), ((0,), (-166,293)), ((0,1), (56,110))]:
        effects = set()
        for base in gaps[1]:
            pre = base - gaps[0]
            if 0 in cr:
                pre = pre[(pre >= rows[0][0]-1e-9) & (pre <= rows[0][1]+1e-9)]
            if not len(pre):
                continue
            allowance = max(abs(pre))
            for y, post in zip(yy[2], gaps[2]):
                d = post - base
                if 1 in cr and not rows[1][0]-1e-9 <= d <= rows[1][1]+1e-9:
                    continue
                effects.update(range(max(y-400,d-allowance), min(y,d+allowance)+1))
        assert effects == set(range(expected[0], expected[1]+1))
        print("Four-precinct integer RM, CR changes", cr, ": endpoints", expected, "/ 400")


if __name__ == "__main__":
    pair_loading()
    shared_baseline()
    causal_paths()
    disconnected_rm()
    finite_electorate()
    degree_three()
    finite_sample()
    worked_examples()
    integer_worked_example()
