import sys; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from phase2e_code_review_round2 import LP
x = np.array([.25, .5, .75]); Y0 = np.array([1, 41, 81]) / 200; Y1 = np.array([51, 31, 1]) / 200; w = x / x.sum()
def ci(eps, valid=True, breakdown=False):
    L = LP(); m0 = {(g, i): L.var("m0", 0 if valid else None, 1 if valid else None) for g in "VN" for i in range(3)}
    d = {(g, i): L.var("d") for g in "VN" for i in range(3)}; ab = {g: (L.var("a"), L.var("b")) for g in "VN"}; ch = {g: (L.var("c"), L.var("e")) for g in "VN"}
    for i in range(3):
        L.eq({m0["V", i]: x[i], m0["N", i]: 1 - x[i]}, Y0[i]); L.eq({d["V", i]: x[i], d["N", i]: 1 - x[i]}, Y1[i] - Y0[i])
        for g in "VN":
            L.eq({d[g, i]: 1, ch[g][0]: -1, ch[g][1]: -x[i]}, 0)
            if valid: L.le({m0[g, i]: 1, d[g, i]: 1}, 1); L.le({m0[g, i]: -1, d[g, i]: -1}, 0)
            r = {m0[g, i]: 1, ab[g][0]: -1, ab[g][1]: -x[i]}; L.le(r, eps); L.le({k: -v for k, v in r.items()}, eps)
    o = {}
    for i in range(3): o[d["V", i]] = o.get(d["V", i], 0) + w[i]; o[d["N", i]] = o.get(d["N", i], 0) - w[i]
    return L.opt(o, 1), L.opt(o, -1)
def causal(c):
    aV = x * 200; aN = 200 - aV; L = LP()
    v0 = [L.var("v0", max(0, 200 * Y0[i] - aN[i]), min(aV[i], 200 * Y0[i])) for i in range(3)]
    v1 = [L.var("v1", max(0, 200 * Y1[i] - aN[i]), min(aV[i], 200 * Y1[i])) for i in range(3)]
    dw = {}; k0 = 0
    for i in range(3):
        # theta_i = v0/aV + ((Y1-v1)-(Y0-v0))/aN in [0,1]
        r = {v0[i]: 1 / aV[i] + 1 / aN[i], v1[i]: -1 / aN[i]}; k = 200 * (Y1[i] - Y0[i]) / aN[i]
        L.le(r, 1 - k); L.le({a: -b for a, b in r.items()}, k)
        dw[v1[i]] = w[i] * (1 / aV[i] + 1 / aN[i]); dw[v0[i]] = -w[i] * (1 / aV[i] + 1 / aN[i]); k0 -= w[i] * k
    L.le(dw, c[1] - k0); L.le({a: -b for a, b in dw.items()}, -(c[0] - k0))
    return L.opt(dw, 1) + k0, L.opt(dw, -1) + k0
for e in [0, .001, .01, .5]:
    c = ci(e); print("eps %-6g CI_LCE [%.4f, %.4f]  causal [%.4f, %.4f]" % ((e,) + c + causal(c)))
lo, hi = 0, .5
for _ in range(60):
    mid = (lo + hi) / 2
    if causal(ci(mid))[0] > 1e-12: lo = mid
    else: hi = mid
print("eps* = %.8f (1/600 = %.8f)" % (hi, 1 / 600))
