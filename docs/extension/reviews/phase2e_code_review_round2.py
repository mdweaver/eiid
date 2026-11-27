"""Round-2 confirmation (code review of PHASE2e): independent LPs, written from the theory text only,
for Prop. E30(c) (three-precinct table, epsilon*) and Prop. E35's worked example (realised WPT-tau).
Usage: python3.12 -I phase2e_code_review_round2.py
"""
import numpy as np
from scipy.optimize import linprog


class LP:
    def __init__(self):
        self.names, self.lb, self.ub, self.Aeq, self.beq, self.Aub, self.bub = [], [], [], [], [], [], []

    def var(self, name, lo=None, hi=None):
        self.names.append(name); self.lb.append(lo); self.ub.append(hi); return len(self.names) - 1

    def eq(self, coef, rhs): self.Aeq.append(coef); self.beq.append(rhs)

    def le(self, coef, rhs): self.Aub.append(coef); self.bub.append(rhs)

    def _mat(self, rows):
        M = np.zeros((len(rows), len(self.names)))
        for k, r in enumerate(rows):
            for j, v in r.items(): M[k, j] += v
        return M

    def opt(self, obj, sense=1):
        c = np.zeros(len(self.names))
        for j, v in obj.items(): c[j] = sense * v
        r = linprog(c, A_ub=self._mat(self.Aub) if self.Aub else None, b_ub=self.bub or None,
                    A_eq=self._mat(self.Aeq) if self.Aeq else None, b_eq=self.beq or None,
                    bounds=list(zip(self.lb, self.ub)), method="highs")
        return sense * r.fun if r.status == 0 else None


def e30(change, base, post, eps=0.0, breakdown=False):
    x = np.array([.25, .5, .75]); Y0 = np.array([.005, .205, .405]); Y1 = np.array([.255, .155, .005])
    w = x / x.sum()
    L = LP()
    m0 = {(g, i): L.var("m0%s%d" % (g, i), 0, 1) for g in "VN" for i in range(3)}
    d = {(g, i): L.var("d%s%d" % (g, i)) for g in "VN" for i in range(3)}
    ab = {g: (L.var("a" + g), L.var("b" + g)) for g in "VN"}
    ch = {g: (L.var("al" + g), L.var("be" + g)) for g in "VN"}
    po = {g: (L.var("c" + g), L.var("e" + g)) for g in "VN"}
    E = L.var("eps", 0, None) if breakdown else None
    for i in range(3):
        L.eq({m0["V", i]: x[i], m0["N", i]: 1 - x[i]}, Y0[i])
        L.eq({m0["V", i]: x[i], m0["N", i]: 1 - x[i], d["V", i]: x[i], d["N", i]: 1 - x[i]}, Y1[i])
        for g in "VN":                                       # post rates in [0, 1]
            L.le({m0[g, i]: 1, d[g, i]: 1}, 1); L.le({m0[g, i]: -1, d[g, i]: -1}, 0)
            if change: L.eq({d[g, i]: 1, ch[g][0]: -1, ch[g][1]: -x[i]}, 0)
            if post: L.eq({m0[g, i]: 1, d[g, i]: 1, po[g][0]: -1, po[g][1]: -x[i]}, 0)
            if base:
                r = {m0[g, i]: 1, ab[g][0]: -1, ab[g][1]: -x[i]}
                if breakdown:
                    L.le({**r, E: -1}, 0); L.le({**{k: -v for k, v in r.items()}, E: -1}, 0)
                else:
                    L.le(r, eps); L.le({k: -v for k, v in r.items()}, eps)
        L.le({m0["V", i]: 1, d["N", i]: 1}, 1); L.le({m0["V", i]: -1, d["N", i]: -1}, 0)   # WPT counterfactual a rate
    att = {}
    for i in range(3):
        att[d["V", i]] = att.get(d["V", i], 0) + w[i]; att[d["N", i]] = att.get(d["N", i], 0) - w[i]
    if breakdown:
        L.le(att, 0); return L.opt({E: 1})
    return L.opt(att, 1), L.opt(att, -1)


def e35(M=0.0, audit=False, breakdown=False, dw=False):
    aV = np.array([10, 50, 150, 190.]); aN = np.array([190, 150, 50, 10.])
    Y0 = np.array([62, 66, 63, 61.]); Y1 = np.array([74, 86, 113, 120.])
    w = aV / aV.sum(); tau = np.sqrt(1 / aV + 1 / aN)
    L = LP()
    v0 = [L.var("v0%d" % i, max(0, Y0[i] - aN[i]), min(aV[i], Y0[i])) for i in range(4)]
    v1 = [L.var("v1%d" % i, max(0, Y1[i] - aN[i]), min(aV[i], Y1[i])) for i in range(4)]
    th = [L.var("th%d" % i, 0, 1) for i in range(4)]
    Mv = L.var("M", 0, None) if breakdown else None
    if audit:   # non-veteran yes: P3 20 -> 22, P4 4 -> 5
        for i, (a, b) in {2: (20, 22), 3: (4, 5)}.items():
            L.eq({v0[i]: 1}, Y0[i] - a); L.eq({v1[i]: 1}, Y1[i] - b)
    for i in range(4):
        # imputed = v0/aV + ((Y1 - v1) - (Y0 - v0))/aN ;  |theta - imputed| <= M tau
        r = {th[i]: 1, v0[i]: -(1 / aV[i] + 1 / aN[i]), v1[i]: 1 / aN[i]}
        k = (Y1[i] - Y0[i]) / aN[i]
        if dw: continue
        if breakdown:
            L.le({**r, Mv: -tau[i]}, k); L.le({**{a: -b for a, b in r.items()}, Mv: -tau[i]}, -k)
        else:
            L.le(r, k + M * tau[i]); L.le({a: -b for a, b in r.items()}, -k + M * tau[i])
    if dw:   # D_W = sum_i w_i (dB_iV - dB_iN), no rate limits on the imputation
        obj = {}
        for i in range(4):
            obj[v1[i]] = w[i] * (1 / aV[i] + 1 / aN[i]); obj[v0[i]] = -w[i] * (1 / aV[i] + 1 / aN[i])
        c0 = -sum(w[i] * (Y1[i] - Y0[i]) / aN[i] for i in range(4))
        return L.opt(obj, 1) + c0, L.opt(obj, -1) + c0
    att = {}
    for i in range(4):
        att[v1[i]] = w[i] / aV[i]; att[th[i]] = -w[i]
    if breakdown:
        L.le(att, 0); return L.opt({Mv: 1})
    return L.opt(att, 1), L.opt(att, -1)


def e28(size, vet, pre, post_states, kind):
    """Prop. E28(c) with a known law: post totals equally likely per state, baseline fixed; Aumann mean of the
    per-state Duncan-Davis segments (same line slope, so averaged ends); FD-LCE on changes; WPT or PPT."""
    x = np.array(vet, float) / size; Y0 = np.array(pre, float) / size; S = np.array(post_states, float) / size
    p = len(x); wV = x / x.sum(); wN = (1 - x) / (1 - x).sum()
    dd = lambda Y, xi: (max(0, (Y - (1 - xi)) / xi), min(1, Y / xi))
    L = LP()
    mV0 = [L.var("mV0%d" % i, *dd(Y0[i], x[i])) for i in range(p)]
    mV1 = [L.var("mV1%d" % i, np.mean([dd(y, x[i])[0] for y in S[i]]), np.mean([dd(y, x[i])[1] for y in S[i]])) for i in range(p)]
    al = {g: (L.var("al" + g), L.var("be" + g)) for g in "VN"}
    att, c0 = {}, 0.0
    for i in range(p):
        Y1 = S[i].mean()
        # d_V = mV1 - mV0 ; d_N = (Y1 - x mV1)/(1-x) - (Y0 - x mV0)/(1-x)
        dV = {mV1[i]: 1, mV0[i]: -1}; dN = {mV1[i]: -x[i] / (1 - x[i]), mV0[i]: x[i] / (1 - x[i])}; kN = (Y1 - Y0[i]) / (1 - x[i])
        L.eq({**dV, al["V"][0]: -1, al["V"][1]: -x[i]}, 0)
        L.eq({**dN, al["N"][0]: -1, al["N"][1]: -x[i]}, -kN)
        if kind == "wpt":   # m_V0 + d_N in [0, 1]; ATT += w (d_V - d_N)
            r = {mV0[i]: 1 + x[i] / (1 - x[i]), mV1[i]: -x[i] / (1 - x[i])}
            L.le(r, 1 - kN); L.le({a: -b for a, b in r.items()}, kN)
            for a, b in dV.items(): att[a] = att.get(a, 0) + wV[i] * b
            for a, b in dN.items(): att[a] = att.get(a, 0) - wV[i] * b
            c0 -= wV[i] * kN
        else:               # ATT = sum wV d_V - sum wN d_N
            for a, b in dV.items(): att[a] = att.get(a, 0) + wV[i] * b
            for a, b in dN.items(): att[a] = att.get(a, 0) - wN[i] * b
            c0 -= wN[i] * kN
    if kind == "ppt":       # p0 = sum wV mV0 + sum wN d_N in [0, 1]
        r, k = {}, 0.0
        for i in range(p):
            Y1 = S[i].mean(); r[mV0[i]] = r.get(mV0[i], 0) + wV[i] + wN[i] * x[i] / (1 - x[i])
            r[mV1[i]] = r.get(mV1[i], 0) - wN[i] * x[i] / (1 - x[i]); k += wN[i] * (Y1 - Y0[i]) / (1 - x[i])
        L.le(r, 1 - k); L.le({a: -b for a, b in r.items()}, k)
    return 100 * (L.opt(att, 1) + c0), 100 * (L.opt(att, -1) + c0)


if __name__ == "__main__":
    f = lambda t: "[%.6f, %.6f]" % t
    print("E30(c) three-precinct table (expected ATT), theory: no LCE [-29/50, 41/150]; change [-37/75, 11/150]; any two [2/75, 11/150]")
    print("  no LCE      ", f(e30(False, False, False)), " target", f((-29 / 50, 41 / 150)))
    print("  change only ", f(e30(True, False, False)), " target", f((-37 / 75, 11 / 150)))
    for lab, a in [("base+change", (True, True, False)), ("base+post", (False, True, True)), ("post+change", (True, False, True)),
                   ("all three", (True, True, True))]:
        print("  %-12s" % lab, f(e30(*a)), " target", f((2 / 75, 11 / 150)))
    print("  epsilon* (base+change) %.8f   target 1/600 = %.8f" % (e30(True, True, False, breakdown=True), 1 / 600))
    print("  eps = 1/600 - 1e-6 lower %.3e ; eps = 1/600 lower %.3e" % (e30(True, True, False, 1 / 600 - 1e-6)[0], e30(True, True, False, 1 / 600)[0]))
    print("  eps = 1/2  ", f(e30(True, True, False, .5)), " = change only", f(e30(True, False, False)))
    print("E35 worked example (realised WPT-tau, Sec. 9.2 data, base -> post)")
    print("  margins only, M=0 ", f(e35()), " theory [-0.5675, 0.7322]")
    print("  audited, M=0      ", f(e35(0, True)), " theory [0.0525, 0.3522]")
    print("  audited, M=0.1    ", f(e35(.1, True)), " theory [0.0310, 0.3740]")
    print("  audited, M*       %.6f   theory 0.2438" % e35(audit=True, breakdown=True))
    print("  audited, D_W count set", f(e35(audit=True, dw=True)), " theory [-0.0087, 0.3772]")
    print("E28 five-precinct table (points; theory: rare W [10.39, 12.25] P [-3.96, -2.22]; ordinary W [-67.33, 43.78] P [-80.56, 52.13])")
    rare = (20000, [5000, 8000, 11000, 14000, 17000], [40, 12000, 12000, 12000, 12000],
            [[50, 30], [9964, 9956], [8188, 8192], [6682, 6698], [5446, 5474]])
    ordinary = (1000, [100, 300, 500, 700, 900], [500] * 5, [[682, 650], [606, 590], [530, 530], [454, 470], [378, 410]])
    for nm, d in [("rare", rare), ("ordinary", ordinary)]:
        print("  %-9s WPT [%.4f, %.4f]  PPT [%.4f, %.4f]" % ((nm,) + e28(*d, "wpt") + e28(*d, "ppt")))
