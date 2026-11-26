#!/usr/bin/env python3
"""Independent oracle for the code review of eiChange() / eiATT() (PHASE2e, step 2 DiD).

Written from docs/extension/PHASE2e_did_general.md (2e) and 2c Definition 1, without the
package's row builders. Every row is built here from its definition as an affine function of
the stacked T-election table n[t, i, r, c] (Lemma E1), with realised rates beta = n / a:

  * delta_B (W1, one-vs-rest): sum_{s != r} E_p[X_r X_s (beta_r - beta_s)];
  * delta_B (pair r, s):        E_p[X_r X_s (beta_r - beta_s)];
  * delta_W (W1): mu_r Cov_{r}(beta_r, X_r) + (1 - mu_r) Cov_{not r}(beta_notr, X_r);
  * delta_W (W2): sum_g mu_g Cov_g(beta_g, X_r);
  * a row "on changes" is the association at election t minus the one at u (Lemma E4(a));
  * order_change: Delta beta_{irc} - Delta beta_{isc} (s a row, or the rest) in every unit;
  * contrast_change: 2e Prop. E6 (mix: lambda * Dbar >= 0; rate: rho >= 0);
  * CR: {dB >= 0, dW >= 0} U {dB <= 0, dW <= 0} (a union; slack relaxes each side).

For each instance:
  * brute force (tiny instances): every integer table of every unit-election, product over
    unit-elections, exact integer min / max of each target subject to the rows;
  * continuous ends: scipy.optimize.linprog (HiGHS) per union branch;
  * causal sets (Thm E16): an LP in (n, theta_t, S_t, b_t) per (trend branch x count branch),
    with the identity of Lemma E14, rate limits, |S| <= Lambda, the dial, RM (Def. 5.2, joint,
    consecutive differences, union of 2L polyhedra), SD (Def. 5.3), MTR/ND/PS (Prop. E21).

Usage:
  python3.12 -I phase2e_code_review_oracle.py gen   INSTANCES.json [seed]
  python3.12 -I phase2e_code_review_oracle.py check INSTANCES.json RESULTS.json
The R side (phase2e_code_review_checks.R) runs the package on each instance and writes RESULTS.json.
"""
import itertools, json, random, sys
import numpy as np
from scipy.optimize import linprog

TOL = 1e-9

# ---------------------------------------------------------------------------
# design and unit tables
# ---------------------------------------------------------------------------
def unit_tables(a, b):
    R, C = len(a), len(b)
    out = []
    def rec(r, colrem, acc):
        if r == R - 1:
            if sum(colrem) == a[r] and all(v >= 0 for v in colrem):
                out.append(acc + [list(colrem)])
            return
        def comp(c, rem, cur):
            if c == C - 1:
                if rem <= colrem[c]:
                    yield cur + [rem]
                return
            for v in range(0, min(rem, colrem[c]) + 1):
                yield from comp(c + 1, rem - v, cur + [v])
        for row in comp(0, a[r], []):
            rec(r + 1, [colrem[c] - row[c] for c in range(C)], acc + [row])
    rec(0, list(b), [])
    return np.array(out, dtype=float)


class Inst:
    def __init__(self, d):
        self.d = d
        self.A = np.array(d["A"], dtype=float)            # p x R
        self.B = np.array(d["B"], dtype=float)            # T x p x C
        self.T, self.p, self.C = self.B.shape
        self.R = self.A.shape[1]
        self.N = self.A.sum(1)
        self.shape = (self.T, self.p, self.R, self.C)
        self.nv = int(np.prod(self.shape))
        self.pw = self.N / self.N.sum()
        self.X = self.A / self.N[:, None]

    def zero(self):
        return np.zeros(self.shape)

    # ---- targets ---------------------------------------------------------------
    def rate(self, t, r, c, S=None):
        S = range(self.p) if S is None else S
        v = self.zero(); Ar = sum(self.A[i, r] for i in S)
        for i in S: v[t, i, r, c] = 1.0 / Ar
        return v, 0.0

    def gap(self, t, r, s, c, S=None):
        S = list(range(self.p)) if S is None else S
        v, k = self.rate(t, r, c, S)
        if s == "rest":
            den = sum(self.N[i] for i in S) - sum(self.A[i, r] for i in S)
            for i in S:
                for q in range(self.R):
                    if q != r: v[t, i, q, c] -= 1.0 / den
        else:
            w, _ = self.rate(t, s, c, S); v = v - w
        return v, 0.0

    def did(self, r, s, c, u, t, S=None):
        a, _ = self.gap(t, r, s, c, S); b, _ = self.gap(u, r, s, c, S)
        return a - b, 0.0

    def bins(self, c, u, t, L, H):
        """2e Prop. E6: Dbar coefficients, lambda, Dtau_H - Dtau_L (R = 2)."""
        v = self.zero(); xb = {}; dtau = {}
        for name, b in (("L", L), ("H", H)):
            A1 = sum(self.A[i, 0] for i in b); A2 = sum(self.A[i, 1] for i in b); Nb = sum(self.N[i] for i in b)
            for i in b:
                v[t, i, 0, c] += 0.5 / A1; v[u, i, 0, c] -= 0.5 / A1
                v[t, i, 1, c] -= 0.5 / A2; v[u, i, 1, c] += 0.5 / A2
            xb[name] = A1 / Nb
            dtau[name] = sum(self.B[t, i, c] - self.B[u, i, c] for i in b) / Nb
        return v, xb["H"] - xb["L"], dtau["H"] - dtau["L"]

    # ---- associations at one election (2c Definition 1) ------------------------
    def _beta_coef(self, W, t, c):
        """sum_i sum_s W[i, s] beta_{isct} as coefficients on n (beta = n / a)."""
        v = self.zero()
        for i in range(self.p):
            for s in range(self.R):
                if self.A[i, s] > 0 and W[i, s] != 0:
                    v[t, i, s, c] += W[i, s] / self.A[i, s]
        return v

    def assoc_W(self, typ, ver, r, other=None):
        p, X, R = self.pw, self.X, self.R
        mu = (p[:, None] * X).sum(0); xr = X[:, r]
        W = np.zeros((self.p, R))
        if typ == "dB":
            oth = [other] if ver == "pair" else [s for s in range(R) if s != r]
            for s in oth:
                w = p * xr * X[:, s]; W[:, r] += w; W[:, s] -= w
        elif typ == "dW" and ver == "W1":
            xbr = (p * xr * xr).sum() / mu[r]
            xbn = (p * (1 - xr) * xr).sum() / (1 - mu[r])
            W[:, r] = p * xr * (xr - xbr)
            for s in range(R):
                if s != r: W[:, s] = p * X[:, s] * (xr - xbn)
        elif typ == "dW" and ver == "W2":
            for g in range(R):
                if mu[g] <= 0: continue
                xb = (p * X[:, g] * xr).sum() / mu[g]
                W[:, g] = p * X[:, g] * (xr - xb)
        else:
            raise ValueError(typ + ver)
        return W

    def scale(self, typ, ver, r, other=None):
        p, X = self.pw, self.X; xr = X[:, r]; mu = (p[:, None] * X).sum(0)
        def wvar(w):
            w = w / w.sum(); return (w * xr * xr).sum() - (w * xr).sum() ** 2
        if typ == "dB":
            oth = [other] if ver == "pair" else [s for s in range(self.R) if s != r]
            return float((p * xr * X[:, oth].sum(1)).sum())
        if ver == "W1":
            return float(mu[r] * wvar(p * xr) + (1 - mu[r]) * wvar(p * (1 - xr)))
        return float(sum(mu[g] * wvar(p * X[:, g]) for g in range(self.R) if mu[g] > 0))

    def assoc_change(self, typ, ver, r, c, u, t, other=None):
        W = self.assoc_W(typ, ver, r, other)
        return self._beta_coef(W, t, c) - self._beta_coef(W, u, c)

    # ---- rows: list of (coef, lo, hi) with lo <= coef . n <= hi ------------------
    def rows_of(self, spec):
        """returns ('common', [rows]) or ('union', [[rows plus], [rows minus]])"""
        typ = spec["type"]
        if typ in ("sign_change", "cr_change"):
            r = spec["row"]; c = spec["col"]; u = spec["from"]; t = spec["to"]
            sl = spec.get("slack", 0.0); oth = spec.get("other")
            if oth is None and self.R == 2: oth = 1 - r
            if typ == "sign_change":
                st, ver, sg = spec["stype"], spec.get("version", "W1"), spec.get("sign", 1)
                v = self.assoc_change(st, ver, r, c, u, t, oth)
                eps = sl * self.scale(st, ver, r, oth)
                return "common", [(v, -eps, np.inf) if sg > 0 else (v, -np.inf, eps)]
            ver = spec.get("version", "W1")
            vB = "pair" if ver == "pair" else "W1"; vW = "W1" if ver == "W1" else "W2"
            fB = self.assoc_change("dB", vB, r, c, u, t, oth); fW = self.assoc_change("dW", vW, r, c, u, t, oth)
            eB = sl * self.scale("dB", vB, r, oth); eW = sl * self.scale("dW", vW, r, oth)
            return "union", [[(fB, -eB, np.inf), (fW, -eW, np.inf)], [(fB, -np.inf, eB), (fW, -np.inf, eW)]]
        if typ == "order_change":
            r = spec["row"]; c = spec["col"]; u = spec["from"]; t = spec["to"]; sg = spec.get("sign", 1)
            q = spec.get("over", "rest"); sl = spec.get("slack", 0.0)
            if self.R == 2 and q != "rest": q = "rest"
            out = []
            for i in range(self.p):
                a = self.A[i, r]
                if q == "rest":
                    if a == 0 or a == self.N[i]: continue
                    v = self.zero()
                    v[t, i, r, c] += 1 / a; v[u, i, r, c] -= 1 / a
                    for s in range(self.R):
                        if s != r:
                            v[t, i, s, c] -= 1 / (self.N[i] - a); v[u, i, s, c] += 1 / (self.N[i] - a)
                else:
                    if a == 0 or self.A[i, q] == 0: continue
                    v = self.zero(); aq = self.A[i, q]
                    v[t, i, r, c] += 1 / a; v[u, i, r, c] -= 1 / a; v[t, i, q, c] -= 1 / aq; v[u, i, q, c] += 1 / aq
                out.append((v, -sl, np.inf) if sg > 0 else (v, -np.inf, sl))
            return "common", out
        if typ == "contrast_change":
            v, lam, dtau = self.bins(spec["col"], spec["from"], spec["to"], spec["low"], spec["high"])
            bel = spec.get("beliefs", ["mix", "rate"]); out = []
            if "mix" in bel: out.append((v, 0.0, np.inf))
            if "rate" in bel: out.append((v, -np.inf, dtau / lam))
            return "common", out
        if typ == "known":
            v = self.zero(); v[spec["election"], spec["unit"], spec["row"], spec["col"]] = 1
            return "common", [(v, spec["value"], spec["value"])]
        if typ == "total":
            v = self.zero()
            for i in spec["units"]: v[spec["election"], i, spec["row"], spec["col"]] = 1
            return "common", [(v, spec.get("lower", -np.inf), spec.get("upper", np.inf))]
        if typ == "linear":
            v = self.zero()
            for (i, r, c, t, cf) in spec["cells"]: v[t, i, r, c] += cf
            return "common", [(v, spec.get("lower", -np.inf), spec.get("upper", np.inf))]
        raise ValueError(typ)

    def branches(self):
        common, unions = [], []
        for s in self.d.get("rows", []):
            k, rows = self.rows_of(s)
            if k == "common": common += rows
            else: unions.append(rows)
        out = []
        for combo in itertools.product(*[range(2)] * len(unions)):
            rows = list(common)
            for j, side in enumerate(combo): rows += unions[j][side]
            out.append(rows)
        return out

    def target(self, tg):
        k = tg["kind"]
        if k == "did": return self.did(tg["row"], tg["versus"], tg["col"], tg["from"], tg["to"])[0]
        if k == "gap_at": return self.gap(tg["election"], tg["row"], tg["versus"], tg["col"])[0]
        if k == "rate_at": return self.rate(tg["election"], tg["row"], tg["col"])[0]
        if k == "bin_did": return self.bins(tg["col"], tg["from"], tg["to"], tg["low"], tg["high"])[0]
        raise ValueError(k)

    # ---- margins ---------------------------------------------------------------
    def margin_eq(self):
        rows, rhs = [], []
        T, p, R, C = self.shape
        for t in range(T):
            for i in range(p):
                for r in range(R):
                    v = self.zero(); v[t, i, r, :] = 1; rows.append(v.ravel()); rhs.append(self.A[i, r])
                for c in range(C):
                    v = self.zero(); v[t, i, :, c] = 1; rows.append(v.ravel()); rhs.append(self.B[t, i, c])
        return np.array(rows), np.array(rhs)

    # ---- continuous LP ------------------------------------------------------------
    def lp_range(self, obj, rows, extra_nx=0, extra=None, xbounds=None, obj_x=None):
        """min / max of obj . n (+ obj_x . x) over margins, rows (on n) and extra rows (on (n, x))."""
        Aeq, beq = self.margin_eq()
        nx = extra_nx
        Aeq = np.hstack([Aeq, np.zeros((Aeq.shape[0], nx))])
        Aub, bub = [], []
        def add(v, lo, hi):
            if np.isfinite(hi): Aub.append(v); bub.append(hi)
            if np.isfinite(lo): Aub.append(-v); bub.append(-lo)
        eq_extra, beq_extra = [], []
        for (v, lo, hi) in rows:
            vv = np.concatenate([v.ravel(), np.zeros(nx)])
            if np.isfinite(lo) and lo == hi: eq_extra.append(vv); beq_extra.append(lo)
            else: add(vv, lo, hi)
        for (vn, vx, lo, hi) in (extra or []):
            vv = np.concatenate([vn.ravel() if vn is not None else np.zeros(self.nv), vx])
            if np.isfinite(lo) and lo == hi: eq_extra.append(vv); beq_extra.append(lo)
            else: add(vv, lo, hi)
        if eq_extra:
            Aeq = np.vstack([Aeq, np.array(eq_extra)]); beq = np.concatenate([beq, beq_extra])
        bnds = [(0, None)] * self.nv + (xbounds or [])
        cvec = np.concatenate([obj.ravel(), obj_x if obj_x is not None else np.zeros(nx)])
        res = []
        for sg in (1, -1):
            r = linprog(sg * cvec, A_ub=np.array(Aub) if Aub else None, b_ub=np.array(bub) if bub else None,
                        A_eq=Aeq, b_eq=beq, bounds=bnds, method="highs")
            if r.status == 2: return None
            if r.status == 3: res.append(-np.inf * sg)  # unbounded
            elif r.status != 0: raise RuntimeError(r.message)
            else: res.append(sg * r.fun)
        return [res[0], res[1]]

    # ---- brute force ---------------------------------------------------------------
    def brute(self, targets):
        T, p, R, C = self.shape
        tabs = [[unit_tables(self.A[i].astype(int).tolist(), self.B[t, i].astype(int).tolist()) for i in range(p)]
                for t in range(T)]
        sizes = [len(tabs[t][i]) for t in range(T) for i in range(p)]
        M = int(np.prod(sizes))
        if M > 3_000_000: return None
        grids = np.meshgrid(*[np.arange(s) for s in sizes], indexing="ij")
        n = np.zeros((M, T, p, R, C))
        k = 0
        for t in range(T):
            for i in range(p):
                n[:, t, i] = tabs[t][i][grids[k].ravel()]; k += 1
        nf = n.reshape(M, -1)
        def ok_rows(rows):
            ok = np.ones(M, bool)
            for (v, lo, hi) in rows:
                val = nf @ v.ravel()
                sc = 1e-9 * max(1.0, np.abs(v).sum() * self.N.max())
                ok &= (val >= lo - sc) & (val <= hi + sc)
            return ok
        keep = np.zeros(M, bool)
        for rows in self.branches(): keep |= ok_rows(rows)
        out = []
        for tg in targets:
            vals = nf[keep] @ self.target(tg).ravel()
            out.append(None if len(vals) == 0 else [float(vals.min()), float(vals.max())])
        return out

    # ---- continuous ranges per branch ---------------------------------------------------
    def cont(self, targets):
        out = []
        brs = self.branches()
        for tg in targets:
            ivs = []
            for rows in brs:
                iv = self.lp_range(self.target(tg), rows)
                if iv is not None: ivs.append(iv)
            out.append(ivs)
        return out

    # ---- causal LP (Thm E16) -------------------------------------------------------------
    def causal(self, cs):
        """cs: treated r, comparison s ('rest' | int | None), col c, base b0, post (list of t),
        trend {kind, M, Mbar, mtr, nd, ps}, Lambda, rate_limits. Returns per post t a list of
        branch intervals (None for infeasible)."""
        r, s, c, b0 = cs["treated"], cs["comparison"], cs["col"], cs["base"]
        tr = cs["trend"]; Lam = cs["Lambda"]; rl = cs.get("rate_limits", True)
        posts = list(range(b0 + 1, self.T))
        K = len(posts); nx = 3 * K
        ix = lambda k, w: 3 * k + {"theta": 0, "S": 1, "b": 2}[w]
        def ev(pairs):
            v = np.zeros(nx)
            for (j, a) in pairs: v[j] += a
            return v
        Br = lambda t: self.rate(t, r, c)[0]
        if s is not None:
            Bs = (lambda t: self.gap(t, r, s, c)[0] * 0 + (self.rate(t, s, c)[0] if s != "rest" else self._rest_rate(t, r, c)))
            g = lambda t: Br(t) - Bs(t)
            thPT = lambda t: Br(b0) + Bs(t) - Bs(b0)
        xb = []
        for k in range(K):
            xb += [(0, 1) if rl else (None, None), (-Lam, Lam) if s is not None else (0, 0), (None, None)]
        common = []
        for k, t in enumerate(posts):
            if s is not None:
                common.append((thPT(t), ev([(ix(k, "b"), 1), (ix(k, "theta"), -1), (ix(k, "S"), -1)]), 0.0, 0.0))
                if rl: common.append((Bs(t), ev([(ix(k, "S"), -1)]), 0.0, 1.0))
        brs = []
        kind = tr["kind"]
        if kind == "none": brs = [[]]
        elif kind == "dial":
            brs = [[(None, ev([(ix(k, "b"), 1)]), -tr["M"], tr["M"]) for k in range(K)]]
        elif kind == "mtr":
            rr = []
            for k, t in enumerate(posts):
                th = ev([(ix(k, "theta"), 1)])
                if tr["mtr"]: rr.append((-Br(t), th, -np.inf, 0.0))
                if tr["nd"]: rr.append((-Br(b0), th, 0.0, np.inf))
                if tr["ps"]: rr.append((-Bs(t), th + ev([(ix(k, "S"), 1)]), 0.0, np.inf))
            brs = [rr]
        elif kind == "rm":
            pre = list(range(0, b0))                    # pairs (e, e+1), e = 0..b0-1, last touches the base
            for e in pre:
                for sg in (1, -1):
                    d = sg * tr["Mbar"] * (g(e + 1) - g(e))
                    rr = []
                    for k in range(K):
                        dd = ev([(ix(k, "b"), 1)] + ([(ix(k - 1, "b"), -1)] if k > 0 else []))
                        rr.append((d, -dd, 0.0, np.inf)); rr.append((d, dd, 0.0, np.inf))
                    brs.append(rr)
        elif kind == "sd":
            # delta_e: e < b0: D_e = g_e - g_b0; e == b0: 0; e > b0: b
            def dl(e):
                if e < b0: return (g(e) - g(b0), np.zeros(nx))
                if e == b0: return (np.zeros(self.shape), np.zeros(nx))
                return (np.zeros(self.shape), ev([(ix(posts.index(e), "b"), 1)]))
            rr = []
            for j in range(1, self.T - 1):
                a0, a1, a2 = dl(j - 1), dl(j), dl(j + 1)
                rr.append((a2[0] - 2 * a1[0] + a0[0], a2[1] - 2 * a1[1] + a0[1], -tr["M"], tr["M"]))
            brs = [rr]
        out = {}
        fbr = self.branches()
        for k, t in enumerate(posts):
            obj_x = np.zeros(nx); obj_x[ix(k, "theta")] = -1
            ivs = []
            for gb in brs:
                for rows in fbr:
                    iv = self.lp_range(Br(t), rows, nx, common + gb, xb, obj_x)
                    ivs.append(iv)
            out[t] = ivs
        return out

    def _rest_rate(self, t, r, c):
        v = self.zero(); den = self.N.sum() - self.A[:, r].sum()
        for i in range(self.p):
            for q in range(self.R):
                if q != r: v[t, i, q, c] = 1 / den
        return v

    def nrange(self, obj):
        ivs = [self.lp_range(obj, rows) for rows in self.branches()]
        ivs = [v for v in ivs if v is not None]
        return None if not ivs else [min(v[0] for v in ivs), max(v[1] for v in ivs)]


def union_hull(ivs):
    ivs = [v for v in ivs if v is not None]
    if not ivs: return None
    return [min(v[0] for v in ivs), max(v[1] for v in ivs)]

# ---------------------------------------------------------------------------
# instance generation
# ---------------------------------------------------------------------------
def rand_design(rng, p, R, C, T, Nlo, Nhi):
    while True:
        A, Bt = _rand_design(rng, p, R, C, T, Nlo, Nhi)
        x0 = [a[0] / sum(a) for a in A]
        if all(sum(a[r] for a in A) > 0 for r in range(R)) and max(x0) > min(x0):
            return A, Bt


def _rand_design(rng, p, R, C, T, Nlo, Nhi):
    A, Bt = [], [[] for _ in range(T)]
    for i in range(p):
        N = rng.randint(Nlo, Nhi)
        while True:
            a = [rng.randint(0, N) for _ in range(R - 1)]
            if sum(a) <= N: break
        a.append(N - sum(a))
        if R == 2 and (a[0] == 0 or a[0] == N) and rng.random() < 0.8:
            a = [max(1, min(N - 1, a[0])), 0]; a[1] = N - a[0]
        A.append(a)
        for t in range(T):
            b = [rng.randint(0, N) for _ in range(C - 1)]
            while sum(b) > N: b = [rng.randint(0, N) for _ in range(C - 1)]
            b.append(N - sum(b)); Bt[t].append(b)
    return A, Bt


def gen(seed):
    rng = random.Random(seed)
    inst = []
    def add(d): d["id"] = len(inst); inst.append(d)
    # ---- set A: 2x2, T = 2, brute force, single rows and pairs of rows -----------------
    rowlib = [
        [], [dict(type="sign_change", stype="dW", sign=1)], [dict(type="sign_change", stype="dW", sign=-1)],
        [dict(type="sign_change", stype="dB", sign=1)], [dict(type="sign_change", stype="dB", sign=-1)],
        [dict(type="cr_change")], [dict(type="cr_change", slack=0.05)],
        [dict(type="sign_change", stype="dW", version="W2", sign=1)],
        [dict(type="sign_change", stype="dB", version="pair", sign=-1)],
        [dict(type="order_change", sign=1)], [dict(type="order_change", sign=-1, slack=0.1)],
        [dict(type="contrast_change")], [dict(type="contrast_change", beliefs=["rate"])],
        [dict(type="cr_change"), dict(type="known")], [dict(type="sign_change", stype="dW", sign=1), dict(type="total")],
        [dict(type="cr_change", version="W2")],
    ]
    for k in range(80):
        A, Bt = rand_design(rng, 3, 2, 2, 2, 3, 6)
        rows = [dict(x) for x in rowlib[k % len(rowlib)]]
        d = dict(set="A", p=3, R=2, C=2, T=2, A=A, B=Bt, base=0, brute=True)
        tg = [dict(kind="did", row=0, versus="rest", col=0, **{"from": 0, "to": 1}),
              dict(kind="rate_at", row=0, col=0, election=1), dict(kind="gap_at", row=1, versus=0, col=0, election=0)]
        N = [sum(a) for a in A]; x = [A[i][0] / N[i] for i in range(3)]
        order = sorted(range(3), key=lambda i: x[i])
        L, H = [order[0]], order[1:]
        xbL = sum(A[i][0] for i in L) / sum(N[i] for i in L); xbH = sum(A[i][0] for i in H) / sum(N[i] for i in H)
        for rw in rows:
            rw.setdefault("row", 0); rw.setdefault("col", 0); rw["from"] = 0; rw["to"] = 1
            if rw["type"] == "known":
                rw.update(election=1, unit=0, value=int(min(A[0][0], Bt[1][0][0])))
            if rw["type"] == "total":
                lo = sum(max(0, A[i][0] + Bt[0][i][0] - N[i]) for i in range(3)); hi = sum(min(A[i][0], Bt[0][i][0]) for i in range(3))
                rw.update(election=0, units=[0, 1, 2], lower=float((lo + hi) // 2), upper=float((lo + hi) // 2 + 1))
            if rw["type"] == "contrast_change":
                rw.update(low=L, high=H)
        if any(r["type"] == "contrast_change" for r in rows):
            if xbH - xbL <= 1e-12 or any(sum(A[i][g] for i in b) == 0 for b in (L, H) for g in (0, 1)):
                rows = [r for r in rows if r["type"] != "contrast_change"]
            else:
                tg.append(dict(kind="bin_did", col=0, low=L, high=H, **{"from": 0, "to": 1}))
        d["rows"] = rows; d["targets"] = tg
        add(d)
    # ---- set B: 2x2, T = 3, brute force, rows on two pairs --------------------------
    for k in range(24):
        A, Bt = rand_design(rng, 2, 2, 2, 3, 3, 5)
        pairs = [(0, 1), (1, 2)] if k % 3 != 2 else [(0, 1), (0, 2)]
        kind = ["cr_change", "sign_change", "order_change"][k % 3]
        rows = []
        for (u, t) in pairs:
            rw = dict(type=kind, row=0, col=0, **{"from": u, "to": t})
            if kind == "sign_change": rw.update(stype="dW", sign=1 if k % 2 else -1)
            rows.append(rw)
        tg = [dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 0}),
              dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 2}),
              dict(kind="did", row=0, versus=1, col=0, **{"from": 0, "to": 2})]
        add(dict(set="B", p=2, R=2, C=2, T=3, A=A, B=Bt, base=1, brute=True, rows=rows, targets=tg))
    # ---- set C: R = 3 (C = 2) and C = 3 (R = 2), T = 2, brute force -------------------
    lib3 = [[], [dict(type="sign_change", stype="dB", sign=1)], [dict(type="sign_change", stype="dW", sign=-1)],
            [dict(type="sign_change", stype="dB", version="pair", other=2, sign=1)],
            [dict(type="sign_change", stype="dW", version="W2", sign=1)], [dict(type="cr_change")],
            [dict(type="cr_change", version="W2")], [dict(type="cr_change", version="pair", other=1)],
            [dict(type="order_change", over=1, sign=1)], [dict(type="order_change", sign=-1)]]
    for k in range(40):
        R, C = (3, 2) if k % 2 == 0 else (2, 3)
        A, Bt = rand_design(rng, 2, R, C, 2, 3, 4)
        rows = [dict(x) for x in lib3[k // 2 % len(lib3)]]
        if R == 2:
            for rw in rows:
                if rw.get("other") == 2: rw["other"] = 1
                if rw.get("over") is not None: rw["over"] = "rest"
        for rw in rows:
            rw.update(row=0, col=1 if C == 3 else 0, **{"from": 0, "to": 1})
        tg = [dict(kind="did", row=0, versus="rest", col=1 if C == 3 else 0, **{"from": 0, "to": 1}),
              dict(kind="did", row=0, versus=1, col=0, **{"from": 0, "to": 1})]
        add(dict(set="C", p=2, R=R, C=C, T=2, A=A, B=Bt, base=0, brute=True, rows=rows, targets=tg))
    # ---- set D: 3x3, T = 3, LP only ---------------------------------------------------
    libD = [[], [dict(type="sign_change", stype="dW", sign=1)], [dict(type="cr_change")],
            [dict(type="cr_change", version="W2", slack=0.02)],
            [dict(type="sign_change", stype="dB", version="pair", other=2, sign=1)],
            [dict(type="order_change", over=2, sign=1, slack=0.05)],
            [dict(type="cr_change", **{"from": 0, "to": 1}), dict(type="cr_change", **{"from": 1, "to": 2})]]
    for k in range(28):
        A, Bt = rand_design(rng, 5, 3, 3, 3, 20, 60)
        rows = [dict(x) for x in libD[k % len(libD)]]
        for rw in rows:
            rw.setdefault("from", 1); rw.setdefault("to", 2); rw.update(row=0, col=0)
        tg = [dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 0}),
              dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 2}),
              dict(kind="did", row=0, versus=1, col=2, **{"from": 0, "to": 2})]
        add(dict(set="D", p=5, R=3, C=3, T=3, A=A, B=Bt, base=1, brute=False, rows=rows, targets=tg))
    # ---- set E: causal (Thm E16), T = 3 or 4, LP ---------------------------------------
    trends = [dict(kind="none"), dict(kind="dial", M=0.0), dict(kind="dial", M=0.05), dict(kind="dial", M=0.3),
              dict(kind="rm", Mbar=1.0), dict(kind="rm", Mbar=0.3), dict(kind="sd", M=0.05),
              dict(kind="mtr", mtr=True, nd=True, ps=True), dict(kind="mtr", mtr=True, nd=True, ps=False),
              dict(kind="mtr", mtr=True, nd=False, ps=True)]
    for k in range(60):
        R = 3 if k % 4 == 3 else 2
        T = 4 if k % 5 == 4 else 3
        A, Bt = rand_design(rng, 4, R, 2, T, 30, 80)
        tr = dict(trends[k % len(trends)])
        b0 = 1 if T == 3 else (2 if k % 2 else 1)
        rows = []
        if k % 3 == 1: rows = [dict(type="cr_change", row=0, col=0, **{"from": b0, "to": b0 + 1})]
        if k % 3 == 2: rows = [dict(type="cr_change", row=0, col=0, **{"from": b0 - 1, "to": b0})]
        comp = ("rest" if k % 2 else 2) if R == 3 else 1
        Lam = [0.0, 0.0, 0.03][k % 3]
        rl = not (k % 7 == 6 and tr["kind"] in ("dial", "rm", "sd"))
        cs = dict(treated=0, comparison=comp, col=0, base=b0, trend=tr, Lambda=Lam, rate_limits=rl)
        if tr["kind"] == "mtr" and not tr["ps"]: cs["comparison"] = None
        add(dict(set="E", p=4, R=R, C=2, T=T, A=A, B=Bt, base=b0, brute=False, rows=rows,
                 targets=[], causal=cs))
    # causal with a coupling row (review issue 1 style) and known pre-period totals
    for k in range(12):
        A, Bt = rand_design(rng, 3, 2, 2, 3, 30, 60)
        N = [sum(a) for a in A]
        lo = sum(max(0, A[i][0] + Bt[0][i][0] - N[i]) for i in range(3)); hi = sum(min(A[i][0], Bt[0][i][0]) for i in range(3))
        mid = (lo + hi) / 2
        rows = [dict(type="total", row=0, col=0, election=0, units=[0, 1, 2], lower=float(mid - 2), upper=float(mid + 2))]
        if k % 2:
            y = min(A[0][0], Bt[1][0][0], Bt[2][0][0])
            rows.append(dict(type="linear", cells=[[0, 0, 0, 1, 3.0], [0, 0, 0, 2, -2.0]], lower=float(y) - 3, upper=float(y) + 3))
        tr = [dict(kind="dial", M=0.0), dict(kind="rm", Mbar=1.0), dict(kind="dial", M=0.05), dict(kind="sd", M=0.02)][k % 4]
        cs = dict(treated=0, comparison=1, col=0, base=1, trend=tr, Lambda=0.0, rate_limits=True)
        add(dict(set="F", p=3, R=2, C=2, T=3, A=A, B=Bt, base=1, brute=False, rows=rows, targets=[], causal=cs))
    return inst

# ---------------------------------------------------------------------------
# reference instances: Sec. 9.2, N21, review issues 1, 2, 12
# ---------------------------------------------------------------------------
def reference():
    out = []
    aW = [[10, 190], [50, 150], [150, 50], [190, 10]]
    yes = [[60, 62, 61, 56], [62, 66, 63, 61], [74, 86, 113, 120]]
    BW = [[[yes[t][i], 200 - yes[t][i]] for i in range(4)] for t in range(3)]
    crw = lambda a, b: dict(type="cr_change", row=0, col=0, **{"from": a, "to": b})
    def W(rows, tr, Lam=0.0, rl=True):
        return dict(set="W", p=4, R=2, C=2, T=3, A=aW, B=BW, base=1, brute=False, rows=rows, targets=[
            dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 0}),
            dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 2})],
            causal=dict(treated=0, comparison=1, col=0, base=1, trend=tr, Lambda=Lam, rate_limits=rl))
    for rows, tr, Lam in [([], dict(kind="none"), 0), ([], dict(kind="dial", M=0.0), 0),
                          ([crw(1, 2)], dict(kind="dial", M=0.0), 0), ([crw(1, 2)], dict(kind="dial", M=0.05), 0),
                          ([crw(1, 2)], dict(kind="dial", M=0.05), 0.02), ([], dict(kind="rm", Mbar=1.0), 0),
                          ([crw(1, 2)], dict(kind="rm", Mbar=1.0), 0), ([crw(0, 1)], dict(kind="rm", Mbar=1.0), 0),
                          ([crw(1, 2), crw(0, 1)], dict(kind="rm", Mbar=1.0), 0),
                          ([], dict(kind="mtr", mtr=True, nd=True, ps=True), 0),
                          ([crw(1, 2)], dict(kind="mtr", mtr=True, nd=True, ps=True), 0)]:
        d = W(rows, tr, Lam); d["label"] = "W " + json.dumps([r["from"] for r in rows]) + " " + tr["kind"]; out.append(d)
    a21 = [[20, 80], [60, 40]]; Y = [[30, 45], [32, 48], [30, 60]]
    B21 = [[[Y[t][i], 100 - Y[t][i]] for i in range(2)] for t in range(3)]
    for rl in (False, True):
        out.append(dict(set="N21", label="N21 rl=%s" % rl, p=2, R=2, C=2, T=3, A=a21, B=B21, base=1, brute=False, rows=[],
                        targets=[dict(kind="did", row=0, versus="rest", col=0, **{"from": 0, "to": 1}),
                                 dict(kind="did", row=0, versus="rest", col=0, **{"from": 1, "to": 2})],
                        causal=dict(treated=0, comparison=1, col=0, base=1, trend=dict(kind="rm", Mbar=1.0), Lambda=0.0, rate_limits=rl)))
    pins = []
    for t in range(2):
        lo = sum(max(0, a21[i][0] + Y[t][i] - 100) for i in range(2)); hi = sum(min(a21[i][0], Y[t][i]) for i in range(2))
        c = lo + (hi - lo) / 2
        pins.append(dict(type="total", row=0, col=0, election=t, units=[0, 1], lower=c - 2, upper=c + 2))
    for rl in (False, True):
        out.append(dict(set="N21", label="N21 narrow rl=%s" % rl, p=2, R=2, C=2, T=3, A=a21, B=B21, base=1, brute=False, rows=pins,
                        targets=[dict(kind="did", row=0, versus="rest", col=0, **{"from": 0, "to": 1})],
                        causal=dict(treated=0, comparison=1, col=0, base=1, trend=dict(kind="rm", Mbar=1.0), Lambda=0.0, rate_limits=rl)))
    # review issue 1: one unit, 100 / 100, 100 yes at both elections, 3 y0 - 2 y1 = 140
    out.append(dict(set="I1", label="issue 1", p=1, R=2, C=2, T=2, A=[[100, 100]], B=[[[100, 100]], [[100, 100]]], base=0,
                    brute=False, rows=[dict(type="linear", cells=[[0, 0, 0, 0, 3.0], [0, 0, 0, 1, -2.0]], lower=140.0, upper=140.0)],
                    targets=[dict(kind="did", row=0, versus="rest", col=0, **{"from": 0, "to": 1})],
                    causal=dict(treated=0, comparison=1, col=0, base=0, trend=dict(kind="dial", M=0.0), Lambda=0.0, rate_limits=True)))
    # review issue 2: four pure units
    aE = [[20, 0], [0, 80], [80, 0], [0, 20]]; yE = [[10, 40, 40, 4], [8, 40, 72, 16]]
    BE = [[[yE[t][i], sum(aE[i]) - yE[t][i]] for i in range(4)] for t in range(2)]
    out.append(dict(set="I2", label="issue 2", p=4, R=2, C=2, T=2, A=aE, B=BE, base=0, brute=False,
                    rows=[dict(type="contrast_change", col=0, low=[0, 1], high=[2, 3], **{"from": 0, "to": 1})],
                    targets=[dict(kind="bin_did", col=0, low=[0, 1], high=[2, 3], **{"from": 0, "to": 1})]))
    return out


def run_oracle(d):
    I = Inst(d)
    res = {"id": d["id"]}
    if d["targets"]:
        res["cont"] = I.cont(d["targets"])
        res["dd"] = [I.lp_range(I.target(tg), []) for tg in d["targets"]]
        if d.get("brute"):
            res["brute"] = I.brute(d["targets"])
    if d.get("causal"):
        cs = d["causal"]
        cz = I.causal(cs)
        res["causal"] = {str(t): v for t, v in cz.items()}
        r, s, c, b0 = cs["treated"], cs["comparison"], cs["col"], cs["base"]
        t1 = b0 + 1
        res["B_r"] = I.nrange(I.rate(t1, r, c)[0])
        cs0 = dict(cs, trend=dict(kind="none"), rate_limits=True)
        res["none"] = union_hull(I.causal(cs0)[t1])
        if s is not None:
            res["D"] = I.nrange(I.did(r, s, c, b0, t1)[0])
            if cs["trend"]["kind"] == "rm":
                supA = max(max(abs(v) for v in I.nrange(I.gap(e + 1, r, s, c)[0] - I.gap(e, r, s, c)[0])) for e in range(b0))
                M = cs["trend"]["Mbar"] * supA
                res["supA"] = supA
                res["dial_implied"] = union_hull(I.causal(dict(cs, trend=dict(kind="dial", M=M)))[t1])
    return res


def main():
    if sys.argv[1] == "gen":
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261126
        inst = gen(seed)
        for d in reference():
            d["id"] = len(inst); inst.append(d)
        for d in inst:
            d["oracle"] = run_oracle(d)
        json.dump(inst, open(sys.argv[2], "w"), default=float)
        print("instances:", len(inst))
        return
    inst = json.load(open(sys.argv[2])); res = json.load(open(sys.argv[3]))
    check(inst, res)


def close(a, b, tol=1e-7):
    if a is None or b is None: return a is None and b is None
    return all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def check(inst, res):
    byid = {r["id"]: r for r in res}
    fails = []; counts = {}
    def tally(key, ok, msg):
        counts.setdefault(key, [0, 0]); counts[key][0 if ok else 1] += 1
        if not ok: fails.append(msg)
    for d in inst:
        o = d["oracle"]; r = byid.get(d["id"])
        if r is None:
            tally(d["set"] + ":missing", False, "id %d: no R result" % d["id"]); continue
        if "error" in r:
            tally(d["set"] + ":error", False, "id %d (%s): R error %s" % (d["id"], d["set"], r["error"])); continue
        for j, tg in enumerate(d["targets"]):
            oc = union_hull(o["cont"][j]); rc = r["lp"][j]
            tally(d["set"] + ":cont", close(oc, rc), "id %d target %d cont oracle %s R %s rows %s" % (d["id"], j, oc, rc, d["rows"]))
            tally(d["set"] + ":dd", close(o["dd"][j], r["dd"][j]), "id %d target %d dd oracle %s R %s" % (d["id"], j, o["dd"][j], r["dd"][j]))
            if d.get("brute") and o.get("brute") is not None:
                ob = o["brute"][j]; ri = r["int"][j]
                tally(d["set"] + ":int", close(ob, ri), "id %d target %d int oracle %s R %s (status %s) rows %s" % (d["id"], j, ob, ri, r["status"][j], d["rows"]))
        if d.get("causal"):
            for t, ivs in o["causal"].items():
                oc = union_hull(ivs)
                rc = r["causal"].get(t)
                tally(d["set"] + ":causal", close(oc, rc, 1e-6),
                      "id %d post %s causal oracle %s R %s trend %s rows %s rl %s" % (d["id"], t, oc, rc, d["causal"]["trend"], d["rows"], d["causal"].get("rate_limits")))
            t1 = str(d["causal"]["base"] + 1)
            if "none" in r: tally("E15", close(o["none"], r["none"], 1e-6), "id %d none oracle %s R %s" % (d["id"], o["none"], r["none"]))
            if "B_r" in r and r["B_r"] is not None and o["none"] is not None and d["causal"].get("comparison") is not None:
                tally("E15 closed form", close([o["B_r"][0] - 1, o["B_r"][1]], o["none"], 1e-6), "id %d E15" % d["id"])
            tr = d["causal"]["trend"]; Lam = d["causal"]["Lambda"]
            oc = union_hull(o["causal"][t1])
            if tr["kind"] == "dial" and "dial" in r and oc is not None and d["causal"].get("rate_limits", True):
                ML = tr["M"] + Lam
                eq_hi = abs(oc[1] - (o["D"][1] + ML)) < 1e-7; eq_lo = abs(oc[0] - (o["D"][0] - ML)) < 1e-7
                tally("E17(b) criterion", eq_hi == r["dial"]["crit_upper"] and eq_lo == r["dial"]["crit_lower"],
                      "id %d E17(b): oracle eq (%s, %s) R criterion (%s, %s)" % (d["id"], eq_lo, eq_hi, r["dial"]["crit_lower"], r["dial"]["crit_upper"]))
                tally("E17(a) outer", oc[0] >= o["D"][0] - ML - 1e-7 and oc[1] <= o["D"][1] + ML + 1e-7, "id %d E17(a)" % d["id"])
            if tr["kind"] == "rm" and "rm" in r and oc is not None:
                tally("E19 supA", abs(o["supA"] - r["rm"]["supA"]) < 1e-7, "id %d supA %s %s" % (d["id"], o["supA"], r["rm"]["supA"]))
                tally("E19(a) dial set", close(o["dial_implied"], r["rm"]["dial_set"], 1e-6), "id %d dial set %s %s" % (d["id"], o["dial_implied"], r["rm"]["dial_set"]))
                tally("E19(a) inclusion", oc[0] >= o["dial_implied"][0] - 1e-7 and oc[1] <= o["dial_implied"][1] + 1e-7, "id %d E19(a)" % d["id"])
                if r["rm"].get("vac_upper") is not None and d["causal"].get("rate_limits", True):
                    vu = abs(oc[1] - o["none"][1]) < 1e-7; vl = abs(oc[0] - o["none"][0]) < 1e-7
                    tally("E19(b) vacuity", vu == r["rm"]["vac_upper"] and vl == r["rm"]["vac_lower"],
                          "id %d E19(b): oracle (%s, %s) R (%s, %s) set %s none %s" % (d["id"], vl, vu, r["rm"]["vac_lower"], r["rm"]["vac_upper"], oc, o["none"]))
    print("check                         pass  fail")
    for k in sorted(counts): print("%-28s %5d %5d" % (k, counts[k][0], counts[k][1]))
    print("\nTOTAL: %d pass, %d fail" % (sum(v[0] for v in counts.values()), sum(v[1] for v in counts.values())))
    for f in fails[:60]: print("FAIL:", f)


if __name__ == "__main__":
    main()
