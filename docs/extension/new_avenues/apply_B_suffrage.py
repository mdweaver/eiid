"""
Quick application of Avenue B (context at a second scale) to the veterans / Black-suffrage data.
Usage: python3 apply_B_suffrage.py <wi_twp_analysis_cleaned.csv> <ia_twp_analysis_cleaned.csv>
(The data are the owner's and are not committed.)

Per state and year, a 2x3 problem: rows veterans (v; future veterans before the war) and
non-veterans (n); columns yes / no / other (abstain, roll-off, non-voting) on a base of
eligible men.  Unit = township (or the state's stable unit); z = the unit's county (WI: county
cluster) veteran share, base-weighted over the county's units.

Base (simple choice for this quick look, not the DiD session's interval-valued base):
  WI: N = max(elig_YYYY [1860 census], yes + no, vets);   IA: N = max(eev, yes + no, vets).
Units with negative or missing yes/no/base are dropped.

Model per column c:  rate_v = a_v + b_v x + g_v z,  rate_n = a_n + b_n x + g_n z.
OLS (unweighted) of T_c on [1, x, x^2, z, x z] identifies a_n, g_n, g_v, d = b_v - b_n and
a_v - a_n + b_n; k_c := b_n is unknown, with sum_c k_c = 0.
Rungs: DD; M0 (pooled LCE, no z); B0 (with z, no restriction); B1(M): |b_r,c| <= M |g_r,c|.
Constraints: all three rates in [0,1] at every observed (x, z) (a valid, slightly conservative
version of the hull-vertex outer set).  Targets: base-weighted model-implied yes rates of v and n,
and their gap.  Plug-in only: no confidence intervals.
"""
import sys
import numpy as np, pandas as pd
from scipy.optimize import linprog

def prep_wi(df, y):
    yes, no = df[f"suff_yes_{y}"], df[f"suff_no_{y}"]
    ok = df[f"elig_{y}"].notna() & df.vets.notna() & (yes >= 0) & (no >= 0)
    d = df[ok]; yes, no = yes[ok], no[ok]
    N = np.maximum.reduce([d[f"elig_{y}"].values, (yes + no).values, d.vets.values])
    return pd.DataFrame(dict(county=d.County_Cluster.values, N=N, V=d.vets.values,
                             yes=yes.values, no=no.values))

def prep_ia(df, y):
    yes, no, eev = df[f"suff_yes_{y}"], df[f"suff_no_{y}"], df[f"suff_eev_{y}"]
    ok = eev.notna() & df.vets.notna() & (yes >= 0) & (no >= 0)
    d = df[ok]
    N = np.maximum.reduce([eev[ok].values, (yes[ok] + no[ok]).values, d.vets.values])
    return pd.DataFrame(dict(county=d.County.values, N=N, V=d.vets.values,
                             yes=yes[ok].values, no=no[ok].values))

def analyse(d, label):
    d = d[d.N > 0].copy()
    d["x"] = d.V / d.N
    g = d.groupby("county")[["V", "N"]].transform("sum"); d["z"] = g.V / g.N
    d["other"] = d.N - d.yes - d.no
    x, z, N = d.x.values, d.z.values, d.N.values
    T = np.column_stack([d.yes / d.N, d.no / d.N, d.other / d.N])
    out = [f"### {label}: {len(d)} units, {d.county.nunique()} counties; "
           f"veteran share mean {np.average(x, weights=N):.3f}",
           f"separation of own and county share: SD(x) {x.std():.3f}, SD(x - z) {(x - z).std():.3f}, "
           f"corr(x, z) {np.corrcoef(x, z)[0, 1]:.2f}; units alone in their county: "
           f"{(d.groupby('county').x.transform('size') == 1).mean():.2f}"]
    # Duncan-Davis for yes rates
    Y = d.yes.values; V = d.V.values; NV = N - V
    dd_v = (np.sum(np.maximum(0, V - (N - Y))) / V.sum(), np.sum(np.minimum(V, Y)) / V.sum())
    dd_n = (np.sum(np.maximum(0, NV - (N - Y))) / NV.sum(), np.sum(np.minimum(NV, Y)) / NV.sum())
    def fit(with_z):
        cols = [np.ones_like(x), x, x * x] + ([z, x * z] if with_z else [])
        Xd = np.column_stack(cols)
        C = np.linalg.lstsq(Xd, T, rcond=None)[0]            # rows: terms, cols: yes/no/other
        P = []
        for c in range(3):
            a_n, c1, dd = C[0, c], C[1, c], C[2, c]
            g_n = C[3, c] if with_z else 0.0
            g_v = (C[3, c] + C[4, c]) if with_z else 0.0
            P.append(dict(a_n=a_n, c1=c1, d=dd, g_n=g_n, g_v=g_v))
        return P
    def rates(p, k):
        rn = p["a_n"] + k * x + p["g_n"] * z
        rv = p["c1"] + p["a_n"] + p["d"] * x + p["g_v"] * z + k * (x - 1)
        return rv, rn
    def bounds(P, M=None, mcols=(0, 1, 2), pts=None):
        pts = np.arange(len(x)) if pts is None else pts
        # variables k_yes, k_no; k_other = -k_yes - k_no
        A, b = [], []
        coefs = [np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([-1.0, -1.0])]
        for c in range(3):
            p, e = P[c], coefs[c]
            rv0, rn0 = rates(p, 0.0)
            # rv = rv0 + (x-1) k_c >= 0 ; rn = rn0 + x k_c >= 0  (<= 1 implied by the other columns)
            for i in pts:
                A.append(-(x[i] - 1) * e); b.append(rv0[i])
                A.append(-x[i] * e); b.append(rn0[i])
            if M is not None and c in mcols:   # |k_c| <= M|g_n| ; |d + k_c| <= M|g_v|
                A += [e, -e, e, -e]
                b += [M * abs(p["g_n"]), M * abs(p["g_n"]), M * abs(p["g_v"]) - p["d"], M * abs(p["g_v"]) + p["d"]]
        A, b = np.array(A), np.array(b)
        p = P[0]; rv0, rn0 = rates(p, 0.0)
        wv, wn = N * x, N * (1 - x)
        tv = (np.sum(wv * rv0) / wv.sum(), np.sum(wv * (x - 1)) / wv.sum())   # value = const + slope*k_yes
        tn = (np.sum(wn * rn0) / wn.sum(), np.sum(wn * x) / wn.sum())
        res = {}
        for name, (c0, s) in (("vet", tv), ("non", tn), ("gap", (tv[0] - tn[0], tv[1] - tn[1]))):
            lo = linprog([s, 0], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method="highs")
            hi = linprog([-s, 0], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method="highs")
            if lo.status != 0 or hi.status != 0:
                return None
            res[name] = (c0 + lo.fun, c0 - hi.fun)
        return res
    P0, PB = fit(False), fit(True)
    out.append(f"identified county-share effects on the YES rate: veterans g_v = {PB[0]['g_v']:+.3f}, "
               f"non-veterans g_n = {PB[0]['g_n']:+.3f}; gap in own-share slopes d = {PB[0]['d']:+.3f}")
    keep = np.where((N >= 100) & (x >= np.quantile(x, .025)) & (x <= np.quantile(x, .975))
                    & (z >= np.quantile(z, .025)) & (z <= np.quantile(z, .975)))[0]
    fmt = lambda t: f"[{t[0]:.3f}, {t[1]:.3f}]"
    Ms = np.round(np.arange(0.0, 8.01, 0.05), 2)
    for tag, pts in (("all units as constraint points", None),
                     (f"trimmed constraint points (N >= 100, x and z inside 2.5-97.5%; {len(keep)} units)", keep)):
        out.append(f"**{tag}**")
        out.append("| rung | veterans' yes rate | non-veterans' yes rate | gap (v - n) |")
        out.append("|---|---|---|---|")
        out.append(f"| Duncan-Davis | {fmt(dd_v)} | {fmt(dd_n)} | [{dd_v[0]-dd_n[1]:.3f}, {dd_v[1]-dd_n[0]:.3f}] |")
        rows = [("M0 (pooled LCE)", P0, None, (0, 1, 2)), ("B0 (+ county share)", PB, None, (0, 1, 2))] + \
               [(f"B1(M={M}), yes column only", PB, M, (0,)) for M in (2, 1, 0.5)] + \
               [(f"B1(M={M}), all columns", PB, M, (0, 1, 2)) for M in (2, 1)]
        for lab, P, M, mc in rows:
            r = bounds(P, M, mc, pts)
            out.append(f"| {lab} | " + (" | ".join(fmt(r[k]) for k in ("vet", "non", "gap")) if r else "empty | empty | empty") + " |")
        for mc, nm in (((0,), "yes column only"), ((0, 1, 2), "all columns")):
            mmin = next((M for M in Ms if bounds(PB, M, mc, pts) is not None), None)
            out.append(f"M_min ({nm}, 0.05 grid up to 8): {mmin}")
        out.append("")
    store_pts = keep
    return out, (bounds(PB, None, pts=keep), {M: bounds(PB, M, (0,), keep) for M in (2, 1, 0.5)}, bounds(P0, None, pts=keep))

if __name__ == "__main__":
    wi = pd.read_csv(sys.argv[1]); ia = pd.read_csv(sys.argv[2])
    lines = ["# Avenue B applied to the veterans / suffrage data (quick look, plug-in, no CIs)", ""]
    store = {}
    for lab, d in (("Wisconsin 1857", prep_wi(wi, "1857")), ("Wisconsin 1865", prep_wi(wi, "1865")),
                   ("Iowa 1857", prep_ia(ia, "57")), ("Iowa 1868", prep_ia(ia, "68"))):
        o, s = analyse(d, lab); lines += o + [""]; store[lab] = s
    # DiD of the gap with no cross-period link: difference of the per-period intervals
    lines.append("### DiD of the veteran - non-veteran yes gap, no cross-period link (interval difference; trimmed points; B1 on the yes column)")
    lines.append("| state | M0 | B0 | B1(M=2) | B1(M=1) | B1(M=0.5) |")
    lines.append("|---|---|---|---|---|---|")
    for st, a, b in (("Wisconsin", "Wisconsin 1857", "Wisconsin 1865"), ("Iowa", "Iowa 1857", "Iowa 1868")):
        def did(ra, rb):
            if ra is None or rb is None: return "empty"
            return f"[{rb['gap'][0]-ra['gap'][1]:.3f}, {rb['gap'][1]-ra['gap'][0]:.3f}]"
        A, B = store[a], store[b]
        lines.append(f"| {st} | {did(A[2], B[2])} | {did(A[0], B[0])} | {did(A[1][2], B[1][2])} | "
                     f"{did(A[1][1], B[1][1])} | {did(A[1][0.5], B[1][0.5])} |")
    txt = "\n".join(lines); print(txt)
    open(__file__.replace(".py", "_output.md"), "w").write(txt + "\n")
