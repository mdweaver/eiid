# Cross-check of phase3_montecarlo.py against the R package (eiBounds) on one replication of
# design A (10 x 93, M1) and of design S2 (20 x 50, common non-affine shock; M1 and M0).
# Writes the simulated data to CSV, prints the Python intervals, then runs the same data through
# eiPartialID::eiBounds via Rscript. Run: python3 -I phase3_montecarlo_crosscheck.py
import os, sys, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import phase3_montecarlo as P

tmp = tempfile.mkdtemp()
for nm, simargs, seed in (("A", (10, 93), 1000), ("S2", (20, 50, "nonaffine"), 91000)):
    s = P.sim_A(np.random.default_rng(seed), *simargs)
    np.savetxt(os.path.join(tmp, f"{nm}.csv"), np.column_stack([s["x"], s["t"], s["N"], s["d"]]), delimiter=",")
    r = P.run22(seed, P.sim_A, simargs, ("M1", "M0"), False)
    for mo in ("M1", "M0"):
        rr = r[mo]
        f = lambda iv: "empty" if iv is None or iv[0] > iv[1] else f"[{iv[0]:.4f}, {iv[1]:.4f}]"
        print(f"python {nm} {mo}: plug-in {f(rr['plug'])} | CLR B~ {f(rr['clr_Bt'])} | proj B~ {f(rr['proj_Bt'])} "
              f"| CLR B {f(rr['clr_B_bd'])} | proj B {f(rr['proj_B'])}"
              + (f" | pool W {rr['pool_W']:.3f} df {rr['pool_df']}" if "pool_W" in rr else ""))
R = r'''
suppressMessages(library(eiPartialID))
for (nm in c("A", "S2")) {
  d <- read.csv(file.path("%s", paste0(nm, ".csv")), header = FALSE)
  for (mo in c("M1", "M0")) {
    f <- eiBounds(x = d$V1, t = d$V2, n = d$V3, group = if (mo == "M1") d$V4 else NULL, model = mo,
                  targets = list(b = c(1, 1)), seed = 1, ndraw = 200)
    tb <- f$table; g <- function(a, b) if (is.na(a)) "empty" else sprintf("[%%.4f, %%.4f]", a, b)
    cat(sprintf("R      %%s %%s: plug-in %%s | CLR B~ %%s | proj B~ %%s | CLR B %%s | proj B %%s%%s\n", nm, mo,
        g(tb$estimate_lower, tb$estimate_upper), g(tb$ci_model_lower, tb$ci_model_upper),
        g(tb$proj_model_lower, tb$proj_model_upper), g(tb$ci_actual_lower, tb$ci_actual_upper),
        g(tb$proj_actual_lower, tb$proj_actual_upper),
        if (is.null(f$tests$pool)) "" else sprintf(" | pool W %%.3f df %%d", f$tests$pool$statistic, f$tests$pool$df)))
  }
}
''' % tmp
print(subprocess.run(["Rscript", "-e", R], capture_output=True, text=True).stdout)
print("Notes: the package intersects the model-implied CI with [0, 1] and ci_actual/proj_actual with Duncan-Davis;\n"
      "CLR ends differ by up to about 0.01: the Python script uses every constraint as a candidate (the full dual-vertex\n"
      "set in 2x2), the package collects candidates by simulation, and both simulate the critical values.")
