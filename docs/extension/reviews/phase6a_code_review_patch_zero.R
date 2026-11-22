# Review-only workaround for round-2 issue R2-1 (used with PATCH=... to test everything else); not a proposed fix.
# review-only patch (issue R2-1): a component with a zero objective needs no certificate
ns <- asNamespace("eiPartialID")
f0 <- ns$.ec_glpk
g <- function(nv, obj, trip, dir, rhs, sense, time_limit) {
  if (all(obj == 0) && max(trip[, 1]) < length(dir)) {          # certificate row 0 >= 1: infeasible by construction
    return(list(status = 4, x = NULL, elapsed = 0))
  }
  f0(nv, obj, trip, dir, rhs, sense, time_limit)
}
unlockBinding(".ec_glpk", ns); assign(".ec_glpk", g, envir = ns)
