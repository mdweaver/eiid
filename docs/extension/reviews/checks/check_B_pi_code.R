# Checks of referee report B against the package code (no package install needed).
# Run from the repository root:  Rscript docs/extension/reviews/checks/check_B_pi_code.R
# Needs only lpSolve. Sources the internal rxc_*.R files and calls internal helpers.
suppressMessages(library(lpSolve))
for (f in c("rxc_core.R", "rxc_infer.R", "rxc_sign.R", "rxc_modules.R"))
  source(file.path("package/eiPartialID/R", f))
ok <- function(cond, msg) cat(if (isTRUE(cond)) "PASS" else "FAIL", "-", msg, "\n")

# Unit square [0,1]^2 in z, as rows s + A z >= 0
sq <- list(s = c(0, 1, 0, 1), A = rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1)))

# --- R2 step 4 (5a l.169): target with intercept 2, phi = 2 + z1, U = 3 -----------------
tg <- list(cvec = numeric(0), c0 = 2, gamma = c(1, 0))
keep <- list(list(s = 1.2, A = rbind(c(-1, -1))))   # z1 + z2 <= 1.2 meets the face z1 = 1
cut  <- list(list(s = 0.8, A = rbind(c(-1, -1))))   # z1 + z2 <= 0.8 removes the face
ok(!.rxc_face_moves(sq, tg, numeric(0), keep, "upper", 3, 1e-7),
   "face test with target intercept: cut through the facet does not move U (code includes the intercept)")
ok(.rxc_face_moves(sq, tg, numeric(0), cut, "upper", 3, 1e-7),
   "face test with target intercept: cut removing the facet moves U")
# The printed doc row h'k >= U (without subtracting phi_0) would be z1 >= 3: infeasible on F
ok(.rxc_meet(c(sq$s, -3), rbind(sq$A, c(1, 0)), numeric(0), matrix(0, 0, 2)) == -Inf,
   "doc version of the face LP (h'k >= U, intercept omitted) is infeasible: the doc error is real")

# --- R2(e): 'added last' gain zero for z1 but the row is not implied ----------------------
Fo <- list(s = c(sq$s, 0.5), A = rbind(sq$A, c(-1, 0)))    # square and z1 <= .5
N  <- list(list(s = 0.8, A = rbind(c(-1, -1))))            # z1 + z2 <= .8
ok(!.rxc_implied(Fo, N, 1e-7), "row z1+z2<=.8 is not implied by the others (set-level)")
tx <- list(cvec = numeric(0), c0 = 0, gamma = c(1, 0))
ok(!.rxc_face_moves(Fo, tx, numeric(0), N, "upper", 0.5, 1e-7),
   "... and does not move the upper end of z1: code verdict is 'cuts the set, not this target'")
txy <- list(cvec = numeric(0), c0 = 0, gamma = c(1, 1))
ok(.rxc_face_moves(Fo, txy, numeric(0), N, "upper", 1.5, 1e-7),
   "... but moves the upper end of z1+z2 (target-specific zero gain, as the referee says)")

# --- verify_R10 union length: the package merges overlaps before measuring ---------------
u <- .rxc_union(list(c(0, 1), c(0.5, 1.5)))
ok(nrow(u) == 1 && abs(.rxc_union_len(u) - 1.5) < 1e-12,
   "package union length merges overlapping branch intervals (only verify_R10.py double-counts)")

# --- Proof-status labels of the M-SIGN modules ---------------------------------------------
ps <- c(assoc = .rxc_registry$assoc$meta$proof_status, cr = .rxc_registry$cr$meta$proof_status)
print(ps)
ok(all(grepl("^proved", ps)),
   "assoc and cr are labelled 'proved', so print() shows no EXPERIMENTAL flag (NEXT_BUILD says they should be)")
