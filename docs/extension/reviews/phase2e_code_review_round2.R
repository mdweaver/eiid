# Round-2 confirmation (code review of PHASE2e): package values for the E30/E35 checks (compare with
# phase2e_code_review_round2.py), E16(b) on finding 5, E17(c) minimality, E29/E35+LCE printouts and a scan
# of every printed interval for an IDENTIFICATION / DESCRIPTIVE / CONFIDENCE label.
# Usage: Rscript phase2e_code_review_round2.R [lib]
args <- commandArgs(TRUE); if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(eiPartialID))
p8 <- function(lab, v) cat(sprintf("%-40s [%.6f, %.6f]\n", lab, v[1], v[2]))
iv <- function(a, k = 1) c(a$results[[k]]$lower, a$results[[k]]$upper)
outs <- list()
keep <- function(nm, a) { outs[[nm]] <<- capture.output(print(a)); a }

# ---- E30(c) three precincts (expected WPT) ----
mk <- function(vet, size, y) data.frame(u = paste0("P", seq_along(vet)), v = vet, o = size - vet, yes = y, no = size - y)
f3 <- eiChange(list(b = mk(c(200, 400, 600), 800, c(4, 164, 324)), p = mk(c(200, 400, 600), 800, c(204, 124, 4))),
               cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
W <- function(l, ...) eiATT(f3, treated = "v", comparison = "o", population = "fixed", trend = within_pt(target = "expected", ...),
                           spillover = 0, no_anticipation = TRUE, lce = l)
cat("E30(c)\n")
p8("no LCE", iv(W(lce_mean(change = FALSE)))); p8("change only", iv(W(lce_mean(change = TRUE))))
p8("base+change", iv(a30 <- keep("E30", W(lce_mean(base = TRUE))))); p8("base+post", iv(W(lce_mean(base = TRUE, post = TRUE, change = FALSE))))
p8("post+change", iv(W(lce_mean(post = TRUE)))); p8("eps = 1/2", iv(W(lce_mean(base = TRUE, epsilon = .5))))
cat(sprintf("epsilon* = %.8f\n", a30$results[[1]]$breakdown_eps))
a30p <- keep("E29_persons", eiATT(f3, treated = "v", comparison = "o", population = "fixed", trend = within_pt(target = "expected"),
                                 spillover = 0, no_anticipation = TRUE, lce = lce_mean(), independence = "persons"))

# ---- E28 five-precinct table with the known law (points) ----
cat("E28 (points)\n")
for (ex in list(rare = list(20000, c(5000, 8000, 11000, 14000, 17000), c(40, 12000, 12000, 12000, 12000),
                            rbind(c(50, 30), c(9964, 9956), c(8188, 8192), c(6682, 6698), c(5446, 5474))),
                ordinary = list(1000, c(100, 300, 500, 700, 900), rep(500, 5), rbind(c(682, 650), c(606, 590), c(530, 530), c(454, 470), c(378, 410))))) {
  f <- eiChange(list(b = mk(ex[[2]], ex[[1]], ex[[3]]), p = mk(ex[[2]], ex[[1]], rowMeans(ex[[4]]))), cbind(yes, no) ~ cbind(v, o),
                unit = "u", base = "b", row = "v", col = "yes", integer = "lp")
  g <- function(tr) 100 * iv(eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = tr, spillover = 0, no_anticipation = TRUE,
                                   lce = lce_mean(change = TRUE), law = list(p = ex[[4]])))
  cat(sprintf("  WPT [%.4f, %.4f]  PPT [%.4f, %.4f]\n", g(within_pt(target = "expected"))[1], g(within_pt(target = "expected"))[2], g(pooled_pt())[1], g(pooled_pt())[2]))
}

# ---- E35 worked example ----
aW <- cbind(c(10, 50, 150, 190), c(190, 150, 50, 10)); yes <- rbind(c(60, 62, 61, 56), c(62, 66, 63, 61), c(74, 86, 113, 120))
mk4 <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes[t, ], no = 200 - yes[t, ])
kc <- function(u, e, v) list(unit = u, row = "oth", col = "yes", value = v, election = e)
aud <- list(known = list(kc("P3", "base", 20), kc("P3", "post", 22), kc("P4", "base", 4), kc("P4", "post", 5)))
fit4 <- function(rows = NULL, ...) eiChange(list(early = mk4(1), base = mk4(2), post = mk4(3)), cbind(yes, no) ~ cbind(vet, oth),
                                           unit = "p", base = "base", row = "vet", col = "yes", integer = "lp", restrictions = rows, ...)
E <- function(f, ..., lce = NULL, sp = 0) eiATT(f, treated = "vet", comparison = "oth", population = "fixed", trend = within_pt(...),
                                    spillover = sp, no_anticipation = TRUE, lce = lce)
fa <- fit4(aud)
cat("E35\n")
p8("margins only, M = 0", iv(keep("E35_margins", E(fit4()))))
a0 <- keep("E35_M0", E(fa, M = 0)); p8("audited, M = 0", iv(a0)); p8("audited, M = 0.1", iv(E(fa, M = .1)))
cat(sprintf("M* = %.6f\n", a0$results[[1]]$breakdown_M)); p8("D_W count set", a0$results[[1]]$dw)
fm <- tryCatch(fit4(aud, model = fd()), error = function(e) { cat("fd() fit:", conditionMessage(e), "\n"); NULL })
if (!is.null(fm)) { al <- tryCatch(keep("E35_LCE", E(fm, M = 0, lce = TRUE)), error = function(e) conditionMessage(e)); if (is.character(al)) cat("lce:", al, "\n") else p8("audited, M = 0, lce = TRUE", iv(al)) }
keep("E35_spill", E(fa, M = 0, sp = .05)); keep("E35_amin", E(fa, M = 0, a_min = 20))
keep("PPT", eiATT(fa, treated = "vet", comparison = "oth", population = "fixed", trend = pooled_pt(), spillover = 0, no_anticipation = TRUE))

# ---- E16(b) whole-person witness flag on finding 5 (m = 3): integer ATT set {0}, continuous [-1, 1] ----
m <- 3
h <- function(y) data.frame(u = "U", tr = m, co = m + 1, yes = y, no = 2 * m + 1 - y)
f5 <- eiChange(list(b = h(m), p = h(m)), cbind(yes, no) ~ cbind(tr, co), unit = "u", base = "b", row = "tr", col = "yes", integer = "lp")
a5 <- keep("E16b", eiATT(f5, treated = "tr", comparison = "co", population = "fixed", trend = trend_dial(0), spillover = 0, no_anticipation = TRUE))
p8("finding 5 (m = 3), dial M = 0", iv(a5)); cat(grep("E16\\(b\\)", outs$E16b, value = TRUE), sep = "\n")

# ---- E17(c) sharp breakdown: minimal ----
yes2 <- rbind(c(62, 66, 63, 61), c(74, 86, 113, 120))
mk2 <- function(t) data.frame(p = paste0("P", 1:4), vet = aW[, 1], oth = aW[, 2], yes = yes2[t, ], no = 200 - yes2[t, ])
fc <- eiChange(list(base = mk2(1), post = mk2(2)), cbind(yes, no) ~ cbind(vet, oth), unit = "p", base = "base", row = "vet",
               col = "yes", integer = "lp", restrictions = list(cr_change = list(row = "vet", col = "yes", from = "base", to = "post")))
D <- function(M) eiATT(fc, treated = "vet", comparison = "oth", population = "fixed", trend = trend_dial(M), spillover = 0, no_anticipation = TRUE)
d0 <- keep("E17", D(0)); sb <- d0$results[[1]]$dial$sharp_breakdown
cat(sprintf("E17(c): outer %.6f sharp %.6f; lower at 0.999 sharp %.3e, at sharp %.3e\n", d0$results[[1]]$dial$breakdown, sb,
            iv(D(0.999 * sb))[1], iv(D(sb))[1]))
# a few other trends for the label scan
for (tr in list(rel_mag(.5), smoothness(.05), mtr_nd_ps(), no_trend()))
  keep(paste("trend", class(tr)[1], length(outs)), eiATT(fit4(aud), treated = "vet", comparison = "oth", population = "fixed", trend = tr,
                                                       spillover = 0, no_anticipation = TRUE))

# ---- label scan ----
cat("\nLabel scan: printed lines with an interval [a, b] and no label tag on the line\n")
for (nm in names(outs)) {
  L <- outs[[nm]]; hit <- grepl("\\[-?[0-9.]+, -?[0-9.]+\\]|EMPTY", L) & !grepl("\\[(IDENTIFICATION|DESCRIPTIVE|CONFIDENCE)", L)
  if (any(hit)) cat(sprintf("  %s: %s\n", nm, trimws(L[hit])), sep = "")
}
cat("\n---- E30 printout ----\n"); cat(outs$E30, sep = "\n")
cat("\n---- E29 persons printout (support / independence lines) ----\n"); cat(grep("support|independence|CONFIDENCE|coverage|IDENTIFICATION", outs$E29_persons, value = TRUE), sep = "\n")
cat("\n---- E35 M=0 printout ----\n"); cat(outs$E35_M0, sep = "\n")
if (!is.null(outs$E35_LCE)) { cat("\n---- E35 + LCE printout ----\n"); cat(outs$E35_LCE, sep = "\n") }
cat("\n---- E35 a_min = 20 (P4 unsupported) ----\n"); cat(grep("support|unsupported|ATT", outs$E35_amin, value = TRUE), sep = "\n")

# ---- E29 confidence line under a local dial M > 0 or spillover Lambda > 0 (round-2 issue R2-1) ----
cat("\nE29 with M > 0 / Lambda > 0 (D_W pinned by known comparison cells)\n")
local({
  set.seed(2); p <- 60; vet <- sample(300:700, p, TRUE); sz <- 1000
  nv0 <- rbinom(p, vet, .4); nv1 <- rbinom(p, vet, .5); no0 <- rbinom(p, sz - vet, .4); no1 <- rbinom(p, sz - vet, .42)
  mk <- function(y) data.frame(u = paste0("P", 1:p), v = vet, o = sz - vet, yes = y, no = sz - y)
  kn <- c(lapply(1:p, function(i) list(unit = paste0("P", i), row = "o", col = "yes", value = no0[i], election = "b")),
          lapply(1:p, function(i) list(unit = paste0("P", i), row = "o", col = "yes", value = no1[i], election = "p")))
  f <- eiChange(list(b = mk(nv0 + no0), p = mk(nv1 + no1)), cbind(yes, no) ~ cbind(v, o), unit = "u", base = "b", row = "v", col = "yes",
                integer = "lp", restrictions = list(known = kn))
  for (M in c(0, 0.2)) for (L in c(0, 0.1)) {
    e <- eiATT(f, treated = "v", comparison = "o", population = "fixed", trend = within_pt(target = "expected", M = M),
               spillover = L, no_anticipation = TRUE, lce = lce_mean(), independence = "persons")$results[[1]]
    cat(sprintf("M=%.2f Lambda=%.2f  ID bound [%.4f, %.4f]  D_W [%.4f, %.4f]  E29 CI [%.4f, %.4f]\n", M, L, e$lower, e$upper, e$dw$LU[1], e$dw$LU[2], e$dw$ci[1], e$dw$ci[2]))
  }
})
