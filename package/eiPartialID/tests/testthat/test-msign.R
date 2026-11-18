# M-SIGN (build step 1): signs on mean associations, contextual reinforcement (CR),
# redundancy diagnostics. Theory: docs/extension/PHASE5a_sign_module.md (R2, R3, R10).
# Reference numbers: docs/extension/make_msign_reference.py, which uses verify_R3.py's
# assoc_3x2 / assoc_alt and writes msign_example_3x2.csv, msign_gap_3x2.csv and
# msign_reference_3x2.csv (affine forms "c, a1..a3" in k; for interval rows c = lower,
# a1 = upper).

ref <- utils::read.csv(test_path("msign_reference_3x2.csv"))
ex <- utils::read.csv(test_path("msign_example_3x2.csv"))
gp <- utils::read.csv(test_path("msign_gap_3x2.csv"))
tg3 <- list(B1 = c(1, 1), B2 = c(2, 1), B3 = c(3, 1))
Dgap <- list(D = matrix(c(1, -1, 0, 0), 2, 2))

fit3 <- function(d, ...) eiBounds(x = as.matrix(d[, c("X1", "X2", "X3")]), t = as.matrix(d[, c("T1", "T2")]),
                                  n = d$N, model = "M0", targets = tg3, tests = FALSE, ...)
gap_order <- ref[ref$case == "gap" & ref$name == "order", ]
gap_restr <- function(cr = NULL, assoc = NULL)
  c(list(order = list(list(high = gap_order$r, low = gap_order$s, col = 1))),
    if (!is.null(cr)) list(cr = cr), if (!is.null(assoc)) list(assoc = assoc))

# simulated 2x2 data with linear conditional means (LCE holds), as in verify_R10.py
sim2x2 <- function(n = 150, conc = 15, ab = c(0.40, 0.20, 0.55, -0.05), seed = 1, X = NULL, N = NULL) {
  set.seed(seed)
  if (is.null(X)) X <- pmin(pmax(stats::rbeta(n, 2, 2) * 0.8 + 0.1, 0.02), 0.98)
  if (is.null(N)) N <- sample(300:3000, length(X), replace = TRUE)
  m0 <- ab[1] + ab[2] * X; m1 <- ab[3] + ab[4] * X
  b0 <- stats::rbeta(length(X), m0 * conc, (1 - m0) * conc); b1 <- stats::rbeta(length(X), m1 * conc, (1 - m1) * conc)
  list(x = X, t = (1 - X) * b0 + X * b1, n = N, m0 = m0, m1 = m1)
}

test_that("association forms match the Python reference to 1e-9 (W1, W2, pairwise, ER, NM, gamma, D)", {
  for (case in c("main", "gap")) {
    d <- if (case == "main") ex else gp
    st <- fit3(d, ci = "none")$setup
    m <- ref[ref$case == case, ]
    chk <- function(f, nm, r, s = 0) {
      rr <- m[m$name == nm & m$r == r & m$s == s, ]
      expect_equal(c(eiPartialID:::.rxc_val(f, st$beta), f$a), unlist(rr[, c("c", "a1", "a2", "a3")]),
                   tolerance = 1e-9, ignore_attr = TRUE, label = paste(case, nm, r, s))
    }
    for (r in 1:3) {
      chk(eiPartialID:::.rxc_assoc_form(st, "dB", "W1", r, 1), "dB_W1", r)
      chk(eiPartialID:::.rxc_assoc_form(st, "dW", "W1", r, 1), "dW_W1", r)
      chk(eiPartialID:::.rxc_assoc_form(st, "dW", "W2", r, 1), "dW_W2", r)
      chk(eiPartialID:::.rxc_assoc_form(st, "D", "W1", r, 1), "D", r)
      tg <- eiPartialID:::.rxc_target(st, { a <- matrix(0, 3, 2); a[r, 1] <- 1; a })
      chk(list(coef = tg$cvec, const = tg$c0, a = tg$gamma), "Btilde", r)
      dz <- eiPartialID:::.rxc_design_r(st, r)
      ER <- eiPartialID:::.rxc_ER_form(st, r, 1)
      expect_equal(max(abs(ER$a)), 0, tolerance = 1e-12)
      expect_equal(eiPartialID:::.rxc_val(ER, st$beta), m$c[m$name == "ER" & m$r == r], tolerance = 1e-9)
      expect_equal(dz$gamma * eiPartialID:::.rxc_val(ER, st$beta), m$c[m$name == "NM" & m$r == r], tolerance = 1e-9)
      expect_equal(dz$gamma, m$c[m$name == "gamma" & m$r == r], tolerance = 1e-9)
    }
    P <- eiPartialID:::.rxc_pairs(3)
    for (q in seq_len(nrow(P)))
      chk(eiPartialID:::.rxc_assoc_form(st, "dB", "pair", P[q, 1], 1, P[q, 2]), "dB_pair", P[q, 1], P[q, 2])
  }
})

test_that("W1 identities hold for every k: delta_W = Var(ER - D), delta_B = mu(1-mu)(D - NM)", {
  st <- fit3(ex, ci = "none")$setup
  set.seed(3)
  for (r in 1:3) {
    dz <- eiPartialID:::.rxc_design_r(st, r)
    ER <- eiPartialID:::.rxc_val(eiPartialID:::.rxc_ER_form(st, r, 1), st$beta)
    D <- eiPartialID:::.rxc_assoc_form(st, "D", "W1", r, 1)
    dW <- eiPartialID:::.rxc_assoc_form(st, "dW", "W1", r, 1)
    dB <- eiPartialID:::.rxc_assoc_form(st, "dB", "W1", r, 1)
    for (i in 1:5) {
      k <- stats::rnorm(st$nz)
      v <- function(f) eiPartialID:::.rxc_val(f, st$beta) + sum(f$a * k)
      expect_equal(v(dW), dz$var * (ER - v(D)), tolerance = 1e-12)
      expect_equal(v(dB), dz$mu * (1 - dz$mu) * (v(D) - dz$gamma * ER), tolerance = 1e-12)
    }
  }
})

test_that("3x2 CR point bounds equal the Python branch unions (W1 slab = union; W2/pair with a genuine gap)", {
  for (case in c("main", "gap")) {
    d <- if (case == "main") ex else gp
    m <- ref[ref$case == case, ]
    for (v in c("W1", "W2", "pair")) {
      cr <- list(list(row = 1, col = 1, version = v, other = 2))
      f <- fit3(d, ci = "none", restrictions = if (case == "gap") gap_restr(cr) else list(cr = cr))
      for (t in 1:3) {
        ip <- unlist(m[m$name == paste0("CRplus_", v) & m$r == t, c("c", "a1")])
        im <- unlist(m[m$name == paste0("CRminus_", v) & m$r == t, c("c", "a1")])
        u <- eiPartialID:::.rxc_union(list(ip, im))
        got <- f$pieces[[paste0("B", t)]]$estimate
        expect_equal(unname(got), unname(u), tolerance = 1e-8, label = paste(case, v, t))
        # LCE alone
        l0 <- unlist(m[m$name == paste0("LCE_", v) & m$r == t, c("c", "a1")])
        f0 <- if (case == "gap") fit3(d, ci = "none", restrictions = gap_restr()) else fit3(d, ci = "none")
        expect_equal(unname(unlist(f0$table[t, c("estimate_lower", "estimate_upper")])), unname(l0), tolerance = 1e-8)
      }
      # W1 is one slab LP: the R branch union (CR table) equals the slab estimate
      if (v == "W1") {
        tb <- f$sign$cr
        for (t in 1:3) {
          u <- eiPartialID:::.rxc_union(list(unlist(tb[t, c("plus_lower", "plus_upper")]),
                                            unlist(tb[t, c("minus_lower", "minus_upper")])))
          expect_equal(unname(u[1, ]), unname(unlist(f$table[t, c("estimate_lower", "estimate_upper")])), tolerance = 1e-8)
        }
      }
    }
  }
  # the gap case really has two pieces, and the hull option reports the hull
  f <- fit3(gp, ci = "none", restrictions = gap_restr(list(list(row = 1, col = 1, version = "W2"))))
  expect_equal(f$table$n_pieces[2], 2L)
  expect_true(f$table$estimate_gap_upper[2] - f$table$estimate_gap_lower[2] > 0.05)
  expect_false(f$sign$cr$certificate[2])
  fh <- fit3(gp, ci = "none", union = "hull", restrictions = gap_restr(list(list(row = 1, col = 1, version = "W2"))))
  expect_true(is.na(fh$table$estimate_gap_lower[2]))
  expect_equal(fh$table$estimate_lower, f$table$estimate_lower)
  expect_equal(fh$table$estimate_upper, f$table$estimate_upper)
  expect_output(print(f), "U \\[")
})

test_that("noise-free gap case: both CI versions return the two pieces (union of branch CIs)", {
  for (ci in c("conservative", "clr")) {
    f <- fit3(gp, ci = ci, seed = 1, restrictions = gap_restr(list(list(row = 1, col = 1, version = "W2"))))
    pc <- f$pieces$B2$ci_model
    expect_equal(nrow(pc), 2L)
    expect_equal(unname(pc), unname(f$pieces$B2$estimate), tolerance = 1e-6)
  }
})

test_that("2x2 closed forms: delta_W / delta_B sign rows cut D at ER / NM (R2(d))", {
  v <- utils::read.csv(test_path("voter_vax_MonotoneEI.csv"))
  s <- sim2x2(n = 200, seed = 7)
  dats <- list(vv = list(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes),
               sim = s, sim2 = sim2x2(n = 200, seed = 8, ab = c(0.2, 0.5, 0.6, -0.3)))
  for (nm in names(dats)) for (w in c("none", "n")) {
    d <- dats[[nm]]
    fb <- function(res = NULL) eiBounds(x = d$x, t = d$t, n = d$n, weights = w, targets = Dgap, ci = "none",
                                        tests = FALSE, restrictions = res, diagnostics = FALSE)
    f0 <- fb()
    J <- unlist(f0$table[1, c("estimate_lower", "estimate_upper")])
    ds <- eiPartialID:::.rxc_descriptive(f0$setup, 1, f0$rows, f0$cols)
    ER <- ds$ER; NM <- ds$NM
    expect_equal(NM, ds$gamma * ER)
    cases <- list(list("dW", 1, c(J[1], min(J[2], ER))), list("dW", -1, c(max(J[1], ER), J[2])),
                  list("dB", 1, c(max(J[1], NM), J[2])), list("dB", -1, c(J[1], min(J[2], NM))))
    for (cs in cases) {
      f <- fb(list(assoc = list(list(type = cs[[1]], row = 1, col = 1, sign = cs[[2]]))))
      if (cs[[3]][1] > cs[[3]][2] + 1e-12) expect_true(f$table$empty) else
        expect_equal(unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])), unname(cs[[3]]),
                     tolerance = 1e-9, label = paste(nm, w, cs[[1]], cs[[2]]))
    }
    # 2x2 CR is the slab D in hull{NM, ER}, equal to the single non-empty branch
    f <- fb(list(cr = list(list(row = 1, col = 1))))
    ex2 <- c(max(J[1], min(ER, NM)), min(J[2], max(ER, NM)))
    expect_equal(unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])), unname(ex2), tolerance = 1e-9)
    fd <- eiBounds(x = d$x, t = d$t, n = d$n, weights = w, targets = Dgap, ci = "none", tests = FALSE,
                   restrictions = list(cr = list(list(row = 1, col = 1))))
    tb <- fd$sign$cr
    ne <- c(!is.na(tb$plus_lower), !is.na(tb$minus_lower))
    if (sum(ne) == 2) {      # only possible when both reduce to the point D = ER = NM
      expect_equal(tb$plus_lower, tb$minus_upper, tolerance = 1e-9)
    } else {
      br <- if (ne[1]) c(tb$plus_lower, tb$plus_upper) else c(tb$minus_lower, tb$minus_upper)
      expect_equal(br, unname(ex2), tolerance = 1e-9)
    }
  }
})

test_that("MonotoneEI voter_vax: 2x2 CR bounds for the gap are about [-0.470, -0.050]", {
  # data: github.com/reglab/MonotoneEI (MIT), data/voter_vax.rda; see voter_vax_MonotoneEI.README
  v <- utils::read.csv(test_path("voter_vax_MonotoneEI.csv"))
  f <- eiBounds(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes, weights = "n",
                targets = Dgap, ci = "none", tests = FALSE, restrictions = list(cr = list(list(row = 1, col = 1))))
  expect_equal(unname(unlist(f$table[1, c("estimate_lower", "estimate_upper")])), c(-0.470, -0.050), tolerance = 1e-3)
  ds <- f$sign$descriptive
  expect_lt(abs(ds$ER + 0.470), 5e-4); expect_lt(abs(ds$NM + 0.050), 5e-4)
  # tau-based ER equals the T-based EGGHM ER under N-weighted OLS
  p <- v$total_votes / sum(v$total_votes); x <- v$prop_republican
  ERT <- sum(p * (x - sum(p * x)) * v$prop_covid_vax) / sum(p * (x - sum(p * x))^2)
  expect_equal(ds$ER, ERT, tolerance = 1e-10)
  expect_output(print(f$sign), "DESCRIPTIVE ONLY")
})

test_that("redundancy verdicts (face test) agree with brute force", {
  f0 <- fit3(ex, ci = "none")
  st0 <- f0$setup
  cons0 <- eiPartialID:::.rxc_constraints(st0)
  specs <- list()
  for (r in 1:3) for (sg in c(1, -1)) {
    specs[[length(specs) + 1]] <- list(type = "dB", row = r, col = 1, sign = sg)
    specs[[length(specs) + 1]] <- list(type = "dW", row = r, col = 1, sign = sg)
    specs[[length(specs) + 1]] <- list(type = "dW", row = r, col = 1, sign = sg, version = "W2")
    specs[[length(specs) + 1]] <- list(type = "dB", row = r, col = 1, sign = sg, version = "pair", other = r %% 3 + 1)
  }
  seen <- character(0)
  for (sp in specs) {
    f <- fit3(ex, ci = "none", restrictions = list(assoc = list(sp)))
    if (any(f$table$empty)) next
    vd <- f$sign$verdicts
    # brute force: implied iff min over the rung of sign * delta >= 0; tightens iff the interval moves
    a <- eiPartialID:::.rxc_assoc_form(st0, sp$type, if (is.null(sp$version)) "W1" else sp$version, sp$row, 1, sp$other)
    pl <- eiPartialID:::.rxc_plugin(cons0, list(cvec = sp$sign * a$coef, c0 = sp$sign * a$const, gamma = sp$sign * a$a), st0$beta)
    implied <- pl$lower >= -1e-9
    for (t in 1:3) {
      moved <- max(abs(unlist(f$table[t, c("estimate_lower", "estimate_upper")]) -
                       unlist(f0$table[t, c("estimate_lower", "estimate_upper")]))) > 1e-7
      brute <- if (implied) "implied (buys nothing)" else if (moved) "tightens this target" else "cuts the set, not this target"
      expect_equal(vd$verdict_first[t], brute)
      seen <- c(seen, brute)
    }
  }
  expect_true(all(c("tightens this target", "cuts the set, not this target") %in% seen))
  # duals are never used: a row through an optimal vertex that does not cut the face is not "tightens"
  # (the 2x2 case cannot produce the middle verdict, R2(d))
  v <- utils::read.csv(test_path("voter_vax_MonotoneEI.csv"))
  for (sp in list(list(type = "dB", row = 1, col = 1, sign = 1), list(type = "dW", row = 1, col = 1, sign = -1))) {
    f <- eiBounds(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes, weights = "n", targets = Dgap,
                  ci = "none", tests = FALSE, restrictions = list(assoc = list(sp)))
    expect_false(f$sign$verdicts$verdict_first == "cuts the set, not this target")
  }
})

test_that("added-first and added-last gains for a stack of restrictions", {
  f <- fit3(gp, ci = "none", restrictions = gap_restr(list(list(row = 1, col = 1, version = "W2"))))
  vd <- f$sign$verdicts
  expect_setequal(unique(vd$module), c("order[1]", "cr[1]"))
  # added last = full stack vs stack without the module
  f_no_cr <- fit3(gp, ci = "none", restrictions = gap_restr())
  len <- function(u) if (is.null(u)) 0 else sum(u[, 2] - u[, 1])
  for (t in 1:3) {
    g <- len(f_no_cr$pieces[[t]]$estimate) - len(f$pieces[[t]]$estimate)
    expect_equal(vd$gain_last[vd$module == "cr[1]"][t], g, tolerance = 1e-8)
  }
})

test_that("refutation is joint: empty sets name the combination", {
  # voter_vax: ER < NM, so delta_W >= 0 (D <= ER) and delta_B >= 0 (D >= NM) cannot both hold under LCE
  v <- utils::read.csv(test_path("voter_vax_MonotoneEI.csv"))
  f <- eiBounds(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes, weights = "n", targets = Dgap,
                ci = "clr", seed = 1,
                restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1),
                                                 list(type = "dB", row = 1, col = 1, sign = 1))))
  expect_true(f$table$empty)
  expect_true(any(grepl("the rung \\(M0\\) \\+ delta_W\\(X vs rest\\) for 'T' >= 0 \\+ delta_B", f$notes)))
  expect_true(any(grepl("never a restriction on its own", f$notes)))
  expect_true(f$tests$empty_rejected)
  expect_false(f$tests$empty_detail$rung_alone_rejected)
  expect_true(any(grepl("rung alone is not rejected", f$notes)))
  expect_true(any(grepl("not reassurance", f$notes)))
  expect_equal(f$sign$verdicts$verdict_last, rep("empty: refuted jointly with the base", 2))
  # each alone is fine
  for (sp in list(list(type = "dW", row = 1, col = 1, sign = 1), list(type = "dB", row = 1, col = 1, sign = 1))) {
    g <- eiBounds(x = v$prop_republican, t = v$prop_covid_vax, n = v$total_votes, weights = "n", targets = Dgap,
                  ci = "none", tests = FALSE, restrictions = list(assoc = list(sp)), diagnostics = FALSE)
    expect_false(g$table$empty)
  }
})

test_that("CIs with CR: conservative slab and CLR branch union contain the estimate; CLR never falls back", {
  s <- sim2x2(n = 200, seed = 11)
  for (ci in c("conservative", "clr")) {
    f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = ci, seed = 2,
                  restrictions = list(cr = list(list(row = 1, col = 1))))
    tb <- f$table
    expect_false(tb$empty)
    expect_true(tb$ci_model_lower <= tb$estimate_lower && tb$ci_model_upper >= tb$estimate_upper)
    expect_true(tb$ci_actual_lower <= tb$ci_model_lower + 1e-9)
    j <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = ci, seed = 2)$table
    expect_lt(tb$ci_model_upper - tb$ci_model_lower, j$ci_model_upper - j$ci_model_lower)
  }
  # an empty branch under CLR is dropped (ends cross), not replaced by the conservative CI
  st <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = "none", tests = FALSE,
                 restrictions = list(cr = list(list(row = 1, col = 1))), diagnostics = FALSE)$setup
  cons <- eiPartialID:::.rxc_constraints(st)
  sys <- eiPartialID:::.rxc_systems(cons, "branch")
  tg <- eiPartialID:::.rxc_target(st, Dgap$D)
  set.seed(1); Zs <- eiPartialID:::.rxc_zdraws(st$V, 2000)
  emp <- vapply(sys, function(sy) is.null(eiPartialID:::.rxc_plugin(sy, tg, st$beta)$lower), TRUE)
  expect_equal(sum(emp), 1L)
  cc <- eiPartialID:::.rxc_ci_clr(st, sys[[which(emp)]], tg, 0.025, 0.025, Zs, 50)
  expect_true(all(!is.na(cc)))          # ends computed from relaxed candidates
  expect_error(eiPartialID:::.rxc_ci_clr(st, eiPartialID:::.rxc_systems(cons, "slab")[[1]], tg, 0.025, 0.025, Zs, 10),
               "not valid with disjunctive")
})

test_that("sensitivity(), breakdown() and ladder() work with sign and CR rows", {
  s <- sim2x2(n = 150, seed = 5)
  f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = "conservative", seed = 1,
                restrictions = list(cr = list(list(row = 1, col = 1)),
                                    assoc = list(list(type = "dW", row = 1, col = 1, sign = 1))))
  se <- sensitivity(f, eta = c(0, 0.02))
  expect_equal(nrow(se), 2L)
  expect_true(se$estimate_lower[2] <= se$estimate_lower[1] && se$estimate_upper[2] >= se$estimate_upper[1])
  b <- breakdown(f, "D", threshold = -1, direction = "above", which = "estimate", max_eta = 0.2, tol = 0.01)
  expect_true(is.numeric(b) && b >= 0)
  ld <- ladder(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = "none", rungs = c("DD", "M0"),
               sign_rungs = list("M0 + CR" = list(restrictions = list(cr = list(list(row = 1, col = 1)))),
                                 "M0 + dW>=0" = list(restrictions = list(assoc = list(list(type = "dW", row = 1, col = 1, sign = 1))))))
  expect_equal(ld$rung, c("DD", "M0", "M0 + CR", "M0 + dW>=0"))
  expect_output(print(ld), "M0 \\+ CR")
  w <- ld$estimate_upper - ld$estimate_lower
  expect_true(w[3] <= w[2] + 1e-12 && w[4] <= w[2] + 1e-12)
})

test_that("input checks for assoc and cr", {
  s <- sim2x2(n = 60, seed = 2)
  run <- function(res) eiBounds(x = s$x, t = s$t, n = s$n, ci = "none", tests = FALSE, restrictions = res)
  expect_error(run(list(assoc = list(list(type = "dX", row = 1)))), "type")
  expect_error(run(list(assoc = list(list(type = "dW", row = 1, version = "pair")))), "version")
  expect_error(run(list(cr = list(list(row = 1, version = "W3")))), "version")
  expect_warning(run(list(foo = 1)), "unknown restriction")
  f <- run(list(assoc = list(list(type = "dB", row = 1, version = "pair"))))   # 2 rows: other filled in
  expect_equal(f$setup$extra$assoc[[1]]$other, 2L)
})

test_that("coverage simulation (slow; set EIPARTIALID_SLOW=1)", {
  skip_on_cran()
  skip_if(Sys.getenv("EIPARTIALID_SLOW") == "", "set EIPARTIALID_SLOW=1 to run")
  # truth on the CR edge (delta_W = 0), ER far from 0: verify_R10.py design 1
  set.seed(20261105)
  n <- 150
  X <- pmin(pmax(stats::rbeta(n, 2, 2) * 0.8 + 0.1, 0.02), 0.98); N <- sample(300:3000, n, replace = TRUE)
  p <- N / sum(N); mu <- sum(p * X); q1 <- p * X / mu; q0 <- p * (1 - X) / (1 - mu)
  vr <- function(q) sum(q * X^2) - sum(q * X)^2
  a0 <- 0.40; b0 <- 0.20; b1 <- -b0 * (1 - mu) * vr(q0) / (mu * vr(q1))
  a1 <- 0.15 + a0 + b0 * sum(q0 * X) - b1 * sum(q1 * X)
  m0 <- a0 + b0 * X; m1 <- a1 + b1 * X
  D0 <- sum(p * X * m1) / mu - sum(p * (1 - X) * m0) / (1 - mu)
  reps <- as.integer(Sys.getenv("EIPARTIALID_REPS", "200"))
  cov <- c(conservative = 0, clr = 0); wid <- list(conservative = NULL, clr = NULL)
  for (i in seq_len(reps)) {
    s <- sim2x2(X = X, N = N, ab = c(a0, b0, a1, b1), seed = 1000 + i)
    for (ci in names(cov)) {
      f <- eiBounds(x = s$x, t = s$t, n = s$n, targets = Dgap, ci = ci, seed = i, tests = FALSE,
                    diagnostics = FALSE, nz_draw = 2000, ndraw = 50,
                    restrictions = list(cr = list(list(row = 1, col = 1))))
      pc <- f$pieces$D$ci_model
      cov[ci] <- cov[ci] + (!is.null(pc) && any(pc[, 1] - 1e-12 <= D0 & D0 <= pc[, 2] + 1e-12))
      wid[[ci]] <- c(wid[[ci]], if (is.null(pc)) NA else sum(pc[, 2] - pc[, 1]))
    }
  }
  message(sprintf("coverage at the CR edge (n = %d, %d reps): conservative %.3f (median width %.3f), CLR %.3f (median width %.3f)",
                  n, reps, cov[1] / reps, stats::median(wid$conservative, na.rm = TRUE),
                  cov[2] / reps, stats::median(wid$clr, na.rm = TRUE)))
  expect_gte(cov[["conservative"]] / reps, 0.92)
  expect_gte(cov[["clr"]] / reps, 0.90)
})
