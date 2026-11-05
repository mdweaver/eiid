# M-SIGN: signs on mean associations (delta_B, delta_W), contextual reinforcement
# (CR), redundancy diagnostics and descriptive EGGHM quantities (internal).
# Theory: docs/extension/PHASE5a_sign_module.md (R3, R2, R10) and PHASE2c_monotone.md.
#
# Every association is a design-weighted (p_i = N_i / sum N) linear functional of the
# model-implied unit rates m_s(X_i) = e_s' H_c X_i, hence affine in (beta, z):
#   delta = coef' beta + const + a' z.
# A sign restriction is the single linear row sign * delta >= 0 (R10(a)).
# CR = {delta_B >= 0, delta_W >= 0} U {delta_B <= 0, delta_W <= 0} is a union of two
# linear branches (R3(c)); for W1 (and 2x2) it is the slab
#   D^r in hull{gamma_r ER^r, ER^r}
# i.e. two "disjunctive" rows max_l (s_l(beta)) + a' z >= 0 (R3(d), R10(b)).

# sum_i sum_s W[i, s] * m_s(X_i) for column c (unit i uses its own group's H)
.rxc_mlin <- function(st, c, W) {
  out <- list(coef = numeric(st$nbeta), const = 0, a = numeric(st$nz))
  for (g in seq_len(st$G)) {
    sel <- st$g == g
    if (!any(sel)) next
    for (s in seq_len(st$R)) {
      v <- colSums(st$X[sel, , drop = FALSE] * W[sel, s])
      if (all(v == 0)) next
      rt <- .rxc_rate(st, g, c, s, v)
      out$coef <- out$coef + rt$coef; out$const <- out$const + rt$const; out$a <- out$a + rt$a
    }
  }
  out
}

# weight matrices W[i, s] defining each functional (PHASE2c Definition 1, with m for beta)
.rxc_assoc_W <- function(st, type, version = "W1", r, other = NULL) {
  X <- st$X; p <- st$N / sum(st$N); R <- st$R
  mu <- colSums(p * X)
  W <- matrix(0, nrow(X), R)
  xr <- X[, r]
  if (type == "dB") {
    oth <- if (version == "pair") other else setdiff(seq_len(R), r)
    for (s in oth) { w <- p * xr * X[, s]; W[, r] <- W[, r] + w; W[, s] <- W[, s] - w }
  } else if (type == "dW" && version == "W1") {
    xb_r <- sum(p * xr * xr) / mu[r]
    xb_n <- sum(p * (1 - xr) * xr) / (1 - mu[r])
    W[, r] <- p * xr * (xr - xb_r)
    for (s in setdiff(seq_len(R), r)) W[, s] <- p * X[, s] * (xr - xb_n)
  } else if (type == "dW" && version == "W2") {
    for (g in seq_len(R)) {
      xb <- sum(p * X[, g] * xr) / mu[g]
      W[, g] <- p * X[, g] * (xr - xb)
    }
  } else if (type == "D") {                 # D^r = B~_r - B~_{not r}
    W[, r] <- p * xr / mu[r]
    for (s in setdiff(seq_len(R), r)) W[, s] <- -p * X[, s] / (1 - mu[r])
  } else if (type == "covtau") {            # Cov_p(X_r, tau), tau_i = sum_s X_is m_s(X_i)
    for (s in seq_len(R)) W[, s] <- p * (xr - mu[r]) * X[, s]
  } else stop("unknown association type/version: ", type, "/", version)
  W
}

.rxc_assoc_form <- function(st, type, version = "W1", r, c, other = NULL) {
  .rxc_mlin(st, c, .rxc_assoc_W(st, type, version, r, other))
}

# Positive design-only scale that turns an association into probability units:
# delta_B / scale = average within-unit rate gap (weights p_i X_ir X_is);
# delta_W / scale = average within-group slope of the rate on X_r (W1: groups r and
# not-r; W2: every group), i.e. sum_g mu_g Cov_g(m_g, X_r) / sum_g mu_g Var_g(X_r).
# The slack dial (S9) relaxes sign * delta >= 0 to sign * delta >= -slack * scale.
.rxc_assoc_scale <- function(st, type, version = "W1", r, other = NULL) {
  X <- st$X; p <- st$N / sum(st$N); R <- st$R; xr <- unname(X[, r])
  mu <- unname(colSums(p * X))
  wvar <- function(w) { w <- w / sum(w); sum(w * xr^2) - sum(w * xr)^2 }
  if (type == "dB") {
    oth <- if (version == "pair") other else setdiff(seq_len(R), r)
    return(sum(p * xr * rowSums(X[, oth, drop = FALSE])))
  }
  if (version == "W1")
    return(mu[r] * wvar(p * xr) + (1 - mu[r]) * wvar(p * (1 - xr)))
  sum(vapply(seq_len(R), function(g) mu[g] * wvar(p * X[, g]), 0))
}

# design quantities for row r: mu_r, Var_p(X_r), gamma_r
.rxc_design_r <- function(st, r) {
  p <- st$N / sum(st$N); xr <- st$X[, r]
  mu <- sum(p * xr); v <- sum(p * (xr - mu)^2)
  list(mu = mu, var = v, gamma = v / (mu * (1 - mu)))
}

# ER^r (tau-based) as a linear form in beta (its z-part is identically zero)
.rxc_ER_form <- function(st, r, c) {
  f <- .rxc_assoc_form(st, "covtau", "W1", r, c)
  dz <- .rxc_design_r(st, r)
  list(coef = f$coef / dz$var, const = f$const / dz$var, a = f$a / dz$var)
}

.rxc_val <- function(f, beta) sum(f$coef * beta) + f$const

# the pair (g_B, g_W) of a CR restriction
.rxc_cr_pair <- function(st, cr) {
  if (cr$version == "W1") {
    gB <- .rxc_assoc_form(st, "dB", "W1", cr$row, cr$col)
    gW <- .rxc_assoc_form(st, "dW", "W1", cr$row, cr$col)
  } else if (cr$version == "W2") {
    gB <- .rxc_assoc_form(st, "dB", "W1", cr$row, cr$col)
    gW <- .rxc_assoc_form(st, "dW", "W2", cr$row, cr$col)
  } else {
    gB <- .rxc_assoc_form(st, "dB", "pair", cr$row, cr$col, cr$other)
    gW <- .rxc_assoc_form(st, "dW", "W2", cr$row, cr$col)
  }
  list(gB = gB, gW = gW)
}

.rxc_rows <- function(fs, sg = 1) {
  list(S = do.call(rbind, lapply(fs, function(f) sg * f$coef)),
       s0 = sg * vapply(fs, function(f) f$const, 0),
       A = do.call(rbind, lapply(fs, function(f) sg * f$a)))
}

# CR objects for cons$cr (empty list if none)
.rxc_cr_parts <- function(st) {
  ex <- st$extra
  if (is.null(ex$cr) || st$nz == 0) return(list())
  lapply(seq_along(ex$cr), function(i) {
    cr <- ex$cr[[i]]
    pr <- .rxc_cr_pair(st, cr)
    sl <- if (is.null(cr$slack)) 0 else cr$slack
    # slack (S9): branch + is {g_B >= -eta s_B, g_W >= -eta s_W}, branch - is {g_B <= eta s_B, g_W <= eta s_W}
    sB <- .rxc_assoc_scale(st, "dB", if (cr$version == "pair") "pair" else "W1", cr$row, cr$other)
    sW <- .rxc_assoc_scale(st, "dW", if (cr$version == "W1") "W1" else "W2", cr$row)
    plus <- .rxc_rows(list(pr$gB, pr$gW), 1); plus$s0 <- plus$s0 + sl * c(sB, sW)
    minus <- .rxc_rows(list(pr$gB, pr$gW), -1); minus$s0 <- minus$s0 + sl * c(sB, sW)
    out <- list(module = paste0("cr[", i, "]"), spec = cr,
                gB = pr$gB, gW = pr$gW, plus = plus, minus = minus, slab = NULL)
    if (cr$version == "W1") {
      # slab D in hull{gamma ER, ER}: lower max(D - ER, D - gamma ER) >= 0,
      # upper max(ER - D, gamma ER - D) >= 0 (R10(b)); pieces share the k-normal.
      # With slack: branch + is D in [NM - a, ER + b], branch - is D in [ER - b, NM + a],
      # a = eta s_B / (mu(1-mu)), b = eta s_W / Var(X_r); two non-empty branches always
      # intersect, so the union is [min(NM - a, ER - b), max(ER + b, NM + a)] (PHASE5a §5.2)
      dz <- .rxc_design_r(st, cr$row)
      D <- .rxc_assoc_form(st, "D", "W1", cr$row, cr$col)
      ER <- .rxc_ER_form(st, cr$row, cr$col)
      gm <- dz$gamma
      aa <- sl * sB / (dz$mu * (1 - dz$mu)); bb <- sl * sW / dz$var
      lo <- list(S = rbind(D$coef - ER$coef, D$coef - gm * ER$coef),
                 s0 = c(D$const - ER$const + bb, D$const - gm * ER$const + aa), a = D$a)
      up <- list(S = rbind(ER$coef - D$coef, gm * ER$coef - D$coef),
                 s0 = c(ER$const - D$const + bb, gm * ER$const - D$const + aa), a = -D$a)
      out$slab <- list(lo, up)
      out$D <- D; out$ER <- ER; out$gamma <- gm
    }
    out
  })
}

# evaluate disjunctive rows at beta (optionally inflated piecewise) -> ordinary rows
.rxc_fix_disj <- function(cons, beta, infl = NULL) {
  for (j in seq_along(cons$disj)) {
    d <- cons$disj[[j]]
    v <- as.vector(d$S %*% beta) + d$s0
    if (!is.null(infl)) v <- v + infl[[j]]
    cons$S <- rbind(cons$S, numeric(length(beta)))
    cons$s0 <- c(cons$s0, max(v))
    cons$A <- rbind(cons$A, d$a)
    cons$kind <- c(cons$kind, "cr_slab")
    cons$module <- c(cons$module, d$module)
  }
  cons$disj <- NULL
  cons
}

.rxc_append_rows <- function(cons, rows, kind, module) {
  cons$S <- rbind(cons$S, rows$S); cons$s0 <- c(cons$s0, rows$s0)
  cons$A <- rbind(cons$A, rows$A)
  cons$kind <- c(cons$kind, rep(kind, length(rows$s0)))
  cons$module <- c(cons$module, rep(module, length(rows$s0)))
  cons
}

# Expand CR restrictions into a list of systems whose union is the restricted set.
# mode "slab": W1 CRs as disjunctive slab rows (point bounds, conservative CI);
# mode "branch": every CR as its two linear branches (CLR CI, R10(c)).
# 'keep' (module names) restricts which optional modules are included.
.rxc_systems <- function(cons, mode = c("slab", "branch"), keep = NULL) {
  mode <- match.arg(mode)
  crs <- cons$cr
  base <- cons; base$cr <- NULL; base$branch <- ""
  if (!is.null(keep)) {
    ok <- is.na(base$module) | base$module %in% keep
    base$S <- base$S[ok, , drop = FALSE]; base$s0 <- base$s0[ok]
    base$A <- base$A[ok, , drop = FALSE]; base$kind <- base$kind[ok]; base$module <- base$module[ok]
    crs <- Filter(function(cr) cr$module %in% keep, crs)
  }
  sys <- list(base)
  for (cr in crs) {
    if (mode == "slab" && !is.null(cr$slab)) {
      sl <- lapply(cr$slab, function(d) { d$module <- cr$module; d })
      sys <- lapply(sys, function(s) { s$disj <- c(s$disj, sl); s })
    } else {
      new <- list()
      for (s in sys) for (br in c("plus", "minus")) {
        s2 <- .rxc_append_rows(s, cr[[br]], "cr", cr$module)
        s2$branch <- paste0(s$branch, if (nzchar(s$branch)) ", " else "", cr$module,
                            if (br == "plus") " (+)" else " (-)")
        new[[length(new) + 1]] <- s2
      }
      sys <- new
    }
  }
  sys
}

# union of closed intervals (list of c(lo, hi)); NULL if all empty
.rxc_union <- function(iv, tol = 1e-9) {
  iv <- Filter(function(x) length(x) == 2 && !any(is.na(x)) && x[1] <= x[2] + tol, iv)
  if (!length(iv)) return(NULL)
  m <- do.call(rbind, iv); m[, 2] <- pmax(m[, 1], m[, 2])
  m <- m[order(m[, 1]), , drop = FALSE]
  out <- m[1, , drop = FALSE]
  for (i in seq_len(nrow(m))[-1]) {
    k <- nrow(out)
    if (m[i, 1] <= out[k, 2] + tol) out[k, 2] <- max(out[k, 2], m[i, 2]) else out <- rbind(out, m[i, ])
  }
  colnames(out) <- c("lower", "upper"); rownames(out) <- NULL
  out
}

# summary of a union: hull, number of pieces, widest gap
.rxc_union_summary <- function(u) {
  if (is.null(u)) return(list(lower = NA_real_, upper = NA_real_, pieces = 0L,
                              gap_lower = NA_real_, gap_upper = NA_real_))
  gl <- gu <- NA_real_
  if (nrow(u) > 1) {
    gaps <- u[-1, 1] - u[-nrow(u), 2]
    j <- which.max(gaps); gl <- u[j, 2]; gu <- u[j + 1, 1]
  }
  list(lower = min(u[, 1]), upper = max(u[, 2]), pieces = nrow(u), gap_lower = gl, gap_upper = gu)
}

# ---------------------------------------------------------------------------
# Primal LPs in z (free variables, split as z = u - v for lpSolve)
# rows: s + A z >= 0.  Returns status 0 optimal, 2 infeasible, 3 unbounded.
# ---------------------------------------------------------------------------
.rxc_lp_z <- function(obj, s, A, sense = "max") {
  nv <- ncol(A)
  if (nv == 0) return(list(status = if (all(s >= -1e-9)) 0 else 2, value = 0, z = numeric(0)))
  sol <- lpSolve::lp(sense, objective.in = c(obj, -obj), const.mat = cbind(A, -A),
                     const.dir = rep(">=", nrow(A)), const.rhs = -s)
  z <- sol$solution[seq_len(nv)] - sol$solution[nv + seq_len(nv)]
  list(status = sol$status, value = sum(obj * z), z = z)
}

# largest t such that F (rows sF + AF z >= 0) meets {sN + AN z >= t} (t capped at 1).
# t >= -tol  <=>  F meets N (one feasibility LP, R2(a)).
.rxc_meet <- function(sF, AF, sN, AN) {
  nz <- ncol(AF)
  A <- rbind(cbind(AF, rep(0, nrow(AF))), cbind(AN, rep(-1, nrow(AN))), c(numeric(nz), -1))
  s <- c(sF, sN, 1)
  r <- .rxc_lp_z(c(numeric(nz), 1), s, A)
  if (r$status == 2) return(-Inf)
  r$value
}

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
.rxc_assoc_label <- function(a, rl, cl) {
  nm <- if (a$type == "dB") {
    if (a$version == "pair") sprintf("delta_B(%s, %s)", rl[a$row], rl[a$other])
    else sprintf("delta_B(%s vs rest)", rl[a$row])
  } else {
    if (a$version == "W2") sprintf("Delta_W(%s) [W2]", rl[a$row])
    else sprintf("delta_W(%s vs rest)", rl[a$row])
  }
  sl <- if (is.null(a$slack) || a$slack == 0) "" else sprintf(" (slack %g)", a$slack)
  sprintf("%s for '%s' %s 0%s", nm, cl[a$col], if (a$sign > 0) ">=" else "<=", sl)
}

.rxc_cr_label <- function(cr, rl, cl) {
  pr <- switch(cr$version,
               W1 = sprintf("delta_B(%s vs rest), delta_W(%s vs rest)", rl[cr$row], rl[cr$row]),
               W2 = sprintf("delta_B(%s vs rest), Delta_W(%s) [W2]", rl[cr$row], rl[cr$row]),
               pair = sprintf("delta_B(%s, %s), Delta_W(%s) [W2]", rl[cr$row], rl[cr$other], rl[cr$row]))
  sl <- if (is.null(cr$slack) || cr$slack == 0) "" else sprintf(" (slack %g)", cr$slack)
  sprintf("CR for '%s': same sign of %s%s", cl[cr$col], pr, sl)
}

# ---------------------------------------------------------------------------
# Point-level helpers on systems (at beta-hat)
# ---------------------------------------------------------------------------
# interval of a target over a list of systems (union); NULL if all empty
.rxc_sys_interval <- function(sys, tg, beta) {
  iv <- list(); unb <- FALSE
  for (s in sys) {
    pl <- .rxc_plugin(s, tg, beta)
    if (isTRUE(pl$empty)) next
    if (isTRUE(pl$unbounded)) { unb <- TRUE; next }
    iv[[length(iv) + 1]] <- c(pl$lower, pl$upper)
  }
  list(union = .rxc_union(iv), unbounded = unb, per_system = iv)
}

.rxc_fixed <- function(s, beta) {
  if (length(s$disj)) s <- .rxc_fix_disj(s, beta)
  list(s = as.vector(s$S %*% beta) + s$s0, A = s$A)
}

# the rows of module 'md' evaluated at beta, as a list of alternatives (a union):
# linear rows / W1 slab -> one alternative; other CR -> two (its branches)
.rxc_module_alts <- function(cons, md, beta) {
  cr <- Filter(function(x) x$module == md, cons$cr)
  if (length(cr)) {
    cr <- cr[[1]]
    if (!is.null(cr$slab)) {
      tmp <- list(S = matrix(0, 0, length(beta)), s0 = numeric(0), A = matrix(0, 0, ncol(cons$A)),
                  disj = cr$slab)
      f <- .rxc_fixed(tmp, beta)
      return(list(f))
    }
    return(lapply(c("plus", "minus"), function(br)
      list(s = as.vector(cr[[br]]$S %*% beta) + cr[[br]]$s0, A = cr[[br]]$A)))
  }
  ok <- !is.na(cons$module) & cons$module == md
  list(list(s = as.vector(cons$S[ok, , drop = FALSE] %*% beta) + cons$s0[ok],
            A = cons$A[ok, , drop = FALSE]))
}

# is module N (alternatives) implied by system F, i.e. F subset of N?  (R2(c))
.rxc_implied <- function(Ff, alts, tol) {
  if (length(alts) == 1) {
    N <- alts[[1]]
    for (j in seq_along(N$s)) {
      r <- .rxc_lp_z(N$A[j, ], Ff$s, Ff$A, sense = "min")
      if (r$status == 3) return(FALSE)
      if (r$status == 0 && r$value + N$s[j] < -tol) return(FALSE)
    }
    return(TRUE)
  }
  # CR union {gB >= 0, gW >= 0} U {gB <= 0, gW <= 0}: not implied iff F has a point
  # with gB > 0 > gW or gB < 0 < gW
  P <- alts[[1]]
  for (sg in c(1, -1)) {
    t <- .rxc_meet(Ff$s, Ff$A, c(sg * P$s[1], -sg * P$s[2]), rbind(sg * P$A[1, ], -sg * P$A[2, ]))
    if (t > tol) return(FALSE)
  }
  TRUE
}

# R2(a) face test for one endpoint of a target on a single (convex) system F:
# TRUE if the endpoint moves when N is added (the optimal face misses N)
.rxc_face_moves <- function(Ff, tg, beta, alts, side, value, tol) {
  h <- if (side == "upper") tg$gamma else -tg$gamma
  cb <- sum(tg$cvec * beta) + tg$c0
  # face: h'z >= value - cb (upper) or -h'z >= -(value - cb) (lower)
  lev <- if (side == "upper") value - cb else -(value - cb)
  sface <- c(Ff$s, -lev + tol); Aface <- rbind(Ff$A, h)
  meets <- vapply(alts, function(N) .rxc_meet(sface, Aface, N$s, N$A) >= -tol, logical(1))
  !any(meets)
}

# ---------------------------------------------------------------------------
# Diagnostics: redundancy verdicts (R2), "already implied" table, CR branches /
# certificate (R3(c)), descriptive EGGHM quantities
# ---------------------------------------------------------------------------
.rxc_union_len <- function(u) if (is.null(u)) 0 else sum(u[, 2] - u[, 1])

.rxc_verdict <- function(cons, md, keepF, tg, beta, tol) {
  Fsys <- .rxc_systems(cons, "slab", keep = keepF)
  Nsys <- .rxc_systems(cons, "slab", keep = c(keepF, md))
  IF <- .rxc_sys_interval(Fsys, tg, beta); IN <- .rxc_sys_interval(Nsys, tg, beta)
  hF <- .rxc_union_summary(IF$union); hN <- .rxc_union_summary(IN$union)
  res <- list(verdict = NA_character_, lower_gain = NA_real_, upper_gain = NA_real_,
              width_gain = NA_real_, method = "face")
  if (is.null(IF$union)) { res$verdict <- "base set empty"; return(res) }
  if (is.null(IN$union)) { res$verdict <- "empty: refuted jointly with the base"; return(res) }
  if (IF$unbounded || IN$unbounded) { res$verdict <- "unbounded"; return(res) }
  res$lower_gain <- hN$lower - hF$lower; res$upper_gain <- hF$upper - hN$upper
  res$width_gain <- .rxc_union_len(IF$union) - .rxc_union_len(IN$union)
  alts <- .rxc_module_alts(cons, md, beta)
  Ffix <- lapply(Fsys, .rxc_fixed, beta = beta)
  if (all(vapply(Ffix, function(Ff) .rxc_implied(Ff, alts, tol), logical(1)))) {
    res$verdict <- "implied (buys nothing)"; return(res)
  }
  if (length(Fsys) == 1) {
    # R2(a): an endpoint moves iff the restriction removes the whole optimal face
    up <- .rxc_face_moves(Ffix[[1]], tg, beta, alts, "upper", hF$upper, tol)
    lo <- .rxc_face_moves(Ffix[[1]], tg, beta, alts, "lower", hF$lower, tol)
    moves <- up || lo
  } else {
    res$method <- "brute (base is a union)"
    moves <- res$width_gain > tol * max(1, abs(hF$upper - hF$lower))
  }
  res$verdict <- if (moves) "tightens this target" else "cuts the set, not this target"
  res
}

.rxc_all_modules <- function(cons) {
  unique(c(stats::na.omit(cons$module), vapply(cons$cr, function(x) x$module, "")))
}

# redundancy table for a setup and list of targets (point level, at beta-hat)
.rxc_redundancy <- function(st, targets, labels, tol = 1e-7) {
  cons <- .rxc_constraints(st)
  mods <- .rxc_all_modules(cons)
  if (!length(mods) || st$nz == 0) return(NULL)
  out <- list()
  for (md in mods) for (nm in names(targets)) {
    tg <- .rxc_target(st, targets[[nm]])
    f1 <- .rxc_verdict(cons, md, character(0), tg, st$beta, tol)
    f2 <- if (length(mods) > 1) .rxc_verdict(cons, md, setdiff(mods, md), tg, st$beta, tol) else f1
    out[[length(out) + 1]] <- data.frame(
      module = md, restriction = unname(labels[md]), target = nm,
      verdict_first = f1$verdict, gain_first = f1$width_gain,
      lower_gain_first = f1$lower_gain, upper_gain_first = f1$upper_gain,
      verdict_last = f2$verdict, gain_last = f2$width_gain, method = f1$method,
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}

# delta's under the rung alone: point interval and CI, flagged where they exclude 0
.rxc_implied_table <- function(st, cols, rl, cl, ci = "clr", alpha = 0.05, ndraw = 100,
                               nz_draw = 4000) {
  if (st$nz == 0) return(NULL)
  st0 <- st; st0$extra <- NULL
  cons0 <- .rxc_constraints(st0)
  Zs <- if (ci != "none") .rxc_zdraws(st$V, nz_draw) else NULL
  R <- st$R
  items <- list()
  for (c in cols) {
    for (r in if (R == 2) 1 else seq_len(R)) {
      items[[length(items) + 1]] <- list(nm = sprintf("delta_B(%s vs rest)", rl[r]), type = "dB", v = "W1", r = r, c = c)
      items[[length(items) + 1]] <- list(nm = sprintf("delta_W(%s vs rest)", rl[r]), type = "dW", v = "W1", r = r, c = c)
      if (R > 2) items[[length(items) + 1]] <- list(nm = sprintf("Delta_W(%s) [W2]", rl[r]), type = "dW", v = "W2", r = r, c = c)
    }
    if (R > 2) {
      P <- .rxc_pairs(R)
      for (q in seq_len(nrow(P))) items[[length(items) + 1]] <-
        list(nm = sprintf("delta_B(%s, %s)", rl[P[q, 1]], rl[P[q, 2]]), type = "dB", v = "pair",
             r = P[q, 1], c = c, other = P[q, 2])
    }
  }
  rows <- lapply(items, function(it) {
    f <- .rxc_assoc_form(st0, it$type, it$v, it$r, it$c, it$other)
    sc <- .rxc_assoc_scale(st0, it$type, it$v, it$r, it$other)      # report in probability units
    tg <- list(cvec = f$coef / sc, c0 = f$const / sc, gamma = f$a / sc)
    pl <- .rxc_plugin(cons0, tg, st0$beta)
    lo <- hi <- NA_real_
    if (!isTRUE(pl$empty) && !isTRUE(pl$unbounded)) { lo <- pl$lower; hi <- pl$upper }
    cil <- ciu <- NA_real_
    if (ci == "clr") { cc <- .rxc_ci_clr(st0, cons0, tg, alpha / 2, alpha / 2, Zs, ndraw); cil <- cc[1]; ciu <- cc[2] }
    if (ci == "conservative") { cc <- .rxc_ci_conservative(st0, cons0, tg, alpha, Zs)$ci; cil <- cc[1]; ciu <- cc[2] }
    sgn <- function(l, u) if (is.na(l) || is.na(u) || l > u) NA_character_ else
      if (l > 0) "positive" else if (u < 0) "negative" else "undetermined"
    data.frame(association = it$nm, column = cl[it$c],
               scale = if (it$type == "dB") "average gap" else "average slope", lower = lo, upper = hi,
               ci_lower = cil, ci_upper = ciu, sign_point = sgn(lo, hi), sign_ci = sgn(cil, ciu),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# descriptive EGGHM quantities (tau-based ER, NM = gamma ER, bracket) -- not restrictions
.rxc_descriptive <- function(st, cols, rl, cl) {
  R <- st$R
  rows <- list()
  st0 <- st; st0$extra <- NULL
  cons0 <- if (st$nz > 0) .rxc_constraints(st0) else NULL
  for (c in cols) for (r in if (R == 2) 1 else seq_len(R)) {
    dz <- .rxc_design_r(st, r)
    ER <- .rxc_val(.rxc_ER_form(st, r, c), st$beta)
    NM <- dz$gamma * ER
    D <- .rxc_assoc_form(st, "D", "W1", r, c)
    dl <- du <- NA_real_
    if (!is.null(cons0)) {
      pl <- .rxc_plugin(cons0, list(cvec = D$coef, c0 = D$const, gamma = D$a), st$beta)
      if (!isTRUE(pl$empty) && !isTRUE(pl$unbounded)) { dl <- pl$lower; du <- pl$upper }
    } else dl <- du <- .rxc_val(D, st$beta)
    rows[[length(rows) + 1]] <- data.frame(
      row = rl[r], column = cl[c], mu = dz$mu, var = dz$var, gamma = dz$gamma,
      ER = ER, NM = NM, bracket_lower = min(ER, NM), bracket_upper = max(ER, NM),
      D_lower = dl, D_upper = du, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# CR branches (with all other restrictions) and the connectivity certificate
.rxc_cr_table <- function(st, targets, labels, tol = 1e-7) {
  cons <- .rxc_constraints(st)
  if (!length(cons$cr)) return(NULL)
  mods <- .rxc_all_modules(cons)
  out <- list()
  for (cr in cons$cr) {
    others <- .rxc_systems(cons, "slab", keep = setdiff(mods, cr$module))
    br <- lapply(c("plus", "minus"), function(b)
      lapply(others, function(s) .rxc_append_rows(s, cr[[b]], "cr", cr$module)))
    # certificate: Theta (others) meets {g_B = 0, g_W = 0}
    eq <- list(s = c(cr$plus$s0, cr$minus$s0), A = rbind(cr$plus$A, cr$minus$A),
               S = rbind(cr$plus$S, cr$minus$S))
    cert <- any(vapply(others, function(s) {
      Ff <- .rxc_fixed(s, st$beta)
      .rxc_meet(Ff$s, Ff$A, as.vector(eq$S %*% st$beta) + eq$s, eq$A) >= -tol
    }, logical(1)))
    for (nm in names(targets)) {
      tg <- .rxc_target(st, targets[[nm]])
      ip <- .rxc_sys_interval(br[[1]], tg, st$beta)$union
      im <- .rxc_sys_interval(br[[2]], tg, st$beta)$union
      u <- .rxc_union(c(if (!is.null(ip)) split(ip, row(ip)) else list(),
                        if (!is.null(im)) split(im, row(im)) else list()))
      sm <- .rxc_union_summary(u)
      out[[length(out) + 1]] <- data.frame(
        module = cr$module, restriction = unname(labels[cr$module]), target = nm,
        form = if (is.null(cr$slab)) "union of branches" else "slab (one LP)",
        plus_lower = if (is.null(ip)) NA else min(ip[, 1]), plus_upper = if (is.null(ip)) NA else max(ip[, 2]),
        minus_lower = if (is.null(im)) NA else min(im[, 1]), minus_upper = if (is.null(im)) NA else max(im[, 2]),
        certificate = cert, pieces = sm$pieces, stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out); rownames(res) <- NULL
  res
}

#' Diagnostics for sign and contextual-reinforcement restrictions
#'
#' For an \code{eiRxC} fit, \code{eiRedundancy()} reports what each optional
#' restriction buys (PHASE5a, R2). For each restriction and target it gives a
#' verdict, computed by the face test (one feasibility LP per endpoint; solver
#' dual values are never used for the verdict):
#' \code{"implied (buys nothing)"} (the data and the rung already imply it, for
#' every target), \code{"cuts the set, not this target"}, or
#' \code{"tightens this target"}. Verdicts and width gains are reported with
#' the restriction added first (to the rung alone) and added last (to all the
#' other restrictions).
#'
#' It also returns the "already implied by the data" table (each association's
#' interval and CI under the rung alone, flagged where it excludes 0), the CR
#' table (both branch intervals and the connectivity certificate) and the
#' descriptive EGGHM quantities: the tau-based \eqn{\widehat{ER}},
#' \eqn{\widehat{NM} = \gamma\widehat{ER}}, \eqn{\gamma} and the bracket
#' \eqn{[\gamma\widehat{ER}, \widehat{ER}]}. The bracket is descriptive only:
#' never intersect it with a confidence interval (that under-covers, PHASE5a
#' R10(e)). Our \eqn{\widehat{ER}} uses the fitted \eqn{E[T|X]}; it equals the
#' EGGHM \eqn{T}-based estimate when the regression is \eqn{N}-weighted
#' (\code{weights = "n"}).
#'
#' All verdicts are estimated, sharp-edged classifications at the point
#' estimate; near a tie compare the width change with its uncertainty.
#'
#' @param object an \code{eiRxC} object.
#' @param targets names of targets (default: all).
#' @param cols columns for the association tables (default: those named in the
#'   restrictions, else the first).
#' @param ci compute CIs in the "already implied" table (uses the fit's CI method).
#' @param tol numerical tolerance for verdicts.
#' @return a list of class \code{"eiRedundancy"} with data frames
#'   \code{verdicts}, \code{implied}, \code{cr} and \code{descriptive}.
#' @export
eiRedundancy <- function(object, targets = NULL, cols = NULL, ci = TRUE, tol = 1e-7) {
  st <- object$setup
  tl <- if (is.null(targets)) object$targets else object$targets[targets]
  ex <- st$extra
  if (is.null(cols)) {
    cols <- unique(c(vapply(ex$assoc, function(a) a$col, 0), vapply(ex$cr, function(a) a$col, 0)))
    if (!length(cols)) cols <- 1
  }
  labels <- .rxc_module_labels(ex, object$rows, object$cols)
  structure(list(
    verdicts = .rxc_redundancy(st, tl, labels, tol),
    implied = .rxc_implied_table(st, cols, object$rows, object$cols,
                                 ci = if (ci) object$ci else "none", alpha = object$alpha,
                                 ndraw = object$settings$ndraw, nz_draw = object$settings$nz_draw),
    cr = .rxc_cr_table(st, tl, labels, tol),
    descriptive = .rxc_descriptive(st, cols, object$rows, object$cols),
    weights = object$settings$weights, model = object$model), class = "eiRedundancy")
}

#' @export
print.eiRedundancy <- function(x, digits = 3, ...) {
  r <- function(d) { if (is.null(d)) return(d); n <- vapply(d, is.numeric, TRUE); d[n] <- lapply(d[n], round, digits); d }
  if (!is.null(x$verdicts)) {
    cat("What each restriction buys (face test at the point estimate; PHASE5a R2)\n")
    print(r(x$verdicts[, c("restriction", "target", "verdict_first", "gain_first", "verdict_last", "gain_last")]),
          row.names = FALSE, right = FALSE)
    cat("  gain = reduction in the width of the estimated set; first = added to the rung alone,\n")
    cat("  last = added to all other restrictions\n\n")
  }
  if (!is.null(x$cr)) {
    cat("Contextual reinforcement: branch intervals and connectivity certificate\n")
    print(r(x$cr[, c("restriction", "target", "form", "plus_lower", "plus_upper", "minus_lower",
                     "minus_upper", "certificate", "pieces")]), row.names = FALSE, right = FALSE)
    cat("\n")
  }
  if (!is.null(x$implied)) {
    cat("Associations under the rung alone (signs already implied by the data + rung),\n")
    cat("in probability units: average within-unit gap (delta_B) or average slope in own share (delta_W)\n")
    print(r(x$implied), row.names = FALSE, right = FALSE)
    cat("\n")
  }
  if (!is.null(x$descriptive)) {
    cat("DESCRIPTIVE ONLY (not restrictions, not confidence intervals): EGGHM quantities\n")
    print(r(x$descriptive), row.names = FALSE, right = FALSE)
    cat("  ER, NM are tau-based (fitted E[T|X]); they equal the EGGHM T-based estimates under\n")
    cat("  N-weighted regression", if (identical(x$weights, "n")) " (used here).\n" else " (not used here: weights = \"none\").\n", sep = "")
    cat("  Never intersect the bracket [gamma ER, ER] with a confidence interval: it under-covers.\n")
  }
  invisible(x)
}
