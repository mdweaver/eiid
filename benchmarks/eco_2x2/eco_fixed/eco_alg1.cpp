// eco's DP-mixture Gibbs sampler for 2x2 ecological inference (ecoNP, CAR and
// NCAR; Imai, Lu & Strauss 2008, Appendix B; eco src/gibbsDP.c), re-implemented
// with efficient code and with eco's two departures from its stated model fixed:
//   B1: covariance draws use the inverse-Wishart with the stated degrees of
//       freedom (eco's rWish draws Wishart(df - 1));
//   B2: the Escobar-West mixture weight for alpha is a probability, not the odds.
// Everything else follows eco, step by step, per sweep:
//   1. each latent unit's rates: one independence Metropolis-Hastings step with a
//      uniform proposal on the tomography line (eco's rMH, the default), given the
//      unit's current (mu, Sigma); optionally a slice step instead (wstep = 1);
//   2. partially observed units: missing coordinate from its conditional normal;
//   3. configuration (Neal 2000 Algorithm 1): for each unit, choose another unit's
//      parameters with weight N(W_i* | theta_j), or new parameters with weight
//      alpha * t_0(W_i*), drawn from the NIW posterior given W_i* alone;
//   4. "remixing" (Bush & MacEachern 1996): redraw each cluster's (mu, Sigma)
//      from its NIW posterior given its members;
//   5. alpha by Escobar & West (1995).
// Units sharing parameters have identical weights in step 3, so the sum over
// units j != i equals a sum over clusters weighted by their sizes. grouped = 1
// uses the cluster form (cost n * k); grouped = 0 loops over units as eco does
// (cost n^2). Both target the same distribution.
// The function signature and output match eiDPM's dpm_collapsed_cpp.
//
// eco_param_cpp (end of file) does the same for eco's parametric Bayesian models
// (eco(), CAR and NCAR; src/gibbsBase.c, src/gibbsXBase.c): rate step as above,
// partially observed units, then (mu, Sigma) from the NIW posterior given all
// units, with B1 fixed. eco's other inefficiencies are gone in both: no memory
// allocation or Cholesky factorisation per density evaluation (eco's rMH even
// computes the determinant, which cancels in the MH ratio).
// Two further eco bugs, both for units with X exactly 0 or 1, do not arise here:
//   B3 (ecoNP CAR, gibbsDP.c line 295): the draw of the missing W1* of the i-th
//      X = 0 unit is written into slot n_samp + i (an X = 1 unit's observed W1*),
//      so X = 0 units' W1 is never updated and X = 1 units' data are overwritten;
//   B4 (eco and ecoNP, R/checkdata.R): units are ordered [regular, X=0, X=1] for
//      relabelling but the C code stores [regular, X=1, X=0], so the returned
//      draws of the X = 0 and X = 1 units are attached to the wrong units.
// Here unit bookkeeping is done by eiDPM's dpmData, which keeps them apart.
#include <Rcpp.h>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

namespace {

inline double logitc(double w) { return std::log(w) - std::log1p(-w); }

struct Par {        // cluster parameters
  double mu[3];
  double B[3][3];   // lower triangular, Sigma^{-1} = B B'
  double Sig[3][3];
  double cst;       // -p/2 log(2 pi) + sum log B_kk
};

struct Prior { int p; double mu0[3], tau0, nu0, S0[3][3]; };

struct T0 { double df, cst, loc[3], Li[3][3]; };  // prior predictive t

void chol(const double A[3][3], double L[3][3], int p) {
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) L[i][j] = 0;
  for (int j = 0; j < p; j++) {
    double v = A[j][j];
    for (int k = 0; k < j; k++) v -= L[j][k] * L[j][k];
    L[j][j] = std::sqrt(v);
    for (int i = j + 1; i < p; i++) {
      double w = A[i][j];
      for (int k = 0; k < j; k++) w -= L[i][k] * L[j][k];
      L[i][j] = w / L[j][j];
    }
  }
}
void inv_lower(const double L[3][3], double Li[3][3], int p) {
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Li[i][j] = 0;
  for (int j = 0; j < p; j++) {
    Li[j][j] = 1.0 / L[j][j];
    for (int i = j + 1; i < p; i++) {
      double w = 0;
      for (int k = j; k < i; k++) w -= L[i][k] * Li[k][j];
      Li[i][j] = w / L[i][i];
    }
  }
}

T0 make_t0(const Prior &pr) {
  const int p = pr.p;
  T0 t;
  t.df = pr.nu0 - p + 1;
  double sc = (pr.tau0 + 1) / (pr.tau0 * t.df), A[3][3], L[3][3];
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) A[i][j] = sc * pr.S0[i][j];
  chol(A, L, p);
  inv_lower(L, t.Li, p);
  double ldet = 0;
  for (int j = 0; j < p; j++) ldet += 2 * std::log(L[j][j]);
  for (int i = 0; i < 3; i++) t.loc[i] = pr.mu0[i];
  t.cst = std::lgamma(0.5 * (t.df + p)) - std::lgamma(0.5 * t.df) - 0.5 * p * std::log(t.df * M_PI) - 0.5 * ldet;
  return t;
}
inline double lt0(const double *x, const T0 &t, int p) {
  double q = 0;
  for (int i = 0; i < p; i++) {
    double z = 0;
    for (int j = 0; j <= i; j++) z += t.Li[i][j] * (x[j] - t.loc[j]);
    q += z * z;
  }
  return t.cst - 0.5 * (t.df + p) * std::log1p(q / t.df);
}
inline double lnorm(const double *x, const Par &c, int p) {
  double q = 0;
  for (int k = 0; k < p; k++) {     // (B'(x - mu))_k
    double z = 0;
    for (int i = k; i < p; i++) z += c.B[i][k] * (x[i] - c.mu[i]);
    q += z * z;
  }
  return c.cst - 0.5 * q;
}

// draw (mu, Sigma) from the NIW posterior given m points with sum s and sum of outer products Q
void draw_niw(Par &c, const Prior &pr, int m, const double *s, const double Q[3][3]) {
  const int p = pr.p;
  const double tm = pr.tau0 + m, nm = pr.nu0 + m;
  double mm[3], Sm[3][3];
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Sm[i][j] = (i < p && j < p) ? pr.S0[i][j] : 0;
  for (int i = 0; i < p; i++) mm[i] = (pr.tau0 * pr.mu0[i] + s[i]) / tm;
  if (m > 0) {
    double xb[3], d[3], f = pr.tau0 * m / tm;
    for (int i = 0; i < p; i++) { xb[i] = s[i] / m; d[i] = xb[i] - pr.mu0[i]; }
    for (int i = 0; i < p; i++) for (int j = 0; j < p; j++)
      Sm[i][j] += Q[i][j] - m * xb[i] * xb[j] + f * d[i] * d[j];
  }
  // Sigma^{-1} ~ Wishart(nm, Sm^{-1}): C = chol(Sm^{-1}), Bartlett A, B = C A
  double L[3][3], Li[3][3], Sinv[3][3], C[3][3], A[3][3];
  chol(Sm, L, p); inv_lower(L, Li, p);
  for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) {   // Sm^{-1} = Li' Li
    double v = 0;
    for (int k = std::max(i, j); k < p; k++) v += Li[k][i] * Li[k][j];
    Sinv[i][j] = v;
  }
  chol(Sinv, C, p);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) A[i][j] = 0;
  for (int j = 0; j < p; j++) {
    A[j][j] = std::sqrt(R::rchisq(nm - j));   // correct Bartlett: chi^2(nu - j), j = 0..p-1
    for (int i = j + 1; i < p; i++) A[i][j] = norm_rand();
  }
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
    double v = 0;
    if (i < p && j < p) for (int k = j; k <= i; k++) v += C[i][k] * A[k][j];
    c.B[i][j] = v;
  }
  double Bi[3][3];
  inv_lower(c.B, Bi, p);
  // Sigma = Bi' Bi ; mu = mm + Bi' z / sqrt(tm)
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
    double v = 0;
    if (i < p && j < p) for (int k = std::max(i, j); k < p; k++) v += Bi[k][i] * Bi[k][j];
    c.Sig[i][j] = v;
  }
  double z[3];
  for (int k = 0; k < p; k++) z[k] = norm_rand() / std::sqrt(tm);
  for (int i = 0; i < p; i++) {
    double v = 0;
    for (int k = i; k < p; k++) v += Bi[k][i] * z[k];
    c.mu[i] = mm[i] + v;
  }
  double ld = 0;
  for (int k = 0; k < p; k++) ld += std::log(c.B[k][k]);
  c.cst = -0.5 * p * std::log(2 * M_PI) + ld;
}

inline double lline(double w, double x, double y, double lx, int p, const Par &c) {
  double w2 = (y - x * w) / (1.0 - x);
  if (!(w > 0 && w < 1 && w2 > 0 && w2 < 1)) return R_NegInf;
  double s[3] = {logitc(w), logitc(w2), lx};
  return lnorm(s, c, p) - std::log(w) - std::log1p(-w) - std::log(w2) - std::log1p(-w2);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List eco_alg1_cpp(int p, Rcpp::NumericVector X, Rcpp::NumericVector Y,
                        Rcpp::NumericVector L, Rcpp::NumericVector U,
                        Rcpp::NumericMatrix Wobs, Rcpp::NumericVector w1obs,
                        Rcpp::NumericVector w2obs,
                        Rcpp::NumericVector mu0, double tau0, double nu0,
                        Rcpp::NumericMatrix S0, double a0, double b0,
                        int n_iter, int burnin, int thin,
                        Rcpp::NumericVector W1init, int wstep, int grouped,
                        double init_alpha) {
  Rcpp::RNGScope scope;
  Prior pr; pr.p = p; pr.tau0 = tau0; pr.nu0 = nu0;
  for (int i = 0; i < 3; i++) {
    pr.mu0[i] = i < p ? mu0[i] : 0;
    for (int j = 0; j < 3; j++) pr.S0[i][j] = (i < p && j < p) ? S0(i, j) : 0;
  }
  const T0 t0 = make_t0(pr);
  const int n = X.size(), no = Wobs.ncol(), n1 = w1obs.size(), n0 = w2obs.size();
  const int nt = n + no + n1 + n0;
  std::vector<double> W1(n), lx(n);
  std::vector<std::array<double, 3>> Ws(nt);
  for (int i = 0; i < n; i++) {
    W1[i] = W1init[i]; lx[i] = logitc(X[i]);
    double w2 = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
    Ws[i] = {logitc(W1[i]), logitc(w2), p == 3 ? lx[i] : 0.0};
  }
  for (int j = 0; j < no; j++) Ws[n + j] = {Wobs(0, j), Wobs(1, j), p == 3 ? Wobs(2, j) : 0.0};
  for (int j = 0; j < n1; j++) Ws[n + no + j] = {w1obs[j], 0.0, 0.0};
  for (int j = 0; j < n0; j++) Ws[n + no + n1 + j] = {0.0, w2obs[j], 0.0};

  // start: one cluster drawn from the prior (all units share it)
  std::vector<Par> cl(1);
  std::vector<int> m(1, nt), z(nt, 0);
  { double s[3] = {0, 0, 0}, Q[3][3] = {{0}}; draw_niw(cl[0], pr, 0, s, Q); }
  double alpha = init_alpha;

  int n_save = (n_iter - burnin) / thin;
  Rcpp::NumericMatrix sW1(n_save, n), sW2(n_save, n);
  Rcpp::NumericVector salpha(n_save);
  Rcpp::IntegerVector sk(n_save), smax(n_save);
  int s_ = 0;
  std::vector<double> lp(nt + 1);

  for (int it = 0; it < n_iter; it++) {
    // 1. latent rates given the unit's current parameters
    for (int i = 0; i < n; i++) {
      const Par &c = cl[z[i]];
      double w = W1[i];
      if (wstep == 0) {            // eco's rMH: uniform independence proposal
        double prop = L[i] + unif_rand() * (U[i] - L[i]);
        double r = lline(prop, X[i], Y[i], lx[i], p, c) - lline(w, X[i], Y[i], lx[i], p, c);
        if (std::log(unif_rand()) < r) w = prop;
      } else {                     // slice sampler (shrinkage on [L, U])
        double h = lline(w, X[i], Y[i], lx[i], p, c) - exp_rand();
        double a = L[i], b = U[i];
        for (int k = 0; k < 200; k++) {
          double prop = a + unif_rand() * (b - a);
          if (lline(prop, X[i], Y[i], lx[i], p, c) > h) { w = prop; break; }
          if (prop < w) a = prop; else b = prop;
        }
      }
      W1[i] = w;
      double w2 = (Y[i] - X[i] * w) / (1.0 - X[i]);
      Ws[i][0] = logitc(w); Ws[i][1] = logitc(w2);
    }
    // 2. partially observed units: conditional normal
    for (int j = 0; j < n1 + n0; j++) {
      int i = n + no + j, o = j < n1 ? 0 : 1, mi = 1 - o;
      const Par &c = cl[z[i]];
      double mean = c.mu[mi] + c.Sig[mi][o] / c.Sig[o][o] * (Ws[i][o] - c.mu[o]);
      double var = c.Sig[mi][mi] - c.Sig[mi][o] * c.Sig[mi][o] / c.Sig[o][o];
      Ws[i][mi] = mean + std::sqrt(var) * norm_rand();
    }
    // 3. configuration (Algorithm 1)
    for (int i = 0; i < nt; i++) {
      const double *x = Ws[i].data();
      int c = z[i];
      m[c]--;
      if (m[c] == 0) {               // i held these parameters alone: discard them
        int last = cl.size() - 1;
        if (c != last) {
          cl[c] = cl[last]; m[c] = m[last];
          for (int j = 0; j < nt; j++) if (z[j] == last) z[j] = c;
        }
        cl.pop_back(); m.pop_back();
      }
      z[i] = -1;
      int K = cl.size(), pick;
      double lnew = std::log(alpha) + lt0(x, t0, p);
      if (grouped) {
        for (int k = 0; k < K; k++) lp[k] = std::log((double) m[k]) + lnorm(x, cl[k], p);
        lp[K] = lnew;
        double mx = *std::max_element(lp.begin(), lp.begin() + K + 1), tot = 0;
        for (int k = 0; k <= K; k++) { lp[k] = std::exp(lp[k] - mx); tot += lp[k]; }
        double u = unif_rand() * tot, acc = 0; pick = K;
        for (int k = 0; k <= K; k++) { acc += lp[k]; if (u <= acc) { pick = k; break; } }
      } else {                        // eco's loop over units j != i
        double mx = lnew;
        for (int j = 0; j < nt; j++) {
          if (j == i) { lp[j] = R_NegInf; continue; }
          lp[j] = lnorm(x, cl[z[j]], p);
          if (lp[j] > mx) mx = lp[j];
        }
        double tot = std::exp(lnew - mx);
        double pnew = tot;
        for (int j = 0; j < nt; j++) { lp[j] = (j == i) ? 0 : std::exp(lp[j] - mx); tot += lp[j]; }
        double u = unif_rand() * tot; pick = K;
        if (u > pnew) {
          double acc = pnew;
          for (int j = 0; j < nt; j++) { acc += lp[j]; if (u <= acc && j != i) { pick = z[j]; break; } }
        }
      }
      if (pick == K) {                // new parameters from the posterior given x alone
        Par nc;
        double s[3] = {x[0], x[1], x[2]}, Q[3][3];
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) Q[a][b] = x[a] * x[b];
        draw_niw(nc, pr, 1, s, Q);
        cl.push_back(nc); m.push_back(0);
      }
      z[i] = pick; m[pick]++;
    }
    // 4. remixing: each cluster's parameters given its members
    const int K = cl.size();
    std::vector<std::array<double, 3>> S(K, {0, 0, 0});
    std::vector<std::array<double, 9>> QQ(K);
    for (auto &q : QQ) q.fill(0);
    for (int i = 0; i < nt; i++) {
      int c = z[i];
      for (int a = 0; a < p; a++) {
        S[c][a] += Ws[i][a];
        for (int b = 0; b < p; b++) QQ[c][3 * a + b] += Ws[i][a] * Ws[i][b];
      }
    }
    for (int c = 0; c < K; c++) {
      double Q[3][3];
      for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) Q[a][b] = QQ[c][3 * a + b];
      draw_niw(cl[c], pr, m[c], S[c].data(), Q);
    }
    // 5. alpha (Escobar & West 1995), mixture weight as a probability
    double eta = R::rbeta(alpha + 1.0, (double) nt);
    double rate = b0 - std::log(eta);
    double odds = (a0 + K - 1.0) / (nt * rate);
    double pmix = odds / (1.0 + odds);
    alpha = (unif_rand() < pmix) ? R::rgamma(a0 + K, 1.0 / rate) : R::rgamma(a0 + K - 1.0, 1.0 / rate);

    if (it >= burnin && (it - burnin) % thin == 0 && s_ < n_save) {
      for (int i = 0; i < n; i++) {
        sW1(s_, i) = W1[i];
        sW2(s_, i) = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
      }
      salpha[s_] = alpha; sk[s_] = K;
      smax[s_] = *std::max_element(m.begin(), m.end());
      s_++;
    }
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("W1") = sW1, Rcpp::Named("W2") = sW2,
                            Rcpp::Named("alpha") = salpha, Rcpp::Named("n_clusters") = sk,
                            Rcpp::Named("largest_cluster") = smax);
}

// [[Rcpp::export]]
Rcpp::List eco_param_cpp(int p, Rcpp::NumericVector X, Rcpp::NumericVector Y,
                         Rcpp::NumericVector L, Rcpp::NumericVector U,
                         Rcpp::NumericMatrix Wobs, Rcpp::NumericVector w1obs,
                         Rcpp::NumericVector w2obs,
                         Rcpp::NumericVector mu0, double tau0, double nu0,
                         Rcpp::NumericMatrix S0, int n_iter, int burnin, int thin,
                         Rcpp::NumericVector W1init, int wstep) {
  Rcpp::RNGScope scope;
  Prior pr; pr.p = p; pr.tau0 = tau0; pr.nu0 = nu0;
  for (int i = 0; i < 3; i++) {
    pr.mu0[i] = i < p ? mu0[i] : 0;
    for (int j = 0; j < 3; j++) pr.S0[i][j] = (i < p && j < p) ? S0(i, j) : 0;
  }
  const int n = X.size(), no = Wobs.ncol(), n1 = w1obs.size(), n0 = w2obs.size();
  const int nt = n + no + n1 + n0;
  std::vector<double> W1(n), lx(n);
  std::vector<std::array<double, 3>> Ws(nt);
  for (int i = 0; i < n; i++) {
    W1[i] = W1init[i]; lx[i] = logitc(X[i]);
    double w2 = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
    Ws[i] = {logitc(W1[i]), logitc(w2), p == 3 ? lx[i] : 0.0};
  }
  for (int j = 0; j < no; j++) Ws[n + j] = {Wobs(0, j), Wobs(1, j), p == 3 ? Wobs(2, j) : 0.0};
  for (int j = 0; j < n1; j++) Ws[n + no + j] = {w1obs[j], w1obs[j], 0.0};   // missing coordinate started
  for (int j = 0; j < n0; j++) Ws[n + no + n1 + j] = {w2obs[j], w2obs[j], 0.0}; // at the observed one

  auto update_par = [&](Par &c) {
    double s[3] = {0, 0, 0}, Q[3][3] = {{0}};
    for (int i = 0; i < nt; i++)
      for (int a = 0; a < p; a++) {
        s[a] += Ws[i][a];
        for (int b = 0; b < p; b++) Q[a][b] += Ws[i][a] * Ws[i][b];
      }
    draw_niw(c, pr, nt, s, Q);
  };
  Par c;
  update_par(c);

  int n_save = (n_iter - burnin) / thin;
  Rcpp::NumericMatrix sW1(n_save, n), sW2(n_save, n), smu(n_save, p), sSig(n_save, p * p);
  int s_ = 0;
  for (int it = 0; it < n_iter; it++) {
    for (int i = 0; i < n; i++) {
      double w = W1[i];
      if (wstep == 0) {
        double prop = L[i] + unif_rand() * (U[i] - L[i]);
        double r = lline(prop, X[i], Y[i], lx[i], p, c) - lline(w, X[i], Y[i], lx[i], p, c);
        if (std::log(unif_rand()) < r) w = prop;
      } else {
        double h = lline(w, X[i], Y[i], lx[i], p, c) - exp_rand();
        double a = L[i], b = U[i];
        for (int k = 0; k < 200; k++) {
          double prop = a + unif_rand() * (b - a);
          if (lline(prop, X[i], Y[i], lx[i], p, c) > h) { w = prop; break; }
          if (prop < w) a = prop; else b = prop;
        }
      }
      W1[i] = w;
      double w2 = (Y[i] - X[i] * w) / (1.0 - X[i]);
      Ws[i][0] = logitc(w); Ws[i][1] = logitc(w2);
    }
    for (int j = 0; j < n1 + n0; j++) {
      int i = n + no + j, o = j < n1 ? 0 : 1, mi = 1 - o;
      double mean = c.mu[mi] + c.Sig[mi][o] / c.Sig[o][o] * (Ws[i][o] - c.mu[o]);
      double var = c.Sig[mi][mi] - c.Sig[mi][o] * c.Sig[mi][o] / c.Sig[o][o];
      Ws[i][mi] = mean + std::sqrt(var) * norm_rand();
    }
    update_par(c);
    if (it >= burnin && (it - burnin) % thin == 0 && s_ < n_save) {
      for (int i = 0; i < n; i++) {
        sW1(s_, i) = W1[i];
        sW2(s_, i) = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
      }
      for (int a = 0; a < p; a++) {
        smu(s_, a) = c.mu[a];
        for (int b = 0; b < p; b++) sSig(s_, a * p + b) = c.Sig[a][b];
      }
      s_++;
    }
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("W1") = sW1, Rcpp::Named("W2") = sW2,
                            Rcpp::Named("mu") = smu, Rcpp::Named("Sigma") = sSig);
}
