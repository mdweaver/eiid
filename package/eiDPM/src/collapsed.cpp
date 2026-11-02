// Collapsed Gibbs sampler for the (untruncated) Dirichlet-process mixture of
// normals in 2x2 ecological inference (Imai, Lu & Strauss 2008):
// Neal (2000), Algorithm 3, with the conjugate normal-inverse-Wishart base
// measure integrated out, and the Escobar & West (1995) update for alpha.
//
// Base measure (eco's code conventions): mu | Sigma ~ N(mu0, Sigma / tau0),
// Sigma ~ InvWishart(nu0, S0) (scale S0). A cluster holding m points has
// posterior (tau_m, mu_m, nu_m, S_m) and posterior predictive
//   t_{nu_m - p + 1}( mu_m, S_m (tau_m + 1) / (tau_m (nu_m - p + 1)) ).
//
// One sweep, for each unit i: remove i from its cluster; for a latent unit,
// update W_i on its tomography line given the predictive of its cluster
// without i (slice sampling); then draw the cluster of i given W_i from
//   P(existing c) ~ n_c t_c(W_i*),  P(new) ~ alpha t_0(W_i*);
// add i back. Then update alpha. Labels carry no order, so there is no
// truncation and no label-ordering problem.
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

inline double logitc(double w) { return std::log(w) - std::log1p(-w); }

struct Clust {
  int m;
  double s[3];     // sum of points
  double Q[3][3];  // sum of outer products
  // predictive t cache
  double df, cst;
  double loc[3];
  double Li[3][3];  // inverse of chol(scale), lower
};

struct Prior {
  int p;
  double mu0[3];
  double tau0, nu0;
  double S0[3][3];
  std::vector<double> lgam_ratio;  // lgamma((df+p)/2) - lgamma(df/2) - p/2 log(df pi), by m
};

void cache(Clust &c, const Prior &pr) {
  const int p = pr.p;
  const double tm = pr.tau0 + c.m;
  double mm[3], Sm[3][3];
  for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) Sm[i][j] = pr.S0[i][j];
  if (c.m > 0) {
    double xb[3], d[3];
    for (int i = 0; i < p; i++) {
      xb[i] = c.s[i] / c.m;
      mm[i] = (pr.tau0 * pr.mu0[i] + c.s[i]) / tm;
      d[i] = xb[i] - pr.mu0[i];
    }
    double f = pr.tau0 * c.m / tm;
    for (int i = 0; i < p; i++) for (int j = 0; j < p; j++)
      Sm[i][j] += c.Q[i][j] - c.m * xb[i] * xb[j] + f * d[i] * d[j];
  } else {
    for (int i = 0; i < p; i++) mm[i] = pr.mu0[i];
  }
  c.df = pr.nu0 + c.m - p + 1;
  const double sc = (tm + 1) / (tm * c.df);
  // Cholesky of sc * Sm (p <= 3)
  double L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  double ldet = 0;
  for (int j = 0; j < p; j++) {
    double v = sc * Sm[j][j];
    for (int k = 0; k < j; k++) v -= L[j][k] * L[j][k];
    L[j][j] = std::sqrt(v);
    ldet += 2.0 * std::log(L[j][j]);
    for (int i = j + 1; i < p; i++) {
      double w = sc * Sm[i][j];
      for (int k = 0; k < j; k++) w -= L[i][k] * L[j][k];
      L[i][j] = w / L[j][j];
    }
  }
  // inverse of the lower-triangular factor
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) c.Li[i][j] = 0;
  for (int j = 0; j < p; j++) {
    c.Li[j][j] = 1.0 / L[j][j];
    for (int i = j + 1; i < p; i++) {
      double w = 0;
      for (int k = j; k < i; k++) w -= L[i][k] * c.Li[k][j];
      c.Li[i][j] = w / L[i][i];
    }
  }
  for (int i = 0; i < p; i++) c.loc[i] = mm[i];
  c.cst = pr.lgam_ratio[c.m] - 0.5 * ldet;
}

inline double lt(const double *x, const Clust &c, int p) {
  double q = 0;
  for (int i = 0; i < p; i++) {
    double zi = 0;
    for (int j = 0; j <= i; j++) zi += c.Li[i][j] * (x[j] - c.loc[j]);
    q += zi * zi;
  }
  return c.cst - 0.5 * (c.df + p) * std::log1p(q / c.df);
}

void add_pt(Clust &c, const double *x, int p, int sign) {
  c.m += sign;
  for (int i = 0; i < p; i++) {
    c.s[i] += sign * x[i];
    for (int j = 0; j < p; j++) c.Q[i][j] += sign * x[i] * x[j];
  }
}

Clust empty_clust(const Prior &pr) {
  Clust c;
  c.m = 0;
  for (int i = 0; i < 3; i++) { c.s[i] = 0; for (int j = 0; j < 3; j++) c.Q[i][j] = 0; }
  cache(c, pr);
  return c;
}

inline double lline(double w, double x, double y, double lx, int p, const Clust &c) {
  double w2 = (y - x * w) / (1.0 - x);
  if (!(w > 0 && w < 1 && w2 > 0 && w2 < 1)) return R_NegInf;
  double s[3] = {logitc(w), logitc(w2), lx};
  return lt(s, c, p) - std::log(w) - std::log1p(-w) - std::log(w2) - std::log1p(-w2);
}

double slice1(double w, double lo, double hi, double x, double y, double lx,
              int p, const Clust &c) {
  double f0 = lline(w, x, y, lx, p, c);
  double h = f0 - exp_rand();
  double a = lo, b = hi;
  for (int it = 0; it < 200; it++) {
    double prop = a + unif_rand() * (b - a);
    if (lline(prop, x, y, lx, p, c) > h) return prop;
    if (prop < w) a = prop; else b = prop;
  }
  return w;
}

// draw the missing coordinate of a bivariate t given the observed one
double rcond_t(const Clust &c, int obs, double xo) {
  // rebuild the 2x2 scale from its inverse Cholesky factor
  MatrixXd Li(2, 2);
  Li << c.Li[0][0], 0, c.Li[1][0], c.Li[1][1];
  MatrixXd L = Li.triangularView<Eigen::Lower>().solve(MatrixXd::Identity(2, 2));
  MatrixXd S = L * L.transpose();
  int mis = 1 - obs;
  double d = xo - c.loc[obs];
  double m = c.loc[mis] + S(mis, obs) / S(obs, obs) * d;
  double s2 = S(mis, mis) - S(mis, obs) * S(mis, obs) / S(obs, obs);
  double dfc = c.df + 1;
  double sc2 = (c.df + d * d / S(obs, obs)) / dfc * s2;
  return m + std::sqrt(sc2) * norm_rand() / std::sqrt(R::rchisq(dfc) / dfc);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List dpm_collapsed_cpp(int p, Rcpp::NumericVector X, Rcpp::NumericVector Y,
                             Rcpp::NumericVector L, Rcpp::NumericVector U,
                             Rcpp::NumericMatrix Wobs, Rcpp::NumericVector w1obs,
                             Rcpp::NumericVector w2obs,
                             Rcpp::NumericVector mu0, double tau0, double nu0,
                             Rcpp::NumericMatrix S0, double a0, double b0,
                             int n_iter, int burnin, int thin,
                             Rcpp::NumericVector W1init, int slice_steps,
                             double init_alpha) {
  Rcpp::RNGScope scope;
  Prior pr;
  pr.p = p; pr.tau0 = tau0; pr.nu0 = nu0;
  for (int i = 0; i < 3; i++) {
    pr.mu0[i] = i < p ? mu0[i] : 0;
    for (int j = 0; j < 3; j++) pr.S0[i][j] = (i < p && j < p) ? S0(i, j) : 0;
  }
  const int n = X.size(), no = Wobs.ncol(), n1 = w1obs.size(), n0 = w2obs.size();
  const int nt = n + no + n1 + n0;
  pr.lgam_ratio.resize(nt + 2);
  for (int m = 0; m <= nt + 1; m++) {
    double df = nu0 + m - p + 1;
    pr.lgam_ratio[m] = std::lgamma(0.5 * (df + p)) - std::lgamma(0.5 * df) - 0.5 * p * std::log(df * M_PI);
  }
  // logit-scale points (columns)
  MatrixXd Ws = MatrixXd::Zero(p, nt);
  std::vector<double> W1(n), lx(n);
  for (int i = 0; i < n; i++) {
    W1[i] = W1init[i]; lx[i] = logitc(X[i]);
    double w2 = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
    Ws(0, i) = logitc(W1[i]); Ws(1, i) = logitc(w2);
    if (p == 3) Ws(2, i) = lx[i];
  }
  for (int j = 0; j < no; j++) for (int r = 0; r < p; r++) Ws(r, n + j) = Wobs(r, j);
  for (int j = 0; j < n1; j++) { Ws(0, n + no + j) = w1obs[j]; Ws(1, n + no + j) = 0; }
  for (int j = 0; j < n0; j++) { Ws(0, n + no + n1 + j) = 0; Ws(1, n + no + n1 + j) = w2obs[j]; }

  Clust prior_c = empty_clust(pr);
  // start: every unit in one cluster
  std::vector<Clust> cl;
  cl.push_back(empty_clust(pr));
  std::vector<int> z(nt, 0);
  // partially observed units: impute from the prior predictive first
  for (int j = 0; j < n1 + n0; j++) {
    int col = n + no + j, obs = j < n1 ? 0 : 1;
    Ws(1 - obs, col) = rcond_t(prior_c, obs, Ws(obs, col));
  }
  for (int i = 0; i < nt; i++) { double x[3] = {Ws(0, i), Ws(1, i), p == 3 ? Ws(2, i) : 0}; add_pt(cl[0], x, p, 1); }
  cache(cl[0], pr);
  double alpha = init_alpha;

  int n_save = (n_iter - burnin) / thin;
  Rcpp::NumericMatrix sW1(n_save, n), sW2(n_save, n);
  Rcpp::NumericVector salpha(n_save);
  Rcpp::IntegerVector sk(n_save), smax(n_save);
  int s = 0;
  std::vector<double> lp;

  for (int it = 0; it < n_iter; it++) {
    for (int i = 0; i < nt; i++) {
      int c = z[i];
      double x[3] = {Ws(0, i), Ws(1, i), p == 3 ? Ws(2, i) : 0};
      add_pt(cl[c], x, p, -1);
      if (cl[c].m == 0) {                    // drop the empty cluster
        int last = cl.size() - 1;
        if (c != last) {
          cl[c] = cl[last];
          for (int j = 0; j < nt; j++) if (z[j] == last) z[j] = c;
        }
        cl.pop_back();
        c = -1;
      } else {
        cache(cl[c], pr);
      }
      // latent unit: update W_i given the predictive of its (former) cluster
      if (i < n) {
        const Clust &cc = c >= 0 ? cl[c] : prior_c;
        for (int r = 0; r < slice_steps; r++)
          W1[i] = slice1(W1[i], L[i], U[i], X[i], Y[i], lx[i], p, cc);
        double w2 = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
        x[0] = logitc(W1[i]); x[1] = logitc(w2);
        Ws(0, i) = x[0]; Ws(1, i) = x[1];
      } else if (i >= n + no) {               // partially observed (CAR only)
        int obs = (i < n + no + n1) ? 0 : 1;
        const Clust &cc = c >= 0 ? cl[c] : prior_c;
        x[1 - obs] = rcond_t(cc, obs, x[obs]);
        Ws(1 - obs, i) = x[1 - obs];
      }
      // reassign
      int K = cl.size();
      lp.assign(K + 1, 0.0);
      for (int k = 0; k < K; k++) lp[k] = std::log((double) cl[k].m) + lt(x, cl[k], p);
      lp[K] = std::log(alpha) + lt(x, prior_c, p);
      double mx = *std::max_element(lp.begin(), lp.end()), tot = 0;
      for (int k = 0; k <= K; k++) { lp[k] = std::exp(lp[k] - mx); tot += lp[k]; }
      double u = unif_rand() * tot, acc = 0;
      int kk = K;
      for (int k = 0; k <= K; k++) { acc += lp[k]; if (u <= acc) { kk = k; break; } }
      if (kk == K) cl.push_back(empty_clust(pr));
      add_pt(cl[kk], x, p, 1);
      cache(cl[kk], pr);
      z[i] = kk;
    }
    // alpha: Escobar & West (1995), auxiliary eta
    int k = cl.size();
    double eta = R::rbeta(alpha + 1.0, (double) nt);
    double rate = b0 - std::log(eta);
    double odds = (a0 + k - 1.0) / (nt * rate);
    double pmix = odds / (1.0 + odds);
    alpha = (unif_rand() < pmix) ? R::rgamma(a0 + k, 1.0 / rate) : R::rgamma(a0 + k - 1.0, 1.0 / rate);

    if (it >= burnin && (it - burnin) % thin == 0 && s < n_save) {
      for (int i = 0; i < n; i++) {
        sW1(s, i) = W1[i];
        sW2(s, i) = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
      }
      salpha[s] = alpha; sk[s] = k;
      int mx = 0;
      for (int j = 0; j < k; j++) mx = std::max(mx, cl[j].m);
      smax[s] = mx;
      s++;
    }
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("W1") = sW1, Rcpp::Named("W2") = sW2,
                            Rcpp::Named("alpha") = salpha, Rcpp::Named("n_clusters") = sk,
                            Rcpp::Named("largest_cluster") = smax);
}
