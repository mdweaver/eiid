// Blocked Gibbs sampler (Ishwaran & James 2001) for the truncated
// Dirichlet-process mixture of normals in 2x2 ecological inference
// (Imai, Lu & Strauss 2008). p = 2 (CAR) or p = 3 (NCAR, third row logit X).
// K = 1 gives the parametric model.
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

const double LOG2PI = std::log(2.0 * M_PI);

inline double logit(double w) { return std::log(w) - std::log1p(-w); }

struct Comp {
  VectorXd mu;
  MatrixXd Sigma;
  double Li[3][3];  // inverse of chol(Sigma), lower triangular
  double m[3];
  double cst;       // -0.5 (p log 2 pi + log det Sigma)
  int p;
};

void refresh(Comp &c) {
  c.p = c.mu.size();
  Eigen::LLT<MatrixXd> llt(c.Sigma);
  MatrixXd L = llt.matrixL();
  MatrixXd Li = L.triangularView<Eigen::Lower>().solve(MatrixXd::Identity(c.p, c.p));
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    c.Li[i][j] = (i < c.p && j < c.p) ? Li(i, j) : 0.0;
  for (int i = 0; i < c.p; i++) c.m[i] = c.mu(i);
  double ldet = 2.0 * L.diagonal().array().log().sum();
  c.cst = -0.5 * (c.p * LOG2PI + ldet);
}

inline double lnormv(const double *x, const Comp &c) {
  double q = 0;
  for (int i = 0; i < c.p; i++) {
    double zi = 0;
    for (int j = 0; j <= i; j++) zi += c.Li[i][j] * (x[j] - c.m[j]);
    q += zi * zi;
  }
  return c.cst - 0.5 * q;
}

// Sigma ~ InvWishart(nu, S): Sigma^{-1} ~ Wishart(nu, S^{-1}), Bartlett.
MatrixXd rinvwish(double nu, const MatrixXd &S) {
  int p = S.rows();
  MatrixXd Sinv = S.inverse();
  Eigen::LLT<MatrixXd> llt(Sinv);
  MatrixXd L = llt.matrixL();
  MatrixXd A = MatrixXd::Zero(p, p);
  for (int j = 0; j < p; j++) {
    A(j, j) = std::sqrt(R::rchisq(nu - j));
    for (int i = j + 1; i < p; i++) A(i, j) = norm_rand();
  }
  MatrixXd LA = L * A;
  MatrixXd Prec = LA * LA.transpose();
  return Prec.inverse();
}

VectorXd rmvn(const VectorXd &mu, const MatrixXd &Sigma) {
  Eigen::LLT<MatrixXd> llt(Sigma);
  MatrixXd L = llt.matrixL();
  VectorXd z(mu.size());
  for (int j = 0; j < mu.size(); j++) z(j) = norm_rand();
  return mu + L * z;
}

// log density of W1 = w on the tomography line (up to a constant)
inline double lline(double w, double x, double y, double lx, int p, const Comp &c) {
  double w2 = (y - x * w) / (1.0 - x);
  if (!(w > 0 && w < 1 && w2 > 0 && w2 < 1)) return R_NegInf;
  double s[3] = {logit(w), logit(w2), lx};
  return lnormv(s, c) - std::log(w) - std::log1p(-w) - std::log(w2) - std::log1p(-w2);
}

// one slice-sampling update on (lo, hi) with shrinkage from the full interval
double slice1(double w, double lo, double hi, double x, double y, double lx,
              int p, const Comp &c) {
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

}  // namespace

// [[Rcpp::export]]
Rcpp::List dpm_gibbs_cpp(int p, int K,
                         Rcpp::NumericVector X, Rcpp::NumericVector Y,
                         Rcpp::NumericVector L, Rcpp::NumericVector U,
                         Rcpp::NumericMatrix Wobs,   // p x n_obs, logit scale
                         Rcpp::NumericVector w1obs,  // X = 1 units (CAR)
                         Rcpp::NumericVector w2obs,  // X = 0 units (CAR)
                         Rcpp::NumericVector mu0, double tau0, double nu0,
                         Rcpp::NumericMatrix S0, double a0, double b0,
                         int n_iter, int burnin, int thin,
                         Rcpp::NumericVector W1init, int slice_steps,
                         int init_k = 1, double init_alpha = 1.0) {
  Rcpp::RNGScope scope;
  const int n = X.size(), no = Wobs.ncol(), n1 = w1obs.size(), n0 = w2obs.size();
  const int nt = n + no + n1 + n0;
  VectorXd m0 = Rcpp::as<VectorXd>(mu0);
  MatrixXd S0m = Rcpp::as<MatrixXd>(S0);

  // logit-scale data matrix for all units (columns), filled each sweep
  MatrixXd Ws(p, nt);
  std::vector<double> W1(n), lx(n);
  for (int i = 0; i < n; i++) {
    W1[i] = W1init[i];
    lx[i] = logit(X[i]);
  }
  for (int j = 0; j < no; j++) for (int r = 0; r < p; r++) Ws(r, n + j) = Wobs(r, j);
  for (int j = 0; j < n1; j++) { Ws(0, n + no + j) = w1obs[j]; Ws(1, n + no + j) = 0; }
  for (int j = 0; j < n0; j++) { Ws(0, n + no + n1 + j) = 0; Ws(1, n + no + n1 + j) = w2obs[j]; }

  // initial components: prior draws; all units start in component 1 with
  // alpha = 1 (starting at large alpha with many occupied components can leave
  // the truncated sampler in a metastable large-alpha region where the last
  // component absorbs most units)
  std::vector<Comp> comp(K);
  for (int k = 0; k < K; k++) {
    comp[k].Sigma = rinvwish(nu0, S0m);
    comp[k].mu = rmvn(m0, comp[k].Sigma / tau0);
    refresh(comp[k]);
  }
  std::vector<int> z(nt);
  for (int i = 0; i < nt; i++) z[i] = init_k <= 1 ? 0 : (int) std::floor(unif_rand() * std::min(K, init_k));
  VectorXd logpi(K);
  {
    double acc = 0;
    for (int k = 0; k < K - 1; k++) { logpi(k) = std::log(0.5) + acc; acc += std::log(0.5); }
    logpi(K - 1) = acc;
  }
  std::vector<double> v(K, 0.5);
  double alpha = init_alpha;

  int n_save = (n_iter - burnin) / thin;
  Rcpp::NumericMatrix sW1(n_save, n), sW2(n_save, n);
  Rcpp::NumericVector salpha(n_save), slast(n_save);
  Rcpp::IntegerVector socc(n_save), slastn(n_save);
  int s = 0;

  for (int it = 0; it < n_iter; it++) {
    // 1. latent rates given components
    for (int i = 0; i < n; i++) {
      const Comp &c = comp[z[i]];
      for (int r = 0; r < slice_steps; r++)
        W1[i] = slice1(W1[i], L[i], U[i], X[i], Y[i], lx[i], p, c);
      double w2 = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
      Ws(0, i) = logit(W1[i]); Ws(1, i) = logit(w2);
      if (p == 3) Ws(2, i) = lx[i];
    }
    // partially observed units (CAR): impute the missing coordinate
    for (int j = 0; j < n1 + n0; j++) {
      int col = n + no + j;
      const Comp &c = comp[z[col]];
      int obs = j < n1 ? 0 : 1, mis = 1 - obs;
      double m = c.mu(mis) + c.Sigma(mis, obs) / c.Sigma(obs, obs) * (Ws(obs, col) - c.mu(obs));
      double sd = std::sqrt(c.Sigma(mis, mis) - c.Sigma(mis, obs) * c.Sigma(mis, obs) / c.Sigma(obs, obs));
      Ws(mis, col) = m + sd * norm_rand();
    }
    // 2. assignments
    std::vector<int> cnt(K, 0);
    VectorXd lp(K);
    for (int i = 0; i < nt; i++) {
      double x[3] = {Ws(0, i), Ws(1, i), p == 3 ? Ws(2, i) : 0.0};
      for (int k = 0; k < K; k++) lp(k) = logpi(k) + lnormv(x, comp[k]);
      double mx = lp.maxCoeff(), tot = 0;
      for (int k = 0; k < K; k++) { lp(k) = std::exp(lp(k) - mx); tot += lp(k); }
      double u = unif_rand() * tot, acc = 0;
      int kk = K - 1;
      for (int k = 0; k < K; k++) { acc += lp(k); if (u <= acc) { kk = k; break; } }
      z[i] = kk; cnt[kk]++;
    }
    // 3. stick-breaking weights
    if (K > 1) {
      int above = nt;
      double acc = 0, sumlog = 0;
      for (int k = 0; k < K - 1; k++) {
        above -= cnt[k];
        v[k] = R::rbeta(1.0 + cnt[k], alpha + above);
        v[k] = std::min(std::max(v[k], 1e-300), 1.0 - 1e-16);
        logpi(k) = std::log(v[k]) + acc;
        acc += std::log1p(-v[k]);
        sumlog += std::log1p(-v[k]);
      }
      logpi(K - 1) = acc;
      // 5. concentration
      alpha = R::rgamma(a0 + K - 1, 1.0 / (b0 - sumlog));
    } else {
      logpi(0) = 0;
    }
    // 4. component parameters (conjugate NIW)
    for (int k = 0; k < K; k++) {
      VectorXd xbar = VectorXd::Zero(p);
      MatrixXd Sc = MatrixXd::Zero(p, p);
      if (cnt[k] > 0) {
        for (int i = 0; i < nt; i++) if (z[i] == k) xbar += Ws.col(i);
        xbar /= cnt[k];
        for (int i = 0; i < nt; i++) if (z[i] == k) {
          VectorXd d = Ws.col(i) - xbar; Sc += d * d.transpose();
        }
      }
      double tn = tau0 + cnt[k];
      VectorXd mn = (tau0 * m0 + cnt[k] * xbar) / tn;
      VectorXd dm = xbar - m0;
      MatrixXd Sn = S0m + Sc + (tau0 * cnt[k] / tn) * dm * dm.transpose();
      comp[k].Sigma = rinvwish(nu0 + cnt[k], Sn);
      comp[k].mu = rmvn(mn, comp[k].Sigma / tn);
      refresh(comp[k]);
    }
    if (it >= burnin && (it - burnin) % thin == 0 && s < n_save) {
      for (int i = 0; i < n; i++) {
        sW1(s, i) = W1[i];
        sW2(s, i) = (Y[i] - X[i] * W1[i]) / (1.0 - X[i]);
      }
      salpha[s] = alpha;
      slast[s] = std::exp(logpi(K - 1));
      int occ = 0;
      for (int k = 0; k < K; k++) occ += cnt[k] > 0;
      socc[s] = occ;
      slastn[s] = cnt[K - 1];
      s++;
    }
    if (it % 100 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("W1") = sW1, Rcpp::Named("W2") = sW2,
                            Rcpp::Named("alpha") = salpha, Rcpp::Named("pi_last") = slast,
                            Rcpp::Named("n_clusters") = socc,
                            Rcpp::Named("n_last") = slastn);
}
