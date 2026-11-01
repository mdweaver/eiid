// Truncated Dirichlet-process (stick-breaking) mixture of normals for 2x2
// ecological inference (Imai, Lu & Strauss 2008), latent-rate formulation.
//   p = 2: CAR,  rows (logit W1, logit W2)
//   p = 3: NCAR, rows (logit W1, logit W2, logit X)
//   K = 1: the parametric model.
// Each latent unit's W1 is a parameter bounded by its Duncan-Davis interval;
// W2 is implied by the accounting identity. The density of W1 on the line,
// w.r.t. Lebesgue measure in w1, is f_W(w1, w2(w1)) up to a data constant,
// with f_W(w) = phi(logit w) / prod_j w_j (1 - w_j).
// Components use a Bartlett (non-centred) parametrisation of the
// normal-inverse-Wishart base measure: Sigma^{-1} = M M', M = chol(S0^{-1}) A.
data {
  int<lower=2, upper=3> p;
  int<lower=1> K;
  int<lower=0> n;                 // latent units, 0 < X < 1, 0 < Y < 1
  vector<lower=0, upper=1>[n] X;
  vector<lower=0, upper=1>[n] Y;
  vector<lower=0, upper=1>[n] L;  // bounds of W1 (shrunk slightly)
  vector<lower=0, upper=1>[n] U;
  int<lower=0> n_obs;             // fully observed units (logit scale, p rows)
  matrix[p, n_obs] Wobs;
  int<lower=0> n_x1;              // CAR only: X = 1, logit W1 observed
  vector[n_x1] w1obs;
  int<lower=0> n_x0;              // CAR only: X = 0, logit W2 observed
  vector[n_x0] w2obs;
  vector[p] mu0;
  real<lower=0> tau0;
  real<lower=p - 1> nu0;
  cov_matrix[p] S0;
  real<lower=0> a0;
  real<lower=0> b0;
}
transformed data {
  matrix[p, p] LS = cholesky_decompose(inverse_spd(S0));
  int n_off = (p * (p - 1)) %/% 2;
  row_vector[n] lx;
  for (i in 1:n) lx[i] = logit(X[i]);
}
parameters {
  vector<lower=L, upper=U>[n] W1;
  real<lower=0> alpha;
  vector<lower=0, upper=1>[K - 1] v;
  array[K] vector<lower=0>[p] csq;   // chi-square diagonals of the Bartlett factor
  array[K] vector[n_off] zoff;       // standard-normal off-diagonals
  array[K] vector[p] zmu;            // non-centred component means
}
transformed parameters {
  vector[n] W2 = (Y - X .* W1) ./ (1 - X);
}
model {
  vector[K] logpi;
  matrix[K, n] lp;
  matrix[K, n_obs] lpo;
  matrix[K, n_x1] lp1;
  matrix[K, n_x0] lp0;
  matrix[p, n] Ws;
  Ws[1] = to_row_vector(logit(W1));
  Ws[2] = to_row_vector(logit(W2));
  if (p == 3) Ws[3] = lx;

  alpha ~ gamma(a0, b0);
  v ~ beta(1, alpha);
  if (K == 1) {
    logpi[1] = 0;
  } else {
    real acc = 0;
    for (k in 1:(K - 1)) {
      logpi[k] = log(v[k]) + acc;
      acc += log1m(v[k]);
    }
    logpi[K] = acc;
  }

  for (k in 1:K) {
    matrix[p, p] A = rep_matrix(0, p, p);
    matrix[p, p] M;
    vector[p] mu;
    int idx = 1;
    for (j in 1:p) {
      csq[k][j] ~ chi_square(nu0 - j + 1);
      A[j, j] = sqrt(csq[k][j]);
      for (i in (j + 1):p) {
        A[i, j] = zoff[k][idx];
        idx += 1;
      }
    }
    zoff[k] ~ std_normal();
    zmu[k] ~ std_normal();
    M = LS * A;                       // Sigma^{-1} = M M'
    // mu = mu0 + Sigma^{1/2} zmu / sqrt(tau0), with Sigma^{1/2} = M^{-T}
    mu = mu0 + to_vector(mdivide_right_tri_low(to_row_vector(zmu[k]), M)) / sqrt(tau0);
    {
      real ldet = sum(log(diagonal(M)));   // = -0.5 log det Sigma
      real c = -0.5 * p * log(2 * pi()) + ldet;
      matrix[p, n] Z = M' * (Ws - rep_matrix(mu, n));
      lp[k] = logpi[k] + c - 0.5 * columns_dot_self(Z);
      if (n_obs > 0) {
        matrix[p, n_obs] Zo = M' * (Wobs - rep_matrix(mu, n_obs));
        lpo[k] = logpi[k] + c - 0.5 * columns_dot_self(Zo);
      }
      if (n_x1 + n_x0 > 0) {
        matrix[p, p] Minv = mdivide_left_tri_low(M, diag_matrix(rep_vector(1, p)));
        matrix[p, p] Sig = crossprod(Minv);
        if (n_x1 > 0)
          for (i in 1:n_x1) lp1[k, i] = logpi[k] + normal_lpdf(w1obs[i] | mu[1], sqrt(Sig[1, 1]));
        if (n_x0 > 0)
          for (i in 1:n_x0) lp0[k, i] = logpi[k] + normal_lpdf(w2obs[i] | mu[2], sqrt(Sig[2, 2]));
      }
    }
  }
  for (i in 1:n) target += log_sum_exp(col(lp, i));
  for (i in 1:n_obs) target += log_sum_exp(col(lpo, i));
  for (i in 1:n_x1) target += log_sum_exp(col(lp1, i));
  for (i in 1:n_x0) target += log_sum_exp(col(lp0, i));
  // Jacobian of the logit transform of the latent rates
  target += -sum(log(W1) + log1m(W1)) - sum(log(W2) + log1m(W2));
}
generated quantities {
  real pi_last = K == 1 ? 1 : exp(sum(log1m(v)));
}
