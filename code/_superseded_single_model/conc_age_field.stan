// =============================================================================
// conc_age_field.stan
// -----------------------------------------------------------------------------
// Same likelihood as conc_age_sim.stan (see that file for full notation), but
// with wider priors on alpha, beta, and sigma. We use this version for all
// field-data analyses, where the detection threshold and measurement noise
// are less well-constrained than in the simulation or carboy settings.
// =============================================================================

data {
  int<lower=1> Nt;
  int<lower=1> Nloci;
  vector[Nloci] r;

  int<lower=0> N_obs;
  array[N_obs] int<lower=1, upper=Nt>    obs_i;
  array[N_obs] int<lower=1, upper=Nloci> obs_j;
  array[N_obs] int<lower=0, upper=1>     z;

  int<lower=0> N_y;
  array[N_y] int<lower=1, upper=Nt>    y_i;
  array[N_y] int<lower=1, upper=Nloci> y_j;
  vector[N_y] y_obs;

  vector[Nloci] p;

  real C0_mean;
  real C0_sd;
  real<lower=0> t_mean;
  real<lower=0> t_sd;
}

parameters {
  vector<lower=0>[Nt] t;
  vector<lower=0>[Nt] C;
  real alpha;
  real<lower=0> beta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[Nt, Nloci] mu;
  for (i in 1:Nt) {
    for (j in 1:Nloci) {
      mu[i, j] = C[i] + p[j] + (r[j] * t[i]);
    }
  }
}

model {
  // Priors (wide; field regime)
  t     ~ normal(t_mean, t_sd);
  C     ~ normal(C0_mean, C0_sd);
  alpha ~ normal(-2, 3);
  beta  ~ normal(4, 2);
  sigma ~ normal(0, 1);

  // Likelihood
  for (n in 1:N_obs) {
    z[n] ~ bernoulli_logit(alpha + beta * mu[obs_i[n], obs_j[n]]);
  }
  for (n in 1:N_y) {
    y_obs[n] ~ normal(mu[y_i[n], y_j[n]], sigma);
  }
}
