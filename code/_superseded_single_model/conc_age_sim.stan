// =============================================================================
// conc_age_sim.stan
// -----------------------------------------------------------------------------
// Joint age-estimation model used for the simulation study and the
// leave-one-out carboy validation. This variant uses tight priors on the
// detection parameters (alpha, beta) and measurement noise (sigma), which is
// appropriate when the data-generating process (simulation) or the
// calibration context (carboy) is well-characterised.
//
// For the field-data analysis we use conc_age_field.stan, which has the same
// likelihood but wider priors on alpha / beta / sigma.
//
// Notation
//   t[i]   : age (hours) of sample i
//   C[i]   : baseline log eDNA concentration of sample i at t = 0
//   p[j]   : log-concentration offset of marker j (relative to reference marker)
//   r[j]   : first-order decay rate of marker j (1/hour; negative)
//   mu[i,j] = C[i] + p[j] + r[j] * t[i]
//
// Two likelihood components per observation:
//   (1) Bernoulli detection: z ~ bernoulli_logit(alpha + beta * mu)
//   (2) Gaussian concentration (when detected): y_obs ~ normal(mu, sigma)
// =============================================================================

data {
  int<lower=1> Nt;                                  // number of samples
  int<lower=1> Nloci;                               // number of markers
  vector[Nloci] r;                                  // decay rates (known)

  // Detection / non-detection observations (all PCR replicates)
  int<lower=0> N_obs;
  array[N_obs] int<lower=1, upper=Nt>    obs_i;
  array[N_obs] int<lower=1, upper=Nloci> obs_j;
  array[N_obs] int<lower=0, upper=1>     z;

  // Quantitative observations (detected only)
  int<lower=0> N_y;
  array[N_y] int<lower=1, upper=Nt>    y_i;
  array[N_y] int<lower=1, upper=Nloci> y_j;
  vector[N_y] y_obs;

  // Marker proportions at t = 0 (known)
  vector[Nloci] p;

  // Priors
  real C0_mean;
  real C0_sd;
  real<lower=0> t_mean;
  real<lower=0> t_sd;
}

parameters {
  vector<lower=0>[Nt] t;       // sample age
  vector<lower=0>[Nt] C;       // baseline log-concentration at t=0
  real alpha;                  // detection intercept
  real<lower=0> beta;          // detection slope on log-concentration
  real<lower=0> sigma;         // measurement SD of log-concentration
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
  // Priors (tight; simulation/LOO regime)
  t     ~ normal(t_mean, t_sd);
  C     ~ normal(C0_mean, C0_sd);
  alpha ~ normal(-3, 0.5);
  beta  ~ normal(4, 1);
  sigma ~ normal(0, 0.5);

  // Likelihood
  for (n in 1:N_obs) {
    z[n] ~ bernoulli_logit(alpha + beta * mu[obs_i[n], obs_j[n]]);
  }
  for (n in 1:N_y) {
    y_obs[n] ~ normal(mu[y_i[n], y_j[n]], sigma);
  }
}
