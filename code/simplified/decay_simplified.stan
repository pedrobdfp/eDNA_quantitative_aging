// =============================================================================
// decay_simplified.stan
// -----------------------------------------------------------------------------
// Marker-specific first-order decay rates from a controlled decay experiment,
// using the same concentration observation model as conc_age_simplified.stan
// so that the rates are estimated under the likelihood they are later used
// with.
//
// PROCESS
//     log C_ikt = C0[i,k] - lambda[k] * t                   (log copies / L)
//
//   for vessel i, marker k and elapsed time t. Each vessel-marker pair has its
//   own starting concentration; lambda[k] is constrained positive and
//   subtracted, so concentration declines.
//
// OBSERVATION
//     z_n ~ Bernoulli(logit^-1(beta * (log C_n - logC50)))
//     y_n ~ Normal(log C_n, sigma_obs)                      where z = 1
//
//   Inputs are concentrations, not droplet counts, so this file works for qPCR
//   and ddPCR alike. A decay experiment usually has technical replication only
//   -- one vessel sampled repeatedly, several PCR wells per sample -- so there
//   is a single observation SD rather than the nested pair used in the field.
//
// OUTPUTS USED DOWNSTREAM
//   r      the decay rates written negative, the convention the age model uses
//   p_log  each marker's t = 0 log-concentration offset relative to marker 1,
//          so p_log[1] is exactly 0
// =============================================================================

data {
  int<lower=1> N;                                  // all replicate observations
  int<lower=1> N_ik;                               // vessel x marker combinations
  int<lower=1> N_k;                                // markers

  array[N] int<lower=1, upper=N_ik> ik;
  array[N] int<lower=1, upper=N_k>  k;
  vector[N] time;                                  // hours since t = 0
  array[N] int<lower=0, upper=1> z;                // 1 = detected

  int<lower=0> N_y;                                // detected observations
  array[N_y] int<lower=1, upper=N> y_row;
  vector[N_y] y_obs;                               // log copies / L

  array[N_ik] int<lower=1, upper=N_k> ik_to_k;

  real C0_mean;
  real<lower=0> C0_sd;
  real logC50_mean;
  real<lower=0> logC50_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real<lower=0> sigma_sd;
  real<lower=0> lambda_sd;

  int<lower=0> N_time_sim;
  vector[N_time_sim] time_sim;
}

parameters {
  vector[N_ik] C0;                 // log-concentration at t = 0, per vessel x marker
  vector<lower=0>[N_k] lambda;     // decay rate per marker (1/hour, positive)
  real logC50;                     // log-conc at 50% detection
  real<lower=0> beta;              // detection slope
  real<lower=0> sigma_obs;         // observation SD of log-concentration
}

transformed parameters {
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = C0[ik[n]] - lambda[k[n]] * time[n];
  }
}

model {
  C0        ~ normal(C0_mean, C0_sd);
  lambda    ~ normal(0, lambda_sd);
  logC50    ~ normal(logC50_mean, logC50_sd);
  beta      ~ normal(beta_mean, beta_sd);
  sigma_obs ~ normal(0, sigma_sd);

  z     ~ bernoulli_logit(beta * (mu - logC50));
  y_obs ~ normal(mu[y_row], sigma_obs);
}

generated quantities {
  vector[N_k] r;                   // decay rate, negative (age-model convention)
  vector[N_k] half_life;           // hours
  vector[N_k] C0_bar;              // mean t = 0 log-concentration per marker
  vector[N_k] p_log;               // t = 0 log offset relative to marker 1
  vector[N_k] p_ratio;             // the same thing as a concentration ratio
  matrix[N_time_sim, N_k] C_sim;   // mean decay curve per marker
  vector[N] log_lik;
  real alpha = -beta * logC50;

  r = -lambda;

  for (kk in 1:N_k) {
    half_life[kk] = log(2) / lambda[kk];
    {
      real s = 0;
      int n_ik = 0;
      for (i in 1:N_ik) {
        if (ik_to_k[i] == kk) {
          s += C0[i];
          n_ik += 1;
        }
      }
      C0_bar[kk] = n_ik == 0 ? negative_infinity() : s / n_ik;
    }
  }

  for (kk in 1:N_k) {
    p_log[kk]   = C0_bar[kk] - C0_bar[1];
    p_ratio[kk] = exp(p_log[kk]);
    for (tt in 1:N_time_sim) {
      C_sim[tt, kk] = C0_bar[kk] - lambda[kk] * time_sim[tt];
    }
  }

  {
    vector[N] ll;
    for (n in 1:N) ll[n] = bernoulli_logit_lpmf(z[n] | beta * (mu[n] - logC50));
    for (n in 1:N_y) ll[y_row[n]] += normal_lpdf(y_obs[n] | mu[y_row[n]], sigma_obs);
    log_lik = ll;
  }
}

