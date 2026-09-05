// =============================================================================
// decay_monophasic_ddPCR.stan
// -----------------------------------------------------------------------------
// Marker-specific first-order decay rates estimated from the carboy
// experiment, using the same observation model as droplet_age_ddPCR.stan so
// that the rates are estimated under the likelihood they are later used with.
//
// PROCESS
//     log C_ikt = C_0[i,k] - lambda[k] * t
//
//   for carboy i, marker k and elapsed time t. Each carboy and marker has its
//   own starting concentration C_0, and lambda[k] is the decay rate of marker
//   k, constrained positive and subtracted so that concentration declines.
//
// OBSERVATION
//     log C_obs = log C_ikt + eps - sigma_obs^2 / 2
//                 eps ~ Normal(0, sigma_obs)
//     omega     = log C_obs - log(20) + log(1.818) + log(Dilution)
//                 - log(50 / Filtered) - 7.07
//     W ~ Binomial(U, 1 - exp(-exp(omega)))
//
//   The offset terms convert a seawater concentration into expected copies per
//   droplet given the reaction, elution, dilution and filtered volumes; they
//   are algebraically identical to the offset used for the field samples.
//   sigma_obs is the replicate-level SD of log concentration, shared across
//   markers so that the mean correction is constant and cannot displace lambda.
// =============================================================================

data {
  int N_obs_0;
  int N_obs_1;
  int N_ijk_1;
  int N_ik_0;
  int N_k_1;

  array[N_obs_0] int ik_idx_0;
  array[N_obs_1] int ijk_idx_1;
  array[N_ijk_1] int s_ik_idx_1;
  array[N_ijk_1] int s_k_idx_1;
  array[N_ik_0] int ik_to_k;

  array[N_ijk_1] real time;
  array[N_obs_0] real Dilution_0;
  array[N_obs_1] real Dilution_1;
  array[N_obs_0] real Filtered_0;
  array[N_obs_1] real Filtered_1;
  array[N_obs_0] int U_0;
  array[N_obs_1] int U_1;
  array[N_obs_0] int W_0;
  array[N_obs_1] int W_1;

  int N_time_sim;
  array[N_time_sim] real time_sim;

  real<lower=0> sigma_obs_sd;      // half-normal prior scale
}

transformed data {
  // marker index of every observation, from the existing index maps
  array[N_obs_0] int k_of_0;
  array[N_obs_1] int k_of_1;

  for (i in 1:N_obs_0) k_of_0[i] = ik_to_k[ik_idx_0[i]];
  for (i in 1:N_obs_1) k_of_1[i] = s_k_idx_1[ijk_idx_1[i]];
}

parameters {
  vector[N_ik_0]         C_0;        // log-concentration at t=0, per carboy x locus
  vector<lower=0>[N_k_1] lambda;     // monophasic decay, per locus (positive)
  real<lower=0> sigma_obs;           // observation noise, shared across markers
  vector[N_obs_0] eps0_raw;
  vector[N_obs_1] eps1_raw;
}

transformed parameters {
  vector[N_ijk_1] C;
  vector[N_obs_0] omega_0;
  vector[N_obs_1] omega_1;

  // decay curve
  for (j in 1:N_ijk_1) {
    C[j] = C_0[s_ik_idx_1[j]] - lambda[s_k_idx_1[j]] * time[j];
  }

  // t = 0 droplets
  for (i in 1:N_obs_0) {
    real s = sigma_obs;
    omega_0[i] =
        C_0[ik_idx_0[i]]
      - log(20)
      + log(1.818)
      + log(Dilution_0[i])
      - log(50 / Filtered_0[i])
      - 7.07
      + s * eps0_raw[i] - 0.5 * square(s);
  }

  // droplets, t > 0
  for (i in 1:N_obs_1) {
    real s = sigma_obs;
    omega_1[i] =
        C[ijk_idx_1[i]]
      - log(20)
      + log(1.818)
      + log(Dilution_1[i])
      - log(50 / Filtered_1[i])
      - 7.07
      + s * eps1_raw[i] - 0.5 * square(s);
  }
}

model {
  W_0 ~ binomial(U_0, inv_cloglog(omega_0));
  W_1 ~ binomial(U_1, inv_cloglog(omega_1));

  C_0       ~ normal(0, 4);
  lambda    ~ normal(0, 1);
  sigma_obs ~ normal(0, sigma_obs_sd);
  eps0_raw  ~ std_normal();
  eps1_raw  ~ std_normal();
}

generated quantities {
  matrix[N_time_sim, N_k_1] C_sim;
  array[N_obs_0] int W_pred_0;
  array[N_obs_1] int W_pred_1;
  vector[N_k_1] half_life;

  for (k in 1:N_k_1) half_life[k] = log(2) / lambda[k];

  for (i in 1:N_obs_0) {
    W_pred_0[i] = binomial_rng(U_0[i], inv_cloglog(omega_0[i]));
  }
  for (i in 1:N_obs_1) {
    W_pred_1[i] = binomial_rng(U_1[i], inv_cloglog(omega_1[i]));
  }

  // mean decay curve per marker
  for (t in 1:N_time_sim) {
    for (k in 1:N_k_1) {
      real sum_c0 = 0;
      int n_c0 = 0;
      for (ik in 1:N_ik_0) {
        if (ik_to_k[ik] == k) {
          sum_c0 += C_0[ik];
          n_c0 += 1;
        }
      }
      if (n_c0 == 0) {
        C_sim[t, k] = negative_infinity();
      } else {
        real C0_bar = sum_c0 / n_c0;
        C_sim[t, k] = C0_bar - lambda[k] * time_sim[t];
      }
    }
  }
}

