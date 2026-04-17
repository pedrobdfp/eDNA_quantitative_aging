// Decay model: Monophasic exponential decay – estimate across carboys

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
}

parameters {
  vector[N_ik_0]        C_0;       // log‐concentration at t=0
  vector<lower=0>[N_k_1] lambda;   // monophasic decay (per component)
}

transformed parameters {
  vector[N_ijk_1] C;
  vector[N_obs_0] omega_0;
  vector[N_obs_1] omega_1;

  // t=0 droplets
  for (i in 1:N_obs_0) {
    omega_0[i] =
        C_0[ik_idx_0[i]]
      - log(20)
      + log(1.818)
      + log(Dilution_0[i])
      - log(50 / Filtered_0[i])
      - 7.07;
  }

  // decay curve
  for (j in 1:N_ijk_1) {
    C[j] = C_0[s_ik_idx_1[j]] - lambda[s_k_idx_1[j]] * time[j];
  }

  // droplets >0
  for (i in 1:N_obs_1) {
    omega_1[i] =
        C[ijk_idx_1[i]]
      - log(20)
      + log(1.818)
      + log(Dilution_1[i])
      - log(50 / Filtered_1[i])
      - 7.07;
  }
}

model {
  W_0 ~ binomial(U_0, inv_cloglog(omega_0));
  W_1 ~ binomial(U_1, inv_cloglog(omega_1));

  C_0    ~ normal(0, 4);
  lambda ~ normal(0, 1);
}

generated quantities {
  matrix[N_time_sim, N_k_1] C_sim;
  int W_pred_0[N_obs_0];
  int W_pred_1[N_obs_1];
  vector[N_obs_0] p_pred_0;
  vector[N_obs_1] p_pred_1;

  // posterior predictive for observed droplets
  for (i in 1:N_obs_0) {
    real p0 = inv_cloglog(omega_0[i]);
    W_pred_0[i] = binomial_rng(U_0[i], p0);
    p_pred_0[i] = W_pred_0[i] / U_0[i];
  }

  for (i in 1:N_obs_1) {
    real p1 = inv_cloglog(omega_1[i]);
    W_pred_1[i] = binomial_rng(U_1[i], p1);
    p_pred_1[i] = W_pred_1[i] / U_1[i];
  }

  // simulate mean decay curve per k (marker/component combo)
  for (t in 1:N_time_sim) {
    for (k in 1:N_k_1) {

      // mean C0 across all ik that belong to this k
      real sum_c0 = 0;
      int n_c0 = 0;
      for (ik in 1:N_ik_0) {
        if (ik_to_k[ik] == k) {
          sum_c0 += C_0[ik];
          n_c0 += 1;
        }
      }

      // safety: if mapping is wrong, avoid division by zero
      // (better to crash, but this keeps it defined)
      if (n_c0 == 0) {
        C_sim[t, k] = negative_infinity();
      } else {
        real C0_bar = sum_c0 / n_c0;
        C_sim[t, k] = C0_bar - lambda[k] * time_sim[t];
      }
    }
  }
}

