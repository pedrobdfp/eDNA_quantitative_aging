// =============================================================================
// decay_simplified.stan
//
// Marker-specific first-order decay rates, estimated from a controlled
// experiment using concentration measurements.
// -----------------------------------------------------------------------------
//
// PURPOSE
//
//   The age model (conc_age_simplified.stan) needs to know two things about
//   each marker before it can date anything: how quickly the marker
//   disappears, and where it starts relative to the others. This file measures
//   both, from an experiment in which water is held in a sealed vessel and
//   subsampled over time.
//
//   Here the elapsed time is known and the decay rate is estimated. In the age
//   model it is the other way round: the rate is known and the time is
//   estimated. The same observation model is used in both, so the rates are
//   measured under the same likelihood that later consumes them.
//
//   To characterise a new system: hold water containing the target in a
//   sealed vessel under conditions resembling the field and subsample it
//   over a period long enough to observe a clear decline. Denser sampling
//   early is more informative than a few distant timepoints.
//
// THE PROCESS MODEL
//
//   Each vessel-and-marker combination has its own starting level. Within one,
//   concentration falls exponentially, which on the log scale is a straight
//   line:
//
//       log C[n] = C0[ik[n]] - lambda[k[n]] * time[n]
//
//   in log copies per litre, where
//
//       C0      the log concentration at the start, one per vessel and marker,
//               estimated as the intercept of that line
//       lambda  the decay rate of a marker, per hour, estimated. It is held
//               positive and subtracted, so a larger lambda means faster loss.
//       time    hours since the start of the experiment, known
//
// OBSERVATION MODEL
//
//       z[n] ~ Bernoulli(logit^-1(beta * (log C[n] - logC50)))
//       y[n] ~ Normal(log C[n], sigma_obs)                    where z = 1
//
//   Every replicate contributes a detection outcome; those that detected
//   something contribute a measurement as well. Keep your non-detections: near
//   the end of a decay experiment they are the observations that pin down how
//   far the concentration has fallen.
//
//   sigma_obs is the spread between replicates measuring the same sample. A
//   decay experiment usually has technical replication only -- one vessel
//   sampled repeatedly, several PCR replicates per sample -- so a single
//   observation spread is enough here, unlike the field model which separates
//   variation between water samples from variation between replicates.
//
//   logC50 is the log concentration at which half of replicates amplify, the
//   effective limit of detection, and beta is how sharply detection switches on
//   there. The equivalent intercept of the more familiar alpha + beta * log C
//   form is alpha = -beta * logC50 and is reported below.
//
// QUANTITIES PASSED TO THE AGE MODEL
//
//   r      the decay rates written as negative numbers, which is the sign
//          convention the age model expects
//   p_log  each marker's starting level relative to the first marker, taken
//          from the fitted intercepts. This is the right quantity to pass on,
//          because the age model assumes a straight line of slope r, so the
//          offset it needs is the intercept of that same line.
// =============================================================================


data {
  // ---- Dimensions ---------------------------------------------------------
  int<lower=1> N;                        // number of PCR replicates
  int<lower=1> N_ik;                     // vessel x marker combinations
  int<lower=1> N_k;                      // number of markers

  // ---- Observations: one row per PCR replicate ----------------------------
  array[N] int<lower=1, upper=N_ik> ik;  // which vessel-and-marker this belongs to
  array[N] int<lower=1, upper=N_k>  k;   // which marker it measured
  vector[N] time;                        // hours since the start of the experiment
  array[N] int<lower=0, upper=1> z;      // 1 if it amplified, 0 if not

  // ---- Detected replicates only -------------------------------------------
  int<lower=0> N_y;                      // how many replicates detected something
  array[N_y] int<lower=1, upper=N> y_row;   // their position in the list above
  vector[N_y] y_obs;                        // their measurements, log copies per litre

  // ---- Marker of each vessel-by-marker combination ------------------------
  array[N_ik] int<lower=1, upper=N_k> ik_to_k;

  // ---- Prior hyperparameters ----------------------------------------------
  real C0_mean;                          // expected starting log concentration
  real<lower=0> C0_sd;                   // how uncertain that expectation is
  real logC50_mean;                      // expected limit of detection, log scale
  real<lower=0> logC50_sd;               // how uncertain that expectation is
  real beta_mean;                        // expected sharpness of the detection curve
  real<lower=0> beta_sd;                 // how uncertain that expectation is
  real<lower=0> sigma_sd;                // scale of the prior on the observation spread
  real<lower=0> lambda_sd;               // scale of the prior on the decay rates

  // ---- Prediction grid for the fitted decay curve -------------------------
  int<lower=0> N_time_sim;
  vector[N_time_sim] time_sim;
}


parameters {
  vector[N_ik] C0;                       // starting log concentration, per vessel and marker
  vector<lower=0>[N_k] lambda;           // decay rate of each marker, per hour
  real logC50;                           // log concentration at 50% detection
  real<lower=0> beta;                    // steepness of the detection curve
  real<lower=0> sigma_obs;               // spread between replicates of one sample
}


transformed parameters {
  vector[N] mu;                          // expected log concentration of each replicate

  for (n in 1:N) {
    mu[n] = C0[ik[n]] - lambda[k[n]] * time[n];
  }
}


model {
  // ---- Priors -------------------------------------------------------------
  C0        ~ normal(C0_mean, C0_sd);
  lambda    ~ normal(0, lambda_sd);      // half-normal: lambda cannot be negative
  logC50    ~ normal(logC50_mean, logC50_sd);
  beta      ~ normal(beta_mean, beta_sd);   // half-normal
  sigma_obs ~ normal(0, sigma_sd);          // half-normal

  // ---- Likelihood ---------------------------------------------------------
  z     ~ bernoulli_logit(beta * (mu - logC50));   // every replicate
  y_obs ~ normal(mu[y_row], sigma_obs);            // those that detected something
}


generated quantities {
  vector[N_k] r;                         // decay rates in the age model's convention
  vector[N_k] half_life;                 // hours for a marker to halve
  vector[N_k] C0_bar;                    // average starting level of each marker
  vector[N_k] p_log;                     // starting level relative to marker 1
  vector[N_k] p_ratio;                   // the same thing as a plain ratio
  matrix[N_time_sim, N_k] C_sim;         // fitted decay curve, for plotting
  vector[N] log_lik;                     // fit of each replicate, for model comparison
  real alpha = -beta * logC50;           // intercept of the same detection curve

  r = -lambda;

  for (kk in 1:N_k) {
    half_life[kk] = log(2) / lambda[kk];

    // Average the intercepts of every vessel measured for this marker.
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
    // Offsets are relative, so the first marker is the reference and is zero
    // by construction.
    p_log[kk]   = C0_bar[kk] - C0_bar[1];
    p_ratio[kk] = exp(p_log[kk]);

    for (tt in 1:N_time_sim) {
      C_sim[tt, kk] = C0_bar[kk] - lambda[kk] * time_sim[tt];
    }
  }

  {
    vector[N] ll;
    for (n in 1:N) ll[n] = bernoulli_logit_lpmf(z[n] | beta * (mu[n] - logC50));
    // Replicates that detected something contribute their measurement too.
    for (n in 1:N_y) ll[y_row[n]] += normal_lpdf(y_obs[n] | mu[y_row[n]], sigma_obs);
    log_lik = ll;
  }
}
