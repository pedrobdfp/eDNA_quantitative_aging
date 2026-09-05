// =============================================================================
// conc_age_simplified.stan
// -----------------------------------------------------------------------------
// Age estimation from replicate CONCENTRATION measurements.
//
// This is the general-purpose form of the model. It takes whatever quantity
// your instrument reports -- a qPCR standard-curve concentration, a ddPCR
// concentration, anything on a copies-per-volume scale -- together with a
// detected / not-detected flag for every replicate. It never sees droplets, so
// it is not tied to ddPCR.
//
// PROCESS  (identical to the ddPCR version)
//     log C_ij = C[i] + p[j] + r[j] * t[i]                  (log copies / L)
//
//   for unit i (a sample, or a grab of several samples) and marker j. t[i] is
//   the time since shedding and C[i] the baseline log-concentration at t = 0.
//   The marker offsets p[j] and decay rates r[j] are supplied as data.
//
// OBSERVATION
//     log C_ijs = log C_ij + eta_sj                 eta ~ Normal(0, sigma_bio)
//
//     z_ijrs    ~ Bernoulli(logit^-1(beta * (log C_ijs - logC50)))
//     y_ijrs    ~ Normal(log C_ijs, sigma_tech)              where z = 1
//
//   s indexes a biological replicate (a separate water sample from the same
//   unit) and r a technical replicate (a separate PCR well from the same
//   extract). sigma_bio is the SD between water samples; sigma_tech the SD
//   between wells of one sample.
//
//   Every replicate contributes a detection term, so wells that amplified
//   nothing are informative rather than missing: a non-detection says the
//   concentration was low, which is exactly what a long-decayed marker looks
//   like. Only detected replicates contribute a quantitative term.
//
// WHY THERE IS NO MEAN CORRECTION HERE
//   The ddPCR version needs a -sigma^2/2 term because its likelihood applies
//   the exponential link to a lognormal quantity. Here the likelihood is
//   written directly on the log scale and is symmetric, so no correction is
//   needed and both sigmas may be interpreted directly.
//
// DETECTION PARAMETERS
//   Detection is parameterised by its 50% point rather than by an intercept.
//   logC50 is the log-concentration at which half of replicates amplify -- the
//   limit of detection of the assay -- and beta is the slope there. The usual
//   intercept, alpha = -beta * logC50, is reported in generated quantities.
//
//   Written this way the two parameters are close to independent. In the
//   intercept-and-slope form they are strongly correlated, because moving the
//   slope moves the whole curve sideways unless the intercept compensates, and
//   the sampler spends most of its time following that ridge.
//
//   Set the priors from R, so that the same file serves a well-characterised
//   calibration setting and a field setting.
//
// use_bio = 0 removes the biological level entirely, for designs with
// technical replication only. N_bio and bio_idx are then ignored.
// =============================================================================

data {
  int<lower=1> Nt;                                  // units being aged
  int<lower=1> Nloci;                               // markers
  vector[Nloci] r;                                  // decay rates (known)
  vector[Nloci] p;                                  // t = 0 log offsets (known)

  int<lower=0> N;                                   // replicates, all of them
  array[N] int<lower=1, upper=Nt>    obs_i;
  array[N] int<lower=1, upper=Nloci> obs_j;
  array[N] int<lower=0, upper=1>     z;             // 1 = detected

  int<lower=0> N_y;                                 // detected replicates
  array[N_y] int<lower=1> y_row;                    // their row in 1..N
  vector[N_y] y_obs;                                // log copies / L

  int<lower=0, upper=1> use_bio;                    // include biological level
  int<lower=1> N_bio;                               // biological replicates
  array[N] int<lower=1> bio_idx;                    // replicate of each row

  real C0_mean;
  real<lower=0> C0_sd;
  real<lower=0> t_mean;
  real<lower=0> t_sd;
  real logC50_mean;
  real<lower=0> logC50_sd;
  real beta_mean;
  real<lower=0> beta_sd;
  real<lower=0> sigma_sd;                           // half-normal prior scale
}

transformed data {
  int N_bio_eff = use_bio ? N_bio : 0;
  int Nloci_eff = use_bio ? Nloci : 0;
}

parameters {
  vector<lower=0>[Nt] t;                            // age of each unit
  vector<lower=0>[Nt] C;                            // baseline log-conc at t = 0
  real logC50;                                      // log-conc at 50% detection
  real<lower=0> beta;                               // detection slope
  real<lower=0> sigma_tech;                         // well-to-well SD
  vector<lower=0>[use_bio ? 1 : 0] sigma_bio;       // sample-to-sample SD
  matrix[N_bio_eff, Nloci_eff] eta_raw;             // biological deviates
}

transformed parameters {
  matrix[Nt, Nloci] mu;
  vector[N] level;                                  // log C of each replicate's sample
  real sb = use_bio ? sigma_bio[1] : 0;

  for (i in 1:Nt) {
    for (j in 1:Nloci) {
      mu[i, j] = C[i] + p[j] + (r[j] * t[i]);
    }
  }

  for (n in 1:N) {
    level[n] = mu[obs_i[n], obs_j[n]];
    if (use_bio) {
      level[n] += sb * eta_raw[bio_idx[n], obs_j[n]];
    }
  }
}

model {
  t          ~ normal(t_mean, t_sd);
  C          ~ normal(C0_mean, C0_sd);
  logC50     ~ normal(logC50_mean, logC50_sd);
  beta       ~ normal(beta_mean, beta_sd);
  sigma_tech ~ normal(0, sigma_sd);

  if (use_bio) {
    sigma_bio          ~ normal(0, sigma_sd);
    to_vector(eta_raw) ~ std_normal();
  }

  z     ~ bernoulli_logit(beta * (level - logC50));
  y_obs ~ normal(level[y_row], sigma_tech);
}

generated quantities {
  vector[N] log_lik;
  vector[N] p_detect;
  real sigma_total = sqrt(square(sigma_tech) + square(sb));
  real alpha = -beta * logC50;                      // intercept of the same function

  {
    vector[N] ll;
    for (n in 1:N) {
      p_detect[n] = inv_logit(beta * (level[n] - logC50));
      ll[n] = bernoulli_logit_lpmf(z[n] | beta * (level[n] - logC50));
    }
    for (n in 1:N_y) {
      ll[y_row[n]] += normal_lpdf(y_obs[n] | level[y_row[n]], sigma_tech);
    }
    log_lik = ll;
  }
}

