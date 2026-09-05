// =============================================================================
// droplet_age_ddPCR.stan
// -----------------------------------------------------------------------------
// Age estimation from ddPCR droplet counts, with observation error split into
// nested biological and technical components, each shared across markers.
//
// PROCESS
//     log C_ij = C[i] + p[j] + r[j] * t[i]                  (log copies / L)
//
//   for grab i and marker j. t[i] is the time since shedding and C[i] the
//   baseline log-concentration at t = 0. The marker offsets p[j] and decay
//   rates r[j] are supplied as data.
//
// OBSERVATION
//     log C_ijs  = log C_ij  + eta_sj                eta ~ Normal(0, sigma_bio)
//     log C_ijrs = log C_ijs + eps_ijrs - correction eps ~ Normal(0, sigma_tech)
//     omega      = log C_ijrs + log_offset
//     W ~ Binomial(U, 1 - exp(-exp(omega)))
//
//   s indexes a biological replicate (a separate water sample from the same
//   grab) and r a technical replicate (a separate well from the same extract).
//   sigma_bio is the SD between water samples of one grab; sigma_tech the SD
//   between wells of one sample.
//
// WHY BOTH SIGMAS ARE SHARED ACROSS MARKERS
//   Because the noise is lognormal, the value of C that the data identify
//   depends on sigma through a term sigma^2/2. If sigma differed by marker,
//   that term would differ by marker, which is a shift in the marker pattern --
//   and the marker pattern is precisely what identifies the age. With the
//   fitted per-marker values this moved ages by about four hours. Sharing the
//   two sigmas makes the correction a single constant, which C[i] absorbs, so
//   observation error widens the age posterior without displacing it.
//
// WHY THE NESTING MATTERS
//   Wells of one sample and samples of one grab are not exchangeable. Treating
//   them as a single flat set lets technical replication buy precision that has
//   not been earned, because variation shared by all wells of a sample does not
//   average away. Separating the two also lets each be reported.
//
// use_bio = 0 removes the biological level entirely, for designs with technical
// replication only. In that case N_bio and bio_idx are ignored.
// =============================================================================

data {
  int<lower=1> Nt;                                  // grabs
  int<lower=1> Nloci;                               // markers
  vector[Nloci] r;                                  // decay rates (known)
  vector[Nloci] p;                                  // t = 0 log offsets (known)

  int<lower=0> N;                                   // wells
  array[N] int<lower=1, upper=Nt>    obs_i;         // grab of each well
  array[N] int<lower=1, upper=Nloci> obs_j;         // marker of each well
  array[N] int<lower=0> W;                          // positive droplets
  array[N] int<lower=1> U;                          // accepted droplets
  vector[N] log_offset;

  int<lower=0, upper=1> use_bio;                    // include biological level
  int<lower=1> N_bio;                               // biological replicates
  array[N] int<lower=1> bio_idx;                    // replicate of each well

  real C0_mean;
  real<lower=0> C0_sd;
  real<lower=0> t_mean;
  real<lower=0> t_sd;
  real<lower=0> sigma_sd;                           // half-normal prior scale
}

transformed data {
  int N_bio_eff = use_bio ? N_bio : 0;
  int Nloci_eff = use_bio ? Nloci : 0;
}

parameters {
  vector<lower=0>[Nt] t;                            // age of each grab
  vector<lower=0>[Nt] C;                            // baseline log-conc at t = 0
  real<lower=0> sigma_tech;                         // well-to-well SD
  vector<lower=0>[use_bio ? 1 : 0] sigma_bio;       // sample-to-sample SD
  matrix[N_bio_eff, Nloci_eff] eta_raw;             // biological deviates
  vector[N] eps_raw;                                // technical deviates
}

transformed parameters {
  matrix[Nt, Nloci] mu;
  vector[N] omega;
  real sb = use_bio ? sigma_bio[1] : 0;

  // Constant across markers, so it relabels C without moving t.
  real mean_correction = 0.5 * (square(sigma_tech) + square(sb));

  for (i in 1:Nt) {
    for (j in 1:Nloci) {
      mu[i, j] = C[i] + p[j] + (r[j] * t[i]);
    }
  }

  for (n in 1:N) {
    real level = mu[obs_i[n], obs_j[n]];
    if (use_bio) {
      level += sb * eta_raw[bio_idx[n], obs_j[n]];
    }
    omega[n] = level + log_offset[n]
               + sigma_tech * eps_raw[n] - mean_correction;
  }
}

model {
  t          ~ normal(t_mean, t_sd);
  C          ~ normal(C0_mean, C0_sd);
  sigma_tech ~ normal(0, sigma_sd);
  eps_raw    ~ std_normal();

  if (use_bio) {
    sigma_bio         ~ normal(0, sigma_sd);
    to_vector(eta_raw) ~ std_normal();
  }

  W ~ binomial(U, inv_cloglog(omega));
}

generated quantities {
  vector[N] log_lik;
  array[N] int W_rep;
  vector[N] p_detect;
  real sigma_total = sqrt(square(sigma_tech) + square(sb));

  for (n in 1:N) {
    real pi_n = inv_cloglog(omega[n]);
    log_lik[n]  = binomial_lpmf(W[n] | U[n], pi_n);
    W_rep[n]    = binomial_rng(U[n], pi_n);
    p_detect[n] = 1 - exp(-U[n] * exp(omega[n]));
  }
}

