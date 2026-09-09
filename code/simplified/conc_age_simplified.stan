// =============================================================================
// conc_age_simplified.stan
//
// HOW OLD IS THIS eDNA? Estimating time since shedding from concentration
// measurements. Works with qPCR or digital PCR.
// -----------------------------------------------------------------------------
//
// START HERE IF YOU ARE NEW TO THE METHOD. This is the general-purpose version.
// It reads whatever your instrument reports on a copies-per-volume scale,
// together with a yes/no detection flag for each replicate. It never sees
// droplets, so nothing about it is specific to ddPCR.
//
// THE IDEA IN ONE PARAGRAPH
//
//   eDNA starts breaking down the moment it is released. If you measure several
//   genetic markers from the same organism, and those markers break down at
//   different speeds, then the balance between them tells you how long the
//   breakdown has been going on. A sample where the fragile markers have
//   already faded is old; a sample where all markers remain in their original
//   proportions is fresh. That balance does not depend on how much DNA was
//   released to begin with, which is a quantity you almost never know.
//
// THE PROCESS MODEL
//
//   For unit i (one water sample, or a group collected together) and marker j:
//
//       log C[i,j] = C[i] + p[j] + r[j] * t[i]
//
//   in log copies per litre of water, where
//
//       t[i]  the age of the sample in hours          <- ESTIMATED, the answer
//       C[i]  its log concentration when shed         <- estimated
//       p[j]  marker j's starting level relative to
//             the first marker (so p[1] = 0)          <- supplied as data
//       r[j]  marker j's decay rate per hour,
//             a negative number                       <- supplied as data
//
//   Age and concentration do not get confused with one another because they
//   leave different fingerprints. C[i] raises or lowers every marker by the
//   same amount. t[i] is multiplied by r[j], which differs between markers, so
//   age tilts the pattern rather than shifting it. If all markers decayed at
//   the same rate the two would be indistinguishable and the method would not
//   work.
//
//   r and p come from a decay experiment, fitted by decay_simplified.stan, and
//   are treated as known here.
//
// WHAT THE MODEL DOES WITH YOUR MEASUREMENTS
//
//   Two things vary between replicates, and they are kept separate:
//
//     1. Water samples from the same unit differ, because eDNA is patchy.
//        That is eta, with standard deviation sigma_bio.
//
//     2. PCR replicates of the same water sample differ. That is sigma_tech.
//
//   Writing s for a water sample and r for a replicate:
//
//       log C[i,j,s] = log C[i,j] + eta[s,j]        eta ~ Normal(0, sigma_bio)
//
//   Every replicate then contributes a detection outcome, and those that
//   detected something contribute a measurement as well:
//
//       z[i,j,r,s] ~ Bernoulli(logit^-1(beta * (log C[i,j,s] - logC50)))
//       y[i,j,r,s] ~ Normal(log C[i,j,s], sigma_tech)          where z = 1
//
//   KEEP YOUR NON-DETECTIONS. A replicate that amplified nothing is not
//   missing data. It says the concentration was below what your assay can see,
//   which is exactly what a marker that has been decaying a long time looks
//   like. Every replicate contributes a detection term whether or not it
//   amplified; only those that did contribute a measurement. Deleting them
//   biases ages young.
//
// THE DETECTION FUNCTION
//
//   Detection is described by two numbers with direct physical meaning:
//
//       logC50  the log concentration at which half of replicates amplify,
//               that is, the effective limit of detection of your assay
//       beta    how sharply detection switches on around that point
//
//   The equivalent intercept of the more familiar form alpha + beta * log C is
//   alpha = -beta * logC50, and is reported in generated quantities. Both
//   describe the same curve. Priors for logC50 and beta are set in R, so this
//   one file serves a well-characterised assay and a poorly-characterised one.
//
// WHY THE TWO STANDARD DEVIATIONS ARE SHARED BETWEEN MARKERS
//
//   sigma_bio and sigma_tech are single numbers applied to every marker rather
//   than one per marker. A marker-specific spread would shift each marker's
//   estimated concentration by a different amount, which tilts the marker
//   pattern -- and the marker pattern is exactly what carries the age.
//
// THERE IS NO MEAN CORRECTION IN THIS FILE
//
//   The ddPCR version of this model needs a -sigma^2/2 term because its
//   likelihood applies an exponential to a quantity carrying log-scale noise.
//   Here the likelihood is written directly on the log scale and is symmetric,
//   so no correction is needed and both standard deviations can be read
//   straight off.
//
// USING THIS FILE WITHOUT BIOLOGICAL REPLICATION
//
//   Set use_bio = 0 when each unit is a single water sample. The biological
//   level then disappears completely and only sigma_tech is estimated. N_bio
//   and bio_idx are ignored.
// =============================================================================


data {
  // ---- What is being aged --------------------------------------------------
  int<lower=1> Nt;                       // number of units to date
  int<lower=1> Nloci;                    // number of markers

  // ---- Known from the decay experiment -------------------------------------
  vector[Nloci] r;                       // decay rate per marker, per hour, negative
  vector[Nloci] p;                       // starting level per marker, relative to
                                         //   marker 1, on the log scale; p[1] = 0

  // ---- Every replicate, detected or not ------------------------------------
  int<lower=0> N;                                   // number of PCR replicates
  array[N] int<lower=1, upper=Nt>    obs_i;         // which unit this replicate belongs to
  array[N] int<lower=1, upper=Nloci> obs_j;         // which marker it measured
  array[N] int<lower=0, upper=1>     z;             // 1 if it amplified, 0 if not

  // ---- The subset that produced a number -----------------------------------
  int<lower=0> N_y;                      // how many replicates detected something
  array[N_y] int<lower=1> y_row;         // their position in the list above
  vector[N_y] y_obs;                     // their measurements, log copies per litre

  // ---- Grouping of replicates into water samples ---------------------------
  int<lower=0, upper=1> use_bio;         // 1 = model between-sample variation
  int<lower=1> N_bio;                    // number of water samples
  array[N] int<lower=1> bio_idx;         // which water sample each replicate came from

  // ---- Prior settings, chosen in R -----------------------------------------
  real C0_mean;                          // expected log concentration when shed
  real<lower=0> C0_sd;                   // how uncertain that expectation is
  real<lower=0> t_mean;                  // expected age in hours
  real<lower=0> t_sd;                    // how uncertain that expectation is
  real logC50_mean;                      // expected limit of detection, log scale
  real<lower=0> logC50_sd;               // how uncertain that expectation is
  real beta_mean;                        // expected sharpness of the detection curve
  real<lower=0> beta_sd;                 // how uncertain that expectation is
  real<lower=0> sigma_sd;                // scale of the prior on both noise terms
}


transformed data {
  // Stan needs array sizes fixed before sampling begins. Setting these to zero
  // when use_bio = 0 makes the biological deviations vanish entirely: they
  // still exist in the code, but as empty containers that cost nothing and
  // contribute nothing.
  int N_bio_eff = use_bio ? N_bio : 0;
  int Nloci_eff = use_bio ? Nloci : 0;
}


parameters {
  // ---- The quantities of interest ------------------------------------------
  vector<lower=0>[Nt] t;                 // age of each unit, hours, cannot be negative
  vector<lower=0>[Nt] C;                 // log concentration of each unit when shed

  // ---- The detection curve -------------------------------------------------
  real logC50;                           // log concentration at 50% detection
  real<lower=0> beta;                    // steepness of the curve there

  // ---- How much things vary ------------------------------------------------
  real<lower=0> sigma_tech;              // SD between replicates of one water sample
  vector<lower=0>[use_bio ? 1 : 0] sigma_bio;
                                         // SD between water samples of one unit.
                                         //   Declared as a vector of length 1 (or 0
                                         //   when switched off) so it can be removed
                                         //   cleanly; sb below unwraps it.

  // ---- The individual departures -------------------------------------------
  // Held in standard units, mean 0 and SD 1, and multiplied by sigma_bio below.
  // Written this way rather than drawn directly at scale sigma_bio, which gives
  // the sampler an evenly shaped space to explore.
  matrix[N_bio_eff, Nloci_eff] eta_raw;  // one per water sample and marker
}


transformed parameters {
  matrix[Nt, Nloci] mu;                  // expected log concentration per unit and marker
  vector[N] level;                       // log concentration of each replicate's
                                         //   own water sample

  // sigma_bio as a plain number, or zero when the biological level is off.
  real sb = use_bio ? sigma_bio[1] : 0;

  // The process model, evaluated for every unit and marker.
  for (i in 1:Nt) {
    for (j in 1:Nloci) {
      mu[i, j] = C[i] + p[j] + (r[j] * t[i]);
    }
  }

  // Step out to the individual water sample behind each replicate.
  for (n in 1:N) {
    level[n] = mu[obs_i[n], obs_j[n]];
    if (use_bio) {
      level[n] += sb * eta_raw[bio_idx[n], obs_j[n]];
    }
  }
}


model {
  // ---- What we believe before seeing the data ------------------------------
  t ~ normal(t_mean, t_sd);              // truncated at 0 by the declaration above
  C ~ normal(C0_mean, C0_sd);            // likewise

  logC50 ~ normal(logC50_mean, logC50_sd);
  beta   ~ normal(beta_mean, beta_sd);   // half-normal: beta cannot be negative

  // Variation between replicates of one water sample.
  sigma_tech ~ normal(0, sigma_sd);      // half-normal

  // Variation between water samples of one unit. Both lines do nothing at all
  // when use_bio = 0, because the containers are then empty.
  sigma_bio          ~ normal(0, sigma_sd);
  to_vector(eta_raw) ~ std_normal();     // to_vector flattens the matrix so the
                                         //   same prior applies to every entry

  // ---- What the data say ---------------------------------------------------
  z     ~ bernoulli_logit(beta * (level - logC50));   // every replicate
  y_obs ~ normal(level[y_row], sigma_tech);           // those that detected something
}


generated quantities {
  vector[N] log_lik;                     // fit of each replicate, for model comparison
  vector[N] p_detect;                    // chance this replicate detects anything
  real sigma_total = sqrt(square(sigma_tech) + square(sb));   // both sources combined
  real alpha = -beta * logC50;           // intercept of the same detection curve

  {
    vector[N] ll;
    for (n in 1:N) {
      p_detect[n] = inv_logit(beta * (level[n] - logC50));
      ll[n] = bernoulli_logit_lpmf(z[n] | beta * (level[n] - logC50));
    }
    // Replicates that detected something contribute their measurement too.
    for (n in 1:N_y) {
      ll[y_row[n]] += normal_lpdf(y_obs[n] | level[y_row[n]], sigma_tech);
    }
    log_lik = ll;
  }
}
