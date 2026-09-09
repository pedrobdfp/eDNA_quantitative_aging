// =============================================================================
// droplet_age_ddPCR.stan
//
// Time since shedding, estimated from ddPCR droplet counts.
//
// -----------------------------------------------------------------------------
//
// BACKGROUND
//
//   eDNA starts breaking down the moment it is released. If you measure several
//   genetic markers from the same organism, and those markers break down at
//   different speeds, then the balance between them tells you how long the
//   breakdown has been going on. A sample where the fragile markers have
//   already faded is old; a sample where all markers are still in their
//   original proportions is fresh. Crucially, that balance does not depend on
//   how much DNA was released to begin with, which is a quantity you almost
//   never know.
//
// THE PROCESS MODEL
//
//   For unit i (one water sample, or a group of samples collected together)
//   and marker j:
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
//   Notice why age and concentration do not get confused with one another.
//   C[i] raises or lowers all markers together by the same amount. t[i] is
//   multiplied by r[j], which differs between markers, so age tilts the
//   pattern instead of shifting it. Those are different fingerprints, which is
//   what lets the model tell them apart. If every marker decayed at the same
//   rate the two would be indistinguishable and the method would not work.
//
//   r and p come from a separate decay experiment, fitted by
//   decay_monophasic_ddPCR.stan, and are treated as known here.
//
// OBSERVATION MODEL
//
//   Droplet digital PCR splits one reaction into roughly twenty thousand
//   droplets and reports how many of them contained the target. So the data
//   are counts, not concentrations, and the model has to travel from a
//   concentration in seawater to an expected count.
//
//   Three things happen on the way, each adding to the log concentration:
//
//     1. Water samples from the same unit differ from one another, because
//        eDNA is patchy. That is eta, with standard deviation sigma_bio.
//
//     2. PCR wells from the same water sample differ from one another. That is
//        eps, with standard deviation sigma_tech.
//
//     3. A fixed volumetric conversion, log_offset, turns copies per litre of
//        seawater into expected copies per droplet, accounting for how much
//        water was filtered, any dilution, the elution volume and the size of
//        a droplet. It is computed in R and passed in.
//
//   Writing s for a water sample and r for a well:
//
//       log C[i,j,s]   = log C[i,j]   + eta[s,j]        eta ~ Normal(0, sigma_bio)
//       log C[i,j,r,s] = log C[i,j,s] + eps[r]          eps ~ Normal(0, sigma_tech)
//       omega          = log C[i,j,r,s] + log_offset    (log copies per droplet)
//
//       W ~ Binomial(U, 1 - exp(-exp(omega)))
//
//   W is the number of droplets that lit up and U the number the machine
//   accepted. The expression 1 - exp(-exp(omega)) is the chance that any one
//   droplet contains at least one copy, assuming copies scatter randomly among
//   droplets. In Stan it is written inv_cloglog(omega).
//
//   Wells that produced no positive droplets are kept and are informative: a
//   zero says the concentration was low, which is what an old marker looks
//   like. Discarding them would bias ages young.
//
// SHARED VARIANCE COMPONENTS
//
//   Both sigma_bio and sigma_tech are single numbers applied to every marker,
//   rather than one per marker. The reason is that noise on a log scale
//   interacts with the exponential in the likelihood: the concentration the
//   data imply shifts by an amount proportional to sigma squared (see
//   mean_correction below). A marker-specific sigma would shift each marker by
//   a different amount, which tilts the marker pattern -- and the marker
//   pattern is exactly what carries the age. Sharing them makes that shift one
//   constant, absorbed harmlessly by C[i], so measurement noise widens the age
//   estimate without moving it.
//
// NESTED REPLICATE STRUCTURE
//
//   Wells from one water sample are not independent observations of the unit.
//   Whatever makes one bottle differ from another is shared by all of its
//   wells and does not average away no matter how many wells you run. Keeping
//   the levels separate stops technical replication from buying precision it
//   has not earned, and lets each source of variation be reported on its own.
//
// DESIGNS WITHOUT BIOLOGICAL REPLICATION
//
//   Set use_bio = 0 when each unit is a single water sample, as in the carboy
//   experiment. The biological level then disappears completely and only
//   sigma_tech is estimated. N_bio and bio_idx are ignored.
// =============================================================================


data {
  // ---- Dimensions ---------------------------------------------------------
  int<lower=1> Nt;                       // number of units to date
  int<lower=1> Nloci;                    // number of markers

  // ---- Fixed from the decay experiment ------------------------------------
  vector[Nloci] r;                       // decay rate per marker, per hour, negative
  vector[Nloci] p;                       // starting level per marker, relative to
                                         //   marker 1, on the log scale; p[1] = 0

  // ---- Observations: one row per ddPCR well -------------------------------
  int<lower=0> N;                                   // number of wells
  array[N] int<lower=1, upper=Nt>    obs_i;         // which unit this well belongs to
  array[N] int<lower=1, upper=Nloci> obs_j;         // which marker this well measured
  array[N] int<lower=0> W;                          // droplets that were positive
  array[N] int<lower=1> U;                          // droplets the reader accepted
  vector[N] log_offset;                             // volumetric conversion, see header

  // ---- Replicate structure ------------------------------------------------
  int<lower=0, upper=1> use_bio;         // 1 = model between-sample variation
  int<lower=1> N_bio;                    // number of water samples
  array[N] int<lower=1> bio_idx;         // which water sample each well came from

  // ---- Prior hyperparameters ----------------------------------------------
  real C0_mean;                          // expected log concentration when shed
  real<lower=0> C0_sd;                   // how uncertain that expectation is
  real<lower=0> t_mean;                  // expected age in hours
  real<lower=0> t_sd;                    // how uncertain that expectation is
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
  // ---- Parameters of interest ---------------------------------------------
  vector<lower=0>[Nt] t;                 // age of each unit, hours, cannot be negative
  vector<lower=0>[Nt] C;                 // log concentration of each unit when shed

  // ---- Variance components ------------------------------------------------
  real<lower=0> sigma_tech;              // SD between wells of one water sample
  vector<lower=0>[use_bio ? 1 : 0] sigma_bio;
                                         // SD between water samples of one unit.
                                         //   Declared as a vector of length 1 (or 0
                                         //   when switched off) so that it can be
                                         //   removed cleanly; sb below unwraps it.

  // ---- Standardised random effects ----------------------------------------
  // These are held in standard units, mean 0 and SD 1, and are multiplied by
  // the standard deviations above to give the actual departures. Writing them
  // this way rather than drawing them directly at scale sigma leaves the
  // sampler an evenly shaped space to explore, which it handles far better.
  matrix[N_bio_eff, Nloci_eff] eta_raw;  // one per water sample and marker
  vector[N] eps_raw;                     // one per well
}


transformed parameters {
  matrix[Nt, Nloci] mu;                  // expected log concentration per unit and marker
  vector[N] omega;                       // log expected copies per droplet, per well

  // sigma_bio as a plain number, or zero when the biological level is off.
  real sb = use_bio ? sigma_bio[1] : 0;

  // Exponentiating a quantity that carries symmetric noise on the log scale
  // inflates its average by exp(sigma^2 / 2). Subtracting this term cancels
  // that, so C keeps its plain meaning as the concentration when shed. Because
  // both standard deviations are shared across markers this is a single
  // constant: it relabels C and cannot move t.
  real mean_correction = 0.5 * (square(sigma_tech) + square(sb));

  // The process model, evaluated for every unit and marker.
  for (i in 1:Nt) {
    for (j in 1:Nloci) {
      mu[i, j] = C[i] + p[j] + (r[j] * t[i]);
    }
  }

  // Walk out to each individual well.
  for (n in 1:N) {
    real level = mu[obs_i[n], obs_j[n]];        // where this unit and marker sit
    if (use_bio) {
      level += sb * eta_raw[bio_idx[n], obs_j[n]];   // this water sample's departure
    }
    omega[n] = level + log_offset[n]                 // convert to copies per droplet
               + sigma_tech * eps_raw[n]             // this well's departure
               - mean_correction;
  }
}


model {
  // ---- Priors -------------------------------------------------------------
  t ~ normal(t_mean, t_sd);              // truncated at 0 by the declaration above
  C ~ normal(C0_mean, C0_sd);            // half-normal

  // Variation between wells of one water sample.
  sigma_tech ~ normal(0, sigma_sd);      // half-normal: sigma_tech cannot be negative
  eps_raw    ~ std_normal();

  // Variation between water samples of one unit. Both lines do nothing at all
  // when use_bio = 0, because the containers are then empty.
  sigma_bio          ~ normal(0, sigma_sd);
  to_vector(eta_raw) ~ std_normal();     // to_vector flattens the matrix so the
                                         //   same prior applies to every entry

  // ---- Likelihood ---------------------------------------------------------
  W ~ binomial(U, inv_cloglog(omega));
}


generated quantities {
  vector[N] log_lik;                     // fit of each well, for model comparison
  array[N] int W_rep;                    // counts the fitted model would produce,
                                         //   for checking it can reproduce the data
  vector[N] p_detect;                    // chance this well detects anything at all
  real sigma_total = sqrt(square(sigma_tech) + square(sb));   // both sources combined

  for (n in 1:N) {
    real pi_n = inv_cloglog(omega[n]);   // chance a single droplet is positive

    log_lik[n]  = binomial_lpmf(W[n] | U[n], pi_n);
    W_rep[n]    = binomial_rng(U[n], pi_n);

    // Not the same as pi_n: this is the chance that at least one of the well's
    // U droplets comes up positive, which is the quantity comparable to a
    // detection rate reported from qPCR.
    p_detect[n] = 1 - exp(-U[n] * exp(omega[n]));
  }
}
