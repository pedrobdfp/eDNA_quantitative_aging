// =============================================================================
// decay_monophasic_ddPCR.stan
//
// Marker-specific first-order decay rates, estimated from a controlled
// experiment using ddPCR droplet counts.
// -----------------------------------------------------------------------------
//
// PURPOSE
//
//   The age model (droplet_age_ddPCR.stan) needs to know two things about each
//   marker before it can date anything: how quickly the marker disappears, and
//   where it starts relative to the others. This file measures both, from an
//   experiment in which seawater was held in sealed carboys and subsampled
//   repeatedly over time.
//
//   Here the elapsed time is known and the decay rate is estimated. In the age
//   model it is the other way round: the rate is known and the time is
//   estimated. The same observation model is used in both, so the rates are
//   measured under the same likelihood that later consumes them.
//
// THE PROCESS MODEL
//
//   Each carboy-and-marker combination is called a series and has its own
//   starting level. Within a series, concentration falls exponentially, which
//   on the log scale is a straight line:
//
//       log C[n] = C0[series[n]] - lambda[marker[n]] * time[n]
//
//   in log copies per litre, where
//
//       C0      the log concentration at the start of the experiment, one per
//               series, estimated as the intercept of that line
//       lambda  the decay rate of a marker, per hour, estimated. It is held
//               positive and subtracted, so a larger lambda means faster loss.
//       time    hours since the carboys were filled, known
//
//   The samples taken at the start are ordinary observations with time = 0.
//   They inform C0 through the same likelihood as every other observation and
//   are given no special treatment, so C0 really is just the intercept of the
//   fitted line.
//
// OBSERVATION MODEL
//
//       omega[n] = log C[n] + log_offset[n] + sigma_obs * eps[n] - sigma_obs^2/2
//       W[n] ~ Binomial(U[n], 1 - exp(-exp(omega[n])))
//
//   W is the number of droplets that came up positive out of U accepted. The
//   expression 1 - exp(-exp(omega)), written inv_cloglog in Stan, is the chance
//   that any one droplet holds at least one copy.
//
//   log_offset converts copies per litre of seawater into expected copies per
//   droplet, given how much water was filtered, any dilution, the elution
//   volume and the size of a droplet. It is computed in R by the same function
//   the field analysis uses, so both parts of the study share one definition.
//
//   sigma_obs is the spread between PCR wells measuring the same sample. The
//   final term corrects for the fact that exponentiating symmetric noise on a
//   log scale would otherwise inflate the average; because sigma_obs is shared
//   across markers this correction is one constant and cannot bend the decay
//   rates.
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
  int<lower=1> N;                        // number of PCR wells
  int<lower=1> N_series;                 // carboy x marker combinations
  int<lower=1> N_marker;                 // number of markers

  // ---- Observations: one row per ddPCR well -------------------------------
  array[N] int<lower=1, upper=N_series> series;   // which series this well belongs to
  array[N] int<lower=1, upper=N_marker> marker;   // which marker it measured
  vector<lower=0>[N] time;                        // hours since the carboys were filled
  vector[N] log_offset;                           // volumetric conversion, see header
  array[N] int<lower=0> W;                        // droplets that were positive
  array[N] int<lower=1> U;                        // droplets the reader accepted

  // ---- Marker of each series, used for the per-marker summaries -----------
  array[N_series] int<lower=1, upper=N_marker> series_marker;

  // ---- Prior hyperparameters ----------------------------------------------
  real<lower=0> sigma_obs_sd;            // scale of the prior on well-to-well spread
  real<lower=0> lambda_sd;               // scale of the prior on the decay rates
  real C0_mean;                          // expected starting log concentration
  real<lower=0> C0_sd;                   // how uncertain that expectation is
}


parameters {
  vector[N_series] C0;                   // starting log concentration of each series
  vector<lower=0>[N_marker] lambda;      // decay rate of each marker, per hour
  real<lower=0> sigma_obs;               // spread between wells of one sample

  // Departures of individual wells, held in standard units and multiplied by
  // sigma_obs below. Written this way rather than drawn directly at scale
  // sigma_obs, which gives the sampler an evenly shaped space to explore.
  vector[N] eps_raw;
}


transformed parameters {
  vector[N] omega;                       // log expected copies per droplet, per well

  {
    // Cancels the inflation that exponentiating log-scale noise would cause,
    // so C0 keeps its plain meaning. One constant, shared by every marker.
    real correction = 0.5 * square(sigma_obs);

    for (n in 1:N) {
      omega[n] = C0[series[n]] - lambda[marker[n]] * time[n]   // the decay line
                 + log_offset[n]                               // to copies per droplet
                 + sigma_obs * eps_raw[n]                      // this well's departure
                 - correction;
    }
  }
}


model {
  // ---- Priors -------------------------------------------------------------
  C0        ~ normal(C0_mean, C0_sd);
  lambda    ~ normal(0, lambda_sd);      // half-normal: lambda cannot be negative
  sigma_obs ~ normal(0, sigma_obs_sd);   // half-normal
  eps_raw   ~ std_normal();

  // ---- Likelihood ---------------------------------------------------------
  W ~ binomial(U, inv_cloglog(omega));
}


generated quantities {
  vector[N_marker] r = -lambda;          // decay rates in the age model's convention
  vector[N_marker] half_life;            // hours for a marker to halve
  vector[N_marker] C0_bar;               // average starting level of each marker
  vector[N_marker] p_log;                // starting level relative to marker 1
  vector[N] log_lik;                     // fit of each well, for model comparison

  for (k in 1:N_marker) {
    half_life[k] = log(2) / lambda[k];

    // Average the intercepts of every series belonging to this marker, that is,
    // across carboys.
    real total = 0;
    int  n_k   = 0;
    for (s in 1:N_series) {
      if (series_marker[s] == k) { total += C0[s]; n_k += 1; }
    }
    C0_bar[k] = n_k == 0 ? negative_infinity() : total / n_k;
  }

  // Offsets are relative, so the first marker is the reference and is zero by
  // construction.
  for (k in 1:N_marker) p_log[k] = C0_bar[k] - C0_bar[1];

  for (n in 1:N) log_lik[n] = binomial_lpmf(W[n] | U[n], inv_cloglog(omega[n]));
}
