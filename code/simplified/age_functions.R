# =============================================================================
# age_functions.R
# -----------------------------------------------------------------------------
# A small set of functions for applying the simplified model to your own data.
# Four steps: make or load a dataset, fit, summarise, plot.
#
#     source(here::here("code", "simplified", "age_functions.R"))
#
#     toy <- simulate_toy_edna()
#     fit <- fit_edna_age(toy$obs, r = toy$r, p = toy$p)
#     ages <- age_table(fit, truth = toy$truth)
#     plot_ages(ages)
#
# quickstart.R at the repository root runs exactly that.
# =============================================================================

source(here::here("code", "simplified", "00_setup_simplified.R"))

# -----------------------------------------------------------------------------
# simulate_toy_edna()
# -----------------------------------------------------------------------------
# A small dataset in the format the model expects, for trying the workflow out
# before you have your own data ready.
#
#   n_units   how many things you want an age for. Each is a water body, a
#             station visit, or a pooled grab -- whatever you are dating.
#   n_bio     water samples collected per unit. Set to 1 if you took one.
#   n_reps    PCR wells run per water sample and marker.
#   r         decay rate of each marker, per hour, negative.
#   p         each marker's log concentration at t = 0 relative to the first
#             marker, so p[1] is 0. A marker starting at half the concentration
#             of the reference has p = log(0.5) = -0.69.
#   C0        log concentration of the reference marker at t = 0, copies/L.
#   sigma_bio SD between water samples of one unit, on the log scale.
#   sigma_tech SD between wells of one water sample.
#
# Returns the observation table, the r and p it used, and the true ages.

simulate_toy_edna <- function(n_units   = 12,
                              n_bio     = 2,
                              n_reps    = 3,
                              r         = c(-0.05, -0.12, -0.25),
                              p         = c(0, -0.5, -1.5),
                              C0        = 9,
                              t_max     = 40,
                              sigma_bio = 0.4,
                              sigma_tech = 0.15,
                              logC50    = log(10),
                              slope     = 2,
                              seed      = 1) {
  stopifnot(length(r) == length(p), p[1] == 0)
  set.seed(seed)

  t_true <- seq(1, t_max, length.out = n_units)

  obs <- expand.grid(unit   = seq_len(n_units),
                     bio    = seq_len(n_bio),
                     marker = seq_along(r),
                     rep    = seq_len(n_reps))
  obs$sample <- paste0("unit", obs$unit, "_s", obs$bio)

  # true log concentration of each water sample
  mu  <- C0 + p[obs$marker] + r[obs$marker] * t_true[obs$unit]
  eta <- matrix(rnorm(n_units * n_bio * length(r), 0, sigma_bio),
                nrow = n_units * n_bio, ncol = length(r),
                dimnames = list(unique(obs$sample), NULL))
  level <- mu + eta[cbind(match(obs$sample, rownames(eta)), obs$marker)]

  # detection, then measurement where detected
  obs$detected <- rbinom(nrow(obs), 1, plogis(slope * (level - logC50)))
  obs$logC     <- ifelse(obs$detected == 1,
                         rnorm(nrow(obs), level, sigma_tech), NA_real_)

  obs$obs_i <- obs$unit
  obs$obs_j <- obs$marker

  list(obs   = obs[, c("obs_i", "obs_j", "sample", "rep", "detected", "logC")],
       r     = r,
       p     = p,
       truth = data.frame(obs_i = seq_len(n_units), t_true = t_true))
}

# -----------------------------------------------------------------------------
# fit_edna_age()
# -----------------------------------------------------------------------------
# Fit the simplified model.
#
#   obs   one row per PCR replicate, with columns
#           obs_i     which unit is being aged (integer, 1..Nt)
#           obs_j     which marker (integer, matching the order of r and p)
#           detected  0 or 1
#           logC      natural log of copies per litre; ignored where detected
#                     is 0, so NA there is fine
#           sample    optional, which water sample the replicate came from.
#                     Supply it when a unit contains several water samples;
#                     leave it out when replication is technical only.
#   r, p  decay rates and t = 0 offsets, one per marker. Estimate them from a
#         decay experiment with 01_decay_simplified.R, or supply published
#         values.
#   C0_mean, C0_sd  prior on the log concentration at t = 0, copies/L. Use what
#         you know about your system; a previous survey at the same site is
#         ideal.
#   t_mean, t_sd    prior on age in hours, truncated at zero.
#   control         passed to Stan. The default caps the treedepth at 10,
#                   because the per-marker biological deviations open a long,
#                   nearly flat ridge against age and deeper trajectories cost
#                   a great deal for very little. Raise it if you see poor
#                   effective sample sizes rather than warnings.
#
# Returns the stanfit together with the data list it was given.

fit_edna_age <- function(obs, r, p,
                         C0_mean = 9, C0_sd = 2,
                         t_mean  = 12, t_sd = 24,
                         sample_col = "sample",
                         chains = 4, iter = 2000, warmup = 1000,
                         seed = 1,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         ...) {
  stopifnot(length(r) == length(p))
  bio <- if (!is.null(sample_col) && sample_col %in% names(obs))
           obs[[sample_col]] else NULL

  stan_data <- build_conc_stan_data(
    obs     = obs,
    Nt      = max(obs$obs_i),
    Nloci   = length(r),
    r_vec   = r,
    p_vec   = p,
    C0_mean = C0_mean,
    C0_sd   = C0_sd,
    t_mean  = t_mean,
    t_sd    = t_sd,
    bio     = bio)

  fit <- stan(file = stan_conc_age, data = stan_data,
              chains = chains, iter = iter, warmup = warmup, seed = seed,
              control = control, ...)

  structure(list(fit = fit, data = stan_data, obs = obs),
            class = "edna_age_fit")
}

# -----------------------------------------------------------------------------
# age_table()
# -----------------------------------------------------------------------------
# Posterior age of each unit: median, 50% and 95% intervals. Pass `truth` (a
# data frame with obs_i and t_true) when you know the answer, as in a
# simulation, and the error is added.

age_table <- function(x, truth = NULL) {
  fit <- if (inherits(x, "edna_age_fit")) x$fit else x
  ts  <- summary(fit, pars = "t",
                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  out <- data.frame(obs_i = seq_len(nrow(ts)),
                    t_med = ts[, "50%"],
                    lo95  = ts[, "2.5%"],
                    lo50  = ts[, "25%"],
                    hi50  = ts[, "75%"],
                    hi95  = ts[, "97.5%"],
                    row.names = NULL)

  if (!is.null(truth)) {
    out <- merge(out, truth, by = "obs_i", all.x = TRUE)
    out$error <- out$t_med - out$t_true
  }
  out
}

# -----------------------------------------------------------------------------
# plot_ages()
# -----------------------------------------------------------------------------
# If the table has a t_true column, estimates are plotted against it with a
# one-to-one line. Otherwise each unit gets an interval.

plot_ages <- function(ages) {
  if ("t_true" %in% names(ages)) {
    ggplot(ages, aes(x = t_true, y = t_med)) +
      geom_abline(slope = 1, linetype = "dashed", colour = "red") +
      geom_linerange(aes(ymin = lo95, ymax = hi95),
                     colour = "#2166ac", linewidth = 0.6, alpha = 0.5) +
      geom_linerange(aes(ymin = lo50, ymax = hi50),
                     colour = "#2166ac", linewidth = 1.4) +
      geom_point(size = 2.5, colour = "#2166ac") +
      labs(x = "True age (hours)", y = "Estimated age (hours)") +
      theme_bw(base_size = 12)
  } else {
    ggplot(ages, aes(x = factor(obs_i), y = t_med)) +
      geom_linerange(aes(ymin = lo95, ymax = hi95),
                     colour = "#2166ac", linewidth = 0.6, alpha = 0.5) +
      geom_linerange(aes(ymin = lo50, ymax = hi50),
                     colour = "#2166ac", linewidth = 1.4) +
      geom_point(size = 2.5, colour = "#2166ac") +
      coord_flip() +
      labs(x = "Unit", y = "Estimated age (hours)") +
      theme_bw(base_size = 12)
  }
}

# -----------------------------------------------------------------------------
# observation_table()
# -----------------------------------------------------------------------------
# The fitted observation model: the two noise components and the detection
# function. logC50 is reported both on the log scale and as a concentration.

observation_table <- function(x) {
  fit  <- if (inherits(x, "edna_age_fit")) x$fit else x
  pars <- c("sigma_tech", "alpha", "beta", "logC50")
  if ("sigma_bio[1]" %in% names(fit)) pars <- c("sigma_bio", pars)

  ss <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  out <- data.frame(parameter = rownames(ss),
                    mean = ss[, "mean"], sd = ss[, "sd"],
                    lo95 = ss[, "2.5%"], hi95 = ss[, "97.5%"],
                    Rhat = ss[, "Rhat"], row.names = NULL)
  i <- match("logC50", out$parameter)
  if (!is.na(i)) {
    out$copies_per_L <- NA_real_
    out$copies_per_L[i] <- exp(out$mean[i])
  }
  out
}
