# =============================================================================
# quickstart.R
# -----------------------------------------------------------------------------
# A five-minute tour of the method on a toy dataset.
#
# Open eDNA_quantitative_aging.Rproj, then run this file top to bottom. It
# makes a small simulated dataset, fits the simplified model to it, and checks
# the ages that come back against the ages that went in. Nothing here touches
# the manuscript data, so it is safe to edit and re-run.
#
# Runtime: about a minute, plus a one-off Stan compilation the first time.
# =============================================================================

source(here::here("code", "simplified", "age_functions.R"))


# --- 1. Make a toy dataset ---------------------------------------------------
# Twelve units to be aged. Each was sampled twice, and each water sample was
# run in three PCR wells for each of three markers, so 12 x 2 x 3 x 3 = 216
# replicates. The markers decay at -0.05, -0.12 and -0.25 per hour.

toy <- simulate_toy_edna(n_units = 12, n_bio = 2, n_reps = 3,
                         r = c(-0.05, -0.12, -0.25),
                         p = c(0, -0.5, -1.5),
                         seed = 1)

head(toy$obs, 10)

# Five columns is all the model needs:
#
#   obs_i     which unit is being aged
#   obs_j     which marker
#   sample    which water sample the replicate came from
#   detected  1 if the marker amplified, 0 if it did not
#   logC      log copies per litre, where it amplified
#
# Wells that did not amplify are kept. A non-detection is not missing data: it
# says the concentration was below the limit of detection, which is exactly
# what a marker that has been decaying for a long time looks like.

cat("replicates:", nrow(toy$obs),
    " non-detections:", sum(toy$obs$detected == 0), "\n")


# --- 2. Fit ------------------------------------------------------------------
# The decay rates r and the starting offsets p are treated as known. Here they
# come from the simulation; with real data they come from a decay experiment
# (see code/simplified/01_decay_simplified.R) or from published values.
#
# C0_mean and C0_sd are your prior for the log concentration at the moment of
# shedding, in copies per litre. t_mean and t_sd are your prior for age in
# hours, truncated at zero.

fit <- fit_edna_age(toy$obs,
                    r = toy$r, p = toy$p,
                    C0_mean = 9, C0_sd = 2,
                    t_mean  = 12, t_sd = 24,
                    chains = 4, iter = 2000, warmup = 1000, seed = 1)


# --- 3. Read the ages --------------------------------------------------------

ages <- age_table(fit, truth = toy$truth)
print(round(ages, 2))

cat("\nmean error:", round(mean(ages$error), 2), "hours\n")
cat("mean absolute error:", round(mean(abs(ages$error)), 2), "hours\n")
cat("true age inside the 95% interval:",
    round(mean(ages$t_true >= ages$lo95 & ages$t_true <= ages$hi95), 2), "\n")


# --- 4. Look at the observation model ----------------------------------------
# sigma_bio is the spread between water samples of one unit and sigma_tech the
# spread between wells of one sample, both on the log scale. logC50 is the
# concentration at which half of replicates amplify: the limit of detection the
# model inferred from your non-detections.

print(observation_table(fit))


# --- 5. Plot -----------------------------------------------------------------

plot_ages(ages)


# =============================================================================
# Using your own data
# =============================================================================
#
# Build the same five columns from your own table and call fit_edna_age().
# Three things to get right:
#
#   * one row per PCR replicate, non-detections included
#   * logC on a single volumetric scale. Everything here is log copies per
#     litre of water, so divide by the volume you filtered before you take the
#     log. If your instrument reports copies per mL, multiply by 1000 first.
#   * r and p from the same markers you measured, on the same scale. p[1] must
#     be exactly 0, since every offset is relative to the first marker.
#
# Leave the `sample` column out if you have no biological replication -- one
# water sample per unit, several wells per sample. The nested level is then
# switched off automatically and only sigma_tech is estimated.
#
# You need at least two markers with different decay rates. The further apart
# the rates, the sharper the age estimate; the supplementary simulation figure
# shows how much difference this makes.
# =============================================================================
