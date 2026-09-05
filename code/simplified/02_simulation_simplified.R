# =============================================================================
# 02_simulation_simplified.R
# -----------------------------------------------------------------------------
# Simulation study for conc_age_simplified.stan: the concentration model with
# observation error split into nested biological and technical components, both
# shared across markers.
#
# DATA-GENERATING PROCESS
#
#   1. Decay curve, per unit i and marker j
#          mu[i,j] = C0 + p[j] + r[j] * t[i]                (log copies / L)
#
#   2. Biological replicate s, a separate water sample from the same unit
#          log C[i,j,s] = mu[i,j] + eta[s,j]        eta ~ Normal(0, sigma_bio)
#
#   3. Detection, one draw per PCR replicate
#          z ~ Bernoulli(logit^-1(alpha + beta * log C[i,j,s]))
#
#   4. Measurement, only where the marker was detected
#          y = log C[i,j,s] + eps                   eps ~ Normal(0, sigma_tech)
#
#   Only z and y reach the model. Concentration panels show layer 4, the
#   measurements an instrument would return.
#
# WHAT THIS CHECKS
#   Whether sigma_bio and sigma_tech are separately recovered, whether the
#   detection function is recovered, and whether age posteriors are calibrated.
#   Non-detections are not discarded: each one tells the model the
#   concentration was below the limit of detection, which is what a
#   long-decayed marker looks like.
#
# FIGURES
#   plots/simplified/simulation_main.png        distinct vs observed decay rates
#   plots/simplified/simulation_supplement.png  noise, marker number, rate spread
#
# TABLES
#   outputs/simplified/simulation_recovery_main.csv   per-unit age posteriors
#   outputs/simplified/simulation_calibration.csv     bias, RMSE, coverage
#   outputs/simplified/simulation_sigma_recovery.csv  variance components
# =============================================================================

source(here::here("code", "simplified", "00_setup_simplified.R"))

set.seed(666)

# =============================================================================
# Simulation settings
# =============================================================================

Nt_main <- 30           # units, each with its own age
Nt_supp <- 20           # units per supplementary scenario
N_BIO   <- 3            # water samples per unit
K_REPS  <- 3            # PCR wells per water sample and marker
t_max   <- 50

C0_true <- 9.0          # baseline log concentration at t = 0 (copies/L)

r_actual <- unname(r_by_marker[locus_levels])
p_true   <- unname(p_by_marker[locus_levels])

# Noise levels measured in the field samples, defined in 00_setup_simplified.R.
# sigma_tech is the SD between wells of one water sample, sigma_bio the SD
# between water samples of one unit.
SIGMA_TECH_MAIN <- SIGMA_LEVELS$tech
SIGMA_BIO_MAIN  <- SIGMA_LEVELS$bio_field

# True detection function: half of replicates amplify at 10 copies per litre.
BETA_TRUE   <- 2
LOGC50_TRUE <- log(10)
ALPHA_TRUE  <- -BETA_TRUE * LOGC50_TRUE

message("Simulation settings:")
message("  ", Nt_main, " units x ", N_BIO, " water samples x ", K_REPS,
        " wells x ", length(r_actual), " markers")
message("  sigma_bio = ", SIGMA_BIO_MAIN, " | sigma_tech = ", SIGMA_TECH_MAIN)
message("  detection: 50% at ", round(exp(LOGC50_TRUE), 1),
        " copies/L, slope ", BETA_TRUE)

# =============================================================================
# Generative model
# =============================================================================

simulate_data <- function(r_true, p_local, sigma_bio, sigma_tech,
                          C0 = C0_true, Nt_local = Nt_main,
                          n_bio = N_BIO, n_tech = K_REPS,
                          t_max_local = t_max, seed_sim = 666) {
  set.seed(seed_sim)
  K_local <- length(r_true)
  t_local <- seq(0.5, t_max_local, length.out = Nt_local)

  obs <- expand.grid(obs_i = seq_len(Nt_local),
                     bio   = seq_len(n_bio),
                     obs_j = seq_len(K_local),
                     rep   = seq_len(n_tech))
  obs$t_true <- t_local[obs$obs_i]

  # a biological replicate is identified by unit and sample number
  obs$sample_id <- paste(obs$obs_i, obs$bio, sep = "_")

  # 1. decay curve
  obs$mu <- C0 + p_local[obs$obs_j] + r_true[obs$obs_j] * obs$t_true

  # 2. biological replicate: one deviate per (water sample, marker)
  eta <- matrix(rnorm(Nt_local * n_bio * K_local, 0, sigma_bio),
                nrow = Nt_local * n_bio, ncol = K_local,
                dimnames = list(unique(obs$sample_id), NULL))
  obs$logC_bio <- obs$mu + eta[cbind(match(obs$sample_id, rownames(eta)),
                                     obs$obs_j)]

  # 3. detection, one draw per well
  obs$detected <- rbinom(nrow(obs), 1,
                         plogis(ALPHA_TRUE + BETA_TRUE * obs$logC_bio))

  # 4. measurement, only where detected
  obs$logC <- ifelse(obs$detected == 1,
                     rnorm(nrow(obs), obs$logC_bio, sigma_tech), NA_real_)

  list(obs = obs, t_local = t_local,
       sigma_bio = sigma_bio, sigma_tech = sigma_tech)
}

simulate_and_fit <- function(r_true, sigma_bio = SIGMA_BIO_MAIN,
                             sigma_tech = SIGMA_TECH_MAIN, p_local = NULL,
                             seed_sim = 666, seed_stan = 42,
                             t_mean = 24, t_sd = 24, Nt_local = Nt_main) {
  K_local <- length(r_true)
  if (is.null(p_local)) {
    p_local <- c(p_true, rep(0, max(0, K_local - length(p_true))))[seq_len(K_local)]
  }
  sim <- simulate_data(r_true, p_local, sigma_bio, sigma_tech,
                       Nt_local = Nt_local, seed_sim = seed_sim)

  stan_data <- build_conc_stan_data(
    obs = sim$obs, Nt = length(sim$t_local), Nloci = K_local,
    r_vec = r_true, p_vec = p_local,
    C0_mean = C0_true, C0_sd = 1.8, t_mean = t_mean, t_sd = t_sd,
    bio = sim$obs$sample_id)

  # The biological deviations are indexed by marker, so a set of them can
  # partly mimic the marker pattern that identifies age: raising t[i] shifts
  # each marker by r[j] * dt, and eta[s,j] can absorb that shift at a cost of
  # well under one prior SD. The resulting ridge is long and nearly flat, so
  # deep trajectories buy almost nothing. Capping the treedepth at 10 rather
  # than 12 quarters the work per iteration for a similar effective sample
  # size; divergences are reported in outputs/simplified/simulation_calibration.csv
  # and are near zero at this adapt_delta.
  fit <- stan(file = stan_conc_age, data = stan_data,
              chains = 4, iter = 2000, warmup = 1000, seed = seed_stan,
              control = list(adapt_delta = 0.9, max_treedepth = 10),
              refresh = 200)

  list(fit = fit, sim = sim, t_local = sim$t_local,
       p_local = p_local, r_true = r_true)
}

# =============================================================================
# Summaries
# =============================================================================

recovery_table <- function(res, label) {
  ts <- summary(res$fit, pars = "t",
                probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
  data.frame(scenario = label, true = res$t_local,
             med = ts[, "50%"], mean = ts[, "mean"], sd = ts[, "sd"],
             lo95 = ts[, "2.5%"], hi95 = ts[, "97.5%"],
             lo50 = ts[, "25%"],  hi50 = ts[, "75%"],
             row.names = NULL)
}

# Variance components and detection function against the values simulated
sigma_table <- function(res, label) {
  pars <- c("sigma_bio", "sigma_tech", "beta", "logC50")
  ss   <- summary(res$fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  true <- c(res$sim$sigma_bio, res$sim$sigma_tech, BETA_TRUE, LOGC50_TRUE)
  data.frame(scenario  = label,
             component = pars,
             true      = true,
             estimate  = ss[, "mean"],
             lo95      = ss[, "2.5%"],
             hi95      = ss[, "97.5%"],
             covered   = true >= ss[, "2.5%"] & true <= ss[, "97.5%"],
             row.names = NULL)
}

calib_row <- function(rt, res, label) {
  div <- sum(sapply(rstan::get_sampler_params(res$fit, inc_warmup = FALSE),
                    function(x) sum(x[, "divergent__"])))
  data.frame(scenario = label, n = nrow(rt),
             non_detections = mean(res$sim$obs$detected == 0),
             bias = mean(rt$med - rt$true),
             mae  = mean(abs(rt$med - rt$true)),
             rmse = sqrt(mean((rt$med - rt$true)^2)),
             cover50 = mean(rt$true >= rt$lo50 & rt$true <= rt$hi50),
             cover95 = mean(rt$true >= rt$lo95 & rt$true <= rt$hi95),
             ci95_med = median(rt$hi95 - rt$lo95),
             divergences = div,
             max_rhat = max(summary(res$fit, pars = "t")$summary[, "Rhat"],
                            na.rm = TRUE))
}

regime_check <- function(sim_obs) {
  f <- here("data", "field_droplets.csv")
  if (!file.exists(f)) return(invisible(NULL))
  d <- read.csv(f, stringsAsFactors = FALSE)
  message("\nSimulated replicates vs field replicates:")
  message("  field     : non-detection fraction ", round(mean(d$positives == 0), 3))
  message("  simulated : non-detection fraction ",
          round(mean(sim_obs$detected == 0), 3))
  invisible(NULL)
}

# =============================================================================
# Figure panels
# =============================================================================

plot_concentrations <- function(res, component_labels, title_str) {
  K_local <- length(component_labels)
  res$sim$obs %>%
    mutate(Components = factor(obs_j, levels = seq_len(K_local),
                               labels = component_labels),
           y_plot = ifelse(detected == 1, logC, NA_real_)) %>%
    ggplot(aes(x = t_true, y = y_plot, color = Components)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
    labs(title = title_str, x = "Time elapsed",
         y = "Log simulated eDNA concentration", color = "Components") +
    theme_bw(base_size = 12) +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14))
}

plot_t_recovery <- function(res, title_str) {
  rt <- recovery_table(res, "")
  ggplot(rt, aes(x = true, y = med)) +
    geom_abline(slope = 1, linetype = "dashed", color = "red") +
    geom_linerange(aes(ymin = lo95, ymax = hi95),
                   colour = "#2166ac", linewidth = 0.6, alpha = 0.5) +
    geom_linerange(aes(ymin = lo50, ymax = hi50),
                   colour = "#2166ac", linewidth = 1.4) +
    geom_point(size = 2.5, colour = "#2166ac") +
    labs(title = title_str, x = "True time elapsed",
         y = "Estimated time elapsed") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14))
}

rate_labels_for <- function(r_vec, digits = 3) {
  paste0(seq_along(r_vec), " (λ = ",
         sprintf(paste0("%.", digits, "f"), r_vec), ")")
}

# =============================================================================
# Main figure
# =============================================================================

message("\nRun 1: distinct decay rates")
r_distinct <- c(-0.05, -0.10, -0.20, -0.40)
res1 <- simulate_and_fit(r_distinct, seed_sim = 66, seed_stan = 42)
message("  non-detection rate: ", round(mean(res1$sim$obs$detected == 0), 3))

message("\nRun 2: decay rates observed in the carboy experiment")
res2 <- simulate_and_fit(r_actual, seed_sim = 66, seed_stan = 123)
message("  non-detection rate: ", round(mean(res2$sim$obs$detected == 0), 3))
regime_check(res2$sim$obs)

simulation_main <- plot_grid(
  plot_concentrations(res2, rate_labels_for(r_actual, 3),
                      "Components (actual rates)"),
  plot_t_recovery(res2, "Model age estimation (actual rates)"),
  plot_concentrations(res1, rate_labels_for(r_distinct, 2),
                      "Components (distinct rates)"),
  plot_t_recovery(res1, "Model age estimation (distinct rates)"),
  labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 16)

ggsave(plot_path("simulation_main.png"), plot = simulation_main,
       width = 12, height = 9, dpi = 300)
message("Saved: plots/simplified/simulation_main.png")

rec_main <- bind_rows(recovery_table(res1, "distinct rates"),
                      recovery_table(res2, "actual rates"))
write.csv(rec_main, out_path("simulation_recovery_main.csv"), row.names = FALSE)

calib <- bind_rows(calib_row(recovery_table(res1, "distinct"), res1, "distinct rates"),
                   calib_row(recovery_table(res2, "actual"),   res2, "actual rates"))
sig   <- bind_rows(sigma_table(res1, "distinct rates"),
                   sigma_table(res2, "actual rates"))

message("\nObservation parameters, simulated versus estimated:")
print(as.data.frame(sig %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

# =============================================================================
# Supplementary figure
# =============================================================================
# Panels 1-3 hold the decay rates at their observed values and vary the amount
# of biological variation. Panels 4-6 hold noise fixed and vary the number of
# markers and the spread of their decay rates.

supp_scenarios <- list(
  list(r = r_actual, sb = SIGMA_LEVELS$bio_low, st = SIGMA_TECH_MAIN,
       lab = "low biological variation"),
  list(r = r_actual, sb = SIGMA_LEVELS$bio_field, st = SIGMA_TECH_MAIN,
       lab = "observed biological variation"),
  list(r = r_actual, sb = SIGMA_LEVELS$bio_high, st = SIGMA_TECH_MAIN,
       lab = "high biological variation"),
  list(r = c(-0.10, -0.12), sb = SIGMA_BIO_MAIN, st = SIGMA_TECH_MAIN,
       lab = "2 markers, similar rates"),
  list(r = c(-0.10, -0.20, -0.50), sb = SIGMA_BIO_MAIN, st = SIGMA_TECH_MAIN,
       lab = "3 markers"),
  list(r = c(-0.05, -0.10, -0.15, -0.20, -0.40, -0.80),
       sb = SIGMA_BIO_MAIN, st = SIGMA_TECH_MAIN, lab = "6 markers")
)

run_supp <- function(s, seed = 13) {
  res <- simulate_and_fit(s$r, sigma_bio = s$sb, sigma_tech = s$st,
                          seed_sim = seed, seed_stan = seed,
                          t_mean = 20, t_sd = 15, Nt_local = Nt_supp)
  labs <- rate_labels_for(s$r, 2)
  list(res    = res,
       top    = plot_concentrations(res, labs, paste0("Components (", s$lab, ")")),
       bottom = plot_t_recovery(res, paste0("Age estimation (", s$lab, ")")),
       calib  = calib_row(recovery_table(res, s$lab), res, s$lab),
       sigma  = sigma_table(res, s$lab))
}

supp_res <- lapply(supp_scenarios, run_supp)

simulation_supp <- plot_grid(
  plotlist = list(supp_res[[1]]$top,    supp_res[[2]]$top,    supp_res[[3]]$top,
                  supp_res[[1]]$bottom, supp_res[[2]]$bottom, supp_res[[3]]$bottom,
                  supp_res[[4]]$top,    supp_res[[5]]$top,    supp_res[[6]]$top,
                  supp_res[[4]]$bottom, supp_res[[5]]$bottom, supp_res[[6]]$bottom),
  ncol = 3, labels = paste0("(", letters[1:12], ")"), label_size = 16)

ggsave(plot_path("simulation_supplement.png"), plot = simulation_supp,
       width = 14, height = 16, dpi = 300)
message("Saved: plots/simplified/simulation_supplement.png")

calib <- bind_rows(calib, bind_rows(lapply(supp_res, `[[`, "calib")))
sig   <- bind_rows(sig,   bind_rows(lapply(supp_res, `[[`, "sigma")))

# =============================================================================
# Report
# =============================================================================

write.csv(calib, out_path("simulation_calibration.csv"),    row.names = FALSE)
write.csv(sig,   out_path("simulation_sigma_recovery.csv"), row.names = FALSE)

message("\nObservation parameters across all scenarios:")
print(as.data.frame(sig %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
message("95% interval coverage of the true values: ", round(mean(sig$covered), 3))

message("\nAge recovery across all scenarios:")
print(as.data.frame(calib %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

message("\nSaved: outputs/simplified/simulation_calibration.csv")
message("       outputs/simplified/simulation_sigma_recovery.csv")
message("       outputs/simplified/simulation_recovery_main.csv")
