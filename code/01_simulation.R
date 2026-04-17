# =============================================================================
# 01_simulation.R
# -----------------------------------------------------------------------------
# Simulation study for the joint age-estimation model.
#
# Produces two figures:
#   (1) Main: 4-panel grid comparing "actual" (observed in carboys) versus
#       "distinct" (artificially spread out) decay rates.
#       -> plots/simulation_main.png
#   (2) Supplement: 12-panel grid varying the number of loci, the spread of
#       decay rates, and the measurement noise sigma.
#       -> plots/simulation_supplement.png
#
# The model used here is conc_age_sim.stan (tight priors appropriate for
# simulation). See 02_decay.R onwards for the field variant.
# =============================================================================

source(here::here("code", "00_setup.R"))

# =============================================================================
# Shared simulation parameters
# =============================================================================

set.seed(666)

Nt_main   <- 50                       # number of "samples" (time points)
Nloci     <- 4                        # number of markers
K_reps    <- 3                        # PCR replicates per sample * marker

# Marker proportions at t = 0 (matches the carboy marker set)
p_true    <- c(1, 0.7957, 0.4096, 0.1029)

t_true    <- seq(0.5, 50, length.out = Nt_main)
C0_true   <- 6
alpha_true <- -3.0
beta_true  <-  4.0
sigma_true <-  0.5

# =============================================================================
# Helpers
# =============================================================================

# ---- Simulate one dataset + fit the model ----------------------------------

simulate_and_fit <- function(r_true,
                             seed_sim  = 666,
                             seed_stan = 42,
                             t_mean    = 24,
                             t_sd      = 24,
                             Nt_local  = Nt_main,
                             sigma     = sigma_true,
                             p_local   = NULL) {

  K_local <- length(r_true)
  if (is.null(p_local)) {
    # Pad / truncate p_true to match the number of markers in r_true
    p_local <- c(p_true, rep(0, max(0, K_local - length(p_true))))[seq_len(K_local)]
  }
  t_local <- seq(0.5, 50, length.out = Nt_local)

  set.seed(seed_sim)

  # True latent log-concentrations
  mu_true <- sweep(outer(t_local, r_true), 2, p_local + C0_true, "+")

  # Build observation list: all (sample, marker, replicate) combinations
  obs_list <- expand.grid(i = seq_len(Nt_local),
                          j = seq_len(K_local),
                          k = seq_len(K_reps))
  obs_list$mu_ij <- mu_true[cbind(obs_list$i, obs_list$j)]
  obs_list$p     <- plogis(alpha_true + beta_true * obs_list$mu_ij)
  obs_list$z     <- rbinom(nrow(obs_list), 1, prob = obs_list$p)
  obs_list$y     <- ifelse(obs_list$z == 1,
                           rnorm(nrow(obs_list), obs_list$mu_ij, sigma),
                           NA_real_)
  obs_list$t_true_val <- t_local[obs_list$i]

  det_rows <- obs_list[obs_list$z == 1, ]

  stan_data <- list(
    Nt    = Nt_local,
    Nloci = K_local,
    r     = r_true,
    p     = p_local,

    N_obs = nrow(obs_list),
    obs_i = obs_list$i,
    obs_j = obs_list$j,
    z     = obs_list$z,

    N_y   = nrow(det_rows),
    y_i   = det_rows$i,
    y_j   = det_rows$j,
    y_obs = det_rows$y,

    C0_mean = 4.4,
    C0_sd   = 1.8,
    t_mean  = t_mean,
    t_sd    = t_sd
  )

  fit <- stan(
    file    = stan_conc_sim,
    data    = stan_data,
    chains  = 4,
    iter    = 2000,
    warmup  = 1000,
    seed    = seed_stan,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 200
  )

  list(fit       = fit,
       obs_list  = obs_list,
       det_rows  = det_rows,
       t_local   = t_local,
       stan_data = stan_data)
}

# ---- Concentration scatter (observed detections, by marker) -----------------

plot_concentrations <- function(obs_list, component_labels, title_str) {
  K_local <- length(component_labels)
  obs_list %>%
    mutate(
      Components = factor(j, levels = seq_len(K_local), labels = component_labels),
      y_plot     = ifelse(z == 1, y, NA_real_)
    ) %>%
    ggplot(aes(x = t_true_val, y = y_plot, color = Components)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
    labs(
      title = title_str,
      x     = "Time elapsed",
      y     = "Log simulated eDNA concentration",
      color = "Components"
    ) +
    theme_bw(base_size = 12) +
    scale_color_brewer(palette = "Set2") +
    theme(
      legend.position = "right",
      plot.title      = element_text(hjust = 0.5, size = 16),
      axis.title      = element_text(size = 14)
    )
}

# ---- Age-recovery plot (estimated t vs true t) ------------------------------

plot_t_recovery <- function(fit, title_str, t_vec) {
  t_sum <- summary(fit, pars = "t",
                   probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  tibble(t = t_vec) %>%
    bind_cols(as.data.frame(t_sum)) %>%
    ggplot(aes(x = t, y = `50%`)) +
    geom_abline(slope = 1, linetype = "dashed", color = "red") +
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`),
                   colour = "#2166ac", linewidth = 0.6, alpha = 0.5) +
    geom_linerange(aes(ymin = `25%`, ymax = `75%`),
                   colour = "#2166ac", linewidth = 1.4) +
    geom_point(size = 2.5, colour = "#2166ac") +
    labs(title = title_str,
         x     = "True time elapsed",
         y     = "Estimated time elapsed") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = 16),
          axis.title      = element_text(size = 14))
}

# =============================================================================
# Main figure: distinct vs actual decay rates
# =============================================================================

# ---- Run 1: distinct rates ---------------------------------------------------

message("\n=== Simulation run 1: distinct decay rates ===")

r_distinct <- c(-0.05, -0.10, -0.20, -0.40)
res1       <- simulate_and_fit(r_true = r_distinct, seed_sim = 66, seed_stan = 42)

message("Detection rate (distinct):", round(mean(res1$obs_list$z), 3))

labels_distinct <- c("1 (\u03bb = -0.05)",
                     "2 (\u03bb = -0.1)",
                     "3 (\u03bb = -0.2)",
                     "4 (\u03bb = -0.4)")

p_conc1 <- plot_concentrations(res1$obs_list, labels_distinct,
                               "Components (distinct rates)")
p_trec1 <- plot_t_recovery(res1$fit,
                           "Model age estimation (distinct rates)",
                           t_vec = res1$t_local)
p_trec1
# ---- Run 2: actual (carboy-observed) rates ----------------------------------

message("\n=== Simulation run 2: actual decay rates ===")

r_actual <- c(-0.114, -0.165, -0.167, -0.192)
res2     <- simulate_and_fit(r_true = r_actual, seed_sim = 66, seed_stan = 123)

message("Detection rate (actual):", round(mean(res2$obs_list$z), 3))

labels_actual <- c("1 (\u03bb = -0.114)",
                   "2 (\u03bb = -0.165)",
                   "3 (\u03bb = -0.167)",
                   "4 (\u03bb = -0.192)")

p_conc2 <- plot_concentrations(res2$obs_list, labels_actual,
                               "Components (actual rates)")
p_trec2 <- plot_t_recovery(res2$fit,
                           "Model age estimation (actual rates)",
                           t_vec = res2$t_local)
p_trec2 
# ---- Assemble 4-panel grid --------------------------------------------------

simulation_main <- plot_grid(
  p_conc2, p_trec2,
  p_conc1, p_trec1,
  labels     = c("(a)", "(b)", "(c)", "(d)"),
  label_size = 16
)

ggsave(here("plots", "simulation_main.png"),
       plot   = simulation_main,
       width  = 12, height = 9, dpi = 300)

message("Saved: plots/simulation_main.png")

# ---- Summary statistics for the text ----------------------------------------

for (run_label in c("distinct rates", "actual rates")) {
  res_obj <- if (run_label == "distinct rates") res1 else res2
  t_sum   <- summary(res_obj$fit, pars = "t",
                     probs = c(0.025, 0.5, 0.975))$summary
  recovery <- data.frame(
    true  = res_obj$t_local,
    est   = t_sum[, "50%"],
    lo95  = t_sum[, "2.5%"],
    hi95  = t_sum[, "97.5%"]
  )
  recovery$covered <- recovery$true >= recovery$lo95 &
                      recovery$true <= recovery$hi95
  message("\n--- t recovery: ", run_label, " ---")
  message("95% CI coverage: ", round(mean(recovery$covered), 2))
  print(summary(res_obj$fit,
                pars  = c("beta", "sigma", "alpha"),
                probs = c(0.025, 0.5, 0.975))$summary)
}

# =============================================================================
# Supplementary figure: varying N markers, rate spread, and sigma
# =============================================================================

# Use 25 time points per panel for speed
Nt_supp <- 25

supp_scenarios <- list(
  list(r = c(-0.10, -0.20),                              sigma = 0.2, seed = 13),
  list(r = c(-0.10, -0.20),                              sigma = 0.5, seed = 13),
  list(r = c(-0.10, -0.20),                              sigma = 1.0, seed = 13),
  list(r = c(-0.10, -0.12),                              sigma = 0.5, seed = 13),
  list(r = c(-0.10, -0.20, -0.50),                       sigma = 0.5, seed = 13),
  list(r = c(-0.05, -0.10, -0.15, -0.20, -0.40, -0.80), sigma = 0.5, seed = 13)
)

run_supp_scenario <- function(s, t_mean = 20, t_sd = 15) {
  res <- simulate_and_fit(
    r_true    = s$r,
    seed_sim  = s$seed,
    seed_stan = s$seed,
    t_mean    = t_mean,
    t_sd      = t_sd,
    Nt_local  = Nt_supp,
    sigma     = s$sigma
  )

  K_local     <- length(s$r)
  rate_labels <- paste0(seq_len(K_local), " (\u03bb = ",
                        sub("^0\\.", ".", sprintf("%.2f", abs(s$r))), ")")
  title_suffix <- paste0("N = ", K_local, ", \u03c3 = ", s$sigma)

  list(
    top    = plot_concentrations(res$obs_list, rate_labels,
                                 paste0("Components (", title_suffix, ")")),
    bottom = plot_t_recovery(res$fit,
                             paste0("Model age estimation (", title_suffix, ")"),
                             t_vec = res$t_local)
  )
}

supp_res <- lapply(supp_scenarios, run_supp_scenario)

# Arrange: row 1 = scenarios 1-3 concentrations, row 2 = their recoveries,
# row 3 = scenarios 4-6 concentrations, row 4 = their recoveries
panels_ordered <- list(
  supp_res[[1]]$top,    supp_res[[2]]$top,    supp_res[[3]]$top,
  supp_res[[1]]$bottom, supp_res[[2]]$bottom, supp_res[[3]]$bottom,
  supp_res[[4]]$top,    supp_res[[5]]$top,    supp_res[[6]]$top,
  supp_res[[4]]$bottom, supp_res[[5]]$bottom, supp_res[[6]]$bottom
)

simulation_supp <- plot_grid(
  plotlist   = panels_ordered,
  ncol       = 3,
  labels     = paste0("(", letters[1:12], ")"),
  label_size = 16
)
simulation_supp
ggsave(here("plots", "simulation_supplement.png"),
       plot = simulation_supp, width = 14, height = 16, dpi = 300)

message("Saved: plots/simulation_supplement.png")
