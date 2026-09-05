# =============================================================================
# 01_decay_simplified.R
# -----------------------------------------------------------------------------
# Carboy decay experiment under the simplified concentration model. Produces
# the marker-specific decay rates and t = 0 offsets that the field analyses
# treat as known, and validates the age model by leave-one-carboy-out
# cross-validation.
#
#   Part A  Fits all three carboys jointly and writes r and p to
#           outputs/simplified/carboy_rates.csv.
#
#   Part B  For each carboy in turn, estimates decay rates from the other two
#           and uses them to recover the elapsed time of the held-out carboy's
#           samples. Held-out ages are compared against the known sampling
#           times in plots/simplified/carboy_loo_grid.png.
#
# Both steps use the concentration observation model, so the decay rates are
# estimated under the likelihood that later consumes them. r is reported as a
# negative rate; the Stan model parameterises the positive magnitude lambda and
# subtracts it.
#
# The carboys have technical replication only -- one water sample per carboy
# per timepoint, several PCR wells per sample -- so the biological level of the
# age model is switched off here and only sigma_tech is estimated.
# =============================================================================

source(here::here("code", "simplified", "00_setup_simplified.R"))

# --- Input -------------------------------------------------------------------

decay_raw <- read.csv(data_path("Final_decay_ddPCR_datasheet.csv"),
                      stringsAsFactors = FALSE)

# --- Settings ----------------------------------------------------------------

components_use  <- "DNA"
hours_max_decay <- 30
hours_max_time  <- 30
t_mean_prior    <- 12
t_sd_prior      <- 24

locus_levels_master <- c("DNA_Cytb", "DNA_16S", "DNA_Dloop", "DNA_Bridge")

locus_to_marker <- c(DNA_Cytb   = "cytb",
                     DNA_16S    = "Tt_16S",
                     DNA_Dloop  = "Tt_DLL1",
                     DNA_Bridge = "Tt_longFrag")

# =============================================================================
# Helpers
# =============================================================================

# One row per PCR replicate, with the two columns the observation model needs:
# a detection flag and, where detected, log copies per litre.
prepare_decay <- function(dat, locus_levels_ord) {
  dat %>%
    mutate(locus      = factor(paste(Component, Marker, sep = "_"),
                               levels = locus_levels_ord),
           k_idx      = as.integer(locus),
           Hours_base = as.numeric(Hours_base),
           copies_mL  = as.numeric(copies_mL),
           detected   = as.integer(!is.na(copies_mL) & copies_mL > 0),
           logC       = ifelse(detected == 1, log(copies_mL * 1000), NA_real_)) %>%
    filter(!is.na(k_idx), !is.na(Hours_base))
}

build_decay_data <- function(dat, locus_levels_ord) {
  d <- prepare_decay(dat, locus_levels_ord) %>%
    mutate(ik_key = paste(Carboy, k_idx, sep = "|"))

  ik_lev  <- unique(d$ik_key)
  d$ik    <- match(d$ik_key, ik_lev)
  ik_to_k <- d %>% distinct(ik, k_idx) %>% arrange(ik) %>% pull(k_idx)

  det <- which(d$detected == 1)

  list(
    N          = nrow(d),
    N_ik       = length(ik_lev),
    N_k        = length(locus_levels_ord),
    ik         = as.integer(d$ik),
    k          = as.integer(d$k_idx),
    time       = as.numeric(d$Hours_base),
    z          = as.integer(d$detected),
    N_y        = length(det),
    y_row      = as.integer(det),
    y_obs      = as.numeric(d$logC[det]),
    ik_to_k    = as.integer(ik_to_k),
    C0_mean    = 0,
    C0_sd      = 4,
    logC50_mean = detection_priors()$logC50_mean,
    logC50_sd   = detection_priors()$logC50_sd,
    beta_mean   = detection_priors()$beta_mean,
    beta_sd     = detection_priors()$beta_sd,
    sigma_sd   = SIGMA_SD,
    lambda_sd  = 1,
    N_time_sim = length(unique(d$Hours_base)),
    time_sim   = sort(unique(as.numeric(d$Hours_base)))
  )
}

# (p, C0) from the t = 0 samples, as log offsets on the copies/L scale
compute_p_and_C0_L <- function(dat_t0, locus_levels_ord) {
  t0_summ <- prepare_decay(dat_t0, locus_levels_ord) %>%
    filter(detected == 1) %>%
    group_by(locus) %>%
    summarise(mean_logC = mean(logC), sd_logC = sd(logC), n = n(),
              .groups = "drop") %>%
    arrange(locus)

  ref <- t0_summ %>% filter(locus == locus_levels_ord[1])

  list(p       = t0_summ$mean_logC - ref$mean_logC,   # p[ref] == 0
       C0_mean = ref$mean_logC,
       C0_sd   = max(ref$sd_logC, 1.5, na.rm = TRUE),
       table   = t0_summ)
}

# Concentration observations for a held-out carboy, ready for the age model
build_age_conc_data <- function(dat_test, locus_levels_ord) {
  d <- prepare_decay(dat_test, locus_levels_ord) %>%
    mutate(obs_j = k_idx)

  tp_map <- d %>%
    group_by(Timepoint) %>%
    summarise(True_Hours = median(Hours_base, na.rm = TRUE), .groups = "drop") %>%
    arrange(True_Hours) %>%
    mutate(obs_i = row_number())

  d <- d %>% left_join(tp_map %>% dplyr::select(Timepoint, obs_i),
                       by = "Timepoint")

  list(obs = d, tp_map = tp_map)
}

# =============================================================================
# Part A. Decay fit on all three carboys
# =============================================================================

message("\n=== Part A: decay fit, concentration observation model ===")

dat_all <- decay_raw %>%
  filter(!Control,
         Component %in% components_use,
         Hours_base <= hours_max_decay,
         Carboy %in% c(1, 2, 3))

decay_sd_all <- build_decay_data(dat_all, locus_levels_master)

message("Observations: ", decay_sd_all$N,
        " | detected: ", decay_sd_all$N_y,
        " | non-detections: ", decay_sd_all$N - decay_sd_all$N_y)

fit_decay_all <- stan(
  file    = stan_decay,
  data    = decay_sd_all,
  init    = function() list(C0     = rep(9,    decay_sd_all$N_ik),
                            lambda = rep(0.15, decay_sd_all$N_k)),
  chains  = 4, iter = 3000, warmup = 1000, seed = 99,
  control = list(adapt_delta = 0.95), refresh = 300
)

lambda_summ  <- summary(fit_decay_all, pars = "lambda",
                        probs = c(0.025, 0.5, 0.975))$summary
r_raw        <- -as.numeric(lambda_summ[, "mean"])
names(r_raw) <- locus_levels_master

message("\nDecay rates:")
print(round(as.data.frame(lambda_summ[, c("mean", "sd", "2.5%", "97.5%",
                                          "n_eff", "Rhat")]), 4))

message("\nObservation model parameters:")
.obs <- summary(fit_decay_all, pars = c("sigma_obs", "alpha", "beta", "logC50"),
                probs = c(0.025, 0.5, 0.975))$summary
print(round(.obs[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 4))
write.csv(data.frame(parameter = rownames(.obs), .obs, row.names = NULL),
          out_path("decay_observation_parameters.csv"), row.names = FALSE)

pc_all <- compute_p_and_C0_L(dat_all %>% filter(Hours_base == 0),
                             locus_levels_master)
p_raw        <- pc_all$p
names(p_raw) <- locus_levels_master

message("\nt = 0 summary per locus (log copies/L):")
print(as.data.frame(pc_all$table))
message("C0_mean (reference locus, log copies/L): ", round(pc_all$C0_mean, 4))

r_fit        <- r_raw[names(locus_to_marker)]
names(r_fit) <- unname(locus_to_marker)
p_fit        <- p_raw[names(locus_to_marker)]
names(p_fit) <- unname(locus_to_marker)

carboy_rates <- data.frame(
  marker  = names(r_fit),
  r       = as.numeric(r_fit),
  p       = as.numeric(p_fit),
  p_ratio = exp(as.numeric(p_fit)),
  r_sd    = as.numeric(lambda_summ[, "sd"])
)

write.csv(carboy_rates, out_path("carboy_rates.csv"), row.names = FALSE)
message("\nSaved: outputs/simplified/carboy_rates.csv")
print(carboy_rates %>% mutate(across(where(is.numeric), ~round(.x, 4))))

# =============================================================================
# Part B. Leave-one-carboy-out validation
# =============================================================================

message("\n=== Part B: leave-one-out carboy validation ===")
message("sigma prior: Half-Normal(0, ", SIGMA_SD, ")")

folds <- list(
  list(train = c(1, 2), test = 3, label = "1+2"),
  list(train = c(1, 3), test = 2, label = "1+3"),
  list(train = c(2, 3), test = 1, label = "2+3")
)

pred_list  <- list()
sigma_list <- list()

for (f in folds) {

  message("\n--- Fold: train = ", paste(f$train, collapse = "+"),
          " | test = ", f$test, " ---")

  dat_train <- decay_raw %>%
    filter(!Control, Component %in% components_use,
           Hours_base <= hours_max_decay, Carboy %in% f$train)

  decay_sd <- build_decay_data(dat_train, locus_levels_master)

  fit_decay <- stan(
    file   = stan_decay,
    data   = decay_sd,
    init   = function() list(C0     = rep(9,    decay_sd$N_ik),
                             lambda = rep(0.15, decay_sd$N_k)),
    chains = 2, iter = 2500, warmup = 500, seed = 42 + f$test, refresh = 0
  )

  r_vec <- -as.numeric(summary(fit_decay, pars = "lambda")$summary[, "mean"])
  pc    <- compute_p_and_C0_L(dat_train %>% filter(Hours_base == 0),
                              locus_levels_master)

  dat_test <- decay_raw %>%
    filter(!Control, Component %in% components_use,
           Carboy == f$test, Hours_base < hours_max_time)

  built <- build_age_conc_data(dat_test, locus_levels_master)

  message("  test replicates: ", nrow(built$obs),
          " | non-detections: ", sum(built$obs$detected == 0),
          " | timepoints: ", nrow(built$tp_map))

  stan_data <- build_conc_stan_data(
    obs     = built$obs,
    Nt      = nrow(built$tp_map),
    Nloci   = length(locus_levels_master),
    r_vec   = r_vec,
    p_vec   = pc$p,
    C0_mean = pc$C0_mean,
    C0_sd   = pc$C0_sd,
    t_mean  = t_mean_prior,
    t_sd    = t_sd_prior
  )

  fit_time <- stan(
    file    = stan_conc_age,
    data    = stan_data,
    chains  = 4, iter = 4000, warmup = 2000, seed = 100 + f$test,
    control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 200
  )

  t_sum <- summary(fit_time, pars = "t",
                   probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  pred_list[[f$label]] <- tibble(obs_i = seq_len(nrow(t_sum))) %>%
    bind_cols(as.data.frame(t_sum)) %>%
    rename(lo95 = `2.5%`, lo50 = `25%`, med = `50%`,
           hi50 = `75%`, hi95 = `97.5%`) %>%
    left_join(built$tp_map, by = "obs_i") %>%
    mutate(Carboy = f$test, Fold = f$label)

  .so <- summary(fit_time, pars = "sigma_tech",
                 probs = c(0.025, 0.975))$summary
  sigma_list[[f$label]] <- tibble(Fold       = f$label,
                                  sigma_tech = .so[1, "mean"],
                                  lo95       = .so[1, "2.5%"],
                                  hi95       = .so[1, "97.5%"])
  message("  ", sprintf("sigma_tech = %.3f", .so[1, "mean"]))
}

pred_all <- bind_rows(pred_list) %>% mutate(Timepoint = as.numeric(Timepoint))

loo_skill <- pred_all %>%
  summarise(n        = n(),
            bias     = mean(med - True_Hours),
            mae      = mean(abs(med - True_Hours)),
            rmse     = sqrt(mean((med - True_Hours)^2)),
            cover95  = mean(True_Hours >= lo95 & True_Hours <= hi95),
            cover50  = mean(True_Hours >= lo50 & True_Hours <= hi50),
            ci95_med = median(hi95 - lo95))

message("\nLeave-one-out skill (simplified model):")
print(as.data.frame(round(loo_skill, 3)))

write.csv(pred_all,  out_path("carboy_loo_predictions.csv"), row.names = FALSE)
write.csv(loo_skill, out_path("carboy_loo_skill.csv"),       row.names = FALSE)
write.csv(bind_rows(sigma_list), out_path("carboy_loo_sigma.csv"),
          row.names = FALSE)

# --- Panel (a): observed decay ----------------------------------------------

decay_plotdat <- decay_raw %>%
  filter(!Control, Component %in% components_use, Hours_base <= 30,
         Marker %in% c("Cytb", "16S", "Dloop", "Bridge")) %>%
  mutate(copies_mL = as.numeric(copies_mL),
         locus     = paste(Marker, Component, sep = "_"),
         Carboy    = factor(Carboy))

set2_cols        <- brewer.pal(4, "Set2")[c(1, 2, 4, 3)]
names(set2_cols) <- c("16S_DNA", "Cytb_DNA", "Bridge_DNA", "Dloop_DNA")
set2_cols_use    <- set2_cols[names(set2_cols) %in% unique(decay_plotdat$locus)]

plot_locus_order <- c("Cytb_DNA", "16S_DNA", "Dloop_DNA", "Bridge_DNA")
plot_locus_names <- c("Cytb", "16S", "D-loop", "Bridge")
lambda_labels    <- setNames(
  sprintf("%s (λ = %.3f)", plot_locus_names, as.numeric(r_raw)),
  plot_locus_order)

p_carboy_decay <- ggplot(decay_plotdat,
                         aes(x = Hours_base, y = copies_mL * 1000,
                             color = locus, shape = Carboy)) +
  geom_point(alpha = 0.6) +
  geom_smooth(aes(group = locus), method = "lm", se = FALSE, na.rm = TRUE) +
  scale_color_manual(values = set2_cols_use, breaks = plot_locus_order,
                     labels = lambda_labels, name = "Marker") +
  scale_shape_manual(values = c(16, 17, 15), name = "Carboy") +
  labs(title = "Observed eDNA decay across carboys",
       x = "Time elapsed (hours)", y = "eDNA concentration (copies/L)") +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 30)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# --- Panel (b): LOO ages ----------------------------------------------------

p_loo_time <- ggplot(pred_all, aes(x = True_Hours, color = factor(Carboy))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_linerange(aes(y = med, ymin = lo95, ymax = hi95),
                 linewidth = 0.6, alpha = 0.5) +
  geom_linerange(aes(y = med, ymin = lo50, ymax = hi50), linewidth = 1.4) +
  geom_point(aes(y = med), size = 2.5) +
  scale_color_manual(values = c("1" = "#C6DBEF", "2" = "#6BAED6",
                                "3" = "#08519C"), name = "Carboy") +
  ylim(0, max(pred_all$hi95, na.rm = TRUE)) +
  xlim(0, 30) +
  labs(title = "Leave-one-out age estimations",
       x = "Actual elapsed time (hours)",
       y = "Estimated elapsed time (hours)", color = "Carboy") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

carboy_loo_grid <- plot_grid(p_carboy_decay, p_loo_time,
                             labels = c("(a)", "(b)"), label_size = 16,
                             ncol = 2, rel_widths = c(0.55, 0.45))

ggsave(plot_path("carboy_loo_grid.png"), plot = carboy_loo_grid,
       width = 12, height = 5, dpi = 300)

message("Saved: plots/simplified/carboy_loo_grid.png")
