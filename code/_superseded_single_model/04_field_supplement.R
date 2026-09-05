# =============================================================================
# 04_field_supplement.R
# -----------------------------------------------------------------------------
# Field-data age estimation -- supplementary figure.
#
# Every individual sample receives its own age estimate; biological replicates
# are NOT pooled. This version appears in the supplementary material; the
# main-text analysis (grab-averaged) is in 03_field_main.R.
#
# Outputs:
#   - plots/field_individual_samples.png
#   - outputs/age_estimates_individual_samples.csv
#   - outputs/all_parameters_individual_samples.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

set.seed(44)

# --- Load and filter raw data -----------------------------------------------

d <- readRDS(here("data", "ddpcr_combined_all.rds")) %>%
  filter(
    station   %in% c("NOAA Boat", "Husbandry Area"),
    component == "DNA1 (column 1)",
    haplotype %in% c("mixed", "mixed0", "unknown"),
    !is.na(time_collected),
    time_collected <= 1730
  ) %>%
  group_by(sample) %>%
  filter(sum(conc_copies_L, na.rm = TRUE) > 0) %>%
  ungroup()

sample_levels <- sort(unique(d$sample))

# --- Build Stan observation arrays (obs_i indexes individual sample) --------

obs       <- d[d$marker %in% locus_levels & !is.na(d$binary), ]
obs$obs_i <- match(obs$sample, unique(obs$sample))
obs$obs_j <- match(obs$marker, locus_levels)
obs$z     <- as.integer(obs$binary)

det       <- obs[obs$z == 1 & obs$conc_copies_L > 0, ]
det$y_obs <- log(det$conc_copies_L)

# --- C0 prior from ESP reference data (Brasseale et al.) --------------------

ESP_data  <- read.csv(here("data", "ESP_timestamps_mLseawater.csv"),
                      stringsAsFactors = FALSE)
C0_ref_m  <- mean(log(ESP_data$copies_per_mLseawater * 1000))
C0_ref_sd <- sd(  log(ESP_data$copies_per_mLseawater * 1000))

# --- Assemble stan_data ------------------------------------------------------
# Note: here we use a broader prior on t (t_mean = 0, t_sd = 40) to avoid a
# t-C0 degeneracy that arose with a tighter prior in this individual-sample
# regime.

Nt    <- length(unique(obs$obs_i))
Nloci <- length(locus_levels)

field_stan_data <- list(
  Nt    = Nt,
  Nloci = Nloci,
  r     = unname(r_by_marker[locus_levels]),

  N_obs = nrow(obs),
  obs_i = obs$obs_i,
  obs_j = obs$obs_j,
  z     = obs$z,

  N_y   = nrow(det),
  y_i   = det$obs_i,
  y_j   = det$obs_j,
  y_obs = det$y_obs,

  p       = unname(p_by_marker[locus_levels]),
  C0_mean = C0_ref_m,
  C0_sd   = C0_ref_sd,
  t_mean  = 0,
  t_sd    = 40
)

# --- Sanity checks -----------------------------------------------------------

cat("Samples (Nt):", Nt,    "\n")
cat("Loci (Nloci):", Nloci, "\n")
cat("Locus order and r values:\n")
print(data.frame(locus = locus_levels, r = field_stan_data$r))
cat("\nTotal observations (N_obs):", field_stan_data$N_obs, "\n")
cat("Detected observations (N_y):", field_stan_data$N_y,   "\n")
cat("Overall detection rate:",
    round(field_stan_data$N_y / field_stan_data$N_obs, 3), "\n")

cat("\nDetection rate by locus:\n")
print(tapply(obs$z, obs$obs_j, mean) |> setNames(locus_levels) |> round(3))

cat("\nlog(conc_copies_L) summary for detected obs:\n")
print(summary(field_stan_data$y_obs))

cat("\nObservations per sample (min/median/max):\n")
print(summary(tabulate(obs$obs_i, nbins = Nt)))

# --- Fit ---------------------------------------------------------------------

fit_real <- stan(
  file    = stan_conc_field,
  data    = field_stan_data,
  chains  = 3,
  iter    = 4000,
  warmup  = 2000,
  seed    = 44,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  refresh = 200
)

# =============================================================================
# Build `out` -- per-sample posterior summaries and marker composition
# =============================================================================

t_mat <- rstan::extract(fit_real, pars = "t")$t

out_t <- data.frame(idx = obs$obs_i) %>%
  distinct() %>%
  arrange(idx) %>%
  mutate(
    med  = apply(t_mat, 2, quantile, probs = 0.50),
    lo95 = apply(t_mat, 2, quantile, probs = 0.025),
    lo50 = apply(t_mat, 2, quantile, probs = 0.25),
    hi50 = apply(t_mat, 2, quantile, probs = 0.75),
    hi95 = apply(t_mat, 2, quantile, probs = 0.975)
  )

out <- obs %>%
  arrange(obs_i) %>%
  dplyr::select(sample, station, Instrument, depth, marker,
                conc_copies_ul, obs_i, id, time_collected) %>%
  group_by(sample, station, Instrument, depth, marker, obs_i, id, time_collected) %>%
  summarise(conc = mean(conc_copies_ul, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample, station, Instrument, depth, obs_i, id, time_collected) %>%
  mutate(marker_prop = conc / sum(conc, na.rm = TRUE)) %>%
  dplyr::select(-conc) %>%
  ungroup() %>%
  pivot_wider(names_from  = marker,
              values_from = marker_prop,
              values_fill = 0) %>%
  arrange(obs_i) %>%
  left_join(out_t, by = join_by(obs_i == idx)) %>%
  arrange(station, time_collected) %>%
  mutate(
    time_str = sprintf("%02d:%02d",
                       as.integer(time_collected) %/% 100L,
                       as.integer(time_collected) %% 100L),
    sample_label = factor(
      paste0(id, " - ", time_str),
      levels = unique(paste0(id, " - ", time_str))
    ),
    station = factor(station, levels = c("Husbandry Area", "NOAA Boat"))
  )

# =============================================================================
# Plots
# =============================================================================

# ---- Panel (a): marker composition per sample ------------------------------

prop_long <- out %>%
  dplyr::select(sample_label, all_of(locus_levels)) %>%
  pivot_longer(cols      = all_of(locus_levels),
               names_to  = "marker",
               values_to = "prop") %>%
  mutate(marker = factor(marker, levels = locus_levels))

gg_prop <- ggplot(prop_long,
                  aes(x = sample_label, y = prop, fill = marker)) +
  geom_col(width = 0.85) +
  coord_flip() +
  scale_fill_manual(values = marker_colors, name = "Marker",
                    labels = marker_display) +
  labs(x     = "Sample ID - Collection time",
       y     = "Observed proportion",
       title = "Observed marker composition") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 8), legend.position = "top")
gg_prop
# ---- Panel (a2): total eDNA concentration per sample (slim lollipop) -------

conc_df <- obs %>%
  group_by(obs_i) %>%
  summarise(total_conc = sum(conc_copies_L, na.rm = TRUE), .groups = "drop") %>%
  left_join(out %>% dplyr::select(obs_i, sample_label), by = "obs_i")

gg_conc <- ggplot(conc_df, aes(x = sample_label, y = total_conc)) +
  geom_segment(aes(xend = sample_label, y = 0, yend = total_conc),
               color = "black", linewidth = 0.6) +
  geom_point(size = 1.2, color = "black") +
  coord_flip() +
  labs(x = NULL, y = "Total eDNA\n(copies/L)", title = "") +
  geom_point(aes(color = " "), alpha = 0) +
  scale_color_manual(values = c(" " = "white")) +
  guides(color = guide_legend(title = " ", override.aes = list(alpha = 0))) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(size = 7),
    legend.position = "top",
    legend.text     = element_text(color = "white"),
    legend.title    = element_text(color = "white"),
    legend.key      = element_rect(fill = "white", color = "white"),
    plot.margin     = margin(5, 2, 5, 2)
  )

# ---- Panel (b): estimated age per sample -----------------------------------

gg_t <- ggplot(out, aes(x = sample_label, y = med, color = station)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95),
                 linewidth = 0.6, alpha = 0.5) +
  geom_linerange(aes(ymin = lo50, ymax = hi50),
                 linewidth = 1.4) +
  geom_point(size = 2.5) +
  coord_flip() +
  scale_color_manual(values = station_colors, labels = station_labels) +
  labs(x     = "Sample ID - Collection time",
       y     = "Estimated age t (hours)",
       color = "Station",
       title = "Model estimated eDNA age") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 8), legend.position = "top")
gg_t 
# ---- Compose: (a) | (a2 slim) | (b) ----------------------------------------

field_grid <- plot_grid(
  gg_prop, gg_conc, gg_t,
  labels     = c("(a)", "", "(b)"),
  label_size = 16,
  ncol       = 3,
  rel_widths = c(2, 0.5, 2),
  align      = "h",
  axis       = "tb"
)

ggsave(here("plots", "field_individual_samples.png"),
       plot = field_grid, width = 14, height = 8, dpi = 300)

message("Saved: plots/field_individual_samples.png")

# --- Summary by station (for the text) --------------------------------------

out %>%
  group_by(station) %>%
  summarise(mean_age = mean(med, na.rm = TRUE)) %>%
  print()

# =============================================================================
# Post-fit exports
# =============================================================================

age_table_samples <- tibble(
  obs_i  = seq_len(ncol(t_mat)),
  t_med  = apply(t_mat, 2, quantile, probs = 0.50),
  t_lo95 = apply(t_mat, 2, quantile, probs = 0.025),
  t_hi95 = apply(t_mat, 2, quantile, probs = 0.975)
) %>%
  left_join(obs %>% distinct(obs_i, sample) %>% arrange(obs_i), by = "obs_i") %>%
  dplyr::select(sample, obs_i, t_med, t_lo95, t_hi95)
age_table_samples
write.csv(age_table_samples,
          here("outputs", "age_estimates_individual_samples.csv"),
          row.names = FALSE)

all_params <- rstan::summary(fit_real)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column("parameter") %>%
  filter(!grepl("^t\\[", parameter))

write.csv(all_params,
          here("outputs", "all_parameters_individual_samples.csv"),
          row.names = FALSE)

message("Saved: outputs/age_estimates_individual_samples.csv",
        " + all_parameters_individual_samples.csv")

