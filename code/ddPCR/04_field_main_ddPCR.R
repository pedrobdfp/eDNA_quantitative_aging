# =============================================================================
# 04_field_main_ddPCR.R
# -----------------------------------------------------------------------------
# Age estimation for the field samples: the main-text analysis.
#
# Biological replicates taken at the same station, depth and instrument within
# 10 minutes of one another are pooled into a single grab and given one joint
# age estimate.
#
# Input is the droplet table from 01_build_field_droplets.R, one row per well
# with positive and accepted droplet counts. The prior on baseline
# concentration comes from the ESP reference data of Brasseale et al.
#
# Outputs:
#   plots/ddPCR/field_grabs.png
#   outputs/ddPCR/age_estimates_grabs.csv
#   outputs/ddPCR/all_parameters_grabs.csv
# =============================================================================

source(here::here("code", "ddPCR", "00_setup_ddPCR.R"))

set.seed(44)

# --- Load droplet-level field data ------------------------------------------

d <- load_field_droplets() %>%
  mutate(marker = as.character(marker)) %>%
  filter(marker %in% locus_levels)

message("Droplet wells loaded: ", nrow(d),
        " | samples: ", n_distinct(d$sample))

# --- Assign grab IDs ---------------------------------------------------------
# A grab is the set of water samples taken at the same station, depth and
# instrument within 10 minutes of one another. 06_field_removal_ddPCR.R applies
# the same rule.

grab_map <- d %>%
  distinct(sample, station, Instrument, depth, time_min) %>%
  arrange(station, Instrument, depth, time_min) %>%
  group_by(station, Instrument, depth) %>%
  mutate(
    new_grab = cumsum(c(TRUE, diff(time_min) > 10)),
    grab_id  = paste(station, Instrument, depth, new_grab, sep = "_")
  ) %>%
  ungroup() %>%
  dplyr::select(sample, grab_id)

d <- d %>% left_join(grab_map, by = "sample")

message("Unique grabs: ", n_distinct(d$grab_id))

# --- Build Stan observation arrays (obs_i indexes grab) ---------------------

obs <- d %>%
  mutate(obs_i = match(grab_id, unique(grab_id)),
         obs_j = match(marker,  locus_levels)) %>%
  arrange(obs_i, obs_j)

# --- C0 prior from ESP reference data (Brasseale et al.) --------------------
# log copies/L, the same scale as C in droplet_age.stan.

ESP_data  <- read.csv(data_path("ESP_timestamps_mLseawater.csv"),
                      stringsAsFactors = FALSE)
C0_ref_m  <- mean(log(ESP_data$copies_per_mLseawater * 1000))
C0_ref_sd <- sd(  log(ESP_data$copies_per_mLseawater * 1000))

# --- Assemble stan_data ------------------------------------------------------

Nt    <- n_distinct(obs$obs_i)
Nloci <- length(locus_levels)

field_stan_data <- build_droplet_stan_data(
  obs     = obs,
  Nt      = Nt,
  Nloci   = Nloci,
  r_vec   = unname(r_by_marker[locus_levels]),
  p_vec   = unname(p_by_marker[locus_levels]),
  C0_mean = C0_ref_m,
  C0_sd   = C0_ref_sd,
  t_mean  = 12,
  t_sd    = 24,
  bio     = obs$sample
)

# --- Sanity checks -----------------------------------------------------------

cat("Grabs (Nt):",   Nt,    "\n")
cat("Loci (Nloci):", Nloci, "\n")
cat("Locus order, r and p:\n")
print(data.frame(locus = locus_levels,
                 r = field_stan_data$r, p = field_stan_data$p))
cat("\nWells (N):", field_stan_data$N, "\n")
cat("Wells with at least one positive droplet:", sum(obs$positives > 0),
    "(", round(mean(obs$positives > 0), 3), ")\n")
cat("Zero wells (the informative non-detections):", sum(obs$positives == 0), "\n")

cat("\nDetection rate by locus:\n")
print(tapply(obs$positives > 0, obs$obs_j, mean) |>
        setNames(locus_levels) |> round(3))

cat("\nPositive droplets per well (min/median/max):\n")
print(summary(obs$positives))
cat("Accepted droplets per well:\n")
print(summary(obs$accepted))
cat("\nWells per grab (min/median/max):\n")
print(summary(tabulate(obs$obs_i, nbins = Nt)))
cat("\nlog_offset:\n"); print(summary(obs$log_offset))
cat("\nsigma prior: Half-Normal(0,", SIGMA_SD, ")\n")

# --- Fit ---------------------------------------------------------------------

fit_real <- stan(
  file    = stan_droplet_age,
  data    = field_stan_data,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 44,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 500
)

cat("\nsigma_bio = between water samples, sigma_tech = between wells:\n")
.so <- summary(fit_real, pars = c("sigma_bio", "sigma_tech", "sigma_total"),
                probs = c(0.025, 0.5, 0.975))$summary

print(round(.so[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 4))

# =============================================================================
# Build `out` -- per-grab posterior summaries and marker composition
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

grab_time <- obs %>%
  group_by(grab_id, obs_i) %>%
  summarise(mean_time_min = mean(time_min, na.rm = TRUE), .groups = "drop")

# Marker composition, from the machine concentrations carried in the droplet
# table (copies/uL), exactly as the concentration version plotted it.
out <- obs %>%
  arrange(obs_i) %>%
  dplyr::select(grab_id, station, Instrument, depth, marker,
                conc_machine, obs_i) %>%
  group_by(grab_id, station, Instrument, depth, marker, obs_i) %>%
  summarise(conc = mean(conc_machine, na.rm = TRUE), .groups = "drop") %>%
  group_by(grab_id, station, Instrument, depth, obs_i) %>%
  mutate(marker_prop = conc / sum(conc, na.rm = TRUE)) %>%
  dplyr::select(-conc) %>%
  ungroup() %>%
  pivot_wider(names_from = marker, values_from = marker_prop, values_fill = 0) %>%
  arrange(obs_i) %>%
  left_join(out_t,     by = join_by(obs_i == idx)) %>%
  left_join(grab_time, by = c("grab_id", "obs_i")) %>%
  arrange(station, mean_time_min) %>%
  mutate(
    time_str = sprintf("%02d:%02d",
                       as.integer(mean_time_min) %/% 60L,
                       as.integer(mean_time_min) %% 60L),
    station_short = case_when(
      station == "Husbandry Area" ~ "HUSB",
      station == "NOAA Boat"      ~ "BOAT",
      TRUE                        ~ as.character(station)
    ),
    sample_label = factor(
      paste0(station_short, " ", depth, "m - ", time_str),
      levels = unique(paste0(station_short, " ", depth, "m - ", time_str))
    ),
    station = factor(station, levels = c("Husbandry Area", "NOAA Boat"))
  )

# =============================================================================
# Plots
# =============================================================================

prop_long <- out %>%
  dplyr::select(sample_label, all_of(locus_levels)) %>%
  pivot_longer(cols = all_of(locus_levels), names_to = "marker",
               values_to = "prop") %>%
  mutate(marker = factor(marker, levels = locus_levels))

gg_prop <- ggplot(prop_long, aes(x = sample_label, y = prop, fill = marker)) +
  geom_col(width = 0.85) +
  coord_flip() +
  scale_fill_manual(values = marker_colors, name = "Marker",
                    labels = marker_display) +
  labs(x = "Station / depth / collection time", y = "Observed proportion",
       title = "Observed marker composition") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 8), legend.position = "top")

# Total eDNA per grab, back on the copies/L scale
conc_df <- obs %>%
  mutate(conc_L = conc_machine * K_CONC / vol) %>%
  group_by(obs_i, sample) %>%
  summarise(sample_total = sum(conc_L, na.rm = TRUE), .groups = "drop") %>%
  group_by(obs_i) %>%
  summarise(total_conc = mean(sample_total), .groups = "drop") %>%
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
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), axis.text.x = element_text(size = 7),
        legend.position = "top",
        legend.text  = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.key   = element_rect(fill = "white", color = "white"),
        plot.margin  = margin(5, 2, 5, 2))

gg_t <- ggplot(out, aes(x = sample_label, y = med, color = station)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6, alpha = 0.5) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.4) +
  geom_point(size = 2.5) +
  coord_flip() +
  scale_color_manual(values = station_colors, labels = station_labels) +
  labs(x = "Station / depth / collection time",
       y = "Estimated age t (hours)", color = "Station",
       title = "Model estimated eDNA age (droplet model)") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 8), legend.position = "top")

field_grid <- plot_grid(gg_prop, gg_conc, gg_t,
                        labels = c("(a)", "", "(b)"), label_size = 16,
                        ncol = 3, rel_widths = c(2, 0.5, 2),
                        align = "h", axis = "tb")

ggsave(plot_path("field_grabs.png"), plot = field_grid,
       width = 14, height = 8, dpi = 300)

message("Saved: plots/ddPCR/field_grabs.png")

# --- Summary by station ------------------------------------------------------

out %>%
  group_by(station) %>%
  summarise(mean_age  = mean(med, na.rm = TRUE),
            min_age   = min(med,  na.rm = TRUE),
            max_age   = max(med,  na.rm = TRUE),
            mean_lo95 = mean(lo95, na.rm = TRUE),
            mean_hi95 = mean(hi95, na.rm = TRUE)) %>%
  print()

# =============================================================================
# Post-fit exports
# =============================================================================

grab_sample_lookup <- obs %>% distinct(obs_i, grab_id, sample) %>%
  arrange(obs_i, sample)

age_table_grabs <- tibble(
  obs_i  = seq_len(ncol(t_mat)),
  t_med  = apply(t_mat, 2, quantile, probs = 0.50),
  t_lo95 = apply(t_mat, 2, quantile, probs = 0.025),
  t_hi95 = apply(t_mat, 2, quantile, probs = 0.975)
) %>%
  left_join(
    grab_sample_lookup %>%
      group_by(obs_i, grab_id) %>%
      summarise(samples_in_grab = paste(sort(unique(sample)), collapse = "; "),
                n_samples       = n_distinct(sample), .groups = "drop"),
    by = "obs_i") %>%
  left_join(out %>% dplyr::select(obs_i, station, depth, mean_time_min),
            by = "obs_i") %>%
  dplyr::select(grab_id, obs_i, station, depth, mean_time_min,
                n_samples, samples_in_grab, t_med, t_lo95, t_hi95)

write.csv(age_table_grabs, out_path("age_estimates_grabs.csv"),
          row.names = FALSE)

all_params_grabs <- rstan::summary(
  fit_real, pars = c("C", "sigma_bio", "sigma_tech"), probs = c(0.025, 0.5, 0.975))$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column("parameter")

write.csv(all_params_grabs, out_path("all_parameters_grabs.csv"),
          row.names = FALSE)

message("Saved: outputs/ddPCR/age_estimates_grabs.csv + all_parameters_grabs.csv")
