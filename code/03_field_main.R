# =============================================================================
# 03_field_main.R
# -----------------------------------------------------------------------------
# Field-data age estimation -- main-text figure.
#
# Biological replicates collected at the same station, depth, and instrument
# within 10 minutes of each other are pooled into a single "grab" and receive
# a single joint age estimate. This is the main-text analysis.
#
# Outputs:
#   - plots/field_grabs.png                 (main-text figure)
#   - outputs/age_estimates_grabs.csv       (per-grab age posteriors)
#   - outputs/all_parameters_grabs.csv      (all other parameters)
#
# The individual-sample version (no grab pooling) is the supplementary
# figure and is produced by 04_field_supplement.R.
# =============================================================================

source(here::here("code", "00_setup.R"))

set.seed(44)

# --- Load and filter raw data -----------------------------------------------

d <- readRDS(here("data", "ddpcr_combined_all.rds")) %>%
  filter(
    station   %in% c("NOAA Boat", "Husbandry Area"),
    component == "DNA1 (column 1)",                 # remove carryover column
    haplotype %in% c("mixed", "mixed0", "unknown"), # drop cytb haplotype duplicates
    !is.na(time_collected),
    time_collected <= 1730
  ) %>%
  group_by(sample) %>%
  filter(sum(conc_copies_L, na.rm = TRUE) > 0) %>%  # drop all-zero samples
  ungroup()

# --- Convert HHMM integer -> minutes since midnight -------------------------

d <- d %>%
  mutate(time_min = (as.integer(time_collected) %/% 100L) * 60L +
                    (as.integer(time_collected) %% 100L))

# --- Assign grab IDs ---------------------------------------------------------
# A "grab" = same station + depth + instrument, samples within 10 min of
# each other. Within each (station, instrument, depth) group, sort by time
# and start a new grab whenever the gap exceeds 10 min.

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

# --- Build Stan observation arrays (obs_i indexes grab, not sample) ---------

obs       <- d[d$marker %in% locus_levels & !is.na(d$binary), ]
obs$obs_i <- match(obs$grab_id, unique(obs$grab_id))
obs$obs_j <- match(obs$marker,  locus_levels)
obs$z     <- as.integer(obs$binary)

det       <- obs[obs$z == 1 & obs$conc_copies_L > 0, ]
det$y_obs <- log(det$conc_copies_L)

# --- C0 prior from ESP reference data (Brasseale et al.) --------------------

ESP_data  <- read.csv(here("data", "ESP_timestamps_mLseawater.csv"),
                      stringsAsFactors = FALSE)
C0_ref_m  <- mean(log(ESP_data$copies_per_mLseawater * 1000))
C0_ref_sd <- sd(  log(ESP_data$copies_per_mLseawater * 1000))

# --- Assemble stan_data ------------------------------------------------------

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
  t_mean  = 12,
  t_sd    = 24
)

# --- Sanity checks -----------------------------------------------------------

cat("Grabs (Nt):",  Nt,    "\n")
cat("Loci (Nloci):", Nloci, "\n")
cat("Locus order and r values:\n")
print(data.frame(locus = locus_levels, r = field_stan_data$r))
cat("\nTotal observations (N_obs):", field_stan_data$N_obs, "\n")
cat("Detected observations (N_y):", field_stan_data$N_y,   "\n")
cat("Overall detection rate:",
    round(field_stan_data$N_y / field_stan_data$N_obs, 3), "\n")

cat("\nDetection rate by locus:\n")
print(tapply(obs$z, obs$obs_j, mean) |> setNames(locus_levels) |> round(3))

cat("\nObservations per grab (min/median/max):\n")
print(summary(tabulate(obs$obs_i, nbins = Nt)))

# --- Fit ---------------------------------------------------------------------

fit_real <- stan(
  file    = stan_conc_field,
  data    = field_stan_data,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  seed    = 44,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 500
)

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

# Mean collection time per grab (for labels and ordering)
grab_time <- obs %>%
  group_by(grab_id, obs_i) %>%
  summarise(mean_time_min = mean(time_min, na.rm = TRUE), .groups = "drop")

# Per-grab marker proportions: mean concentration across replicates/markers,
# then renormalise within grab.
out <- obs %>%
  arrange(obs_i) %>%
  dplyr::select(grab_id, station, Instrument, depth, marker,
                conc_copies_ul, obs_i) %>%
  group_by(grab_id, station, Instrument, depth, marker, obs_i) %>%
  summarise(conc = mean(conc_copies_ul, na.rm = TRUE), .groups = "drop") %>%
  group_by(grab_id, station, Instrument, depth, obs_i) %>%
  mutate(marker_prop = conc / sum(conc, na.rm = TRUE)) %>%
  dplyr::select(-conc) %>%
  ungroup() %>%
  pivot_wider(names_from  = marker,
              values_from = marker_prop,
              values_fill = 0) %>%
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

# ---- Panel (a): marker composition per grab --------------------------------

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
  labs(x     = "Station / depth / collection time",
       y     = "Observed proportion",
       title = "Observed marker composition") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 8), legend.position = "top")
gg_prop
# ---- Panel (a2): total eDNA concentration per grab (slim lollipop) ---------

conc_df <- obs %>%
  group_by(obs_i, sample) %>%
  summarise(sample_total = sum(conc_copies_L, na.rm = TRUE), .groups = "drop") %>%
  group_by(obs_i) %>%
  summarise(total_conc = mean(sample_total), .groups = "drop") %>%
  left_join(out %>% dplyr::select(obs_i, sample_label), by = "obs_i")

gg_conc <- ggplot(conc_df, aes(x = sample_label, y = total_conc)) +
  geom_segment(aes(xend = sample_label, y = 0, yend = total_conc),
               color = "black", linewidth = 0.6) +
  geom_point(size = 1.2, color = "black") +
  coord_flip() +
  labs(x = NULL, y = "Total eDNA\n(copies/L)", title = "") +
  # Invisible legend point to keep vertical alignment with the other panels
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

# ---- Panel (b): estimated age per grab --------------------------------------

gg_t <- ggplot(out, aes(x = sample_label, y = med, color = station)) +
  geom_linerange(aes(ymin = lo95, ymax = hi95),
                 linewidth = 0.6, alpha = 0.5) +
  geom_linerange(aes(ymin = lo50, ymax = hi50),
                 linewidth = 1.4) +
  geom_point(size = 2.5) +
  coord_flip() +
  scale_color_manual(values = station_colors, labels = station_labels) +
  labs(x     = "Station / depth / collection time",
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
field_grid
ggsave(here("plots", "field_grabs.png"),
       plot = field_grid, width = 14, height = 8, dpi = 300)

message("Saved: plots/field_grabs.png")

# --- Summary by station (for the text) --------------------------------------

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

# ---- Per-grab age posteriors -----------------------------------------------

grab_sample_lookup <- obs %>%
  distinct(obs_i, grab_id, sample) %>%
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
                n_samples       = n_distinct(sample),
                .groups         = "drop"),
    by = "obs_i"
  ) %>%
  left_join(out %>% dplyr::select(obs_i, station, depth, mean_time_min),
            by = "obs_i") %>%
  dplyr::select(grab_id, obs_i, station, depth, mean_time_min,
                n_samples, samples_in_grab,
                t_med, t_lo95, t_hi95)

write.csv(age_table_grabs,
          here("outputs", "age_estimates_grabs.csv"),
          row.names = FALSE)

# ---- All non-t parameters (supplementary table) ----------------------------

all_params_grabs <- rstan::summary(fit_real)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column("parameter") %>%
  filter(!grepl("^t\\[", parameter))

write.csv(all_params_grabs,
          here("outputs", "all_parameters_grabs.csv"),
          row.names = FALSE)

message("Saved: outputs/age_estimates_grabs.csv + all_parameters_grabs.csv")
