# =============================================================================
# 05_field_removal.R
# -----------------------------------------------------------------------------
# Field-data age estimation with physical-removal sensitivity.
#
# We re-fit the grab-level field model after adding a physical-removal term
# to each marker's effective decay rate:
#
#     r_effective[j] = r_biological[j] + rho,    rho = -0.45 /hr
#
# rho = -0.45 /hr is the physical transport component reported by
# Brasseale et al. (2025, JGR Oceans) for a tracer-dispersion model
# fit to the same Hood Canal site where the field samples were collected.
#
# The resulting age estimates are plotted side-by-side with the baseline
# (biological decay only) estimates from 03_field_main.R to show how
# accounting for physical dilution compresses inferred ages.
#
# Requires: 03_field_main.R has already been run
#           (outputs/age_estimates_grabs.csv must exist).
# =============================================================================

source(here::here("code", "00_setup.R"))

set.seed(44)

# --- Load baseline age estimates (from 03_field_main.R) ---------------------

baseline_file <- here("outputs", "age_estimates_grabs.csv")
if (!file.exists(baseline_file)) {
  stop("Run 03_field_main.R first to produce ", baseline_file)
}

baseline_ages <- read.csv(baseline_file, stringsAsFactors = FALSE)

# =============================================================================
# Rebuild field stan_data and re-fit with r adjusted for physical removal
# (Same data pipeline as 03_field_main.R.)
# =============================================================================

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
  ungroup() %>%
  mutate(time_min = (as.integer(time_collected) %/% 100L) * 60L +
                    (as.integer(time_collected) %% 100L))

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

obs       <- d[d$marker %in% locus_levels & !is.na(d$binary), ]
obs$obs_i <- match(obs$grab_id, unique(obs$grab_id))
obs$obs_j <- match(obs$marker,  locus_levels)
obs$z     <- as.integer(obs$binary)

det       <- obs[obs$z == 1 & obs$conc_copies_L > 0, ]
det$y_obs <- log(det$conc_copies_L)

ESP_data  <- read.csv(here("data", "ESP_timestamps_mLseawater.csv"),
                      stringsAsFactors = FALSE)
C0_ref_m  <- mean(log(ESP_data$copies_per_mLseawater * 1000))
C0_ref_sd <- sd(  log(ESP_data$copies_per_mLseawater * 1000))

# --- Apply physical removal to the decay rates ------------------------------

rho_removal <- -0.45
r_removal   <- r_by_marker + rho_removal

message("Adjusted decay rates (biological + physical removal):")
print(round(r_removal, 4))

# --- Assemble stan_data (baseline shape from 03, but with r_removal) --------

field_stan_data <- list(
  Nt    = length(unique(obs$obs_i)),
  Nloci = length(locus_levels),
  r     = unname(r_removal[locus_levels]),

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

# --- Fit ---------------------------------------------------------------------

fit_removal <- stan(
  file    = stan_conc_field,
  data    = field_stan_data,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 44,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 500
)

# --- Extract posterior-median age per grab + attach station -----------------

t_mat_removal <- rstan::extract(fit_removal, pars = "t")$t

# Grab-to-station lookup (from the same obs we just built)
grab_station <- obs %>%
  distinct(obs_i, station) %>%
  arrange(obs_i)

removal_ages <- grab_station %>%
  mutate(med_t = apply(t_mat_removal, 2, median))

# =============================================================================
# Compose comparison plot: baseline vs with-removal, by station
# =============================================================================

boxplot_data <- bind_rows(
  baseline_ages %>%
    transmute(station, med_t = t_med,
              scenario = "No removal\n(biological decay only)"),
  removal_ages %>%
    transmute(station,
              med_t,
              scenario = "With physical removal\n(\u03c1 = -0.45 /hr)")
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("No removal\n(biological decay only)",
                                      "With physical removal\n(\u03c1 = -0.45 /hr)")))

y_max <- max(boxplot_data$med_t, na.rm = TRUE) * 1.08

make_box_panel <- function(df_s, show_y = TRUE) {
  ggplot(df_s, aes(x = station, y = med_t, fill = station)) +
    geom_boxplot(width = 0.55, linewidth = 0.5,
                 outlier.size = 1.2, outlier.alpha = 0.6) +
    scale_fill_manual(values = station_colors, labels = station_labels) +
    scale_x_discrete(labels = station_labels) +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = seq(0, ceiling(y_max / 10) * 10, by = 5)) +
    labs(title = unique(as.character(df_s$scenario)),
         x     = NULL,
         y     = if (show_y) "Estimated age t (hours)" else NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text.x     = element_text(size = 10),
      axis.text.y     = element_text(size = 10),
      axis.title.y    = element_text(size = 11),
      plot.margin     = margin(5, 10, 5, 10)
    )
}

scenario_levels <- levels(boxplot_data$scenario)
panels <- lapply(seq_along(scenario_levels), function(i) {
  make_box_panel(
    df_s   = boxplot_data %>% filter(scenario == scenario_levels[i]),
    show_y = (i == 1)
  )
})

removal_figure <- plot_grid(
  panels[[1]], panels[[2]],
  labels     = c("(a)", "(b)"),
  label_size = 12,
  ncol       = 2,
  align      = "h",
  axis       = "tb"
)
removal_figure
ggsave(here("plots", "field_removal.png"),
       plot = removal_figure, width = 8, height = 4, dpi = 300)

message("Saved: plots/field_removal.png")

# --- Summaries for the text --------------------------------------------------

message("\nMedian age per station x scenario:")
boxplot_data %>%
  group_by(scenario, station) %>%
  summarise(median_med_t = median(med_t, na.rm = TRUE),
            min_med_t    = min(med_t,    na.rm = TRUE),
            max_med_t    = max(med_t,    na.rm = TRUE),
            .groups      = "drop") %>%
  print()

# --- Export removal-scenario age table --------------------------------------

removal_export <- removal_ages %>%
  mutate(
    t_med  = med_t,
    t_lo95 = apply(t_mat_removal, 2, quantile, probs = 0.025),
    t_hi95 = apply(t_mat_removal, 2, quantile, probs = 0.975),
    rho    = rho_removal
  ) %>%
  dplyr::select(obs_i, station, rho, t_med, t_lo95, t_hi95)

write.csv(removal_export,
          here("outputs", "age_estimates_grabs_removal.csv"),
          row.names = FALSE)

message("Saved: outputs/age_estimates_grabs_removal.csv")
