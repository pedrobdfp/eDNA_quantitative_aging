# =============================================================================
# 06_field_removal_ddPCR.R
# -----------------------------------------------------------------------------
# Field age estimation with a physical-removal term added to each marker's
# effective decay rate:
#
#     r_effective[j] = r_biological[j] + rho,     rho = -0.45 /hr
#
# rho is the physical transport component reported by Brasseale et al. (2025,
# JGR Oceans) for a tracer-dispersion model fitted at the same Hood Canal site.
# Ages are plotted beside the biological-decay-only estimates to show how
# accounting for dilution compresses inferred ages.
#
# Requires 04_field_main_ddPCR.R to have been run
# (outputs/ddPCR/age_estimates_grabs.csv must exist).
#
# Outputs:
#   plots/ddPCR/field_removal.png
#   outputs/ddPCR/age_estimates_grabs_removal.csv
# =============================================================================

source(here::here("code", "ddPCR", "00_setup_ddPCR.R"))

set.seed(44)

# --- Load baseline age estimates (from 04_field_main_ddPCR.R) ------------------

baseline_file <- out_path("age_estimates_grabs.csv")
if (!file.exists(baseline_file)) {
  stop("Run 04_field_main_ddPCR.R first to produce ", baseline_file)
}
baseline_ages <- read.csv(baseline_file, stringsAsFactors = FALSE)

# =============================================================================
# Rebuild the grab-level droplet data (same pipeline as 04_field_main_ddPCR.R)
# =============================================================================

d <- load_field_droplets() %>%
  mutate(marker = as.character(marker)) %>%
  filter(marker %in% locus_levels)

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

obs <- d %>%
  mutate(obs_i = match(grab_id, unique(grab_id)),
         obs_j = match(marker,  locus_levels)) %>%
  arrange(obs_i, obs_j)

ESP_data  <- read.csv(data_path("ESP_timestamps_mLseawater.csv"),
                      stringsAsFactors = FALSE)
C0_ref_m  <- mean(log(ESP_data$copies_per_mLseawater * 1000))
C0_ref_sd <- sd(  log(ESP_data$copies_per_mLseawater * 1000))

# --- Apply physical removal to the decay rates ------------------------------

rho_removal <- -0.45
r_removal   <- r_by_marker + rho_removal

message("Adjusted decay rates (biological + physical removal):")
print(round(r_removal, 4))

# --- Assemble stan_data ------------------------------------------------------

Nt    <- n_distinct(obs$obs_i)
Nloci <- length(locus_levels)

field_stan_data <- build_droplet_stan_data(
  obs     = obs,
  Nt      = Nt,
  Nloci   = Nloci,
  r_vec   = unname(r_removal[locus_levels]),
  p_vec   = unname(p_by_marker[locus_levels]),
  C0_mean = C0_ref_m,
  C0_sd   = C0_ref_sd,
  t_mean  = 12,
  t_sd    = 24,
  bio     = obs$sample
)

cat("Grabs (Nt):", Nt, " | wells (N):", field_stan_data$N, "\n")
cat("Zero wells:", sum(obs$positives == 0), "\n")
cat("sigma prior: Half-Normal(0,", SIGMA_SD, ")\n")

# --- Fit ---------------------------------------------------------------------

fit_removal <- stan(
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
.so <- summary(fit_removal, pars = c("sigma_bio", "sigma_tech", "sigma_total"),
                probs = c(0.025, 0.5, 0.975))$summary

print(round(.so[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 4))

# --- Extract posterior-median age per grab ----------------------------------

t_mat_removal <- rstan::extract(fit_removal, pars = "t")$t

grab_station <- obs %>% distinct(obs_i, station) %>% arrange(obs_i)

removal_ages <- grab_station %>%
  mutate(med_t = apply(t_mat_removal, 2, median))

# =============================================================================
# Comparison plot: baseline vs with-removal, by station
# =============================================================================

boxplot_data <- bind_rows(
  baseline_ages %>%
    transmute(station, med_t = t_med,
              scenario = "No removal\n(biological decay only)"),
  removal_ages %>%
    transmute(station, med_t,
              scenario = "With physical removal\n(ρ = -0.45 /hr)")
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("No removal\n(biological decay only)",
                                      "With physical removal\n(ρ = -0.45 /hr)")))

y_max <- max(boxplot_data$med_t, na.rm = TRUE) * 1.08

make_box_panel <- function(df_s, show_y = TRUE) {
  ggplot(df_s, aes(x = station, y = med_t, fill = station)) +
    geom_boxplot(width = 0.55, linewidth = 0.5,
                 outlier.size = 1.2, outlier.alpha = 0.6) +
    scale_fill_manual(values = station_colors, labels = station_labels) +
    scale_x_discrete(labels = station_labels) +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = seq(0, ceiling(y_max / 10) * 10, by = 5)) +
    labs(title = unique(as.character(df_s$scenario)), x = NULL,
         y = if (show_y) "Estimated age t (hours)" else NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10),
          axis.title.y = element_text(size = 11),
          plot.margin  = margin(5, 10, 5, 10))
}

scenario_levels <- levels(boxplot_data$scenario)
panels <- lapply(seq_along(scenario_levels), function(i) {
  make_box_panel(boxplot_data %>% filter(scenario == scenario_levels[i]),
                 show_y = (i == 1))
})

removal_figure <- plot_grid(panels[[1]], panels[[2]],
                            labels = c("(a)", "(b)"), label_size = 12,
                            ncol = 2, align = "h", axis = "tb")

ggsave(plot_path("field_removal.png"), plot = removal_figure,
       width = 8, height = 4, dpi = 300)

message("Saved: plots/ddPCR/field_removal.png")

# --- Summaries ---------------------------------------------------------------

message("\nMedian age per station x scenario:")
boxplot_data %>%
  group_by(scenario, station) %>%
  summarise(median_med_t = median(med_t, na.rm = TRUE),
            min_med_t    = min(med_t,    na.rm = TRUE),
            max_med_t    = max(med_t,    na.rm = TRUE), .groups = "drop") %>%
  print()

# --- Export ------------------------------------------------------------------

removal_export <- removal_ages %>%
  mutate(t_med  = med_t,
         t_lo95 = apply(t_mat_removal, 2, quantile, probs = 0.025),
         t_hi95 = apply(t_mat_removal, 2, quantile, probs = 0.975),
         rho    = rho_removal) %>%
  dplyr::select(obs_i, station, rho, t_med, t_lo95, t_hi95)

write.csv(removal_export, out_path("age_estimates_grabs_removal.csv"),
          row.names = FALSE)

message("Saved: outputs/ddPCR/age_estimates_grabs_removal.csv")
