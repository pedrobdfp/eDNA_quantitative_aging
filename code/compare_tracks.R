# =============================================================================
# compare_tracks.R
# -----------------------------------------------------------------------------
# Side-by-side summary of the two observation models, for the supplement.
#
# Both tracks are fitted to the same wells with the same process model, the
# same fixed decay rates, the same priors on age and baseline concentration,
# and the same nested error structure. They differ only in what they read from
# the instrument, so every difference below is attributable to that.
#
# Run code/ddPCR/ and code/simplified/ first. Writes
# outputs/model_comparison.csv and prints the same table.
# =============================================================================

library(here)
suppressPackageStartupMessages(library(dplyr))

rd <- function(track, f) {
  p <- here("outputs", track, f)
  if (!file.exists(p)) {
    stop(f, " not found for the ", track, " track. Run that track first.")
  }
  read.csv(p, stringsAsFactors = FALSE)
}

tracks <- c(ddPCR = "ddPCR", simplified = "simplified")

# --- Decay rates -------------------------------------------------------------

rates <- bind_rows(lapply(names(tracks), function(tk) {
  rd(tracks[tk], "carboy_rates.csv") %>% mutate(track = tk)
}))

rates_wide <- rates %>%
  dplyr::select(track, marker, r) %>%
  tidyr::pivot_wider(names_from = track, values_from = r) %>%
  mutate(difference = simplified - ddPCR)

# --- Leave-one-carboy-out validation ----------------------------------------

loo <- bind_rows(lapply(names(tracks), function(tk) {
  rd(tracks[tk], "carboy_loo_skill.csv") %>% mutate(track = tk)
})) %>%
  dplyr::select(track, n, bias, mae, rmse, cover50, cover95, ci95_med)

# --- Field ages, grab level --------------------------------------------------

ages <- bind_rows(lapply(names(tracks), function(tk) {
  rd(tracks[tk], "age_estimates_grabs.csv") %>% mutate(track = tk)
}))

ages_by_station <- ages %>%
  group_by(track, station) %>%
  summarise(n_grabs      = n(),
            median_age   = median(t_med),
            min_age      = min(t_med),
            max_age      = max(t_med),
            median_ci95  = median(t_hi95 - t_lo95),
            .groups = "drop")

ages_paired <- ages %>%
  dplyr::select(track, grab_id, t_med, t_lo95, t_hi95) %>%
  tidyr::pivot_wider(names_from = track,
                     values_from = c(t_med, t_lo95, t_hi95))

age_agreement <- ages_paired %>%
  summarise(n                 = n(),
            mean_difference   = mean(t_med_simplified - t_med_ddPCR),
            max_abs_difference = max(abs(t_med_simplified - t_med_ddPCR)),
            correlation       = cor(t_med_ddPCR, t_med_simplified),
            median_ci95_ddPCR      = median(t_hi95_ddPCR - t_lo95_ddPCR),
            median_ci95_simplified = median(t_hi95_simplified - t_lo95_simplified))

# --- Variance components -----------------------------------------------------

sigmas <- bind_rows(lapply(names(tracks), function(tk) {
  rd(tracks[tk], "all_parameters_grabs.csv") %>%
    filter(grepl("^sigma", parameter)) %>%
    transmute(track = tk, parameter, mean,
              lo95 = X2.5., hi95 = X97.5.)
}))

# --- Report ------------------------------------------------------------------

r3 <- function(d) d %>% mutate(across(where(is.numeric), ~round(.x, 3)))

cat("\n=== Decay rates (1/hour) ===\n");        print(as.data.frame(r3(rates_wide)))
cat("\n=== Leave-one-carboy-out validation ===\n"); print(as.data.frame(r3(loo)))
cat("\n=== Field grabs, by station ===\n");     print(as.data.frame(r3(ages_by_station)))
cat("\n=== Field grabs, paired agreement ===\n"); print(as.data.frame(r3(age_agreement)))
cat("\n=== Observation error ===\n");           print(as.data.frame(r3(sigmas)))

comparison <- bind_rows(
  rates_wide      %>% mutate(section = "decay rates")   %>% mutate(across(everything(), as.character)),
  loo             %>% mutate(section = "carboy LOO")    %>% mutate(across(everything(), as.character)),
  ages_by_station %>% mutate(section = "field grabs")   %>% mutate(across(everything(), as.character)),
  age_agreement   %>% mutate(section = "field agreement") %>% mutate(across(everything(), as.character)),
  sigmas          %>% mutate(section = "observation error") %>% mutate(across(everything(), as.character))
) %>% dplyr::select(section, everything())

write.csv(comparison, here("outputs", "model_comparison.csv"), row.names = FALSE)
message("\nSaved: outputs/model_comparison.csv")
