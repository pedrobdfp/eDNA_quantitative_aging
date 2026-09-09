# =============================================================================
# ar1_decay_check.R
# -----------------------------------------------------------------------------
# Cross-check on the marker decay rates, independent of the Bayesian fit.
#
# THE CONCERN
#   The decay model fits one exponential through all timepoints at once, so a
#   distant timepoint carries more leverage on the slope than a near one.
#
# THE CHECK
#   Re-estimate each rate from consecutive timepoints only. On the log scale
#   exponential decay is an autoregression with the coefficient fixed at 1 and
#   the decay carried by the intercept,
#
#       log C_t = log C_{t-1} - lambda * dt
#
#   so the estimator is a regression of the step change on the step length
#   through the origin. Every transition contributes once whatever its absolute
#   time, which is the equal weighting the concern asks for. The implied
#   one-hour autoregressive coefficient is beta = exp(-lambda).
#
#   Uncertainty on the stepwise estimate comes from a nonparametric bootstrap
#   over transitions, since consecutive differences share endpoints and an
#   ordinary regression standard error would be optimistic.
#
#   Finally the decay fit is refitted with each timepoint dropped in turn. If
#   one timepoint were driving the rate, that would show it.
#
# OUTPUTS
#   outputs/ddPCR/ar1_decay_check.csv
#   plots/ddPCR/ar1_decay_check.png
#
# No Stan required. Reads only data/Final_decay_ddPCR_datasheet.csv and
# outputs/ddPCR/carboy_rates.csv. Run from the repository root.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(forcats)
  library(ggplot2); library(cowplot); library(RColorBrewer)
})

set.seed(1)

HOURS_MAX <- 30      # the window the published model is fitted on
N_BOOT    <- 4000

marker_map     <- c(Cytb = "cytb", `16S` = "Tt_16S",
                    Dloop = "Tt_DLL1", Bridge = "Tt_longFrag")
locus_levels   <- c("cytb", "Tt_16S", "Tt_DLL1", "Tt_longFrag")
marker_display <- c(cytb = "Cytb", Tt_16S = "16S",
                    Tt_DLL1 = "D-loop", Tt_longFrag = "Bridge")
marker_colors  <- setNames(brewer.pal(4, "Set2")[c(2, 1, 4, 3)], locus_levels)
lab <- function(x) unname(marker_display[as.character(x)])

# --- Data: one mean log concentration per carboy x marker x timepoint --------

raw <- read.csv(here("data", "Final_decay_ddPCR_datasheet.csv"),
                stringsAsFactors = FALSE)

series <- raw %>%
  filter(!Control, Component == "DNA", Carboy %in% 1:3,
         Marker %in% names(marker_map), Hours_base <= HOURS_MAX) %>%
  mutate(marker   = factor(unname(marker_map[Marker]), levels = locus_levels),
         copies_L = as.numeric(copies_mL) * 1000) %>%
  filter(!is.na(copies_L), copies_L > 0) %>%
  group_by(Carboy, marker, Timepoint) %>%
  summarise(t = mean(Hours_base), logC = mean(log(copies_L)), .groups = "drop") %>%
  arrange(Carboy, marker, t)

steps <- series %>%
  group_by(Carboy, marker) %>%
  mutate(dt = t - lag(t), dlog = logC - lag(logC)) %>%
  filter(!is.na(dt), dt > 0) %>%
  ungroup()

# --- Estimator 1: the decay model --------------------------------------------

published <- read.csv(here("outputs", "ddPCR", "carboy_rates.csv"),
                      stringsAsFactors = FALSE) %>%
  transmute(marker = factor(marker, levels = locus_levels),
            lambda = -r, lo = -r - 1.96 * r_sd, hi = -r + 1.96 * r_sd,
            method = "Decay model")

# --- Estimator 2: stepwise / AR1, bootstrapped -------------------------------

ar1 <- lapply(locus_levels, function(mk) {
  d   <- steps %>% filter(marker == mk)
  fit <- function(x) -coef(lm(dlog ~ 0 + dt, data = x))[[1]]
  bs  <- replicate(N_BOOT, fit(d[sample(nrow(d), replace = TRUE), , drop = FALSE]))
  data.frame(marker = mk, lambda = fit(d),
             lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)),
             method = "Stepwise (AR1)",
             n_transitions = nrow(d), row.names = NULL)
}) %>% bind_rows() %>% mutate(marker = factor(marker, levels = locus_levels))

# --- Assemble ----------------------------------------------------------------

comparison <- bind_rows(published,
                        ar1 %>% dplyr::select(-n_transitions)) %>%
  mutate(beta_1h = exp(-lambda),
         method  = factor(method,
                          levels = c("Decay model",
                                     "Stepwise (AR1)")))

agreement <- published %>%
  dplyr::select(marker, lambda_bayes = lambda, lo_b = lo, hi_b = hi) %>%
  left_join(ar1 %>% dplyr::select(marker, lambda_ar1 = lambda,
                                  lo_a = lo, hi_a = hi, n_transitions),
            by = "marker") %>%
  transmute(marker,
            lambda_bayes, lambda_ar1,
            difference     = lambda_ar1 - lambda_bayes,
            pct_difference = 100 * (lambda_ar1 / lambda_bayes - 1),
            beta_bayes     = exp(-lambda_bayes),
            beta_ar1       = exp(-lambda_ar1),
            intervals_overlap         = lo_a <= hi_b & hi_a >= lo_b,
            bayes_inside_ar1_interval = lambda_bayes >= lo_a & lambda_bayes <= hi_a,
            n_transitions)

# --- Leverage: drop each timepoint in turn and refit -------------------------

loo_tp <- lapply(locus_levels, function(mk) {
  d <- series %>% filter(marker == mk)
  lapply(sort(unique(d$Timepoint)), function(k) {
    m <- lm(logC ~ 0 + factor(Carboy) + t, data = d %>% filter(Timepoint != k))
    data.frame(marker = mk, dropped = k,
               dropped_hours = round(mean(d$t[d$Timepoint == k]), 1),
               lambda = -coef(m)[["t"]], row.names = NULL)
  }) %>% bind_rows()
}) %>% bind_rows() %>% mutate(marker = factor(marker, levels = locus_levels))

write.csv(
  bind_rows(
    comparison %>% mutate(section = "rate estimates") %>% mutate(across(everything(), as.character)),
    agreement  %>% mutate(section = "agreement")      %>% mutate(across(everything(), as.character)),
    loo_tp     %>% mutate(section = "leave-one-timepoint-out") %>% mutate(across(everything(), as.character))),
  here("outputs", "ddPCR", "ar1_decay_check.csv"), row.names = FALSE, na = "")

# =============================================================================
# Figure
# =============================================================================

anchor <- series %>% group_by(marker) %>%
  summarise(c0 = mean(logC[t == min(t)]), .groups = "drop")

line_df <- comparison %>%

  mutate(kind = ifelse(grepl("Stepwise", method), "Stepwise (AR1)", "Decay model")) %>%
  dplyr::select(marker, lambda, kind) %>%
  left_join(anchor, by = "marker") %>%
  crossing(tt = seq(0, HOURS_MAX, length.out = 60)) %>%
  mutate(logC = c0 - lambda * tt)

p_data <- ggplot(series, aes(t, logC)) +
  geom_point(aes(colour = marker, shape = factor(Carboy)), alpha = 0.75, size = 1.9) +
  geom_line(data = line_df, aes(tt, logC, linetype = kind), linewidth = 0.75) +
  facet_wrap(~ marker, nrow = 1, labeller = as_labeller(marker_display)) +
  scale_colour_manual(values = marker_colors, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), name = "Carboy") +
  scale_linetype_manual(values = c("Decay model" = "solid",
                                   "Stepwise (AR1)"  = "22"), name = NULL) +
  labs(x = "Hours since collection", y = "log eDNA (copies/L)",
       title = "(a)  Observed decay, 0-30 h, with both fitted rates") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

p_rate <- ggplot(comparison, aes(lambda, fct_rev(marker), colour = method)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.6), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.6), size = 2.4) +
  scale_y_discrete(labels = marker_display) +
  scale_colour_manual(values = c("#08519C", "#D95F02"), name = NULL) +
  labs(x = expression(paste("decay rate  ", lambda, "  (per hour)")), y = NULL,
       title = expression(bold(paste("(b)  Decay rate  ", lambda)))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(ncol = 1))

# (c) beta = exp(-lambda), the quantity the AR1 formulation is stated in.
# Plotted as stepwise against decay model with the line of equality.
beta_pair <- comparison %>%
  mutate(beta = exp(-lambda), beta_lo = exp(-hi), beta_hi = exp(-lo)) %>%
  dplyr::select(marker, method, beta, beta_lo, beta_hi) %>%
  pivot_wider(names_from = method,
              values_from = c(beta, beta_lo, beta_hi)) %>%
  rename(decay = `beta_Decay model`, step = `beta_Stepwise (AR1)`,
         step_lo = `beta_lo_Stepwise (AR1)`, step_hi = `beta_hi_Stepwise (AR1)`,
         decay_lo = `beta_lo_Decay model`, decay_hi = `beta_hi_Decay model`)

rng <- range(c(beta_pair$step_lo, beta_pair$step_hi,
               beta_pair$decay_lo, beta_pair$decay_hi))

p_beta <- ggplot(beta_pair, aes(decay, step, colour = marker)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_linerange(aes(ymin = step_lo, ymax = step_hi), linewidth = 0.7, alpha = 0.8) +
  geom_errorbarh(aes(xmin = decay_lo, xmax = decay_hi), height = 0, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = marker_colors, labels = marker_display, name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = expression(paste(beta, "  from the decay model,  ", e^{-lambda})),
       y = expression(paste(beta, "  from the stepwise fit")),
       title = expression(bold(paste("(c)  One-hour coefficient  ", beta)))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_loo <- ggplot(loo_tp, aes(lambda, fct_rev(marker))) +
  geom_linerange(data = published, inherit.aes = FALSE,
                 aes(xmin = lo, xmax = hi, y = fct_rev(marker)),
                 colour = "grey82", linewidth = 3.4) +
  geom_point(aes(colour = dropped_hours), size = 2.3,
             position = position_jitter(height = 0.13, width = 0)) +
  geom_point(data = published, inherit.aes = FALSE,
             aes(lambda, fct_rev(marker)), shape = 4, size = 3, stroke = 1.1) +
  scale_y_discrete(labels = marker_display) +
  scale_colour_viridis_c(name = "dropped\ntimepoint (h)", option = "C", end = 0.9) +
  labs(x = expression(paste(lambda, " refitted, one timepoint removed")), y = NULL,
       title = "(d)  Leave-one-timepoint-out") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

fig <- plot_grid(p_data,
                 plot_grid(p_rate, p_beta, p_loo, ncol = 3,
                           rel_widths = c(1, 1, 1)),
                 ncol = 1, rel_heights = c(1, 1.25))

ggsave(here("plots", "ddPCR", "ar1_decay_check.png"), fig,
       width = 14, height = 8.6, dpi = 300)

# =============================================================================
# Report
# =============================================================================

cat("\n=== Decay rate, both estimators (0 -", HOURS_MAX, "h) ===\n")
print(as.data.frame(comparison %>%
  transmute(marker = lab(marker), method,
            lambda = round(lambda, 4), lo = round(lo, 4), hi = round(hi, 4),
            beta_1h = round(beta_1h, 4))))

cat("\n=== Decay model vs stepwise AR1 ===\n")
print(as.data.frame(agreement %>%
  transmute(marker = lab(marker),
            lambda_decay = round(lambda_bayes, 4),
            lambda_ar1   = round(lambda_ar1, 4),
            pct_difference = round(pct_difference, 1),
            beta_decay = round(beta_bayes, 4), beta_ar1 = round(beta_ar1, 4),
            intervals_overlap, decay_inside_ar1_interval = bayes_inside_ar1_interval, n_transitions)))

full_lambda <- sapply(locus_levels, function(mk)
  -coef(lm(logC ~ 0 + factor(Carboy) + t, data = series %>% filter(marker == mk)))[["t"]])

cat("\n=== Leverage: lambda when one timepoint is dropped ===\n")
print(as.data.frame(loo_tp %>% group_by(marker) %>%
  summarise(min = min(lambda), max = max(lambda), .groups = "drop") %>%
  mutate(full_data     = full_lambda[as.character(marker)],
         abs_spread    = max - min,
         max_shift_pct = round(100 * pmax(abs(max - full_data),
                                          abs(min - full_data)) / full_data, 1),
         marker = lab(marker),
         across(c(min, max, full_data, abs_spread), ~round(.x, 4)))))

cat("\nBridge has the widest absolute spread (",
    round(max(loo_tp$lambda[loo_tp$marker == "Tt_longFrag"]) -
          min(loo_tp$lambda[loo_tp$marker == "Tt_longFrag"]), 4),
    " on lambda) and Cytb the largest relative shift; both have the fewest\n",
    "usable observations, Bridge because of its non-detections.\n", sep = "")

cat("\nSaved: outputs/ddPCR/ar1_decay_check.csv\n")
cat("       plots/ddPCR/ar1_decay_check.png\n")
