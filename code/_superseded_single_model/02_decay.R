# =============================================================================
# 02_decay.R
# -----------------------------------------------------------------------------
# Carboy decay experiment: fits the mono-phasic decay model to the ddPCR
# data and uses it two ways:
#
#   Part A. Fit to ALL THREE carboys jointly to obtain the locus-specific
#           decay rates (r) and log-concentration offsets at t = 0 (p) that
#           the field analyses require.
#           -> writes outputs/carboy_rates.csv
#
#   Part B. Leave-one-out validation: for each carboy in turn, fit the
#           decay model on the other two, then use those rates with the
#           age-estimation model (conc_age_sim.stan) to recover the
#           elapsed time of the held-out carboy's samples.
#           -> plots/carboy_loo_grid.png
# =============================================================================

source(here::here("code", "00_setup.R"))

# --- Input -------------------------------------------------------------------

Final_decay_ddPCR_datasheet <- read.csv(
  here("data", "Final_decay_ddPCR_datasheet.csv"),
  stringsAsFactors = FALSE
)

# --- Settings (shared across A and B) ----------------------------------------

components_use      <- "DNA"
hours_max_decay     <- 30
hours_max_time      <- 30
t_mean_prior        <- 12
t_sd_prior          <- 24

locus_levels_master <- c("DNA_Cytb", "DNA_16S", "DNA_Dloop", "DNA_Bridge")
ref_locus           <- locus_levels_master[1]

# Mapping: carboy locus names -> field-script marker names, in the order
# expected by the field analyses (must match names(r_by_marker) in 00_setup.R).
locus_to_marker <- c(
  DNA_Cytb   = "cytb",
  DNA_16S    = "Tt_16S",
  DNA_Dloop  = "Tt_DLL1",
  DNA_Bridge = "Tt_longFrag"
)

# =============================================================================
# Helpers (shared by Part A and Part B)
# =============================================================================

# Assign 1..K integer indices to the sorted unique values of a column.
make_index <- function(df, col, new_col) {
  vals <- sort(unique(df[[col]]))
  df[[new_col]] <- match(df[[col]], vals)
  df
}

# Combine multiple index columns into a single sequential index.
combine_index <- function(df, cols, new_col) {
  key <- df %>%
    transmute(.key = do.call(paste, c(across(all_of(cols)), sep = "||")))
  lev <- unique(key$.key)
  map <- tibble(.key = lev, !!as.name(new_col) := seq_along(lev))
  df %>% bind_cols(key) %>% left_join(map, by = ".key") %>% select(-.key)
}

# Build the data list that decay_monophasic_3.stan consumes.
build_decay_data <- function(dat, locus_levels_ord) {
  dat <- dat %>%
    mutate(
      locus          = factor(paste(Component, Marker, sep = "_"),
                              levels = locus_levels_ord),
      k_idx          = as.integer(locus),
      Hours_base     = as.numeric(Hours_base),
      Total_droplets = Positive_droplets + Negative_droplets
    ) %>%
    filter(!is.na(k_idx), !is.na(Negative_droplets)) %>%
    make_index("Carboy",    "i_idx") %>%
    make_index("Timepoint", "j_idx") %>%
    combine_index(c("k_idx", "i_idx"), "ik_idx")

  df0 <- dat %>% filter(Hours_base == 0)
  df1 <- dat %>%
    filter(Hours_base > 0) %>%
    combine_index(c("k_idx", "i_idx", "j_idx"), "ijk_idx") %>%
    arrange(ijk_idx)

  secondary_idx <- df1 %>% distinct(ijk_idx, k_idx, ik_idx, Hours_base)
  ik_map        <- dat %>% distinct(ik_idx, k_idx) %>% arrange(ik_idx)

  list(
    N_obs_0    = nrow(df0),
    N_obs_1    = nrow(df1),
    N_ijk_1    = n_distinct(df1$ijk_idx),
    N_ik_0     = n_distinct(df0$ik_idx),
    N_k_1      = n_distinct(df1$k_idx),
    ik_idx_0   = as.integer(df0$ik_idx),
    ijk_idx_1  = as.integer(df1$ijk_idx),
    s_ik_idx_1 = as.integer(secondary_idx$ik_idx),
    s_k_idx_1  = as.integer(secondary_idx$k_idx),
    ik_to_k    = as.integer(ik_map$k_idx),
    time       = as.numeric(secondary_idx$Hours_base),
    Dilution_0 = as.numeric(df0$Dilution),
    Dilution_1 = as.numeric(df1$Dilution),
    Filtered_0 = as.numeric(df0$Filt_mL) / 1000,
    Filtered_1 = as.numeric(df1$Filt_mL) / 1000,
    W_0        = as.integer(df0$Positive_droplets),
    W_1        = as.integer(df1$Positive_droplets),
    U_0        = as.integer(df0$Total_droplets),
    U_1        = as.integer(df1$Total_droplets),
    N_time_sim = length(unique(df1$Hours_base)),
    time_sim   = sort(unique(as.numeric(df1$Hours_base)))
  )
}

# Compute (p, C0_mean, C0_sd) from the training-set t=0 observations.
compute_p_and_C0 <- function(dat_train_t0, locus_levels_ord) {
  t0_summ <- dat_train_t0 %>%
    mutate(
      locus     = paste(Component, Marker, sep = "_"),
      copies_mL = as.numeric(copies_mL),
      logC      = if_else(!is.na(copies_mL) & copies_mL > 0,
                          log(copies_mL), NA_real_)
    ) %>%
    filter(locus %in% locus_levels_ord, !is.na(logC)) %>%
    mutate(locus = factor(locus, levels = locus_levels_ord)) %>%
    group_by(locus) %>%
    summarise(mean_logC = mean(logC, na.rm = TRUE),
              sd_logC   = sd(logC,   na.rm = TRUE),
              .groups   = "drop") %>%
    arrange(locus)

  ref_row  <- t0_summ %>% filter(locus == locus_levels_ord[1])
  logC_ref <- ref_row$mean_logC

  list(
    p       = t0_summ$mean_logC - logC_ref,
    C0_mean = logC_ref,
    C0_sd   = max(ref_row$sd_logC, 1.5, na.rm = TRUE)
  )
}

# Build the data list for conc_age_sim.stan from a held-out carboy.
build_age_stan_data <- function(dat_test, locus_levels_ord, r_vec, p_vec,
                                C0_mean, C0_sd, t_mean, t_sd) {
  Nloci <- length(locus_levels_ord)

  dat_prep <- dat_test %>%
    mutate(
      locus      = factor(paste(Component, Marker, sep = "_"),
                          levels = locus_levels_ord),
      j          = as.integer(locus),
      copies_mL  = as.numeric(copies_mL),
      Hours_base = as.numeric(Hours_base),
      y          = if_else(!is.na(copies_mL) & copies_mL > 0,
                           log(copies_mL), NA_real_),
      z          = as.integer(!is.na(y))
    ) %>%
    filter(!is.na(j))

  tp_map <- dat_prep %>%
    group_by(Timepoint) %>%
    summarise(True_Hours = median(Hours_base, na.rm = TRUE), .groups = "drop") %>%
    arrange(True_Hours) %>%
    mutate(i = row_number())

  dat_prep <- dat_prep %>%
    left_join(tp_map %>% select(Timepoint, i), by = "Timepoint")

  det_rows <- dat_prep %>% filter(z == 1)

  list(
    stan_data = list(
      Nt      = as.integer(nrow(tp_map)),
      Nloci   = as.integer(Nloci),
      r       = r_vec,
      p       = p_vec,
      N_obs   = nrow(dat_prep),
      obs_i   = as.integer(dat_prep$i),
      obs_j   = as.integer(dat_prep$j),
      z       = as.integer(dat_prep$z),
      N_y     = nrow(det_rows),
      y_i     = as.integer(det_rows$i),
      y_j     = as.integer(det_rows$j),
      y_obs   = as.numeric(det_rows$y),
      C0_mean = C0_mean,
      C0_sd   = C0_sd,
      t_mean  = t_mean,
      t_sd    = t_sd
    ),
    tp_map = tp_map
  )
}

# =============================================================================
# Part A. Full fit on all three carboys -> r and p for field analyses
# =============================================================================

message("\n=== Part A: decay fit on all three carboys ===")

dat_all <- Final_decay_ddPCR_datasheet %>%
  filter(!Control,
         Component %in% components_use,
         Hours_base <= hours_max_decay,
         Carboy %in% c(1, 2, 3))

decay_sd_all <- build_decay_data(dat_all, locus_levels_master)

init_decay_all <- function() list(
  C_0    = rep(10,   decay_sd_all$N_ik_0),
  lambda = rep(0.01, decay_sd_all$N_k_1)
)

fit_decay_all <- stan(
  file    = stan_decay,
  data    = decay_sd_all,
  init    = init_decay_all,
  chains  = 4,
  iter    = 3000,
  warmup  = 1000,
  seed    = 99,
  control = list(adapt_delta = 0.95),
  refresh = 300
)

# --- Extract r (decay rates, expressed as negative for downstream use) ------

lambda_hat   <- summary(fit_decay_all, pars = "lambda")$summary[, "mean"]
r_raw        <- -as.numeric(lambda_hat)
names(r_raw) <- locus_levels_master

# Rename to field-marker names, preserving order
r_by_marker_fit        <- r_raw[names(locus_to_marker)]
names(r_by_marker_fit) <- locus_to_marker[names(locus_to_marker)]

# --- Extract p (t=0 mean concentration, normalised to the reference locus) ---

dat_t0 <- dat_all %>%
  filter(Hours_base == 0) %>%
  mutate(locus     = paste(Component, Marker, sep = "_"),
         copies_mL = as.numeric(copies_mL)) %>%
  filter(locus %in% locus_levels_master, !is.na(copies_mL), copies_mL > 0)

mean_conc_by_locus <- dat_t0 %>%
  mutate(locus = factor(locus, levels = locus_levels_master)) %>%
  group_by(locus) %>%
  summarise(mean_conc = mean(copies_mL, na.rm = TRUE), .groups = "drop") %>%
  arrange(locus)

conc_cytb          <- mean_conc_by_locus$mean_conc[1]
p_raw              <- mean_conc_by_locus$mean_conc / conc_cytb
names(p_raw)       <- locus_levels_master

p_by_marker_fit        <- p_raw[names(locus_to_marker)]
names(p_by_marker_fit) <- locus_to_marker[names(locus_to_marker)]

# --- Save to /outputs/carboy_rates.csv (consumed by 00_setup.R) -------------

carboy_rates <- data.frame(
  marker = names(r_by_marker_fit),
  r      = as.numeric(r_by_marker_fit),
  p      = as.numeric(p_by_marker_fit)
)

write.csv(carboy_rates, here("outputs", "carboy_rates.csv"), row.names = FALSE)
message("Saved: outputs/carboy_rates.csv")
print(round(as.data.frame(carboy_rates[-1]), 4))

# =============================================================================
# Part B. Leave-one-out validation
# =============================================================================

message("\n=== Part B: leave-one-out carboy validation ===")

folds <- list(
  list(train = c(1, 2), test = 3, label = "1+2"),
  list(train = c(1, 3), test = 2, label = "1+3"),
  list(train = c(2, 3), test = 1, label = "2+3")
)

pred_list <- list()

for (f in folds) {

  message("\n--- Fold: train = ", paste(f$train, collapse = "+"),
          " | test = ", f$test, " ---")

  # Training data (two carboys) -> r vector and (p, C0)
  dat_train <- Final_decay_ddPCR_datasheet %>%
    filter(!Control,
           Component %in% components_use,
           Hours_base <= hours_max_decay,
           Carboy %in% f$train)

  decay_sd <- build_decay_data(dat_train, locus_levels_master)

  fit_decay <- stan(
    file   = stan_decay,
    data   = decay_sd,
    init   = function() list(C_0    = rep(10,   decay_sd$N_ik_0),
                             lambda = rep(0.01, decay_sd$N_k_1)),
    chains = 2,
    iter   = 2500,
    warmup = 500,
    seed   = 42 + f$test
  )

  lambda_hat <- summary(fit_decay, pars = "lambda")$summary[, "mean"]
  r_vec      <- -as.numeric(lambda_hat)

  pc <- compute_p_and_C0(dat_train %>% filter(Hours_base == 0),
                         locus_levels_master)

  # Test data (held-out carboy) -> age estimation
  dat_test <- Final_decay_ddPCR_datasheet %>%
    filter(!Control,
           Component %in% components_use,
           Carboy     == f$test,
           Hours_base <  hours_max_time)

  built <- build_age_stan_data(
    dat_test         = dat_test,
    locus_levels_ord = locus_levels_master,
    r_vec            = r_vec,
    p_vec            = pc$p,
    C0_mean          = pc$C0_mean,
    C0_sd            = pc$C0_sd,
    t_mean           = t_mean_prior,
    t_sd             = t_sd_prior
  )

  fit_time <- stan(
    file    = stan_conc_sim,
    data    = built$stan_data,
    chains  = 4,
    iter    = 2000,
    warmup  = 1000,
    seed    = 100 + f$test,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 200
  )

  t_sum <- summary(fit_time, pars = "t",
                   probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary

  pred_list[[f$label]] <- tibble(i = seq_len(nrow(t_sum))) %>%
    bind_cols(as.data.frame(t_sum)) %>%
    rename(lo95 = `2.5%`, lo50 = `25%`, med = `50%`,
           hi50 = `75%`, hi95 = `97.5%`) %>%
    left_join(built$tp_map, by = "i") %>%
    mutate(Carboy = f$test, Fold = f$label)
}

pred_all <- bind_rows(pred_list) %>%
  mutate(Timepoint = as.numeric(Timepoint))

# --- Panel (a): observed eDNA decay across carboys ---------------------------

wanted_loci <- tibble::tribble(
  ~Marker,  ~Component,
  "Cytb",   "DNA",
  "16S",    "DNA",
  "Dloop",  "DNA",
  "Bridge", "DNA"
)

decay_plotdat <- Final_decay_ddPCR_datasheet %>%
  filter(!Control, Component %in% components_use, Hours_base <= 30) %>%
  inner_join(wanted_loci, by = c("Marker", "Component")) %>%
  mutate(
    copies_mL     = as.numeric(copies_mL),
    locus         = paste(Marker, Component, sep = "_"),
    Carboy        = factor(Carboy),
    log_copies_mL = if_else(!is.na(copies_mL) & copies_mL > 0,
                            log(copies_mL), NA_real_)
  )

set2_cols        <- brewer.pal(4, "Set2")[c(1, 2, 4, 3)]
names(set2_cols) <- c("16S_DNA", "Cytb_DNA", "Bridge_DNA", "Dloop_DNA")
set2_cols_use    <- set2_cols[names(set2_cols) %in% unique(decay_plotdat$locus)]

p_carboy_decay <- ggplot(decay_plotdat,
                         aes(x = Hours_base, y = copies_mL * 1000,
                             color = locus, shape = Carboy)) +
  geom_point(alpha = 0.6) +
  geom_smooth(aes(group = locus), method = "lm", se = FALSE, na.rm = TRUE) +
  scale_color_manual(
    values = set2_cols_use,
    breaks = c("Cytb_DNA", "16S_DNA", "Dloop_DNA", "Bridge_DNA"),
    labels = c(
      "Cytb_DNA"   = "Cytb (\u03bb = -0.114)",
      "16S_DNA"    = "16S (\u03bb = -0.165)",
      "Dloop_DNA"  = "D-loop (\u03bb = -0.167)",
      "Bridge_DNA" = "Bridge (\u03bb = -0.192)"
    ),
    name = "Marker"
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Carboy") +
  labs(title = "Observed eDNA decay across carboys",
       x     = "Time elapsed (hours)",
       y     = "eDNA concentration (copies/L)") +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 30)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# --- Panel (b): leave-one-out time estimation -------------------------------

p_loo_time <- ggplot(pred_all, aes(x = True_Hours, color = factor(Carboy))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_linerange(aes(y = med, ymin = lo95, ymax = hi95),
                 linewidth = 0.6, alpha = 0.5) +
  geom_linerange(aes(y = med, ymin = lo50, ymax = hi50),
                 linewidth = 1.4) +
  geom_point(aes(y = med), size = 2.5) +
  scale_color_manual(
    values = c("1" = "#C6DBEF", "2" = "#6BAED6", "3" = "#08519C"),
    name   = "Carboy"
  ) +
  ylim(0, max(pred_all$hi95, na.rm = TRUE)) +
  xlim(0, 30) +
  labs(title = "Leave-one-out age estimations",
       x     = "Actual elapsed time (hours)",
       y     = "Estimated elapsed time (hours)",
       color = "Carboy") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# --- Compose and save --------------------------------------------------------

carboy_loo_grid <- plot_grid(
  p_carboy_decay, p_loo_time,
  labels     = c("(a)", "(b)"),
  label_size = 16,
  ncol       = 2,
  rel_widths = c(0.55, 0.45)
)
carboy_loo_grid 
ggsave(here("plots", "carboy_loo_grid.png"),
       plot = carboy_loo_grid, width = 12, height = 5, dpi = 300)

message("Saved: plots/carboy_loo_grid.png")
