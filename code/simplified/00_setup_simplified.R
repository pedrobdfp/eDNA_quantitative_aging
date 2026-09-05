# =============================================================================
# 00_setup_simplified.R
# -----------------------------------------------------------------------------
# Shared setup for the simplified (concentration) analysis presented in the
# supplement. Every other script in code/simplified/ sources this file first.
#
# The simplified model reads concentrations and detection flags rather than
# droplet counts, so it applies to qPCR as readily as to ddPCR. Everything
# above the observation model -- the decay process, the marker offsets, the
# priors on age and baseline concentration -- is identical to the ddPCR
# version in code/ddPCR/.
# =============================================================================

library(here)
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(RColorBrewer)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# --- Output locations --------------------------------------------------------

MODEL_TAG <- "simplified"

dir.create(here("plots",   MODEL_TAG), showWarnings = FALSE, recursive = TRUE)
dir.create(here("outputs", MODEL_TAG), showWarnings = FALSE, recursive = TRUE)

out_path  <- function(...) here("outputs", MODEL_TAG, ...)
plot_path <- function(...) here("plots",   MODEL_TAG, ...)

data_path <- function(fname) {
  f <- here("data", fname)
  if (!file.exists(f)) stop("Could not find '", fname, "' in data/.")
  f
}

# --- Concentration scale -----------------------------------------------------
# Every concentration in this analysis is copies per LITRE of water, and every
# model quantity is its natural log. Choosing one volumetric scale and staying
# on it is the only convention the simplified model imposes; if your data are
# copies per mL, multiply by 1000 before you start.

CONC_UNITS <- "copies per litre of water"

# --- Observation noise -------------------------------------------------------
# Observation error is split into two nested components, both shared across
# markers:
#
#   sigma_bio    SD between water samples of the same unit
#   sigma_tech   SD between PCR wells of the same water sample
#
# Both are shared rather than marker-specific. A marker-specific SD would shift
# the marker pattern, and the marker pattern is what identifies the age.

SIGMA_SD <- 1        # half-normal prior scale for both components

# Noise levels used by 02_simulation_simplified.R, measured directly in the
# field samples: tech is the SD between wells of one extract, bio_field the SD
# between water samples of the same grab, and bio_low / bio_high bracket it.
SIGMA_LEVELS <- list(tech      = 0.11,
                     bio_low   = 0.10,
                     bio_field = 0.60,
                     bio_high  = 1.20)

# --- Detection function ------------------------------------------------------
# Detection is modelled as
#
#     P(detect) = logit^-1(beta * (log C - logC50))
#
# parameterised by the two quantities worth having an opinion about: logC50,
# the log-concentration at which half of replicates amplify, and beta, the
# slope there. The equivalent intercept alpha = -beta * logC50 is reported by
# the model. Written this way the two parameters are close to independent; in
# the intercept-and-slope form they are strongly correlated and the sampler is
# several times slower.
#
# For these data one positive droplet in a 20,000-droplet well of a 3 L filtrate
# corresponds to about 10 copies per litre, so the limit of detection is put
# there, loosely.

DETECT_LOGC50    <- log(10)      # prior centre for the 50% detection point
DETECT_LOGC50_SD <- 1.5          # how loosely, on the log scale
DETECT_SLOPE     <- c(mean = 2, sd = 1)

detection_priors <- function(logC50    = DETECT_LOGC50,
                             logC50_sd = DETECT_LOGC50_SD,
                             slope     = DETECT_SLOPE) {
  list(logC50_mean = logC50,
       logC50_sd   = logC50_sd,
       beta_mean   = slope[["mean"]],
       beta_sd     = slope[["sd"]])
}

# --- Stan models -------------------------------------------------------------

stan_conc_age <- here("code", "simplified", "conc_age_simplified.stan")
stan_decay    <- here("code", "simplified", "decay_simplified.stan")

for (f in c(stan_conc_age, stan_decay)) {
  if (!file.exists(f)) stop(basename(f), " not found in code/simplified/.")
}

# --- Field data --------------------------------------------------------------
# The simplified analysis uses exactly the same wells as the ddPCR analysis, so
# that any difference between the two is attributable to the observation model
# and nothing else. Each well is reduced to two numbers: whether the marker was
# detected, and, if it was, the concentration reported by the instrument.
#
# data/field_droplets.csv reports conc_machine, the concentration in the
# reaction (copies/uL). Converting it to seawater concentration uses the same
# volumetric offset the droplet model uses:
#
#   log C[copies/L] = log(conc_machine) - log_offset + log(v_droplet)
#
# For a qPCR dataset there is no such step: the standard curve already reports
# a concentration, and only the filtered volume has to be divided out.

V_DROPLET_UL <- 0.00085

load_field_concentrations <- function() {
  d <- read.csv(data_path("field_droplets.csv"), stringsAsFactors = FALSE)

  d$detected <- as.integer(d$positives > 0)
  d$logC     <- ifelse(d$detected == 1,
                       log(d$conc_machine) - d$log_offset + log(V_DROPLET_UL),
                       NA_real_)
  d
}

# --- Marker decay rates and t = 0 offsets ------------------------------------
# Estimated by 01_decay_simplified.R under the simplified observation model and
# treated as known by the field analyses. The ddPCR track has its own estimates;
# the two are close but not identical, and each track uses its own throughout.

carboy_rates_file <- out_path("carboy_rates.csv")

if (file.exists(carboy_rates_file)) {
  .cr         <- read.csv(carboy_rates_file, stringsAsFactors = FALSE)
  r_by_marker <- setNames(.cr$r, .cr$marker)
  p_by_marker <- setNames(.cr$p, .cr$marker)
  rm(.cr)
} else {
  message("00_setup_simplified.R: ", basename(carboy_rates_file), " not found. ",
          "Run 01_decay_simplified.R to generate it; ddPCR values used below.")

  r_by_marker <- c(cytb        = -0.114074,
                   Tt_16S      = -0.164945,
                   Tt_DLL1     = -0.166744,
                   Tt_longFrag = -0.191465)

  p_by_marker <- c(cytb        =  0.000000,
                   Tt_16S      = -0.220540,
                   Tt_DLL1     = -0.886862,
                   Tt_longFrag = -2.343509)
}

locus_levels <- names(r_by_marker)

if (!isTRUE(all.equal(unname(p_by_marker[1]), 0))) {
  stop("p_by_marker[1] must be 0 on the log scale but is ", p_by_marker[1], ".")
}

# --- Display names and palette -----------------------------------------------

marker_display <- c(cytb        = "Cytb",
                    Tt_16S      = "16S",
                    Tt_DLL1     = "D-loop",
                    Tt_longFrag = "Bridge")

marker_colors        <- RColorBrewer::brewer.pal(4, "Set2")[c(2, 1, 4, 3)]
names(marker_colors) <- locus_levels

station_colors <- c("Husbandry Area" = "#53CFDA",
                    "NOAA Boat"      = "#FF7994")

station_labels <- c("Husbandry Area" = "Nearby station (HUSB)",
                    "NOAA Boat"      = "Further station (BOAT)")

# --- Stan data assembly ------------------------------------------------------
# This is the entry point for your own data. `obs` needs one row per PCR
# replicate and these five columns:
#
#   obs_i     integer, which unit is being aged (a sample, or a pooled grab)
#   obs_j     integer, which marker, indexed the same way as r_vec and p_vec
#   detected  0 or 1
#   logC      log concentration of the water, natural log of copies per litre;
#             ignored where detected == 0, so NA there is fine
#   bio       optional; which water sample the replicate came from. Supply it
#             when a unit contains several water samples, leave it NULL when
#             replication is technical only.

build_conc_stan_data <- function(obs, Nt, Nloci, r_vec, p_vec,
                                 C0_mean, C0_sd, t_mean, t_sd,
                                 bio = NULL, sigma_sd = SIGMA_SD,
                                 detect = detection_priors()) {
  stopifnot(all(c("obs_i", "obs_j", "detected", "logC") %in% names(obs)))
  stopifnot(all(obs$detected %in% c(0L, 1L)))
  det_rows <- which(obs$detected == 1)
  stopifnot(all(is.finite(obs$logC[det_rows])))

  use_bio <- !is.null(bio)
  if (use_bio) {
    bio_idx <- match(bio, unique(bio))
    N_bio   <- max(bio_idx)
  } else {
    bio_idx <- rep(1L, nrow(obs))
    N_bio   <- 1L
  }

  list(
    Nt         = as.integer(Nt),
    Nloci      = as.integer(Nloci),
    r          = as.numeric(r_vec),
    p          = as.numeric(p_vec),
    N          = nrow(obs),
    obs_i      = as.integer(obs$obs_i),
    obs_j      = as.integer(obs$obs_j),
    z          = as.integer(obs$detected),
    N_y        = length(det_rows),
    y_row      = as.integer(det_rows),
    y_obs      = as.numeric(obs$logC[det_rows]),
    use_bio    = as.integer(use_bio),
    N_bio      = as.integer(N_bio),
    bio_idx    = as.integer(bio_idx),
    C0_mean    = C0_mean,
    C0_sd      = C0_sd,
    t_mean     = t_mean,
    t_sd       = t_sd,
    logC50_mean = detect$logC50_mean,
    logC50_sd   = detect$logC50_sd,
    beta_mean   = detect$beta_mean,
    beta_sd     = detect$beta_sd,
    sigma_sd   = sigma_sd
  )
}

message("Setup loaded: simplified concentration model, nested observation error.")
print(round(data.frame(r = r_by_marker, p = p_by_marker), 4))
