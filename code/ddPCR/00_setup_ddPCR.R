# =============================================================================
# 00_setup_ddPCR.R
# -----------------------------------------------------------------------------
# Shared setup for the ddPCR droplet analysis presented in the main text.
# Every other script in code/ddPCR/ sources this file first.
#
# It loads packages, resolves paths, defines the constants of the ddPCR
# observation model, and supplies the marker-specific decay rates and t = 0
# offsets that the field analyses treat as known.
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

MODEL_TAG <- "ddPCR"

dir.create(here("plots",   MODEL_TAG), showWarnings = FALSE, recursive = TRUE)
dir.create(here("outputs", MODEL_TAG), showWarnings = FALSE, recursive = TRUE)

out_path  <- function(...) here("outputs", MODEL_TAG, ...)
plot_path <- function(...) here("plots",   MODEL_TAG, ...)

data_path <- function(fname) {
  f <- here("data", fname)
  if (!file.exists(f)) stop("Could not find '", fname, "' in data/.")
  f
}

# --- ddPCR observation-model constants ---------------------------------------
# A well reports W positive droplets out of U accepted droplets. The expected
# number of target copies per droplet is
#
#     lambda = C * vol * dilution * v_droplet / K
#
# where C is the seawater concentration (copies/L), vol the volume filtered (L),
# v_droplet the effective droplet volume (uL), and K the factor relating the
# concentration in the reaction to the concentration in seawater. On the log
# scale this is an additive offset, so that
#
#     log(lambda) = log(C) + droplet_log_offset(vol, dilution).
#
# Both constants were verified against the plate exports in
# 01_build_field_droplets.R: K is exactly 550 for every record, and the droplet
# volume implied by -log(1 - W/U) / conc has a median of 0.00085 uL.

K_CONC       <- 550
V_DROPLET_UL <- 0.00085

droplet_log_offset <- function(vol_L, dilution = 1) {
  log(vol_L) + log(dilution) - log(K_CONC) + log(V_DROPLET_UL)
}

# --- Observation noise -------------------------------------------------------
# Observation error is split into two nested components, both shared across
# markers:
#
#   sigma_bio    SD between water samples of the same grab
#   sigma_tech   SD between wells of the same water sample
#
# Both are shared rather than marker-specific. Because the noise is lognormal,
# the concentration the data identify depends on sigma through a term
# sigma^2 / 2. A marker-specific sigma would make that term differ by marker,
# which is a shift in the marker pattern -- and the marker pattern is what
# identifies the age. A shared sigma makes the correction a single constant
# that the baseline concentration absorbs, so observation error widens the age
# posterior without displacing it.

SIGMA_SD <- 1        # half-normal prior scale for both components

# --- Stan models -------------------------------------------------------------

stan_droplet_age <- here("code", "ddPCR", "droplet_age_ddPCR.stan")
stan_decay       <- here("code", "ddPCR", "decay_monophasic_ddPCR.stan")

for (f in c(stan_droplet_age, stan_decay)) {
  if (!file.exists(f)) stop(basename(f), " not found in code/ddPCR/.")
}

# --- Droplet-level field data ------------------------------------------------
# data/field_droplets.csv is the analysed well-level table: one row per ddPCR
# well, with positive and accepted droplet counts, the volume filtered and the
# sample metadata. It is written by 01_build_field_droplets.R from the raw
# plate exports and is shipped here so the field analyses can be reproduced
# without them.

load_field_droplets <- function() {
  d <- read.csv(data_path("field_droplets.csv"), stringsAsFactors = FALSE)
  stopifnot(all(c("sample", "marker", "positives", "accepted",
                  "log_offset") %in% names(d)))
  d
}

# --- Marker decay rates and t = 0 offsets ------------------------------------
# r[j] is the first-order decay rate of marker j (1/hour, negative) and p[j] its
# log-concentration offset at t = 0 relative to the reference marker (cytb), so
# p[cytb] = 0. Both are estimated from the carboy experiment by 02_decay_ddPCR.R
# and treated as known by the field analyses. The values below are the published
# ones and are used if that script has not been run.

carboy_rates_file <- out_path("carboy_rates.csv")

if (file.exists(carboy_rates_file)) {
  .cr         <- read.csv(carboy_rates_file, stringsAsFactors = FALSE)
  r_by_marker <- setNames(.cr$r, .cr$marker)
  p_by_marker <- setNames(.cr$p, .cr$marker)
  rm(.cr)
} else {
  message("00_setup_ddPCR.R: ", basename(carboy_rates_file), " not found. ",
          "Run 02_decay_ddPCR.R to regenerate it; published values used below.")

  r_by_marker <- c(cytb        = -0.118916,
                   Tt_16S      = -0.161128,
                   Tt_DLL1     = -0.177387,
                   Tt_longFrag = -0.189234)

  p_by_marker <- c(cytb        =  0.000000,
                   Tt_16S      = -0.220540,
                   Tt_DLL1     = -0.886862,
                   Tt_longFrag = -2.343509)
}

locus_levels <- names(r_by_marker)

# p is a log offset, so the reference marker must be exactly zero.
if (!isTRUE(all.equal(unname(p_by_marker[1]), 0))) {
  stop("p_by_marker[1] must be 0 on the log scale but is ", p_by_marker[1], ".")
}

# --- Observation-noise levels used by the simulation -------------------------
# 03_simulation_ddPCR.R generates data at these values, all taken from the real
# samples so that the simulated regime is the one the method is applied in.
#
#   tech        SD between wells of one extract on a single plate, measured
#               directly in the field data; ddPCR precision alone
#   bio_low     a design in which replicate water samples agree closely
#   bio_field   SD between water samples of the same grab, measured directly in
#               the field data. The posterior sigma_bio from the grab-level fit
#               is somewhat larger (about 0.81) because it also absorbs
#               departures of real samples from a single-source decay curve;
#               bio_high brackets that.
#   bio_high    roughly the spread between wells pooled across the two plates
#               on which 15 samples were re-run, i.e. sample-to-sample
#               variation with run-to-run variation added

SIGMA_LEVELS <- list(tech      = 0.11,
                     bio_low   = 0.10,
                     bio_field = 0.60,
                     bio_high  = 1.20)

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
# obs must contain obs_i (unit index), obs_j (marker index), positives,
# accepted and log_offset.

build_droplet_stan_data <- function(obs, Nt, Nloci, r_vec, p_vec,
                                    C0_mean, C0_sd, t_mean, t_sd,
                                    bio = NULL, sigma_sd = SIGMA_SD) {
  stopifnot(all(c("obs_i", "obs_j", "positives", "accepted", "log_offset")
                %in% names(obs)))
  stopifnot(all(obs$positives <= obs$accepted), all(obs$accepted > 0))
  stopifnot(!any(is.na(obs$log_offset)))

  # bio identifies the biological replicate (water sample) each well came from.
  # Supplying it enables the nested level; omitting it fits technical error
  # only, which is correct for designs without biological replication such as
  # the carboy experiment.
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
    W          = as.integer(obs$positives),
    U          = as.integer(obs$accepted),
    log_offset = as.numeric(obs$log_offset),
    use_bio    = as.integer(use_bio),
    N_bio      = as.integer(N_bio),
    bio_idx    = as.integer(bio_idx),
    C0_mean    = C0_mean,
    C0_sd      = C0_sd,
    t_mean     = t_mean,
    t_sd       = t_sd,
    sigma_sd   = sigma_sd
  )
}

message("Setup loaded: ddPCR droplet model, nested observation error.")
print(round(data.frame(r = r_by_marker, p = p_by_marker), 4))
