# =============================================================================
# 00_setup.R
# -----------------------------------------------------------------------------
# Shared project setup: load packages, set global options, resolve paths,
# and expose reference parameters used across multiple scripts.
#
# Every analysis script sources this file at the top.
# =============================================================================

# --- Packages ----------------------------------------------------------------

library(here)
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(RColorBrewer)

# --- Global options ----------------------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# --- Directory layout --------------------------------------------------------
# here() anchors to the project root (identified by a .Rproj, .here, or
# .git file at the repo root). Outputs go to /plots and /outputs.

dir.create(here("plots"),   showWarnings = FALSE, recursive = TRUE)
dir.create(here("outputs"), showWarnings = FALSE, recursive = TRUE)

# --- Stan model files --------------------------------------------------------

stan_conc_sim   <- here("code", "conc_age_sim.stan")    # tight priors  (sim + LOO)
stan_conc_field <- here("code", "conc_age_field.stan")  # wider priors (field data)
stan_decay      <- here("code", "decay_monophasic_3.stan")

# --- Marker reference values -------------------------------------------------
# Locus-specific decay rates (r) and log-concentration offsets at t = 0 (p).
# These are the posterior means from fitting the mono-phasic decay model to
# all three carboys jointly (see 02_decay.R). Field scripts pull these
# from outputs/carboy_rates.csv once 02_decay.R has been run; if the file
# is not yet present we fall back to the published values.

carboy_rates_file <- here("outputs", "carboy_rates.csv")

if (file.exists(carboy_rates_file)) {
  .cr         <- read.csv(carboy_rates_file, stringsAsFactors = FALSE)
  r_by_marker <- setNames(.cr$r, .cr$marker)
  p_by_marker <- setNames(.cr$p, .cr$marker)
  rm(.cr)
} else {
  # Published values (reproduce 02_decay.R to regenerate).
  r_by_marker <- c(cytb        = -0.1141,
                   Tt_16S      = -0.1649,
                   Tt_DLL1     = -0.1667,
                   Tt_longFrag = -0.1915)
  p_by_marker <- c(cytb        =  1.0000,
                   Tt_16S      =  0.7957,
                   Tt_DLL1     =  0.4096,
                   Tt_longFrag =  0.1029)
}

locus_levels <- names(r_by_marker)

# --- Display names and palette for field markers -----------------------------

marker_display <- c(
  cytb        = "Cytb",
  Tt_16S      = "16S",
  Tt_DLL1     = "D-loop",
  Tt_longFrag = "Bridge"
)

marker_colors        <- RColorBrewer::brewer.pal(4, "Set2")[c(2, 1, 4, 3)]
names(marker_colors) <- locus_levels

station_colors <- c("Husbandry Area" = "#53CFDA",
                    "NOAA Boat"      = "#FF7994")

station_labels <- c("Husbandry Area" = "Nearby station (HUSB)",
                    "NOAA Boat"      = "Further station (BOAT)")
