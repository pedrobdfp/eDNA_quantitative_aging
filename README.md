# Quantitative aging of environmental DNA

Code and data accompanying the manuscript *[Quantitative aging of environmental DNA using multiple components]* submitted to
*Methods in Ecology and Evolution*.

This repository reproduces every figure and parameter estimate in the main
text and supplementary material.

## Repository layout

```
.
├── code/                       analysis scripts and Stan models
│   ├── 00_setup.R              packages, paths, shared constants
│   ├── 01_simulation.R         simulation study (Figures: simulation_main, simulation_supplement)
│   ├── 02_decay.R              carboy decay fit (all 3) + leave-one-out validation
│   ├── 03_field_main.R         field data, grab-averaged -- MAIN field figure
│   ├── 04_field_supplement.R   field data, individual samples -- SUPPLEMENTARY
│   ├── 05_field_removal.R      field data with physical removal (rho = -0.45 /hr)
│   ├── 06_map.R                study-area map (satellite imagery)
│   ├── conc_age_sim.stan       age-estimation model (tight priors; sim + LOO)
│   ├── conc_age_field.stan     age-estimation model (wider priors; field data)
│   └── decay_monophasic_3.stan mono-phasic decay model (carboy experiment)
├── data/                       input data (not redistributed here; see below)
├── plots/                      figure outputs (written by scripts)
└── outputs/                    tabular outputs (written by scripts)
```

## Dependencies

All scripts run under R (tested on R ≥ 4.3). Required packages:

```r
install.packages(c(
  "here", "rstan", "ggplot2", "dplyr", "tidyr", "tibble",
  "cowplot", "RColorBrewer",
  # map only:
  "sf", "terra", "maptiles", "ggspatial"
))
```

A working `rstan` / Stan toolchain is required. See
<https://mc-stan.org/users/interfaces/rstan> for installation.

## Data files

Place these files in `data/`:

| File | Produced by | Used in |
|---|---|---|
| `Final_decay_ddPCR_datasheet.csv` | carboy experiment | `02_decay.R` |
| `ddpcr_combined_all.rds` | field experiment | `03_field_main.R`, `04_field_supplement.R`, `05_field_removal.R` |
| `ESP_timestamps_mLseawater.csv` | Brasseale et al. (2025) | `03`, `04`, `05` (C0 prior) |

## Reproducing the analysis

Open the project at its root (in RStudio, via the `.Rproj` file; otherwise
make sure `here::here()` resolves to the repository root). Then run the
scripts in order:

```r
source(here::here("code", "01_simulation.R"))         # Figure: simulation
source(here::here("code", "02_decay.R"))              # Figure: carboy LOO
                                                      # + outputs/carboy_rates.csv
source(here::here("code", "03_field_main.R"))         # Figure: main field
source(here::here("code", "04_field_supplement.R"))   # Figure: supplementary field
source(here::here("code", "05_field_removal.R"))      # Figure: removal sensitivity
source(here::here("code", "06_map.R"))                # Figure: study-area map (standalone)
```

`02_decay.R` writes `outputs/carboy_rates.csv`, which `00_setup.R` loads for
the subsequent field scripts. If you skip `02_decay.R`, the published values
are used as a fallback.

## Outputs

Figures are written to `plots/`:

* `simulation_main.png`, `simulation_supplement.png`
* `carboy_loo_grid.png`
* `field_grabs.png` (main text)
* `field_individual_samples.png` (supplementary)
* `field_removal.png`
* `map_study_area.tiff`

Per-fit tables are written to `outputs/`:

* `carboy_rates.csv` -- decay rates and t=0 offsets from the full carboy fit
* `age_estimates_grabs.csv`, `all_parameters_grabs.csv`
* `age_estimates_individual_samples.csv`, `all_parameters_individual_samples.csv`
* `age_estimates_grabs_removal.csv`

## Model overview

The joint age-estimation model (`conc_age_sim.stan` / `conc_age_field.stan`)
infers the elapsed time `t` for each sample by combining information across
markers with different decay rates `r` and known relative initial
concentrations `p`. It uses both the quantitative log-concentrations of
detected markers and the Bernoulli detection pattern of non-detections, so
non-detections carry information via a logistic detection function of the
latent mu:

```
mu[i,j] = C[i] + p[j] + r[j] * t[i]
y_obs | z=1   ~ Normal(mu, sigma)
z             ~ Bernoulli( logit^-1( alpha + beta * mu ) )
```

The two Stan files share this likelihood; they differ only in the priors
on `alpha`, `beta`, and `sigma` (tight for the simulation / carboy
calibration regime, wider for the field regime).

## Contact

[Corresponding author name + email redacted for peer review]
