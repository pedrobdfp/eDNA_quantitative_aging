# Quantitative aging of environmental DNA

Code and data accompanying the manuscript *[Quantitative Aging of Environmental DNA Using Multiple Fragment Sizes]*.

This repository reproduces the figures and parameter estimates in the manuscript
and provides the model in a form that can be applied to other systems.

---

## The method

Environmental DNA begins to degrade from the moment it is shed. Fragments of
different length degrade at different rates, so the relative abundance of
several markers from the same organism changes systematically with time since
shedding. That relative abundance is informative about age independently of the
absolute quantity released, which is rarely known.

The model has a single process equation. For unit *i*, a water sample or a set
of water samples collected together, and marker *j*:

```
log C[i,j]  =  C[i]  +  p[j]  +  r[j] * t[i]
```

| Term | Definition | Status |
|---|---|---|
| `t[i]` | time since shedding, hours | estimated |
| `C[i]` | log concentration at the time of shedding | estimated |
| `r[j]` | first-order decay rate of marker *j*, per hour, negative | fixed, from a decay experiment |
| `p[j]` | log concentration of marker *j* at *t* = 0 relative to the reference marker, so `p[1] = 0` | fixed, from the same experiment |

Age and initial concentration are separately identifiable because they act
differently on the marker set: `C[i]` shifts all markers by the same amount,
whereas `r[j] * t[i]` shifts each marker in proportion to its own decay rate.
Identification therefore depends on the decay rates differing between markers,
and the precision of an age estimate increases with the spread of those rates.

Observation error is partitioned into a biological component, between water
samples of the same unit, and a technical component, between PCR replicates of
the same water sample. Both are shared across markers, so that measurement error
widens the age posterior without displacing it.

---

## Two observation models

The process equation and the error structure are identical in both. They differ
in what format the data are supplied.

**`code/ddPCR/`** — droplet model, used in the manuscript. Takes the number of
positive droplets out of the number accepted, per well, and models them as
binomial with a complementary log-log link. Non-detections enter as exact zeros.
Applicable to droplet digital PCR only.

**`code/simplified/`** — concentration model, general. Takes an estimated
concentration and a detection indicator per replicate. Detection is modelled
with a logistic function whose 50% point is estimated, so non-detections
contribute information about the limit of detection rather than being discarded.
Applicable to quantitative PCR and digital PCR alike, and the appropriate
starting point for a new system.

Both are fitted here to the same wells, so results are directly comparable;
`code/compare_tracks.R` summarises them side by side.

---

## Quickstart

Open `eDNA_quantitative_aging.Rproj` and run [`quickstart.R`](quickstart.R). It
simulates a small dataset, fits the concentration model, and compares the
estimated ages with the values used to generate them.

The interface is four functions, defined in
[`code/simplified/age_functions.R`](code/simplified/age_functions.R):

```r
source(here::here("code", "simplified", "age_functions.R"))

toy  <- simulate_toy_edna(n_units = 12, n_bio = 2, n_reps = 3,
                          r = c(-0.05, -0.12, -0.25),
                          p = c(0, -0.5, -1.5))

fit  <- fit_edna_age(toy$obs, r = toy$r, p = toy$p,
                     C0_mean = 9, C0_sd = 2,
                     t_mean = 12, t_sd = 24)

ages <- age_table(fit, truth = toy$truth)
plot_ages(ages)
observation_table(fit)     # variance components and detection function
```

### Input format

One row per PCR replicate:

| `obs_i` | `obs_j` | `sample` | `detected` | `logC` |
|---|---|---|---|---|
| 1 | 1 | `unit1_s1` | 1 | 8.71 |
| 1 | 1 | `unit1_s1` | 1 | 8.55 |
| 1 | 2 | `unit1_s1` | 0 | `NA` |
| 1 | 1 | `unit1_s2` | 1 | 8.44 |

* `obs_i` — index of the unit being dated
* `obs_j` — index of the marker, ordered as in `r` and `p`
* `sample` — identifier of the water sample the replicate came from. Omit the
  column when each unit is a single water sample; the biological level is then
  switched off automatically.
* `detected` — 1 or 0
* `logC` — natural log of concentration, `NA` where the marker was not detected

Non-detections must be retained. A replicate that did not amplify constrains the
concentration to lie below the limit of detection, which is the expected result
for a marker that has been degrading for some time. Removing these rows biases
age estimates downward.

### Concentration units and the volumetric correction

The model does not impose a unit. `p[j]` is a difference of logarithms and
`r[j]` has units of inverse time, so neither carries a concentration scale;
whatever unit the concentrations are supplied in is the unit in which `C[i]` is
reported. The only requirement is that the prior on `C` (`C0_mean`, `C0_sd`) and
the prior on the detection threshold are expressed on the same scale as the
data.

What does matter is that concentrations refer to a **fixed volume of water
sampled**, not to a volume of PCR reaction or of extract. Converting from the
instrument reading to a concentration in water requires the dilution applied,
the volume of template added to the reaction, the elution volume, and above all
the **volume of water filtered**, which typically varies between samples.

Within a single sample this correction is common to all markers and is absorbed
by `C[i]`, so it does not bias that sample's age. It becomes important when
samples that filtered different volumes are compared, when several water samples
are pooled into one unit, and whenever the detection threshold is interpreted,
since that threshold is defined on the concentration scale.

In this repository the conversion is applied by `droplet_log_offset()` in each
setup file, which assembles the filtered volume, the dilution, the elution and
reaction volumes and the droplet volume into a single additive offset on the log
scale.

### Obtaining `r` and `p`

Both come from a decay experiment: water from the system of interest is held in
a sealed vessel under conditions resembling the field and subsampled over a
period long enough to observe a clear decline, with denser sampling early.
[`code/simplified/01_decay_simplified.R`](code/simplified/01_decay_simplified.R)
fits that experiment and writes `r` and `p` in the form the age model expects.

The values reported here apply to *Tursiops truncatus* eDNA in Hood Canal
seawater and should not be transferred to another system without validation.

---

## Repository layout

```
.
├── quickstart.R                     worked example on simulated data
├── code/
│   ├── simplified/                  concentration model
│   │   ├── 00_setup_simplified.R    packages, paths, constants, decay rates
│   │   ├── age_functions.R          four-function interface for new data
│   │   ├── 01_decay_simplified.R    decay rates and leave-one-out validation
│   │   ├── 02_simulation_simplified.R
│   │   ├── 03_field_main_simplified.R         field ages, pooled to grabs
│   │   ├── 04_field_supplement_simplified.R   field ages, per water sample
│   │   ├── conc_age_simplified.stan age model
│   │   └── decay_simplified.stan    decay model
│   ├── ddPCR/                       droplet model
│   │   ├── 00_setup_ddPCR.R
│   │   ├── 01_build_field_droplets.R          plate exports to well-level table
│   │   ├── 02_decay_ddPCR.R         decay rates and leave-one-out validation
│   │   ├── 03_simulation_ddPCR.R
│   │   ├── 04_field_main_ddPCR.R              field ages, pooled to grabs
│   │   ├── 05_field_supplement_ddPCR.R        field ages, per water sample
│   │   ├── 06_field_removal_ddPCR.R           field ages including physical removal
│   │   ├── droplet_age_ddPCR.stan   age model
│   │   └── decay_monophasic_ddPCR.stan   decay model
│   ├── map_study_area.R             study-area map
│   ├── ar1_decay_check.R            decay rates re-estimated from consecutive
│   │                                timepoints as an independent check
│   ├── compare_tracks.R             side-by-side summary of the two models
│   └── build_faire_package.R        writes data/faire/
├── data/
├── outputs/{simplified,ddPCR}/      tabular results
└── plots/{simplified,ddPCR}/        figures
```

---

## Data

| File | Contents |
|---|---|
| `data/field_droplets.csv` | Field dataset: 747 ddPCR wells, being 66 water samples × 4 markers × 3 replicates, with positive and accepted droplet counts, filtered volume and sample metadata. 146 wells returned no positive droplets. |
| `data/Final_decay_ddPCR_datasheet.csv` | Decay experiment: 234 wells across 3 carboys and 6 timepoints within the first 30 hours, from which `r` and `p` are estimated. |
| `data/ESP_timestamps_mLseawater.csv` | 42 reference concentrations from an autonomous sampler at the same site (Brasseale et al. 2025), used to set the prior on concentration at the time of shedding. |
| `data/ddpcr_combined_all.rds` | Sample-level field metadata, used when rebuilding the well-level table. |
| `data/faire/` | Both datasets in the FAIRe metadata standard, checklist v1.0.2. |

`data/faire/` is the archival form of the same data, written by
`code/build_faire_package.R`: project metadata, sample metadata, and one
amplification-data file per assay. See
[`data/faire/README_FAIRe.md`](data/faire/README_FAIRe.md) for the structure, the
terms added beyond the checklist to carry digital PCR partition counts, and the
laboratory metadata still to be supplied.

`code/ddPCR/01_build_field_droplets.R` regenerates `data/field_droplets.csv` from
the raw plate exports, applying the quality control: control wells excluded,
column-2 carryover extracts excluded, and, where a sample and assay were run on
more than one plate, only the plate with the higher mean concentration retained.
The raw exports are archived separately; the resulting table is included here, so
that script need not be run.

---

## Reproducing the analysis

R (≥ 4.3) and a working Stan toolchain via **rstan**
(<https://mc-stan.org/users/interfaces/rstan>).

```r
install.packages(c(
  "here", "rstan", "ggplot2", "dplyr", "tidyr", "tibble",
  "cowplot", "RColorBrewer", "forcats",
  # map only:
  "sf", "terra", "maptiles", "ggspatial"
))
```

Open the project at its root so that `here::here()` resolves there, then:

```r
# droplet model
source(here::here("code", "ddPCR", "02_decay_ddPCR.R"))
source(here::here("code", "ddPCR", "03_simulation_ddPCR.R"))
source(here::here("code", "ddPCR", "04_field_main_ddPCR.R"))
source(here::here("code", "ddPCR", "05_field_supplement_ddPCR.R"))
source(here::here("code", "ddPCR", "06_field_removal_ddPCR.R"))

# concentration model
source(here::here("code", "simplified", "01_decay_simplified.R"))
source(here::here("code", "simplified", "02_simulation_simplified.R"))
source(here::here("code", "simplified", "03_field_main_simplified.R"))
source(here::here("code", "simplified", "04_field_supplement_simplified.R"))

# map, comparison, and independent check on the decay rates
source(here::here("code", "map_study_area.R"))
source(here::here("code", "compare_tracks.R"))
source(here::here("code", "ar1_decay_check.R"))
```

The decay script of each track must be run first: it writes `carboy_rates.csv`,
which that track's setup file loads for the field analyses. Each track retains
its own estimates, since each is fitted under its own observation model. If the
decay script is not run, published values are used instead.

The full set takes several hours, most of it in the two simulation scripts.

---

## Citation and contact

Please cite the manuscript above. Corresponding author details are redacted
while the manuscript is under peer review.
