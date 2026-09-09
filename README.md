# Quantitative aging of environmental DNA

Code and data accompanying the manuscript *[Quantitative Aging of Environmental DNA Using Multiple Fragment Sizes]*.

This repository reproduces every figure and parameter estimate in the manuscript,
and provides the model in a form you can apply to your own data.

---

## What the method does

eDNA starts breaking down the moment it is shed. If you measure several genetic
markers from the same organism, and those markers break down at different
speeds, then the **balance between them** tells you how long the breakdown has
been going on.

A sample in which the fragile markers have already faded is old. A sample in
which all markers remain in their original proportions is fresh. Crucially, that
balance does not depend on how much DNA was released to begin with — which is a
quantity you almost never know.

The whole method is one equation. For unit *i* (a water sample, or a group of
samples collected together) and marker *j*:

```
log C[i,j]  =  C[i]  +  p[j]  +  r[j] * t[i]
```

| Symbol | Meaning | Where it comes from |
|---|---|---|
| `t[i]` | age of the sample, in hours | **estimated — this is the answer** |
| `C[i]` | log concentration at the moment of shedding | estimated |
| `r[j]` | decay rate of marker *j*, per hour, negative | known, from a decay experiment |
| `p[j]` | marker *j*'s starting level relative to the first marker | known, from the same experiment |

Age and concentration do not get confused with one another because they leave
different fingerprints. `C[i]` raises or lowers every marker by the same amount;
`t[i]` is multiplied by `r[j]`, which differs between markers, so age *tilts*
the pattern rather than shifting it. If all markers decayed at the same rate the
two would be indistinguishable and the method would not work.

Everything else in this repository connects that equation to what a PCR machine
actually reports.

---

## Two versions of the model

The equation above is identical in both. So is the treatment of measurement
error, which is split into a **biological** component (variation between water
samples of the same unit) and a **technical** component (variation between PCR
replicates of the same water sample).

They differ in one place only: **what the model reads from the instrument.**

|  | **Simplified version** | **ddPCR version** |
|---|---|---|
| Directory | `code/simplified/` | `code/ddPCR/` |
| Used for | the supplement, and for reuse | the main text |
| Input per replicate | a concentration and a detected / not-detected flag | positive droplets out of accepted droplets |
| Instrument | any qPCR or digital PCR | droplet digital PCR only |
| How non-detections are used | a detection curve whose 50% point is estimated from the data | a zero is an exact zero of a binomial |
| Extra parameters | two, describing the detection curve | none |

**Which should you use?**

* **Starting out, or working with qPCR** → the **simplified version**. It is the
  general form of the method and the one to build on. Begin with
  [`quickstart.R`](quickstart.R).
* **You have ddPCR droplet counts** → the **ddPCR version** uses information a
  summarised concentration throws away: how many droplets amplified out of how
  many, and the fact that a zero is exactly zero rather than "below some
  threshold".

Both are run here on the same wells, so any difference between them is
attributable to the observation model and nothing else.

### How the two compare on these data

Decay rates agree closely except for the Bridge marker, which has the lowest
concentrations and the most non-detections:

| Marker | ddPCR | simplified |
|---|---|---|
| Cytb | −0.119 | −0.117 |
| 16S | −0.161 | −0.161 |
| D-loop | −0.178 | −0.172 |
| Bridge | −0.189 | −0.171 |

Both recover the same variation between water samples (σ_bio 0.76 and 0.77). They
differ in the technical term (σ_tech 0.089 and 0.345) for an understandable
reason: the droplet model represents counting noise explicitly through the
binomial, so its σ_tech measures only what is left over, whereas a summarised
concentration already contains counting noise and its σ_tech must carry both.

On the leave-one-carboy-out validation (18 held-out timepoints):

| | ddPCR | simplified |
|---|---|---|
| bias (h) | 0.30 | 1.59 |
| mean absolute error (h) | 4.31 | 3.80 |
| RMSE (h) | 5.42 | 4.61 |
| 50% interval coverage | 0.11 | 0.44 |
| 95% interval coverage | 0.61 | 0.83 |
| median 95% interval width (h) | 12.0 | 17.9 |

The droplet model returns intervals about a third narrower. On this small
validation set it is nonetheless the less well calibrated of the two, with half
intervals covering the truth in 2 of 18 cases against a nominal 9. Eighteen
points cannot settle which observation model is preferable in general, and the
carboys are a more forgiving setting than the field — a single vessel, no
biological replication, known elapsed times. The table is reported as a caution
against reading narrower intervals as straightforward gains in precision.

On the field data the two order the grabs almost identically (correlation 0.97
between per-grab medians) and agree that eDNA at the further station is
substantially older, though they disagree about how old the oldest grabs are.

---

## Quickstart

Open `eDNA_quantitative_aging.Rproj`, then run [`quickstart.R`](quickstart.R).
It builds a small simulated dataset, fits the simplified model to it, and checks
the recovered ages against the ages that went in. About a minute, plus a one-off
Stan compilation the first time.

The whole workflow is four functions, defined in
[`code/simplified/age_functions.R`](code/simplified/age_functions.R):

```r
source(here::here("code", "simplified", "age_functions.R"))

# 1. a toy dataset: 12 units, 2 water samples each, 3 replicates per sample and
#    marker, 3 markers decaying at -0.05, -0.12 and -0.25 per hour
toy <- simulate_toy_edna(n_units = 12, n_bio = 2, n_reps = 3,
                         r = c(-0.05, -0.12, -0.25),
                         p = c(0, -0.5, -1.5))

# 2. fit
fit <- fit_edna_age(toy$obs, r = toy$r, p = toy$p,
                    C0_mean = 9, C0_sd = 2,
                    t_mean = 12, t_sd = 24)

# 3. read the ages
ages <- age_table(fit, truth = toy$truth)

# 4. look at them
plot_ages(ages)
observation_table(fit)   # the two noise terms and the detection curve
```

### The data format

`fit_edna_age()` wants one row per PCR replicate:

| `obs_i` | `obs_j` | `sample` | `detected` | `logC` |
|---|---|---|---|---|
| 1 | 1 | `unit1_s1` | 1 | 8.71 |
| 1 | 1 | `unit1_s1` | 1 | 8.55 |
| 1 | 2 | `unit1_s1` | 0 | `NA` |
| 1 | 1 | `unit1_s2` | 1 | 8.44 |
| … | | | | |

* `obs_i` — which unit you want an age for
* `obs_j` — which marker, indexed in the same order as `r` and `p`
* `sample` — which water sample the replicate came from. Leave the column out
  if you have one water sample per unit; the biological level then switches off
  on its own.
* `detected` — 1 or 0
* `logC` — natural log of copies per litre of water, `NA` where nothing amplified

**Keep your non-detections.** A replicate that amplified nothing is not missing
data. It says the concentration was below what your assay can see, which is
exactly what a marker that has been decaying a long time looks like. Both
versions of the model use that information; deleting those rows biases ages
young.

### Three things to get right

1. **One volumetric scale throughout.** Everything here is log copies per
   **litre of water**. Divide by the volume you filtered before taking the log.
   If your instrument reports copies per mL, multiply by 1000 first.
2. **`p[1]` must be exactly 0.** Every offset is relative to the first marker.
3. **At least two markers with different decay rates**, and the further apart
   the rates, the sharper the age estimate. The supplementary simulation figure
   shows how much difference this makes.

### Getting `r` and `p` for your own markers

They come from a decay experiment: hold water from your system in a sealed
vessel, keep it under conditions resembling the field, and subsample it over a
period long enough for a clear decline. Denser sampling early is worth more than
a few distant points. Then fit
[`code/simplified/01_decay_simplified.R`](code/simplified/01_decay_simplified.R),
which writes `r` and `p` in the format the age model expects.

The values published here are specific to *Tursiops truncatus* eDNA in Hood Canal
seawater and should not be transferred to another system without checking.

---

## Repository layout

```
.
├── quickstart.R                     five-minute tour on a toy dataset
├── code/
│   ├── simplified/                  concentration model: qPCR or digital PCR
│   │   ├── 00_setup_simplified.R    packages, paths, constants, decay rates
│   │   ├── age_functions.R          the four-function interface
│   │   ├── 01_decay_simplified.R    decay rates + leave-one-out validation
│   │   ├── 02_simulation_simplified.R
│   │   ├── 03_field_main_simplified.R    field ages, pooled to grabs
│   │   ├── 04_field_supplement_simplified.R   field ages, per water sample
│   │   ├── conc_age_simplified.stan age model
│   │   └── decay_simplified.stan    decay model
│   ├── ddPCR/                       droplet model: digital PCR counts
│   │   ├── 00_setup_ddPCR.R
│   │   ├── 01_build_field_droplets.R   raw plate exports -> well-level table
│   │   ├── 02_decay_ddPCR.R         decay rates + leave-one-out validation
│   │   ├── 03_simulation_ddPCR.R
│   │   ├── 04_field_main_ddPCR.R    field ages, pooled to grabs
│   │   ├── 05_field_supplement_ddPCR.R   field ages, per water sample
│   │   ├── 06_field_removal_ddPCR.R field ages with physical removal added
│   │   ├── droplet_age_ddPCR.stan   age model
│   │   └── decay_monophasic_ddPCR.stan   decay model
│   ├── map_study_area.R             study-area map
│   ├── ar1_decay_check.R            decay rates re-estimated from consecutive
│   │                                timepoints, as an independent check
│   ├── compare_tracks.R             side-by-side summary of the two models
│   └── build_faire_package.R        writes data/faire/
├── data/                            see below
├── outputs/{simplified,ddPCR}/      tables written by the scripts
└── plots/{simplified,ddPCR}/        figures written by the scripts
```

---

## Data

| File | Contents |
|---|---|
| `data/field_droplets.csv` | The analysed field dataset: 747 ddPCR wells, being 66 water samples × 4 markers × 3 replicates, with positive and accepted droplet counts, filtered volume and sample metadata. 146 wells produced no positive droplets. |
| `data/Final_decay_ddPCR_datasheet.csv` | The carboy decay experiment, from which `r` and `p` are estimated. 234 wells across 3 carboys and 6 timepoints within the first 30 hours. |
| `data/ESP_timestamps_mLseawater.csv` | 42 reference concentrations from an autonomous sampler at the same site (Brasseale et al. 2025), used only to set the prior on concentration at the moment of shedding. |
| `data/ddpcr_combined_all.rds` | Sample-level field metadata, used by the droplet builder. |
| `data/faire/` | Both datasets in the FAIRe (FAIR environmental DNA) metadata standard, checklist v1.0.2. |

`data/faire/` is the archival form of the same data, written by
`code/build_faire_package.R`. It contains project metadata, sample metadata and
one amplification-data file per assay. See
[`data/faire/README_FAIRe.md`](data/faire/README_FAIRe.md) for the file
structure, the terms added beyond the checklist to carry digital PCR partition
counts, and the laboratory metadata still to be supplied.

`code/ddPCR/01_build_field_droplets.R` regenerates `data/field_droplets.csv`
from the raw plate exports, including the quality control: control wells
removed, column-2 carryover extracts removed, and, where a sample and assay were
run on more than one plate, only the plate with the higher mean concentration
retained. The raw exports are archived separately; the table it produces is
shipped here, so that script does not need to be run.

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

**Main text — ddPCR model**

```r
source(here::here("code", "ddPCR", "02_decay_ddPCR.R"))            # decay rates + validation
source(here::here("code", "ddPCR", "03_simulation_ddPCR.R"))       # simulation figures
source(here::here("code", "ddPCR", "04_field_main_ddPCR.R"))       # main field figure
source(here::here("code", "ddPCR", "05_field_supplement_ddPCR.R")) # per-sample ages
source(here::here("code", "ddPCR", "06_field_removal_ddPCR.R"))    # removal sensitivity
source(here::here("code", "map_study_area.R"))                     # study-area map
```

**Supplement — simplified model**

```r
source(here::here("code", "simplified", "01_decay_simplified.R"))
source(here::here("code", "simplified", "02_simulation_simplified.R"))
source(here::here("code", "simplified", "03_field_main_simplified.R"))
source(here::here("code", "simplified", "04_field_supplement_simplified.R"))
```

**Comparisons and checks**

```r
source(here::here("code", "compare_tracks.R"))     # outputs/model_comparison.csv
source(here::here("code", "ar1_decay_check.R"))    # independent check on the decay rates
```

Run the decay script of a track first: it writes `carboy_rates.csv`, which that
track's setup file loads for the field analyses. Each track keeps its own
estimates, since each is fitted under its own observation model. If the decay
script is skipped, published values are used instead.

The full set takes a few hours, most of it in the two simulation scripts.

---

## The models in brief

### `code/simplified/conc_age_simplified.stan`

```
log C[i,j,s] = C[i] + p[j] + r[j]*t[i] + eta[s,j]     eta ~ Normal(0, sigma_bio)
z            ~ Bernoulli(logit^-1(beta * (log C[i,j,s] - logC50)))
y            ~ Normal(log C[i,j,s], sigma_tech)       where z = 1
```

`z` is detected / not detected for every replicate, `y` the measured log
concentration where detected. The detection curve is described by `logC50`, the
concentration at which half of replicates amplify — the effective limit of
detection — and `beta`, how sharply detection switches on there. The equivalent
intercept of the more familiar `alpha + beta * log C` form is
`alpha = -beta * logC50` and is reported alongside.

### `code/ddPCR/droplet_age_ddPCR.stan`

```
log C[i,j,s]   = C[i] + p[j] + r[j]*t[i] + eta[s,j]   eta ~ Normal(0, sigma_bio)
omega          = log C[i,j,s] + eps + log_offset - correction
                                                      eps ~ Normal(0, sigma_tech)
W              ~ Binomial(U, 1 - exp(-exp(omega)))
```

`W` positive droplets out of `U` accepted. `log_offset` converts a seawater
concentration into expected copies per droplet, given the volume filtered, the
dilution and the droplet volume. `correction` is `(sigma_tech² + sigma_bio²)/2`,
which the exponential link requires; because both standard deviations are shared
across markers it is a single constant absorbed by `C[i]`, so measurement error
widens the age estimate without moving it.

Both files switch the biological level off when no biological replication is
supplied, which is the right choice for a decay experiment where one vessel is
sampled repeatedly.

---

## Citation and contact

Please cite the manuscript above. Corresponding author details are redacted
while the manuscript is under peer review.
