# Quantitative aging of environmental DNA

Code and data accompanying the manuscript *[Quantitative aging of environmental DNA using multiple components]*.

The repository reproduces every figure and parameter estimate in the manuscript,
and provides the model in a form you can apply to your own data.

---

## What the method does

eDNA in water starts decaying the moment it is shed. If you measure several
genetic markers from the same organism, and those markers decay at different
rates, then the **ratio between them** tells you how long the decay has been
going on. The ratio does not depend on how much DNA was released in the first
place, which is the quantity you almost never know.

A marker that decays quickly disappears first. So a sample dominated by the
slow-decaying marker is old; a sample where all markers are still in their
original proportions is fresh.

The model has one process equation. For unit *i* (a water sample, or a group of
samples collected together) and marker *j*:

```
log C[i,j]  =  C[i]  +  p[j]  +  r[j] * t[i]
```

| Symbol | Meaning | Where it comes from |
|---|---|---|
| `t[i]` | age of the sample, in hours | **estimated** — this is the answer |
| `C[i]` | log concentration at the moment of shedding | estimated |
| `r[j]` | decay rate of marker *j*, per hour, negative | known, from a decay experiment |
| `p[j]` | marker *j*'s log concentration at t = 0, relative to the first marker | known, from the same experiment |

Everything else in the repository is about connecting that equation to what a
PCR machine actually reports.

---

## Two versions of the model

The process equation above is identical in both. So is the observation noise:
error is split into a **biological** component (spread between water samples of
the same unit) and a **technical** component (spread between PCR wells of the
same water sample), both shared across markers.

They differ in one place only — **what the model reads from the instrument.**

|  | **ddPCR version** | **Simplified version** |
|---|---|---|
| Directory | `code/ddPCR/` | `code/simplified/` |
| Used for | the main text | the supplement, and for reuse |
| Input per replicate | positive droplets out of accepted droplets | a concentration and a detected / not-detected flag |
| Instrument | droplet digital PCR only | any qPCR or dPCR |
| How non-detections are handled | a well with zero positives is an exact zero of a binomial | a logistic detection function whose 50% point is estimated from the data |
| Extra parameters | none | two: the slope and intercept of that detection function |
| Age intervals | tighter | wider |

**Which one should you use?**

* You have ddPCR droplet counts and want the most out of them → **ddPCR
  version**. Counting droplets uses information a summarised concentration
  throws away: how many droplets amplified, out of how many, and the fact that
  a zero is exactly zero rather than "below some threshold".
* You have qPCR, or ddPCR concentrations without the raw counts, or you want a
  starting point → **simplified version**. It is the general form of the method
  and the one to build on.

The two are run side by side in this repository on the same wells, so that the
difference between them is attributable to the observation model and nothing
else. On the leave-one-carboy-out validation (18 held-out timepoints) the ddPCR
version gives intervals about a third narrower, but on this small test set it
is the less accurate and the less honest of the two:

| | ddPCR | simplified |
|---|---|---|
| bias (h) | −1.48 | 0.02 |
| RMSE (h) | 5.44 | 4.33 |
| 50% interval coverage | 0.28 | 0.50 |
| 95% interval coverage | 0.67 | 0.78 |
| median 95% interval width (h) | 12.7 | 18.4 |

Eighteen points is far too few to settle which observation model is better, and
the carboys are a more forgiving setting than the field — one vessel, no
biological replication, known times. Read the table as a caution against taking
the ddPCR version's narrower intervals at face value, not as evidence that the
simplified version is superior.

---

## Quickstart

Open `eDNA_quantitative_aging.Rproj`, then run
[`quickstart.R`](quickstart.R). It builds a small simulated dataset, fits the
simplified model to it, and checks the recovered ages against the true ones. It
takes about a minute, plus a one-off Stan compilation the first time.

The whole workflow is four functions, defined in
[`code/simplified/age_functions.R`](code/simplified/age_functions.R):

```r
source(here::here("code", "simplified", "age_functions.R"))

# 1. a toy dataset: 12 units, 2 water samples each, 3 wells per sample and
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
observation_table(fit)   # the two noise terms and the detection function
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

* `obs_i` — which unit you want an age for.
* `obs_j` — which marker, indexed in the same order as `r` and `p`.
* `sample` — which water sample the replicate came from. Leave the column out
  if you have one water sample per unit; the biological level of the model then
  switches off on its own.
* `detected` — 1 or 0.
* `logC` — natural log of copies per litre of water. `NA` where nothing
  amplified.

**Keep your non-detections.** A well that amplified nothing is not missing data.
It says the concentration was below the limit of detection, which is exactly
what a marker that has been decaying for a long time looks like. Both versions
of the model use that information; deleting those rows biases ages young.

### Three things to get right

1. **One volumetric scale, all the way through.** Everything here is log copies
   per **litre of water**. Divide by the volume you filtered before taking the
   log. If your instrument reports copies per mL, multiply by 1000 first.
2. **`p[1]` must be exactly 0.** Every offset is relative to the first marker.
3. **At least two markers with different decay rates**, and the further apart
   the rates, the sharper the age estimate. Panels (j)–(l) of the supplementary
   simulation figure show how much difference this makes.

### Getting `r` and `p` for your own markers

They come from a decay experiment: hold water from your system in a sealed
vessel, subsample it over time, and fit
[`code/simplified/01_decay_simplified.R`](code/simplified/01_decay_simplified.R),
which writes `r` and `p` in the format the age model expects. The values
published here are specific to *Tursiops truncatus* eDNA in Hood Canal seawater
and should not be transferred to another system without checking.

---

## Repository layout

```
.
├── quickstart.R                     five-minute tour on a toy dataset
├── code/
│   ├── ddPCR/                       MAIN TEXT: droplet observation model
│   │   ├── 00_setup_ddPCR.R         packages, paths, constants, decay rates
│   │   ├── 01_build_field_droplets.R  raw plate exports -> well-level table
│   │   ├── 02_decay_ddPCR.R         carboy decay fit + leave-one-out validation
│   │   ├── 03_simulation_ddPCR.R    simulation study
│   │   ├── 04_field_main_ddPCR.R    field ages, pooled to grabs
│   │   ├── 05_field_supplement_ddPCR.R  field ages, one per water sample
│   │   ├── 06_field_removal_ddPCR.R field ages with physical removal added
│   │   ├── droplet_age_ddPCR.stan   age model
│   │   └── decay_monophasic_ddPCR.stan  decay model
│   ├── simplified/                  SUPPLEMENT: concentration model
│   │   ├── 00_setup_simplified.R
│   │   ├── age_functions.R          the four-function interface
│   │   ├── 01_decay_simplified.R    carboy decay fit + leave-one-out validation
│   │   ├── 02_simulation_simplified.R
│   │   ├── 03_field_main_simplified.R    field ages, pooled to grabs
│   │   ├── 04_field_supplement_simplified.R  field ages, per water sample
│   │   ├── conc_age_simplified.stan age model
│   │   └── decay_simplified.stan    decay model
│   ├── map_study_area.R             study-area map (standalone)
│   ├── compare_tracks.R             side-by-side summary of the two models
│   └── build_faire_package.R        writes data/faire/
├── data/
│   ├── field_droplets.csv           the analysed field wells
│   ├── Final_decay_ddPCR_datasheet.csv   the carboy decay experiment
│   ├── ESP_timestamps_mLseawater.csv     reference concentrations for the prior
│   ├── ddpcr_combined_all.rds       field sample metadata
│   └── faire/                       both datasets in the FAIRe standard
├── outputs/{ddPCR,simplified}/      tables written by the scripts
└── plots/{ddPCR,simplified}/        figures written by the scripts
```

---

## Reproducing the manuscript

Open the project at its root so that `here::here()` resolves there. Then:

**Main text — ddPCR model**

```r
source(here::here("code", "ddPCR", "02_decay_ddPCR.R"))            # decay rates + LOO figure
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

**Comparing the two**

```r
source(here::here("code", "compare_tracks.R"))   # outputs/model_comparison.csv
```

Run the decay script of a track first: it writes `carboy_rates.csv`, which the
setup file of that track loads for the field analyses. Each track keeps its own
estimates, since each is fitted under its own observation model. If you skip it,
the published ddPCR values are used as a fallback.

`code/ddPCR/01_build_field_droplets.R` regenerates `data/field_droplets.csv`
from the raw plate exports. Those exports are archived separately; the table it
produces is shipped here, so you do not need to run it.

Every fit uses four chains of 2,000–4,000 iterations. The full set takes a few
hours on a laptop, most of it in the two simulation scripts.

---

## Dependencies

R ≥ 4.3 and a working Stan toolchain
(<https://mc-stan.org/users/interfaces/rstan>).

```r
install.packages(c(
  "here", "rstan", "ggplot2", "dplyr", "tidyr", "tibble",
  "cowplot", "RColorBrewer",
  # map only:
  "sf", "terra", "maptiles", "ggspatial"
))
```

---

## Data

| File | What it is |
|---|---|
| `data/field_droplets.csv` | 747 ddPCR wells: 66 water samples × 4 markers × 3 replicates, with positive and accepted droplet counts, filtered volume and sample metadata. This is the field dataset the manuscript analyses. |
| `data/Final_decay_ddPCR_datasheet.csv` | the carboy decay experiment, from which `r` and `p` are estimated |
| `data/ESP_timestamps_mLseawater.csv` | reference concentrations from an autonomous sampler at the same site (Brasseale et al. 2025), used only to set the prior on concentration at t = 0 |
| `data/ddpcr_combined_all.rds` | sample-level field metadata, used by the droplet builder |

`data/faire/` holds the same two datasets — decay experiment and field survey —
in the FAIRe (FAIR environmental DNA) metadata format, checklist v1.0.2, written
by `code/build_faire_package.R`. See
[`data/faire/README_FAIRe.md`](data/faire/README_FAIRe.md) for the file
structure, the terms added beyond the checklist for digital PCR partition
counts, and the laboratory metadata that still has to be filled in before
archiving.

---

## Model files in detail

### `code/ddPCR/droplet_age_ddPCR.stan`

```
log C[i,j,s]   = C[i] + p[j] + r[j]*t[i] + eta[s,j]      eta ~ Normal(0, sigma_bio)
omega[i,j,r,s] = log C[i,j,s] + eps + log_offset - correction
                                                          eps ~ Normal(0, sigma_tech)
W              ~ Binomial(U, 1 - exp(-exp(omega)))
```

`W` positive droplets out of `U` accepted; `log_offset` converts a seawater
concentration into expected copies per droplet given the volume filtered, the
dilution and the droplet volume. `correction` is `(sigma_tech² + sigma_bio²)/2`,
which the exponential link requires; because both sigmas are shared across
markers it is a single constant that `C[i]` absorbs, so observation error widens
the age posterior without displacing it.

### `code/simplified/conc_age_simplified.stan`

```
log C[i,j,s] = C[i] + p[j] + r[j]*t[i] + eta[s,j]        eta ~ Normal(0, sigma_bio)
z            ~ Bernoulli(logit^-1(beta * (log C[i,j,s] - logC50)))
y            ~ Normal(log C[i,j,s], sigma_tech)          where z = 1
```

`z` is detected / not detected for every replicate, `y` the measured log
concentration where detected. No mean correction is needed here: the likelihood
is written directly on the log scale and is symmetric.

The detection function is parameterised by its own 50% point, `logC50` — the
concentration at which half of replicates amplify, which is the limit of
detection the model infers from your non-detections — rather than by an
intercept. The equivalent intercept `alpha = -beta * logC50` is reported
alongside it. This is the same likelihood written two ways, but `logC50` and
`beta` are close to independent whereas an intercept and a slope are strongly
correlated, and the sampler is several times faster for it.

Both files switch the biological level off when no biological replication is
supplied, which is the right choice for a decay experiment where one vessel is
sampled repeatedly.

---

## Contact

[Corresponding author name and email redacted for peer review]
