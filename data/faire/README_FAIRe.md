# FAIRe-formatted data

The data analysed in *Quantitative aging of environmental DNA using multiple
components* are archived here in the FAIRe (FAIR environmental DNA) metadata
format, checklist version 1.0.2.

* Standard: <https://fair-edna.github.io>
* Reference: Takahashi et al. (2025) *Environmental DNA*,
  [doi:10.1002/edn3.70100](https://doi.org/10.1002/edn3.70100)
* Written by [`code/build_faire_package.R`](../../code/build_faire_package.R)

Project identifier: `ednaAging`

## Files

| File | Rows | Contents |
|---|---|---|
| `projectMetadata_ednaAging.csv` | 40 terms | project- and assay-level terms, one row per term, one column per assay |
| `sampleMetadata_ednaAging.csv` | 84 samples | one row per water sample: 66 field, 18 decay-experiment |
| `ampData_ednaAging_Tt_Cytb.csv` | 234 | one row per PCR replicate |
| `ampData_ednaAging_Tt_16S.csv` | 273 | " |
| `ampData_ednaAging_Tt_DLL1.csv` | 237 | " |
| `ampData_ednaAging_Tt_longFrag.csv` | 237 | " |

`samp_name` links `ampData` to `sampleMetadata`; `assay_name` links both to the
assay columns of `projectMetadata`.

## The two datasets

A `dataset` column in `sampleMetadata` and `ampData` separates them.

**Decay experiment** — three carboys of seawater, subsampled repeatedly, used to
estimate each marker's decay rate. Archived here are the records the decay model
consumes: the DNA fraction of carboys 1–3 within the first 30 hours, controls
excluded. Sample names are `carboy<n>_<hours>h`.

**Field survey** — two stations in Hood Canal, Washington, sampled over one day
at several depths, used to estimate the age of the eDNA in each water sample.
Archived here are the DNA1 extracts that pass the plate-level quality control in
[`code/ddPCR/01_build_field_droplets.R`](../../code/ddPCR/01_build_field_droplets.R):
control wells removed, column-2 carryover extracts removed, and, where a sample
and assay were run on more than one plate, only the plate with the higher mean
concentration retained. Every sample and marker then has exactly three
replicates.

Samples taken at the same station, depth and instrument within ten minutes of
one another are biological replicates of one grab; `biological_rep_relation`
lists the members of each grab.

## Four assays, three reactions

`Tt_longFrag` is not a separate PCR. It is the double-positive channel of the
`Tt_16S` / `Tt_DLL1` duplex reaction: droplets that carry both probes, which
indicates a template molecule long enough to span both target regions. Its
replicates therefore share `pcr_plate_id` and `well_id` with the `Tt_16S` and
`Tt_DLL1` rows of the same well, and `assay_reaction_format` in
`projectMetadata` records this.

The positive-droplet counts for `Tt_longFrag` are the raw double positives and
are **not** corrected for double positives expected by chance. For these
concentrations the correction is negligible — single-marker counts are single
digits, so chance doubles are of order 10⁻³ droplets — but the columns needed to
apply it (`positivePartitions` of the two single markers, `acceptedPartitions`)
are all present should a reuser want to.

## Terms added beyond the checklist

FAIRe v1.0.2 has no slots for digital PCR partition counts, which are the raw
observation this study models. Following the guideline that new terms may be
added with clear, descriptive names, `ampData` carries these additional columns:

| Column | Meaning |
|---|---|
| `detected_notDetected` | `detected` if the replicate had at least one positive partition |
| `positivePartitions` | positive droplets in the well |
| `acceptedPartitions` | total droplets accepted by the reader |
| `partitionVolumeMicroliters` | effective droplet volume, 0.00085 µL |
| `concentrationInReaction` | concentration in the reaction, copies µL⁻¹ |
| `concentrationInReaction_unit` | `copies/uL` |
| `volumeFilteredLiters` | volume of water filtered for that sample |
| `dilutionFactor` | dilution applied to the extract before the reaction |
| `well_id` | plate well, so duplex channels of one reaction can be identified |
| `dataset` | `field survey` or `decay experiment` |

`quantificationCycle` and `rfu` are qPCR terms and are set to
`not applicable: dPCR`.

`estimatedNumberOfCopies` is the seawater concentration in copies per litre —
the reaction concentration converted using the volume filtered — because that is
the quantity the models are written on. `concentrationInReaction` keeps the
instrument's own output alongside it.

## Source of the laboratory metadata

Assay chemistry, thermocycling, extraction, decontamination and controls are as
published in the decay study these datasets were generated alongside
(*Differential decay of environmental nucleic acid components from a marine
mammal*, Scientific Reports; Table 1 and Table S1), and are transcribed into
`projectMetadata` by `build_faire_package.R`.

Both datasets share the assays and the extraction protocol. They do **not**
share filtration: the decay samples were filtered in the laboratory through a
tandem stack of 5.0, 1.0 and 0.45 µm MCE filters, while the field samples were
filtered in situ by three different autonomous samplers (Smith-Root, Pufferfish,
Ascension G). Filter terms are therefore recorded per dataset, and are left
missing for the field samples.

Both experiments share a sampling date, **2024-10-16**: the water for the decay
carboys was collected during the same field day and the incubation began at the
laboratory that evening. Elapsed time for the decay samples is measured from
collection and is carried by `Hours_base` in the source datasheet.

## Still to be completed before archiving

Seven terms remain at the INSDC value `missing: not provided`. None are
derivable from the analysed data or from the decay study.

**Project level** — `project_contact`, `institution`, `license`, `recordedBy`,
`bibliographicCitation`. Administrative; fill in at submission.

**Assay level** — `pcr_assay_lod`, `pcr_assay_loq`. The decay study does not
report a formal limit of detection or quantification per assay. Note that the
simplified model in `code/simplified/` estimates an effective 50% detection
point from the non-detections (about 12.5 copies L⁻¹ for the field samples),
which is not the same quantity as a validated assay LOD and should not be
entered here as one.

**Sample level** — `filter_name` and `filter_material` for the field samples
only (see above); `samp_store_dur` for the field samples; `date_ext` for both;
`tidal_stage` for the field samples, which is available from the NOAA tide
record for station 9445478 (Union, Hood Canal) if wanted.

## Validating

Once the terms above are filled, the files can be checked against the checklist
with the FAIRe-fier web validator, <https://shiny.csiro.au/FAIRe-fier/>.
