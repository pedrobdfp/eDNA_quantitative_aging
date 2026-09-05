# =============================================================================
# build_faire_package.R
# -----------------------------------------------------------------------------
# Writes the data underlying the manuscript into data/faire/ in the FAIRe
# (FAIR environmental DNA) metadata format, version 1.0.2.
#
#   https://fair-edna.github.io
#   Takahashi et al. (2025) Environmental DNA, doi:10.1002/edn3.70100
#
# Two datasets are covered, and only the records the manuscript actually
# analyses:
#
#   decay experiment   three carboys sampled repeatedly over the first 30 hours,
#                      DNA fraction, controls excluded
#   field survey       two stations in Hood Canal, one day, four markers,
#                      DNA1 extracts, after the plate-level quality control in
#                      code/ddPCR/01_build_field_droplets.R
#
# FILES WRITTEN, following the FAIRe convention [datatype]_[project_id]_[assay]:
#
#   projectMetadata_ednaAging.csv          project- and assay-level terms
#   sampleMetadata_ednaAging.csv           one row per water sample
#   ampData_ednaAging_<assay>.csv          one row per PCR replicate, per assay
#   README_FAIRe.md                        structure, deviations, and the terms
#                                          that still need to be filled in
#
# TERMS ADDED BEYOND THE CHECKLIST
#   FAIRe v1.0.2 has no slots for droplet digital PCR partition counts, which
#   are the raw observation this manuscript models. They are added to ampData
#   with descriptive names and documented in README_FAIRe.md:
#     positivePartitions, acceptedPartitions, partitionVolumeMicroliters,
#     concentrationInReaction, concentrationInReaction_unit,
#     volumeFilteredLiters, detected_notDetected
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
})

out_dir <- here("data", "faire")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

PROJECT_ID <- "ednaAging"

# -----------------------------------------------------------------------------
# Sampling dates
# -----------------------------------------------------------------------------
# Neither source table carries a calendar date, only clock times, so both are
# supplied here. The two experiments share a date: the water for the decay
# carboys was collected during the same field day and transported to the
# laboratory, where the incubation began that evening. Elapsed time for the
# decay samples is given by Hours_base, measured from the moment of collection.

FIELD_EVENT_DATE <- "2024-10-16"
DECAY_EVENT_DATE <- "2024-10-16"

MISSING <- "missing: not provided"

# -----------------------------------------------------------------------------
# Laboratory metadata
# -----------------------------------------------------------------------------
# Assay chemistry, thermocycling, extraction and filtration are as published in
# the decay study that these two datasets were generated alongside:
#
#   Differential decay of environmental nucleic acid components from a marine
#   mammal (Scientific Reports), Table 1 and Table S1.
#
# The assays and the extraction protocol are shared by both datasets. Filtration
# is not: the decay samples were filtered through a fixed tandem MCE stack in
# the laboratory, whereas the field samples were filtered in situ by three
# different autonomous samplers, so filter terms are recorded per dataset.

DECAY_REF <- paste("Differential decay of environmental nucleic acid components",
                   "from a marine mammal, Scientific Reports")

# Oligonucleotides (Table S1). Bridge has none of its own: it is the
# double-positive channel of the 16S / D-loop duplex.
OLIGO <- list(
  Tt_Cytb = list(
    fwd_name = "Ttru-cytb_F", fwd = "TTATTCTTCCATTCATCATCAC",
    rev_name = "Ttru-cytb_R", rev = "GTGGGGTTGTTGGATCCTGT",
    probe    = "/5FAM/AATAGTAGG/ZEN/TGAACGGCTGCCA/3IABkFQ/",
    reporter = "FAM", primer_conc = "900 nM each",
    anneal   = "60", cycles = "45", amplicon = "79",
    ref      = "primers from Xiong et al. 2025; probe designed in the decay study"),
  Tt_16S = list(
    fwd_name = "Ttru-16S_F", fwd = "AGGGTTTTACTGTCTCTTACTCTT",
    rev_name = "Ttru-16S_R", rev = "CTTTGTTATCCCTTGGTGGTATTA",
    probe    = "/5-HEX/ATAATGCAA/ZEN/TAAGACGAGAAGACCCTATGG/3IABkFQ/",
    reporter = "HEX", primer_conc = "600 nM each",
    anneal   = "56", cycles = "45", amplicon = "146",
    ref      = "designed in the decay study"),
  Tt_DLL1 = list(
    fwd_name = "DL1-HL_F", fwd = "CACCCAAAGCTGRARTTCTDYATAAACT",
    rev_name = "Oordlp4",  rev = "GCGGGTTGCTGGTTTCACG",
    probe    = "/5FAM/ACTGACGTA/ZEN/GTACTGTGATGTTGTGGTAACTGTAC/3IABkFQ/",
    reporter = "FAM", primer_conc = "900 nM each",
    anneal   = "56", cycles = "45", amplicon = "390",
    ref      = paste("forward designed in the decay study;",
                     "reverse Oordlp4 from Baker et al. 2018")),
  Tt_longFrag = list(
    fwd_name = "not applicable: derived channel",
    fwd      = "not applicable: derived channel",
    rev_name = "not applicable: derived channel",
    rev      = "not applicable: derived channel",
    probe    = "not applicable: derived channel",
    reporter = "FAM and HEX", primer_conc = "not applicable: derived channel",
    anneal   = "56", cycles = "45", amplicon = ">=2746",
    ref      = "double-positive channel of the 16S / D-loop duplex")
)

og <- function(a, f) OLIGO[[a]][[f]]

NUCL_ACID_EXT <- paste(
  "Filters thawed in DNA/RNA Shield, incubated 15 min at 37 C with agitation,",
  "buffer concentrated from 2 mL to ~400 uL on an Amicon Ultra-15 30 kDa unit",
  "(Millipore-Sigma), then coextracted with the Quick-DNA/RNA Miniprep Kit",
  "(Zymo Research) with a 30 min proteinase K incubation before binding buffer.",
  "The protocol yields two sequential extracts; only the first (DNA1) is",
  "analysed here.")

STERILISE <- paste("Forceps, filter holders, adapters and tubing soaked in 5%",
                   "bleach, rinsed with deionised water and dried before reuse;",
                   "carboys and collection vessels cleaned overnight in 1%",
                   "bleach and rinsed with deionised water.")

# Station coordinates, WGS84 decimal degrees (as used in code/map_study_area.R)
station_coords <- tibble::tribble(
  ~station,          ~lat,       ~lon,
  "Husbandry Area",  47.742236,  -122.729975,
  "NOAA Boat",       47.736441,  -122.743206
)

# Marker to assay name. Tt_longFrag is not a separate PCR: it is the
# double-positive channel of the 16S / D-loop duplex reaction, so its
# replicates share a plate and well with Tt_16S and Tt_DLL1.
assay_of_marker <- c(cytb        = "Tt_Cytb",
                     Tt_16S      = "Tt_16S",
                     Tt_DLL1     = "Tt_DLL1",
                     Tt_longFrag = "Tt_longFrag")

assay_of_decay_marker <- c(Cytb   = "Tt_Cytb",
                           `16S`  = "Tt_16S",
                           Dloop  = "Tt_DLL1",
                           Bridge = "Tt_longFrag")

# =============================================================================
# 1. Field survey
# =============================================================================

field <- read.csv(here("data", "field_droplets.csv"), stringsAsFactors = FALSE)

field_samples <- field %>%
  group_by(sample) %>%
  summarise(id                = first(id),
            bio_rep           = first(bio_rep),
            station           = first(station),
            Instrument        = first(Instrument),
            depth             = first(depth),
            time_collected    = first(time_collected),
            filter_start_time = first(filter_start_time),
            vol               = first(vol),
            .groups = "drop") %>%
  left_join(station_coords, by = "station") %>%
  mutate(time_str = sprintf("%02d:%02d",
                            time_collected %/% 100L, time_collected %% 100L),
         grab_key = paste(station, Instrument, depth, sep = "|"))

# Water samples taken at the same station, depth and instrument within ten
# minutes of one another are biological replicates of the same grab. This is
# the grouping the grab-level analysis uses.
field_samples <- field_samples %>%
  arrange(station, Instrument, depth, time_collected) %>%
  group_by(station, Instrument, depth) %>%
  mutate(time_min = (time_collected %/% 100L) * 60L + time_collected %% 100L,
         grab_n   = cumsum(c(TRUE, diff(time_min) > 10)),
         grab_id  = paste(station, Instrument, depth, grab_n, sep = "_")) %>%
  ungroup()

grab_members <- field_samples %>%
  group_by(grab_id) %>%
  summarise(members = paste(sort(sample), collapse = " | "), .groups = "drop")

field_samples <- field_samples %>% left_join(grab_members, by = "grab_id")

field_sample_meta <- field_samples %>%
  transmute(
    samp_name                     = sample,
    materialSampleID              = sample,
    samp_category                 = "sample",
    eventDate                     = FIELD_EVENT_DATE,
    verbatimEventTime             = time_str,
    decimalLatitude               = lat,
    decimalLongitude              = lon,
    verbatimCoordinateSystem      = "decimal degrees",
    verbatimSRS                   = "WGS84",
    geo_loc_name                  = "USA: Washington, Hood Canal",
    env_broad_scale               = "marine biome [ENVO:00000447]",
    env_local_scale               = "estuarine fjord [ENVO:00000406]",
    env_medium                    = "sea water [ENVO:00002149]",
    habitat_natural_artificial_0_1 = 0,
    minimumDepthInMeters          = depth,
    maximumDepthInMeters          = depth,
    samp_collect_device           = Instrument,
    samp_collect_method           = "in situ filtration of seawater",
    samp_size                     = vol,
    samp_size_unit                = "L",
    biological_rep_relation       = members,
    sample_derived_from           = "not applicable: collected in situ",
    # Filtration differs by sampler and is not recorded in the decay study.
    filter_name                   = MISSING,
    filter_material               = MISSING,
    filter_diameter               = NA_real_,
    filter_passive_active_0_1     = 1,
    samp_store_sol                = "DNA/RNA Shield (Zymo)",
    samp_store_temp               = -80,
    samp_store_dur                = MISSING,
    nucl_acid_ext                 = NUCL_ACID_EXT,
    nucl_acid_ext_kit             = "Quick-DNA/RNA Miniprep Kit (Zymo Research)",
    date_ext                      = MISSING,
    samp_vol_we_dna_ext           = 2,
    samp_vol_we_dna_ext_unit      = "mL of preservation buffer",
    tidal_stage                   = MISSING,
    temp                          = NA_real_,
    salinity                      = NA_real_,
    dataset                       = "field survey"
  )

field_amp <- field %>%
  mutate(assay_name        = unname(assay_of_marker[marker]),
         samp_name         = sample,
         technical_rep_id  = rep,
         pcr_plate_id      = plate,
         well_id           = well,
         detected_notDetected = ifelse(positives > 0, "detected", "notDetected"),
         # seawater concentration, from the reaction concentration and the
         # volume filtered; the conversion used throughout the analysis
         estimatedNumberOfCopies = ifelse(
           positives > 0,
           exp(log(conc_machine) - log_offset + log(0.00085)),
           0),
         estimatedNumberOfCopies_unit   = "copies/L",
         estimatedNumberOfCopies_method = paste(
           "ddPCR Poisson estimate from partition counts, converted to",
           "copies per litre of water using the volume filtered"),
         estimatedNumberOfCopies_error      = MISSING,
         estimatedNumberOfCopies_error_type = MISSING,
         quantificationCycle = "not applicable: dPCR",
         rfu                 = MISSING,
         positivePartitions           = positives,
         acceptedPartitions           = accepted,
         partitionVolumeMicroliters   = 0.00085,
         concentrationInReaction      = conc_machine,
         concentrationInReaction_unit = "copies/uL",
         volumeFilteredLiters         = vol,
         dilutionFactor               = dilution,
         dataset                      = "field survey") %>%
  dplyr::select(assay_name, samp_name, technical_rep_id, pcr_plate_id, well_id,
                detected_notDetected,
                estimatedNumberOfCopies, estimatedNumberOfCopies_unit,
                estimatedNumberOfCopies_method, estimatedNumberOfCopies_error,
                estimatedNumberOfCopies_error_type,
                quantificationCycle, rfu,
                positivePartitions, acceptedPartitions,
                partitionVolumeMicroliters,
                concentrationInReaction, concentrationInReaction_unit,
                volumeFilteredLiters, dilutionFactor, dataset)

# =============================================================================
# 2. Decay experiment
# =============================================================================
# Only the records the decay model consumes: the DNA fraction of carboys 1-3
# over the first 30 hours, controls excluded.

decay_raw <- read.csv(here("data", "Final_decay_ddPCR_datasheet.csv"),
                      stringsAsFactors = FALSE)

decay <- decay_raw %>%
  filter(!Control, Component == "DNA", Hours_base <= 30, Carboy %in% 1:3,
         Marker %in% names(assay_of_decay_marker)) %>%
  mutate(samp_name  = sprintf("carboy%d_%.1fh", Carboy, Hours_base),
         assay_name = unname(assay_of_decay_marker[Marker]),
         copies_mL  = as.numeric(copies_mL))

decay_sample_meta <- decay %>%
  group_by(samp_name) %>%
  summarise(Carboy      = first(Carboy),
            Hours_base  = first(Hours_base),
            Time_sampled = first(Time_sampled),
            Filt_mL     = first(Filt_mL),
            Extraction  = first(Extraction),
            .groups = "drop") %>%
  arrange(Carboy, Hours_base) %>%
  transmute(
    samp_name                      = samp_name,
    materialSampleID               = Extraction,
    samp_category                  = "sample",
    eventDate                      = DECAY_EVENT_DATE,
    verbatimEventTime              = Time_sampled,
    decimalLatitude                = station_coords$lat[1],
    decimalLongitude               = station_coords$lon[1],
    verbatimCoordinateSystem       = "decimal degrees",
    verbatimSRS                    = "WGS84",
    geo_loc_name                   = "USA: Washington, Hood Canal",
    env_broad_scale                = "marine biome [ENVO:00000447]",
    env_local_scale                = "mesocosm [ENVO:00003924]",
    env_medium                     = "sea water [ENVO:00002149]",
    habitat_natural_artificial_0_1 = 1,
    minimumDepthInMeters           = NA_real_,
    maximumDepthInMeters           = NA_real_,
    samp_collect_device            = "carboy subsample",
    samp_collect_method            = paste(
      "seawater collected from a netted enclosure into a 25 L carboy,",
      "transported at ambient temperature and held in an environmental chamber",
      "at 15 C on a 16 h light cycle with aeration, then subsampled;",
      "elapsed time since collection is",
      sprintf("%.1f h", Hours_base)),
    samp_size                      = Filt_mL / 1000,
    samp_size_unit                 = "L",
    biological_rep_relation        = "the three carboys are the biological replicates",
    sample_derived_from            = sprintf("carboy %d", Carboy),
    filter_name                    = paste("Advantec mixed cellulose ester, in a",
                                           "tandem Smith-Root housing stack of",
                                           "5.0, 1.0 and 0.45 um; only the",
                                           "analysed fraction is reported here"),
    filter_material                = "mixed cellulose ester (MCE)",
    filter_diameter                = 47,
    filter_passive_active_0_1      = 1,
    samp_store_sol                 = "DNA/RNA Shield (Zymo), 2 mL",
    samp_store_temp                = -80,
    samp_store_dur                 = "within one month of collection",
    nucl_acid_ext                  = NUCL_ACID_EXT,
    nucl_acid_ext_kit              = "Quick-DNA/RNA Miniprep Kit (Zymo Research)",
    date_ext                       = MISSING,
    samp_vol_we_dna_ext            = 2,
    samp_vol_we_dna_ext_unit       = "mL of preservation buffer",
    tidal_stage                    = "not applicable: mesocosm",
    temp                           = 15,
    salinity                       = NA_real_,
    dataset                        = "decay experiment"
  )

decay_amp <- decay %>%
  arrange(Carboy, Hours_base, assay_name, PCR_rep) %>%
  transmute(
    assay_name              = assay_name,
    samp_name               = samp_name,
    technical_rep_id        = PCR_rep,
    pcr_plate_id            = MISSING,
    well_id                 = MISSING,
    detected_notDetected    = ifelse(!is.na(copies_mL) & copies_mL > 0,
                                     "detected", "notDetected"),
    estimatedNumberOfCopies = copies_mL * 1000,
    estimatedNumberOfCopies_unit   = "copies/L",
    estimatedNumberOfCopies_method = paste(
      "ddPCR Poisson estimate from partition counts, converted to copies per",
      "litre of water using the volume filtered"),
    estimatedNumberOfCopies_error      = MISSING,
    estimatedNumberOfCopies_error_type = MISSING,
    quantificationCycle     = "not applicable: dPCR",
    rfu                     = MISSING,
    positivePartitions      = Positive_droplets,
    acceptedPartitions      = Positive_droplets + Negative_droplets,
    partitionVolumeMicroliters   = 0.00085,
    concentrationInReaction      = ddPCR_copies_ul_well,
    concentrationInReaction_unit = "copies/uL",
    volumeFilteredLiters         = Filt_mL / 1000,
    dilutionFactor               = Dilution,
    dataset                      = "decay experiment"
  )

# =============================================================================
# 3. Write sampleMetadata and ampData
# =============================================================================

sample_meta <- bind_rows(field_sample_meta, decay_sample_meta)

write.csv(sample_meta,
          file.path(out_dir, sprintf("sampleMetadata_%s.csv", PROJECT_ID)),
          row.names = FALSE, na = "")

amp_all <- bind_rows(field_amp, decay_amp)

for (a in sort(unique(amp_all$assay_name))) {
  f <- file.path(out_dir, sprintf("ampData_%s_%s.csv", PROJECT_ID, a))
  write.csv(amp_all %>% filter(assay_name == a), f, row.names = FALSE, na = "")
  cat("wrote", basename(f), "-", sum(amp_all$assay_name == a), "rows\n")
}

# =============================================================================
# 4. projectMetadata
# =============================================================================
# FAIRe lays projectMetadata out with one row per term, a project_level column
# for values shared across assays, and one column per assay for values that
# differ between them.

assays <- sort(unique(amp_all$assay_name))

pm <- function(term, project_level = "", ...) {
  vals <- list(...)
  row <- c(term_name = term, project_level = as.character(project_level))
  for (a in assays) row[a] <- if (!is.null(vals[[a]])) as.character(vals[[a]]) else ""
  row
}

project_rows <- list(
  pm("project_id",   PROJECT_ID),
  pm("project_name", "Quantitative aging of environmental DNA using multiple components"),
  pm("project_contact", MISSING),
  pm("institution",     MISSING),
  pm("license",         MISSING),
  pm("recordedBy",      MISSING),
  pm("bibliographicCitation", MISSING),
  pm("associated_resource",
     "https://github.com/pedrobdfp/eDNA_quantitative_aging"),
  pm("checkls_ver",  "FAIRe checklist v1.0.2"),
  pm("study_factor",
     paste("time since shedding, estimated from the relative concentrations of",
           "four mitochondrial markers with different decay rates")),
  pm("detection_type", "targeted taxon detection"),
  pm("assay_type",     "targeted"),
  pm("amp_vis_method", "dPCR"),
  pm("platform",       "Bio-Rad QX200 Droplet Digital PCR"),
  pm("instrument",     "Bio-Rad QX200 droplet reader"),
  pm("thermocycler",   "Bio-Rad C1000 Touch"),
  pm("droplet_generator", "Bio-Rad AutoDG"),
  pm("pcr_analysis_software", "Bio-Rad QuantaSoft"),
  pm("target_taxon",   "Tursiops truncatus"),
  pm("targetTaxonomicScope", "Tursiops truncatus"),

  pm("target_gene", "",
     Tt_Cytb     = "cytb",
     Tt_16S      = "16S rRNA",
     Tt_DLL1     = "D-loop",
     Tt_longFrag = "16S rRNA to D-loop"),

  pm("assay_validation", "",
     Tt_Cytb     = DECAY_REF,
     Tt_16S      = DECAY_REF,
     Tt_DLL1     = DECAY_REF,
     Tt_longFrag = paste("derived channel: droplets positive for both the 16S",
                         "and D-loop probes in the duplex reaction, indicating",
                         "a template molecule spanning both target regions.",
                         DECAY_REF)),

  pm("pcr_primer_reference", "",
     Tt_Cytb = og("Tt_Cytb", "ref"), Tt_16S = og("Tt_16S", "ref"),
     Tt_DLL1 = og("Tt_DLL1", "ref"), Tt_longFrag = og("Tt_longFrag", "ref")),

  pm("target_length_bp", "",
     Tt_Cytb = og("Tt_Cytb", "amplicon"), Tt_16S = og("Tt_16S", "amplicon"),
     Tt_DLL1 = og("Tt_DLL1", "amplicon"),
     Tt_longFrag = og("Tt_longFrag", "amplicon")),

  pm("assay_reaction_format", "",
     Tt_Cytb     = "monoplex",
     Tt_16S      = "duplex with Tt_DLL1",
     Tt_DLL1     = "duplex with Tt_16S",
     Tt_longFrag = "double-positive channel of the Tt_16S / Tt_DLL1 duplex"),

  pm("probeReporter", "",
     Tt_Cytb = og("Tt_Cytb", "reporter"), Tt_16S = og("Tt_16S", "reporter"),
     Tt_DLL1 = og("Tt_DLL1", "reporter"),
     Tt_longFrag = og("Tt_longFrag", "reporter")),

  pm("pcr_primer_forward", "",
     Tt_Cytb = og("Tt_Cytb", "fwd"), Tt_16S = og("Tt_16S", "fwd"),
     Tt_DLL1 = og("Tt_DLL1", "fwd"), Tt_longFrag = og("Tt_longFrag", "fwd")),
  pm("pcr_primer_reverse", "",
     Tt_Cytb = og("Tt_Cytb", "rev"), Tt_16S = og("Tt_16S", "rev"),
     Tt_DLL1 = og("Tt_DLL1", "rev"), Tt_longFrag = og("Tt_longFrag", "rev")),
  pm("pcr_primer_name_forward", "",
     Tt_Cytb = og("Tt_Cytb", "fwd_name"), Tt_16S = og("Tt_16S", "fwd_name"),
     Tt_DLL1 = og("Tt_DLL1", "fwd_name"),
     Tt_longFrag = og("Tt_longFrag", "fwd_name")),
  pm("pcr_primer_name_reverse", "",
     Tt_Cytb = og("Tt_Cytb", "rev_name"), Tt_16S = og("Tt_16S", "rev_name"),
     Tt_DLL1 = og("Tt_DLL1", "rev_name"),
     Tt_longFrag = og("Tt_longFrag", "rev_name")),
  pm("pcr_primer_conc_forward", "",
     Tt_Cytb = og("Tt_Cytb", "primer_conc"),
     Tt_16S  = og("Tt_16S", "primer_conc"),
     Tt_DLL1 = og("Tt_DLL1", "primer_conc"),
     Tt_longFrag = og("Tt_longFrag", "primer_conc")),
  pm("pcr_primer_conc_reverse", "",
     Tt_Cytb = og("Tt_Cytb", "primer_conc"),
     Tt_16S  = og("Tt_16S", "primer_conc"),
     Tt_DLL1 = og("Tt_DLL1", "primer_conc"),
     Tt_longFrag = og("Tt_longFrag", "primer_conc")),
  pm("probe_seq", "",
     Tt_Cytb = og("Tt_Cytb", "probe"), Tt_16S = og("Tt_16S", "probe"),
     Tt_DLL1 = og("Tt_DLL1", "probe"),
     Tt_longFrag = og("Tt_longFrag", "probe")),
  pm("probe_conc", "250 nM"),
  pm("probeQuencher", "ZEN internal quencher with 3' Iowa Black FQ (IDT)"),
  pm("annealingTemp", "",
     Tt_Cytb = og("Tt_Cytb", "anneal"), Tt_16S = og("Tt_16S", "anneal"),
     Tt_DLL1 = og("Tt_DLL1", "anneal"),
     Tt_longFrag = og("Tt_longFrag", "anneal")),
  pm("pcr_cycles", "45"),
  pm("thermocycling_program", "",
     Tt_Cytb = paste("4 C 10 min; 95 C 10 min; 45 cycles of 94 C 30 s and",
                     "60 C 60 s (annealing/extension); 98 C 10 min; 4 C hold"),
     Tt_16S  = paste("4 C 10 min; 95 C 10 min; 45 cycles of 94 C 30 s,",
                     "56 C 30 s and 72 C 120 s; 98 C 10 min; 4 C hold"),
     Tt_DLL1 = paste("4 C 10 min; 95 C 10 min; 45 cycles of 94 C 30 s,",
                     "56 C 30 s and 72 C 120 s; 98 C 10 min; 4 C hold"),
     Tt_longFrag = paste("4 C 10 min; 95 C 10 min; 45 cycles of 94 C 30 s,",
                         "56 C 30 s and 72 C 120 s; 98 C 10 min; 4 C hold")),
  pm("amplificationReactionVolume", "22 uL"),
  pm("templateVolume", "2 uL"),
  pm("mastermix", "ddPCR Supermix for Probes (Bio-Rad), 11 uL per reaction"),
  pm("nucleic_acid_elution_volume", "50 uL"),
  pm("pcr_assay_lod",  "", Tt_Cytb = MISSING, Tt_16S = MISSING,
     Tt_DLL1 = MISSING, Tt_longFrag = MISSING),
  pm("pcr_assay_loq",  "", Tt_Cytb = MISSING, Tt_16S = MISSING,
     Tt_DLL1 = MISSING, Tt_longFrag = MISSING),

  pm("neg_cont_type",
     paste("three no-template controls per plate; extraction blanks; and, for",
           "the decay experiment, 2 L of post-bleach rinse water filtered at",
           "each timepoint. All excluded from the archived records.")),
  pm("pos_cont_type",
     paste("one per plate, DNA extracted from the target species, used to set",
           "the positive threshold. Excluded from the archived records.")),
  pm("sterilise_method", STERILISE),
  pm("nucl_acid_ext",    NUCL_ACID_EXT),
  pm("nucl_acid_ext_kit","Quick-DNA/RNA Miniprep Kit (Zymo Research)"),
  pm("dataGeneralization",
     paste("Only the records analysed in the manuscript are included: the DNA",
           "fraction of carboys 1-3 within the first 30 hours for the decay",
           "experiment, and the DNA1 extracts of field samples that passed the",
           "plate-level quality control described in",
           "code/ddPCR/01_build_field_droplets.R.")),
  pm("informationWithheld",
     "Controls and unused nucleic-acid fractions are retained in the raw plate exports.")
)

project_meta <- as.data.frame(do.call(rbind, project_rows),
                              stringsAsFactors = FALSE)

write.csv(project_meta,
          file.path(out_dir, sprintf("projectMetadata_%s.csv", PROJECT_ID)),
          row.names = FALSE, na = "")

# =============================================================================
# 5. Report
# =============================================================================

cat("\nwrote projectMetadata_", PROJECT_ID, ".csv - ", nrow(project_meta),
    " terms\n", sep = "")
cat("wrote sampleMetadata_", PROJECT_ID, ".csv - ", nrow(sample_meta),
    " samples (", sum(sample_meta$dataset == "field survey"), " field, ",
    sum(sample_meta$dataset == "decay experiment"), " decay)\n", sep = "")
cat("ampData rows total:", nrow(amp_all), "\n")

has_missing <- function(x) any(!is.na(x) & x == MISSING)

todo <- c(project_meta$term_name[apply(project_meta[, -1], 1, has_missing)],
          paste0("sampleMetadata: ",
                 paste(names(sample_meta)[sapply(sample_meta, has_missing)],
                       collapse = ", ")))

cat("\nTerms still set to '", MISSING, "':\n", sep = "")
cat(paste0("  ", todo, collapse = "\n"), "\n")
if (DECAY_EVENT_DATE == MISSING) {
  cat("\nDECAY_EVENT_DATE is still unset; edit it at the top of this script.\n")
}
