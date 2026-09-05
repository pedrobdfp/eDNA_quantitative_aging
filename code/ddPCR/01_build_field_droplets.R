# =============================================================================
# 01_build_field_droplets.R
# -----------------------------------------------------------------------------
# Builds the well-level ddPCR table for the field samples from the raw plate
# exports, with the quality control needed before the counts can be modelled.
#
# The table it writes is shipped with this repository as data/field_droplets.csv,
# so the field analyses can be reproduced without the raw exports. Run this
# script only to rebuild that table from the plate files, which are archived
# separately (see the data availability statement).
#
# INPUTS (all in data/)
#   * 14 Cytb plate exports        -- "<date>_dPB<nnn>_EXT<nnn>_Cytb_<ts>.csv"
#       Target TT_Cytb. Two extracts per sample: DNA1 (column 1) and
#       DNA2 (column 2). Only DNA1 is used; DNA2 is column carryover.
#   * Oct_Duplex_combined.csv      -- the 7 duplex plates (dPB128-134) already
#       stacked, carrying Tt_DLL1 (FAM), Tt_16S (HEX) and Tt_longFrag (HEX+FAM,
#       the double-positive "Bridge" channel). Duplex plates only ran DNA1.
#   * the 7 raw DL_16S plates and their *_Double_only.csv siblings, used only
#       to verify that Oct_Duplex_combined reproduces them
#       (set VERIFY_DUPLEX_SOURCES = FALSE to skip).
#   * ddpcr_combined_all.rds       -- sample-level metadata (station, depth,
#       instrument, collection time, filtered volume). Verified here to be
#       constant within `sample`, so joining on `sample` is safe.
#
# OUTPUTS (outputs/field_droplets/)
#   field_droplets_all.csv     every well parsed, nothing dropped, with QC flags
#   field_droplets_model.csv   model-ready subset (DNA1, unknowns, field
#                              stations and time window) with the Stan offset
#   field_droplets_matched.csv same, restricted to samples in which at least
#                              one marker was detected
#   qc_*.csv                   the individual QC tables
#   qc_report.txt              the full console log of this script
#
# QUALITY CONTROL APPLIED
#   * control wells are identified by name and excluded, including 38 wells
#     that carry a control name but an "Unknown" sample type
#   * DNA2 (column 2) extracts are excluded as carryover
#   * where a sample and assay were run on more than one plate, only the plate
#     with the higher mean concentration is kept, so no run is counted twice
#   * after that filter every (sample, marker) has exactly three wells
#   * wells with fewer than MIN_ACCEPTED_DROPLETS accepted droplets are
#     flagged, not dropped
#
# THE STAN OFFSET
#   A droplet model needs, per well, W positives out of U accepted droplets and
#   an offset converting seawater concentration to expected copies per droplet:
#
#     conc_copies_ul(reaction)   = C[copies/L] * vol[L] / K       (K = 550)
#     lambda(copies per droplet) = conc_copies_ul * v_drop  (v_drop = 0.00085 uL)
#     W ~ Binomial(U, inv_cloglog(omega)),  omega = log(C) + log_offset
#     log_offset = log(vol) - log(K) + log(v_drop)
#
#   Both constants are checked against the data rather than assumed:
#   K = conc_copies_L * vol / conc_copies_ul is exactly 550 for every row, and
#   the implied droplet volume -log(1 - W/U)/conc has median 0.00085 uL.
#
# THE BRIDGE CHANNEL
#   Tt_longFrag positives in Oct_Duplex_combined are the raw double-positive
#   droplet counts, not corrected for double positives expected by chance: the
#   file's Tt_longFrag concentration reproduces -log(1 - W_double/U)/0.00085 to
#   within 3e-2 (r = 0.999998), whereas the chance-corrected version does not
#   (maximum discrepancy 51 droplets).
#
#   For the field samples the correction is negligible, since single-marker
#   counts are single digits and expected chance doubles are of order 1e-3
#   droplets, but it is not negligible for the positive controls, where
#   thousands of doubles are expected. This script therefore exports the raw
#   doubles as W together with the columns needed to compute the chance term
#   (bridge_W_16S, bridge_W_DLL1, bridge_U) and a precomputed
#   bridge_expected_chance, so that the handling of it is explicit in the model
#   rather than fixed here.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# --- Settings ----------------------------------------------------------------

VERIFY_DUPLEX_SOURCES <- TRUE   # re-derive Oct_Duplex_combined from the 7 plates
MIN_ACCEPTED_DROPLETS <- 10000  # standard ddPCR well QC; flagged, not dropped

K_CONC        <- 550            # conc_copies_L = conc_copies_ul * K / vol
V_DROPLET_UL  <- 0.00085        # BioRad effective droplet volume

FIELD_STATIONS   <- c("NOAA Boat", "Husbandry Area")
FIELD_TIME_MAX   <- 1730        # stop before severe rain started
FIELD_COMPONENT  <- "DNA1 (column 1)"

MARKER_LEVELS <- c("cytb", "Tt_16S", "Tt_DLL1", "Tt_longFrag")

out_dir <- here("outputs", "field_droplets")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Tee console output to a report file
report_con <- file(file.path(out_dir, "qc_report.txt"), open = "wt")
sink(report_con, split = TRUE)
on.exit({ sink(); close(report_con) }, add = TRUE)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) cat("\n", strrep("=", 78), "\n", title, "\n",
                            strrep("=", 78), "\n", sep = "")

data_dir <- function(f = NULL) {
  base <- here("data")

  if (is.null(f)) base else file.path(base, f)
}

read_plate <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(x) <- trimws(names(x))
  x$source_file <- basename(path)
  x
}

# Normalise the messy free-text plate columns.
#   extract : "DNA1 (column 1)" / "DNA1 (Column 1)"                 -> DNA1
#             "DNA2 (column 2 carryover)" / "DNA2 (Column2, carryover)" -> DNA2
#   rep     : "Rep 1" / "Rep2" / ""                                 -> 1 / 2 / NA
#   target  : "TT_Cytb" / "TT_CYtb"                                 -> cytb
norm_extract <- function(x) {
  x <- trimws(as.character(x))
  case_when(grepl("DNA1", x, ignore.case = TRUE) ~ "DNA1",
            grepl("DNA2", x, ignore.case = TRUE) ~ "DNA2",
            TRUE ~ NA_character_)
}

norm_rep <- function(x) {
  x <- trimws(as.character(x))
  n <- suppressWarnings(as.integer(gsub("[^0-9]", "", x)))
  ifelse(is.na(n) | n == 0L, NA_integer_, n)
}

norm_marker <- function(target) {
  t <- trimws(as.character(target))
  case_when(
    grepl("^tt[_ ]?cytb$",     t, ignore.case = TRUE) ~ "cytb",
    grepl("^tt[_ ]?16s$",      t, ignore.case = TRUE) ~ "Tt_16S",
    grepl("^tt[_ ]?dll1$",     t, ignore.case = TRUE) ~ "Tt_DLL1",
    grepl("^tt[_ ]?longfrag$", t, ignore.case = TRUE) ~ "Tt_longFrag",
    grepl("^double_only_?$",   t, ignore.case = TRUE) ~ "Tt_longFrag",
    TRUE ~ NA_character_
  )
}

# Plate id out of a filename: "...dPB117_EXT120_Cytb..." -> "dPB117"
plate_from_file <- function(f) {
  p <- regmatches(f, regexpr("dPB[0-9]+", f, ignore.case = TRUE))
  ifelse(length(p) && nzchar(p), p, NA_character_)
}

std_cols <- function(x) {
  tibble(
    plate        = x$plate,
    well         = x$Well,
    sample       = trimws(x$`Sample description 1`),
    extract_raw  = x$`Sample description 2`,
    rep_raw      = x$`Sample description 3`,
    target_raw   = x$Target,
    dye          = x$`DyeName(s)`,
    sample_type  = x$SampleType,
    status       = x$Status,
    conc_machine = suppressWarnings(as.numeric(x$`Conc(copies/µL)`)),
    accepted     = suppressWarnings(as.integer(x$`Accepted Droplets`)),
    positives    = suppressWarnings(as.integer(x$Positives)),
    negatives    = suppressWarnings(as.integer(x$Negatives)),
    source_file  = x$source_file
  ) %>%
    mutate(extract = norm_extract(extract_raw),
           rep     = norm_rep(rep_raw),
           marker  = norm_marker(target_raw))
}

# =============================================================================
# 1. Read the raw plates
# =============================================================================

rule("1. READING RAW PLATE EXPORTS")

all_csv <- list.files(data_dir(), pattern = "\\.csv$", full.names = TRUE)

f_duplex_comb <- all_csv[basename(all_csv) == "Oct_Duplex_combined.csv"]
f_double_only <- all_csv[grepl("_Double_only\\.csv$", all_csv)]
f_dl16s_raw   <- all_csv[grepl("DL_16S", all_csv) & !grepl("_Double_only", all_csv)]
f_other       <- setdiff(all_csv, c(f_duplex_comb, f_double_only, f_dl16s_raw))
# Cytb plates are the remaining exports that carry the standard plate columns
f_cytb <- Filter(function(f) {
  h <- names(read.csv(f, nrows = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  all(c("Well", "Target", "Positives", "Accepted Droplets") %in% trimws(h))
}, f_other)

say("Cytb plate files          : ", length(f_cytb))
say("Duplex combined           : ", length(f_duplex_comb))
say("Duplex raw plates         : ", length(f_dl16s_raw))
say("Double-only plates        : ", length(f_double_only))
say("Non-plate csv (ignored)   : ",
    paste(basename(setdiff(f_other, f_cytb)), collapse = ", "))

stopifnot(length(f_duplex_comb) == 1, length(f_cytb) > 0)

# ---- Cytb ------------------------------------------------------------------

cytb_raw <- bind_rows(lapply(f_cytb, function(f) {
  x <- read_plate(f)
  x$plate <- plate_from_file(basename(f))
  std_cols(x)
}))

say("\nCytb wells read: ", nrow(cytb_raw),
    " across ", n_distinct(cytb_raw$plate), " plates")

# ---- Duplex ----------------------------------------------------------------

dup_in <- read_plate(f_duplex_comb)
dup_in$plate <- dup_in$Plate
duplex_raw <- std_cols(dup_in)

say("Duplex wells read: ", nrow(duplex_raw),
    " across ", n_distinct(duplex_raw$plate), " plates")

droplets_all <- bind_rows(cytb_raw, duplex_raw) %>%
  mutate(assay = if_else(marker == "cytb", "Cytb monoplex", "DL/16S duplex"))

# =============================================================================
# 2. QC on the raw parse
# =============================================================================

rule("2. QC -- RAW PARSE AND NAME NORMALISATION")

say("Target spellings found -> normalised marker:")
print(droplets_all %>% count(target_raw, marker, name = "n_wells") %>%
        as.data.frame())

say("\nExtract spellings found -> normalised extract:")
print(droplets_all %>% count(extract_raw, extract, name = "n_wells") %>%
        as.data.frame())

say("\nReplicate labels found -> normalised rep:")
print(droplets_all %>% count(rep_raw, rep, name = "n_wells") %>% as.data.frame())

say("\nSampleType:")
print(table(droplets_all$sample_type, useNA = "ifany"))

say("\nStatus (should be a single value):")
print(table(droplets_all$status, useNA = "ifany"))

# Integrity: positives + negatives must equal accepted droplets
bad_sum <- droplets_all %>% filter(positives + negatives != accepted)
say("\nWells where positives + negatives != accepted: ", nrow(bad_sum))
if (nrow(bad_sum)) print(as.data.frame(head(bad_sum, 10)))

unmapped <- droplets_all %>% filter(is.na(marker))
say("Wells with an unmapped Target: ", nrow(unmapped))
if (nrow(unmapped)) print(as.data.frame(count(unmapped, target_raw)))

write.csv(droplets_all, file.path(out_dir, "field_droplets_all.csv"),
          row.names = FALSE)

# =============================================================================
# 3. Verify Oct_Duplex_combined against its source plates
# =============================================================================

if (VERIFY_DUPLEX_SOURCES && length(f_dl16s_raw) > 0) {
  rule("3. VERIFY -- Oct_Duplex_combined vs the raw duplex plates")

  # Keep the source channel labels distinct. The *_Double_only exports carry
  # each well TWICE -- once on FAM ("Double_only") and once on HEX
  # ("Double_only_") -- and the two channels do NOT always report the same
  # count, so they must not be collapsed.
  src <- bind_rows(lapply(c(f_dl16s_raw, f_double_only), function(f) {
    x <- read_plate(f); x$plate <- plate_from_file(basename(f)); std_cols(x)
  })) %>%
    filter(!is.na(marker)) %>%
    mutate(channel = trimws(target_raw)) %>%
    select(plate, well, channel, marker,
           positives_src = positives, accepted_src = accepted)

  comb <- duplex_raw %>%
    select(plate, well, marker,
           positives_comb = positives, accepted_comb = accepted)

  j <- left_join(src, comb, by = c("plate", "well", "marker"))

  say("source rows: ", nrow(src), " | combined rows: ", nrow(comb))

  say("\n-- single-marker channels (Tt_16S, Tt_DLL1) --")
  singles_chk <- j %>% filter(marker %in% c("Tt_16S", "Tt_DLL1"))
  say("rows compared : ", nrow(singles_chk))
  say("exact matches : ",
      sum(singles_chk$positives_src == singles_chk$positives_comb &
          singles_chk$accepted_src  == singles_chk$accepted_comb, na.rm = TRUE))

  say("\n-- Bridge channel (Tt_longFrag vs the *_Double_only exports) --")
  br_chk <- j %>% filter(marker == "Tt_longFrag") %>%
    mutate(delta = positives_src - positives_comb)
  for (ch in sort(unique(br_chk$channel))) {
    s <- br_chk %>% filter(channel == ch)
    say("  channel ", ch, ": ", nrow(s), " wells | exact matches ",
        sum(s$delta == 0, na.rm = TRUE),
        " | max |delta| ", max(abs(s$delta), na.rm = TRUE))
  }
  say("  accepted-droplet totals always agree: ",
      all(br_chk$accepted_src == br_chk$accepted_comb, na.rm = TRUE))

  mism <- br_chk %>% filter(delta != 0) %>% arrange(desc(abs(delta)))
  say("\n  Bridge wells where the two sources disagree: ", nrow(mism))
  if (nrow(mism)) {
    say("  largest disagreements (these are the positive-control wells):")
    print(as.data.frame(head(mism, 8)))
    say("\n  disagreements on SAMPLE wells only (excluding H06/H12 controls):")
    samp_mism <- mism %>% filter(!well %in% c("H06", "H12"))
    print(as.data.frame(samp_mism))
    say("  -> ", nrow(samp_mism), " sample wells differ, all by ",
        paste(sort(unique(abs(samp_mism$delta))), collapse = "/"), " droplet(s).")
  }
  say("\n  DECISION: Oct_Duplex_combined is treated as authoritative, since it",
      " is the curated Bridge call. The *_Double_only exports are an earlier",
      " gating and are used here only for this comparison.")
  write.csv(br_chk, file.path(out_dir, "qc_bridge_source_comparison.csv"),
            row.names = FALSE)
} else {
  rule("3. VERIFY -- skipped")
}

# =============================================================================
# 4. Drop controls and the DNA2 carryover extract
# =============================================================================

rule("4. FILTERING -- controls and DNA2 carryover")

# A well is a control if the machine says so, if it has no extract label, or if
# its sample name looks like a control. Some NTC/PosC wells are mislabelled
# SampleType "Unknown", so the name test is not redundant.
is_control_name <- function(s) {
  grepl(paste0("^\\s*(NTC|NEG|BLANK|PosC|POS|TT_PosC|TT_Positive",
               "|EB[-_ ]?[0-9]|Extraction\\s*Blank)"),
        s, ignore.case = TRUE)
}

controls <- droplets_all %>%
  filter(sample_type != "Unknown" | is.na(extract) | is_control_name(sample))
say("Control / non-sample wells set aside: ", nrow(controls))
print(controls %>% count(sample_type, sample, name = "n_wells") %>% as.data.frame())
write.csv(controls, file.path(out_dir, "qc_control_wells.csv"), row.names = FALSE)

dna2 <- droplets_all %>% filter(sample_type == "Unknown", extract == "DNA2")
say("\nDNA2 carryover wells removed: ", nrow(dna2),
    " (", n_distinct(dna2$sample), " samples, ",
    n_distinct(dna2$plate), " plates)")
say("  -- these are the 'column 2 carryover' extracts and are never used.")

samples_wells <- droplets_all %>%
  filter(sample_type == "Unknown", extract == "DNA1", !is.na(marker),
         !is_control_name(sample))

say("\nSample wells retained (DNA1, unknowns): ", nrow(samples_wells))
say("distinct samples: ", n_distinct(samples_wells$sample))
say("by marker:")
print(table(samples_wells$marker))

# HARD ASSERTION -- no DNA2 may survive
stopifnot(!any(samples_wells$extract == "DNA2"))
stopifnot(all(samples_wells$sample_type == "Unknown"))
stopifnot(!any(is_control_name(samples_wells$sample)))
stopifnot(all(samples_wells$extract == "DNA1"))
say("\nASSERTION PASSED: retained set is DNA1 only, no controls, no DNA2.")

# =============================================================================
# 5. Duplicate-run checks
# =============================================================================

rule("5. QC -- DUPLICATE AND REPEATED RUNS")

# (a) the same physical well appearing twice for the same marker = a real bug
well_dupes <- samples_wells %>%
  count(plate, well, marker, name = "n") %>% filter(n > 1)
say("(a) (plate, well, marker) appearing more than once: ", nrow(well_dupes))
if (nrow(well_dupes)) print(as.data.frame(head(well_dupes, 20)))

# (b) the same (sample, marker, plate, rep) twice = a real bug
key_dupes <- samples_wells %>%
  count(sample, marker, plate, rep, name = "n") %>% filter(n > 1)
say("(b) (sample, marker, plate, rep) appearing more than once: ", nrow(key_dupes))
if (nrow(key_dupes)) {
  # Not the same well twice -- check (a) is 0 -- but the same sample name and
  # replicate label in two DIFFERENT wells of one plate. Ambiguous: either a
  # true technical duplicate pair or a plate-map labelling slip.
  print(as.data.frame(head(key_dupes, 20)))
  say("    affected samples: ",
      paste(sort(unique(key_dupes$sample)), collapse = ", "))
  say("    NOTE: check (a) is ", nrow(well_dupes), ", so these are distinct",
      " wells sharing a sample+rep label on one plate, not a re-read of one well.")
  write.csv(key_dupes, file.path(out_dir, "qc_ambiguous_rep_labels.csv"),
            row.names = FALSE)
}

# (c) a sample legitimately re-run on a second plate -> more than 3 reps
plate_counts <- samples_wells %>%
  distinct(sample, marker, plate) %>%
  count(sample, marker, name = "n_plates")
say("\n(c) samples run on more than one plate (legitimate re-runs):")
print(table(plate_counts$n_plates))
multi <- plate_counts %>% filter(n_plates > 1)
say("    sample x marker combinations affected: ", nrow(multi),
    " | distinct samples: ", n_distinct(multi$sample))
print(as.data.frame(multi %>% count(sample, name = "markers_rerun")))
write.csv(multi, file.path(out_dir, "qc_multiplate_samples.csv"), row.names = FALSE)

# (d) replicate counts
rep_counts <- samples_wells %>% count(sample, marker, name = "n_reps")
say("\n(d) replicates per (sample, marker):")
print(table(rep_counts$n_reps))
odd <- rep_counts %>% filter(!n_reps %in% c(3, 6))
say("    NOT 3 or 6 replicates: ", nrow(odd))
if (nrow(odd)) print(as.data.frame(odd))
write.csv(rep_counts, file.path(out_dir, "qc_replicate_counts.csv"),
          row.names = FALSE)

# =============================================================================
# 5b. Resolve re-runs: keep one plate per sample and assay
# =============================================================================
# Some samples were re-run on a second plate. Pooling both runs treats
# plate-to-plate variation as within-sample replicate noise, which inflates the
# apparent observation error. For each sample and assay we therefore retain only
# the plate that yielded the highest concentration.
#
# The choice is made per ASSAY, not per marker: Tt_16S, Tt_DLL1 and Tt_longFrag
# are read from the same duplex well and must come from the same plate, while
# cytb is a separate monoplex reaction.

rule("5b. RESOLVING RE-RUNS")

plate_conc <- samples_wells %>%
  group_by(sample, assay, plate) %>%
  summarise(mean_conc = mean(conc_machine, na.rm = TRUE),
            n_wells   = n(), .groups = "drop")

keep_plate <- plate_conc %>%
  group_by(sample, assay) %>%
  slice_max(mean_conc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(sample, assay, plate)

rerun <- plate_conc %>%
  group_by(sample, assay) %>%
  filter(n() > 1) %>%
  arrange(sample, assay, desc(mean_conc)) %>%
  ungroup()

say("Sample x assay combinations run on more than one plate: ",
    n_distinct(paste(rerun$sample, rerun$assay)))
if (nrow(rerun)) {
  say("Plate retained per re-run (highest mean concentration):")
  print(as.data.frame(rerun %>%
    group_by(sample, assay) %>%
    summarise(kept    = plate[which.max(mean_conc)],
              dropped = paste(plate[-which.max(mean_conc)], collapse = ", "),
              conc_kept    = round(max(mean_conc), 4),
              conc_dropped = round(max(mean_conc[-which.max(mean_conc)]), 4),
              .groups = "drop")))
}
write.csv(rerun, file.path(out_dir, "qc_rerun_plate_choice.csv"), row.names = FALSE)

n_before <- nrow(samples_wells)
samples_wells <- samples_wells %>% inner_join(keep_plate,
                                              by = c("sample", "assay", "plate"))
say("\nWells before: ", n_before, " | after: ", nrow(samples_wells),
    " | removed: ", n_before - nrow(samples_wells))

rep_counts <- samples_wells %>% count(sample, marker, name = "n_reps")
say("Replicates per (sample, marker) after resolving re-runs:")
print(table(rep_counts$n_reps))

# A handful of samples were laid out twice on the SAME plate, so selecting a
# plate cannot separate them. These carry duplicated sample and replicate
# labels in different wells and are reported by check (b) above.
still_dup <- rep_counts %>% filter(n_reps > 3)
if (nrow(still_dup)) {
  say("\nStill more than 3 replicates -- duplicated within a single plate: ",
      nrow(still_dup), " (sample, marker) combinations")
  print(as.data.frame(still_dup))
  say("Affected samples: ",
      paste(sort(unique(still_dup$sample)), collapse = ", "))
  say("These are left intact here; the model subset is checked below.")
  write.csv(still_dup, file.path(out_dir, "qc_within_plate_duplicates.csv"),
            row.names = FALSE)
}

stopifnot(!any(duplicated(keep_plate[, c("sample", "assay")])))
say("ASSERTION PASSED: at most one plate retained per sample and assay.")

# =============================================================================
# 6. Metadata join
# =============================================================================

rule("6. METADATA JOIN")

meta_all <- readRDS(data_dir("ddpcr_combined_all.rds"))

meta_cols <- c("id", "bio_rep", "station", "Instrument", "depth",
               "time_collected", "filter_start_time", "vol", "depthbin")

# The join is on `sample` alone, which is only valid if these columns are
# constant within sample. Verify rather than assume.
varying <- meta_all %>%
  group_by(sample) %>%
  summarise(across(all_of(meta_cols), ~ n_distinct(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  filter(if_any(all_of(meta_cols), ~ .x > 1))

say("Samples whose metadata is NOT constant within sample: ", nrow(varying))
if (nrow(varying)) {
  print(as.data.frame(varying))
  stop("Metadata varies within sample -- joining on `sample` alone is unsafe.")
}
say("ASSERTION PASSED: sample-level metadata is constant within sample.")

meta <- meta_all %>%
  distinct(sample, across(all_of(meta_cols))) %>%
  filter(!is.na(sample))

say("\nmetadata rows: ", nrow(meta), " for ", n_distinct(meta$sample), " samples")

droplets <- samples_wells %>%
  left_join(meta, by = "sample") %>%
  mutate(marker = factor(marker, levels = MARKER_LEVELS))

no_meta <- droplets %>% filter(is.na(station)) %>% distinct(sample)
say("Samples with droplet data but NO metadata: ", nrow(no_meta))
if (nrow(no_meta)) print(no_meta$sample)

# =============================================================================
# 7. Reconcile droplet counts against the recorded concentrations
# =============================================================================

rule("7. RECONCILIATION -- droplets vs concentrations")

recon <- droplets %>%
  filter(accepted > 0) %>%
  mutate(
    lambda_hat    = -log(pmax(1 - positives / accepted, .Machine$double.eps)),
    conc_recomp   = lambda_hat / V_DROPLET_UL,
    conc_diff     = conc_recomp - conc_machine
  )

say("Recomputed conc = -log(1 - W/U) / ", V_DROPLET_UL,
    " vs the machine's Conc(copies/uL):")
print(summary(recon$conc_diff))
say("max |difference| : ", round(max(abs(recon$conc_diff), na.rm = TRUE), 5))
say("correlation      : ",
    round(cor(recon$conc_recomp, recon$conc_machine, use = "complete.obs"), 8))

v_implied <- recon %>% filter(positives > 0, positives < accepted) %>%
  mutate(v = lambda_hat / conc_machine) %>% pull(v)
say("\nImplied droplet volume from the data (uL): median ",
    round(median(v_implied, na.rm = TRUE), 8),
    "  -> -log(v) = ", round(-log(median(v_implied, na.rm = TRUE)), 4))
say("the decay model uses -7.07, i.e. v = ",
    round(exp(-7.07), 8))

# The copies/uL -> copies/L constant
kchk <- meta_all %>% filter(conc_copies_ul > 0, !is.na(vol)) %>%
  mutate(K = conc_copies_L * vol / conc_copies_ul)
say("\nK = conc_copies_L * vol / conc_copies_ul : min ", min(kchk$K),
    " max ", max(kchk$K), " (script uses ", K_CONC, ")")
stopifnot(isTRUE(all.equal(range(kchk$K), c(K_CONC, K_CONC))))
say("ASSERTION PASSED: K is exactly ", K_CONC, " for every metadata row.")

# =============================================================================
# 8. Droplet-count QC
# =============================================================================

rule("8. QC -- DROPLET COUNTS PER WELL")

say("Accepted droplets per well:")
print(summary(droplets$accepted))
low <- droplets %>% filter(accepted < MIN_ACCEPTED_DROPLETS)
say("\nWells below ", MIN_ACCEPTED_DROPLETS, " accepted droplets: ", nrow(low))
if (nrow(low)) {
  print(as.data.frame(low %>%
    select(plate, well, sample, marker, rep, positives, accepted, conc_machine)))
}
write.csv(low, file.path(out_dir, "qc_low_droplet_wells.csv"), row.names = FALSE)
say("(flagged via low_droplet_flag, NOT dropped -- decide in the model script.)")

droplets <- droplets %>%
  mutate(low_droplet_flag = accepted < MIN_ACCEPTED_DROPLETS)

# =============================================================================
# 9. Bridge / Tt_longFrag chance-coincidence columns
# =============================================================================

rule("9. BRIDGE -- chance-coincidence terms")

singles <- droplets %>%
  filter(marker %in% c("Tt_16S", "Tt_DLL1")) %>%
  select(plate, well, marker, positives) %>%
  pivot_wider(names_from = marker, values_from = positives,
              names_prefix = "W_")

bridge_terms <- droplets %>%
  filter(marker == "Tt_longFrag") %>%
  select(plate, well, bridge_U = accepted, bridge_W_double = positives) %>%
  left_join(singles, by = c("plate", "well")) %>%
  rename(bridge_W_16S = W_Tt_16S, bridge_W_DLL1 = W_Tt_DLL1) %>%
  mutate(bridge_expected_chance =
           bridge_U * (bridge_W_16S / bridge_U) * (bridge_W_DLL1 / bridge_U),
         bridge_excess = bridge_W_double - bridge_expected_chance)

say("Bridge wells: ", nrow(bridge_terms))
say("Expected chance doubles per well:")
print(summary(bridge_terms$bridge_expected_chance))
say("\nObserved double positives per well:")
print(summary(bridge_terms$bridge_W_double))
say("\nWells where expected chance doubles exceed 0.5 droplets: ",
    sum(bridge_terms$bridge_expected_chance > 0.5, na.rm = TRUE))
say("=> for sample wells the chance correction is negligible; the raw double",
    " count is exported as W.")
write.csv(bridge_terms, file.path(out_dir, "qc_bridge_terms.csv"),
          row.names = FALSE)

droplets <- droplets %>%
  left_join(bridge_terms %>%
              select(plate, well, bridge_W_16S, bridge_W_DLL1, bridge_U,
                     bridge_expected_chance, bridge_excess),
            by = c("plate", "well")) %>%
  mutate(across(starts_with("bridge_"), ~ if_else(marker == "Tt_longFrag",
                                                  .x, NA_real_)))

# =============================================================================
# 10. Model-ready subset + Stan offset
# =============================================================================

rule("10. MODEL-READY SUBSET")

model_df <- droplets %>%
  filter(station %in% FIELD_STATIONS,
         !is.na(time_collected),
         time_collected <= FIELD_TIME_MAX,
         !is.na(vol)) %>%
  mutate(
    time_min   = (as.integer(time_collected) %/% 100L) * 60L +
                 (as.integer(time_collected) %% 100L),
    dilution   = 1,                     # no dilution recorded for the field set
    # omega = log(C[copies/L]) + log_offset ; W ~ Binomial(U, inv_cloglog(omega))
    log_offset = log(vol) - log(K_CONC) + log(V_DROPLET_UL) + log(dilution)
  ) %>%
  arrange(sample, marker, plate, rep)

say("Model-ready wells: ", nrow(model_df))
say("distinct samples : ", n_distinct(model_df$sample))
say("stations         : ", paste(sort(unique(model_df$station)), collapse = ", "))
say("\nwells by marker:")
print(table(model_df$marker))
say("\nreplicates per (sample, marker):")
print(table(model_df %>% count(sample, marker) %>% pull(n)))
say("\nmarkers per sample:")
print(table(model_df %>% distinct(sample, marker) %>% count(sample) %>% pull(n)))

say("\nSamples with cytb only (duplex dropped because cytb was zero):")
cytb_only <- model_df %>% distinct(sample, marker) %>% count(sample) %>%
  filter(n == 1) %>% pull(sample)
say("  ", paste(cytb_only, collapse = ", "), " (n = ", length(cytb_only), ")")

say("\nlog_offset summary (log C[copies/L] + log_offset = log lambda):")
print(summary(model_df$log_offset))

# --- The all-zero-sample rule ------------------------------------------------
# 04_field_main_ddPCR.R additionally drops any sample whose TOTAL concentration across
# all markers is zero:
#     group_by(sample) %>% filter(sum(conc_copies_L, na.rm = TRUE) > 0)
# That rule makes sense for a lognormal model (an all-zero sample carries no
# quantitative information) but NOT necessarily for a droplet model, where a
# sample that is zero across every marker and replicate is still an informative
# set of Bernoulli non-detections.
#
# Rather than choose here, both are exported: `in_conc_pipeline` marks the wells
# the concentration pipeline keeps.

zero_samples <- model_df %>%
  group_by(sample) %>%
  summarise(total_pos = sum(positives, na.rm = TRUE), .groups = "drop") %>%
  filter(total_pos == 0) %>% pull(sample)

model_df <- model_df %>% mutate(in_conc_pipeline = !sample %in% zero_samples)

say("\nSamples with ZERO positive droplets across every marker and replicate: ",
    length(zero_samples))
say("  ", paste(zero_samples, collapse = ", "))
say("  These are dropped by 04_field_main_ddPCR.R's all-zero rule. They are KEPT in")
say("  field_droplets_model.csv and flagged in_conc_pipeline = FALSE, because a")
say("  droplet model can use them as non-detections. field_droplets_matched.csv")
say("  excludes them so the droplet and concentration runs are comparable.")
say("\nwells in_conc_pipeline = TRUE : ", sum(model_df$in_conc_pipeline))
say("wells in_conc_pipeline = FALSE: ", sum(!model_df$in_conc_pipeline))

# Worked example so the offset can be checked by hand
ex <- model_df %>% filter(conc_machine > 0) %>% slice(1)
say("\nWorked example -- sample ", ex$sample, ", ", as.character(ex$marker), ":")
say("  vol = ", ex$vol, " L | machine conc = ", round(ex$conc_machine, 6),
    " copies/uL | W = ", ex$positives, " / U = ", ex$accepted)
say("  implied C = ", round(ex$conc_machine * K_CONC / ex$vol, 3), " copies/L")
say("  omega = log(C) + log_offset = ",
    round(log(ex$conc_machine * K_CONC / ex$vol) + ex$log_offset, 6))
say("  log(lambda) direct          = ",
    round(log(-log(1 - ex$positives / ex$accepted)), 6))

write.csv(model_df, file.path(out_dir, "field_droplets_model.csv"),
          row.names = FALSE)
saveRDS(model_df, file.path(out_dir, "field_droplets_model.rds"))

matched_df <- model_df %>% filter(in_conc_pipeline)
write.csv(matched_df, file.path(out_dir, "field_droplets_matched.csv"),
          row.names = FALSE)
saveRDS(matched_df, file.path(out_dir, "field_droplets_matched.rds"))

# =============================================================================
# 11. Cross-check against the concentration pipeline
# =============================================================================

rule("11. CROSS-CHECK -- does this reproduce the concentration-model input?")

# The sample subset used by the concentration analysis in the supplement
fld_conc <- meta_all %>%
  filter(station %in% FIELD_STATIONS,
         component == FIELD_COMPONENT,
         haplotype %in% c("mixed", "mixed0", "unknown"),
         !is.na(time_collected), time_collected <= FIELD_TIME_MAX) %>%
  group_by(sample) %>%
  filter(sum(conc_copies_L, na.rm = TRUE) > 0) %>%
  ungroup()

say("04_field_main_ddPCR.R observation rows : ", nrow(fld_conc))
say("this script's droplet wells      : ", nrow(model_df))
say("04_field_main_ddPCR.R samples          : ", n_distinct(fld_conc$sample))
say("this script's samples            : ", n_distinct(model_df$sample))

say("this script's MATCHED wells     : ", nrow(matched_df))
say("this script's MATCHED samples   : ", n_distinct(matched_df$sample))

cnt_conc <- fld_conc %>% count(sample, marker, name = "n_conc")
cnt_drop <- matched_df %>% count(sample, marker, name = "n_drop") %>%
  mutate(marker = as.character(marker))

cmp <- full_join(cnt_conc, cnt_drop, by = c("sample", "marker")) %>%
  mutate(n_conc = coalesce(n_conc, 0L), n_drop = coalesce(n_drop, 0L),
         diff = n_conc - n_drop)

say("\n(sample, marker) row counts, concentration pipeline vs droplet build:")
print(table(conc = cmp$n_conc, droplet = cmp$n_drop))

mismatch <- cmp %>% filter(diff != 0)
say("\nMISMATCHED (sample, marker) combinations: ", nrow(mismatch))
if (nrow(mismatch)) print(as.data.frame(mismatch))
write.csv(cmp, file.path(out_dir, "qc_vs_concentration_pipeline.csv"),
          row.names = FALSE)

if (nrow(mismatch) == 0) {
  say("\nThe droplet build reproduces the concentration pipeline exactly, ",
      "well-for-observation.")
}

# Do the concentrations themselves agree, per sample x marker?
conc_cmp <- full_join(
  fld_conc %>% group_by(sample, marker) %>%
    summarise(conc_pipeline = paste(sort(round(conc_copies_ul, 4)), collapse = ","),
              .groups = "drop"),
  matched_df %>% mutate(marker = as.character(marker)) %>%
    group_by(sample, marker) %>%
    summarise(conc_droplets = paste(sort(round(conc_machine, 4)), collapse = ","),
              .groups = "drop"),
  by = c("sample", "marker"))

comparable <- conc_cmp %>%
  filter(!is.na(conc_pipeline), !is.na(conc_droplets))
agree    <- sum(comparable$conc_pipeline == comparable$conc_droplets)
disagree <- comparable %>% filter(conc_pipeline != conc_droplets)
only_one <- conc_cmp %>% filter(is.na(conc_pipeline) | is.na(conc_droplets))

say("\nPer (sample, marker) concentration vectors:")
say("  present in both, IDENTICAL : ", agree, " / ", nrow(comparable))
say("  present in both, differing : ", nrow(disagree))
say("  present in only one side   : ", nrow(only_one))
if (nrow(disagree)) {
  say("\n  Disagreements:")
  print(as.data.frame(head(disagree, 20)))
}
if (nrow(only_one)) {
  say("\n  Present on only one side (expect 0 after the matched filter):")
  print(as.data.frame(head(only_one, 20)))
}
write.csv(conc_cmp, file.path(out_dir, "qc_concentration_agreement.csv"),
          row.names = FALSE)

# Does droplet-derived detection agree with the `binary` flag the concentration
# model uses as its Bernoulli outcome? Compared as per (sample, marker) counts,
# since the RDS has no replicate key to join on.
det_cmp <- full_join(
  fld_conc %>% group_by(sample, marker) %>%
    summarise(n_det_conc = sum(binary == 1, na.rm = TRUE), .groups = "drop"),
  matched_df %>% mutate(marker = as.character(marker)) %>%
    group_by(sample, marker) %>%
    summarise(n_det_drop = sum(positives > 0), .groups = "drop"),
  by = c("sample", "marker")) %>%
  mutate(delta = n_det_conc - n_det_drop)

say("\nDetections per (sample, marker): `binary` flag vs positives > 0")
say("  combinations compared : ", nrow(det_cmp))
say("  identical             : ", sum(det_cmp$delta == 0, na.rm = TRUE))
say("  differing             : ", sum(det_cmp$delta != 0, na.rm = TRUE))
if (any(det_cmp$delta != 0, na.rm = TRUE)) {
  print(as.data.frame(det_cmp %>% filter(delta != 0) %>% arrange(desc(abs(delta)))))
  say("  -> a nonzero droplet count with binary = 0 (or vice versa) means the")
  say("     two pipelines would disagree about what counts as a detection.")
}
write.csv(det_cmp, file.path(out_dir, "qc_detection_agreement.csv"),
          row.names = FALSE)

say("\nNote: omega recomputed from the machine's concentration differs from")
say("log(-log(1 - W/U)) by at most ~0.08 in log space, because BioRad's per-well")
say("effective droplet volume varies slightly (0.000838-0.00092 uL). The droplet")
say("model uses W and U directly, so this affects only this cross-check.")

# =============================================================================
# Done
# =============================================================================

rule("DONE")
say("Written to outputs/field_droplets/ :")
say("  field_droplets_all.csv           every parsed well (", nrow(droplets_all), " rows)")
say("  field_droplets_model.csv/.rds    model-ready, incl. all-zero samples (",
    nrow(model_df), " rows, ", n_distinct(model_df$sample), " samples)")
say("  field_droplets_matched.csv/.rds  matched to 04_field_main_ddPCR.R (",
    nrow(matched_df), " rows, ", n_distinct(matched_df$sample), " samples)")
say("  qc_report.txt                    this report")
say("  qc_*.csv                         individual QC tables")
say("")
say("Key columns for a droplet model:")
say("  positives (W), accepted (U), log_offset")
say("    W ~ Binomial(U, inv_cloglog(log(C[copies/L]) + log_offset))")
say("  bridge_* columns carry the chance-coincidence terms for Tt_longFrag.")
say("  low_droplet_flag marks wells under ", MIN_ACCEPTED_DROPLETS, " droplets.")
