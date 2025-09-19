# =============================================================================
# Script: 01_data_preparation.R
# Purpose: Data preparation and preprocessing for nanoparticle SRFA ratio analysis
# Author: Michel Gad
# Date: 2025-09-19
# Description: 
#   - Load pre-processed formulas_avg and formulas_clean datasets
#   - Apply normalization and magnitude filtering
#   - Calculate reproducibility metrics
#   - Export processed data for downstream analysis
# =============================================================================

# Print script header information
cat("=============================================================================\n")
cat("Script: 01_data_preparation.R\n")
cat("Purpose: Data preparation and preprocessing for nanoparticle SRFA ratio analysis\n")
cat("Author: Michel Gad\n")
cat("Date: 2025-09-19\n")
cat("=============================================================================\n\n")

# --- Load Required Libraries ---
message("Loading required libraries...")
library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(htmlwidgets)

# --- Set Working Directory and Create Output Folders ---
message("Setting up directories...")
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}
if (!dir.exists("output/01_data_preparation")) {
  dir.create("output/01_data_preparation", recursive = TRUE)
}

# --- Project Initialization ---
message("Initializing project...")
project_dir <- getwd()
project_name <- "NP_SRFA_Ratio"

message(paste("Project directory:", project_dir))
message(paste("Project name:", project_name))

# --- Load Pre-processed Data Files ---
message("Loading pre-processed data files...")

# Load the averaged and clean formula datasets
formulas_avg_file <- "input/formulas_avg.csv"
formulas_clean_file <- "input/formulas_clean.csv"

if (file.exists(formulas_avg_file)) {
  formulas.avg <- read_delim(formulas_avg_file, delim = ";", show_col_types = FALSE, 
                            locale = locale(decimal_mark = ","))
  message(paste("Loaded averaged formulas:", nrow(formulas.avg), "records"))
} else {
  stop(paste("Averaged formulas file not found:", formulas_avg_file))
}

if (file.exists(formulas_clean_file)) {
  formulas.clean <- read_delim(formulas_clean_file, delim = ";", show_col_types = FALSE,
                              locale = locale(decimal_mark = ","))
  message(paste("Loaded clean formulas:", nrow(formulas.clean), "records"))
} else {
  stop(paste("Clean formulas file not found:", formulas_clean_file))
}

# --- Data Validation ---
message("Validating input data...")

# Check required columns
required_cols <- c("measurement_name", "formula_string", "peak_intensity", "group")
missing_cols <- setdiff(required_cols, names(formulas.avg))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns in formulas_avg:", paste(missing_cols, collapse = ", ")))
}

message("Data validation complete. All required columns present.")

# --- Define Sample Groups ---
message("Defining sample groups...")

# Define common prefixes for group assignment
common_prefixes <- c(
  "AuNP-HH", "P25-HH", 
  "AuNP-MM", "P25-MM",
  "AuNP-LL", "P25-LL",
  "AuNP-HM", "P25-HM", 
  "AuNP-HL", "P25-HL",
  "AuNP-MH", "P25-MH", 
  "AuNP-LH", "P25-LH",
  "BLK_P25", "BLK_AuNP"
)

# Ensure group column exists in both datasets
if (!"group" %in% names(formulas.avg)) {
  formulas.avg <- formulas.avg %>%
  mutate(
    # Extract group prefix from measurement names
    group = str_extract(
      measurement_name, 
      pattern = paste0("^(", paste(common_prefixes, collapse = "|"), ")")
    ),
    # Assign "SRFA-STD" to any remaining NA values
    group = if_else(is.na(group), "SRFA-STD", group)
  )
}

if (!"group" %in% names(formulas.clean)) {
  formulas.clean <- formulas.clean %>%
  mutate(
      # Extract group prefix from measurement names
    group = str_extract(
      measurement_name,
      pattern = paste0("^(", paste(common_prefixes, collapse = "|"), ")")
    ),
      # Assign "SRFA-STD" to any remaining NA values
    group = if_else(is.na(group), "SRFA-STD", group)
  )
}

message("Sample groups assigned successfully:")
message(paste("-", length(unique(formulas.avg$group)), "groups in averaged formulas"))
message(paste("-", length(unique(formulas.clean$group)), "groups in clean formulas"))

# --- Exclude Pristine NPs ---
message("Excluding pristine nanoparticles...")

# Remove blank samples from both datasets
formulas.avg.blk <- formulas.avg %>%
  filter(!stringr::str_detect(measurement_name, "BLK"))

formulas.clean.blk <- formulas.clean %>%
  filter(!stringr::str_detect(measurement_name, "BLK"))

message(paste("Pristine NPs excluded. Remaining formulas:"))
message(paste("- Averaged:", nrow(formulas.avg.blk)))
message(paste("- Clean:", nrow(formulas.clean.blk)))

# --- Calculate Replicate Reproducibility Metrics ---
message("Calculating replicate reproducibility metrics...")

formulas.clean <- formulas.clean %>%
  group_by(group, formula_string) %>%
  mutate(
    formula_count = n(),
    occurrence_type = if_else(formula_count > 2, "Common", "Unique"),
    RIdiff = if_else(
      occurrence_type == "Common",
      (max(peak_relint_bp, na.rm = TRUE) - min(peak_relint_bp, na.rm = TRUE)) / 
        min(peak_relint_bp, na.rm = TRUE),
      NA_real_
    ),
    Prec_RIdiff = if_else(
      occurrence_type == "Common",
      RIdiff * 100,
      NA_real_
    )
  ) %>%
  ungroup()

# Calculate reproducibility thresholds
reproducibility_threshold_window <- formulas.clean %>%
  filter(occurrence_type == "Common") %>%
  group_by(group) %>%
  summarise(
    threshold_min = quantile(peak_relint_bp, 0.025, na.rm = TRUE),
    threshold_max = quantile(peak_relint_bp, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

message(paste("Replicate analysis complete. Calculated RIdiff for", 
              sum(formulas.clean$occurrence_type == "Common", na.rm = TRUE),
              "common formulas"))

# --- Generate Summary Table ---
message("Generating summary table...")

summary_table <- formulas.clean %>%
  group_by(group) %>%
  summarise(
    Common = sum(occurrence_type == "Common", na.rm = TRUE),
    Unique = sum(occurrence_type == "Unique", na.rm = TRUE),
    Total = n(),
    Common_pct = (Common / Total) * 100,
    Unique_pct = (Unique / Total) * 100,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c(common_prefixes, "SRFA-STD"))) %>%
  arrange(group)

message(paste("Summary statistics calculated for", nrow(summary_table), "sample groups"))

# --- Process Common Formulas Across Replicates ---
message("Processing common formulas across replicates...")

# Filter for formulas that appear in all replicates
common_formulas <- formulas.clean %>% 
  filter(occurrence_type == "Common")

message(paste("Found", nrow(common_formulas), "common formulas across all replicates"))

# Aggregate intensity values across technical replicates
grouped_formulas <- common_formulas %>%
  # Create base measurement name by removing replicate suffix
  mutate(
    measurement_name_base = str_remove(measurement_name, "_\\d+$")
  ) %>%
  # Group by compound identity and sample group
  group_by(group, formula_string, measurement_name_base) %>%
  # Sum intensity values across replicates
  summarise(
    peak_intensity = sum(peak_intensity, na.rm = TRUE),
    peak_relint_tic = sum(peak_relint_tic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Ensure proper ordering of results
  arrange(group, formula_string)

message(paste("Aggregated intensities for", nrow(grouped_formulas),
              "unique formula-sample combinations"))
message("Average intensity values calculated across replicates")


# --- Normalizing to Total Sum (RI (‰)) ---
message("Normalizing to total sum (RI (‰))...")

# Step 1: Normalizing to Total Sum (RI (‰))
data_normalized <- formulas.avg %>%
  group_by(measurement_name) %>%
  mutate(
    total_intensity = sum(peak_intensity, na.rm = TRUE),
    # Calculate the RI (‰) (per thousand) RI
    rel_intensity_permille = (peak_intensity / total_intensity) * 1000
  ) %>%
  ungroup() %>%
  select(-total_intensity) # Clean up the intermediate column

message("Step 1: Normalization complete. New column 'rel_intensity_permille' created.")

# Step 2: Applying Magnitude Cutoff of 0.05 ‰
magnitude_cutoff <- 0.05
rows_before <- nrow(data_normalized)

data_mag_filtered <- data_normalized %>%
  filter(rel_intensity_permille >= magnitude_cutoff)

rows_after <- nrow(data_mag_filtered)
message(paste("Step 2: Removed", rows_before - rows_after, "peaks based on the magnitude cutoff of", magnitude_cutoff, "‰."))
message(paste("Remaining peaks:", rows_after))

formulas.avg.processed <- data_mag_filtered

message("--- FINAL PROCESSED DATASET ---")
message("Workflow complete. The final data is stored in the 'formulas.avg.processed' dataframe.")
message(paste("Total peaks remaining in the final dataset:", nrow(formulas.avg.processed)))

# --- Calculate Replicate Reproducibility Metrics ---
message("Calculating replicate reproducibility metrics...")

formulas.clean <- formulas.clean %>%
  group_by(group, formula_string) %>%
  mutate(
    formula_count = n(),
    occurrence_type = if_else(formula_count > 2, "Common", "Unique"),
    RIdiff = if_else(
      occurrence_type == "Common",
      (max(peak_relint_bp, na.rm = TRUE) - min(peak_relint_bp, na.rm = TRUE)) / 
        min(peak_relint_bp, na.rm = TRUE),
      NA_real_
    ),
    Prec_RIdiff = if_else(
      occurrence_type == "Common",
      RIdiff * 100,
      NA_real_
    )
  ) %>%
  ungroup()

# Calculate reproducibility thresholds
reproducibility_threshold_window <- formulas.clean %>%
  filter(occurrence_type == "Common") %>%
  group_by(group) %>%
  summarise(
    threshold_min = quantile(peak_relint_bp, 0.025, na.rm = TRUE),
    threshold_max = quantile(peak_relint_bp, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

message(paste("Replicate analysis complete. Calculated RIdiff for", 
              sum(formulas.clean$occurrence_type == "Common", na.rm = TRUE),
              "common formulas"))

# --- Generate Summary Table ---
message("Generating summary table...")

summary_table <- formulas.clean %>%
  group_by(group) %>%
  summarise(
    Common = sum(occurrence_type == "Common", na.rm = TRUE),
    Unique = sum(occurrence_type == "Unique", na.rm = TRUE),
    Total = n(),
    Common_pct = (Common / Total) * 100,
    Unique_pct = (Unique / Total) * 100,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c(common_prefixes, "SRFA-STD"))) %>%
  arrange(group)

message(paste("Summary statistics calculated for", nrow(summary_table), "sample groups"))

# --- Process Common Formulas Across Replicates ---
message("Processing common formulas across replicates...")

# Filter for formulas that appear in all replicates
common_formulas <- formulas.clean %>% 
  filter(occurrence_type == "Common")

message(paste("Found", nrow(common_formulas), "common formulas across all replicates"))

# Aggregate intensity values across technical replicates
grouped_formulas <- common_formulas %>%
  # Create base measurement name by removing replicate suffix
  mutate(
    measurement_name_base = str_remove(measurement_name, "_\\d+$")
  ) %>%
  # Group by compound identity and sample group
  group_by(group, formula_string, measurement_name_base) %>%
  # Sum intensity values across replicates
  summarise(
    peak_intensity = sum(peak_intensity, na.rm = TRUE),
    peak_relint_tic = sum(peak_relint_tic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Ensure proper ordering of results
  arrange(group, formula_string)

message(paste("Aggregated intensities for", nrow(grouped_formulas),
              "unique formula-sample combinations"))
message("Average intensity values calculated across replicates")

# --- Export Processed Data ---
message("Exporting processed data...")

# Save main processed datasets
write_csv(formulas.avg.processed, "output/01_data_preparation/formulas_processed.csv")
write_csv(summary_table, "output/01_data_preparation/summary_statistics.csv")
write_csv(grouped_formulas, "output/01_data_preparation/grouped_formulas.csv")
write_csv(reproducibility_threshold_window, "output/01_data_preparation/reproducibility_thresholds.csv")

# Save intermediate datasets
write_csv(formulas.clean, "output/01_data_preparation/formulas_clean.csv")
write_csv(formulas.avg, "output/01_data_preparation/formulas_averaged.csv")

message("Data preparation complete!")
message("Files saved to output/01_data_preparation/")
message("- formulas_processed.csv: Final processed dataset")
message("- summary_statistics.csv: Summary statistics by group")
message("- grouped_formulas.csv: Aggregated formulas across replicates")
message("- reproducibility_thresholds.csv: Reproducibility thresholds by group")
message("- formulas_clean.csv: Clean dataset after blank subtraction")
message("- formulas_averaged.csv: Averaged dataset across replicates")
