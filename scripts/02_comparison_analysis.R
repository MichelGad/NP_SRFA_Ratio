# =============================================================================
# Script: 02_comparison_analysis.R
# Purpose: Comparative analysis of AuNP and P25 nanoparticle formulations
# Author: Michel Gad
# Date: 2025-09-19
# Description: 
#   - Load processed molecular formula data
#   - Create separate datasets for each nanoparticle type
#   - Perform blank subtraction and comparative analysis
#   - Calculate intensity-weighted averages and generate summary statistics
# =============================================================================

# Print script header information
cat("=============================================================================\n")
cat("Script: 02_comparison_analysis.R\n")
cat("Purpose: Comparative analysis of AuNP and P25 nanoparticle formulations\n")
cat("Author: Michel Gad\n")
cat("Date: 2025-09-19\n")
cat("=============================================================================\n\n")

# --- Load Required Libraries ---
message("Loading required libraries...")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(scales)

# --- Set Working Directory and Create Output Folders ---
message("Setting up directories...")
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}
if (!dir.exists("output/02_comparison_analysis")) {
  dir.create("output/02_comparison_analysis", recursive = TRUE)
}

# --- Load Processed Data ---
message("Loading processed molecular formula data...")

# Import and clean raw data
formulas_processed <- read_csv("output/01_data_preparation/formulas_processed.csv", show_col_types = FALSE)

message(paste("Data loaded and cleaned:"))
message(paste("- File: formulas_processed.csv"))
message(paste("- Rows:", nrow(formulas_processed)))
message(paste("- Columns:", ncol(formulas_processed)))
message(paste("- Samples:", length(unique(formulas_processed$measurement_name))))

# --- Extract Nanoparticle Subsets ---
message("Creating separate datasets for each nanoparticle type...")

STD_SRFA <- formulas_processed %>% filter(str_detect(measurement_name, "STD"))
formulas_P25 <- formulas_processed %>% filter(str_detect(measurement_name, "P25"))
formulas_AuNP <- formulas_processed %>% filter(str_detect(measurement_name, "AuNP"))
BLK_P25 <- formulas_processed %>% filter(str_detect(measurement_name, "BLK_P25"))
BLK_AuNP <- formulas_processed %>% filter(str_detect(measurement_name, "BLK_AuNP"))

message("Data subsets created:")
message(paste("- SRFA standard samples:", nrow(STD_SRFA), "rows"))
message(paste("- BLK P25:", nrow(BLK_P25), "rows"))
message(paste("- BLK AuNP:", nrow(BLK_AuNP), "rows"))
message(paste("- P25 samples:", nrow(formulas_P25), "rows"))
message(paste("- AuNP samples:", nrow(formulas_AuNP), "rows"))

# --- Helper Functions for Blank Subtraction ---
message("Defining helper functions for blank subtraction...")

# Helper Function 1: Report Current Row Counts by Measurement Name
report_rows_by_measurement <- function(data_df, df_label) {
  message(paste0("--- Current Rows in ", df_label, " by measurement_name ---"))
  if (nrow(data_df) > 0) {
    # Count rows for each measurement_name
    row_counts <- data_df %>%
      group_by(measurement_name) %>%
      summarise(count = n(), .groups = 'drop')
    print(as.data.frame(row_counts))
  } else {
    message(paste0("No rows found in ", df_label, "."))
  }
  message("---")
}

# Helper Function 2: Report Rows That Will Be Removed (Pre-AntiJoin)
report_pre_anti_join_removed_rows <- function(original_data_df, block_data_df, df_label) {
  # Use semi_join to find rows in original_data_df that have matches in block_data_df
  rows_to_be_removed <- original_data_df %>%
    semi_join(block_data_df, by = "formula_string")

  message(paste0("--- Rows to be REMOVED from ", df_label, " by measurement_name (Pre-AntiJoin) ---"))
  if (nrow(rows_to_be_removed) > 0) {
    # Count how many of these "removed" rows belong to each measurement_name
    removed_summary <- rows_to_be_removed %>%
      group_by(measurement_name) %>%
      summarise(removed_count = n(), .groups = 'drop')
    print(as.data.frame(removed_summary))
  } else {
    message(paste0("No rows will be removed from ", df_label, " based on the BLK data provided."))
  }
  message("---")
}

# --- Process P25 Data ---
message("Processing P25 data...")

# Report current state of ORIGINAL formulas_P25 BEFORE anti_join
message(paste("Total rows in formulas_P25 (original) BEFORE anti_join:", nrow(formulas_P25), "rows"))
report_rows_by_measurement(formulas_P25, "formulas_P25 (ORIGINAL)")

# Report specifically WHICH rows will be removed by the upcoming anti_join
report_pre_anti_join_removed_rows(formulas_P25, BLK_P25, "formulas_P25")

# Perform the actual anti_join to remove BLK_P25 formula_strings
formulas_P25_filt <- formulas_P25 %>%
  dplyr::anti_join(BLK_P25, by = "formula_string")

# Report state of FILTERED formulas_P25_filt AFTER anti_join
message(paste("Total rows in formulas_P25_filt AFTER anti_join:", nrow(formulas_P25_filt), "rows"))
report_rows_by_measurement(formulas_P25_filt, "formulas_P25_filt (FILTERED)")

# --- Process AuNP Data ---
message("Processing AuNP data...")

# Report current state of ORIGINAL formulas_AuNP BEFORE anti_join
message(paste("Total rows in formulas_AuNP (original) BEFORE anti_join:", nrow(formulas_AuNP), "rows"))
report_rows_by_measurement(formulas_AuNP, "formulas_AuNP (ORIGINAL)")

# Report specifically WHICH rows will be removed by the upcoming anti_join
report_pre_anti_join_removed_rows(formulas_AuNP, BLK_AuNP, "formulas_AuNP")

# Perform the actual anti_join to remove BLK_AuNP formula_strings
formulas_AuNP_filt <- formulas_AuNP %>%
  dplyr::anti_join(BLK_AuNP, by = "formula_string")

# Report state of FILTERED formulas_AuNP_filt AFTER anti_join
message(paste("Total rows in formulas_AuNP_filt AFTER anti_join:", nrow(formulas_AuNP_filt), "rows"))
report_rows_by_measurement(formulas_AuNP_filt, "formulas_AuNP_filt (FILTERED)")

# --- Final Overall Summary ---
message("--- Overall Data Subsets Summary ---")
message(paste("- P25 samples (filtered):", nrow(formulas_P25_filt), "rows"))
message(paste("- AuNP samples (filtered):", nrow(formulas_AuNP_filt), "rows"))

# --- Categorize Nanoparticle Ratios ---
message("Categorizing nanoparticle ratios...")

# Identify and extract ratio patterns from measurement names
categorize_ratios <- function(data) {
  ratios <- c("-HH", "-MM", "-LL", "-MH", "-LH", "-HM", "-HL")
  
  data %>%
    mutate(
      ratio = map_chr(measurement_name, function(name) {
        # Find which ratio pattern exists in the name
        matches <- ratios[str_detect(name, fixed(ratios))]
        if (length(matches) > 0) matches[1] else NA_character_
      }),
      new_name = ifelse(!is.na(ratio), paste0("NP", ratio), NA_character_)
    )
}

# Apply to both nanoparticle types
formulas_P25 <- categorize_ratios(formulas_P25) %>% mutate(NP = "P25")
formulas_AuNP <- categorize_ratios(formulas_AuNP) %>% mutate(NP = "AuNP")

# Verify the results
message("Nanoparticles categorized by ratio:")
message(paste("- P25 ratios:", toString(na.omit(unique(formulas_P25$new_name)))))
message(paste("- AuNP ratios:", toString(na.omit(unique(formulas_AuNP$new_name)))))

# --- Combine and Analyze Nanoparticle Data ---
message("Combining and analyzing nanoparticle data...")

formulas_NP <- bind_rows(formulas_AuNP, formulas_P25) %>%
  group_by(new_name, formula_string) %>%
  mutate(
    common = case_when(
      n_distinct(NP) == 2 ~ "Common",
      n_distinct(NP) == 1 ~ NP,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(is_SRFA = formula_string %in% STD_SRFA$formula_string)

message("Combined dataset created:")
message(paste("- Total formulas:", nrow(formulas_NP)))
message(paste("- Common formulas:", sum(formulas_NP$common == "Common", na.rm = TRUE)))
message(paste("- Unique to AuNP:", sum(formulas_NP$common == "AuNP", na.rm = TRUE)))
message(paste("- Unique to P25:", sum(formulas_NP$common == "P25", na.rm = TRUE)))
message(paste("- Matching SRFA:", sum(formulas_NP$is_SRFA)))

# --- Define Visualization Parameters ---
message("Setting up visualization parameters...")

color_palette <- c(
  "Common" = "#999999",       # Common: appears in all NP types
  "AuNP" = "red",             # AuNP: for formulas that appear only in AuNP
  "P25" = "blue",             # P25: for formulas that appear only in P25
  "Unique" = "#984EA3"        # Unique: for formulas that appear only in P25
)

theme_custom <- function(base_size = 14) {
  theme_minimal(base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

message("Visualization parameters set:")
message("- Color palette defined for 4 categories")
message("- Custom ggplot theme configured")

# --- Generate Comprehensive Summary Statistics ---
message("Generating comprehensive summary statistics...")

# 1. Summary by nanoparticle ratio (Common vs Unique)
summary_common <- formulas_NP %>%
  group_by(new_name) %>%
  summarise(
    Common = sum(common == "Common", na.rm = TRUE),
    AuNP_unique = sum(common == "AuNP", na.rm = TRUE),
    P25_unique = sum(common == "P25", na.rm = TRUE),
    Total = n(),
    .groups = "drop"
  ) %>%
  mutate(across(c(Common, AuNP_unique, P25_unique), 
                ~ round(.x / Total * 100, 1), 
                .names = "{.col}_pct"))

message("Common formula summary by ratio group:")
print(summary_common, n = Inf)

# 2. SRFA comparison by ratio group
summary_SRFA <- formulas_NP %>%
  filter(!is.na(new_name)) %>%
  filter(new_name != "NA") %>%
  filter(!str_detect(new_name, "BLK")) %>%
  group_by(new_name) %>%
  summarise(
    SRFA_match = sum(is_SRFA),
    Non_SRFA = sum(!is_SRFA),
    Total = n(),
    Increased = sum(is_SRFA & rel_intensity_permille > 
                    STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                  na.rm = TRUE),
    Decreased = sum(is_SRFA & rel_intensity_permille < 
                    STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                  na.rm = TRUE),
    Unchanged = sum(is_SRFA & rel_intensity_permille == 
                    STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), 
                ~ round(.x / Total * 100, 1), 
                .names = "{.col}_pct"))

message("SRFA comparison by ratio group:")
print(summary_SRFA, n = Inf)

# 3. SRFA comparison by individual measurement
summary_SRFA_measurement <- formulas_NP %>%
  filter(!is.na(measurement_name)) %>%
  filter(measurement_name != "NA") %>%
  filter(!str_detect(measurement_name, "BLK")) %>%
  filter(!measurement_name %in% c("BLK_AuNP_1", "BLK_P25_1")) %>%
  group_by(measurement_name) %>%
  summarise(
    SRFA_match = sum(is_SRFA),
    Non_SRFA = sum(!is_SRFA),
    Total = n(),
    Increased = sum(is_SRFA & rel_intensity_permille > 
                              STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                            na.rm = TRUE),
    Decreased = sum(is_SRFA & rel_intensity_permille < 
                              STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                            na.rm = TRUE),
    Intensity_unchanged = sum(is_SRFA & rel_intensity_permille == 
                              STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)], 
                            na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), 
                ~ round(.x / Total * 100, 1), 
                .names = "{.col}_pct")) %>%
  # Add NP type for easier filtering
  left_join(distinct(formulas_NP, measurement_name, NP), by = "measurement_name")

message("SRFA comparison by individual measurement:")
print(summary_SRFA_measurement, n = Inf)

# 4. STD SRFA coverage statistics
STD_comparison_summary <- formulas_NP %>%
  group_by(measurement_name) %>%
  summarise(
    SRFA_coverage = sum(STD_SRFA$formula_string %in% formula_string),
    SRFA_missing = sum(!STD_SRFA$formula_string %in% formula_string),
    Total_SRFA = n_distinct(STD_SRFA$formula_string),
    .groups = "drop"
  ) %>%
  mutate(
    Coverage_pct = round(SRFA_coverage / Total_SRFA * 100, 1),
    Missing_pct = round(SRFA_missing / Total_SRFA * 100, 1)
  ) %>%
  # Add NP type for easier filtering
  left_join(distinct(formulas_NP, measurement_name, NP), by = "measurement_name")

message("SRFA standard coverage statistics:")
print(STD_comparison_summary, n = Inf)

# --- Create Specialized Datasets ---
message("Creating specialized datasets...")

unique_AuNP <- formulas_NP %>% filter(common == "AuNP")
unique_P25 <- formulas_NP %>% filter(common == "P25")
common_MFs <- formulas_NP %>% filter(common == "Common")

message("Specialized datasets created:")
message(paste("- Unique AuNP formulas:", nrow(unique_AuNP)))
message(paste("- Unique P25 formulas:", nrow(unique_P25)))
message(paste("- Common formulas:", nrow(common_MFs)))

# --- Calculate Intensity-Weighted Averages ---
message("Calculating intensity-weighted averages...")

# Comprehensive function to calculate weighted statistics for all specified columns
weighted_stats <- function(data) {
  # Define all columns to calculate weighted averages for
  chem_columns <- c(
    "formula_mass", "formula_hc", "formula_oc", "formula_nc", 
    "formula_sc", "formula_pc", "formula_dbe", "formula_dbe_o",
    "peak_kmd_ch2", "peak_kmd_co2", "peak_kmd_h2", "peak_kmd_nh2",
    "formula_ai", "formula_ai_mod", "formula_ai_con", "formula_xc",
    "formula_nosc", "formula_f_star", "formula_dbe_div_c",
    "formula_dbe_div_h", "formula_dbe_div_o", "peak_kmd_div_zstar",
    "peak_zstar", "C", "H", "N", "O"
  )
  
  # Calculate weighted averages with safety checks
  data %>%
    group_by(measurement_name, new_name, NP) %>%
    summarise(
      across(
        all_of(chem_columns),
        ~ {
          total_intensity <- sum(rel_intensity_permille, na.rm = TRUE)
          if (total_intensity == 0 || all(is.na(.x))) {
            NA_real_
          } else {
            sum(.x * rel_intensity_permille, na.rm = TRUE) / total_intensity
          }
        },
        .names = "weighted_{.col}"
      ),
      n_formulas = n(),
      total_intensity = sum(rel_intensity_permille, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove any columns that ended up with all NA values
    select(where(~!all(is.na(.x))))
}

# Calculate for both nanoparticle types
IWA_AuNP <- weighted_stats(formulas_AuNP)
IWA_P25 <- weighted_stats(formulas_P25)

# Verify results
message("Intensity-weighted averages calculated:")
message(paste("- AuNP samples:", nrow(IWA_AuNP), "measurements with", 
              ncol(IWA_AuNP) - 3, "weighted metrics"))
message(paste("- P25 samples:", nrow(IWA_P25), "measurements with", 
              ncol(IWA_P25) - 3, "weighted metrics"))
message(paste("- Total intensity range (AuNP):", 
              round(range(IWA_AuNP$total_intensity, na.rm = TRUE), 2)))
message(paste("- Total intensity range (P25):", 
              round(range(IWA_P25$total_intensity, na.rm = TRUE), 2)))

# --- Export Results ---
message("Exporting comparison analysis results...")

# Save summary tables
write_csv(summary_common, "output/02_comparison_analysis/summary_common.csv")
write_csv(summary_SRFA, "output/02_comparison_analysis/summary_SRFA.csv")
write_csv(summary_SRFA_measurement, "output/02_comparison_analysis/summary_SRFA_measurement.csv")
write_csv(STD_comparison_summary, "output/02_comparison_analysis/STD_SRFA_vs_measurements_summary.csv")

# Save intensity-weighted averages
write_csv(IWA_AuNP, "output/02_comparison_analysis/IWA_AuNP.csv")
write_csv(IWA_P25, "output/02_comparison_analysis/IWA_P25.csv")

# Save specialized datasets
write_csv(unique_AuNP, "output/02_comparison_analysis/unique_AuNP_formulas.csv")
write_csv(unique_P25, "output/02_comparison_analysis/unique_P25_formulas.csv")
write_csv(common_MFs, "output/02_comparison_analysis/common_formulas.csv")

# Save combined dataset
write_csv(formulas_NP, "output/02_comparison_analysis/combined_nanoparticle_data.csv")

message("Comparison analysis complete!")
message("Results saved to output/02_comparison_analysis/")
message("- summary_common.csv: Common formula summary by ratio group")
message("- summary_SRFA.csv: SRFA comparison by ratio group")
message("- summary_SRFA_measurement.csv: SRFA comparison by individual measurement")
message("- STD_SRFA_vs_measurements_summary.csv: SRFA standard coverage statistics")
message("- IWA_AuNP.csv: Intensity-weighted averages for AuNP")
message("- IWA_P25.csv: Intensity-weighted averages for P25")
message("- unique_AuNP_formulas.csv: Formulas unique to AuNP")
message("- unique_P25_formulas.csv: Formulas unique to P25")
message("- common_formulas.csv: Formulas common to both nanoparticle types")
message("- combined_nanoparticle_data.csv: Complete combined dataset")
