# =============================================================================
# Script: 04_export_results.R
# Purpose: Export final results and create comprehensive summary reports for nanoparticle analysis
# Author: Michel Gad
# Date: 2025-09-19
# Description: 
#   - Compile all analysis results into final export files
#   - Create summary reports and metadata files
#   - Export data in multiple formats for downstream analysis
#   - Generate analysis documentation
# =============================================================================

# Print script header information
cat("=============================================================================\n")
cat("Script: 04_export_results.R\n")
cat("Purpose: Export final results and create comprehensive summary reports for nanoparticle analysis\n")
cat("Author: Michel Gad\n")
cat("Date: 2025-09-19\n")
cat("=============================================================================\n\n")

# --- Load Required Libraries ---
message("Loading required libraries...")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(knitr)
library(rmarkdown)

# --- Set Working Directory and Create Output Folders ---
message("Setting up directories...")
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}
if (!dir.exists("output/04_export_results")) {
  dir.create("output/04_export_results", recursive = TRUE)
}
if (!dir.exists("output/04_export_results/data")) {
  dir.create("output/04_export_results/data", recursive = TRUE)
}
if (!dir.exists("output/04_export_results/reports")) {
  dir.create("output/04_export_results/reports", recursive = TRUE)
}

# --- Load Analysis Results for Summary Creation ---
message("Loading analysis results for summary creation...")

# Load only the data needed for creating new summary files
formulas_NP <- read_csv("output/02_comparison_analysis/combined_nanoparticle_data.csv", show_col_types = FALSE)
summary_statistics <- read_csv("output/01_data_preparation/summary_statistics.csv", show_col_types = FALSE)

message("Analysis results loaded for summary creation")

# --- Create Analysis Summary ---
message("Creating analysis summary...")

# Basic statistics
total_formulas <- nrow(formulas_NP)
total_measurements <- length(unique(formulas_NP$measurement_name))
n_nanoparticle_types <- 2  # AuNP and P25
n_ratio_combinations <- length(unique(formulas_NP$new_name))

# Create summary table
analysis_summary <- data.frame(
  Parameter = c(
    "Total Molecular Formulas",
    "Total Measurements",
    "Nanoparticle Types",
    "Ratio Combinations",
    "Analysis Date",
    "Scripts Used"
  ),
  Value = c(
    as.character(total_formulas),
    as.character(total_measurements),
    as.character(n_nanoparticle_types),
    as.character(n_ratio_combinations),
    as.character(Sys.Date()),
    "01_data_preparation.R, 02_comparison_analysis.R, 03_visualization.R, 04_export_results.R"
  )
)

write_csv(analysis_summary, "output/04_export_results/analysis_summary.csv")
message("Analysis summary saved")

# --- Copy Core Datasets (avoiding duplicate exports) ---
message("Copying core datasets from previous scripts...")

# Copy files from data_prep output
file.copy("output/01_data_preparation/formulas_processed.csv", "output/04_export_results/data/formulas_processed.csv", overwrite = TRUE)
file.copy("output/01_data_preparation/reproducibility_thresholds.csv", "output/04_export_results/data/reproducibility_thresholds.csv", overwrite = TRUE)
file.copy("output/01_data_preparation/grouped_formulas.csv", "output/04_export_results/data/grouped_formulas.csv", overwrite = TRUE)
message("Data preparation files copied")

# Copy files from comparison_analysis output
file.copy("output/02_comparison_analysis/combined_nanoparticle_data.csv", "output/04_export_results/data/combined_nanoparticle_data.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/IWA_AuNP.csv", "output/04_export_results/data/IWA_AuNP.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/IWA_P25.csv", "output/04_export_results/data/IWA_P25.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/summary_common.csv", "output/04_export_results/data/summary_common.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/summary_SRFA.csv", "output/04_export_results/data/summary_SRFA.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/summary_SRFA_measurement.csv", "output/04_export_results/data/summary_SRFA_measurement.csv", overwrite = TRUE)
file.copy("output/02_comparison_analysis/STD_SRFA_vs_measurements_summary.csv", "output/04_export_results/data/STD_comparison_summary.csv", overwrite = TRUE)
message("Comparison analysis files copied")

# --- Create Nanoparticle Comparison Summary ---
message("Creating nanoparticle comparison summary...")

# Detailed comparison information
comparison_summary <- formulas_NP %>%
  group_by(NP, new_name) %>%
  summarise(
    n_formulas = n(),
    common_formulas = sum(common == "Common", na.rm = TRUE),
    unique_formulas = sum(common != "Common", na.rm = TRUE),
    srfa_matches = sum(is_SRFA, na.rm = TRUE),
    mean_intensity = mean(rel_intensity_permille, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(NP, new_name)

write_csv(comparison_summary, "output/04_export_results/data/nanoparticle_comparison_summary.csv")
message("Nanoparticle comparison summary exported")

# --- Create SRFA Analysis Summary ---
message("Creating SRFA analysis summary...")

# SRFA coverage analysis
srfa_analysis <- formulas_NP %>%
  group_by(measurement_name, NP) %>%
  summarise(
    total_formulas = n(),
    srfa_formulas = sum(is_SRFA, na.rm = TRUE),
    unique_formulas = sum(!is_SRFA, na.rm = TRUE),
    srfa_coverage_pct = round(srfa_formulas / total_formulas * 100, 1),
    unique_pct = round(unique_formulas / total_formulas * 100, 1),
    .groups = 'drop'
  ) %>%
  arrange(NP, measurement_name)

write_csv(srfa_analysis, "output/04_export_results/data/srfa_analysis_summary.csv")
message("SRFA analysis summary exported")

# --- Create Intensity-Weighted Averages Summary ---
message("Creating intensity-weighted averages summary...")

# Load IWA data for summary creation
IWA_AuNP <- read_csv("output/02_comparison_analysis/IWA_AuNP.csv", show_col_types = FALSE)
IWA_P25 <- read_csv("output/02_comparison_analysis/IWA_P25.csv", show_col_types = FALSE)

# Combine IWA data for comparison
IWA_combined <- bind_rows(
  IWA_AuNP %>% mutate(NP_type = "AuNP"),
  IWA_P25 %>% mutate(NP_type = "P25")
)

# Create IWA comparison summary
IWA_summary <- IWA_combined %>%
  group_by(NP_type, new_name) %>%
  summarise(
    n_measurements = n(),
    mean_hc = mean(weighted_formula_hc, na.rm = TRUE),
    mean_oc = mean(weighted_formula_oc, na.rm = TRUE),
    mean_dbe = mean(weighted_formula_dbe, na.rm = TRUE),
    mean_ai = mean(weighted_formula_ai, na.rm = TRUE),
    total_intensity = mean(total_intensity, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(NP_type, new_name)

write_csv(IWA_summary, "output/04_export_results/data/IWA_summary.csv")
message("Intensity-weighted averages summary exported")

# --- Create Formula Class Analysis ---
message("Creating formula class analysis...")

# Analyze formula classes if available
if ("formula_class" %in% names(formulas_NP)) {
  formula_class_analysis <- formulas_NP %>%
    group_by(NP, new_name, formula_class) %>%
    summarise(
      n_formulas = n(),
      mean_intensity = mean(rel_intensity_permille, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = formula_class,
      values_from = n_formulas,
      values_fill = 0
    ) %>%
    arrange(NP, new_name)
  
  write_csv(formula_class_analysis, "output/04_export_results/data/formula_class_analysis.csv")
  message("Formula class analysis exported")
}

# --- Create Analysis Documentation ---
message("Creating analysis documentation...")

# Create a comprehensive text report
report_content <- paste0(
  "NANOPARTICLE SRFA RATIO ANALYSIS REPORT\n",
  "=====================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Total Molecular Formulas: ", total_formulas, "\n",
  "Total Measurements: ", total_measurements, "\n",
  "Nanoparticle Types: ", n_nanoparticle_types, " (AuNP, P25)\n",
  "Ratio Combinations: ", n_ratio_combinations, "\n\n",
  "ANALYSIS PIPELINE:\n",
  "01_data_preparation.R - Data preparation, cleaning, and replicate analysis\n",
  "02_comparison_analysis.R - Nanoparticle comparison and SRFA analysis\n",
  "03_visualization.R - Comprehensive visualizations including:\n",
  "  - Van Krevelen diagrams (SRFA comparison, ΔRI, unique formulas)\n",
  "  - Van Krevelen variants (clean_bp_RI, avg_MFs_bp, avg_MFs_RI, rep_MFs)\n",
  "  - Intensity-weighted averages comparison (AuNP vs P25)\n",
  "  - Stacked plots (common vs unique, by measurement)\n",
  "  - DOC analysis plots\n",
  "  - Reproducibility analysis plots\n",
  "04_export_results.R - Final data export and summary reports\n\n",
  "OUTPUT FILES:\n",
  "- formulas_processed.csv: Complete processed dataset\n",
  "- combined_nanoparticle_data.csv: Combined AuNP and P25 data\n",
  "- IWA_AuNP.csv: Intensity-weighted averages for AuNP\n",
  "- IWA_P25.csv: Intensity-weighted averages for P25\n",
  "- summary_common.csv: Common formula summary by ratio group\n",
  "- summary_SRFA.csv: SRFA comparison by ratio group\n",
  "- summary_SRFA_measurement.csv: SRFA comparison by individual measurement\n",
  "- nanoparticle_comparison_summary.csv: Detailed nanoparticle comparison\n",
  "- srfa_analysis_summary.csv: SRFA coverage analysis\n",
  "- IWA_summary.csv: Intensity-weighted averages summary\n",
  "- reproducibility_thresholds.csv: Reproducibility thresholds by group\n",
  "- grouped_formulas.csv: Aggregated formulas across replicates\n\n",
  "KEY FINDINGS:\n"
)

# Add key findings to report
if (nrow(comparison_summary) > 0) {
  total_common <- sum(comparison_summary$common_formulas, na.rm = TRUE)
  total_unique <- sum(comparison_summary$unique_formulas, na.rm = TRUE)
  total_srfa <- sum(comparison_summary$srfa_matches, na.rm = TRUE)
  
  report_content <- paste0(
    report_content,
    "- Total common formulas between AuNP and P25: ", total_common, "\n",
    "- Total unique formulas: ", total_unique, "\n",
    "- Total SRFA matches: ", total_srfa, "\n",
    "- Average SRFA coverage: ", round(mean(srfa_analysis$srfa_coverage_pct, na.rm = TRUE), 1), "%\n"
  )
}

# Write report to file
writeLines(report_content, "output/04_export_results/reports/analysis_report.txt")
message("Analysis report saved")

# --- Create Data Dictionary ---
message("Creating data dictionary...")

# Define column descriptions
data_dictionary <- data.frame(
  Column_Name = c(
    "measurement_name", "formula_string", "formula_mass", "formula_hc", "formula_oc", 
    "formula_nc", "formula_sc", "formula_dbe", "formula_dbe_o", "formula_ai", 
    "formula_ai_mod", "formula_ai_con", "formula_xc", "formula_nosc", "formula_f_star",
    "peak_kmd_ch2", "peak_kmd_co2", "peak_kmd_h2", "peak_kmd_nh2", "peak_kmd_div_zstar", 
    "peak_zstar", "formula_dbe_div_c", "formula_dbe_div_h", "formula_dbe_div_o", 
    "C", "H", "N", "O", "S", "rel_intensity_permille", "group", "new_name", "NP", 
    "common", "is_SRFA", "occurrence_type", "RIdiff", "Prec_RIdiff"
  ),
  Description = c(
    "Sample measurement identifier",
    "Molecular formula string identifier",
    "Molecular formula mass (Da)",
    "Hydrogen to carbon ratio",
    "Oxygen to carbon ratio", 
    "Nitrogen to carbon ratio",
    "Sulfur to carbon ratio",
    "Double bond equivalent",
    "Double bond equivalent minus oxygen",
    "Aromaticity index",
    "Modified aromaticity index",
    "Condensed aromaticity index",
    "Carbon oxidation state",
    "Nominal oxidation state of carbon",
    "Modified aromaticity index (f*)",
    "Kendrick mass defect for CH2",
    "Kendrick mass defect for CO2",
    "Kendrick mass defect for H2",
    "Kendrick mass defect for NH2",
    "Kendrick mass defect divided by z*",
    "Charge state (z*)",
    "Double bond equivalent to carbon ratio",
    "Double bond equivalent to hydrogen ratio",
    "Double bond equivalent to oxygen ratio",
    "Number of carbon atoms",
    "Number of hydrogen atoms",
    "Number of nitrogen atoms",
    "Number of oxygen atoms",
    "Number of sulfur atoms",
    "Relative intensity normalized to total sum (‰)",
    "Sample group identifier",
    "Nanoparticle ratio identifier",
    "Nanoparticle type (AuNP or P25)",
    "Commonality classification (Common, AuNP, P25)",
    "Whether formula matches SRFA standard",
    "Reproducibility classification (Common, Unique)",
    "Relative intensity difference between replicates",
    "Percentage relative intensity difference"
  ),
  Data_Type = c(
    "Character", "Character", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric",
    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric",
    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric",
    "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric", "Numeric",
    "Numeric", "Character", "Character", "Character", "Character", "Logical", "Character", "Numeric", "Numeric"
  ),
  stringsAsFactors = FALSE
)

write_csv(data_dictionary, "output/04_export_results/data/data_dictionary.csv")
message("Data dictionary saved")

# --- Create Final Summary Statistics ---
message("Creating final summary statistics...")

final_summary <- list(
  analysis_date = Sys.Date(),
  total_formulas = total_formulas,
  total_measurements = total_measurements,
  n_nanoparticle_types = n_nanoparticle_types,
  n_ratio_combinations = n_ratio_combinations,
  mean_reproducibility = round(mean(summary_statistics$Common_pct, na.rm = TRUE), 1)
)

# Convert to dataframe and save
final_summary_df <- data.frame(
  Parameter = names(final_summary),
  Value = unlist(final_summary)
)

write_csv(final_summary_df, "output/04_export_results/final_summary_statistics.csv")
message("Final summary statistics saved")

# --- Copy Visualizations Tree ---
message("Copying visualization outputs...")

# Ensure plots root exists
if (!dir.exists("output/04_export_results/plots")) {
  dir.create("output/04_export_results/plots", recursive = TRUE)
}

# Directories to copy recursively (if present)
viz_dirs <- c(
  "output/03_visualization/comparison_plots",
  "output/03_visualization/stacked_plots", 
  "output/03_visualization/reproducibility",
  "output/03_visualization/van_krevelen"
)

# Copy main visualization directories
for (src_dir in viz_dirs) {
  if (dir.exists(src_dir)) {
    dest_dir <- file.path("output/04_export_results/plots", basename(src_dir))
    if (dir.exists(dest_dir)) {
      unlink(dest_dir, recursive = TRUE, force = TRUE)
    }
    ok <- file.copy(from = src_dir, to = "output/04_export_results/plots", recursive = TRUE)
    if (isTRUE(ok)) {
      message(paste("Copied directory:", src_dir, "->", dest_dir))
    } else {
      warning(paste("Failed to copy directory:", src_dir))
    }
  } else {
    message(paste("Source directory not found (skipping):", src_dir))
  }
}

# Check for specific Van Krevelen subdirectories and report
vk_subdirs <- c(
  "SRFA_comparison", "delta_RI_comparison", "unique_molecular_formulas",
  "clean_bp_RI", "avg_MFs_bp", "avg_MFs_RI", "rep_MFs"
)

message("Checking Van Krevelen subdirectories...")
for (subdir in vk_subdirs) {
  src_path <- file.path("output/03_visualization/van_krevelen", subdir)
  if (dir.exists(src_path)) {
    message(paste("  Found:", subdir))
  } else {
    message(paste("  Missing:", subdir))
  }
}

# --- Create README File ---
message("Creating README file...")

readme_content <- paste0(
  "# Nanoparticle SRFA Ratio Analysis Results\n\n",
  "This directory contains the complete results from the nanoparticle SRFA ratio analysis.\n\n",
  "## Directory Structure\n\n",
  "- `data/`: All analysis datasets and results\n",
  "- `plots/`: Visualizations and figures copied from visualization output\n",
  "- `reports/`: Analysis reports and documentation\n\n",
  "## Key Files\n\n",
  "### Core Datasets\n",
  "- `formulas_processed.csv`: Complete processed dataset with all molecular formulas\n",
  "- `combined_nanoparticle_data.csv`: Combined AuNP and P25 data with comparisons\n",
  "- `IWA_AuNP.csv`: Intensity-weighted averages for AuNP samples\n",
  "- `IWA_P25.csv`: Intensity-weighted averages for P25 samples\n\n",
  "### Summary Files\n",
  "- `analysis_summary.csv`: Basic analysis parameters\n",
  "- `nanoparticle_comparison_summary.csv`: Detailed nanoparticle comparison\n",
  "- `srfa_analysis_summary.csv`: SRFA coverage analysis\n",
  "- `IWA_summary.csv`: Intensity-weighted averages summary\n",
  "- `data_dictionary.csv`: Column descriptions and data types\n\n",
  "### Comparison Analysis (Data)\n",
  "- `summary_common.csv`: Common formula summary by ratio group\n",
  "- `summary_SRFA.csv`: SRFA comparison by ratio group\n",
  "- `summary_SRFA_measurement.csv`: SRFA comparison by individual measurement\n",
  "- `STD_comparison_summary.csv`: SRFA standard coverage statistics\n\n",
  "### Reproducibility Data\n",
  "- `reproducibility_thresholds.csv`: Reproducibility thresholds by group\n",
  "- `grouped_formulas.csv`: Aggregated formulas across replicates\n\n",
  "### Visualization Outputs\n",
  "- `plots/comparison_plots/intensity_weighted_averages/`: IWA comparison (AuNP vs P25)\n",
  "- `plots/comparison_plots/DOC_analysis/`: DOC analysis plots for each sample type\n",
  "- `plots/stacked_plots/common_vs_unique_formulas/`: Common vs unique stacked plots\n",
  "- `plots/stacked_plots/measurement_comparison/`: Stacked plots by measurement\n",
  "- `plots/reproducibility/`: Reproducibility violin plots (log10 RIdiff, %RIdiff)\n",
  "- `plots/van_krevelen/SRFA_comparison/`: Van Krevelen diagrams comparing with SRFA standard\n",
  "- `plots/van_krevelen/delta_RI_comparison/`: Van Krevelen diagrams showing ΔRI changes\n",
  "- `plots/van_krevelen/unique_molecular_formulas/`: Van Krevelen diagrams for unique formulas\n",
  "- `plots/van_krevelen/clean_bp_RI/`: Van Krevelen diagrams colored by RI(bp)\n",
  "- `plots/van_krevelen/avg_MFs_bp/`: Van Krevelen diagrams for averaged formulas (RI bp)\n",
  "- `plots/van_krevelen/avg_MFs_RI/`: Van Krevelen diagrams for averaged formulas (RI ‰)\n",
  "- `plots/van_krevelen/rep_MFs/`: Van Krevelen diagrams with reproducibility thresholds\n\n",
  "## Analysis Pipeline\n\n",
  "1. **Data Preparation** (`01_data_preparation.R`): Load, clean, and process molecular formula data\n",
  "2. **Comparison Analysis** (`02_comparison_analysis.R`): Compare AuNP and P25 nanoparticles\n",
  "3. **Visualization** (`03_visualization.R`): Create Van Krevelen diagrams and comparison plots\n",
  "4. **Export Results** (`04_export_results.R`): Compile final results and documentation\n\n",
  "## Usage\n\n",
  "The exported datasets can be used for:\n",
  "- Nanoparticle comparison studies\n",
  "- SRFA standard analysis\n",
  "- Intensity-weighted average calculations\n",
  "- Reproducibility assessments\n",
  "- Publication figures and tables\n\n",
  "For questions about the analysis, refer to the individual script files or the analysis report.\n"
)

writeLines(readme_content, "output/04_export_results/README.md")
message("README file created")

message("\nFinal export complete!")
message("All results saved to output/04_export_results/")
message("\nKey outputs:")
message("- data/: Complete datasets and summaries")
message("- plots/: Comprehensive visualizations including:")
message("  - Van Krevelen diagrams (SRFA comparison, ΔRI, unique formulas)")
message("  - Van Krevelen variants (clean_bp_RI, avg_MFs_bp, avg_MFs_RI, rep_MFs)")
message("  - Intensity-weighted averages comparison (AuNP vs P25)")
message("  - Stacked plots (common vs unique, by measurement)")
message("  - DOC analysis plots")
message("  - Reproducibility analysis plots")
message("- reports/: Analysis documentation")
message("- README.md: Usage guide and file descriptions")
message("\nNanoparticle SRFA ratio analysis completed successfully!")
