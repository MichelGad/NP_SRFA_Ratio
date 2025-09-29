###############################################################################
# Script: 03_visualization.R
# Purpose: Unified, consistent visualizations for nanoparticle SRFA ratio analysis
# Author: Michel Gad
# Date: 2025-09-19
# Description: 
#   - Robust helpers and data loading
#   - Van Krevelen diagrams (SRFA comparison, ΔRI, unique, clean_bp_RI,
#     avg_MFs by RI(bp) and RI(‰), reproducibility thresholds)
#   - Intensity-weighted averages comparison (AuNP vs P25)
#   - Stacked plots (common vs unique; by measurement)
#   - DOC analysis plots
#   - Reproducibility analysis plots
###############################################################################

cat("=============================================================================\n")
cat("Script: 03_visualization.R\n")
cat("Purpose: Unified visualizations for nanoparticle SRFA ratio analysis\n")
cat("Author: Michel Gad\n")
cat("Date: 2025-09-19\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
message("Loading required libraries...")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(scales)
library(ggpattern)
library(ggrepel)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------
output_base <- "output/03_visualization"

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

ensure_dir("output")
ensure_dir(output_base)
ensure_dir(file.path(output_base, "van_krevelen"))

message("Setting up plot saving function...")
save_plot <- function(plot_object, 
                      filename, 
                      width = 7, 
                      height = 7, 
                      dpi = 300,
                      device = cairo_pdf,
                      plot_exports = output_base) {
  if (!inherits(plot_object, "ggplot")) stop("plot_object must be a ggplot object")
  filepath <- file.path(plot_exports, filename)
  tryCatch({
    pdf(filepath, width = width, height = height)
    print(plot_object)
    dev.off()
    message("Successfully saved plot to: ", filepath)
  }, error = function(e) {
    dev.off()
    stop("Failed to save plot: ", e$message)
  })
}

# Helper function for consistent Van Krevelen plot styling
create_vk_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 0.5, 1, 1, "cm"),
      plot.caption = element_text(size = 10, color = "gray40", hjust = 0)
    )
}

# Helper function for consistent section headers
create_section_header <- function(section_name, description = "") {
  cat("\n")
  cat("=============================================================================\n")
  cat("Section:", section_name, "\n")
  if (description != "") cat("Description:", description, "\n")
  cat("=============================================================================\n")
  message("Creating", section_name, "...")
}

# Helper function for optional data validation
validate_optional_data <- function(data, required_cols, data_name) {
  if (is.null(data)) {
    message("  ", data_name, " not available - skipping related plots")
    return(FALSE)
  }
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    message("  Missing required columns in", data_name, ":", paste(missing_cols, collapse = ", "))
    message("  Skipping related plots")
    return(FALSE)
  }
  return(TRUE)
}

# Try multiple candidate paths and read_csv/read_csv2
safe_read_any <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) {
      message("  Reading: ", p)
      dat <- tryCatch(
        readr::read_csv(p, show_col_types = FALSE),
        error = function(e1) {
          message("    read_csv failed for ", p, " - trying read_csv2()...")
          tryCatch(
            readr::read_csv2(p, show_col_types = FALSE),
            error = function(e2) {
              message("    read_csv2 also failed for ", p, " - skipping file.")
              NULL
            }
          )
        }
      )
      if (!is.null(dat)) return(as_tibble(dat))
    }
  }
  return(NULL)
}

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
message("Loading data for visualization...")

# Mandatory inputs
formulas_NP <- read_csv("output/02_comparison_analysis/combined_nanoparticle_data.csv", show_col_types = FALSE)
STD_SRFA <- read_csv("output/01_data_preparation/formulas_processed.csv", show_col_types = FALSE) %>%
  filter(str_detect(measurement_name, "STD"))
IWA_AuNP <- read_csv("output/02_comparison_analysis/IWA_AuNP.csv", show_col_types = FALSE)
IWA_P25 <- read_csv("output/02_comparison_analysis/IWA_P25.csv", show_col_types = FALSE)
summary_SRFA <- read_csv("output/02_comparison_analysis/summary_SRFA.csv", show_col_types = FALSE)
summary_SRFA_measurement <- read_csv("output/02_comparison_analysis/summary_SRFA_measurement.csv", show_col_types = FALSE)

# Optional inputs (multiple candidate paths)
paths_formulas_avg <- c(
  "output/02_comparison_analysis/formulas_avg.csv",
  "output/01_data_preparation/formulas_avg.csv",
  "processed_output/formulas_avg.csv"
)
paths_formulas_processed <- c(
  "output/01_data_preparation/formulas_processed.csv",
  "output/02_comparison_analysis/formulas_processed.csv",
  "processed_output/formulas_processed.csv"
)
paths_formulas_clean <- c(
  "output/01_data_preparation/formulas_clean.csv",
  "output/02_comparison_analysis/formulas_clean.csv",
  "processed_output/formulas_clean.csv"
)
paths_summary_table <- c(
  "output/01_data_preparation/summary_statistics.csv",
  "output/02_comparison_analysis/summary_statistics.csv",
  "processed_output/summary_statistics.csv"
)

formulas_avg <- safe_read_any(paths_formulas_avg)
formulas_processed_all <- safe_read_any(paths_formulas_processed)
formulas_clean <- safe_read_any(paths_formulas_clean)
summary_table <- safe_read_any(paths_summary_table)

message("Data loaded successfully for visualization")

# -----------------------------------------------------------------------------
# Section: Van Krevelen — SRFA comparison
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen SRFA Comparison", "Comparison with SRFA Standard")

vk_srfa_dir <- "van_krevelen/SRFA_comparison"
ensure_dir(file.path(output_base, vk_srfa_dir))

required_cols <- c("measurement_name", "formula_oc", "formula_hc", 
                   "rel_intensity_permille", "is_SRFA", "NP")
if (!all(required_cols %in% names(formulas_NP))) {
  stop("Required columns missing in formulas_NP: ",
       paste(setdiff(required_cols, names(formulas_NP)), collapse = ", "))
}

unique_measurements <- unique(formulas_NP$measurement_name)
for (measurement in unique_measurements) {
  plot_data <- formulas_NP %>%
    filter(measurement_name == measurement) %>%
    mutate(category = ifelse(is_SRFA, "SRFA", "NP")) %>%
    arrange(rel_intensity_permille)
  if (nrow(plot_data) == 0) next

  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    geom_point(data = filter(plot_data, category == "SRFA"), aes(color = "SRFA"), size = 2.5, alpha = 0.7) +
    geom_point(data = filter(plot_data, category == "NP"), aes(color = "NP"), size = 3, alpha = 0.5) +
    scale_color_manual(name = "Formula", values = c("SRFA" = "#E41A1C", "NP" = "#377EB8")) +
    labs(title = paste("Van Krevelen Diagram:", measurement),
         subtitle = "Comparison with SRFA Standard",
         x = "O/C Ratio", y = "H/C Ratio",
         caption = "Unique NP formulas (blue) are those not found in SRFA standard") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    create_vk_theme()
  
  filename <- paste0("vK_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(vk_plot, file.path(vk_srfa_dir, filename))
}
message(paste("Completed generating", length(unique_measurements), 
              "SRFA comparison diagrams in:", file.path(output_base, vk_srfa_dir)))

# -----------------------------------------------------------------------------
# Section: Van Krevelen — ΔRI comparison (common vs unique)
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen ΔRI Comparison", "Common vs unique formulas with intensity changes")

vk_deltari_dir <- "van_krevelen/delta_RI_comparison"
ensure_dir(file.path(output_base, vk_deltari_dir))

for (measurement in unique_measurements) {
  plot_data <- formulas_NP %>%
    filter(measurement_name == measurement) %>%
    mutate(
      δRI = ifelse(
        is_SRFA,
        (rel_intensity_permille - STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)]) / 
          STD_SRFA$rel_intensity_permille[match(formula_string, STD_SRFA$formula_string)] * 100,
        NA_real_
      )
    ) %>%
    arrange(rel_intensity_permille)
  if (nrow(plot_data) == 0) next

  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    geom_point(data = filter(plot_data, is.na(δRI)), aes(shape = "Unique"), color = "#E41A1C", size = 2.5, alpha = 0.5) +
    geom_point(data = filter(plot_data, !is.na(δRI)), aes(color = δRI), size = 3, alpha = 0.7) +
    scale_color_viridis_c(option = "D", name = expression(paste(Delta, "RI(%)")),
                          limits = quantile(plot_data$δRI, c(0.05, 0.95), na.rm = TRUE), oob = scales::squish) +
    scale_shape_manual(values = c("Unique" = 16), name = NULL) +
    labs(title = paste("Van Krevelen Diagram:", measurement), 
         x = "O/C Ratio", y = "H/C Ratio",
         caption = expression(paste("Unique (red) are NP-specific | Color shows ", Delta, "RI(%) for common formulas"))) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    create_vk_theme() +
    theme(legend.box.spacing = unit(0.2, "cm"))
  
  filename <- paste0("vK_deltaRI_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(vk_plot, file.path(vk_deltari_dir, filename))
}
message(paste("Completed generating", length(unique_measurements), 
              "δRI comparison diagrams in:", file.path(output_base, vk_deltari_dir)))

# -----------------------------------------------------------------------------
# Section: Van Krevelen — unique NP molecular formulas (log10 RI%)
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen Unique Formulas", "Unique NP molecular formulas colored by log10 intensity")

vk_unique_dir <- "van_krevelen/unique_molecular_formulas"
ensure_dir(file.path(output_base, vk_unique_dir))

for (measurement in unique_measurements) {
  plot_data <- formulas_NP %>%
    filter(measurement_name == measurement, !is_SRFA) %>%
    mutate(log_intensity = log10(rel_intensity_permille)) %>%
    arrange(rel_intensity_permille)
  if (nrow(plot_data) == 0) next

  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    geom_point(aes(color = log_intensity), size = 3.5, alpha = 0.8) +
    scale_color_viridis_c(option = "D", name = expression(log[10] * "(RI%)"),
                          guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"),
                                                 title.position = "top", title.hjust = 0.5)) +
    labs(title = paste("Unique Molecular Formulas:", measurement),
         subtitle = expression(paste("Colored by ", log[10], " Relative Intensity")),
         x = "O/C Ratio", y = "H/C Ratio",
         caption = "Showing only formulas unique to nanoparticle samples (not found in SRFA standard)") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    create_vk_theme()
  
  filename <- paste0("vK_unique_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(vk_plot, file.path(vk_unique_dir, filename))
}
message(paste("Completed generating", length(unique_measurements), 
              "unique molecular formula diagrams in:", file.path(output_base, vk_unique_dir)))

# -----------------------------------------------------------------------------
# Section: Intensity-weighted averages (AuNP vs P25)
# -----------------------------------------------------------------------------
create_section_header("Intensity-Weighted Averages", "AuNP vs P25 comparison of H/C and O/C ratios")

iwa_dir <- "comparison_plots/intensity_weighted_averages"
ensure_dir(file.path(output_base, iwa_dir))

required_cols_iwa <- c("new_name", "weighted_formula_hc", "weighted_formula_oc")
if (!all(required_cols_iwa %in% names(IWA_AuNP)) || !all(required_cols_iwa %in% names(IWA_P25))) {
  stop("Required columns missing in IWA data frames")
}

selected_ratios <- c("NP-HH", "NP-HM", "NP-HL")
IWA_AuNP_filtered <- IWA_AuNP %>% filter(new_name %in% selected_ratios) %>% mutate(NP_type = "AuNP")
IWA_P25_filtered <- IWA_P25 %>% filter(new_name %in% selected_ratios) %>% mutate(NP_type = "P25")
combined_data <- bind_rows(IWA_AuNP_filtered, IWA_P25_filtered)
if (nrow(combined_data) == 0) stop("No data available for selected ratios")

hc_range <- range(combined_data$weighted_formula_hc, na.rm = TRUE)
oc_range <- range(combined_data$weighted_formula_oc, na.rm = TRUE)
hc_min <- hc_range[1] - 0.05 * diff(hc_range)
hc_max <- hc_range[2] + 0.05 * diff(hc_range)
oc_min <- oc_range[1] - 0.05 * diff(oc_range)
oc_max <- oc_range[2] + 0.05 * diff(oc_range)
  
iwa_plot <- ggplot(combined_data, aes(x = weighted_formula_hc, y = weighted_formula_oc, color = NP_type)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(aes(label = new_name), size = 4, box.padding = 0.5, point.padding = 0.3,
                  max.overlaps = 20, min.segment.length = 0.2, seed = 123) +
  scale_color_manual(name = "Nanoparticle", values = c("AuNP" = "#E41A1C", "P25" = "#377EB8")) +
  labs(title = "Intensity-Weighted Averages: AuNP vs P25", 
       subtitle = "Comparison of H/C and O/C Ratios",
       x = "H/C Ratio", y = "O/C Ratio") +
  coord_cartesian(xlim = c(hc_min, hc_max), ylim = c(oc_min, oc_max)) +
  create_vk_theme() +
  theme(plot.margin = margin(1, 3, 1, 1, "cm"))

filename <- "IWA_comparison_AuNP_vs_P25.pdf"
save_plot(iwa_plot, file.path(iwa_dir, filename), width = 10, height = 7)
message(paste("Successfully saved intensity-weighted averages plot to:", file.path(output_base, iwa_dir, filename)))

# -----------------------------------------------------------------------------
# Section: Stacked plots — common vs unique (by group)
# -----------------------------------------------------------------------------
create_section_header("Stacked Plots - Common vs Unique", "Formula changes relative to SRFA standard by group")

stacked_common_dir <- "stacked_plots/common_vs_unique_formulas"
ensure_dir(file.path(output_base, stacked_common_dir))

custom_x_order_groups <- c("NP-HH", "NP-MM", "NP-LL", "NP-HM", "NP-HL", "NP-MH", "NP-LH")

plot_data_groups <- summary_SRFA %>%
  select(new_name, Increased_pct, Decreased_pct, Non_SRFA_pct) %>%
  pivot_longer(cols = -new_name, names_to = "category", values_to = "percentage") %>%
  mutate(
    fill_group = case_when(
      category %in% c("Increased_pct", "Decreased_pct") ~ "common",
      category == "Non_SRFA_pct" ~ "unique"
    ),
    pattern_type = case_when(
      category == "Increased_pct" ~ "increased",
      category == "Decreased_pct" ~ "decreased",
      TRUE ~ "none"
    ),
    fill_group = factor(fill_group, levels = c("common", "unique")),
    new_name = factor(new_name, levels = custom_x_order_groups)
  )

fill_colors <- c("common" = "#3E8E41", "unique" = "#0072B2")
pattern_types <- c("increased" = "stripe", "decreased" = "circle", "none" = "none")

stacked_plot_groups <- ggplot(plot_data_groups, aes(x = new_name, y = percentage, fill = fill_group, pattern = pattern_type)) +
  geom_col_pattern(position = position_stack(reverse = FALSE), color = "black",
                   pattern_density = 0.05, pattern_spacing = 0.02, pattern_fill = "gray30") +
  geom_text(data = summary_SRFA %>% distinct(new_name, Total), aes(x = new_name, y = -1, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 3.5, color = "black", vjust = 1.5) +
  scale_fill_manual(values = fill_colors, labels = c("common" = "Common with SRFA", "unique" = "Unique to Nanoparticles"), name = "Formula Type") +
  scale_pattern_manual(values = pattern_types, labels = c("increased" = "Increased Intensity", "decreased" = "Decreased Intensity", "none" = "No Change"), name = "Intensity Change", na.translate = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = "Percentage of Molecular Formulas",
    title = "Molecular Formula Changes Relative to SRFA Standard",
    subtitle = "Showing common formulas with intensity changes and Unique",
       caption = "Percentage of common formulas shown at the top of the bar and unique ones at the bottom | Total formulas shown below each bar") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.1), margin = margin(t = 5)),
        axis.title.y = element_text(size = rel(1.2), face = "bold"),
        panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = "right", legend.title = element_text(size = rel(1), face = "bold"), legend.text = element_text(size = rel(0.9)),
        plot.title = element_text(face = "bold", size = rel(1.3)), plot.subtitle = element_text(color = "gray40", size = rel(1.1)), plot.caption = element_text(color = "gray50", size = rel(0.8))) +
  guides(fill = guide_legend(order = 1, override.aes = list(pattern = "none")),
         pattern = guide_legend(order = 2, override.aes = list(fill = "white"))) +
  coord_cartesian(clip = "off")

filename <- "formula_changes_stacked_patterns.pdf"
save_plot(stacked_plot_groups, file.path(stacked_common_dir, filename), width = 12, height = 7)
message(paste("- Plot saved to:", file.path(output_base, stacked_common_dir, filename)))

# -----------------------------------------------------------------------------
# Section: Stacked plots — by measurement
# -----------------------------------------------------------------------------
create_section_header("Stacked Plots - By Measurement", "Formula changes relative to SRFA standard by measurement")

stacked_meas_dir <- "stacked_plots/measurement_comparison"
ensure_dir(file.path(output_base, stacked_meas_dir))

custom_x_order_meas <- c("AuNP-HH", "P25-HH", "AuNP-MM", "P25-MM",
                   "AuNP-LL", "P25-LL", "AuNP-HM", "P25-HM",
                   "AuNP-HL", "P25-HL", "AuNP-MH", "P25-MH",
                   "AuNP-LH", "P25-LH")

plot_data_meas <- summary_SRFA_measurement %>%
  select(measurement_name, Increased_pct, Decreased_pct, Non_SRFA_pct) %>%
  pivot_longer(cols = -measurement_name, names_to = "category", values_to = "percentage") %>%
  mutate(
    fill_group = case_when(
      category %in% c("Increased_pct", "Decreased_pct") ~ "common",
      category == "Non_SRFA_pct" ~ "unique"
    ),
    pattern_type = case_when(
      category == "Increased_pct" ~ "increased",
      category == "Decreased_pct" ~ "decreased",
      TRUE ~ "none"
    ),
    fill_group = factor(fill_group, levels = c("common", "unique")),
    measurement_name = factor(measurement_name, levels = custom_x_order_meas)
  )

stacked_plot_meas <- ggplot(plot_data_meas, aes(x = measurement_name, y = percentage, fill = fill_group, pattern = pattern_type)) +
  geom_col_pattern(position = position_stack(reverse = FALSE), color = "black",
                   pattern_density = 0.05, pattern_spacing = 0.02, pattern_fill = "gray30") +
  geom_text(data = summary_SRFA_measurement %>% distinct(measurement_name, Total), aes(x = measurement_name, y = -1, label = paste0("n=", Total)),
            inherit.aes = FALSE, size = 3.5, color = "black", vjust = 1.5) +
  scale_fill_manual(values = fill_colors, labels = c("common" = "Common with SRFA", "unique" = "Unique to Nanoparticles"), name = "Formula Type") +
  scale_pattern_manual(values = pattern_types, labels = c("increased" = "Increased Intensity", "decreased" = "Decreased Intensity", "none" = "No Change"), name = "Intensity Change", na.translate = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = "Percentage of Molecular Formulas", title = "Molecular Formula Changes by Measurement",
    subtitle = "Showing common formulas with intensity changes and Unique",
       caption = "Percentage of common formulas shown at the top of the bar and unique ones at the bottom | Total formulas shown below each bar") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.1), margin = margin(t = 5)),
        axis.title.y = element_text(size = rel(1.2), face = "bold"),
        panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = "right", legend.title = element_text(size = rel(1), face = "bold"), legend.text = element_text(size = rel(0.9)),
        plot.title = element_text(face = "bold", size = rel(1.3)), plot.subtitle = element_text(color = "gray40", size = rel(1.1)), plot.caption = element_text(color = "gray50", size = rel(0.8))) +
  guides(fill = guide_legend(order = 1, override.aes = list(pattern = "none")),
         pattern = guide_legend(order = 2, override.aes = list(fill = "white"))) +
  coord_cartesian(clip = "off")

filename <- "formula_changes_by_measurement.pdf"
save_plot(stacked_plot_meas, file.path(stacked_meas_dir, filename), width = 14, height = 7)
message(paste("- Plot saved to:", file.path(output_base, stacked_meas_dir, filename)))


# -----------------------------------------------------------------------------
# Section: Reproducibility violins (log10 RIdiff and %RIdiff)
# -----------------------------------------------------------------------------
message("Creating reproducibility violin plots (log10 RIdiff and %RIdiff)...")
reprod_dir <- "reproducibility"
ensure_dir(file.path(output_base, reprod_dir))

if (!is.null(formulas_clean) && all(c("measurement_name", "formula_string") %in% names(formulas_clean))) {
  if (!"occurrence_type" %in% names(formulas_clean)) {
    message("  Deriving occurrence_type and RIdiff metrics from formulas_clean...")
    formulas_clean <- formulas_clean %>%
      group_by(group, formula_string) %>%
      mutate(formula_count = n(), occurrence_type = if_else(formula_count > 2, "Common", "Unique")) %>%
      ungroup()
  }

  if (!"RIdiff" %in% names(formulas_clean)) {
    if (!"peak_relint_bp" %in% names(formulas_clean)) {
      message("  Cannot compute RIdiff: 'peak_relint_bp' missing in formulas_clean")
    } else {
      formulas_clean <- formulas_clean %>%
        group_by(group, formula_string) %>%
        mutate(
          RIdiff = if_else(occurrence_type == "Common",
                           (max(peak_relint_bp, na.rm = TRUE) - min(peak_relint_bp, na.rm = TRUE)) /
                             pmax(min(peak_relint_bp, na.rm = TRUE), .Machine$double.eps),
                           NA_real_),
          Prec_RIdiff = if_else(occurrence_type == "Common", RIdiff * 100, NA_real_)
        ) %>%
        ungroup()
    }
  }

  common_prefixes <- c("AuNP-HH", "P25-HH", "AuNP-MM", "P25-MM", "AuNP-LL", "P25-LL",
                       "AuNP-HM", "P25-HM", "AuNP-HL", "P25-HL", "AuNP-MH", "P25-MH",
                       "AuNP-LH", "P25-LH")

  rep_data <- formulas_clean %>% filter(occurrence_type == "Common", group %in% common_prefixes)
  if (nrow(rep_data) > 0) {
    rep_data_log <- rep_data %>% filter(!is.na(RIdiff) & RIdiff > 0) %>% mutate(group = factor(group, levels = rev(common_prefixes)))
    if (nrow(rep_data_log) > 0) {
      medians_log <- rep_data_log %>% group_by(group) %>% summarise(median_log = median(log10(RIdiff), na.rm = TRUE), .groups = "drop")
      log_plot <- ggplot(rep_data_log, aes(x = log10(RIdiff), y = group)) +
        geom_violin(fill = "#9ecae1", alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.color = "#d62728") +
        geom_text(data = medians_log, aes(x = median_log, y = group, label = sprintf("%.2f", median_log)), color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +
        labs(title = bquote(log[10] * " RIdiffs Between Replicates"), x = expression(log[10]("RIdiff")), y = NULL, caption = "Median values for each group shown as text labels") +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none")
      save_plot(log_plot, file.path(reprod_dir, "reproducibility_logRIdiff.pdf"), width = 10, height = 6)
    } else {
      message("  No positive RIdiff values available for log10 plot — skipping log10 RIdiff violin.")
    }

    rep_data_pct <- rep_data %>% filter(!is.na(Prec_RIdiff)) %>% mutate(group = factor(group, levels = rev(common_prefixes)))
    if (nrow(rep_data_pct) > 0) {
      medians_pct <- rep_data_pct %>% group_by(group) %>% summarise(median_pct = median(Prec_RIdiff, na.rm = TRUE), .groups = "drop")
      pct_plot <- ggplot(rep_data_pct, aes(x = Prec_RIdiff, y = group)) +
        geom_violin(fill = "#a1d99b", alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.color = "#d62728") +
        geom_text(data = medians_pct, aes(x = median_pct, y = group, label = sprintf("%.1f%%", median_pct)), color = "black", size = 4, inherit.aes = FALSE) +
        labs(title = "%RIdiff Between Replicates", x = "%RIdiff", y = NULL, caption = "Median values shown as percentage labels") +
        scale_x_continuous(labels = scales::percent_format(scale = 1)) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
      save_plot(pct_plot, file.path(reprod_dir, "reproducibility_pctRIdiff.pdf"), width = 10, height = 6)
    } else {
      message("  No Prec_RIdiff values available — skipping %RIdiff violin.")
    }
  } else {
    message("  No 'Common' formulas found in the selected groups — skipping reproducibility plots.")
  }
} else {
  message("  formulas_clean not available or missing required columns — skipping reproducibility violin plots.")
}

# -----------------------------------------------------------------------------
# Section: Van Krevelen — clean_bp_RI
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen Clean bp RI", "Van Krevelen diagrams colored by RI(bp) from clean data")

vk_clean_bp_dir <- "van_krevelen/clean_bp_RI"
ensure_dir(file.path(output_base, vk_clean_bp_dir))

if (validate_optional_data(formulas_clean, c("measurement_name", "formula_oc", "formula_hc", "peak_relint_bp"), "formulas_clean")) {
  for (measurement in unique(formulas_clean$measurement_name)) {
    plot_data <- formulas_clean %>% filter(measurement_name == measurement) %>% arrange(peak_relint_bp)
    if (nrow(plot_data) == 0) next
    plot_limits <- c(0, 1)
    if (nrow(plot_data) > 1) {
      qlim <- quantile(plot_data$peak_relint_bp, probs = c(0.01, 0.99), na.rm = TRUE)
      plot_limits <- c(qlim[1], qlim[2])
    }
    vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc, color = peak_relint_bp)) +
      geom_point(alpha = 0.9, size = 3) +
      scale_color_viridis_c(name = "RI(bp)", option = "viridis", direction = 1, limits = plot_limits, oob = scales::squish,
                            guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"), title.position = "top", title.hjust = 0.5)) +
      labs(title = paste("Van Krevelen Diagram:", measurement), 
           x = "O/C Ratio", y = "H/C Ratio", 
           caption = "RI(bp) is the peak intensity normalized to base peak intensity") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
      create_vk_theme()
    
    filename <- paste0("vK_bp_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
    save_plot(vk_plot, file.path(vk_clean_bp_dir, filename))
  }
  message(paste("Completed generating clean bp RI diagrams in:", file.path(output_base, vk_clean_bp_dir)))
}

# -----------------------------------------------------------------------------
# Section: Van Krevelen — avg_MFs_bp
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen Averaged MFs bp", "Van Krevelen diagrams colored by RI(bp) from averaged data")

vk_avg_bp_dir <- "van_krevelen/avg_MFs_bp"
ensure_dir(file.path(output_base, vk_avg_bp_dir))

if (validate_optional_data(formulas_processed_all, c("measurement_name", "formula_oc", "formula_hc", "peak_relint_bp"), "formulas_processed_all")) {
  for (measurement in unique(formulas_processed_all$measurement_name)) {
    plot_data <- formulas_processed_all %>% filter(measurement_name == measurement) %>% arrange(peak_relint_bp)
    if (nrow(plot_data) == 0) next
    plot_limits <- c(0, 1)
    if (nrow(plot_data) > 1) {
      qlim <- quantile(plot_data$peak_relint_bp, probs = c(0.01, 0.99), na.rm = TRUE)
      plot_limits <- c(qlim[1], qlim[2])
    }
    vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc, color = peak_relint_bp)) +
      geom_point(alpha = 0.9, size = 3) +
      scale_color_viridis_c(name = "RI(bp)", option = "viridis", direction = 1, limits = plot_limits, oob = scales::squish,
                            guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"), title.position = "top", title.hjust = 0.5)) +
      labs(title = paste("Van Krevelen Diagram:", measurement), 
           x = "O/C Ratio", y = "H/C Ratio", 
           caption = "RI(bp) is the peak intensity normalized to base peak intensity") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
      create_vk_theme()
    
    filename <- paste0("vK_bp_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
    save_plot(vk_plot, file.path(vk_avg_bp_dir, filename))
  }
  message(paste("Completed generating averaged bp RI diagrams in:", file.path(output_base, vk_avg_bp_dir)))
}

# -----------------------------------------------------------------------------
# Section: Van Krevelen — avg_MFs_RI
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen Averaged MFs RI", "Van Krevelen diagrams colored by RI(‰) from averaged data")

vk_avg_ri_dir <- "van_krevelen/avg_MFs_RI"
ensure_dir(file.path(output_base, vk_avg_ri_dir))

if (validate_optional_data(formulas_processed_all, c("measurement_name", "formula_oc", "formula_hc", "rel_intensity_permille"), "formulas_processed_all")) {
  for (measurement in unique(formulas_processed_all$measurement_name)) {
    plot_data <- formulas_processed_all %>% filter(measurement_name == measurement) %>% arrange(rel_intensity_permille)
    if (nrow(plot_data) == 0) next
    plot_limits <- c(0, 1000)
    if (nrow(plot_data) > 1) {
      qlim <- quantile(plot_data$rel_intensity_permille, probs = c(0.01, 0.99), na.rm = TRUE)
      plot_limits <- c(qlim[1], qlim[2])
    }
    vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc, color = rel_intensity_permille)) +
      geom_point(alpha = 0.9, size = 3) +
      scale_color_viridis_c(name = "RI(‰)", option = "viridis", direction = 1, limits = plot_limits, oob = scales::squish,
                            guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"), title.position = "top", title.hjust = 0.5)) +
      labs(title = paste("Van Krevelen Diagram:", measurement), 
           x = "O/C Ratio", y = "H/C Ratio",
           caption = "RI(‰) is the peak intensities normalized to sum of all peaks present (‰)") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
      create_vk_theme()
    
    filename <- paste0("vK_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
    save_plot(vk_plot, file.path(vk_avg_ri_dir, filename))
  }
  message(paste("Completed generating averaged RI diagrams in:", file.path(output_base, vk_avg_ri_dir)))
}

# -----------------------------------------------------------------------------
# Section: Van Krevelen — rep_MFs
# -----------------------------------------------------------------------------
create_section_header("Van Krevelen Reproducibility MFs", "Van Krevelen diagrams with reproducibility thresholds")

vk_rep_dir <- "van_krevelen/rep_MFs"
ensure_dir(file.path(output_base, vk_rep_dir))

# Calculate reproducibility thresholds
reproducibility_threshold_window <- NULL
if (!is.null(formulas_clean) && all(c("group", "occurrence_type", "peak_relint_bp") %in% names(formulas_clean))) {
  reproducibility_threshold_window <- formulas_clean %>%
    filter(occurrence_type == "Common") %>%
    group_by(group) %>%
    summarise(threshold_min = quantile(peak_relint_bp, 0.025, na.rm = TRUE), 
              threshold_max = quantile(peak_relint_bp, 0.975, na.rm = TRUE), 
              .groups = "drop")
}

if (validate_optional_data(formulas_processed_all, c("measurement_name", "formula_oc", "formula_hc", "rel_intensity_permille", "group"), "formulas_processed_all") && 
    !is.null(reproducibility_threshold_window)) {
  for (measurement in unique(formulas_processed_all$measurement_name)) {
    plot_data <- formulas_processed_all %>% filter(measurement_name == measurement) %>% arrange(rel_intensity_permille)
    if (nrow(plot_data) == 0) next
    group_name <- unique(plot_data$group)
    if (length(group_name) != 1) next
    group_threshold <- reproducibility_threshold_window %>% filter(group == group_name)
    if (nrow(group_threshold) == 0) next

    minv <- min(plot_data$rel_intensity_permille, na.rm = TRUE)
    maxv <- max(plot_data$rel_intensity_permille, na.rm = TRUE)
    plot_data <- plot_data %>%
      mutate(color_value = ifelse(rel_intensity_permille >= group_threshold$threshold_min & 
                                 rel_intensity_permille <= group_threshold$threshold_max, 
                                 NA_real_, rel_intensity_permille))
    values_vec <- scales::rescale(c(minv, group_threshold$threshold_min, group_threshold$threshold_max, maxv), to = c(0, 1))

    vk_plot <- ggplot() +
      geom_point(data = dplyr::filter(plot_data, is.na(color_value)), 
                 aes(x = formula_oc, y = formula_hc), color = "grey", alpha = 0.9, size = 3) +
      geom_point(data = dplyr::filter(plot_data, !is.na(color_value)), 
                 aes(x = formula_oc, y = formula_hc, color = color_value), alpha = 0.9, size = 3) +
      scale_color_gradientn(colours = c("blue", "grey", "red"), values = values_vec, name = "RI(‰)", na.value = "grey",
                            limits = c(minv, maxv), oob = scales::squish,
                            guide = guide_colorbar(barwidth = unit(0.5, "cm"), barheight = unit(5, "cm"), title.position = "top", title.hjust = 0.5)) +
      labs(title = paste("Van Krevelen Diagram:", measurement), 
           x = "O/C Ratio", y = "H/C Ratio",
           caption = "Grey points within reproducibility threshold (‰) — colored points are outside the reproducibility window") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
      create_vk_theme()

    filename <- paste0("vK_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
    save_plot(vk_plot, file.path(vk_rep_dir, filename))
  }
  message(paste("Completed generating reproducibility threshold diagrams in:", file.path(output_base, vk_rep_dir)))
}

# -----------------------------------------------------------------------------
# Section: DOC Analysis Plots
# -----------------------------------------------------------------------------
create_section_header("DOC Analysis Plots", "Dissolved Organic Carbon analysis plots for each sample type")

doc_dir <- "comparison_plots/DOC_analysis"
ensure_dir(file.path(output_base, doc_dir))

# Read DOC data
doc_file <- "input/DOC.csv"
if (file.exists(doc_file)) {
  message("  Reading DOC data from:", doc_file)
  
  # Read DOC data with proper decimal handling (comma as decimal separator)
  doc_data <- read_csv2(doc_file, show_col_types = FALSE) %>%
    # The first column is unnamed, so we need to handle it properly
    rename(sample_type = `...1`, sample_name = `Samples`, doc_mg_L = `DOC (mg/L)`) %>%
    # Convert comma decimal to dot decimal for numeric conversion
    mutate(
      doc_mg_L = as.numeric(gsub(",", ".", doc_mg_L))
    ) %>%
    # Filter out any rows with missing data
    filter(!is.na(doc_mg_L))
  
  if (nrow(doc_data) > 0) {
    # Generate replicate data for SRFA samples (different numbers of replicates)
    if ("SRFA" %in% doc_data$sample_type) {
      srfa_data <- doc_data %>% filter(sample_type == "SRFA")
      
      # Create replicates with different numbers based on sample name
      srfa_replicates_list <- list()
      
      for (i in seq_len(nrow(srfa_data))) {
        sample_name <- srfa_data$sample_name[i]
        n_replicates <- case_when(
          str_detect(sample_name, "SRFA-M|SRFA-L") ~ 2,  # 2 replicates for SRFA-M and SRFA-L
          TRUE ~ 3  # 3 replicates for other samples (e.g., SRFA-H)
        )
        
        sample_replicates <- srfa_data[i, ] %>%
          slice(rep(1, each = n_replicates)) %>%
          mutate(
            replicate = 1:n_replicates,
            sample_name = paste0(sample_name, "_", replicate)
          )
        
        srfa_replicates_list[[i]] <- sample_replicates
      }
      
      srfa_replicates <- bind_rows(srfa_replicates_list)
      
      # Replace original SRFA data with replicates
      doc_data <- doc_data %>%
        filter(sample_type != "SRFA") %>%
        bind_rows(srfa_replicates)
    }
    
    # Create stacked bar plot for each sample type
    sample_types <- unique(doc_data$sample_type)
    
    for (sample_type in sample_types) {
      type_data <- doc_data %>% filter(sample_type == !!sample_type)
      
      if (nrow(type_data) > 0) {
        # Prepare data for plot (only DOC mg/L) with custom ordering
        plot_data <- type_data %>%
          select(sample_name, doc_mg_L) %>%
          rename(value = doc_mg_L) %>%
          # Fix sample names for display (convert subscript 2 to regular 2)
          mutate(display_name = str_replace(sample_name, "₂", "2")) %>%
          mutate(
            # Custom ordering based on sample type
            order_factor = case_when(
              sample_type == "SRFA" ~ case_when(
                str_detect(sample_name, "SRFA-H") ~ 1,
                str_detect(sample_name, "SRFA-M") ~ 2,
                str_detect(sample_name, "SRFA-L") ~ 3,
                TRUE ~ 4
              ),
              TRUE ~ case_when(
                str_extract(sample_name, "[HML]$") == "H" ~ 1,
                str_extract(sample_name, "[HML]$") == "M" ~ 2,
                str_extract(sample_name, "[HML]$") == "L" ~ 3,
                TRUE ~ 4
              )
            ),
            # For SRFA replicates, also order by replicate number
            replicate_order = ifelse(sample_type == "SRFA", 
                                   as.numeric(str_extract(sample_name, "[0-9]+$")), 
                                   0)
          ) %>%
          arrange(order_factor, replicate_order, sample_name) %>%
          mutate(display_name = factor(display_name, levels = display_name))
        
        # Define colors for each sample type
        sample_colors <- c("SRFA" = "#E41A1C", "TiO₂" = "#377EB8", "AuNP" = "#4DAF4A")
        line_color <- sample_colors[sample_type]
        
        # Create proper display name for title
        display_name <- case_when(
          sample_type == "TiO₂" ~ "TiO[2]",
          TRUE ~ sample_type
        )
        
        # Calculate positions for vertical lines to separate H, M, L regions
        # First line between SRFA-H_3 and SRFA-M_1, second line between SRFA-M_2 and SRFA-L_1
        n_points <- nrow(plot_data)
        if (n_points >= 6) {
          # Find the position after SRFA-H_3 (which should be at position 3)
          h_region_end <- 3.5  # Between SRFA-H_3 and SRFA-M_1
          # Find the position after SRFA-M_2 (which should be at position 5)
          m_region_end <- 5.5   # Between SRFA-M_2 and SRFA-L_1
        } else {
          h_region_end <- 1.5
          m_region_end <- 2.5
        }
        
        # Create the plot
        doc_plot <- ggplot(plot_data, aes(x = display_name, y = value, group = 1)) +
          # Add vertical lines to separate regions
          geom_vline(xintercept = h_region_end, linetype = "dashed", color = "gray50", alpha = 0.7, size = 0.8) +
          geom_vline(xintercept = m_region_end, linetype = "dashed", color = "gray50", alpha = 0.7, size = 0.8) +
          geom_line(size = 1.2, alpha = 0.8, color = line_color) +
          geom_point(size = 3, alpha = 0.9, color = line_color) +
          geom_text(aes(label = display_name), size = 3, vjust = -1.2, hjust = 0.5, color = "black", 
                    fontface = "bold", alpha = 0.8) +
          labs(
            title = paste("DOC Analysis:", display_name),
            subtitle = "Dissolved Organic Carbon Measurements",
            x = "Sample Name",
            y = "DOC Concentration (mg/L)",
            caption = "Line plot showing DOC (mg/L) measurements"
          ) +
    theme_bw(base_size = 14) +
          theme(
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold", size = 16),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
            plot.margin = margin(1, 1, 1, 1, "cm"),
            plot.caption = element_text(size = 10, color = "gray40", hjust = 0)
          ) +
          scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(-0.5, NA)) +
          coord_cartesian(ylim = c(-0.5, NA), clip = "off")
        
        # Save the plot
        filename <- paste0("DOC_analysis_", gsub("[^[:alnum:]₂]", "_", sample_type), ".pdf")
        save_plot(doc_plot, file.path(doc_dir, filename), width = 10, height = 7)
        message(paste("  - Saved DOC plot for", sample_type, "to:", file.path(doc_dir, filename)))
      }
    }
    
    # Create a combined comparison plot (only DOC mg/L) with custom ordering
    combined_plot_data <- doc_data %>%
      select(sample_name, doc_mg_L, sample_type) %>%
      rename(value = doc_mg_L) %>%
      # Fix sample names for display (convert subscript 2 to regular 2)
      mutate(display_name = str_replace(sample_name, "₂", "2")) %>%
  mutate(
        # Custom ordering based on sample type
        order_factor = case_when(
          sample_type == "SRFA" ~ case_when(
            str_detect(sample_name, "SRFA-H") ~ 1,
            str_detect(sample_name, "SRFA-M") ~ 2,
            str_detect(sample_name, "SRFA-L") ~ 3,
            TRUE ~ 4
          ),
          TRUE ~ case_when(
            str_extract(sample_name, "[HML]$") == "H" ~ 1,
            str_extract(sample_name, "[HML]$") == "M" ~ 2,
            str_extract(sample_name, "[HML]$") == "L" ~ 3,
            TRUE ~ 4
          )
        ),
        # For SRFA replicates, also order by replicate number
        replicate_order = ifelse(sample_type == "SRFA", 
                               as.numeric(str_extract(sample_name, "[0-9]+$")), 
                               0)
      ) %>%
      arrange(sample_type, order_factor, replicate_order, sample_name) %>%
      # Create a continuous x-axis position for each sample type
      group_by(sample_type) %>%
      mutate(x_position = row_number()) %>%
      ungroup()
    
    # Define colors for each sample type
    sample_colors <- c("SRFA" = "#E41A1C", "TiO₂" = "#377EB8", "AuNP" = "#4DAF4A")
    
    # Create proper labels for the legend with correct mathematical notation
    sample_labels <- c("SRFA" = "SRFA", "TiO₂" = "TiO[2]", "AuNP" = "AuNP")
    
    # Calculate positions for vertical lines to separate H, M, L regions in combined plot
    # First line between SRFA-H_3 and SRFA-M_1, second line between SRFA-M_2 and SRFA-L_1
    n_combined_points <- nrow(combined_plot_data)
    if (n_combined_points >= 6) {
      # Find the position after SRFA-H_3 (which should be at position 3)
      h_region_end_combined <- 3.5  # Between SRFA-H_3 and SRFA-M_1
      # Find the position after SRFA-M_2 (which should be at position 5)
      m_region_end_combined <- 5.5   # Between SRFA-M_2 and SRFA-L_1
    } else {
      h_region_end_combined <- 1.5
      m_region_end_combined <- 2.5
    }
    
    combined_plot <- ggplot(combined_plot_data, aes(x = x_position, y = value, color = sample_type, group = sample_type)) +
      # Add vertical lines to separate regions
      geom_vline(xintercept = h_region_end_combined, linetype = "dashed", color = "gray50", alpha = 0.7, size = 0.8) +
      geom_vline(xintercept = m_region_end_combined, linetype = "dashed", color = "gray50", alpha = 0.7, size = 0.8) +
      geom_line(size = 1.2, alpha = 0.8) +
      geom_point(size = 3, alpha = 0.9) +
      geom_text_repel(aes(label = display_name), size = 3.5, color = "black", 
                      box.padding = 1.2, point.padding = 0.8, 
                      max.overlaps = 100, min.segment.length = 0.2, 
                      force = 4, force_pull = 3, seed = 123,
                      direction = "both", nudge_x = 0.3, nudge_y = 0.3,
                      fontface = "bold", alpha = 0.9) +
      scale_color_manual(values = sample_colors, labels = sample_labels, name = "Sample Type") +
      labs(
        title = "DOC Comparison",
        subtitle = "Comparing DOC concentrations between SRFA solutions and supernatant of each NP mixture",
        x = NULL,
        y = "DOC Concentration (mg/L)"
        ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.caption = element_text(size = 10, color = "gray40", hjust = 0)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(-0.5, NA)) +
      coord_cartesian(ylim = c(-0.5, NA), clip = "off")
    
    # Save the combined plot
    filename <- "DOC_analysis_comparison.pdf"
    save_plot(combined_plot, file.path(doc_dir, filename), width = 12, height = 10)
    message(paste("  - Saved combined DOC comparison plot to:", file.path(doc_dir, filename)))
    
    message(paste("Completed generating DOC analysis plots in:", file.path(output_base, doc_dir)))
} else {
    message("  No valid DOC data found after processing")
  }
} else {
  message("  DOC file not found at:", doc_file)
}

# -----------------------------------------------------------------------------
# Footer
# -----------------------------------------------------------------------------
message("\nVisualization complete!")
message("Results saved to output/03_visualization/")
message("- van_krevelen/: Van Krevelen diagrams for various comparisons")
message("- comparison_plots/: Intensity-weighted average comparison plots")
message("- stacked_plots/: Stacked column plots for formula analysis")


