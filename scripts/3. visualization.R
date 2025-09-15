# =============================================================================
# Script: 3. visualization.R
# Purpose: Create comprehensive visualizations for nanoparticle SRFA ratio analysis
# Author: Michel Gad
# Date: 2025-09-15
# Description: 
#   - Generate Van Krevelen diagrams for various comparisons
#   - Create stacked column plots and intensity-weighted average plots
#   - Export publication-ready figures for nanoparticle analysis
# =============================================================================

# Print script header information
cat("=============================================================================\n")
cat("Script: 3. visualization.R\n")
cat("Purpose: Create comprehensive visualizations for nanoparticle SRFA ratio analysis\n")
cat("Author: Michel Gad\n")
cat("Date: 2025-09-15\n")
cat("=============================================================================\n\n")

# --- Load Required Libraries ---
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

# --- Set Working Directory and Create Output Folders ---
message("Setting up directories...")
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}
if (!dir.exists("output/visualization")) {
  dir.create("output/visualization", recursive = TRUE)
}
if (!dir.exists("output/visualization/van_krevelen")) {
  dir.create("output/visualization/van_krevelen", recursive = TRUE)
}

# --- Load Data ---
message("Loading data for visualization...")

# Load comparison analysis results
formulas_NP <- read_csv("output/comparison_analysis/combined_nanoparticle_data.csv", show_col_types = FALSE)
STD_SRFA <- read_csv("output/data_prep/formulas_processed.csv", show_col_types = FALSE) %>%
  filter(str_detect(measurement_name, "STD"))
IWA_AuNP <- read_csv("output/comparison_analysis/IWA_AuNP.csv", show_col_types = FALSE)
IWA_P25 <- read_csv("output/comparison_analysis/IWA_P25.csv", show_col_types = FALSE)
summary_SRFA <- read_csv("output/comparison_analysis/summary_SRFA.csv", show_col_types = FALSE)
summary_SRFA_measurement <- read_csv("output/comparison_analysis/summary_SRFA_measurement.csv", show_col_types = FALSE)

message("Data loaded successfully for visualization")

# --- Define Plot Saving Function ---
message("Setting up plot saving function...")

# Define a robust plot saving function with error handling
save_plot <- function(plot_object, 
                      filename, 
                      width = 7, 
                      height = 7, 
                      dpi = 300,
                      device = cairo_pdf,
                      plot_exports = "output/visualization") {
  
  # Validate inputs
  if (!inherits(plot_object, "ggplot")) {
    stop("plot_object must be a ggplot object")
  }
  
  # Construct full file path
  filepath <- file.path(plot_exports, filename)
  
  # Save as PDF with error handling
  tryCatch({
    pdf(filepath, width = width, height = height)
    print(plot_object)
    dev.off()
    message("Successfully saved plot to: ", filepath)
  }, error = function(e) {
    dev.off() # Ensure device is closed
    stop("Failed to save plot: ", e$message)
  })
}

# --- Van Krevelen Diagrams: SRFA Comparison ---
message("Creating Van Krevelen diagrams for SRFA comparison...")

# Create directory structure
output_base <- "output/visualization"
subdir <- "van_krevelen/SRFA_comparison"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Check if required data exists
required_cols <- c("measurement_name", "formula_oc", "formula_hc", 
                   "rel_intensity_permille", "is_SRFA", "NP")
if (!all(required_cols %in% names(formulas_NP))) {
  stop("Required columns missing in formulas_NP: ",
       paste(setdiff(required_cols, names(formulas_NP)), collapse = ", "))
}

# Generate plots for each measurement
unique_measurements <- unique(formulas_NP$measurement_name)

for (measurement in unique_measurements) {
  # Filter and prepare data
  plot_data <- formulas_NP %>%
    filter(measurement_name == measurement) %>%
    mutate(
      category = ifelse(is_SRFA, "SRFA Standard", "NP Formulas")
    ) %>%
    arrange(rel_intensity_permille)
  
  if (nrow(plot_data) == 0) {
    warning("No data available for measurement: ", measurement)
    next
  }
  
  # Calculate statistics
  n_srfa <- sum(plot_data$is_SRFA)
  n_np <- sum(!plot_data$is_SRFA)
  
  message(paste("Processing:", measurement))
  message(paste("  SRFA:", n_srfa, paste0("(", round(n_srfa/nrow(plot_data)*100, 1), "%)")))
  message(paste("  NP:", n_np, paste0("(", round(n_np/nrow(plot_data)*100, 1), "%)")))
  
  # Create optimized plot
  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    # SRFA points (red)
    geom_point(
      data = filter(plot_data, is_SRFA),
      aes(color = "SRFA"),
      size = 2.5,
      alpha = 0.7
    ) +
    # NP points (blue)
    geom_point(
      data = filter(plot_data, !is_SRFA),
      aes(color = "NP"),
      size = 3,
      alpha = 0.5
    ) +
    scale_color_manual(
      name = "Formula",
      values = c(
        "SRFA" = "#E41A1C",  # Red
        "NP" = "#377EB8"     # Blue
      )
    ) +
    labs(
      title = paste("Van Krevelen Diagram:", measurement),
      subtitle = "Comparison with SRFA Standard",
      x = "O/C Ratio",
      y = "H/C Ratio",
      caption = "Unique NP formulas (blue) are those not found in SRFA standard"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    theme_bw(base_size = 14) +
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
  
  # Save plot with optimized dimensions
  filename <- paste0("vK_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(
    plot_object = vk_plot,
    filename = file.path(subdir, filename)
  )
  
  message(paste("  Saved plot:", file.path(subdir, filename)))
}

message(paste("Completed generating", length(unique_measurements), 
              "SRFA comparison diagrams in:", file.path(output_base, subdir)))

# --- Van Krevelen Diagrams: δRI Comparison ---
message("Creating Van Krevelen diagrams for δRI comparison...")

# Create directory structure
subdir <- "van_krevelen/delta_RI_comparison"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Generate plots for each measurement
for (measurement in unique_measurements) {
  # Filter and prepare data
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
  
  if (nrow(plot_data) == 0) {
    warning("No data available for measurement: ", measurement)
    next
  }
  
  # Calculate statistics
  n_common <- sum(!is.na(plot_data$δRI))
  n_unique <- sum(is.na(plot_data$δRI))
  
  message(paste("Processing:", measurement))
  message(paste("  Common formulas:", n_common, paste0("(", round(n_common/nrow(plot_data)*100, 1), "%)")))
  message(paste("  Unique:", n_unique, paste0("(", round(n_unique/nrow(plot_data)*100, 1), "%)")))
  message(paste("  δRI range:", paste(round(range(plot_data$δRI, na.rm = TRUE), 2), collapse = " to ")))
  
  # Create optimized plot
  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    # Unique (red)
    geom_point(
      data = filter(plot_data, is.na(δRI)),
      aes(shape = "Unique"),
      color = "#E41A1C",
      size = 2.5,
      alpha = 0.5
    ) +
    # Common formulas (viridis color scale)
    geom_point(
      data = filter(plot_data, !is.na(δRI)),
      aes(color = δRI),
      size = 3,
      alpha = 0.7
    ) +
    # Viridis color scale
    scale_color_viridis_c(
      option = "D",
      name = expression(paste(Delta, "RI(%)")),
      limits = quantile(plot_data$δRI, c(0.05, 0.95), na.rm = TRUE),
      oob = scales::squish
    ) +
    # Shape scale for Unique
    scale_shape_manual(
      values = c("Unique" = 16),
      name = NULL
    ) +
    labs(
      title = paste("Van Krevelen Diagram:", measurement),
      x = "O/C Ratio",
      y = "H/C Ratio",
      caption = expression(paste("Unique (red) are NP-specific | Color shows ", Delta, "RI(%) for common formulas"))
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.box.spacing = unit(0.2, "cm"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 0.5, 1, 1, "cm"),
      plot.caption = element_text(size = 10, color = "gray40", hjust = 0)
    )
  
  # Save plot with optimized dimensions
  filename <- paste0("vK_deltaRI_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(
    plot_object = vk_plot,
    filename = file.path(subdir, filename)
  )
  
  message(paste("  Saved plot:", file.path(subdir, filename)))
}

message(paste("Completed generating", length(unique_measurements), 
              "δRI comparison diagrams in:", file.path(output_base, subdir)))

# --- Van Krevelen Diagrams: Unique Molecular Formulas ---
message("Creating Van Krevelen diagrams for unique molecular formulas...")

# Create directory structure
subdir <- "van_krevelen/unique_molecular_formulas"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Generate plots for each measurement
for (measurement in unique_measurements) {
  # Filter and prepare data
  plot_data <- formulas_NP %>%
    filter(measurement_name == measurement, !is_SRFA) %>%
    mutate(log_intensity = log10(rel_intensity_permille)) %>%
    arrange(rel_intensity_permille)
  
  if (nrow(plot_data) == 0) {
    warning("No Unique available for measurement: ", measurement)
    next
  }
  
  # Calculate statistics
  intensity_range <- round(range(plot_data$rel_intensity_permille, na.rm = TRUE), 2)
  log_range <- round(range(plot_data$log_intensity, na.rm = TRUE), 2)
  
  message(paste("Processing:", measurement))
  message(paste("  Unique:", nrow(plot_data)))
  message(paste("  Intensity range:", paste(intensity_range, collapse = " to "), "%"))
  message(paste("  Log10 intensity range:", paste(log_range, collapse = " to ")))
  
  # Create optimized plot with proper subscript formatting
  vk_plot <- ggplot(plot_data, aes(x = formula_oc, y = formula_hc)) +
    geom_point(
      aes(color = log_intensity),
      size = 3.5,
      alpha = 0.8
    ) +
    scale_color_viridis_c(
      option = "D",
      name = expression(log[10] * "(RI%)"),  # Subscript 10 in log10
      guide = guide_colorbar(
        barwidth = unit(0.5, "cm"),
        barheight = unit(5, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      title = paste("Unique Molecular Formulas:", measurement),
      subtitle = expression(paste("Colored by ", log[10], " Relative Intensity")),  # Subscript 10
      x = "O/C Ratio",
      y = "H/C Ratio",
      caption = "Showing only formulas unique to nanoparticle samples (not found in SRFA standard)"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 0.5, 1, 1, "cm"),
      plot.caption = element_text(size = 10, color = "gray40", hjust = 0)
    )
  
  # Save plot with optimized dimensions
  filename <- paste0("vK_unique_", gsub("[^[:alnum:]]", "_", measurement), ".pdf")
  save_plot(
    plot_object = vk_plot,
    filename = file.path(subdir, filename)
  )
  
  message(paste("  Saved plot:", file.path(subdir, filename)))
}

message(paste("Completed generating", length(unique_measurements), 
              "unique molecular formula diagrams in:", file.path(output_base, subdir)))

# --- Intensity-Weighted Averages Comparison Plot ---
message("Creating intensity-weighted averages comparison plot...")

# Create directory structure
subdir <- "comparison_plots/intensity_weighted_averages"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Check if required data exists
required_cols <- c("new_name", "weighted_formula_hc", "weighted_formula_oc")
if (!all(required_cols %in% names(IWA_AuNP)) || !all(required_cols %in% names(IWA_P25))) {
  stop("Required columns missing in IWA data frames")
}

# Filter data for selected ratios
selected_ratios <- c("NP-HH", "NP-HM", "NP-HL")
IWA_AuNP_filtered <- IWA_AuNP %>% 
  filter(new_name %in% selected_ratios) %>%
  mutate(NP_type = "AuNP")

IWA_P25_filtered <- IWA_P25 %>% 
  filter(new_name %in% selected_ratios) %>%
  mutate(NP_type = "P25")

# Combine data for plotting
combined_data <- bind_rows(IWA_AuNP_filtered, IWA_P25_filtered)

# Calculate axis limits safely
if (nrow(combined_data) > 0) {
  hc_range <- range(combined_data$weighted_formula_hc, na.rm = TRUE)
  oc_range <- range(combined_data$weighted_formula_oc, na.rm = TRUE)
  
  # Add 5% padding to ranges
  hc_min <- hc_range[1] - 0.05 * diff(hc_range)
  hc_max <- hc_range[2] + 0.05 * diff(hc_range)
  oc_min <- oc_range[1] - 0.05 * diff(oc_range)
  oc_max <- oc_range[2] + 0.05 * diff(oc_range)
  
  message("Calculated axis limits:")
  message(paste("  H/C range:", round(hc_min, 3), "to", round(hc_max, 3)))
  message(paste("  O/C range:", round(oc_min, 3), "to", round(oc_max, 3)))
} else {
  stop("No data available for selected ratios")
}

# Create comparison plot
iwa_plot <- ggplot(
  combined_data,
  aes(x = weighted_formula_hc, y = weighted_formula_oc, color = NP_type)
) +
  geom_point(
    size = 4,
    alpha = 0.7
  ) +
  # Add labels only if there's data
  {if (nrow(combined_data) > 0) {
    geom_text_repel(
      aes(label = new_name),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 20,
      min.segment.length = 0.2,
      seed = 123  # For reproducible label positions
    )
  }} +
  # Color scheme
  scale_color_manual(
    name = "Nanoparticle",
    values = c("AuNP" = "#E41A1C", "P25" = "#377EB8"),
    labels = c("AuNP", "P25")
  ) +
  # Labels and titles
  labs(
    title = "Intensity-Weighted Averages: AuNP vs P25",
    subtitle = "Comparison of H/C and O/C Ratios",
    x = "H/C Ratio",
    y = "O/C Ratio"
  ) +
  # Axis limits with slight buffer
  coord_cartesian(
    xlim = c(hc_min, hc_max),
    ylim = c(oc_min, oc_max)
  ) +
  # Theme
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 3, 1, 1, "cm")
  )

# Save plot
filename <- "IWA_comparison_AuNP_vs_P25.pdf"
save_plot(
  plot_object = iwa_plot,
  filename = file.path(subdir, filename),
  width = 10,
  height = 7
)

message(paste("Successfully saved intensity-weighted averages plot to:", 
              file.path(output_base, subdir, filename)))

# --- Stacked Column Plots ---
message("Creating stacked column plots...")

# Create directory structure
subdir <- "stacked_plots/common_vs_unique_formulas"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Define custom order for x-axis
custom_x_order <- c("NP-HH", "NP-MM", "NP-LL",
                   "NP-HM", "NP-HL", "NP-MH", "NP-LH")

# Transform data with proper ordering and patterns
plot_data <- summary_SRFA %>%
  filter(!is.na(new_name)) %>%  # Filter out NA values
  filter(!str_detect(new_name, "BLK")) %>%  # Filter out blank samples
  select(new_name, Increased_pct, Decreased_pct, Non_SRFA_pct) %>%
  pivot_longer(
    cols = -new_name,
    names_to = "category",
    values_to = "percentage"
  ) %>%
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
    new_name = factor(new_name, levels = custom_x_order)
  )

# Define visual properties
fill_colors <- c(
  "common" = "#3E8E41",  # Forest green
  "unique" = "#0072B2"   # Steel blue
)

pattern_types <- c(
  "increased" = "stripe",
  "decreased" = "circle",
  "none" = "none"
)

# Create plot
stacked_plot <- ggplot(
  plot_data,
  aes(x = new_name, y = percentage, fill = fill_group, pattern = pattern_type)
) +
  geom_col_pattern(
    position = position_stack(reverse = FALSE),
    color = "black",
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_fill = "gray30"
  ) +
  # Add sample totals with reduced spacing
  geom_text(
    data = summary_SRFA %>% distinct(new_name, Total),
    aes(x = new_name, y = -1, label = paste0("n=", Total)),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black",
    vjust = 1.5
  ) +
  scale_fill_manual(
    values = fill_colors,
    labels = c(
      "common" = "Common with SRFA",
      "unique" = "Unique to Nanoparticles"
    ),
    name = "Formula Type"
  ) +
  scale_pattern_manual(
    values = pattern_types,
    labels = c(
      "Increased" = "Increased Intensity",
      "Decreased" = "Decreased Intensity",
      "none" = "No Change"
    ),
    name = "Intensity Change",
    na.translate = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05)),
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    x = NULL,
    y = "Percentage of Molecular Formulas",
    title = "Molecular Formula Changes Relative to SRFA Standard",
    subtitle = "Showing common formulas with intensity changes and Unique",
    caption = "Percentage of common formulas shown at the top of the bar and unique ones at the bottom | Total formulas shown below each bar"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = rel(1.1),
      margin = margin(t = 5)
    ),
    axis.title.y = element_text(
      size = rel(1.2),
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = rel(1), face = "bold"),
    legend.text = element_text(size = rel(0.9)),
    plot.title = element_text(face = "bold", size = rel(1.3)),
    plot.subtitle = element_text(color = "gray40", size = rel(1.1)),
    plot.caption = element_text(color = "gray50", size = rel(0.8))
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(pattern = "none")
    ),
    pattern = guide_legend(
      order = 2,
      override.aes = list(fill = "white")
    )
  ) +
  coord_cartesian(clip = "off")

# Save and report results
filename <- "formula_changes_stacked_patterns.pdf"
save_plot(
  plot_object = stacked_plot,
  filename = file.path(subdir, filename),
  width = 12,
  height = 7
)

message("Created stacked column plot with patterns:")
message("- Common formulas (green) with intensity change patterns")
message("- Unique (blue) with no pattern")
message(paste("-", n_distinct(plot_data$new_name), "sample groups analyzed"))
message("- Patterns show increased (stripes) and decreased (dots) intensity changes")
message(paste("- Plot saved to:", file.path(output_base, subdir, filename)))

# --- Stacked Column Plot by Measurement Name ---
message("Creating stacked column plot by measurement name...")

# Create directory structure
subdir <- "stacked_plots/measurement_comparison"
dir.create(file.path(output_base, subdir), recursive = TRUE, showWarnings = FALSE)

# Define custom order for x-axis
custom_x_order <- c("AuNP-HH", "P25-HH", "AuNP-MM", "P25-MM",
                   "AuNP-LL", "P25-LL", "AuNP-HM", "P25-HM",
                   "AuNP-HL", "P25-HL", "AuNP-MH", "P25-MH",
                   "AuNP-LH", "P25-LH")

# Transform data with proper ordering and patterns
plot_data <- summary_SRFA_measurement %>%
  filter(!str_detect(measurement_name, "BLK")) %>%  # Filter out blank samples
  select(measurement_name, Increased_pct, Decreased_pct, Non_SRFA_pct) %>%
  pivot_longer(
    cols = -measurement_name,
    names_to = "category",
    values_to = "percentage"
  ) %>%
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
    measurement_name = factor(measurement_name, levels = custom_x_order)
  )

# Create plot
stacked_plot <- ggplot(
  plot_data,
  aes(x = measurement_name, y = percentage, fill = fill_group, pattern = pattern_type)
) +
  geom_col_pattern(
    position = position_stack(reverse = FALSE),
    color = "black",
    pattern_density = 0.05,
    pattern_spacing = 0.02,
    pattern_fill = "gray30"
  ) +
  # Add sample totals with reduced spacing
  geom_text(
    data = summary_SRFA_measurement %>% distinct(measurement_name, Total),
    aes(x = measurement_name, y = -1, label = paste0("n=", Total)),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black",
    vjust = 1.5
  ) +
  scale_fill_manual(
    values = fill_colors,
    labels = c(
      "common" = "Common with SRFA",
      "unique" = "Unique to Nanoparticles"
    ),
    name = "Formula Type"
  ) +
  scale_pattern_manual(
    values = pattern_types,
    labels = c(
      "increased" = "Increased Intensity",
      "decreased" = "Decreased Intensity",
      "none" = "No Change"
    ),
    name = "Intensity Change",
    na.translate = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05)),
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    x = NULL,
    y = "Percentage of Molecular Formulas",
    title = "Molecular Formula Changes by Measurement",
    subtitle = "Showing common formulas with intensity changes and Unique",
    caption = "Percentage of common formulas shown at the top of the bar and unique ones at the bottom | Total formulas shown below each bar"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = rel(1.1),
      margin = margin(t = 5)
    ),
    axis.title.y = element_text(
      size = rel(1.2),
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = rel(1), face = "bold"),
    legend.text = element_text(size = rel(0.9)),
    plot.title = element_text(face = "bold", size = rel(1.3)),
    plot.subtitle = element_text(color = "gray40", size = rel(1.1)),
    plot.caption = element_text(color = "gray50", size = rel(0.8))
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(pattern = "none")
    ),
    pattern = guide_legend(
      order = 2,
      override.aes = list(fill = "white")
    )
  ) +
  coord_cartesian(clip = "off")

# Save and report results
filename <- "formula_changes_by_measurement.pdf"
save_plot(
  plot_object = stacked_plot,
  filename = file.path(subdir, filename),
  width = 14,  # Wider to accommodate more x-axis labels
  height = 7
)

message("Created stacked column plot by measurement:")
message("- Common formulas (green) with intensity change patterns")
message("- Unique (blue) with no pattern")
message(paste("-", n_distinct(plot_data$measurement_name), "measurements analyzed"))
message("- Patterns show increased (stripes) and decreased (dots) intensity changes")
message(paste("- Plot saved to:", file.path(output_base, subdir, filename)))

message("\nVisualization complete!")
message("Results saved to output/visualization/")
message("- van_krevelen/: Van Krevelen diagrams for various comparisons")
message("- comparison_plots/: Intensity-weighted average comparison plots")
message("- stacked_plots/: Stacked column plots for formula analysis")
