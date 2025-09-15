# =============================================================================
# R Package Installation Script
# Purpose: Install all required R packages for nanoparticle SRFA ratio analysis
# Author: Michel Gad
# Date: 2025-09-15
# Description: 
#   - Check for missing R packages and install them
#   - Verify all packages can be loaded successfully
#   - Used by run_all.sh for automated package management
# =============================================================================

# Pin CRAN snapshot to ensure reproducible package versions across machines
snapshot_date <- "2025-09-15"
cran_snapshot <- paste0("https://packagemanager.posit.co/cran/", snapshot_date)
options(repos = c(CRAN = cran_snapshot))

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing missing R packages:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages, dependencies=TRUE)
  } else {
    cat("All required R packages are already installed.\n")
  }
}

# Define required packages
required_packages <- c(
  "tidyverse",      # Data manipulation and visualization
  "data.table",     # Fast data manipulation
  "ggplot2",        # Grammar of graphics plotting
  "dplyr",          # Data manipulation verbs
  "readr",          # Fast file reading
  "stringr",        # String manipulation utilities
  "htmlwidgets",    # HTML widgets for interactive plots
  "knitr",          # Dynamic report generation
  "rmarkdown"       # R Markdown document generation
)

# Install missing packages
cat("Checking R package requirements...\n")
install_if_missing(required_packages)

# Test loading packages
cat("Testing package loading...\n")
failed_packages <- c()

for(pkg in required_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗ Failed to load", pkg, "\n")
    failed_packages <- c(failed_packages, pkg)
  }
}

# Report results
if(length(failed_packages) == 0) {
  cat("\n✓ All R packages are ready!\n")
  cat("Total packages verified:", length(required_packages), "\n")
} else {
  cat("\n✗ Some packages failed to load:\n")
  for(pkg in failed_packages) {
    cat("  -", pkg, "\n")
  }
  cat("\nPlease check the error messages above and try installing manually:\n")
  cat("install.packages(c(", paste0('"', failed_packages, '"', collapse=", "), "))\n")
  quit(status = 1)
}

cat("\nR package setup complete!\n")
