# Required Packages

This document lists all the R packages required for the nanoparticle SRFA ratio analysis pipeline.

# R Packages

## Core Data Manipulation Packages
- **tidyverse**: Collection of R packages for data science
- **data.table**: Fast data manipulation and aggregation
- **dplyr**: Grammar of data manipulation
- **readr**: Fast file reading utilities
- **stringr**: String manipulation utilities

## Visualization and Interactive Packages
- **ggplot2**: Grammar of graphics plotting system
- **htmlwidgets**: HTML widgets for interactive plots (plotly support)

## Specialized Analysis Packages
None required beyond the core and reporting packages listed.

## Report Generation Packages
- **knitr**: Dynamic report generation
- **rmarkdown**: R Markdown document generation

## Installation

All packages are automatically installed when running `run_all.sh`. The script will:

1. **R Package Setup:**
   - Check which R packages are already installed
   - Install any missing packages from CRAN
   - Verify that all packages can be loaded successfully
   - Report any installation or loading failures

## Manual Installation

### R Packages

Install all required R packages:

```r
# Install all required packages
install.packages(c(
  "tidyverse", "data.table", "ggplot2", "dplyr", "readr",
  "stringr", "htmlwidgets", "knitr", "rmarkdown"
), repos="https://cran.rstudio.com/", dependencies=TRUE)
```

Or use the standalone installation script:

```bash
Rscript install_r_packages.R
```

## Package Usage by Script

### R Scripts

#### 01_data_preparation.R
- **tidyverse**: Data manipulation and visualization
- **data.table**: Fast data operations
- **ggplot2**: Plotting system
- **dplyr**: Data manipulation verbs
- **readr**: File reading utilities
- **stringr**: String manipulation for sample names
- **htmlwidgets**: Interactive plot saving
 

#### 02_comparison_analysis.R
- **tidyverse**: Data manipulation and analysis
- **dplyr**: Data grouping and summarization
- **readr**: File I/O operations

#### 03_visualization.R
- **tidyverse**: Data manipulation
- **ggplot2**: Plotting system
- **dplyr**: Data processing for plots
- **readr**: File reading

#### 04_export_results.R
- **tidyverse**: Data manipulation
- **ggplot2**: Plotting system
- **dplyr**: Data processing
- **readr**: File I/O operations
- **knitr**: Report generation
- **rmarkdown**: Document generation

## Package Dependencies
No additional specialized packages are required.

## Troubleshooting

### Common Issues

1. **General Installation Issues**
   - Ensure you have the latest version of R (>= 4.0.0)
   - Some dependencies may require compilation tools
   - On macOS, you may need Xcode command line tools
   - On Linux, you may need development packages

2. **htmlwidgets Issues**
   - Required for saving interactive plots
   - If installation fails, plots will be saved as static images instead

3. **Package Loading Errors**
   - Restart R session after installation
   - Check for conflicting package versions
   - Update R to the latest version if possible

### Manual Package Installation

If automatic installation fails, try installing packages individually:

```r
# Install core packages first
install.packages("tidyverse")
install.packages("data.table")
install.packages("htmlwidgets")

# Install remaining packages
install.packages(c("knitr", "rmarkdown"))
```

## Version Requirements

- **R**: >= 4.0.0
- **tidyverse**: >= 1.3.0
- **data.table**: >= 1.14.0
- **htmlwidgets**: >= 1.5.0

## Verification

After installation, verify all packages can be loaded:

```r
# Test package loading
packages <- c("tidyverse", "data.table", "ggplot2", "dplyr", "readr", 
              "stringr", "htmlwidgets", "knitr", "rmarkdown")

for(pkg in packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("✓", pkg, "loaded successfully\n")
  } else {
    cat("✗ Failed to load", pkg, "\n")
  }
}
```
