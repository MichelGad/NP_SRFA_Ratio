## Nanoparticle-SRFA Analysis Pipeline (NP_SRFA_Ratio)

### 🎯 Introduction

This repository contains a complete, reproducible R-based pipeline for analyzing molecular formula patterns in Suwannee River Fulvic Acid (SRFA) samples and comparing them across gold nanoparticles (AuNP) and P25 TiO₂ nanoparticles. The workflow ingests molecular formula data, performs rigorous preprocessing and group-wise aggregation, computes intensity-weighted averages, identifies common vs. unique formulas, generates publication-ready visualizations (including Van Krevelen diagrams), and exports comprehensive results and summaries.

The pipeline is the exact computational sequence used to produce the outputs in the `output/` directory and the final export package in `output/04_export_results/`.

---

### 📄 Project Scope

This project serves as the computational framework for:

- Cleaning and harmonizing molecular formula datasets (odd/even electrons, blanks, replicates)
- Computing reproducibility thresholds and summary statistics
- Comparing AuNP vs. P25 samples and SRFA standards
- Calculating intensity-weighted average (IWA) molecular descriptors
- Visualizing compositional differences and class distributions
- Exporting a complete, documented data package for downstream use

---

### 🛠️ How It Works

The pipeline is organized into four modular R scripts executed in sequence. Each script logs its execution, writes intermediate outputs, and produces plots and tables.

#### 1) Data Preparation

- Script: `scripts/01_data_preparation.R`
- Purpose: Import, clean, and preprocess molecular formula data; perform QC; aggregate replicates; compute reproducibility and summary statistics.
- Inputs: CSV files in `input/` (e.g., `formulas_clean.csv`, `formulas_avg.csv`)
- Outputs (to `output/01_data_preparation/`):
  - `formulas_clean.csv`
  - `formulas_averaged.csv`
  - `formulas_processed.csv`
  - `grouped_formulas.csv`
  - `reproducibility_thresholds.csv`
  - `summary_statistics.csv`

#### 2) Comparison Analysis

- Script: `scripts/02_comparison_analysis.R`
- Purpose: Build AuNP and P25 datasets, subtract blanks, compute IWAs, identify common vs. unique formulas, and generate comparison summaries.
- Inputs: Outputs from Data Preparation
- Outputs (to `output/02_comparison_analysis/`):
  - `combined_nanoparticle_data.csv`
  - `IWA_AuNP.csv`, `IWA_P25.csv`
  - `common_formulas.csv`, `unique_AuNP_formulas.csv`, `unique_P25_formulas.csv`
  - `summary_common.csv`, `summary_SRFA.csv`, `summary_SRFA_measurement.csv`, `STD_SRFA_vs_measurements_summary.csv`

#### 3) Visualization

- Script: `scripts/03_visualization.R`
- Purpose: Create comprehensive visualizations including Van Krevelen diagrams, comparison plots, stacked plots, DOC analysis, and reproducibility analysis; export publication-ready figures.
- Inputs: Outputs from Comparison Analysis + optional inputs (formulas_clean, formulas_processed_all, DOC.csv)
- **Enhanced Structure**: Now features consistent section organization with helper functions for standardized plotting and error handling.
- Outputs (to `output/03_visualization/`):
  - `van_krevelen/`
    - `SRFA_comparison/` - Van Krevelen diagrams comparing with SRFA standard
    - `delta_RI_comparison/` - Van Krevelen diagrams showing ΔRI changes
    - `unique_molecular_formulas/` - Van Krevelen diagrams for unique formulas
    - `clean_bp_RI/` - Van Krevelen diagrams colored by RI(bp) from clean data
    - `avg_MFs_bp/` - Van Krevelen diagrams colored by RI(bp) from averaged data
    - `avg_MFs_RI/` - Van Krevelen diagrams colored by RI(‰) from averaged data
    - `rep_MFs/` - Van Krevelen diagrams with reproducibility thresholds
  - `comparison_plots/`
    - `intensity_weighted_averages/IWA_comparison_AuNP_vs_P25.pdf`
    - `DOC_analysis/` - DOC analysis plots for each sample type
  - `stacked_plots/`
    - `common_vs_unique_formulas/formula_changes_stacked_patterns.pdf`
    - `measurement_comparison/formula_changes_by_measurement.pdf`
  - `reproducibility/`
    - `reproducibility_logRIdiff.pdf`
    - `reproducibility_pctRIdiff.pdf`

#### 4) Export Results

- Script: `scripts/04_export_results.R`
- Purpose: Consolidate outputs, create final data package, generate reports and documentation.
- Inputs: All prior outputs
- Outputs (to `output/04_export_results/`):
  - `data/` complete tables (processed, grouped, IWAs, summaries, dictionary)
  - `plots/` visualization folders copied from `output/03_visualization/` (full tree)
  - `reports/analysis_report.txt`
  - `analysis_summary.csv`, `final_summary_statistics.csv`, `README.md`

---

### 🎨 Enhanced Visualization Structure

The visualization script (`03_visualization.R`) has been completely restructured for better organization and consistency:

#### **Consistent Section Organization**
- Each visualization type now has its own dedicated section
- Standardized section headers using `create_section_header()`
- Consistent error handling and validation using `validate_optional_data()`
- Uniform completion messages across all sections

#### **Helper Functions**
- `create_vk_theme()` - Consistent Van Krevelen plot styling
- `validate_optional_data()` - Standardized optional data validation
- `create_section_header()` - Uniform section headers

#### **Van Krevelen Diagram Variants**
Each Van Krevelen variant is now a separate, well-organized section:
- **SRFA Comparison** - Comparison with SRFA standard
- **ΔRI Comparison** - Common vs unique formulas with intensity changes
- **Unique Formulas** - Unique NP molecular formulas colored by log10 intensity
- **Clean bp RI** - Clean data colored by RI(bp)
- **Averaged MFs bp** - Averaged data colored by RI(bp)
- **Averaged MFs RI** - Averaged data colored by RI(‰)
- **Reproducibility MFs** - Diagrams with reproducibility thresholds

#### **Additional Visualizations**
- **DOC Analysis** - Comprehensive DOC analysis plots for each sample type
- **Intensity-Weighted Averages** - AuNP vs P25 comparison
- **Stacked Plots** - Common vs unique and measurement comparison plots
- **Reproducibility Analysis** - Violin plots for reproducibility assessment

---

### 🚀 Getting Started

Follow these steps to set up and run the full Ratio pipeline.

#### 1) Project Setup

Clone the repository and enter the project directory:

```bash
git clone https://github.com/MichelGad/NP_SRFA_Ratio.git
cd NP_SRFA_Ratio
```

Create the required directories if they do not exist:

```bash
mkdir -p input output
```

Place required CSV inputs in `input/`:
- `input/formulas_clean.csv` (required)
- `input/formulas_avg.csv` (required)
- `input/DOC.csv` (optional - for DOC analysis plots)

**Required files** should contain the molecular formula data exported from upstream processing.

**Optional files** enable additional visualizations:
- `DOC.csv` - Enables DOC analysis plots for each sample type
- Additional optional inputs are automatically detected and used when available

#### 2) Automated Environment Setup

Use the provided script to install R packages pinned in `PACKAGES.md` and run the pipeline end-to-end.

```bash
# Make executable and run
chmod +x run_all.sh
./run_all.sh
```

What this does:
- ✅ Verifies R is available
- ✅ Installs required R packages via `install_r_packages.R`
- ✅ Creates output directories
- ✅ Runs all scripts in sequence with logging
- ✅ Produces a completion summary and organized outputs

#### 3) Manual Execution (if preferred)

Install R packages:

```bash
Rscript install_r_packages.R
```

Run scripts individually in order:

```bash
Rscript "scripts/01_data_preparation.R"
Rscript "scripts/02_comparison_analysis.R"
Rscript "scripts/03_visualization.R"
Rscript "scripts/04_export_results.R"
```

---

### 📋 Output and Logs

All key outputs are organized under `output/`. Each run generates a unique log directory at `output/logs/<RUN_ID>/` with one log per script and a `run_metadata.txt` file (environment, git info, input checksums). A high-level summary is written to `output/analysis_completion_summary.txt` and includes the Run ID and log directory path.

#### Output Structure

```
output/
├── analysis_completion_summary.txt
├── 01_data_preparation/
│   ├── formulas_averaged.csv
│   ├── formulas_clean.csv
│   ├── formulas_processed.csv
│   ├── grouped_formulas.csv
│   ├── reproducibility_thresholds.csv
│   └── summary_statistics.csv
├── 02_comparison_analysis/
│   ├── combined_nanoparticle_data.csv
│   ├── common_formulas.csv
│   ├── IWA_AuNP.csv
│   ├── IWA_P25.csv
│   ├── STD_SRFA_vs_measurements_summary.csv
│   ├── summary_common.csv
│   ├── summary_SRFA.csv
│   ├── summary_SRFA_measurement.csv
│   ├── unique_AuNP_formulas.csv
│   └── unique_P25_formulas.csv
├── 03_visualization/
│   ├── comparison_plots/
│   │   ├── intensity_weighted_averages/
│   │   │   └── IWA_comparison_AuNP_vs_P25.pdf
│   │   └── DOC_analysis/
│   │       ├── DOC_analysis_SRFA.pdf
│   │       ├── DOC_analysis_TiO₂.pdf
│   │       ├── DOC_analysis_AuNP.pdf
│   │       └── DOC_analysis_comparison.pdf
│   ├── stacked_plots/
│   │   ├── common_vs_unique_formulas/
│   │   │   └── formula_changes_stacked_patterns.pdf
│   │   └── measurement_comparison/
│   │       └── formula_changes_by_measurement.pdf
│   ├── van_krevelen/
│   │   ├── SRFA_comparison/
│   │   ├── delta_RI_comparison/
│   │   ├── unique_molecular_formulas/
│   │   ├── clean_bp_RI/
│   │   ├── avg_MFs_bp/
│   │   ├── avg_MFs_RI/
│   │   └── rep_MFs/
│   └── reproducibility/
│       ├── reproducibility_logRIdiff.pdf
│       └── reproducibility_pctRIdiff.pdf
├── logs/
│   └── <RUN_ID>/
│       ├── 1.__data_prep.log
│       ├── 2.__comparison_analysis.log
│       ├── 3.__visualization.log
│       ├── 4.__export_results.log
│       └── run_metadata.txt
└── 04_export_results/
    ├── analysis_summary.csv
    ├── data/
    │   ├── combined_nanoparticle_data.csv
    │   ├── data_dictionary.csv
    │   ├── formula_class_analysis.csv
    │   ├── formulas_processed.csv
    │   ├── grouped_formulas.csv
    │   ├── IWA_AuNP.csv
    │   ├── IWA_P25.csv
    │   ├── IWA_summary.csv
    │   ├── nanoparticle_comparison_summary.csv
    │   ├── reproducibility_thresholds.csv
    │   ├── srfa_analysis_summary.csv
    │   ├── STD_comparison_summary.csv
    │   ├── summary_common.csv
    │   ├── summary_SRFA.csv
    │   └── summary_SRFA_measurement.csv
    ├── final_summary_statistics.csv
    ├── plots/
    │   ├── comparison_plots/
    │   ├── stacked_plots/
    │   ├── reproducibility/
    │   └── van_krevelen/
    ├── README.md
    └── reports/
        └── analysis_report.txt
```

#### Monitoring Execution

When running `run_all.sh`, you'll see real-time progress. Each run writes logs to `output/logs/<RUN_ID>/` and the completion summary includes the Run ID and log directory.

Useful commands:

```bash
# View completion summary (includes Run ID and log dir)
cat output/analysis_completion_summary.txt

# List the most recent run's log directory
ls -td output/logs/*/ | head -1

# Inspect a specific script log from the latest run
LATEST_RUN=$(ls -td output/logs/*/ | head -1)
cat "${LATEST_RUN}/1.__data_prep.log" | sed -n '1,120p'

# Search for errors across the latest run's logs
grep -i "error" "${LATEST_RUN}"/*.log || true
```

---

### 📚 Dependencies

R 4.0+ is required. R package versions are fixed via a CRAN snapshot in `install_r_packages.R` to ensure reproducibility across machines. The package list is documented in `PACKAGES.md`.

#### 📊 R Libraries

| Library       | Version | Description                                                                 |
| :------------ | :------ | :-------------------------------------------------------------------------- |
| tidyverse     | 2.0.0   | Collection of core tidyverse packages (dplyr, tidyr, ggplot2, readr, purrr, forcats, stringr, tibble) |
| data.table    | 1.17.8  | Enhanced data frames for efficient large data operations                    |
| FactoMineR    | 2.12    | Multivariate exploratory data analysis                                     |
| factoextra    | 1.0.7   | Visualization of multivariate analysis results                             |
| corrr         | 0.4.4   | Tidy correlation analysis                                                   |
| RColorBrewer  | 1.1-3   | Color palettes for attractive plots                                         |

Install with:

```bash
Rscript install_r_packages.R
```

---

### 🔧 Troubleshooting

#### 1) Permission Denied Error

```bash
# Error: Permission denied: ./run_all.sh
# Solution: Make the script executable
chmod +x run_all.sh
```

#### 2) Missing Input Files

```bash
# Error: Missing required input files
# Solution: Place required CSVs in input/ directory
# Required files:
# - input/formulas_clean.csv
# - input/formulas_avg.csv
```

#### 3) R Not Found

```bash
# Error: Rscript is not installed or not in PATH
# Solution: Install R from https://www.r-project.org/
```

#### 4) Missing R Packages

```bash
# Error: Missing R packages
# Solution: The run_all.sh script installs them automatically. For manual installation:

# Install with the provided script (recommended)
Rscript install_r_packages.R

# Or install manually from CRAN
Rscript -e "install.packages(c('tidyverse','data.table','ggplot2','dplyr','readr','stringr','htmlwidgets','knitr','rmarkdown'), repos='https://cran.rstudio.com/', dependencies=TRUE)"
```

Note: R package versions are fixed via a CRAN snapshot in `install_r_packages.R` to ensure reproducibility.

#### 5) Script Execution Failures

```bash
# View overall completion summary (has Run ID and log dir)
cat output/analysis_completion_summary.txt

# Inspect the latest run's logs
LATEST_RUN=$(ls -td output/logs/*/ | head -1)
ls -l "${LATEST_RUN}"
for f in "${LATEST_RUN}"/*.log; do echo "---- ${f}"; tail -n 50 "$f"; done

# Search for error patterns across the latest run's logs
grep -i "error" "${LATEST_RUN}"/*.log || true
grep -i "failed" "${LATEST_RUN}"/*.log || true
```

#### 6) Memory Issues

```bash
# For large datasets, increase memory limits
# Windows (run before scripts):
Rscript -e "memory.limit(size = 16000)"

# Linux/macOS:
export R_MAX_MEM_SIZE=16G
```

#### 7) Partial Pipeline Execution

```bash
# Run individual scripts if a step fails
Rscript "scripts/01_data_preparation.R"
Rscript "scripts/02_comparison_analysis.R"
Rscript "scripts/03_visualization.R"
Rscript "scripts/04_export_results.R"

# Check which outputs were produced
cat output/analysis_completion_summary.txt 2>/dev/null || true
```

#### 8) Quick Checks

```bash
# Validate expected files exist after each step
ls output/01_data_preparation/ | wc -l
ls output/02_comparison_analysis/ | wc -l

# Confirm core figures were produced
test -f output/03_visualization/stacked_plots/common_vs_unique_formulas/formula_changes_stacked_patterns.pdf && echo "Stacked plots: OK"
test -f output/03_visualization/comparison_plots/intensity_weighted_averages/IWA_comparison_AuNP_vs_P25.pdf && echo "IWA comparison: OK"
test -f output/03_visualization/van_krevelen/SRFA_comparison/ && echo "Van Krevelen diagrams: OK"
test -f output/03_visualization/comparison_plots/DOC_analysis/ && echo "DOC analysis: OK" || echo "DOC analysis: Not available (optional input missing)"
```

---

### 🔬 Scientific Context

This pipeline quantifies and visualizes compositional differences in SRFA-related molecular formula space across AuNP and P25 contexts. It supports:

- Intensity-normalized formula comparisons and descriptor trends
- Identification of common vs. unique molecular formulas
- IWA-based interpretation of compositional shifts
- Publication-ready Van Krevelen and comparison plots

---

### 📖 Citation

If you use this pipeline in your research, please cite appropriately and include a link to this repository and its version/tag used for your analysis.

**To be added**

---

### 🤝 Contributing

Contributions are welcome. You can:
- Report bugs or issues
- Suggest features or enhancements
- Submit pull requests with improvements

Please keep changes modular and documented.

---

### 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

### 👥 Authors

- Michel Gad

### 📞 Contact

For questions, please contact: michel.gad@outlook.com
