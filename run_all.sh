#!/bin/bash

# =============================================================================
# Script: run_all.sh
# Purpose: Execute the complete nanoparticle SRFA ratio analysis pipeline
# Author: Michel Gad
# Date: 2025-09-19
# Description: 
#   - Runs all analysis scripts in sequence
#   - Provides progress updates and error handling
#   - Creates output directories as needed
#   - Generates completion summary
# =============================================================================

# Set script options
set -e  # Exit on any error
set -u  # Exit on undefined variables

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Generate a unique RUN_ID and per-run log directory
RUN_ID=$(date +%Y%m%d_%H%M%S)
LOG_ROOT_DIR="output/logs"
LOG_DIR="${LOG_ROOT_DIR}/${RUN_ID}"

# Capture run metadata
capture_run_metadata() {
    mkdir -p "${LOG_DIR}"
    META_FILE="${LOG_DIR}/run_metadata.txt"

    # Basic metadata
    {
        echo "RUN_ID: ${RUN_ID}"
        echo "START_TIME: $(date -u "+%Y-%m-%dT%H:%M:%SZ")"
        echo "HOSTNAME: $(hostname 2>/dev/null || echo unknown)"
        echo "OS: $(uname -a 2>/dev/null || echo unknown)"
        echo "R_VERSION: $(Rscript -e 'cat(R.version.string)' 2>/dev/null || echo not_found)"
    } > "${META_FILE}"

    # Git metadata (if repo)
    if command -v git >/dev/null 2>&1 && git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        {
            echo "GIT_COMMIT: $(git rev-parse HEAD 2>/dev/null || echo unknown)"
            echo "GIT_BRANCH: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
            git diff --quiet 2>/dev/null; dirty=$?; [ "$dirty" -eq 0 ] && echo "GIT_DIRTY: false" || echo "GIT_DIRTY: true"
        } >> "${META_FILE}"
    fi

    # Input file checksums (simplified to avoid timeout on large files)
    echo "" >> "${META_FILE}"
    echo "INPUT_FILES:" >> "${META_FILE}"
    for f in input/formulas_avg.csv input/formulas_clean.csv; do
        if [ -f "$f" ]; then
            echo "FILE_SIZE=$(stat -f%z "$f" 2>/dev/null || stat -c%s "$f" 2>/dev/null || echo unknown)  $f" >> "${META_FILE}"
        else
            echo "MISSING: $f" >> "${META_FILE}"
        fi
    done
}

# Function to check if R is available
check_r() {
    if ! command -v Rscript &> /dev/null; then
        print_error "Rscript is not installed or not in PATH"
        print_error "Please install R and ensure Rscript is available"
        exit 1
    fi
    print_success "Rscript found: $(which Rscript)"
}

# Function to install R packages
install_r_packages() {
    print_status "Setting up R packages..."
    
    # Check if install script exists
    if [ ! -f "install_r_packages.R" ]; then
        print_error "install_r_packages.R not found"
        print_error "Please ensure the R package installation script exists"
        exit 1
    fi
    
    # Run the R package installation script
    if Rscript install_r_packages.R; then
        print_success "All R packages installed and verified successfully"
    else
        print_error "Failed to install or verify R packages"
        print_error "Please check the error messages above and install packages manually"
        exit 1
    fi
}

# Function to create output directories
create_directories() {
    print_status "Creating output directories..."
    
    directories=(
        "output"
        "${LOG_ROOT_DIR}"
        "${LOG_DIR}"
        "output/01_data_preparation"
        "output/02_comparison_analysis"
        "output/03_visualization"
        "output/03_visualization/van_krevelen"
        "output/03_visualization/van_krevelen/SRFA_comparison"
        "output/03_visualization/van_krevelen/delta_RI_comparison"
        "output/03_visualization/van_krevelen/unique_molecular_formulas"
        "output/03_visualization/van_krevelen/clean_bp_RI"
        "output/03_visualization/van_krevelen/avg_MFs_bp"
        "output/03_visualization/van_krevelen/avg_MFs_RI"
        "output/03_visualization/van_krevelen/rep_MFs"
        "output/03_visualization/comparison_plots"
        "output/03_visualization/comparison_plots/intensity_weighted_averages"
        "output/03_visualization/stacked_plots"
        "output/03_visualization/stacked_plots/common_vs_unique_formulas"
        "output/03_visualization/stacked_plots/measurement_comparison"
        "output/03_visualization/reproducibility"
        "output/04_export_results"
        "output/04_export_results/data"
        "output/04_export_results/plots"
        "output/04_export_results/reports"
    )
    
    for dir in "${directories[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir"
            print_status "Created directory: $dir"
        fi
    done
    
    print_success "All output directories created"
}

# Function to check input files
check_input_files() {
    print_status "Checking input files..."
    
    input_files=(
        "input/formulas_avg.csv"
        "input/formulas_clean.csv"
    )
    
    missing_files=()
    
    for file in "${input_files[@]}"; do
        if [ ! -f "$file" ]; then
            missing_files+=("$file")
        fi
    done
    
    if [ ${#missing_files[@]} -eq 0 ]; then
        print_success "All input files found"
    else
        print_error "Missing input files: ${missing_files[*]}"
        print_error "Please ensure the pre-processed CSV files are in the input/ directory"
        exit 1
    fi
}

# Function to run R script with error handling
run_r_script() {
    local script_name="$1"
    local script_path="$2"
    
    print_status "Running $script_name..."
    
    if [ ! -f "$script_path" ]; then
        print_error "Script not found: $script_path"
        return 1
    fi
    
    # Run the script and capture output into per-run log file
    local base_no_ext="${script_name%.R}"
    local safe_base
    safe_base=$(echo "$base_no_ext" | tr ' /' '__')
    local log_file="${LOG_DIR}/${safe_base}.log"

    if Rscript "$script_path" 2>&1 | tee "$log_file"; then
        print_success "$script_name completed successfully (log: $log_file)"
        return 0
    else
        print_error "$script_name failed"
        print_error "Check the log file: $log_file"
        return 1
    fi
}

# Function to generate completion summary
generate_summary() {
    print_status "Generating completion summary..."
    
    summary_file="output/analysis_completion_summary.txt"
    
    cat > "$summary_file" << EOF
NANOPARTICLE SRFA RATIO ANALYSIS PIPELINE - COMPLETION SUMMARY
============================================================

Analysis Date: $(date)
Pipeline Version: 1.0
Author: Michel Gad
Run ID: ${RUN_ID}
Log Directory: ${LOG_DIR}

EXECUTION SUMMARY:
==================

EOF

    # Report status per script based on run results
    scripts=("01_data_preparation.R" "02_comparison_analysis.R" "03_visualization.R" "04_export_results.R")
    for script in "${scripts[@]}"; do
        base_no_ext="${script%.R}"
        safe_base=$(echo "$base_no_ext" | tr ' /' '__')
        log_path="${LOG_DIR}/${safe_base}.log"
        status="COMPLETED"
        for failed in "${failed_scripts[@]:-}"; do
            if [ "$failed" = "$(basename "$script")" ]; then
                status="FAILED"
                break
            fi
        done
        prefix="✓"
        [ "$status" = "FAILED" ] && prefix="✗"
        echo "${prefix} ${script} - ${status}  (log: ${log_path})" >> "$summary_file"
    done
    
    echo "" >> "$summary_file"
    echo "OUTPUT FILES GENERATED:" >> "$summary_file"
    echo "=======================" >> "$summary_file"
    
    # Count output files
    if [ -d "output" ]; then
        find output -name "*.csv" -o -name "*.pdf" -o -name "*.txt" -o -name "*.md" | wc -l | xargs echo "Total files:" >> "$summary_file"
        echo "" >> "$summary_file"
        echo "Key output directories:" >> "$summary_file"
        echo "- output/logs/: Execution logs for all scripts" >> "$summary_file"
        echo "- output/01_data_preparation/: Data preparation results" >> "$summary_file"
        echo "- output/02_comparison_analysis/: Nanoparticle comparison results" >> "$summary_file"
        echo "- output/03_visualization/: All generated plots and figures" >> "$summary_file"
        echo "- output/04_export_results/: Complete results package" >> "$summary_file"
    fi
    
    echo "" >> "$summary_file"
    echo "NEXT STEPS:" >> "$summary_file"
    echo "===========" >> "$summary_file"
    echo "1. Review the generated plots in output/03_visualization/" >> "$summary_file"
    echo "2. Check the final results in output/04_export_results/" >> "$summary_file"
    echo "3. Read the analysis report in output/04_export_results/reports/" >> "$summary_file"
    echo "4. Use the exported data for further analysis or publication" >> "$summary_file"
    echo "" >> "$summary_file"
    echo "RUN METADATA:" >> "$summary_file"
    echo "=============" >> "$summary_file"
    echo "- See ${LOG_DIR}/run_metadata.txt for environment, git info, and input checksums" >> "$summary_file"
    
    print_success "Completion summary saved to: $summary_file"
}

# Main execution function
main() {
    echo "============================================================================="
    echo "NANOPARTICLE SRFA RATIO ANALYSIS PIPELINE"
    echo "============================================================================="
    echo "Author: Michel Gad"
    echo "Date: $(date)"
    echo "Version: 1.0"
    echo "============================================================================="
    echo ""
    
    # Pre-flight checks
    print_status "Performing pre-flight checks..."
    check_r
    install_r_packages
    create_directories
    check_input_files
    
    echo ""
    print_status "Starting analysis pipeline..."
    echo ""
    
    # Track execution time
    start_time=$(date +%s)
    
    # Run scripts in sequence
    scripts=(
        "scripts/01_data_preparation.R:Data preparation and preprocessing"
        "scripts/02_comparison_analysis.R:Comparative analysis of nanoparticles"
        "scripts/03_visualization.R:Generate plots and visualizations"
        "scripts/04_export_results.R:Export final results and documentation"
    )
    
    failed_scripts=()
    
    for script_info in "${scripts[@]}"; do
        IFS=':' read -r script_path description <<< "$script_info"
        script_name=$(basename "$script_path")
        
        echo "-------------------------------------------------------------------------"
        print_status "Step: $description"
        echo "-------------------------------------------------------------------------"
        
        if run_r_script "$script_name" "$script_path"; then
            print_success "$script_name completed successfully"
        else
            print_error "$script_name failed"
            failed_scripts+=("$script_name")
        fi
        
        echo ""
    done
    
    # Calculate execution time
    end_time=$(date +%s)
    execution_time=$((end_time - start_time))
    hours=$((execution_time / 3600))
    minutes=$(((execution_time % 3600) / 60))
    seconds=$((execution_time % 60))
    
    echo "============================================================================="
    echo "PIPELINE EXECUTION COMPLETE"
    echo "============================================================================="
    
    if [ ${#failed_scripts[@]} -eq 0 ]; then
        print_success "All scripts completed successfully!"
        echo ""
        print_success "Total execution time: ${hours}h ${minutes}m ${seconds}s"
        echo ""
        print_success "Results are available in the output/ directory"
        print_success "Check output/04_export_results/ for the complete results package"
    else
        print_error "Pipeline completed with errors"
        print_error "Failed scripts: ${failed_scripts[*]}"
        echo ""
        print_warning "Partial results may be available in the output/ directory"
        print_warning "Check individual log files for error details"
    fi
    
    # Generate completion summary
    generate_summary
    
    echo ""
    echo "============================================================================="
    echo "For detailed information, see:"
    echo "- README.md: Complete pipeline documentation"
    echo "- output/analysis_completion_summary.txt: Execution summary"
    echo "- output/04_export_results/README.md: Results usage guide"
    echo "============================================================================="
}

# Run main function
main "$@"
