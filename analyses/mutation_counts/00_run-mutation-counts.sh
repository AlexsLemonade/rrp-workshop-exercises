#!/bin/bash

set -euo pipefail

# Navigate to this script's directory to ensure
#  relative paths are correct
cd "$(dirname "${BASH_SOURCE[0]}")"

# Define relative path to data from this script
DATA_DIR=../../data/processed/maf_files
RESULTS_DIR=../../results/mutation_counts

# Create the results directory if it doesn't exist
mkdir -p $RESULTS_DIR

# Run the code over each of the following MAF names using a for loop
for NAME in Ependymoma Ganglioglioma LGAT Medulloblastoma; do

    # Define input file
    MAF_FILE="${DATA_DIR}/${NAME}.maf.tsv.gz"

    # Define output file(s)
    RESULT_FILE="${RESULTS_DIR}/${NAME}_gene-mutations.tsv"

    # Run the script to process the MAF file and export a table of gene counts
    Rscript 01_count-gene-mutations.R \
     --maf "$MAF_FILE" \
     --outfile "$RESULT_FILE"

done


# Run the notebook to visualize LGAT and Medulloblastoma counts
#  Using `Rscript -e` allows us to include a very short script (in quotes) directly on
#  the command line without writing a separate script file first
Rscript -e "rmarkdown::render('02_mutation-count-plots.Rmd')"
