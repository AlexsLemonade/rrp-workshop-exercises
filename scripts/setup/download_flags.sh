#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

FLAGS_URL="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5706417/bin/12920_2017_309_MOESM3_ESM.txt"
DEST_DIR="../../data/reference"
FLAGS_FILE="$DEST_DIR/FLAGS.txt"

# create header line
echo -e "gene\tcount" > "$FLAGS_FILE"
# download and take the top 50 genes
curl -sN "$FLAGS_URL" | head -n 50 >> "$FLAGS_FILE"
