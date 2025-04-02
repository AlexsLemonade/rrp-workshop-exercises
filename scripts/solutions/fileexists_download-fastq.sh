#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Define study ID
STUDY_ID="SRP255885"

# Define file names
FASTQ_R1="SRR11518889_1.fastq.gz"
FASTQ_R2="SRR11518889_2.fastq.gz"
FASTQ_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/089/SRR11518889"

# Define and create destination directories for FASTQ files to live in
FASTQ_DEST="../data/raw/fastq/${STUDY_ID}"
TRIMMED_DIR="../data/trimmed/${STUDY_ID}"
REPORTS_DIR="../reports/fastp"

mkdir -p $FASTQ_DEST $TRIMMED_DIR $REPORTS_DIR

##### Process the R1 file #####

# If the file doesn't exist, then we need to download it
if [ ! -e "$FASTQ_DEST/$FASTQ_R1" ]; then

    # Print an indicator:
    echo "Obtaining $FASTQ_R1"

    # Curl the file (using one of several approaches)
    curl -O $FASTQ_URL/$FASTQ_R1 # this approach preserves the original internet file name

    # Move the file to its destination directory
    mv $FASTQ_R1 $FASTQ_DEST

fi

# Explore: how many lines are in the file?
echo "The number of lines in $FASTQ_R1 is:"
gunzip -c $FASTQ_DEST/$FASTQ_R1 | wc -l



##### Process the R2 file #####

# If the file doesn't exist, then we need to download it
if [ ! -e "$FASTQ_DEST/$FASTQ_R2" ]; then

    # Print an indicator:
    echo "Obtaining $FASTQ_R2"

    # Curl the file (using one of several approaches)
    curl -O $FASTQ_URL/$FASTQ_R2 # this approach preserves the original internet file name

    # Move the file to its destination directory
    mv $FASTQ_R2 $FASTQ_DEST

fi

# Explore: how many lines are in the file?
echo "The number of lines in $FASTQ_R2 is:"
gunzip -c $FASTQ_DEST/$FASTQ_R2 | wc -l


## Trim the files with fastp
fastp \
  --in1 $FASTQ_DEST/$FASTQ_R1 \
  --in2 $FASTQ_DEST/$FASTQ_R2 \
  --out1 $TRIMMED_DIR/$FASTQ_R1 \
  --out2 $TRIMMED_DIR/$FASTQ_R2 \
  --html "$REPORTS_DIR/${STUDY_ID}_report.html" \
  --json "$REPORTS_DIR/${STUDY_ID}_report.json"



############### Downloading files with curl ###################
##### A few alternative approaches for downloading with curl are shown below, with R1 for example.
##### Depending on how the workshop is moving along, multiple versions can be optionally introduced.
##### The versions shown below will place the downloaded file into the destination directory
#####    as part of the same step without needing a separate line to mv.
##### For any version, you can always add `-s` to keep it quiet!

# Approach 2: Download directly to the location you want
#curl $FASTQ_URL/$FASTQ_R1 > $FASTQ_DEST/$FASTQ_R1

# Approach 3: Same as above, but without stdout
#curl $FASTQ_URL/$FASTQ_R1 -o $FASTQ_DEST/$FASTQ_R1
##############################################################
