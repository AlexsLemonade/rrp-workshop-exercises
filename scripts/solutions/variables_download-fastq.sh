#!/bin/bash
set -euo pipefail

# Define study ID
STUDY_ID="SRP255885"

# Define file names
FASTQ_R1="NC16_S1_L004_R1_001.fastq.gz"
FASTQ_R2="NC16_S1_L004_R2_001.fastq.gz"
FASTQ_URL="https://sra-download.ncbi.nlm.nih.gov/traces/sra63/SRZ/011518/SRR11518889"

# Define and create destination directory for FASTQ files to live in
FASTQ_DEST="../data/raw/fastq/$STUDY_ID/"
mkdir -p $FASTQ_DEST 

##### Process the R1 file #####

# Print an indicator:
echo "Obtaining $FASTQ_R1"

# Curl the file (using one of several approaches)
curl -O $FASTQ_URL/$FASTQ_R1 # this approach preserves the original internet file name

# Explore: how many lines are in the file?
echo "The number of lines in $FASTQ_R1 is:"
gunzip -c $FASTQ_R1 | wc -l

  
# Move the file to its destination directory
mv $FASTQ_R1 $FASTQ_DEST


##### Process the R2 file #####

# Print an indicator:
echo "Obtaining $FASTQ_R2"

# Curl the file (using one of several approaches)
curl -O $FASTQ_URL/$FASTQ_R2 # this approach preserves the original internet file name

# Explore: how many lines are in the file?
gunzip -c $FASTQ_R2 | wc -l
  
# Move the file to its destination directory
mv $FASTQ_R2 $FASTQ_DEST



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


