#!/bin/bash
set -euo pipefail


##### Process the R1 file #####

# Print an indicator:
echo "Obtaining NC16_S1_L004_R1_001.fastq.gz"

# Curl the file (using one of several approaches, scroll down for more)
curl -O https://sra-download.ncbi.nlm.nih.gov/traces/sra63/SRZ/011518/SRR11518889/NC16_S1_L004_R1_001.fastq.gz


# Explore: how many lines are in the file?
echo "The number of lines in NC16_S1_L004_R1_001.fastq.gz is:"
gunzip -c NC16_S1_L004_R1_001.fastq.gz | wc -l
  
# Move the file to its destination directory
mv NC16_S1_L004_R1_001.fastq.gz ../data/raw/fastq/SRP255885/


##### Process the R2 file #####

# Print an indicator:
echo "Obtaining NC16_S1_L004_R2_001.fastq.gz"

# Curl the file (using one of several approaches)
curl -O https://sra-download.ncbi.nlm.nih.gov/traces/sra63/SRZ/011518/SRR11518889/NC16_S1_L004_R2_001.fastq.gz

# Explore: how many lines are in the file?
echo "The number of lines in NC16_S1_L004_R2_001.fastq.gz is:"
gunzip -c NC16_S1_L004_R2_001.fastq.gz | wc -l
  
# Move the file to its destination directory
mv NC16_S1_L004_R2_001.fastq.gz ../data/raw/fastq/SRP255885/




############### Downloading files with curl ###################
##### A few alternative approaches for downloading files with curl are shown below, with R1 for example. 
##### Depending on how the workshop is moving along, multiple versions can be optionally introduced.
##### The versions shown below will place the downloaded file into the destination directory 
#####    as part of the same step without needing a separate line to mv.
##### For any version, you can always add `-s` to keep it quiet!

# Approach 2: Download directly to the location you want
#curl https://sra-download.ncbi.nlm.nih.gov/traces/sra63/SRZ/011518/SRR11518889/NC16_S1_L004_R2_001.fastq.gz > ../data/raw/fastq/SRP255885/NC16_S1_L004_R1_001.fastq.gz

# Approach 3: Same as above, but without stdout
#curl https://sra-download.ncbi.nlm.nih.gov/traces/sra63/SRZ/011518/SRR11518889/NC16_S1_L004_R2_001.fastq.gz -o ../data/raw/fastq/SRP255885/NC16_S1_L004_R1_001.fastq.gz
##############################################################

