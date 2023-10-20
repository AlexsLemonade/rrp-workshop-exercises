#!/bin/bash
set -euo pipefail


##### Process the R1 file #####

# Print an indicator:
echo "Obtaining SRR11518889_1.fastq.gz"

# Curl the file (using one of several approaches, scroll down for more)
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/089/SRR11518889/SRR11518889_1.fastq.gz


# Explore: how many lines are in the file?
echo "The number of lines in SRR11518889_1.fastq.gz is:"
gunzip -c SRR11518889_1.fastq.gz | wc -l
  
# Move the file to its destination directory
mv SRR11518889_1.fastq.gz ../data/raw/fastq/SRP255885/


##### Process the R2 file #####

# Print an indicator:
echo "Obtaining SRR11518889_2.fastq.gz"

# Curl the file (using one of several approaches)
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/089/SRR11518889/SRR11518889_2.fastq.gz

# Explore: how many lines are in the file?
echo "The number of lines in SRR11518889_2.fastq.gz is:"
gunzip -c SRR11518889_2.fastq.gz | wc -l
  
# Move the file to its destination directory
mv SRR11518889_2.fastq.gz ../data/raw/fastq/SRP255885/




############### Downloading files with curl ###################
##### A few alternative approaches for downloading files with curl are shown below, with R1 for example. 
##### Depending on how the workshop is moving along, multiple versions can be optionally introduced.
##### The versions shown below will place the downloaded file into the destination directory 
#####    as part of the same step without needing a separate line to mv.
##### For any version, you can always add `-s` to keep it quiet!

# Approach 2: Download directly to the location you want
#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/089/SRR11518889/SRR11518889_1.fastq.gz > ../data/raw/fastq/SRP255885/SRR11518889_1.fastq.g

# Approach 3: Same as above, but without stdout
#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/089/SRR11518889/SRR11518889_1.fastq.gz -o ../data/raw/fastq/SRP255885/SRR11518889_1.fastq.g
##############################################################

