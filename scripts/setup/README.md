# Data download and preparation scripts

- `split_maf.R`: Takes a MAF file and splits it into separate files by cancer type, using provided metadata.
By default, uses data directly downloaded from the OpenPBTA project.

- `download_flags.sh`: Downloads a list of frequently mutated genes in exomes (FLAGS), as published in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/ (corrected). Adds a header and selects the top 50 genes from that list.
