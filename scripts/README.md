The `scripts` directory is where we might store scripts that are not associated with a particular analysis, or may be used across multiple different analyses.

In this case, the `setup` subdirectory contains scripts that were used to download and preprocess some of the data files included with the repository.

The `solutions/` directory contains two versions of a "solution" script for the `download-fastq.sh` that will be written interactively during the workshop.
One version (`variables_download-fastq.sh`) makes use of bash variables, and the other version (`novariables_download-fastq.sh`) does not. 
