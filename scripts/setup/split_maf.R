#!/usr/bin/env Rscript

# Split a maf file by disease types
# Default settings reads sample and MAF files from OpenPBTA and creates
# MAF files for the 5 most common histologies

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(magrittr)
  library(dplyr)
})

### Define functions

#' Filter a maf table to a single histology
#'
#' @param histology The histology to filter to
#' @param sample_df Input sample info table 
#'  with histology info as `short_histology`
#'  and sample ids as `Kids_First_Biospecimen_ID`
#' @param maf_df The maf data to filter
#' @param outdir Output directory to write files to
#' @param suffix Suffix for the output file (to be combined with the histology name)
#'
filter_write_maf <- function(histology, sample_df, maf_df, outdir, suffix = ".maf.tsv.gz"){
  # get samples
  sample_ids <- sample_df %>%
    filter(short_histology == histology) %>%
    pull(Kids_First_Biospecimen_ID)
  
  # set output file path
  maf_path <- file.path(outdir, 
                        paste0(gsub(" ", "_", histology), suffix))
  
  # filter maf and write
  maf_df %>%
    filter(Tumor_Sample_Barcode %in% sample_ids) %>%
    readr::write_tsv(maf_path)
}

# end Functions


# Set up options
option_list <- list(
  make_option(
    opt_str = c("--maf", "-m"),
    default = "https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data/release-v21-20210820/pbta-snv-consensus-mutation.maf.tsv.gz",
    type = "character",
    help = "Path of MAF file to be split. Can be .gz compressed."
  ),
  make_option(
    opt_str = c("--outdir", "-o"),
    type = "character",
    default = here::here("data/processed/maf_files"),
    help = "File path where output table will be placed."
  ),
  make_option(
    opt_str = c("--sample_info", "-s"),
    type = "character",
    default = "https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data/release-v21-20210820/pbta-histologies.tsv",
    help = "Sample information file, 
            with `Kids_First_Biospecimen_ID`, `experimental_strategy`,
            and `short_histology` columns",
  ),
  make_option(
    opt_str = c("--n_splits", "-n"),
    type = "numeric",
    default = 5
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# create output directory if it does not exist
if(!dir.exists(opts$outdir)){
  dir.create(opts$outdir)
}

# read and filter sample info
samples_df <- readr::read_tsv(opts$sample_info) %>%
  select(
    Kids_First_Biospecimen_ID,
    short_histology,
    experimental_strategy
  )

histologies <- samples_df %>%
  filter(
    !is.na(short_histology), 
    experimental_strategy == "WGS"
  ) %>%
  count(short_histology, sort = TRUE) %>%
  head(n = opts$n_splits) %>%
  pull(short_histology)


# read maf file
maf_df <- readr::read_tsv(opts$maf)

# filter and write for each type
histologies %>%
  purrr::walk(
    filter_write_maf,
    sample_df = samples_df,
    maf_df = maf_df,
    outdir = opts$outdir
  )
    
