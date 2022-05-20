#!/usr/bin/env Rscript

# Count the number of samples mutated for each gene in a MAF file.

# This script reads in a MAF file and writes out a table of the genes with
# mutations. The table includes the number of samples that have at least one 
# mutation in the gene, as well as the total number of mutations detected across
# all samples.

# Option descriptions:
#
# --maf :  File path to MAF file to be analyzed. Can be .gz compressed.
# --outfile : The path of the output file to create
# --vaf: Minimum variant allele fraction of mutations to include.
# --min_depth: Minimum sequencing depth to call mutations.
# --include_syn: Flag to include synonymous mutations in counts.

# Example invocation:
#
# Rscript count_gene_mutations.R \
#   --maf mutations.maf.tsv.gz \
#   --outfile gene_counts.tsv

# Load libraries
library(optparse)
library(magrittr)

# Set up options
option_list <- list(
  make_option(
    opt_str = c("--maf", "-m"),
    type = "character",
    default = NA,
    help = "File path of MAF file to be analyzed. Can be .gz compressed."
  ),
  make_option(
    opt_str = c("--outfile", "-o"),
    type = "character",
    default = "gene_counts.tsv",
    help = "File path where output table will be placed."
  ),
  make_option(
    opt_str = c("--vaf", "-v"),
    type = "numeric",
    default = 0.05,
    help = "Minimum variant allele fraction to include (default 0.05)",
    metavar = "numeric"
  ),
  make_option(
    opt_str = c("--min_depth", "-d"),
    type = "numeric",
    default = 0,
    help = "Minimum sequencing depth to include (default 0)",
    metavar = "numeric"
  ),
  make_option(
    opt_str = "--include_syn",
    action = "store_true",
    default = FALSE,
    help = "Include synonymous coding mutations"
  )
)
# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# Check input files
if(!file.exists(opts$maf)){
  stop("The specified MAF file does not exist.")
}

if(!is.null(opts$exclude_genes) && file.exists(opts$exclude_genes)){
  stop("The specified 'excludes_genes' file does not exist.")
}

# Read input MAF file
maf_df <- readr::read_tsv(opts$maf)

# Get the excluded genes list if provided
if(!is.null(opts$exclude_genes)){
  exclude_genes <- readr::read_tsv(opts$exclude_file) %>%
    dplyr::pull("gene")
} else {
  exclude_genes <- c()
}

# Define the MAF `Consequence` values we are interested in and their classification
#  based on definitions in http://asia.ensembl.org/Help/Glossary?id=535
syn_class <- c(
  "Silent",
  "Start_Codon_Ins",
  "Start_Codon_SNP",
  "Stop_Codon_Del",
  "De_novo_Start_InFrame",
  "De_novo_Start_OutOfFrame"
)
nonsyn_class <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

# select mutations to keep
include_class <- nonsyn_class
if(opts$include_syn){
  include_class <- c(include_class, syn_class)
}



# Process the maf table
muts_df <- maf_df %>%
  # select only the fields we need
  dplyr::select(
    Tumor_Sample_Barcode,
    Hugo_Symbol,
    Entrez_Gene_Id,
    Variant_Classification,
    Variant_Type,
    t_depth,
    t_ref_count,
    t_alt_count
  ) %>%
  # calculate VAF
  dplyr::mutate(vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
  # filter by VAF, min depth & Classification
  dplyr::filter(
    vaf >= opts$vaf,
    t_ref_count + t_alt_count >= opts$min_depth,
    Variant_Classification %in% include_class
  )

# count mutations by sample and gene
sample_gene_counts <- muts_df %>%
  dplyr::count(Tumor_Sample_Barcode, Hugo_Symbol, name = "mut_count")

# count mutations by gene
gene_counts <- sample_gene_counts %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(
    mutated_samples = dplyr::n(),
    total_muts = sum(mut_count)
  ) %>%
  # sort genes by sample count, then total (descending)
  dplyr::arrange(desc(mutated_samples), desc(total_muts))

  
# Write output
readr::write_tsv(gene_counts, file = opts$outfile)