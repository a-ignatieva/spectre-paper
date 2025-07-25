# simulate_true_interaction_pheno.R
# Simulates a phenotype with a true interaction effect.

suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))

BETA_INTERACTION <- 0.2
HERITABILITY <- 0.01 # h^2

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript simulate_true_interaction_pheno.R <bcf_file> <output_pheno_file> <snp_A_pos> <snp_B_pos> <chrom_name>", call. = FALSE)
}

bcf_file <- args[1]
output_file <- args[2]
snp_A_pos <- args[3]
snp_B_pos <- args[4]
chrom <- args[5] # Get chromosome name from argument

# Extract the two interacting SNPs
temp_vcf_file <- tempfile(fileext = ".vcf")
system(paste0(
  "bcftools view '", bcf_file, "' -r ", chrom, ":", snp_A_pos, ",", chrom, ":", snp_B_pos, " -Ov -o ", temp_vcf_file
), intern = FALSE)

vcf_pair <- try(read.vcfR(temp_vcf_file, verbose = FALSE), silent = TRUE)
if (inherits(vcf_pair, "try-error") || nrow(getFIX(vcf_pair)) < 2) {
  file.remove(temp_vcf_file)
  quit(save = "no", status = 1)
}

# Get genotypes and create the interaction term
gt_pair <- extract.gt(vcf_pair, element = "GT", as.numeric = TRUE)
fix_block <- getFIX(vcf_pair)
snp_A_raw <- gt_pair[fix_block[, "POS"] == snp_A_pos, ]
snp_B_raw <- gt_pair[fix_block[, "POS"] == snp_B_pos, ]
interaction_term <- snp_A_raw * snp_B_raw

# Simulate phenotype
n_individuals <- ncol(gt_pair)
genetic_effect <- BETA_INTERACTION * interaction_term
genetic_variance <- var(genetic_effect)

if (genetic_variance == 0) {
  environmental_variance <- 1
} else {
  environmental_variance <- (genetic_variance * (1 - HERITABILITY)) / HERITABILITY
}
environmental_effect <- rnorm(n_individuals, mean = 0, sd = sqrt(environmental_variance))

phenotype <- genetic_effect + environmental_effect

pheno_df <- tibble(
  sample_id = colnames(gt_pair),
  phenotype = phenotype
)
write_tsv(pheno_df, output_file)

file.remove(temp_vcf_file)
