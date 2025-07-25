# analyse_interaction_hit.R
# Takes a SNP pair and phenotype, runs regression, and outputs stats if the result is significant.

suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript analyse_interaction_hit.R <bcf_file> <pheno_file> <snp_A_pos> <snp_B_pos> <bonferroni_thresh> <chrom_name> <output_hit_file>", call. = FALSE)
}

bcf_file <- args[1]
pheno_file <- args[2]
snp_A_pos <- args[3]
snp_B_pos <- args[4]
bonferroni_threshold <- as.numeric(args[5])
chrom <- args[6]
output_hit_file <- args[7]

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

# Load data
phenotypes <- read_tsv(pheno_file, show_col_types = FALSE)
gt_pair <- extract.gt(vcf_pair, element = "GT", as.numeric = TRUE)
fix_block <- getFIX(vcf_pair)

snp_A_raw <- gt_pair[fix_block[, "POS"] == snp_A_pos, ]
snp_B_raw <- gt_pair[fix_block[, "POS"] == snp_B_pos, ]
interaction_raw <- snp_A_raw * snp_B_raw

# Check for zero variance
if (var(snp_A_raw) == 0 || var(snp_B_raw) == 0 || var(interaction_raw) == 0) {
  file.remove(temp_vcf_file)
  quit(save = "no", status = 1)
}

# Standardise all variables for regression
pheno_ordered <- phenotypes[match(colnames(gt_pair), phenotypes$sample_id), ]
phenotype_scaled <- scale(pheno_ordered$phenotype)
effectsize <- 0.2/sd(pheno_ordered$phenotype)
snp_A_scaled <- scale(snp_A_raw)
snp_B_scaled <- scale(snp_B_raw)
interaction_scaled <- scale(interaction_raw)

# Run the linear model
model <- lm(phenotype_scaled ~ snp_A_scaled + snp_B_scaled + interaction_scaled)
model_summary <- summary(model)
interaction_term_name <- "interaction_scaled"

if (interaction_term_name %in% rownames(model_summary$coefficients)) {
  p_value <- model_summary$coefficients[interaction_term_name, "Pr(>|t|)"]
  
  # Check against the threshold
  if (!is.na(p_value) && p_value < bonferroni_threshold) {
    t_stat <- model_summary$coefficients[interaction_term_name, "t value"]
    std_beta <- model_summary$coefficients[interaction_term_name, "Estimate"]
    
    # Write to output file
    tibble(
      p_value = p_value,
      t_statistic = t_stat,
      standardised_beta = std_beta,
      effectsize = effectsize
    ) %>% write_tsv(output_hit_file, col_names = FALSE)
    
  } else {
    # If not significant, exit with an error code so the shell script knows
    file.remove(temp_vcf_file)
    quit(save = "no", status = 1)
  }
} else {
  file.remove(temp_vcf_file)
  quit(save = "no", status = 1) # Fail if interaction term is missing
}

file.remove(temp_vcf_file)
