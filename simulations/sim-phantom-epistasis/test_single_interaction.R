# test_single_interaction.R (with explicit chromosome argument)
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr)) # For str_extract

SIGNIFICANCE_THRESHOLD <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript test_single_interaction.R <bcf_file> <pheno_file> <snp_A_pos> <snp_B_pos> <chrom_name>", call. = FALSE)
}

bcf_file <- args[1]
pheno_file <- args[2]
snp_A_pos <- args[3]
snp_B_pos <- args[4]
chrom <- args[5] # UPDATED: Get chromosome name from argument

temp_vcf_file <- tempfile(fileext = ".vcf")

system(paste0(
  "bcftools view '", bcf_file, "' -r ", chrom, ":", snp_A_pos, ",", chrom, ":", snp_B_pos, " -Ov -o ", temp_vcf_file
), intern = FALSE)

vcf_pair <- try(read.vcfR(temp_vcf_file, verbose = FALSE), silent = TRUE)
if (inherits(vcf_pair, "try-error") || nrow(getFIX(vcf_pair)) < 2) {
  file.remove(temp_vcf_file)
  quit(save = "no", status = 1)
}

phenotypes <- read_tsv(pheno_file, show_col_types = FALSE)
gt_pair <- extract.gt(vcf_pair, element = "GT", as.numeric = TRUE)
fix_block <- getFIX(vcf_pair)

snp_A_raw <- gt_pair[fix_block[, "POS"] == snp_A_pos, ]
snp_B_raw <- gt_pair[fix_block[, "POS"] == snp_B_pos, ]
interaction_raw <- snp_A_raw * snp_B_raw

if (var(snp_A_raw) == 0 || var(snp_B_raw) == 0 || var(interaction_raw) == 0) {
  file.remove(temp_vcf_file)
  quit(save = "no", status = 1)
}

pheno_ordered <- phenotypes[match(colnames(gt_pair), phenotypes$sample_id), ]
phenotype_scaled <- scale(pheno_ordered$phenotype)
snp_A_scaled <- scale(snp_A_raw)
snp_B_scaled <- scale(snp_B_raw)
interaction_scaled <- scale(interaction_raw)

model <- lm(phenotype_scaled ~ snp_A_scaled + snp_B_scaled + interaction_scaled)
model_summary <- summary(model)
interaction_term_name <- "interaction_scaled"

if (interaction_term_name %in% rownames(model_summary$coefficients)) {
  p_value <- model_summary$coefficients[interaction_term_name, "Pr(>|t|)"]
  
  if (!is.na(p_value) && p_value < SIGNIFICANCE_THRESHOLD) {
    t_stat <- model_summary$coefficients[interaction_term_name, "t value"]
    std_beta <- model_summary$coefficients[interaction_term_name, "Estimate"]
    cat(paste(p_value, std_beta, t_stat, sep = "\t"))
  }
}

file.remove(temp_vcf_file)
