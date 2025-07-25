# simulate_phenotype.R (with seed argument)
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))

HERITABILITY <- 0.01
BETA <- 0.2

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript simulate_phenotype.R <input_vcf> <output_pheno_file> [optional_seed]", call. = FALSE)
}

vcf_file <- args[1]
output_file <- args[2]

if (length(args) == 3) {
  seed <- as.integer(args[3])
  set.seed(seed)
}

vcf <- read.vcfR(vcf_file, verbose = FALSE)

if (nrow(vcf@gt) == 0) {
  warning("Input VCF is empty. No phenotype will be generated.", call. = FALSE)
  quit(save = "no", status = 1)
}

gt <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
genotypes <- gt[1, ]

n_individuals <- length(genotypes)

# Genetic component
genetic_effect <- BETA * genotypes
genetic_variance <- var(genetic_effect)

# Environmental component
if (genetic_variance == 0) {
  environmental_variance <- 1
} else {
  environmental_variance <- (genetic_variance * (1 - HERITABILITY)) / HERITABILITY
}
environmental_effect <- rnorm(n_individuals, mean = 0, sd = sqrt(environmental_variance))

# Final phenotype
phenotype <- genetic_effect + environmental_effect

pheno_df <- tibble(
  sample_id = colnames(gt),
  phenotype = phenotype
)

write_tsv(pheno_df, output_file)
