suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

COR_THRESHOLD <- 0.1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript test_interactions.R <window_vcf> <pheno_file> <hit_output> <causal_vcf>", call. = FALSE)
}

vcf_file <- args[1]
pheno_file <- args[2]
hit_output_file <- args[3]
causal_vcf_file <- args[4]

# Read data
phenotypes <- read_tsv(pheno_file, show_col_types = FALSE)
vcf_window <- read.vcfR(vcf_file, verbose = FALSE)
gt_window <- extract.gt(vcf_window, element = "GT", as.numeric = TRUE)

vcf_causal <- read.vcfR(causal_vcf_file, verbose = FALSE)
gt_causal <- extract.gt(vcf_causal, element = "GT", as.numeric = TRUE)[1, ]

if (is.null(gt_window) || nrow(gt_window) < 2) {
  quit(save = "no", status = 1)
}

effect_size <- 0.2/sd(phenotypes$phenotype)
cat(paste0("--- Simulated scaled effect size: ", effect_size, " ---\n"))
phenotypes$phenotype_scaled <- scale(phenotypes$phenotype)

num_snps <- nrow(gt_window)
vcf_fix <- getFIX(vcf_window)
snp_ids <- vcf_fix[, "ID"]
snp_pos <- vcf_fix[, "POS"]
snp_af <- extract.info(vcf_window, "AF", as.numeric = TRUE)
causal_af <- extract.info(vcf_causal, "AF", as.numeric = TRUE)

# The Bonferroni correction is based on the total number of pairs
# in the LD-pruned input VCF.
num_total_pairs <- num_snps * (num_snps - 1) / 2
bonferroni_threshold <- 0.05 / num_total_pairs

num_cores <- detectCores() - 1
if (num_cores < 1) num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(paste0("--- Pre-screening ", num_total_pairs, " total pairs on ", num_cores, " CPU cores with r^2 threshold > ", COR_THRESHOLD, " ---\n"))
cat(paste0("--- Using Bonferroni threshold p < ", bonferroni_threshold, " ---\n"))

good_pairs_list <- foreach(i = 1:(num_snps - 1), .combine = 'rbind') %dopar% {
    inner_results <- list()
    for (j in (i + 1):num_snps) {
        snp_A <- gt_window[i, ]
        snp_B <- gt_window[j, ]
        
        carriers_A <- which(snp_A > 0)
        carriers_B <- which(snp_B > 0)
        
        if (all(carriers_A %in% carriers_B) || all(carriers_B %in% carriers_A)) {
            next
        }
        
        interaction_gt <- snp_A * snp_B
        
        if (var(interaction_gt) > 0) {
            r_squared <- stats::cor(interaction_gt, gt_causal)^2
            if (r_squared > COR_THRESHOLD) {
                inner_results[[length(inner_results) + 1]] <- c(i, j)
            }
        }
    }
    
    if (length(inner_results) > 0) {
        do.call(rbind, inner_results)
    } else {
        NULL
    }
}

if (is.null(good_pairs_list) || nrow(good_pairs_list) == 0) {
  cat("--- No pairs passed the pre-screening. ---\n")
  stopCluster(cl)
  quit(save = "no", status = 1)
}

cat(paste0("--- ", nrow(good_pairs_list), " pairs passed pre-screening. Proceeding to test... ---\n"))

pheno_ordered <- phenotypes[match(colnames(gt_window), phenotypes$sample_id), ]

results <- foreach(
  i = 1:nrow(good_pairs_list),
  .combine = 'rbind'
) %dopar% {
  idx1 <- good_pairs_list[i, 1]
  idx2 <- good_pairs_list[i, 2]
  
  snp_A_raw <- gt_window[idx1, ]
  snp_B_raw <- gt_window[idx2, ]

  interaction_raw <- snp_A_raw * snp_B_raw
  snp_A_scaled <- scale(snp_A_raw)
  snp_B_scaled <- scale(snp_B_raw)
  interaction_scaled <- scale(interaction_raw)
  
  model <- stats::lm(pheno_ordered$phenotype_scaled ~ snp_A_scaled + snp_B_scaled + interaction_scaled)
  model_summary <- summary(model)
  
  interaction_term_name <- "interaction_scaled"
  
  if (interaction_term_name %in% rownames(model_summary$coefficients)) {
    p_value <- model_summary$coefficients[interaction_term_name, "Pr(>|t|)"]
    
    if (!is.na(p_value) && p_value < bonferroni_threshold) {
      t_stat <- model_summary$coefficients[interaction_term_name, "t value"]
      std_beta <- model_summary$coefficients[interaction_term_name, "Estimate"]
      neg_log10_p <- -log10(p_value)

      data.frame(
        snp_A_id = snp_ids[idx1],
        pos_A = snp_pos[idx1],
        af_A = snp_af[idx1],
        snp_B_id = snp_ids[idx2],
        pos_B = snp_pos[idx2],
        af_B = snp_af[idx2],
        interaction_p_value = p_value,
        neg_log10_p = neg_log10_p,
        t_statistic = t_stat,
        standardised_beta = std_beta
      )
    } else {
      NULL
    }
  } else {
    NULL
  }
}

stopCluster(cl)
cat("--- All filtered pairs tested. ---\n")

if (!is.null(results) && nrow(results) > 0) {
  first_hit <- results[1, ]
  hit_data <- tibble(
    AF_causal = causal_af,
    interacting_snp_A = first_hit$snp_A_id,
    pos_A = first_hit$pos_A,
    AF_A = first_hit$af_A,
    interacting_snp_B = first_hit$snp_B_id,
    pos_B = first_hit$pos_B,
    AF_B = first_hit$af_B,
    interaction_p_value = first_hit$interaction_p_value,
    neg_log10_p = first_hit$neg_log10_p,
    t_statistic = first_hit$t_statistic,
    standardised_beta = first_hit$standardised_beta,
    standardised_effectsize = effect_size,
    bonferroni_threshold = bonferroni_threshold
  )
  write_tsv(hit_data, hit_output_file, col_names = FALSE)
  quit(save = "no", status = 0)
}

quit(save = "no", status = 1)
