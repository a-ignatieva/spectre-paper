#!/bin/bash

# The BCF file used in the original experiment to find hits
DISCOVERY_BCF="filtered_annotated_chr22.bcf"
# The RAW BCF file for the replication cohort
REPLICATION_BCF_RAW="sim_replicate_chr22.bcf"
# The chromosome name used in the BCF files
CHROM_NAME="1"

# Path to the file containing the significant hits
SIGNIFICANT_HITS_FILE="results/25_pct/significant_hits.tsv"

# Output file for the replication results
REPLICATION_RESULTS_FILE="replication_summary.tsv"

echo "INFO: 🚀 Starting Position-Based Replication Analysis..."

REPLICATION_BCF="sim_replicate_chr22.annotated.bcf"

# Annotate the replication BCF with allele frequencies if it hasn't been done already.
if [ ! -f "${REPLICATION_BCF}" ]; then
    echo "INFO: [SETUP] Annotated replication BCF ('${REPLICATION_BCF}') not found. Creating it from '${REPLICATION_BCF_RAW}'..."
    bcftools +fill-tags "${REPLICATION_BCF_RAW}" -Ob -o "${REPLICATION_BCF}" -- -t AF
    echo "INFO: [SETUP] Indexing the new annotated replication BCF..."
    bcftools index "${REPLICATION_BCF}"
fi

# Check for discovery index
if [ ! -f "${DISCOVERY_BCF}.csi" ]; then
    echo "INFO: [SETUP] Index for discovery BCF not found. Creating index..."
    bcftools index "${DISCOVERY_BCF}"
fi

# Validate that the CHROM_NAME exists in the replication BCF file
AVAILABLE_CHROMS=$(bcftools index -s "${REPLICATION_BCF}" | cut -f 1)
if ! echo "${AVAILABLE_CHROMS}" | grep -q -w "${CHROM_NAME}"; then
    echo "FATAL ERROR: The chromosome name '${CHROM_NAME}' was not found in '${REPLICATION_BCF}'."
    exit 1
fi

echo -e "causal_snp_id\tinteracting_snp_A\tinteracting_snp_B\treplication_p_value\treplication_std_beta\treplication_t_statistic\treplicated_successfully\tstatus" > "${REPLICATION_RESULTS_FILE}"

LOOP_COUNTER=0
tail -n +2 "${SIGNIFICANT_HITS_FILE}" | while IFS=$'\t' read -r causal_snp_id causal_snp_pos AF_causal snp_A_id pos_A AF_A snp_B_id pos_B AF_B p_value neg_log10_p t_statistic std_beta bonferroni_threshold; do
    
    ((LOOP_COUNTER++))
    UNIQUE_SEED=$(( $$ + LOOP_COUNTER ))

    echo "INFO: Testing Hit: Causal POS=${causal_snp_pos}, Pair POS=${pos_A}*${pos_B}"
    TEMP_DIR=$(mktemp -d)

    bcftools view "${REPLICATION_BCF}" -r "${CHROM_NAME}:${causal_snp_pos}" -Ov -o "${TEMP_DIR}/causal_snp_replication.vcf"

    VARIANT_COUNT=$(grep -vc '^#' "${TEMP_DIR}/causal_snp_replication.vcf")
    if [ "${VARIANT_COUNT}" -eq 0 ]; then
        echo "WARN: Causal SNP at POS ${causal_snp_pos} not found in replication BCF. Skipping."
        echo -e "${causal_snp_id}\t${snp_A_id}\t${snp_B_id}\tNA\tNA\tNA\tFALSE\tcausal_snp_not_found" >> "${REPLICATION_RESULTS_FILE}"
        rm -rf "${TEMP_DIR}"
        continue
    fi

    Rscript simulate_phenotype.R "${TEMP_DIR}/causal_snp_replication.vcf" "${TEMP_DIR}/phenotypes_replication.tsv" "${UNIQUE_SEED}"
    
    if [ ! -f "${TEMP_DIR}/phenotypes_replication.tsv" ]; then
        echo "WARN: Phenotype simulation failed for causal SNP at ${causal_snp_pos}. Skipping."
        echo -e "${causal_snp_id}\t${snp_A_id}\t${snp_B_id}\tNA\tNA\tNA\tFALSE\tpheno_sim_failed" >> "${REPLICATION_RESULTS_FILE}"
        rm -rf "${TEMP_DIR}"
        continue
    fi
    
    REPLICATION_OUTPUT=$(Rscript test_single_interaction.R "${REPLICATION_BCF}" "${TEMP_DIR}/phenotypes_replication.tsv" "${pos_A}" "${pos_B}" "${CHROM_NAME}")
    EXIT_CODE=$?

    if [ ${EXIT_CODE} -ne 0 ]; then
        echo "WARN: ❌ R script failed for pair at POS ${pos_A}*${pos_B}. Test not possible (likely missing SNPs or monomorphic)."
        echo -e "${causal_snp_id}\t${snp_A_id}\t${snp_B_id}\tNA\tNA\tNA\tFALSE\tnot_testable" >> "${REPLICATION_RESULTS_FILE}"
    elif [ -n "${REPLICATION_OUTPUT}" ]; then
        echo "INFO: ✅ Replication SUCCESSFUL for pair at POS ${pos_A}*${pos_B}"
        echo -e "${causal_snp_id}\t${snp_A_id}\t${snp_B_id}\t${REPLICATION_OUTPUT}\tTRUE\tsuccess" >> "${REPLICATION_RESULTS_FILE}"
    else
        echo "INFO: ❌ Replication FAILED for pair at POS ${pos_A}*${pos_B}"
        echo -e "${causal_snp_id}\t${snp_A_id}\t${snp_B_id}\tNA\tNA\tNA\tFALSE\tnot_significant" >> "${REPLICATION_RESULTS_FILE}"
    fi

    rm -rf "${TEMP_DIR}"
done

echo "INFO: Done"
echo "INFO: Results saved to ${REPLICATION_RESULTS_FILE}"
