#!/bin/bash

# The annotated BCF file containing all potential SNPs for the simulation.
ANNOTATED_BCF="~/spectre-paper/simulations/filtered_annotated_chr22.bcf"
# The chromosome name as it appears in the BCF file.
CHROM_NAME="1"
# The tree sequence file for Spectre.
CHUNKS_FILE="~/spectre-paper/simulations/trees-relate/chr22_chunknames.txt"
TREEINFO_FILE="~/spectre-paper/simulations/trees-relate/treeinfo_chr22.txt"

# How many successful hits to generate for each parameter combination.
TARGET_HITS_PER_BIN=1
# Allele frequency bins to test.
TARGET_FREQS=(0.10 0.22 0.32 0.50)
# Distance bins (in base pairs) to test.
DISTANCE_BINS_BP=(100000 500000 5000000)
# Tolerance for frequency and distance matching.
FREQ_TOLERANCE=0.02
DISTANCE_TOLERANCE_PCT=0.1 # 10% tolerance for distance

echo "INFO: Simulating true interactions..."
mkdir -p true_interaction_results

for freq_target in "${TARGET_FREQS[@]}"; do
    for dist_target in "${DISTANCE_BINS_BP[@]}"; do

        FREQ_NAME=$(printf "%.0f" $(echo "$freq_target * 100" | bc))pct
        DIST_NAME=$(echo "$dist_target" | awk '{ if ($1 >= 1000000) { printf "%.0fMb", $1/1000000 } else { printf "%.0fkb", $1/1000 } }')
        
        RESULTS_DIR="true_interaction_results/freq_${FREQ_NAME}_dist_${DIST_NAME}"
        mkdir -p "${RESULTS_DIR}"
        
        LOG_FILE="${RESULTS_DIR}/true_interaction_hits.tsv"
        echo -e "snp_A_id\tpos_A\tAF_A\tsnp_B_id\tpos_B\tAF_B\tdistance\tM_snps\tbonferroni_thresh\tp_value\tt_statistic\tstd_beta\tsim_beta" > "${LOG_FILE}"

        echo "------------------------------------------------------------------"
        echo "INFO: Running simulation for freq ~${freq_target} and distance ~${DIST_NAME}"
        
        HIT_COUNTER=0
        TRIAL_COUNTER=0

        # Find suitable SNP pairs
        FREQ_LOWER=$(awk "BEGIN {print ${freq_target} - ${FREQ_TOLERANCE}}")
        FREQ_UPPER=$(awk "BEGIN {print ${freq_target} + ${FREQ_TOLERANCE}}")
        DIST_LOWER=$(awk "BEGIN {print ${dist_target} * (1 - ${DISTANCE_TOLERANCE_PCT})}")
        DIST_UPPER=$(awk "BEGIN {print ${dist_target} * (1 + ${DISTANCE_TOLERANCE_PCT})}")

        # Create a list of candidate SNPs in the target frequency range
        CANDIDATE_SNPS_FILE="${RESULTS_DIR}/candidate_snps.txt"
        bcftools query -i "AF>=${FREQ_LOWER} & AF<=${FREQ_UPPER} & N_ALT==1" -f '%ID\t%POS\n' "${ANNOTATED_BCF}" > "${CANDIDATE_SNPS_FILE}"

        # Loop through candidate SNPs to find pairs
        while read -r snp_A_id snp_A_pos; do
            if [ ${HIT_COUNTER} -ge ${TARGET_HITS_PER_BIN} ]; then break; fi
            ((TRIAL_COUNTER++))

            # Define the search window for SNP B.
            SEARCH_START=$(echo "$snp_A_pos + $DIST_LOWER" | bc | cut -d'.' -f1)
            SEARCH_END=$(echo "$snp_A_pos + $DIST_UPPER" | bc | cut -d'.' -f1)

            # Find a potential SNP B in the window.
            SNP_B_INFO=$(bcftools query -r "${CHROM_NAME}:${SEARCH_START}-${SEARCH_END}" -i "AF>=${FREQ_LOWER} & AF<=${FREQ_UPPER} & N_ALT==1" -f '%ID\t%POS\n' "${ANNOTATED_BCF}" | head -n 1)
            
            if [ -z "${SNP_B_INFO}" ]; then continue; fi # If no pair found, try next SNP A.

            snp_B_id=$(echo "${SNP_B_INFO}" | cut -f1)
            snp_B_pos=$(echo "${SNP_B_INFO}" | cut -f2)

            echo "INFO: [Trial ${TRIAL_COUNTER}] Testing pair: ${snp_A_pos} and ${snp_B_pos}"
            TEMP_DIR=$(mktemp -d)

            BONFERRONI_THRESH=0.000001

            Rscript simulate_true_interaction_pheno.R "${ANNOTATED_BCF}" "${TEMP_DIR}/phenotypes.tsv" "${snp_A_pos}" "${snp_B_pos}" "${CHROM_NAME}"
            if [ $? -ne 0 ]; then echo "WARN: Pheno sim failed. Skipping."; rm -rf "${TEMP_DIR}"; continue; fi

            Rscript analyse_interaction_hit.R "${ANNOTATED_BCF}" "${TEMP_DIR}/phenotypes.tsv" "${snp_A_pos}" "${snp_B_pos}" "${BONFERRONI_THRESH}" "${CHROM_NAME}" "${TEMP_DIR}/hit_details.tsv"
            
            if [ $? -eq 0 ]; then
                ((HIT_COUNTER++))
                HIT_DETAILS=$(cat "${TEMP_DIR}/hit_details.tsv")
                
                AF_A=$(bcftools query -r ${CHROM_NAME}:${snp_A_pos} -f '%AF' ${ANNOTATED_BCF})
                AF_B=$(bcftools query -r ${CHROM_NAME}:${snp_B_pos} -f '%AF' ${ANNOTATED_BCF})
                DISTANCE=$(echo "$snp_B_pos - $snp_A_pos" | bc)

                echo -e "${snp_A_id}\t${snp_A_pos}\t${AF_A}\t${snp_B_id}\t${snp_B_pos}\t${AF_B}\t${DISTANCE}\t${M}\t${BONFERRONI_THRESH}\t${HIT_DETAILS}" >> "${LOG_FILE}"
                echo "INFO: [Hit ${HIT_COUNTER}/${TARGET_HITS_PER_BIN}] ✅ True interaction simulated."

                t_statistic=$(echo "${HIT_DETAILS}" | cut -f2)
                SPECTRE_OUTPUT_DIR="~/spectre-paper/simulations/true-interactions/${RESULTS_DIR}"

                echo "INFO: Running Spectre for this hit..."
                python ~/spectre/spectre.py \
                    --output_dir "${SPECTRE_OUTPUT_DIR}" \
                    --chromosomes 22,22 \
                    --positions "${snp_A_pos},${snp_B_pos}" \
                    --chunknames "${CHUNKS_FILE}" \
                    --treeinfo "${TREEINFO_FILE}" \
                    --teststatistic "${t_statistic}" \
                    --alpha "${BONFERRONI_THRESH}" \
                    --maxeffectsize 2.0 \
                    --overwrite
            fi
            rm -rf "${TEMP_DIR}"
        done < <(sort -R "${CANDIDATE_SNPS_FILE}")
    done
done

echo "INFO: Done."

