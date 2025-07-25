#!/bin/bash

# IMPORTANT: Chromosome name needs to be whatever is in the .bcf file output from stdpopsim simulation
CHROM_NAME="1"

BCF_FILE="sim_data_chr22.bcf"
RECOMB_MAP="HapmapII_GRCh37/genetic_map_GRCh37_chr22.txt"
TARGET_HITS=100
WINDOW_BP=200000
TARGET_FREQS=("0.01" "0.05" "0.10" "0.25")
FREQ_DELTA=0.005

echo "INFO: Starting..."
mkdir -p results

MIN_BP=$(awk 'NR==2 {print $2}' ${RECOMB_MAP})
MAX_BP=$(awk 'END{print $2}' ${RECOMB_MAP})
echo "INFO: [SETUP] Analysis will be constrained to positions ${MIN_BP}-${MAX_BP}."

ANNOTATED_BCF="filtered_annotated_chr22.bcf"
if [ ! -f "${ANNOTATED_BCF}" ]; then
    echo "INFO: [SETUP] Annotating BCF with allele frequencies..."
    bcftools view ${BCF_FILE} | \
    bcftools +fill-tags - -Ob -o ${ANNOTATED_BCF} -- -t all
    bcftools index ${ANNOTATED_BCF}
fi

echo "INFO: [SETUP] Generating lists of candidate SNPs for each frequency..."
for freq in "${TARGET_FREQS[@]}"; do
    LOWER_BOUND=$(echo "${freq} - ${FREQ_DELTA}" | bc)
    UPPER_BOUND=$(echo "${freq} + ${FREQ_DELTA}" | bc)
    FREQ_DIR_NAME=$(printf "%02d" $(echo "scale=0; $freq * 100 / 1" | bc))_pct
    SNP_LIST_FILE="results/${FREQ_DIR_NAME}/candidate_snps.txt"
    mkdir -p "results/${FREQ_DIR_NAME}"

    bcftools query -i "AF>=${LOWER_BOUND} && AF<=${UPPER_BOUND} && POS>=${MIN_BP} && POS<=${MAX_BP} && N_ALT==1" -f '%ID\n' ${ANNOTATED_BCF} > "${SNP_LIST_FILE}"
    echo "INFO: [SETUP] Found $(wc -l < ${SNP_LIST_FILE}) SNPs for ~${freq} AF."
done

for freq in "${TARGET_FREQS[@]}"; do
    FREQ_DIR_NAME=$(printf "%02d" $(echo "scale=0; $freq * 100 / 1" | bc))_pct
    HIT_COUNTER=0
    TRIAL_COUNTER=0
    SNP_LIST_FILE="results/${FREQ_DIR_NAME}/candidate_snps.txt"
    LOG_FILE="results/${FREQ_DIR_NAME}/significant_hits.tsv"
    
    echo -e "causal_snp_id\tcausal_snp_pos\tAF_causal\tinteracting_snp_A\tpos_A\tAF_A\tinteracting_snp_B\tpos_B\tAF_B\tinteraction_p_value\tneg_log10_p\tt_statistic\tbeta_estimated\tbeta_simulated\tbonferroni_threshold" > "${LOG_FILE}"
    echo "INFO: --- Starting frequency category: ${freq} ---"

    for causal_snp_id in $(sort -R "${SNP_LIST_FILE}"); do
        if [ ${HIT_COUNTER} -ge ${TARGET_HITS} ]; then break; fi
        ((TRIAL_COUNTER++))
        echo "INFO: [${freq} AF | Hit ${HIT_COUNTER}/${TARGET_HITS} | Trial ${TRIAL_COUNTER}] Testing Causal SNP: ${causal_snp_id}"

        causal_snp_pos=$(bcftools query -i "ID=='${causal_snp_id}'" -f '%POS' ${ANNOTATED_BCF})
        
        lower_bp=$((causal_snp_pos - WINDOW_BP))
        upper_bp=$((causal_snp_pos + WINDOW_BP))
        
        echo "DEBUG: Defined window ${CHROM_NAME}:${lower_bp}-${upper_bp}"
        TEMP_DIR=$(mktemp -d)
        
        bcftools view ${ANNOTATED_BCF} -i "POS==${causal_snp_pos}" -Ov -o "${TEMP_DIR}/causal_snp.vcf"
        Rscript simulate_phenotype.R "${TEMP_DIR}/causal_snp.vcf" "${TEMP_DIR}/phenotypes.tsv"
        if [ ! -f "${TEMP_DIR}/phenotypes.tsv" ]; then echo "WARN: Pheno sim failed. Skipping."; rm -rf "${TEMP_DIR}"; continue; fi
        
        bcftools view ${ANNOTATED_BCF} -r ${CHROM_NAME}:${lower_bp}-${upper_bp} -i "ID!='${causal_snp_id}' && MAF>=0.01 && N_ALT==1" -Ov -o "${TEMP_DIR}/window_unpruned.vcf"

        # Check if the unpruned window has enough variants for PLINK
        if [ $(bcftools view -H "${TEMP_DIR}/window_unpruned.vcf" | wc -l) -lt 2 ]; then
            echo "INFO: Fewer than 2 variants in window before pruning. Skipping."
            rm -rf "${TEMP_DIR}"
            continue
        fi

        # Run PLINK for LD pruning
        plink --vcf "${TEMP_DIR}/window_unpruned.vcf" --const-fid 0 --indep-pairwise 50 5 0.7 --allow-extra-chr --out "${TEMP_DIR}/ld_prune" > /dev/null 2>&1

        # Create the final VCF with only the pruned set of SNPs
        bcftools view "${TEMP_DIR}/window_unpruned.vcf" -i 'ID=@'"${TEMP_DIR}/ld_prune.prune.in"'' -Ov -o "${TEMP_DIR}/window_final.vcf"

        # Count and print the number of variants after pruning
        pruned_count=$(bcftools view -H "${TEMP_DIR}/window_final.vcf" | wc -l)
        echo "DEBUG: ${pruned_count} variants remain after LD pruning."
        
        if [ $(bcftools view -H "${TEMP_DIR}/window_final.vcf" | wc -l) -lt 2 ]; then
            echo "INFO: Fewer than 2 variants in window for interaction test. Skipping."
            rm -rf "${TEMP_DIR}"
            continue
        fi
            
        Rscript test_interactions.R "${TEMP_DIR}/window_final.vcf" "${TEMP_DIR}/phenotypes.tsv" "${TEMP_DIR}/hit_details.tsv" "${TEMP_DIR}/causal_snp.vcf"
        
        if [ $? -eq 0 ]; then
            ((HIT_COUNTER++))
            HIT_DETAILS=$(cat "${TEMP_DIR}/hit_details.tsv")
            echo -e "${causal_snp_id}\t${causal_snp_pos}\t${HIT_DETAILS}" >> "${LOG_FILE}"
            echo "INFO: [${freq} AF | Hit ${HIT_COUNTER}/${TARGET_HITS}] ✅ SIGNIFICANT INTERACTION FOUND"
        fi
        rm -rf "${TEMP_DIR}"
    done
done
echo "INFO: Done"
