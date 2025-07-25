#!/bin/bash

FREQ_DIR="path/simulations/false-positives/phantom_epistasis_results/05_pct"
HITS_FILE="${FREQ_DIR}/significant_hits_batch2.tsv"
SPECTRE_OUTPUT_DIR="${FREQ_DIR}/relate_trees"
CHUNKS_FILE="path/simulations/false-positives/trees/chr22_chunknames.txt"
TREEINFO_FILE="path/simulations/false-positives/trees/treeinfo_chr22.txt"
CHROM_NAME="22"

echo "INFO: Starting Spectre run on significant hits from ${HITS_FILE}"

if [ ! -f "${HITS_FILE}" ]; then
    echo "ERROR: Hits file not found at '${HITS_FILE}'"
    exit 1
fi

tail -n +2 "${HITS_FILE}" | while IFS=$'\t' read -r causal_snp_id causal_snp_pos AF_causal snp_A_id pos_A AF_A snp_B_id pos_B AF_B p_value neg_log10_p t_statistic est_beta sim_beta bonferroni_threshold; do

    echo "------------------------------------------------------------------"
    echo "INFO: Analyzing hit: Causal=${causal_snp_pos}, Pair=${pos_A}*${pos_B}"

    CHROMS="${CHROM_NAME},${CHROM_NAME},${CHROM_NAME}"
    POSITIONS="${pos_A},${pos_B},${causal_snp_pos}"
    
    echo "INFO: Running Spectre with the following parameters:"
    echo "  - Chromosomes: ${CHROMS}"
    echo "  - Positions: ${POSITIONS}"
    echo "  - Test Statistic: ${t_statistic}"
    echo "  - Alpha: ${bonferroni_threshold}"
    echo "  - Effectsize: ${est_beta}"
    echo "  - Maxeffectsize: 2.0"
    
    python /Users/ignatiev/Dropbox/projects/spectre-paper/spectre.py \
        --output_dir "${SPECTRE_OUTPUT_DIR}" \
        --chromosomes "${CHROMS}" \
        --positions "${POSITIONS}" \
        --chunknames "${CHUNKS_FILE}" \
        --treeinfo "${TREEINFO_FILE}" \
        --teststatistic "${t_statistic}" \
        --alpha "${bonferroni_threshold}" \
        --effectsize "${est_beta}" \
        --maxeffectsize 2.0 \
        --overwrite

    # Check if Spectre ran successfully
    if [ $? -eq 0 ]; then
        echo "INFO: ✅ Spectre analysis finished successfully for this hit."
    else
        echo "WARN: ❌ Spectre analysis failed for this hit."
    fi

done

echo "------------------------------------------------------------------"
echo "INFO: Done"

