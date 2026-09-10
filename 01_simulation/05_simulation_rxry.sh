#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="calculate_rxry_script.py"
SCAFFOLDS="../GRCh38_scaffolds.txt"
OUTPUT_DIR="results"
IDXSTATS_DIR="../cami_idxstats"
RXRY_OUTPUT="${OUTPUT_DIR}/hg38_cami_simulation_rxry_output.txt"
SYSTEM="XY"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --idxstats_dir ${IDXSTATS_DIR} \
    --system ${SYSTEM} \
    --homogametic_id NC_000023.11 \
    --heterogametic_id NC_000024.10 \
    --output ${RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

