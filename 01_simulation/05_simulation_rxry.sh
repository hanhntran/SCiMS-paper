#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/calculate_rxry_script.py"
SCAFFOLDS="./01_simulation/simulated_ref/GRCh38_scaffolds.txt"
OUTPUT_DIR="./01_simulation/results"
IDXSTATS_DIR="./01_simulation/mapped_reads"
RXRY_OUTPUT="${OUTPUT_DIR}/simulation_rxry_output.txt"
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

