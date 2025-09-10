#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
SCAFFOLDS="./04_chicken/chicken_scaffold.txt"
OUTPUT_DIR="./04_chicken/results"
IDXSTATS_DIR="./04_chicken/mapped_reads"
RXRY_OUTPUT="${OUTPUT_DIR}/chicken_rxry_output.txt"
SYSTEM="ZW"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Run RxRy analysis

python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --idxstats_dir ${IDXSTATS_DIR} \
    --system ${SYSTEM} \
    --homogametic_id CM050376.1 \
    --heterogametic_id CM050375.1 \
    --output ${RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

