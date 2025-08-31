#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
SCAFFOLDS="./03_mouse/mouse_scaffold.txt"
OUTPUT_DIR="./03_mouse/results"
IDXSTATS_DIR="./03_mouse/mapped_reads"
RXRY_OUTPUT="${OUTPUT_DIR}/mouse_rxry_output.txt"
SYSTEM="XY"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --idxstats_dir ${IDXSTATS_DIR} \
    --system ${SYSTEM} \
    --homogametic_id CM001013.3 \
    --heterogametic_id CM001014.3 \
    --output ${RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

