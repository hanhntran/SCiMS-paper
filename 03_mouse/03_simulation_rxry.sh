#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
METADATA="./03_mouse/mouse_metadata.tsv"
SCAFFOLDS="./03_mouse/mouse_scaffold.txt"
OUTPUT_DIR="./03_mouse/results"
IDXSTATS_DIR="./03_mouse/mapped_reads"
MASTER_RXRY_OUTPUT="${OUTPUT_DIR}/mouse_rxry_output.txt"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Generate master idxstats file list
ls ${IDXSTATS_DIR}/*.1000x.idxstats > mouse_idxstats.txt

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --metadata ${METADATA} \
    --master_file mouse_idxstats.txt \
    --system XY \
    --homogametic_id CM001013.3 \
    --heterogametic_id CM001014.3 \
    --output ${MASTER_RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

