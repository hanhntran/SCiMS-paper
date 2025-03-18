#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
METADATA="./04_chicken/chicken_metadata.tsv"
SCAFFOLDS="./04_chicken/chicken_scaffold.txt"
OUTPUT_DIR="./04_chicken/results"
IDXSTATS_DIR="./04_chicken/mapped_reads"
MASTER_RXRY_OUTPUT="${OUTPUT_DIR}/chicken_rxry_output.txt"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Generate master idxstats file list
ls ${IDXSTATS_DIR}/*.1000x.idxstats > chicken_idxstats.txt

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --metadata ${METADATA} \
    --master_file chicken_idxstats.txt \
    --system ZW \
    --homogametic_id CM050376.1 \
    --heterogametic_id CM050375.1 \
    --output ${MASTER_RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

