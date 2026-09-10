#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
SCAFFOLDS="./data/ref_genome/GRCh38_scaffold.txt"
OUTPUT_DIR="./02_hmp/results"
IDXSTATS_DIR="./02_hmp/mapped_reads"
RXRY_OUTPUT="${OUTPUT_DIR}/hmp_rxry_output.txt"
SYSTEM="XY"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --idxstats_dir ${IDXSTATS_DIR} \
    --system XY \
    --homogametic_id NC_000023.11 \
    --heterogametic_id NC_000024.10 \
    --output ${RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

