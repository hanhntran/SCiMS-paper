#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="./scripts/rxry_script.py"
METADATA="./02_hmp/dbGap_metadata.txt"
SCAFFOLDS="./data/ref_genome/GRCh38_scaffold.txt"
OUTPUT_DIR="./02_hmp/results"
IDXSTATS_DIR="./02_hmp/mapped_reads"
MASTER_RXRY_OUTPUT="${OUTPUT_DIR}/hmp_rxry_output.txt"
SYSTEM="XY"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Generate master idxstats file list
ls ${IDXSTATS_DIR}/*.1000x.idxstats > hmp_idxstats.txt

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --metadata ${METADATA} \
    --master_file hmp_idxstats.txt \
    --system XY \
    --homogametic_id NC_000023.11 \
    --heterogametic_id NC_000024.10 \
    --output ${MASTER_RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

