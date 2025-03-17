#!/bin/bash

# Variables
THRESHOLD="0.95"
RXRY_SCRIPT="../scripts/rxry_script.py"
METADATA="../01_simulation/simulation_metadata.tsv"
SCAFFOLDS="../data/ref_genome/GRCh38_scaffolds.txt"
OUTPUT_DIR="../01_simulation/results"
IDXSTATS_DIR="../01_simulation/mapped_reads"
MASTER_RXRY_OUTPUT="${OUTPUT_DIR}/simulation_rxry_output.txt"
SYSTEM="XY"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Generate master idxstats file list
ls ${IDXSTATS_DIR}/*.1000x.idxstats > simulation_master.txt

# Run RxRy analysis
python3 ${RXRY_SCRIPT} \
    --scaffolds ${SCAFFOLDS} \
    --metadata ${METADATA} \
    --master_file simulation_master.txt \
    --system XY \
    --homogametic_id NC_000023.11 \
    --heterogametic_id NC_000024.10 \
    --output ${MASTER_RXRY_OUTPUT} \
    --threshold ${THRESHOLD}

