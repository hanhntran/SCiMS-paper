#!/bin/bash

# Set directories consistent with Snakefile outputs
IDXSTATS_DIR="./02_hmp/mapped_reads"
OUTPUT_DIR="./02_hmp/results"
SCAFFOLDS_FILE="./data/ref_genome/GRCh38_scaffolds.txt"
METADATA_FILE="./02_hmp/dbGap_metadata.txt"
X_scaffold="NC_000023.11"
Y_scaffold="NC_000024.10"
id_column="Run"

# Run SCiMS
scims \
  --idxstats_folder ${IDXSTATS_DIR} \
  --scaffolds ${SCAFFOLDS_FILE} \
  --homogametic_scaffold ${X_scaffold} \
  --heterogametic_scaffold ${Y_scaffold} \
  --metadata ${METADATA_FILE} \
  --id_column ${id_column} \
  --output_dir ${OUTPUT_DIR}
