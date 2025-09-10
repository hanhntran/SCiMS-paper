#!/bin/bash
IDXSTATS_DIR="./04_chicken/mapped_reads"
OUTPUT_DIR="./04_chicken/results"
SCAFFOLDS_FILE="./data/ref_genome/chicken_scaffolds.txt"
METADATA_FILE="./04_chicken/chicken_metadata.tsv"
Z_scaffold="CM050376.1"
W_scaffold="CM050375.1"
id_column="Run"

# Run SCiMS
scims \
  --idxstats_folder ${IDXSTATS_DIR} \
  --scaffolds ${SCAFFOLDS_FILE} \
  --homogametic_scaffold ${Z_scaffold} \
  --heterogametic_scaffold ${W_scaffold} \
  --ZW \
  --metadata ${METADATA_FILE} \
  --id_column ${id_column} \
  --output_dir ${OUTPUT_DIR}
