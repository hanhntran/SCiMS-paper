#!/bin/bash
IDXSTATS_DIR="./04_chicken/mapped_reads"
OUTPUT_DIR="./04_chicken/results"
SCAFFOLDS_FILE="./data/ref_genome/chicken_scaffold.txt"
METADATA_FILE="./04_chicken/chicken_metadata.tsv"
Z_scaffold="CM028522.1"
W_scaffold="CM028521.1"
id_colummn="Run"

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
