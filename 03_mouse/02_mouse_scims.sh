#!/bin/bash
IDXSTATS_DIR="./03_mouse/mapped_reads"
OUTPUT_DIR="./03_mouse/results"
SCAFFOLDS_FILE="./data/ref_genome/mouse_scaffold.txt"
METADATA_FILE="./03_mouse/mouse_metadata.tsv"
X_scaffold="NC_000023.11"
Y_scaffold="NC_000024.10"
id_colummn="Run"

# Run SCiMS
scims \
  --idxstats_folder ${IDXSTATS_DIR} \
  --scaffolds ${SCAFFOLDS_FILE} \
  --homogametic_scaffold ${X_scaffold} \
  --heterogametic_scaffold ${Y_scaffold} \
  --metadata ${METADATA_FILE} \
  --id_column ${id_column} \
  --output_dir ${OUTPUT_DIR}
