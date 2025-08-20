#!/bin/bash

IDXSTATS_DIR="./01_simulation/mapped_reads"
OUTPUT_DIR="01_simulation/bexy_output"

mkdir ${OUTPUT_DIR}

ls ${IDXSTATS_DIR}/*.1000x.idxstats > ${OUTPUT_DIR}/simulation_master.txt

# run bexy
scaffolds="data/ref_genome/GRCh38_scaffolds.txt"

./tools/bexy/build/bexy infer --idxstats "${OUTPUT_DIR}/simulation_master.txt" --keepScaffolds $scaffolds

# generate bexy output in R
Rscript ./scripts/bexy.R ${OUTPUT_DIR}/

