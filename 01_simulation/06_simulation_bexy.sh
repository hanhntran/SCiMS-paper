#!/bin/bash

IDXSTATS_DIR="./01_simulation/mapped_reads"
OUTPUT_DIR="01_simulation/bexy_output"

mkdir -p ${OUTPUT_DIR}

ls ${IDXSTATS_DIR}/*.1000x.idxstats > ${OUTPUT_DIR}/simulation_master.txt

# run bexy
scaffolds="data/ref_genome/GRCh38_scaffolds.txt"

bexy infer --idxstats "${OUTPUT_DIR}/simulation_master.txt" --keepScaffolds $scaffolds

