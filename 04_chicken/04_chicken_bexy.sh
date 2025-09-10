#!/bin/bash

IDXSTATS_DIR="./04_chicken/mapped_reads"
OUTPUT_DIR="./04_chicken/results/bexy"

mkdir -p ${OUTPUT_DIR}

ls ${IDXSTATS_DIR}/*.idxstats > ${OUTPUT_DIR}/chicken_idxstats.txt

chicken_scaffolds="data/ref_genome/chicken_scaffolds.txt"

# run bexy
./tools/bexy/build/bexy infer  --idxstats ${OUTPUT_DIR}/chicken_idxstats.txt  --keepScaffolds $chicken_scaffolds --maxNumThreads 16
