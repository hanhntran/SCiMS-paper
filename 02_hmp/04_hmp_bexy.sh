#!/bin/bash

IDXSTATS_DIR="./02_hmp/mapped_reads"
OUTPUT_DIR="./02_hmp/results/bexy"

mkdir -p ${OUTPUT_DIR}

ls ${IDXSTATS_DIR}/*.1000x.idxstats > ${OUTPUT_DIR}/hmp_idxstats.txt

scaffolds="data/ref_genome/GRCh38_scaffolds.txt"

# run bexy
bexy infer  --idxstats ${OUTPUT_DIR}/hmp_idxstats.txt --keepScaffolds $scaffolds --maxNumThreads 16


