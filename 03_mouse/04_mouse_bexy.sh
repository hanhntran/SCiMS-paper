#!/bin/bash

IDXSTATS_DIR="./03_mouse/mapped_reads"
OUTPUT_DIR="./03_mouse/results/bexy"

mkdir -p ${OUTPUT_DIR}

ls ${IDXSTATS_DIR}/*.1000x.idxstats > ${OUTPUT_DIR}/mouse_idxstats.txt

mouse_scaffolds="data/ref_genome/mouse_scaffolds.txt"

# run bexy
bexy infer  --idxstats ${OUTPUT_DIR}/mouse_idxstats.txt  --keepScaffolds $mouse_scaffolds --maxNumThreads 16

