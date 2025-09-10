#!/bin/bash

IDXSTATS_DIR="./04_chicken/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > chicken_idxstats.txt

chicken_scaffolds="data/ref_genome/chicken_scaffolds.txt"

# run bexy
./tools/bexy/build/bexy infer  --idxstats chicken_idxstats.txt  --keepScaffolds $chicken_scaffolds --maxNumThreads 16
