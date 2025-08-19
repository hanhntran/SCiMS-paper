#!/bin/bash

IDXSTATS_DIR="./03_mouse/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > mouse_idxstats.txt

mouse_scaffolds="data/ref_genome/mouse_scaffolds.txt"

# run bexy
./tools/bexy/build/bexy infer  --idxstats mouse_idxstats.txt  --keepScaffolds $mouse_scaffolds --maxNumThreads 16

