#!/bin/bash

IDXSTATS_DIR="./02_hmp/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > hmp_idxstats.txt
scaffolds="data/ref_genome/GRCh38_scaffolds.txt"

# run bexy
infer  --idxstats hmp_idxstats.txt --keepScaffolds $scaffolds --maxNumThreads 16


