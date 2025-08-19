#!/bin/bash

IDXSTATS_DIR="./01_simulation/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > simulation_master.txt

# run bexy
scaffolds="data/ref_genome/GRCh38_scaffolds.txt"

bexy infer --idxstats "simulation_master.txt" --keepScaffolds $scaffolds

# generate bexy output in R
Rscript ./01_simulation/bexy.R

