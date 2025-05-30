#!/bin/bash

IDXSTATS_DIR="./01_simulation/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > simulation_master.txt

# run bexy
./tools/bexy/build/bexy infer  --idxstats simulation_master.txt --keepScaffolds  NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10

# generate bexy output in R
Rscript ./01_simulation/bexy.R

