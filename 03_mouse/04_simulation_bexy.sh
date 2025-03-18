#!/bin/bash

IDXSTATS_DIR="./03_mouse/mapped_reads"

ls ${IDXSTATS_DIR}/*.1000x.idxstats > mouse_idxstats.txt

# run bexy
./tools/bexy/build/bexy infer  --idxstats mouse_idxstats.txt  --keepScaffolds CM000994.3,CM000995.3,CM000996.3,CM000997.3,CM000998.3,CM000999.3,CM001000.3,CM001001.3,CM001002.3,CM001003.3,CM001004.3,CM001005.3,CM001006.3,CM001007.3,CM001008.3,CM001009.3,CM001010.3,CM001011.3,CM001012.3,CM001013.3,CM001014.3

# generate bexy output in R
Rscript ./03_mouse/bexy.R

