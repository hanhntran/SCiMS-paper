#!/bin/bash

scims call --idxstats_folder cami_idxstats --scaffolds GRCh38_scaffolds.txt \
    --homogametic_id NC_000023.11 --heterogametic_id NC_000024.10  \
    --output_dir hg38_cami_scims_out --metadata hg38_cami_simulation_metadata.txt --id_column Sample 
