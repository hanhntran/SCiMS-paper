#!/bin/bash

cd $(pwd)
module load anaconda

conda activate wgsim

# ----- CONFIGURE THESE PATHS -----
MALE_FASTA="reference_genome/hg38_male_XY.fasta"
FEMALE_FASTA="reference_genome/hg38_female_XX.fasta"

OUT_DIR="simulated_reads_raw"

# Simulation Parameters
READ_PAIRS=100000000   # 100 Million pairs = 200 Million total reads
READ_LEN=150          # 150bp Illumina read length
ERROR_RATE=0.01       # 1% flat sequencing error rate
# ----- END CONFIG -----

mkdir -p ${OUT_DIR}

echo "Starting wgsim for MALE (XY) genome..."
wgsim -N ${READ_PAIRS} -1 ${READ_LEN} -2 ${READ_LEN} -e ${ERROR_RATE} \
    ${MALE_FASTA} \
    ${OUT_DIR}/simulated_male_R1.fq \
    ${OUT_DIR}/simulated_male_R2.fq

echo "Male (XY) simulation complete."
echo "----------------------------------------"

echo "Starting wgsim for FEMALE (XX) genome..."
wgsim -N ${READ_PAIRS} -1 ${READ_LEN} -2 ${READ_LEN} -e ${ERROR_RATE} \
    ${FEMALE_FASTA} \
    ${OUT_DIR}/simulated_female_R1.fq \
    ${OUT_DIR}/simulated_female_R2.fq

echo "Female (XX) simulation complete."
echo "----------------------------------------"
echo "All done! Host read datasets are ready in: ${OUT_DIR}"

