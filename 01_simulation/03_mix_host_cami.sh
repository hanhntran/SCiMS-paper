#!/bin/bash 
cd $(pwd)
module load anaconda
conda activate seqtk

# ----- CONFIGURE THESE PATHS -----
CAMI_R1="cami_oral_metagenome/450_genomes/cami_S004_R1.fq.gz"  
CAMI_R2="cami_oral_metagenome/450_genomes/cami_S004_R2.fq.gz"


HOST_DIR="simulated_reads"
OUTPUT_DIR="cami_mixed_reads"

TASK_FILE="${OUTPUT_DIR}/mixing_tasks.txt"
MANIFEST="${OUTPUT_DIR}/mixing_manifest.txt"
 
DEPTHS=(150 250 350 450 1000 10000 20000)
HOST_FRACTIONS=(0.001 0.01 0.10 0.50 0.95 1)
REPLICATES=100

MAX_CONCURRENT_JOBS=10
# ----- END CONFIG -----

mkdir -p ${OUTPUT_DIR}
rm -f ${TASK_FILE}
 
# ==========================================
# RAM DISK SETUP 
# ==========================================
RAM_DIR="/dev/shm/$SLURM_JOB_ID"
mkdir -p ${RAM_DIR}

echo "Loading massive Zymo files into RAM (${RAM_DIR})..."
cp ${CAMI_R1} ${RAM_DIR}/cami_R1_ram.fq.gz
cp ${CAMI_R2} ${RAM_DIR}/cami_R2_ram.fq.gz

RAM_CAMI_R1="${RAM_DIR}/cami_R1_ram.fq.gz"
RAM_CAMI_R2="${RAM_DIR}/cami_R2_ram.fq.gz"
echo "Files loaded to memory successfully!"
# ==========================================

echo -e "sample_id\tsex\thost_depth\thost_fraction\treplicate\thost_reads\tmicrobial_reads\ttotal_reads\tmixed_R1\tmixed_R2" > ${MANIFEST}

TOTAL_SAMPLES=0

echo "Generating empirical task list..."
 
for SEX_LETTER in M F; do
    SEX=$([ "${SEX_LETTER}" == "M" ] && echo "male" || echo "female")
 
    for DEPTH in "${DEPTHS[@]}"; do
        for FRAC in "${HOST_FRACTIONS[@]}"; do
 
            MICRO_READS=$(awk -v h=$DEPTH -v f=$FRAC 'BEGIN {printf "%.0f", (h/f) - h}')
            TOTAL_READS=$(awk -v h=$DEPTH -v f=$FRAC 'BEGIN {printf "%.0f", h/f}')

            for REP in $(seq 1 ${REPLICATES}); do
 
                # Locate the specific Host R1 file
                HOST_FQ_R1="${HOST_DIR}/S${REP}${SEX_LETTER}${DEPTH}_1.fq"
                HOST_FQ_R2="${HOST_DIR}/S${REP}${SEX_LETTER}${DEPTH}_2.fq"
                
                # Check if host files exist, if not skip this iteration
                if [ ! -f "${HOST_FQ_R1}" ] || [ ! -f "${HOST_FQ_R2}" ]; then 
                    continue 
                fi

                SAMPLE_ID="E_S${REP}${SEX_LETTER}${DEPTH}_hf${FRAC}"
                MIXED_R1="${OUTPUT_DIR}/${SAMPLE_ID}_1.fq"
                MIXED_R2="${OUTPUT_DIR}/${SAMPLE_ID}_2.fq"
 
                # Resume protection
                if [ -f "${MIXED_R1}" ] && [ -f "${MIXED_R2}" ]; then continue; fi

                # Unique seed per replicate, ensures identical reads are pulled for both R1 and R2
                SEED=$((REP * 10000 + DEPTH))

                # Safely subset Zymo and concatenate with Host (No Shuffling to preserve pairs)
                CMD_R1="cat ${HOST_FQ_R1} <(seqtk sample -s ${SEED} ${RAM_CAMI_R1} ${MICRO_READS}) > ${MIXED_R1}"
                CMD_R2="cat ${HOST_FQ_R2} <(seqtk sample -s ${SEED} ${RAM_CAMI_R2} ${MICRO_READS}) > ${MIXED_R2}"
                
                # The && ensures R2 only runs if R1 succeeds
                echo "${CMD_R1} && ${CMD_R2}" >> ${TASK_FILE}
                
                echo -e "${SAMPLE_ID}\t${SEX}\t${DEPTH}\t${FRAC}\t${REP}\t${DEPTH}\t${MICRO_READS}\t${TOTAL_READS}\t${MIXED_R1}\t${MIXED_R2}" >> ${MANIFEST}
                TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
 
            done
        done
    done
done

echo "Task list generated. Running ${MAX_CONCURRENT_JOBS} concurrent jobs directly from RAM..."
cat ${TASK_FILE} | xargs -P ${MAX_CONCURRENT_JOBS} -I {} bash -c "{}"
 
# ==========================================
# RAM DISK CLEANUP (Mandatory)
# ==========================================
echo "Cleaning up RAM disk..."
rm -rf ${RAM_DIR}
# ==========================================

echo ""
echo "========================================"
echo "Empirical Mixing complete."
echo "Total paired samples processed: ${TOTAL_SAMPLES}"
echo "========================================"

