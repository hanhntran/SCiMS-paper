#!/bin/bash 
cd $(pwd)
module load anaconda
conda activate seqtk

# ----- CONFIGURE THESE PATHS -----
# Point these to the files wgsim just finished making
MALE_MASTER_R1="simulated_reads_raw/simulated_male_R1.fq"
MALE_MASTER_R2="simulated_reads_raw/simulated_male_R2.fq"

FEMALE_MASTER_R1="simulated_reads_raw/simulated_female_R1.fq"
FEMALE_MASTER_R2="simulated_reads_raw/simulated_female_R2.fq"

OUT_DIR="simulated_reads"
TASK_FILE="${OUT_DIR}/subset_tasks.txt"

DEPTHS=(150 250 350 450 1000 10000 200000)
REPLICATES=100
MAX_CONCURRENT_JOBS=10
# ----- END CONFIG -----

mkdir -p ${OUT_DIR}
rm -f ${TASK_FILE}

echo "Generating task list for Paired-End host replicates..."

for SEX_LETTER in M F; do
    if [ "${SEX_LETTER}" == "M" ]; then
        MASTER_R1="${MALE_MASTER_R1}"
        MASTER_R2="${MALE_MASTER_R2}"
    else
        MASTER_R1="${FEMALE_MASTER_R1}"
        MASTER_R2="${FEMALE_MASTER_R2}"
    fi

    for DEPTH in "${DEPTHS[@]}"; do
        for REP in $(seq 1 ${REPLICATES}); do
            
            # The exact file names the downstream PE script expects
            OUT_FILE_1="${OUT_DIR}/S${REP}${SEX_LETTER}${DEPTH}_1.fq"
            OUT_FILE_2="${OUT_DIR}/S${REP}${SEX_LETTER}${DEPTH}_2.fq"
            
            # Resume protection (checks if BOTH files exist)
            if [ -f "${OUT_FILE_1}" ] && [ -f "${OUT_FILE_2}" ]; then
                continue
            fi
            
            # Unique seed so every replicate gets different host reads
            SEED=$((REP * 10000 + DEPTH))
            
            # Fast random subsampling: Uses the exact same SEED for R1 and R2
            # The '&&' ensures R2 only runs if R1 succeeds
            CMD="seqtk sample -s ${SEED} ${MASTER_R1} ${DEPTH} > ${OUT_FILE_1} && seqtk sample -s ${SEED} ${MASTER_R2} ${DEPTH} > ${OUT_FILE_2}"
            
            echo "$CMD" >> ${TASK_FILE}
            
        done
    done
done

TASK_COUNT=$(wc -l < ${TASK_FILE})

if [ "${TASK_COUNT}" -eq 0 ]; then
    echo "All Paired-End host files already exist!"
    exit 0
fi

echo "Task list created with ${TASK_COUNT} paired samples to generate."
echo "Running ${MAX_CONCURRENT_JOBS} concurrent jobs..."

cat ${TASK_FILE} | xargs -P ${MAX_CONCURRENT_JOBS} -I {} bash -c "{}"

echo "========================================"
echo "Paired-End Host replicate generation complete!"
echo "Files are located in: ${OUT_DIR}"
echo "========================================"
