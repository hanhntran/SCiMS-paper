#!/bin/bash 
cd $(pwd)
module load anaconda
module load bowtie2/2.5.1
module load samtools/1.19.2

# ----- CONFIGURE -----
BT2_INDEX="../hg38_scims_simulation/reference_genome/GRCh38_latest_genomic.fna" 
INPUT_DIR="downloaded_fastqs"
MAPPED_DIR="mapped_reads"
MAX_JOBS=4   
THREADS=3 
# ----- END CONFIG -----

mkdir -p ${MAPPED_DIR}

# Build task list from files on disk
TASK_FILE="${MAPPED_DIR}/map_tasks.txt"
rm -f ${TASK_FILE}

TOTAL=0
SKIPPED=0

# PAIRED-END FIX: Loop ONLY through the _1.fq files
for R1 in ${INPUT_DIR}/*_1.fastq; do
    
    # Safety check to avoid empty directory errors
    if [ ! -f "$R1" ]; then continue; fi
    
    # Dynamically find the matching R2 file
    R2="${R1/_1.fq/_2.fq}"
    
    # Clean the sample name (removes the _1.fq so BAMs are named cleanly)
    SAMPLE=$(basename ${R1} _1.fq)
    
    IDXSTATS="${MAPPED_DIR}/${SAMPLE}.idxstats"

    # Skip if already done
    if [ -f "${IDXSTATS}" ]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # PAIRED-END FIX: Write SAMPLE, R1, and R2 to the task file, separated by tabs
    echo -e "${SAMPLE}\t${R1}\t${R2}" >> ${TASK_FILE}
    TOTAL=$((TOTAL + 1))
done

echo "Paired-end samples to map: ${TOTAL}"
echo "Already done: ${SKIPPED}"

if [ ${TOTAL} -eq 0 ]; then
    echo "Nothing to do."
    exit 0
fi

# Function to map a single paired-end sample
map_sample() {
    LINE="$1"
    SAMPLE=$(echo "${LINE}" | cut -f1)
    R1=$(echo "${LINE}" | cut -f2)
    R2=$(echo "${LINE}" | cut -f3)
    
    BT2_IDX="$2"
    OUT_DIR="$3"
    THR="$4"

    IDXSTATS="${OUT_DIR}/${SAMPLE}.idxstats"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    LOG="${OUT_DIR}/${SAMPLE}.bowtie2.log"

    # PAIRED-END FIX: Replaced -U with -1 and -2
    bowtie2 -x ${BT2_IDX} \
        -1 ${R1} -2 ${R2} \
        --threads ${THR} \
        --no-unal \
    2> ${LOG} \
    | samtools view -b -q 30  \
    | samtools sort -@ 1 -o ${BAM}

    # Index and idxstats
    samtools index ${BAM}
    samtools idxstats ${BAM} > ${IDXSTATS}

    echo "Done: ${SAMPLE}"
}
export -f map_sample

# Run in parallel
if command -v parallel &> /dev/null; then
    cat ${TASK_FILE} | parallel -j ${MAX_JOBS} \
        map_sample "{}" "${BT2_INDEX}" "${MAPPED_DIR}" "${THREADS}"
else
    # xargs fallback 
    cat ${TASK_FILE} | xargs -P ${MAX_JOBS} -I {} bash -c "map_sample \"{}\" \"${BT2_INDEX}\" \"${MAPPED_DIR}\" \"${THREADS}\""
fi

# Build manifest of completed mappings
MAP_MANIFEST="${MAPPED_DIR}/map_manifest.tsv"
echo -e "sample_id\tidxstats_path" > ${MAP_MANIFEST}
for IDX in ${MAPPED_DIR}/*.idxstats; do
    SAMPLE=$(basename ${IDX} .idxstats)
    echo -e "${SAMPLE}\t${IDX}" >> ${MAP_MANIFEST}
done

# Summary
DONE=$(ls ${MAPPED_DIR}/*.idxstats 2>/dev/null | wc -l)

echo ""
echo "========================================"
echo "Paired-End Mapping complete."
echo "idxstats files: ${DONE}"
echo "Output directory: ${MAPPED_DIR}"
