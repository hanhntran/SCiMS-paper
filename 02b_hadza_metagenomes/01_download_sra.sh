#!/bin/bash cd $(pwd)
module load anaconda
conda activate SRA-toolkit

# ----- CONFIGURE THESE PATHS -----
SRR_LIST="srr_list.txt"
OUT_DIR="downloaded_fastqs"
# ----- END CONFIG -----

mkdir -p ${OUT_DIR}

# Check if list exists
if [ ! -f "$SRR_LIST" ]; then
    echo "Error: Cannot find $SRR_LIST"
    exit 1
fi

echo "Starting SRA download pipeline..."

# Read the file line by line
while IFS= read -r SRR; do
    # Skip empty lines
    if [ -z "$SRR" ]; then continue; fi
    
    echo "========================================"
    echo "Processing: $SRR"
    
    # 1. Download the compressed .sra file safely
    # prefetch is smart; if the file is already downloaded, it skips it
    echo "Prefetching $SRR..."
    prefetch ${SRR}
    
    echo "Extracting FASTQ for $SRR..."
    fasterq-dump ${SRR} \
        --split-files \
        --outdir ${OUT_DIR} \
        --threads ${SLURM_NTASKS} \
        --progress
        
    echo "Cleaning up SRA cache for $SRR..."
    rm -rf ${SRR}
    
    echo "$SRR completed successfully."

done < "$SRR_LIST"

echo "========================================"
echo "All downloads complete! Files are in ${OUT_DIR}"
