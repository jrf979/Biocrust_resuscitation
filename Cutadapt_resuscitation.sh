#!/bin/bash
#SBATCH --job-name=cutadapt_resuscitation
#SBATCH --output=cutadapt_resuscitation_%j.log
#SBATCH --error=cutadapt_resuscitation_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G

# =============================================================================
# Cutadapt_resuscitation.sh
# Author: Raul Roman
# =============================================================================

echo "Starting Cutadapt job at $(date)"
echo "Running on node: $(hostname)"

module load anaconda
conda activate cutadapt

cd /storage/work/jxr6215/16S_DADA2_20250923/sequences/ || exit 1

mkdir -p cutadapt_trimmed

for R1 in *_R1_001.fastq.gz; do
    S=${R1%_R1_001.fastq.gz}
    R2=${S}_R2_001.fastq.gz

    if [ ! -e "$R2" ]; then
        echo "Missing $R2 for $S" >&2
        continue
    fi

    echo "Processing sample: $S"

    cutadapt \
        -j 8 \
        -e 0.10 \
        -O 10 \
        -g ^GTGYCAGCMGCCGCGGTAA \
        -G ^GGACTACNVGGGTWTCTAAT \
        -a ATTAGAWACCCBNGTAGTCCX \
        -A TTACCGCGGCKGCTGRCACX \
        --pair-filter=any \
        -m 50 \
        -o cutadapt_trimmed/${S}_R1_001.trim.fastq.gz \
        -p cutadapt_trimmed/${S}_R2_001.trim.fastq.gz \
        "$R1" "$R2" \
        2>&1 | tee cutadapt_trimmed/${S}.cutadapt.log
done

echo "Cutadapt job finished at $(date)"
