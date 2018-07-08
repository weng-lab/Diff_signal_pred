#!/bin/bash
input1=$1
DATASET_PREFIX=${input1%.rep*}
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
# ========================
# Create pooled datasets
# =======================
REP1_TA_FILE="${DATASET_PREFIX}.rep1.tagAlign.gz"
REP2_TA_FILE="${DATASET_PREFIX}.rep2.tagAlign.gz"
POOLED_TA_FILE="${DATASET_PREFIX}.rep0.tagAlign.gz"
#zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -c > ${POOLED_TA_FILE}

bedtools bedtobam -i ${DATASET_PREFIX}.rep0.tagAlign.gz -g ${mm_chrom} > ${DATASET_PREFIX}.rep0.bam
samtools sort ${DATASET_PREFIX}.rep0.bam ${DATASET_PREFIX}.rep0
samtools index ${DATASET_PREFIX}.rep0.bam

