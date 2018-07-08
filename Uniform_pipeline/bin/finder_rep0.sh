#!/bin/bash
input1=$1
DATASET_PREFIX=${input1%.rep*}

REP1_BAM_FILE="${DATASET_PREFIX}.rep1.bam"
REP2_BAM_FILE="${DATASET_PREFIX}.rep2.bam"
POOLED_BAM_FILE="${DATASET_PREFIX}.rep0.bam"
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
samtools merge -f ${POOLED_BAM_FILE} ${REP1_BAM_FILE} ${REP2_BAM_FILE}
samtools sort ${POOLED_BAM_FILE} ${DATASET_PREFIX}.rep0
samtools index ${POOLED_BAM_FILE}

REP1_PR1_BAM_FILE="${DATASET_PREFIX}.rep1.filt.nodup.pr1.bam"
REP1_PR2_BAM_FILE="${DATASET_PREFIX}.rep1.filt.nodup.pr2.bam"
REP2_PR1_BAM_FILE="${DATASET_PREFIX}.rep2.filt.nodup.pr1.bam"
REP2_PR2_BAM_FILE="${DATASET_PREFIX}.rep2.filt.nodup.pr2.bam"
PPR1_BAM_FILE="${DATASET_PREFIX}.rep0.filt.nodup.pr1.bam"
PPR2_BAM_FILE="${DATASET_PREFIX}.rep0.filt.nodup.pr2.bam"

if [ "$filetype" = 'ChIP-seq_Control' ];then
    echo
else
    samtools merge -f ${PPR1_BAM_FILE} ${REP1_PR1_BAM_FILE} ${REP2_PR1_BAM_FILE}
    samtools sort ${PPR1_BAM_FILE} ${DATASET_PREFIX}.rep0.filt.nodup.pr1
    samtools index ${DATASET_PREFIX}.rep0.filt.nodup.pr1.bam

    samtools merge -f ${PPR2_BAM_FILE} ${REP1_PR2_BAM_FILE} ${REP2_PR2_BAM_FILE}
    samtools sort ${PPR2_BAM_FILE} ${DATASET_PREFIX}.rep0.filt.nodup.pr2
    samtools index ${DATASET_PREFIX}.rep0.filt.nodup.pr2.bam

fi
