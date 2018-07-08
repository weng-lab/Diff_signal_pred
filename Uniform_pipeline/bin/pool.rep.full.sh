#!/bin/bash
input1=$1
DATASET_PREFIX=${input1%.rep*}
# ========================
# Create pooled datasets
# =======================
REP1_TA_FILE="${DATASET_PREFIX}.rep1.tagAlign.gz"
REP2_TA_FILE="${DATASET_PREFIX}.rep2.tagAlign.gz"
POOLED_TA_FILE="${DATASET_PREFIX}.rep0.tagAlign.gz"
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
#zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -c > ${POOLED_TA_FILE}

#get bam from tagAlign
#bedtools bedtobam -i ${DATASET_PREFIX}.rep0.tagAlign.gz -g ${mm_chrom} > ${DATASET_PREFIX}.rep0.bam
#samtools sort ${DATASET_PREFIX}.rep0.bam ${DATASET_PREFIX}.rep0
#samtools index ${DATASET_PREFIX}.rep0.bam

# ========================
# Create pooled pseudoreplicates
# =======================
REP1_PR1_TA_FILE="${DATASET_PREFIX}.rep1.filt.nodup.pr1.tagAlign.gz"
REP1_PR2_TA_FILE="${DATASET_PREFIX}.rep1.filt.nodup.pr2.tagAlign.gz"
REP2_PR1_TA_FILE="${DATASET_PREFIX}.rep2.filt.nodup.pr1.tagAlign.gz"
REP2_PR2_TA_FILE="${DATASET_PREFIX}.rep2.filt.nodup.pr2.tagAlign.gz"
PPR1_TA_FILE="${DATASET_PREFIX}.rep0.filt.nodup.pr1.tagAlign.gz"
PPR2_TA_FILE="${DATASET_PREFIX}.rep0.filt.nodup.pr2.tagAlign.gz"

#zcat ${REP1_PR1_TA_FILE} ${REP2_PR1_TA_FILE} | gzip -c > ${PPR1_TA_FILE}
#zcat ${REP1_PR2_TA_FILE} ${REP2_PR2_TA_FILE} | gzip -c > ${PPR2_TA_FILE}

#bedtools bedtobam -i ${PPR1_TA_FILE} -g ${mm_chrom} > ${DATASET_PREFIX}.rep0.filt.nodup.pr1.bam
#samtools sort ${DATASET_PREFIX}.rep0.filt.nodup.pr1.bam ${DATASET_PREFIX}.rep0.filt.nodup.pr1
#samtools index ${DATASET_PREFIX}.rep0.filt.nodup.pr1.bam
#
#bedtools bedtobam -i ${PPR2_TA_FILE} -g ${mm_chrom} > ${DATASET_PREFIX}.rep0.filt.nodup.pr2.bam
#samtools sort ${DATASET_PREFIX}.rep0.filt.nodup.pr2.bam ${DATASET_PREFIX}.rep0.filt.nodup.pr2
#samtools index ${DATASET_PREFIX}.rep0.filt.nodup.pr2.bam


# =============================
# get bed file from tagAlign
# =============================
gunzip -c ${DATASET_PREFIX}.rep0.filt.nodup.pr1.tagAlign.gz | sort > ${DATASET_PREFIX}.rep0.filt.nodup.pr1.bed
gunzip -c ${DATASET_PREFIX}.rep0.filt.nodup.pr2.tagAlign.gz | sort > ${DATASET_PREFIX}.rep0.filt.nodup.pr2.bed
gunzip -c ${DATASET_PREFIX}.rep0.tagAlign.gz | sort > ${DATASET_PREFIX}.rep0.bed
