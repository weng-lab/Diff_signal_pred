#!/bin/bash
# ===================
# Create tagAlign file
# ===================
input1=$1
input2=$2
FINAL_BAM_PREFIX=${input1%.bam}
FINAL_BAM_FILE=${input1}
OFPREFIX=${input2}
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
# Create SE tagAlign file
FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"
bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${FINAL_TA_FILE}

#get bam format and index
#bedtools bedtobam -i ${FINAL_TA_FILE} -g ${mm_chrom} > ${OFPREFIX}.bam
#samtools sort ${OFPREFIX}.bam ${OFPREFIX}
#samtools index ${OFPREFIX}.bam
# =================================
# Subsample tagAlign file
# ================================
NREADS=5000000
SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$((NREADS / 1000000)).tagAlign.gz"
#zcat ${FINAL_TA_FILE} | grep -v "chrM" | shuf -n ${NREADS} | gzip -c > ${SUBSAMPLED_TA_FILE}


# ========================
# Create pseudoReplicates
# =======================
PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_TA_FILE="${PR_PREFIX}.pr1.tagAlign.gz"
PR2_TA_FILE="${PR_PREFIX}.pr2.tagAlign.gz"
# Get total number of read pairs

nlines=$(zcat ${FINAL_TA_FILE} | wc -l)
nlines=$(( (nlines + 1) / 2 ))

# Shuffle and split BED file into 2 equal parts
zcat ${FINAL_TA_FILE} | shuf | split -d -l ${nlines} - ${PR_PREFIX}

# Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01
# Convert reads into standard tagAlign file
gzip -c "${PR_PREFIX}00"> ${PR1_TA_FILE}
rm "${PR_PREFIX}00"
gzip -c "${PR_PREFIX}01" > ${PR2_TA_FILE}
rm "${PR_PREFIX}01"

#get bam format and index
#bedtools bedtobam -i ${PR1_TA_FILE} -g ${mm_chrom} > ${PR_PREFIX}.pr1.bam
#samtools sort  ${PR_PREFIX}.pr1.bam ${PR_PREFIX}.pr1
#samtools index ${PR_PREFIX}.pr1.bam
#bedtools bedtobam -i ${PR2_TA_FILE} -g ${mm_chrom} > ${PR_PREFIX}.pr2.bam
#samtools sort  ${PR_PREFIX}.pr2.bam ${PR_PREFIX}.pr2
#samtools index ${PR_PREFIX}.pr2.bam

# =============================
# run Spp for fragment length
# =============================

SPP="/usr/local/src/phantompeakqualtools/run_spp_nodups.R"
CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.spp.ccscores"
#Rscript ${SPP} -c=${SUBSAMPLED_TA_FILE} -rf -out=${CC_SCORES_FILE} -p=4 -filtchr=chrM
#sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
#mv temp ${CC_SCORES_FILE}

# =============================
# get bed file from tagAlign
# =============================
# gunzip -c ${PR_PREFIX}.pr1.tagAlign.gz | sort > ${PR_PREFIX}.pr1.bed
# gunzip -c ${PR_PREFIX}.pr2.tagAlign.gz | sort > ${PR_PREFIX}.pr2.bed
# gunzip -c ${OFPREFIX}.tagAlign.gz | sort  > ${OFPREFIX}.bed