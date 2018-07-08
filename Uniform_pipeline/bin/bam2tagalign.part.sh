#!/bin/bash
# ===================
# Create tagAlign file
# ===================
input1=$1
input2=$2
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
FINAL_BAM_PREFIX=${input1%.bam}
FINAL_BAM_FILE=${input1}
OFPREFIX=${input2}
# Create SE tagAlign file
FINAL_TA_FILE="${OFPREFIX}.tagAlign.gz"
#bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${FINAL_TA_FILE}

bedtools bedtobam -i ${FINAL_TA_FILE} -g ${mm_chrom} > ${OFPREFIX}.bam
samtools sort ${OFPREFIX}.bam ${OFPREFIX}
samtools index ${OFPREFIX}.bam

# =================================
# Subsample tagAlign file
# ================================
NREADS=5000000
SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$((NREADS / 1000000)).tagAlign.gz"
#zcat ${FINAL_TA_FILE} | grep -v "chrM" | shuf -n ${NREADS} | gzip -c > ${SUBSAMPLED_TA_FILE}




# =============================
# run Spp for fragment length
# =============================

#SPP="/usr/local/src/phantompeakqualtools/run_spp_nodups.R"
#CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.spp.ccscores"
#Rscript ${SPP} -c=${SUBSAMPLED_TA_FILE} -rf -out=${CC_SCORES_FILE} -p=8 -filtchr=chrM
#sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
#mv temp ${CC_SCORES_FILE}