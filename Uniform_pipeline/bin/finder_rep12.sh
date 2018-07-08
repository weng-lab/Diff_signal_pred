#!/usr/bin/env bash
cd /home/adam/dnase/encode/finder_data/
input1=$1
input2=$2
filetype=$3
IN_BAM_PREFIX=${input1%.bam}
IN_BAM_FILE=${input1}
OFPREFIX=${input2}
mm_chrom="/home/public/software/chip-seq-pipeline/dnanexus/shell/resources/home/dnanexus/mm10.chrom.sizes"
mm_ref="/home/public/igenome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa.fai"
FINAL_SAM_FILE="${OFPREFIX}.sam"
FINAL_BAM_FILE="${OFPREFIX}.bam"
#cp $IN_BAM_FILE $FINAL_BAM_FILE
samtools index ${FINAL_BAM_FILE}
PR_PREFIX="${OFPREFIX}.filt.nodup"
PR1_SAM_FILE="${PR_PREFIX}.pr1.sam"
PR2_SAM_FILE="${PR_PREFIX}.pr2.sam"

if [ "$filetype" = 'ChIP-seq_Control' ];then
    echo
else
    nlines=$(samtools view -c ${FINAL_BAM_FILE})
    nlines=$(( (nlines + 1) / 2 ))
    samtools view ${FINAL_BAM_FILE} | shuf > ${FINAL_SAM_FILE}
    head -n${nlines} ${FINAL_SAM_FILE} >${PR1_SAM_FILE}
    tail -n${nlines} ${FINAL_SAM_FILE} >${PR2_SAM_FILE}
    samtools view -bt $mm_ref ${PR1_SAM_FILE} > ${PR_PREFIX}.pr1.bam
    samtools view -bt $mm_ref  ${PR2_SAM_FILE} > ${PR_PREFIX}.pr2.bam
    samtools sort ${PR_PREFIX}.pr1.bam ${PR_PREFIX}.pr1
    samtools sort ${PR_PREFIX}.pr2.bam ${PR_PREFIX}.pr2
    samtools index ${PR_PREFIX}.pr1.bam
    samtools index ${PR_PREFIX}.pr2.bam
#    rm ${FINAL_SAM_FILE}
#    rm ${PR1_SAM_FILE}
#    rm ${PR2_SAM_FILE}
fi