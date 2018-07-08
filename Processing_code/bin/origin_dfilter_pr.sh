#!/usr/bin/env bash
inputs=$1
outpre=$2
tissue=$3

enhancer_dir="/Data/adam/dnase/enhancer/tissue_enhancer/"

awk '{printf "%s\t%d\n",$0,NR}' ${inputs} > ${outpre}.adjust.bed
peaks="${outpre}.adjust.bed"
tissue_enhancer=${enhancer_dir}${tissue}"_enhancer.txt"
negative_enhancer=${enhancer_dir}${tissue}"_negative_enhancer.txt"
total_peaks=$(cat $peaks | wc -l)
total_positive=$( cat $tissue_enhancer |wc -l  )
total_negative=$( cat  $negative_enhancer | wc -l )

output1="${outpre}.positive.bed"
output2="${outpre}.negative.bed"
output3="${outpre}.false.positive.bed"
output4="${outpre}.true.negative.bed"
output="${outpre}.pr.bed"

# true positive/false positive peaks
cut -f1,2,3,4 $peaks |intersectBed -wo -a $tissue_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,4 $peaks |intersectBed -wo -a $negative_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,4 $peaks | intersectBed -v -a $tissue_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3,4 $peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4