#!/usr/bin/env bash
inputs=$1
outpre=$2
tissue=$3
k27bw=$4

enhancer_dir="/Data/adam/dnase/enhancer/tissue_enhancer/"
input_2k=${inputs}.2k
awk '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2-850,$3+850,$4,$5}' ${inputs} > ${input_2k}
bigWigAverageOverBed -bedOut=${outpre}.combine.tmp.adjust.bed ${k27bw} ${input_2k} ${outpre}.combine.tmp.adjust.txt
sort -k 6gr,6gr ${outpre}.combine.tmp.adjust.bed | awk '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.0f\n",$1,$2+850,$3-850,$4,$5,$6,NR,(NR+$5)/2}' \
 > ${outpre}.combine.adjust.bed
rm ${outpre}.combine.tmp.adjust.bed ${outpre}.combine.tmp.adjust.txt ${input_2k}
peaks="${outpre}.combine.adjust.bed"
tissue_enhancer=${enhancer_dir}${tissue}"_enhancer.txt"
negative_enhancer=${enhancer_dir}${tissue}"_negative_enhancer.txt"
total_peaks=$(cat $peaks | wc -l)
total_positive=$( cat $tissue_enhancer |wc -l  )
total_negative=$( cat  $negative_enhancer | wc -l )

output1="${outpre}.positive.bed"
output2="${outpre}.negative.bed"
output3="${outpre}.false.positive.bed"
output4="${outpre}.true.negative.bed"
output="${outpre}.combine.pr.bed"

# true positive/false positive peaks
cut -f1,2,3,8 $peaks |intersectBed -wo -a $tissue_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,8 $peaks |intersectBed -wo -a $negative_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,8 $peaks | intersectBed -v -a $tissue_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3,8 $peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4