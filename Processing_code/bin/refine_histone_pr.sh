#!/usr/bin/env bash
inputs=$1
outpre=$2
tissue=$3
width=$4
dnasepeak=$5
enhancer_dir="/Data/adam/dnase/enhancer/tissue_enhancer/"
awk -v width=${width} '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2-width,$2+width-1,$2,NR}' ${inputs} |  awk '{ if ($2>1){print $0}}'|\
sort -k 1,1 -k 2n,2n | mergeBed -i stdin -c 4,5 -o mean,min | sort -k 5n,5n|\
awk -v width=${width} '{printf("%s\t%s\t%s\t%s\t%s\n"),$1,$4-width,$4+width-1,$5,NR}' |
 awk '{ if ($2>1){print $0}}'> ${outpre}.adjust.bed

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
cut -f1,2,3,5 $peaks |intersectBed -wo -a $tissue_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,5 $peaks |intersectBed -wo -a $negative_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,5 $peaks | intersectBed -v -a $tissue_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3,5 $peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4

combine_peaks="${outpre}.combine.adjust.bed"
cut -f1,2,3,5 $peaks | intersectBed -wa -wb -a stdin -b ${dnasepeak} | sort -k 1,1 -k 2n,2n |\
mergeBed -i stdin -c 4,9 -o min,min| sort -k 4g,4g| awk '{printf "%s\t%s\t%s\t%.0f\n",$1,$2,$3,(NR+$5)/2}'|\
 sort -k 4g,4g | awk '{printf "%s\t%s\n",$0,NR}'>${combine_peaks}

# true positive/false positive peaks
cut -f1,2,3,5 $combine_peaks |intersectBed -wo -a $tissue_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,5 $combine_peaks |intersectBed -wo -a $negative_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,5 $combine_peaks | intersectBed -v -a $tissue_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3,5 $combine_peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > "${outpre}.combine.pr.bed"
rm $output1 $output2 $output3 $output4