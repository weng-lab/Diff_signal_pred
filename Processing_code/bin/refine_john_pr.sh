#!/usr/bin/env bash
inputs=$1
outpre=$2
tissue=$3
k27bw=$4
dnasebw=$5
width=250
enhancer_dir="/Data/adam/dnase/enhancer/tissue_enhancer/"

out_tmp=${outpre}.tmp
awk -v width=${width} '{printf "%s\t%s\t%s\t%s\t%s\n",$1,$2-width,$2+width-1,$2,NR}' ${inputs} |\
  awk '{ if ($2>1){print $0}}'|sort -k 1,1 -k 2n,2n | mergeBed -i stdin -c 4,5 -o mean,min | sort -k 5n,5n|\
awk -v width=${width} '{printf("%s\t%.f\t%.f\t%s\t%s\n"),$1,$4-width,$4+width-1,$5,NR}' |
 awk '{ if ($2>1){print $0}}' > $out_tmp

bigWigAverageOverBed -bedOut=${outpre}.dnase_only.tmp.adjust.bed ${dnasebw} ${out_tmp} ${outpre}.dnase_only.tmp.adjust.txt
sort -k 6gr,6gr ${outpre}.dnase_only.tmp.adjust.bed | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2-750,$3+750,NR}' \
 > ${outpre}.dnase_only.2k.bed
awk '{printf "%s\t%s\t%s\t%s\n",$1,$2+750,$3-750,NR}' ${outpre}.dnase_only.2k.bed |\
 sort -k1,1 -k2g,2g > ${outpre}.adjust.bed

rm ${outpre}.dnase_only.tmp.adjust.bed ${outpre}.dnase_only.tmp.adjust.txt ${out_tmp}

bigWigAverageOverBed -bedOut=${outpre}.combine.tmp.adjust.bed ${k27bw} ${outpre}.dnase_only.2k.bed ${outpre}.combine.tmp.adjust.txt

sort -k 5gr,5gr ${outpre}.combine.tmp.adjust.bed | awk '{printf "%s\t%s\t%s\t%s\t%s\t%.0f\n",$1,$2+750,$3-750,$4,NR,(NR+$4)/2}' \
| sort -k1,1 -k2g,2g > ${outpre}.combine.adjust.bed
rm ${outpre}.combine.tmp.adjust.bed ${outpre}.combine.tmp.adjust.txt ${outpre}.dnase_only.2k.bed


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
cut -f1,2,3,6 $peaks |intersectBed -wo -a $tissue_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,6 $peaks |intersectBed -wo -a $negative_enhancer -b stdin | mergeBed -i stdin -c 7 -o min | cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,6 $peaks | intersectBed -v -a $tissue_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3,6 $peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4

# call pr for origin peaks
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
cut -f1,2,3,5 $peaks | intersectBed -v -a $negative_enhancer -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4