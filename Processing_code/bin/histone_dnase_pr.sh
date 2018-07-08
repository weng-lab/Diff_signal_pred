#!/usr/bin/env bash
histonefile=$1
dnasefile=$2
outpre=$3
mygenome="/Data/Peakcalling/pr_calling/scripts/mm10.chrom.sizes"
shufflenegative="${outpre}.randomdnase.bed"
output1="${outpre}.positive.bed"
output2="${outpre}.negative.bed"
output3="${outpre}.false.positive.bed"
output4="${outpre}.true.negative.bed"
output="${outpre}.dnase_hist.pr.bed"

total_peaks=$(cat $histonefile | wc -l)

shuffleBed -i $dnasefile -g $mygenome -noOverlapping | sort -k1,1 -k2n,2n > $shufflenegative
cut -f1,2,3,5 ${histonefile} | intersectBed -wo -a $dnasefile -b stdin | sort -k1,1 -k2n,2n |mergeBed -i stdin -c 9 -o min | cut -f4 |\
 awk '{printf "%s\t%s\n",$0,1}' > $output1
cut -f1,2,3,5 ${histonefile} | intersectBed -wo -a $shufflenegative -b stdin |sort -k1,1 -k2n,2n | mergeBed -i stdin -c 9 -o min |\
 cut -f4 | awk '{printf "%s\t%s\n",$0,0}' > $output2
cut -f1,2,3,5 ${histonefile} | intersectBed -v  -a $dnasefile -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3

cut -f1,2,3,5 ${histonefile} | intersectBed -v -a $shufflenegative -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4

cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4
