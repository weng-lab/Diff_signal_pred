#!/bin/bash

peak_file=$1
enhancer_file=$2
prefix=$3
outdir=$4

mkdir -p ${outdir}
awk '{if ($4==1){print $0}}' $enhancer_file | sort -k1,1 -k2g,2g |cut -f1-3 > ${outdir}/pos_enhancer.txt
awk '{if ($4==0){print $0}}' $enhancer_file | sort -k1,1 -k2g,2g |cut -f1-3 > ${outdir}/neg_enhancer.txt

output1="${outdir}/${prefix}prauc_positive.bed"
output2="${outdir}/${prefix}prauc_negative.bed"
output3="${outdir}/${prefix}false.positive.bed"
output4="${outdir}/${prefix}true.negative.bed"
output="${outdir}/${prefix}prauc.txt"

total_peaks=$(cat $peak_file | wc -l)
total_positive=$( cat ${outdir}/pos_enhancer.txt | wc -l )
total_negative=$( cat ${outdir}/neg_enhancer.txt | wc -l )

cut -f1,2,3 $peak_file | awk '{printf "%s\t%d\n",$0,NR}' \
|intersectBed -wo -a ${outdir}/pos_enhancer.txt -b stdin | cut -f 1-3,7 |sort -k 1,1 -k 2g,2g -k 4g,4g > ${outdir}/temp1.txt
python ${DIFF_PRED}/pyscript/get_best_peak.py ${outdir}/temp1.txt $output1 1

cut -f1,2,3 $peak_file  | awk '{printf "%s\t%d\n",$0,NR}' \
|intersectBed -wo -a ${outdir}/neg_enhancer.txt -b stdin |  cut -f 1-3,7 |sort -k 1,1 -k 2g,2g -k 4g,4g > ${outdir}/temp2.txt
python ${DIFF_PRED}/pyscript/get_best_peak.py ${outdir}/temp2.txt $output2 0

# false negative/true negative peaks
cut -f1,2,3 $peak_file | awk '{printf "%s\t%d\n",$0,NR}' | intersectBed -v -a ${outdir}/pos_enhancer.txt -b stdin|sort| uniq| awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,1}' > $output3
cut -f1,2,3 $peak_file | awk '{printf "%s\t%d\n",$0,NR}' | intersectBed -v -a ${outdir}/neg_enhancer.txt -b stdin |sort| uniq | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4


cat $output1 $output2 $output3 $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3 $output4 ${outdir}/temp1.txt ${outdir}/temp2.txt