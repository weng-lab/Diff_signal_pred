#!/usr/bin/env bash
inputs=$1
outpre=$2
tissue=$3
types=$4

promoter_dir="/Data/adam/dnase/promoter/mouse_10fpkm_promoter/tissue_promoter/"


#if [ "${types}" = "histone" ]; then
##    echo "running ${inputs}"
#    width=1000
#    awk -v width=${width} '{printf "%s\t%d\t%d\t%d\t%d\n",$1,$2-width,$2+width-1,$2,NR}' ${inputs} |  awk '{if ($2>1){print $0}}'|\
#    sort -k 1,1V -k 2,2n | mergeBed -i stdin -c 4,5 -o mean,min | sort -k 5n,5n|\
#    awk -v width=${width} '{printf("%s\t%d\t%d\t%d\t%d\n"),$1,$4-width,$4+width-1,$5,NR}' |\
#    awk '{if ($2>1){print $0}}'> ${outpre}.adjust.bed
##    head -n1 ${outpre}.adjust.bed
#else
#    cat ${inputs} > ${outpre}.adjust.bed
#
#fi

#awk '{printf "%s\t%d\n",$0,NR}' ${inputs} > ${outpre}.adjust.bed
#peaks="${outpre}.adjust.bed"
peaks=$inputs

tissue_promoter=${promoter_dir}${tissue}"_promoter.txt"
negative_promoter=${promoter_dir}${tissue}"_negative_promoter.txt"
total_peaks=$(cat $peaks | wc -l)
total_positive=$( cat $tissue_promoter |wc -l  )
total_negative=$( cat  $negative_promoter | wc -l )

output1="${outpre}.positive.bed"
output2="${outpre}.negative.bed"
output3="${outpre}.false.positive.bed"
output4="${outpre}.true.negative.bed"
output="${outpre}.pr.bed"

# true positive/false positive peaks
cut -f1,2,3,5 $peaks |intersectBed -wo -a $tissue_promoter -b stdin | mergeBed -i stdin -c 10 -o min | cut -f4 | awk '{printf "%d\t%s\n",$0,1}' > $output1
cut -f1,2,3,5 $peaks |intersectBed -wo -a $negative_promoter -b stdin | mergeBed -i stdin -c 10 -o min | cut -f4 | awk '{printf "%d\t%s\n",$0,0}' > $output2
# false negative/true negative peaks
cut -f1,2,3,5 $peaks | intersectBed -v -a $tissue_promoter -b stdin -s| awk -v totals=$total_peaks '{printf "%d\t%s\n",totals,1}' > $output3
#cut -f1,2,3,5 $peaks | intersectBed -v -a $negative_promoter -b stdin | awk -v totals=$total_peaks '{printf "%s\t%s\n",totals,0}' > $output4
awk -v totals=$total_peaks 'BEGIN{printf "%s\t%s\n",totals,0}' > $output4
cat $output1 $output2 $output3  $output4 |sort -k 1g,1g  > $output
rm $output1 $output2 $output3  $output4