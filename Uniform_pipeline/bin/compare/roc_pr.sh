#!/usr/bin/env bash
tissue_enhancer=$1
negative_enhancer=$2
soft=$3
mark=$4
tissue=$5
peaks="/Data/adam/dnase/enrich_all_merge_bed/${soft}.${mark}.${tissue}.bed"
output1="/Data/adam/dnase/roc_pr_value/${soft}.${mark}.${tissue}.positive.bed"
output2="/Data/adam/dnase/roc_pr_value/${soft}.${mark}.${tissue}.negative.bed"
output="/Data/adam/dnase/roc_pr_value/${soft}.${mark}.${tissue}.bed"





#echo "start"
## tissue enhancer must be 3 columns
#total_peaks=$( wc -l $peaks)
#total_positive=$( wc -l $tissue_enhancer )
#total_negative=$( wc -l $negative_enhancer )
#intersectBed -wo -a $tissue_enhancer -b $peaks | sort -k 1,1 -k 2g,2g -k 8gr,8gr -k 7g,7g | awk '{a[$3]=$0}END{for(i in a)print a[1]}'| cut -f7 | awk -v totalpositive=$total_positive '{printf "%s\t%s\t%s\t%s\t%s\n",$0,NR,NR/totalpositive,0,0}' > $output1
#echo "hello"
#intersectBed -wo -a $negative_enhancer -b $peaks | sort -k 1,1 -k 2g,2g -k 8gr,8gr -k 7g,7g | awk '{a[$3]=$0}END{for(i in a)print a[1]}'| cut -f7 | awk -v totalnegative=$total_negative '{printf "%s\t%s\t%s\t%s\t%s\n",$0,0,0,NR,NR/totalnegative}' >  $output2
#
#cat $output1 $output2 > $output
#rm $output1 $output2
#sort -k 1g,1g $output >$output