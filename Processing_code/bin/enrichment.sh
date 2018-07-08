#!/usr/bin/env bash
peaks=$1
bedfile=$2
soft=$3
tissue=$4
width=$5
export LC_ALL=C
head -n2000 ${peaks}| sort -k1,1 -k2,2g | intersectBed -a stdin -b ${bedfile} -c -sorted | \
awk -v soft=$soft -v tissue=$tissue -v width=$width '{printf "%s\t%s\t%s\t%s\t%s\n",$6,$5,soft,tissue,width}' \
>>/Data/Peakcalling/pr_calling/select_peak/tmp/${soft}.dnase.${tissue}.${width}.txt

