#!/usr/bin/env bash
mark=$1
soft=$2
rep=$3
file=$4
counts=$(zcat $file |wc -l | awk '{print $1}')
echo "${counts} ${mark} ${soft} ${rep}" >> ../rep_basic/counts_rep.txt
zcat $file| awk '{print $3-$2}'|sort -k1,1gr | awk -v mark=$mark -v soft=$soft -v rep=$rep 'NR%500==0{printf("%s\t%s\t%s\t%s\n"),mark,soft,rep,$1}' >>../rep_basic/len_rep.txt
