#!/usr/bin/env bash
PEAK=$1
MARK=$2
TISSUE=$3
SOFT=$4

RESULT="/Data/adam/dnase/enrich_results.txt"
final_temp="/home/adam/dnase/encode/temp/${MARK}.${TISSUE}.rep0.tagAlign.gz"
seq_total=$(intersectBed -a ${final_temp} -b $PEAK -wa -sorted |awk '{sum += $3-$2};END {print sum}')
peak_total=$(awk '{sum += $3-$2};END {print sum}' $PEAK)
enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')
echo "${enrich} ${seq_total}    ${peak_total}   ${SOFT} ${MARK}" >> $RESULT
