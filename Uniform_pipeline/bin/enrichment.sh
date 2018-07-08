#!/usr/bin/env bash
cd /home/adam/dnase/encode/
MARK=$1
TISSUE=$2

final_tagalign="/home/adam/dnase/encode/data/${MARK}/${TISSUE}/${MARK}.${TISSUE}.rep0.tagAlign.gz"
final_bam = "/home/adam/dnase/encode/data/${MARK}/${TISSUE}/${MARK}.${TISSUE}.rep0.bam"
final_temp="/home/adam/dnase/encode/temp/${MARK}.${TISSUE}.rep0.tagAlign.gz"

mkdir ../summit_peak/


#zcat ${final_tagalign} | sort -k1,1 -k2,2n | gzip -c > ${final_temp}
for SOFT in music; do
    for peaktype in broad punctuate; do
        final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.final.${peaktype}.bed.gz"
        seq_total=$(zcat $final_file | sort -k1,1 -k2,2n | intersectBed -a stdin -b $final_temp -wb -sorted| awk '{sum += $3-$2};END {print sum}')
        peak_total=$(zcat $final_file| awk '{sum += $3-$2};END {print sum}')
        enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')
        echo "${enrich} ${seq_total}    ${peak_total}   ${SOFT}.${peaktype} ${MARK}" >> ./enrichment/enrichment.txt
    done
done


for SOFT in macs2; do
    for peaktype in "broad" "gapped" "narrow"; do
        final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.final.${peaktype}Peak.gz"
        seq_total=$(zcat $final_file | sort -k1,1 -k2,2n | intersectBed -a stdin -b $final_temp -wb -sorted| awk '{sum += $3-$2};END {print sum}')
        peak_total=$(zcat $final_file| awk '{sum += $3-$2};END {print sum}')
        enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')

        echo "${enrich} ${seq_total}    ${peak_total}   ${SOFT}.${peaktype}  ${MARK}" >> ./enrichment/enrichment.txt
    done
done


#for SOFT in homer hotspot rseg finder mosaics dfilter;do
for SOFT in homer hotspot dfilter;do
    final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.final.bed.gz"
    seq_total=$(zcat $final_file | sort -k1,1 -k2,2n | intersectBed -a stdin -b $final_temp -wb -sorted| awk '{sum += $3-$2};END {print sum}')
    peak_total=$(zcat $final_file| awk '{sum += $3-$2};END {print sum}')
    enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')
#    echo $seq_total $peak_total $enrich
    echo "${enrich} ${seq_total}    ${peak_total}   ${SOFT}    ${MARK}" >> ./enrichment/enrichment.txt
done



#for SOFT in bcp;do
#    final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.final_results_HM.bed.gz"
#    seq_total=$(zcat $final_file | sort -k1,1 -k2,2n | intersectBed -a stdin -b $final_temp -wb -sorted| awk '{sum += $3-$2};END {print sum}')
#    peak_total=$(zcat $final_file| awk '{sum += $3-$2};END {print sum}')
#    enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')
##    echo $seq_total $peak_total $enrich
#    echo "${enrich} ${seq_total}    ${peak_total}   ${SOFT}    ${MARK}" >> ./enrichment/enrichment.txt
#done

final_tagalign="/home/adam/dnase/encode/data/${MARK}/${TISSUE}/${MARK}.${TISSUE}.rep1.tagAlign.gz"
final_temp="/home/adam/dnase/encode/temp/${MARK}.${TISSUE}.rep1.tagAlign.gz"
#zcat ${final_tagalign} | sort -k1,1 -k2,2n | gzip -c > ${final_temp}
for SOFT in fseq; do
#    final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.final.narrowPeak.gz"
    final_file="/home/adam/dnase/encode/data_out/${SOFT}/${MARK}/${TISSUE}/${MARK}.${TISSUE}.rep1.narrowPeak.gz"
    seq_total=$(zcat $final_file | sort -k1,1 -k2,2g | intersectBed -b stdin -a $final_temp -sorted -wa | awk '{sum += $3-$2};END {print sum}')
    peak_total=$(zcat $final_file| awk '{sum += $3-$2};END {print sum}')
    enrich=$(awk 'BEGIN{printf '${seq_total}'/'${peak_total}'}')
    echo "${enrich}    ${seq_total}   ${peak_total}    ${SOFT} ${MARK}" >> ./enrichment/enrichment.txt
done