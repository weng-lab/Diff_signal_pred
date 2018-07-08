#!/usr/bin/env bash
inputbam=$1
outprefix=$2
samtools depth ${inputbam} | perl -ne 'BEGIN{ print "track type=print wiggle_0 name=fileName description=fileName\n"}; ($c, $start, $depth) = split;if ($c ne $lastC) { print "variableStep chrom=$c span=2\n"; };$lastC=$c;next unless $. % 2 ==0;print "$start\t$depth\n" unless $depth<3' > "/Data/adam/dnase/bigwig/${outprefix}wig"
wigToBigWig "/Data/adam/dnase/bigwig/${outprefix}wig" "/home/public/software/gem/mm10.chrom.sizes" "/Data/adam/dnase/bigwig/${outprefix}bw"