#!/usr/bin/env bash

set -u

machine=`hostname -s|cut -d . -f1`
echo 
echo "Running script:$0 on machine:$machine on date:$(date)"

### The repeat masker annotation is in 3-column bed format, downloaded from UCSC for mm10.

myDir=~/Projects/mouse-piRNA-variation # Project working directory
rmsk=$myDir/data/rmsk_mm10.bed
piRNA=$myDir/data/piRNA.gtf
intersectBedEx=/usr/local/bin/intersectBed # Bedtools v2.29.2
featureCountsEx=/opt/subread-2.0.1-Linux-x86_64/bin/featureCounts # featureCounts v2.0.1

mkdir -p $myDir/out/bowtie/norepeats # Make folder for bam files with alignments overlapping repeats removed
mkdir -p $myDir/out/featureCounts/piRNA/norepeats # Make folder for featureCounts files with alignments overlapping repeats removed
cd $myDir/out/bowtie/ || exit 1

echo `pwd`

for s in `ls sample*.bam|sed "s/.bam//"`; do
    echo
    echo "$(date) : Processing $s"
    if [ -f norepeats/$s\_norepeats.bam ]; then
	echo "Filtered bam file for $s exists. Ignoring."
    else		
	$intersectBedEx -v -abam $s.bam -b $rmsk > norepeats/$s\_norepeats.bam 
    fi
    if [ -f  $myDir/out/featureCounts/piRNA/norepeats/$s.s0.counts ]; then
	echo "Counts for filtered bam file for $s exists. Ignoring."
    else		
	$featureCountsEx -Q 1 -T 7 -s 0 -a $piRNA -F GTF -t gene -g gene_id -o $myDir/out/featureCounts/piRNA/norepeats/$s.s0.counts -O --minOverlap 18 $myDir/out/bowtie/norepeats/$s\_norepeats.bam 2> $myDir/out/featureCounts/piRNA/norepeats/$s.s0.elog	    
    fi
done

