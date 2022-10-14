#!/usr/bin/env bash

set -u

machine=`hostname -s|cut -d . -f1`

echo 
echo "Running script:$0 on machine:$machine on date:$(date)"

##########################################################################
##########################################################################
###
### This script gets and prepares a gtf file with the piRNA coordinates
###
##########################################################################
##########################################################################

myDir=$HOME/Projects/mouse-piRNA-variation
liftOverEx=$HOME/bin/liftOver
bedtoolsEx=/usr/local/bin/bedtools

### Get the chain file for liftOver (needed for conversion of mouse genome coordinates)
if [ ! -f $myDir/data/forLiftOver/mm9ToMm10.over.chain.gz ]; then 
    echo "Necessary chain file not found so getting it from UCSC"
    mkdir -p $myDir/data/forLiftOver/
    wget -P $myDir/data/forLiftOver http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
fi

### Get the piRNA cluster coordinates from GEO (they are in mm9)
if [ ! -f $myDir/data/GSE45049_102-assembled.piRNA.transcripts.bed ]; then 
    echo "piRNA transcript file not found so getting it from GEO"
    wget -P $myDir/data/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45049/suppl/GSE45049_102-assembled.piRNA.transcripts.bed.gz
    gunzip $myDir/data/GSE45049_102-assembled.piRNA.transcripts.bed.gz
fi

# Convert the coordinates of the piRNA clusters from mm9 to mm10
if [ ! -f $myDir/data/GSE45049_102-assembled.piRNA.transcripts.mm10.bed ]; then
    echo "Convert coordinates of piRNA transcripts from mm9 to mm10"
    $liftOverEx -bedPlus=6 $myDir/data/GSE45049_102-assembled.piRNA.transcripts.bed $myDir/data/forLiftOver/mm9ToMm10.over.chain.gz $myDir/data/GSE45049_102-assembled.piRNA.transcripts.mm10.bed $myDir/data/liftOverMm9Mm10Unmapped
fi

# Merge overlapping piRNA transcripts on the same strand to generate piRNA clusters
if [ -f $myDir/data/piRNA.gtf ];
then
    echo "piRNA.gtf file found. Exiting without generating a new one."
else
    echo "Merge overlapping piRNA transcripts and convert to gtf"
    cat $myDir/data/GSE45049_102-assembled.piRNA.transcripts.mm10.bed |sed 's/chr//' | sort -k1,1n -k2,2n | $bedtoolsEx merge -s -c 6 -o distinct | awk 'BEGIN{OFS="\t";n=1} {print $1,"LiMolCell2013\tgene",($2+1),$3,".",$4,".\tgene_id \"piC"n"\"";n++}' > $myDir/data/piRNA.gtf
fi


