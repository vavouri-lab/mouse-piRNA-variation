#!/usr/bin/env bash

set -u

machine=`hostname -s|cut -d . -f1`
echo 
echo "Running script:$0 on machine:$machine on date:$(date)"

###############################################################################
### This script does the basic small RNA processing for our mouse samples. ####
###############################################################################

#######################################################################
############## Set locations of input files ###########################
#######################################################################
refGenome=GRCm38 # Mouse genome version 
refGenomeF=~/Projects/data/Mus_musculus.$refGenome.dna.primary_assembly.fa # Mouse genome fasta file
datasets=(batch1 batch2 batch3) # batch1 and batch2 are ICR data, batch3 are inbred strain data
myDir=~/Projects/mouse-piRNA-variation # Project working directory
#######################################################################

############################################################
############## Set locations of programs ###################
############################################################
bowtieBuildEx=~/miniconda2/bin/bowtie-build # version 1.2.3
bowtieEx=~/miniconda2/bin/bowtie # version 1.2.3
cutadaptEx=~/.local/bin/cutadapt # version 2.10 with Python 3.6.9
fastqQFEx=/usr/local/bin/fastq_quality_filter # version FASTX Toolkit 0.0.14
samtoolsEx=/usr/local/bin/samtools # version 1.10
############################################################


function buildBowtieIndex { # If there is no bowtie index, generate it
    if [ ! -f $myDir/data/bowtie_index/Mus_musculus.GRCm38.dna.primary_assembly.1.ebwt ]; then
	echo "There is no bowtie genome index. Going ahead with generating it." 
	mkdir -p $myDir/data/bowtie_index # Make the folder for the mouse genome bowtie index
	if [ ! -f $refGenomeF ]; then
	    echo "Cannot find the mouse genome fasta file." 
	    echo "You need to dowload it and save it as $refGenomeF"
	    exit 1;
	fi
	echo "Reference mouse genome file $refGenomeF" 
	$bowtieBuildEx --threads 4 $refGenomeF $myDir/data/bowtie_index/Mus_musculus."$refGenome".dna.primary_assembly
	echo "Finished building index"
    fi
}


function mapReads { # Trim, filter and map reads to the genome 
    
    mkdir -p $myDir/out/cutadapt/ # Make folder for trimmed and length filtered reads
    mkdir -p $myDir/out/bowtie/ # Make folder for mapping reads
    
    cd $myDir/data/ || exit 1
    
    echo `pwd`
    
    for s in `ls *.fastq*|sed "s/.fastq.*//"`; do
	echo
	echo "$(date) : Processing $s"
	
	## Trim reads (if necessary)
	if [ ! -f $myDir/out/cutadapt/$s.trimmed.fastq.gz ]; then	    
	    if [ ! -f $myDir/out/cutadapt/$s.trimmed.fastq ]; then
		
		## If the input is gzippped, gunzip the fastq file
		if [ -f $s.fastq.gz ]; then
		    cat $s.fastq.gz | gunzip > $s.fastq 
		fi
		
		echo "$(date) : Running cutadapt"
		$cutadaptEx -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -O 9 -j 0 -m 19 -M 36 --trimmed-only $myDir/data/$s.fastq > $myDir/out/cutadapt/$s.trimmed.fastq
	    fi
	fi
	
	## Filter reads (if necessary)
	if [ ! -f $myDir/out/cutadapt/$s.trimmed.filtered.fastq.gz ]; then
	    if [ ! -f $myDir/out/cutadapt/$s.trimmed.filtered.fastq ]; then
		echo "$(date) : Running fastq quality filter"
		$fastqQFEx -q 30 -p 90 -Q 33 -i $myDir/out/cutadapt/$s.trimmed.fastq -o $myDir/out/cutadapt/$s.trimmed.filtered.fastq
	    fi
	fi
	
	## Map reads to the genome (if necessary)
	if [ ! -f $myDir/out/bowtie/$s.bam ]; then
	    echo "$(date) : Mapping with bowtie"
	    $bowtieEx --threads 4 -t -M 1 --best --strata -v 1 -q -p 7 --seed 666 -S $myDir/data/bowtie_index/Mus_musculus.GRCm38.dna.primary_assembly $myDir/out/cutadapt/$s.trimmed.filtered.fastq | $samtoolsEx view -F 4 -bS - -o $myDir/out/bowtie/$s.unsorted.bam		
	    echo "$(date) : Sorting bam"
	    $samtoolsEx sort -@ 4 -m 4G $myDir/out/bowtie/$s.unsorted.bam > $myDir/out/bowtie/$s.bam
	    echo "$(date) : Indexing bam"
	    $samtoolsEx index $myDir/out/bowtie/$s.bam
	    echo "$(date) : Removing trimmed, filtered fastqs and unsorted bam files"
	    rm $myDir/out/cutadapt/$s.trimmed.fastq  
	    rm $myDir/out/cutadapt/$s.trimmed.filtered.fastq  
	    rm $myDir/out/bowtie/$s.unsorted.bam  
	fi
    done
}

echo "Working directory $myDir"
#buildBowtieIndex
mapReads
echo "End of pipeline"
