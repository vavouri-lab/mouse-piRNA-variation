## Cleanup 
myvars <- ls()
rm(list=myvars)
rm(myvars)

######################
# Script for Fig3B ###
######################

# Load packages
library(DESeq2)
library(beeswarm)

# Filtering of clusters
lowCountThreshold <- 10
 
# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"
 
# Two column text file where the first column contains sample name 
# and the second column T when the sample has IAP in Noct, 
# F when the sample does not have IAP in Noct and NA when not known.
IAP_insertion_f <- "data/sample_NoctIAPstatus.txt"
  
# Cluster ID, coords and raw counts
count_f <- "~/Projects/mouse-piRNA-variation/data/ICR_piRNA_norepeats_featureCounts_s0.tsv"

# Read piRNA counts
countsDF <-  read.table(
  count_f,  
  header = TRUE, sep="\t", 
  row.names = 1, 
  fill=TRUE)
countsDF <- na.omit(countsDF)

# Get the info on the Noct IAP insertion (from PCR)
IAP_insertion <- read.table(paste(dir, IAP_insertion_f, sep = ""), header=F, row.names=1)

# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countsDF,
                              colData = data.frame("IAP_insertion" = IAP_insertion[colnames(countsDF),]), # make this a data.frame
                              design = ~ 1)
dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,]
dds <- estimateSizeFactors(dds)

# Plot piNoct expression, coloring by Noct IAP genotype
piC29NormCounts<-counts(dds, normalize=T)["piC29",]
piC29Color <- rep("grey",length(piC29NormCounts))
piC29data <- data.frame(piC29NormCounts,piC29Color)
piC29data[which(colData(dds)$IAP_insertion == 1),"piC29Color"]  <- "black"
piC29data[which(colData(dds)$IAP_insertion == 0),"piC29Color"] <- "white"
piC29data <- piC29data[which(!is.na(colData(dds)$IAP_insertion)),] # Remove sample with unknown genotype
beeswarm(piC29data[,1], pch=21, pwbg=piC29data[,2], method="compactswarm", 
         ylab = "pi-Noct (counts)")

