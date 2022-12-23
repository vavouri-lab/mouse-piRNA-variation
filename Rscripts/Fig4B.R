## Cleanup before starting
myvars <- ls()
rm(list = myvars)
rm(myvars)
 
####################################################################
# Script to generate Figure 4B
####################################################################
 
##################################
# load required libraries
library(DESeq2)
# library(dplyr)
# library("RColorBrewer") # Needed for the figures
# library("vioplot")
library(beeswarm)
##################################
 
##################################
# Define parameters
##################################
# Filtering of clusters
removeLowCountClusters <- T
lowCountThreshold <- 10
removeOverlapingClusters <- T
# Adj p value threshold
padj_sig <- 0.05 
# log2 fold change threshold
f_sig <- 1 
# Include this flanking region to get overlap with TEs
cflank <- 5000 
###################################
 
##################################
# Drawing settings
##################################
mycex <- .5

##################################
# set up variables
##################################
# Project home dir
dir <- "~/Projects/mouse-piRNA-variation/" 

# Big table with all sample identifiers and other info
t_f <- "data/sampleInfo.txt" 
# File with cluster ID and counts for all samples
count_f <- "out/featureCounts/piRNAPredicted/norepeats/piRNAPredicted.norepeats.featureCounts.s0.tsv" 
# File with coordinates of final set of predicted piRNA clusters
clustercoords_f <- "out/protrac/piRNAPredicted.gtf"

myStrains <- c("BL6", "CAST")
 
# File containing all TEVs from Nellaker et al
tev_f <- "data/NellakerGenomeBiol2012/allTEVs.shared.TEVID.mm9.mm10.tab"
ttypes <- c("LINE", "SINE", "ERV", "IAP")
##################################
 
######################################################################
# Get piRNA counts
######################################################################
countsDF <-  read.table( paste0(dir, count_f),  header = TRUE, row.names = 1)
 
################################################################
# Get cluster coordinates 
################################################################
# Create a GRanges object with cluster coordinates 
clustersDF <- read.table(paste0(dir,clustercoords_f))
clustersDF[clustersDF[, 7] == ".", 7] <- "*" # Bidirectional clusters
clustersGR <- GRanges(
  seqnames = clustersDF[, 1],
  IRanges(start = clustersDF[, 4], end = clustersDF[, 5]),
  strand = clustersDF[, 7],
  clusternames = clustersDF[, 10]
)
 
# Reduce redundancy by removing bidirectional clusters which overlap others
if (removeOverlapingClusters) {
  clustersGR <- clustersGR[!(countOverlaps(clustersGR) > 1 
                             & strand(clustersGR) == "*")]
  countsDF <- countsDF[rownames(countsDF) %in% clustersGR$clusternames, ]
}
 
################################################################
# Get positions of polymorphic transposons
################################################################
# Read TEV coordinates into a dataframe.
tevDF <- read.table(paste0(dir, tev_f), header = T)
# Remove any TEVs that did not convert to mm10
tevDF <- tevDF[!is.na(tevDF[, "start_mm10"]),] 
# Strand is "None" according to Nellaker et al
tevDF[tevDF[, "strand_mm10"] == "N", "strand_mm10"] <- "*"
# Strand is 0 according to Nellaker et al
tevDF[tevDF[, "strand_mm10"] == "0", "strand_mm10"] <- "*" 

# Create a GRanges object with TEV coordinates
tevGR <- GRanges(
  seqnames=gsub('chr', '', tevDF[, "chr_mm10"]),
  IRanges(start = tevDF[, "start_mm10"], end = tevDF[, "end_mm10"]),
  strand = tevDF[, "strand_mm10"],
  TEVtype = tevDF[, "TEtype"],
  BL6 = tevDF[, "C57B6_ref"],
  CAST = tevDF[, "CAST"]
  )

################################################################
#  Get strain names for samples
################################################################
# Get strain names for testis samples
t_names <- read.table(
  paste(dir, t_f, sep = ""), 
  header = T, row.names = "SampleID")[, c("Strain", "Tissue"), drop = F]
# Convert to factors
t_names$Strain <- factor(t_names$Strain) 
# Make sure first level is the control level
t_names$Strain <- relevel(t_names$Strain, "BL6")


################################################################
# Test for differential expression
################################################################
# Generate DESeqDataSet from count data
dds <- DESeqDataSetFromMatrix(countData = countsDF,
                              colData =  t_names[match(names(countsDF), row.names(t_names)), , drop = F],
                              design = ~ Strain)

if(removeLowCountClusters){
  dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,]
}

# Keep only data for BL6 and CAST spermatogonia
dds <- dds[, dds$Strain != "ICR"]
dds <- dds[, dds$Tissue != "testis"]
dds$Strain <- droplevels(dds$Strain)

# Run DESeq
dds <- DESeq(dds)
 
# par(mfrow=c(2,2))
 
strain1 <- myStrains[1]
strain2 <- myStrains[2]

myres <- na.omit(
  results(
    dds,
    contrast = c("Strain", strain1, strain2),
    alpha = padj_sig,
    lfcThreshold = f_sig,
    altHypothesis = "greaterAbs"
  )
)

# Get clusters overlapping any TEVs
NE_withTEV <- subsetByOverlaps(clustersGR,
                               tevGR,
                               maxgap = cflank,
                               ignore.strand = T)$clusternames

# In the supplement of Nellaker et al, TEVs are encoded as follows:
# 1   : BL6 has the TE and also the other strain has it 
# 0   : BL6 does not have the TE and also the other strain does not have it
# INS : BL6 does not have the TE but the other strain has it
# DEL : BL6 has the TE but the other strain does not have it

# Get clusters with SINEs only in BL6
NE_SINE_strain1plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[tevGR$TEVtype == "SINE"
        & elementMetadata(tevGR[, strain2])[, 1] == "DEL", ],
  maxgap = cflank)$clusternames

# Get clusters with LINEs only in BL6
NE_LINE_strain1plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype == "LINE"
         | tevGR$TEVtype == "LINE_frag")
        &  elementMetadata(tevGR[, strain2])[, 1] == "DEL", ],
  maxgap = cflank)$clusternames

# Get clusters with ERVs only in BL6
NE_ERV_strain1plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype != "SINE"
         & tevGR$TEVtype != "LINE"
         & tevGR$TEVtype != "LINE_frag")
        & elementMetadata(tevGR[, strain2])[, 1] == "DEL", ],
  maxgap = cflank)$clusternames

# Get clusters with IAPs only in BL6
NE_IAP_strain1plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype == "IAP-I")
        & elementMetadata(tevGR[, strain2])[, 1] == "DEL", ],
  maxgap = cflank)$clusternames


###############################################
# Get clusters with SINEs only in CAST
NE_SINE_strain2plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[tevGR$TEVtype == "SINE"
        & elementMetadata(tevGR[, strain2])[, 1] == "INS", ],
  maxgap = cflank)$clusternames

# Get clusters with LINEs only in CAST
NE_LINE_strain2plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype == "LINE"
         | tevGR$TEVtype == "LINE_frag")
        & elementMetadata(tevGR[, strain2])[, 1] == "INS", ],
  maxgap = cflank)$clusternames

# Get clusters with ERVs only in CAST
NE_ERV_strain2plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype != "SINE"
         & tevGR$TEVtype != "LINE"
         & tevGR$TEVtype != "LINE_frag")
        & elementMetadata(tevGR[, strain2])[, 1] == "INS", ],
  maxgap = cflank)$clusternames

# Get clusters with IAPs only in CAST
NE_IAP_strain2plus <- subsetByOverlaps(
  ignore.strand = T,
  clustersGR,
  tevGR[(tevGR$TEVtype == "IAP-I")
        & elementMetadata(tevGR[, strain2])[, 1] == "INS", ],
  maxgap = cflank)$clusternames

myListStrain1 = list(
  noTEV = myres[!(row.names(myres) %in% NE_withTEV), "log2FoldChange"],
  LINE = myres[row.names(myres) %in% NE_LINE_strain1plus, "log2FoldChange"],
  SINE = myres[row.names(myres) %in% NE_SINE_strain1plus, "log2FoldChange"],
  ERV = myres[row.names(myres) %in% NE_ERV_strain1plus, "log2FoldChange"],
  IAP = myres[row.names(myres) %in% NE_IAP_strain1plus, "log2FoldChange"]
)

myListStrain2 = list(
  LINE = myres[row.names(myres) %in% NE_LINE_strain2plus, "log2FoldChange"],
  SINE = myres[row.names(myres) %in% NE_SINE_strain2plus, "log2FoldChange"],
  ERV = myres[row.names(myres) %in% NE_ERV_strain2plus, "log2FoldChange"],
  IAP = myres[row.names(myres) %in% NE_IAP_strain2plus, "log2FoldChange"]
)


       beeswarm(myListStrain1, main=paste(strain1," vs ", strain2),
                pch=c(19,24,24,24,24), cex=mycex, corral="random", method="compactswarm",
                col=c("black","red","red","red","red"),bg="grey",
                ylab="Log2FoldChange", ylim=c(-12,12))

       beeswarm(myListStrain2,pch=25, cex=mycex, corral="random", method="compactswarm",
                col="blue", bg="grey",at=c(2:5), add=TRUE)


       lines(x=c(0.8,1.2),y=c(summary(myres[!(row.names(myres) %in% NE_withTEV),"log2FoldChange"])[4],summary(myres[!(row.names(myres) %in% NE_withTEV),"log2FoldChange"])[4]), col="darkgrey", lwd=2)
       lines(x=c(1.8,2.2),y=c(summary(myres[row.names(myres) %in% NE_LINE_strain1plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_LINE_strain1plus,"log2FoldChange"])[4]), col="red", lwd=2)
       lines(x=c(2.8,3.2),y=c(summary(myres[row.names(myres) %in% NE_SINE_strain1plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_SINE_strain1plus,"log2FoldChange"])[4]), col="red", lwd=2)
       lines(x=c(3.8,4.2),y=c(summary(myres[row.names(myres) %in% NE_ERV_strain1plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_ERV_strain1plus,"log2FoldChange"])[4]), col="red", lwd=2)
       lines(x=c(4.8,5.2),y=c(summary(myres[row.names(myres) %in% NE_IAP_strain1plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_IAP_strain1plus,"log2FoldChange"])[4]), col="red", lwd=2)

       lines(x=c(1.8,2.2),y=c(summary(myres[row.names(myres) %in% NE_LINE_strain2plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_LINE_strain2plus,"log2FoldChange"])[4]), col="blue", lwd=2)
       lines(x=c(2.8,3.2),y=c(summary(myres[row.names(myres) %in% NE_SINE_strain2plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_SINE_strain2plus,"log2FoldChange"])[4]), col="blue", lwd=2)
       lines(x=c(3.8,4.2),y=c(summary(myres[row.names(myres) %in% NE_ERV_strain2plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_ERV_strain2plus,"log2FoldChange"])[4]), col="blue", lwd=2)
       lines(x=c(4.8,5.2),y=c(summary(myres[row.names(myres) %in% NE_IAP_strain2plus,"log2FoldChange"])[4],summary(myres[row.names(myres) %in% NE_IAP_strain2plus,"log2FoldChange"])[4]), col="blue", lwd=2)

       print(paste(strain1,"vs",strain2))
       if(length(NE_LINE_strain1plus)>0 & length(NE_LINE_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_LINE_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_LINE_strain2plus,"log2FoldChange"]))
       if(length(NE_SINE_strain1plus)>0 & length(NE_SINE_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_SINE_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_SINE_strain2plus,"log2FoldChange"]))
       if(length(NE_ERV_strain1plus)>0 & length(NE_ERV_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_ERV_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_ERV_strain2plus,"log2FoldChange"]))
       if(length(NE_IAP_strain1plus)>0 & length(NE_IAP_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_IAP_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_IAP_strain2plus,"log2FoldChange"]))

       grid()

#       ## The following is a wilcoxon mann whitney test that can cope with ties. I didn't use these p-values in the paper but they are actually extremely similar to the pvalues of the test that I used (above).
#       library(coin)
#       wilcox_test(fc ~ strainplus, mydf)

#dev.off()
 
