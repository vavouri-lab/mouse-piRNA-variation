## Cleanup 
myvars <- ls()
rm(list=myvars)
rm(myvars)

######################
# Script for Fig4B ###
######################

# Load packages
library(DESeq2)
library(beeswarm)
library(coin)

# Define parameters
lowCountThreshold <- 10
padj_sig <- 0.05 # Adj p value threshold
f_sig <- 1 # log2 fold change threshold
cflank <- 5000 # Max allowed distance to transposon.
mycex <- .5

# In Fig4 we are not interested in the strand of the cluster and the TE
ignoreStrand <- T

# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"

# Table with all sample identifiers and sample info
t_f <- "data/EV1.csv" 

# File with predicted cluster coords 
clustercoords_f <- "data/protrac_merged.mm10.gtf"

# File with predicted cluster counts
count_f <- "data/protrac_merged.rawcounts.tsv"

# File containing all TEVs from Nellaker et al, with coords converted to mm10
tev_f <- "data/NellakerGenomeBiol2012/allTEVs.shared.TEVID.mm9.mm10.tab.gz"  
ttypes <- c("LINE","SINE","ERV","IAP")

# Read piRNA counts
countsDF <-  read.table(paste(dir, count_f, sep = ""),  header = TRUE, row.names = 1)

# Create a GRanges object for piRNA clusters
clustersDF <- read.table(paste(dir,clustercoords_f, sep=""))
clustersDF[clustersDF[,7] == ".",7] <- "*" # Change strand "." to  "*"
clustersGR <- GRanges(
  seqnames = clustersDF[,1],
  IRanges(start=clustersDF[,4], end=clustersDF[,5]),
  strand = clustersDF[,7],
  clusternames = clustersDF[,10]
)


# Remove bidirectional clusters which overlap other clusters to avoid redundancy
clustersGR <- 
  clustersGR[!(countOverlaps(clustersGR) > 1 & strand(clustersGR) == "*")]
countsDF <-
  countsDF[rownames(countsDF) %in% clustersGR$clusternames, ]

# Create a GRanges object for TEVs 
tevDF <- read.table(gzfile(paste0(dir,tev_f)), header=T)
tevDF <- tevDF[!is.na(tevDF[,"start_mm10"]),] # Remove any TEVs that did not convert to mm10
tevDF[tevDF[,"strand_mm10"] == "N","strand_mm10"] <- "*" # Some have strand N
tevDF[tevDF[,"strand_mm10"] == "0","strand_mm10"] <- "*" # Some have strand 0
tevGR <- GRanges(
  seqnames=gsub('chr','',tevDF[,"chr_mm10"]),
  IRanges(start=tevDF[,"start_mm10"],end=tevDF[,"end_mm10"]),
  strand=tevDF[,"strand_mm10"],
  TEVtype=tevDF[,"TEtype"],
  XBL6=tevDF[,"C57B6_ref"],
  XCAST=tevDF[,"CAST"]
)

# Get strain names for samples
t_names <- read.csv(paste(dir, t_f, sep = ""), header=T, row.names="SampleID")[,c("Strain","Tissue"), drop=F]
t_names[,"Strain"]<-paste0("X",t_names[,"Strain"]) #  Adding X to all strain names, due to 129
t_names$Strain <- factor(t_names$Strain)

# Make BL6 the reference
t_names$Strain <- relevel(t_names$Strain, "XBL6")

# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countsDF,
                              colData =  t_names[match(names(countsDF),row.names(t_names)),, drop=F],
                              design = ~ Strain)

dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] # Remove clusters not expr

dds <- dds[,dds$Tissue == "spermatogonia"]
dds$Strain <- droplevels(dds$Strain)

# Run DESeq
dds <- DESeq(dds)

strain1 <- "XBL6"
strain2 <- "XCAST"

myres <- na.omit(results(dds, contrast = c("Strain", strain1, strain2), alpha=padj_sig, lfcThreshold=f_sig,  altHypothesis = "greaterAbs"))

# Get clusters overlapping any TEVs
NE_withTEV <- subsetByOverlaps(clustersGR, tevGR, maxgap = cflank, ignore.strand=T)$clusternames

# Get clusters with TE only in strain1
NE_SINE_strain1plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    tevGR$TEVtype == "SINE"
    & ((elementMetadata(tevGR[,strain1])[,1]==1 & elementMetadata(tevGR[,strain2])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain1])[,1]=="INS" & elementMetadata(tevGR[,strain2])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_LINE_strain1plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype == "LINE" | tevGR$TEVtype == "LINE_frag")
    & ((elementMetadata(tevGR[,strain1])[,1]==1 & elementMetadata(tevGR[,strain2])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain1])[,1]=="INS" & elementMetadata(tevGR[,strain2])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_ERV_strain1plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype != "SINE" & tevGR$TEVtype != "LINE" & tevGR$TEVtype != "LINE_frag"  )
    & ((elementMetadata(tevGR[,strain1])[,1]==1 & elementMetadata(tevGR[,strain2])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain1])[,1]=="INS" & elementMetadata(tevGR[,strain2])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_IAP_strain1plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype == "IAP-I")
    & ((elementMetadata(tevGR[,strain1])[,1]==1 & elementMetadata(tevGR[,strain2])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain1])[,1]=="INS" & elementMetadata(tevGR[,strain2])[,1]==0)),
  ], maxgap = cflank)$clusternames


# Get clusters with TE only in strain2
NE_SINE_strain2plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    tevGR$TEVtype == "SINE"
    & ((elementMetadata(tevGR[,strain2])[,1]==1 & elementMetadata(tevGR[,strain1])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain2])[,1]=="INS" & elementMetadata(tevGR[,strain1])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_LINE_strain2plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype == "LINE" | tevGR$TEVtype == "LINE_frag")
    & ((elementMetadata(tevGR[,strain2])[,1]==1 & elementMetadata(tevGR[,strain1])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain2])[,1]=="INS" & elementMetadata(tevGR[,strain1])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_ERV_strain2plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype != "SINE" & tevGR$TEVtype != "LINE" & tevGR$TEVtype != "LINE_frag"  )
    & ((elementMetadata(tevGR[,strain2])[,1]==1 & elementMetadata(tevGR[,strain1])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain2])[,1]=="INS" & elementMetadata(tevGR[,strain1])[,1]==0)),
  ], maxgap = cflank)$clusternames

NE_IAP_strain2plus <- subsetByOverlaps(
  ignore.strand = ignoreStrand, clustersGR, tevGR[
    (tevGR$TEVtype == "IAP-I")
    & ((elementMetadata(tevGR[,strain2])[,1]==1 & elementMetadata(tevGR[,strain1])[,1]=="DEL") |
         (elementMetadata(tevGR[,strain2])[,1]=="INS" & elementMetadata(tevGR[,strain1])[,1]==0)),
  ], maxgap = cflank)$clusternames

# Put all the data in two lists for plotting
myListStrain1=list(noTEV=myres[!(row.names(myres) %in% NE_withTEV),"log2FoldChange"],
                   LINE=myres[row.names(myres) %in% NE_LINE_strain1plus,"log2FoldChange"],
                   SINE=myres[row.names(myres) %in% NE_SINE_strain1plus,"log2FoldChange"],
                   ERV=myres[row.names(myres) %in% NE_ERV_strain1plus,"log2FoldChange"],
                   IAP=myres[row.names(myres) %in% NE_IAP_strain1plus,"log2FoldChange"]
)
myListStrain2=list(LINE=myres[row.names(myres) %in% NE_LINE_strain2plus,"log2FoldChange"],
                   SINE=myres[row.names(myres) %in% NE_SINE_strain2plus,"log2FoldChange"],
                   ERV=myres[row.names(myres) %in% NE_ERV_strain2plus,"log2FoldChange"],
                   IAP=myres[row.names(myres) %in% NE_IAP_strain2plus,"log2FoldChange"]
)

beeswarm(myListStrain1, main=paste(strain1," vs ", strain2),
         pch=c(19,24,24,24,24), cex=mycex, corral="random",
         method="compactswarm",
         col=c("black","red","red","red","red"),bg="grey",
         ylab="Log2FoldChange", ylim=c(-12,12))

beeswarm(myListStrain2,pch=25, cex=mycex, corral="random",
         method="compactswarm", col="blue", bg="grey",at=c(2:5), add=TRUE)

# Add lines where the means are 
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
# Test whether the two groups come from the same distribution using the Wilcoxon rank sum test
if(length(NE_LINE_strain1plus)>0 & length(NE_LINE_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_LINE_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_LINE_strain2plus,"log2FoldChange"]))
if(length(NE_SINE_strain1plus)>0 & length(NE_SINE_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_SINE_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_SINE_strain2plus,"log2FoldChange"]))
if(length(NE_ERV_strain1plus)>0 & length(NE_ERV_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_ERV_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_ERV_strain2plus,"log2FoldChange"]))
if(length(NE_IAP_strain1plus)>0 & length(NE_IAP_strain2plus)>0) print(wilcox.test(myres[row.names(myres) %in% NE_IAP_strain1plus,"log2FoldChange"],myres[row.names(myres) %in% NE_IAP_strain2plus,"log2FoldChange"]))



