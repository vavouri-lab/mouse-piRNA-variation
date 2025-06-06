# Cleanup
myvars <- ls()
rm(list=myvars)
rm(myvars)

##############################
# Script for Fig2A, 2C, 2D ###
##############################

# Load packages
library(DESeq2)
library("pheatmap")

# Filtering of clusters
lowCountThreshold <- 10
padj_sig <- 0.05 # p adj threshold
f_sig <- 1 # fold change threshold for DESeq

# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"

# Table with all sample identifiers and sample info
t_f <- "data/EV1.csv" 

# File with predicted cluster coords 
clustercoords_f <- "data/protrac_merged.mm10.gtf"

# File with predicted cluster counts
count_f <- "data/protrac_merged.rawcounts.tsv"

# Strain colors (used in the heatmap)
myStrainColors = list(Strain = c(XBL6 = "grey30", X129 = "beige", XC3H = "orange", XNOD = "darkred", XCAST = "lightgoldenrod1")) 

# Read piRNA counts
countsDF <-  read.table(paste(dir, count_f, sep = ""),  header = TRUE, row.names = 1)

# Create a GRanges object
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

dds <- dds[,dds$Strain!="XICR"] # Keep data only for inbred strains
dds$Strain <- droplevels(dds$Strain)
dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] # Remove clusters not expr

# Run DESeq on the testis samples
ddsTestis <- dds[,dds$Tissue=="testis"]
ddsTestis$Strain <- droplevels(ddsTestis$Strain)
ddsTestis <- DESeq(ddsTestis)

# Get the significantly differentially expressed predicted clusters in testes
myStrains <- c("XBL6", "XNOD", "XC3H", "X129")
sigClustersTestis <- c()
numOfSigClusters <-c()
strainPairs <- c()
for(i in 1:(length(myStrains)-1)){
  for(j in (i+1):(length(myStrains))){
    strain1<- myStrains[i]
    strain2<- myStrains[j]
    myres<-na.omit(results(ddsTestis, contrast = c("Strain", strain1, strain2), alpha=padj_sig, lfcThreshold=f_sig, altHypothesis = "greaterAbs"))
    sigClustersTestis <- c(sigClustersTestis,row.names(myres[myres$"padj"<padj_sig,]))
    numOfSigClusters <- c(numOfSigClusters,length(row.names(myres[myres$"padj"<padj_sig,])))
    strainPairs <- c(strainPairs,paste(strain1, strain2))
      }
}
 

# Run DESeq on the spermatogonia samples
ddsSgonia <- dds[,dds$Tissue!="testis"]
ddsSgonia$Strain <- droplevels(ddsSgonia$Strain)
ddsSgonia <- DESeq(ddsSgonia)

# Get the significantly ddifferentially expressed predicted clusters in sgonia
sigClustersSgonia <- c()
strain1<- "XBL6"
strain2<- "XCAST"
myres<-na.omit(results(ddsSgonia, contrast = c("Strain", strain1, strain2), alpha=padj_sig, lfcThreshold=f_sig, altHypothesis = "greaterAbs"))
sigClustersSgonia <- c(sigClustersSgonia,row.names(myres[myres$"padj"<padj_sig,]))
numOfSigClusters <- c(numOfSigClusters,length(row.names(myres[myres$"padj"<padj_sig,])))
strainPairs <- c(strainPairs,paste(strain1,strain2))

dds <- estimateSizeFactors(dds)
ntdds <- normTransform(dds)
 
### Fig 2A heatmap
clusterOrder <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
pheatmap(assay(ntdds)[clusterOrder,], cluster_rows=F, cluster_cols=T, show_rownames=F, show_colnames=F, annotation_col=data.frame(row.names = colnames(dds), 'Strain' = dds$Strain), annotation_colors=myStrainColors)

# Fig 2C barplot
barplot(table(table(sigClustersTestis)), xlab="# of strain pairs showing sig diff in expr in testis", ylab="# of predicted clusters")
 
# Fig 2D barplot
barplot(numOfSigClusters, names=strainPairs, ylab="# of predicted clusters")
