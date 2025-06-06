# Cleanup
myvars <- ls()
rm(list=myvars)
rm(myvars)

##############################
# Script for Fig1A, 1B, 1E ###
##############################

# Load packages
library(DESeq2)
library(pheatmap)
library(tibble)
library(dplyr)
library(tidyr)
library(vioplot)

# Filtering of clusters
lowCountThreshold <- 10
padj_sig <- 0.05 # p adj threshold
f_sig <- 1 # fold change threshold for DESeq

# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"

# Table with all sample identifiers and sample info
t_f <- "data/EV1.csv" 

# Cluster ID, coords and raw counts
count_f <- "~/Projects/mouse-piRNA-variation/data/EV2.csv" 

# Drawing settings
mycex <- .5
mypch <- 20
mylwd <- .5
mySigCol <- "purple"
myNonSigCol <- "black"

# The X before the strain name is because of strain 129
myStrains <- c("XBL6", "XNOD", "XC3H", "X129")

# Strain colors (used in the heatmap)
myStrainColors = list(Strain = c(XBL6 = "grey30", X129 = "beige", XC3H = "orange", XNOD = "darkred")) 

# Read piRNA counts
countsDF <-  read.csv(count_f,  header = TRUE, row.names = 1, fill=TRUE)

# Create a GRanges object
clustersDF <- countsDF[,1:5]
clustersDF[clustersDF[,4] == ".",4] <- "*" 
clustersGR <- GRanges(
  seqnames=clustersDF[,1],
  IRanges(start=clustersDF[,2], end=clustersDF[,3]),
  strand=clustersDF[,4],
  clusternames=row.names(clustersDF)
)

# Remove bidirectional clusters which overlap other clusters to avoid redundancy
clustersGR <- 
  clustersGR[!(countOverlaps(clustersGR) > 1 & strand(clustersGR) == "*")]
countsDF <-
  countsDF[rownames(countsDF) %in% clustersGR$clusternames, ]

# Keep only columns with counts
countsDF <- countsDF[,-c(1:5)]

# Keep only clusters for which there are read counts in all samples
countsDF <- na.omit(countsDF)

# Get strain names for testis samples
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
dds <- dds[,dds$Tissue=="testis"] # Keep data only for testes
dds$Strain <- droplevels(dds$Strain)
dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] # Remove clusters not expr

# Run DESeq on the counts
dds <- DESeq(dds)


##############################
## - Fig1A scatterplots
##############################

# Get mean log2 transformed normalized counts (with 1 pseudocount) and min/max 
ntdds <- normTransform(dds) 
ntCountDF<-data.frame(t(assay(ntdds)),Strain=as.character(ntdds$Strain))
ntCountDFavg<-aggregate(subset(ntCountDF,select = -Strain),by=list(ntCountDF$Strain), FUN=mean)
row.names(ntCountDFavg)<-ntCountDFavg[,1]
ntCountDFavg <- t(ntCountDFavg[,-1])
myCountMin<- min(ntCountDFavg)
myCountMax <-  max(ntCountDFavg)

for(i in 1:(length(myStrains)-1)){
  for(j in (i+1):(length(myStrains))){
    strain1<- myStrains[i]
    strain2<- myStrains[j]
    myres<-
      na.omit(
        results(dds, contrast = c("Strain", strain1, strain2), alpha=padj_sig,
                lfcThreshold = f_sig,  altHypothesis = "greaterAbs")
      )
    sigClus <- 
      row.names(myres[myres$"padj"<padj_sig & abs(myres$log2FoldChange)>f_sig,])
    
    plot(ntCountDFavg[,strain2],ntCountDFavg[,strain1], pch=mypch, cex=mycex, 
         xlim=c(myCountMin,myCountMax), ylim=c(myCountMin,myCountMax), 
         col=myNonSigCol, 
         main = paste(strain1,"vs",strain2), 
         ylab=paste(strain1, "(log2 norm counts)"),
         xlab=paste(strain2, "(log2 norm counts)")
    ) 
    abline(a=1,b=1, col="grey", lty=2, lwd=mylwd)
    abline(a=-1,b=1, col="grey", lty=2, lwd=mylwd)
    abline(a=0,b=1, col="lightgrey", lty=1, lwd=mylwd)
    points(ntCountDFavg[,strain2],ntCountDFavg[,strain1], pch=mypch, col=myNonSigCol, cex=mycex) 
    points(ntCountDFavg[sigClus,strain2],ntCountDFavg[sigClus,strain1], pch=mypch, col=mySigCol, cex=mycex) # replotting so that points appear over lines
    text(ntCountDFavg[sigClus,strain2],ntCountDFavg[sigClus,strain1], labels=c(row.names(ntCountDFavg[sigClus,])))
  }
}


##############################
## - Fig1B Heatmap
##############################
clusterOrder <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
pheatmap(assay(ntdds)[clusterOrder,], cluster_rows=F, cluster_cols=T, show_rownames=F, show_colnames=F, annotation_col=data.frame(row.names = colnames(dds), 'Strain' = dds$Strain), annotation_colors=myStrainColors)


##############################
## - Fig1C Violinplot
##############################
colData(dds)$Sample <- row.names(colData(dds))
vsd <- vst(dds, nsub=nrow(dds), blind=F) # Get variance stabilised expr values
vsd_tibble <- rownames_to_column(as.data.frame(assay(vsd)), var="piC")
mytibble <- pivot_longer(vsd_tibble, cols=starts_with("sample"), names_to="Sample", values_to="expression" )
merged_tibble <- mytibble %>% left_join(as.data.frame(colData(dds)[,c(1,4)]), by=c("Sample" = "Sample"))
all_modelstats_strain <- merged_tibble %>%
  group_by(piC) %>% 
  do({
    model_withStrain <- lm(expression ~ Strain, data = .)
    model_summary <- summary(model_withStrain)
    mypf <- pf(
      q = model_summary$fstatistic[1], 
      df1 = model_summary$fstatistic[2], 
      df2 = model_summary$fstatistic[3], 
      lower.tail = FALSE
    )
    data.frame(piC = unique(.$piC), adj.r.squared = model_summary$adj.r.squared , p_value = mypf)
  })
all_modelstats_strain <- cbind(all_modelstats_strain, p.adjusted = p.adjust(all_modelstats_strain$p_value, method="BH"))
valuesForBarplot <- all_modelstats_strain$adj.r.squared
valuesForBarplot[valuesForBarplot<0] <- 0
vioplot(valuesForBarplot, ylim=c(0,1))


############################################################
## - Fig1E Dot plots of individual piRNA clusters
############################################################
sigClusters <- c("piC29","piC16","piC164","piC134")
for(i in sigClusters){
  plotCounts(dds, gene=i, intgroup="Strain", pch=mypch, transform=F, cex=2)
}
