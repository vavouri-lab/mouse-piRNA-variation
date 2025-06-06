## Cleanup 
myvars <- ls()
rm(list=myvars)
rm(myvars)

######################
# Script for Fig4A ###
######################

# Load packages
library(DESeq2)

# Define parameters
lowCountThreshold <- 10
padj_sig <- 0.05 # Adj p value threshold
f_sig <- 1 # log2 fold change threshold
cflank <- 5000 # Max allowed distance to transposon.

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
  XNOD=tevDF[,"NOD"],
  X129=tevDF[,"X129S1"],
  XC3H=tevDF[,"C3H"],
  XCAST=tevDF[,"CAST"]
)

# Get strain names for samples
t_names <- read.csv(paste(dir, t_f, sep = ""), header=T, row.names="SampleID")[,c("Strain","Tissue"), drop=F]
t_names[,"Strain"]<-paste0("X",t_names[,"Strain"]) #  Adding X to all strain names, due to 129
t_names$Strain <- factor(t_names$Strain)

# Make BL6 the reference
t_names$Strain <- relevel(t_names$Strain, "XBL6")

# In Fig4A there is testis as well as spermatogonia data. Run first for testis
# then for spermatogonia (just BL6 vs CAST)
ors <- c() # Odds ratios
pvals <- c() # pval from fishers test
strainPairs <- c()

for (tissue in c("testis","spermatogonia")){
  if (tissue == "testis"){
    myStrains <- c("XBL6", "XNOD", "XC3H", "X129")
  }else{
    myStrains <- c("XBL6", "XCAST")
  }

  # Generate DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = countsDF,
                                colData =  t_names[match(names(countsDF),row.names(t_names)),, drop=F],
                                design = ~ Strain)
  
  dds <- dds[,dds$Strain!="XICR"] # Keep data only for inbred strains
  dds <- dds[,dds$Tissue== tissue]
  dds$Strain <- droplevels(dds$Strain)
  dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] # Remove clusters not expr
  
  # Run DESeq
  dds <- DESeq(dds)
  

  # Get the numbers of diff expressed clusters overlapping or not different types of TEVs
  # for barplots in Fig4A
  for(i in 1:(length(myStrains)-1)){
    for(j in (i+1):(length(myStrains))){
      strain1 <- myStrains[i]
      strain2 <- myStrains[j]
      if(strain1 == "129" ){
        strain1 <- "X129"
      }else if(strain2 == "129"){
        strain2 <- "X129"
      }
      myres <- na.omit(results(dds, contrast = c("Strain", strain1, strain2), alpha=padj_sig, lfcThreshold=f_sig,  altHypothesis = "greaterAbs"))
      strainPairs <- c(strainPairs,paste(strain1,"vs",strain2))
      
      # Get clusters overlapping TEVs
      NE_LINE_dif <- subsetByOverlaps(
        ignore.strand = ignoreStrand, clustersGR, tevGR[
          (tevGR$TEVtype == "LINE" | tevGR$TEVtype == "LINE_frag")
          & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),
        ], maxgap = cflank)$clusternames
      
      NE_SINE_dif <- subsetByOverlaps(
        ignore.strand = ignoreStrand, clustersGR, tevGR[
          (tevGR$TEVtype == "SINE")
          & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),
        ], maxgap = cflank)$clusternames
      
      NE_ERV_dif <- subsetByOverlaps(
        ignore.strand = ignoreStrand, clustersGR, tevGR[
          (tevGR$TEVtype != "SINE" & tevGR$TEVtype != "LINE" & tevGR$TEVtype != "LINE_frag"  )
          & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),
        ], maxgap = cflank)$clusternames
      
      NE_IAP_dif <- subsetByOverlaps(
        ignore.strand = ignoreStrand, clustersGR, tevGR[
          (tevGR$TEVtype == "IAP-I")
          & (elementMetadata(tevGR[,strain1])[,1] != elementMetadata(tevGR[,strain2])[,1]),
        ], maxgap = cflank)$clusternames
      
      # #       }
      NE_TEV_dif_list <- list(LINE=NE_LINE_dif,SINE=NE_SINE_dif,ERV=NE_ERV_dif,IAP=NE_IAP_dif)
      
      sig <- row.names(myres[which(myres$padj < padj_sig & abs(myres$log2FoldChange)>f_sig),])
      notSig <- row.names(myres[(!row.names(myres) %in% row.names(myres[which(myres$padj < padj_sig & abs(myres$log2FoldChange)>f_sig),])),])
      
      for(k in 1:length(ttypes)){
        ttype <- ttypes[k]
        NE_TEV_dif <- unlist(NE_TEV_dif_list[ttype])
        sigWithTEV <- sig[sig %in% NE_TEV_dif]
        sigWithoutTEV <- sig[!(sig %in% NE_TEV_dif)]
        notSigWithTEV <- notSig[notSig %in% NE_TEV_dif]
        notSigWithoutTEV <- notSig[!(notSig %in% NE_TEV_dif)]
        m <- matrix(c(length(sigWithTEV), length(notSigWithTEV), length(sigWithoutTEV),length(notSigWithoutTEV)),nrow=2, dimnames=list(c("sig","notSig"),c("withTEV","withoutTEV")))
        f <-  fisher.test(m)
        ors<-c(ors,round(f$estimate[[1]],2))
        pvals<-c(pvals,round(f$p.value,4))
        barplot(proportions(m, margin=2), main=ttype, sub=paste(strain1,"vs",strain2))
      }
    }
  }
}  
ors <- data.frame(matrix(ors,nrow=length(strainPairs),ncol=4,byrow=T), row.names=strainPairs)
names(ors) <- c("LINEs","SINEs","ERVs","IAPs")
pvals <- data.frame(matrix(pvals,nrow=length(strainPairs),ncol=4,byrow=T), row.names=strainPairs)
names(pvals) <- c("LINEs","SINEs","ERVs","IAPs")

print("Log odds ratios:")
print(log(ors))

print("P-values:")
print(pvals)


