# Cleanup
myvars <- ls()
rm(list=myvars)
rm(myvars)

######################
# Script for Fig5A ###
# Nocturning is ENSMUSG00000023087
# Genome 1 is 129(IAP-) and Genome 2 is DBA (IAP+)
######################

# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"

# Gene ID, gene length and counts for all samples 
G1_f <- "out/featureCounts/GSE35005/GSE35005.gene.genome1_norepeats.featureCounts.s0.tsv.gz" 
G2_f <- "out/featureCounts/GSE35005/GSE35005.gene.genome2_norepeats.featureCounts.s0.tsv.gz" 
bi_f <- "out/featureCounts/GSE35005/GSE35005.gene.unassigned_norepeats.featureCounts.s0.tsv.gz" 

# Reads mapping uniquely to each of the two genomes of the hybrid
G1 <-  read.table(gzfile(paste0(dir, G1_f)),  header = TRUE, row.names = 1)
G2 <-  read.table(gzfile(paste0(dir, G2_f)),  header = TRUE, row.names = 1)

glengths <- G1[,1] # Gene lengths, needed for RPKM calculation

G1 <- G1[,-1]
G2 <- G2[,-1]

# Total number of reads = reads mapping uniquely to each genome + those mapping to both
total <-  read.table(gzfile(paste0(dir, bi_f)),  header = TRUE, row.names = 1)[,-1] + G1 + G2

totalPerSample <- colSums(total)
rpkm <- (10^9)*t(t(total/glengths)/totalPerSample)

priSGA <- rowMeans(rpkm[,c("SRR398037","SRR398038")])
SGA <- rowMeans(rpkm[,c("SRR398039","SRR398040")])
SGB <- rowMeans(rpkm[,c("SRR398041","SRR398042")])
lepSC <- rowMeans(rpkm[,c("SRR398043","SRR398044")])
pacSC <- rowMeans(rpkm[,c("SRR398045","SRR398046")])
rST <- rowMeans(rpkm[,c("SRR398047","SRR398048")])
eST <- rowMeans(rpkm[,c("SRR398049","SRR398050")])
SE <- rowMeans(rpkm[,c("SRR398051","SRR398052")])

allGeneMeans <- cbind(priSGA, SGA, SGB, lepSC, pacSC, rST, eST, SE)
allGeneMeans<- allGeneMeans[rowSums(allGeneMeans>=1)>1,]

# Fig5 A, upper panel
barplot(matrix(data=rep(1,16),nrow=2,ncol=8, byrow=F), 
        beside=T, col="white", ylim=c(0,1), ylab="Relative allele expression"
)
barplot(
  matrix(data=rep(.5,16),nrow=2,ncol=8, byrow=F), 
  beside=T, col="white", ylim=c(0,1), add=T
)
barplot(
  matrix(
    data=as.numeric(G2["ENSMUSG00000023087",]/(G1["ENSMUSG00000023087",]+G2["ENSMUSG00000023087",])),
    nrow=2,ncol=8, byrow=F
  ), 
  beside=T, col="red", ylim=c(0,1), add=T
)

# Fig5 A, lower panel
boxplot(log(allGeneMeans+1), outline=F, ylab="Gene expression (log2 RPKM)", col="white")
points(allGeneMeans["ENSMUSG00000023087",], pch=19, col="red")
