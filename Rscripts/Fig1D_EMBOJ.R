# Cleanup
myvars <- ls()
rm(list=myvars)
rm(myvars)

######################
# Script for Fig1D ###
######################

# Load packages
library(DESeq2)
library(lme4breeding)
library(sva)
library(pedigreemm)
#library(pedigree)
library(dplyr)
library(vsn)
library(ggplot2) # for the 2d density plot

# Filtering of clusters
lowCountThreshold <- 10

# Parent dir
dir <- "~/Projects/mouse-piRNA-variation/"

# Table with all sample identifiers and sample info
t_f <- "data/EV1.csv" 
 
# Cluster coords
clustercoords_f <- "data/piRNA.gtf"
 
# Cluster ID, coords and raw counts
count_f <- "data/ICR_piRNA_norepeats_featureCounts_s0.tsv"

# Pedigree
pedigree_f <- "data/PedigreeICR.csv"
 
# Read piRNA counts
countsDF <-  read.table(
  paste0(dir,count_f),  
  header = TRUE, sep="\t", 
  row.names = 1, 
  fill=TRUE)
countsDF <- na.omit(countsDF)

# Create GRanges object with cluster coordinates
clustersDF <- read.table(paste0(dir,clustercoords_f))
clustersGR <- GRanges(
seqnames=clustersDF[,1],
   IRanges(start=clustersDF[,4], end=clustersDF[,5]),
   strand=clustersDF[,7],
   clusternames=clustersDF[,10]
 )
 
# Get strain names for testis samples
t_names <- read.csv(
    paste(dir, t_f, sep = ""), 
    header=T, 
    row.names="SampleID"
    )[,c("MouseID","Batch","AncestorOverfed"), drop=F]
t_names$Batch <- factor(t_names$Batch)
t_names$Generation <- nchar(gsub("\\d","",t_names$MouseID))
t_names$Generation <- factor(t_names$Generation)

# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = countsDF,
  colData =  t_names[match(names(countsDF), row.names(t_names)),, drop=F],
  design = ~ 1)
dds <- dds[rowSums(counts(dds)) >= lowCountThreshold,] # Remove clusters not expr

# Run DESeq on the counts
dds <- DESeq(dds)
dds$Batch <- droplevels(dds$Batch)
 
# Get variance stabilised expression values
vsd <- vst(dds, nsub=nrow(dds), blind=T)

# Correct for batches
mymodelmatrix <- model.matrix(~1, dat=colData(dds))
correctedMatrix <- ComBat(
  dat=assay(vsd), 
  batch=colData(dds)$Batch, 
  mod=mymodelmatrix)
assay(vsd) <- correctedMatrix

# Remove "sample" from sample identifier to match the numbers in the pedigree
colnames(vsd) <- sub("sample0","",colnames(vsd))
colnames(vsd) <- sub("sample","",colnames(vsd))
 
### Encoded Pedigree info 
pedigree_data <- read.table(paste0(dir,pedigree_f), header=T, sep=",")
 
mysamples <- 1:39
mypiCs <- rownames(vsd)

# The pedigrees correspond to two sets: 
# One set of pedigrees are from an Ancestor (male F0) that was overnurished
# The second set of pedigrees are from an Ancestor that was standard nurished.
# To account for any possible effects caused by this grouping, we included 
# this as an extra factor.
myformula <- formula(piC ~ AncestralDiet + (1|mouseNum))
 
observed_heritabilities <- c()

probabilities_of_observed_heritabilities <-c()

for(mypiC in mypiCs){
   print(mypiC)
   myvsddata <- assay(vsd)[mypiC,mysamples]
   myvsddf <- data.frame(
     piC = myvsddata, 
     mouseNum = names(myvsddata), 
     Batch = colData(vsd)$Batch, 
     AncestralDiet = colData(vsd)$AncestorOverfed, 
     Generation = colData(vsd)$Generation)
 
   mypedigree_data <- pedigree_data[pedigree_data$ID %in% mysamples,]
 
 
   mypedigree <- pedigree(
     label=as.character(pedigree_data$mouseNum),
     dam = as.integer(pedigree_data$damNum),
     sire = as.integer(pedigree_data$sireNum)
   )
    A <- getA(mypedigree)
    
    Asubset <- A[row.names(A) %in% mysamples,row.names(A) %in% mysamples]
    
    mix1 <- lmebreed(myformula,
                     relmat = list(mouseNum=Asubset),
                     control = lmerControl(
                       check.nobs.vs.nlev = "ignore",
                       check.nobs.vs.rankZ = "ignore",
                       check.nobs.vs.nRE="ignore"
                     ),
                     data=myvsddf)
   
   vc <- VarCorr(mix1)
   genetic_variance <- attr(vc$mouseNum,"stddev")^2
   residual_variance <- attr(vc,"sc")^2
    
   anova_result <- anova(mix1)
   ancestralDietSS <- anova_result$`Sum Sq`[1]
   total_variance <- genetic_variance + residual_variance + ancestralDietSS
    
   piC_heritability  <- 100 *  genetic_variance / total_variance
    
   observed_heritabilities <- c(observed_heritabilities,piC_heritability)
    
   ### Calculate heritability in a permutated dataset
    piC_heritability_byChance_all <- c()
    permuted_myvsddf <- myvsddf
    for(i in 1:100){

      permuted_myvsddf <- myvsddf %>%
        group_by(Generation, AncestralDiet) %>% # Permute values within the same generation and AncestorDiet set
        mutate(mouseNum = sample(mouseNum))
      mix1_permuted <- lmebreed(myformula,
                                relmat = list(mouseNum=Asubset),
                                control = lmerControl(
                                  check.nobs.vs.nlev = "ignore",
                                  check.nobs.vs.rankZ = "ignore",
                                  check.nobs.vs.nRE="ignore"
                                ),
                                data=permuted_myvsddf)
    
      vc <- VarCorr(mix1_permuted)
      genetic_variance <- attr(vc$mouseNum,"stddev")^2
      residual_variance <- attr(vc,"sc")^2
    
      anova_result <- anova(mix1_permuted)
      ancestralDietSS <- anova_result$`Sum Sq`[1]
      total_variance <- genetic_variance + residual_variance + ancestralDietSS
    
      piC_heritability_byChance  <- 100 * genetic_variance / total_variance
      piC_heritability_byChance_all <- c(piC_heritability_byChance_all, piC_heritability_byChance)
    }
    probabilities_of_observed_heritabilities <- c(probabilities_of_observed_heritabilities,sum(piC_heritability_byChance_all>piC_heritability)/100)
}

# Adjust p-values using BH 
probabilities_of_observed_heritabilities_adjusted <- p.adjust(probabilities_of_observed_heritabilities, method="BH")

# Plot Fig 1D
plot(y=observed_heritabilities,x=probabilities_of_observed_heritabilities_adjusted, xlab="Probability of heritability by chance", ylab="Estimated heritability", col = rgb(0, 0, 0, alpha = 0.1), pch=19, )
abline(v=0.05, col="grey", lty=2)
