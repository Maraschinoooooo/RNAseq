## Installing WGCNA
library(WGCNA)
library(GO.db)
library(impute)
library(preprocessCore)
library(tidyverse)
library(magrittr)

# ==== Load and clean data ===== #
#Read the counts data from the file "fcount_all.matrix.txt" into a data frame
countdata<- read.table("fcount_all.matrix.txt", header=TRUE,row.names = 1)
head(countdata)
dim(countdata) #62700    10
#include only 9 samples in new data frame
countdata_1 <- countdata[, c(2:10)]
dim(countdata_1) #62700     9

# Filter countdata 
# Keep rows where the sum of values in each row is greater than or equal to 10
# and all values in each row are greater than or equal to 1
countdata.filter<-countdata_1[rowSums(countdata_1)>=10&apply(countdata_1,1,function(x){all(x>=1)}),]
dim(countdata.filter) 
#16522     9 
# 16522/62700 = 26.35% 

colnames(countdata.filter) <- c("F1","F2","F3", "L1", "L2", "L3", "V1","V2", "V3")
head(countdata.filter)
rownames(countdata.filter) <- lapply(rownames(countdata.filter), sub, pattern = "\\.\\d+$", replacement = "")
head(countdata.filter)

###QC - outlier detection, use goodSampleGenes to detect outlier genes
gsg <- goodSamplesGenes(t(countdata.filter))
summary (gsg)
table(gsg$goodGenes) ## TRUE  14218
table(gsg$goodSamples) ## TRUE  9 

### QC - outlier detection, hierarchical clustering 
htree <- hclust(dist(t(countdata.filter)), method = "average")
plot(htree)  

### QC - outlier detection, PCA
pca <- prcomp(t(countdata.filter))
pca.dat <- pca$x  
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2)) + 
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

countdata_FL <-countdata.filter[,c(1:6)]

#Install ComBat-seq from GitHub
install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
library(sva)

batch <- c(rep(1, 3), rep(2, 3))
countdata.filter_com <- ComBat_seq(countdata_FL, batch=batch, group=NULL)

###normalization, use Deseq 2 creat dds, not spcifying model
library(DESeq2)
condition <- factor(c(rep("F",3),rep("L",3)))
coldata <- data.frame(row.names=colnames(countdata.filter_com), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata.filter_com, colData=coldata, design=~condition) 

##perform variance stabilization 
dds_norm <- vst(dds)
norm.counts <- assay(dds_norm)


#Install ComBat-seq from GitHub
install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
library(sva)

batch <- c(rep(1, 3), rep(2, 3))
countdata.filter_com <- ComBat_seq(countdata.filter, batch=batch, group=NULL)

##network construction 
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 500, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor





