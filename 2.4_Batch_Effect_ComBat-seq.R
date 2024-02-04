
#Read the counts data from the file "fcount_all.matrix.txt" into a data frame
countdata<- read.table("fcount_all.matrix.txt", header=TRUE,row.names = 1)
head(countdata)
dim(countdata) #62700    10

#include only 9 samples in new data frame
countdata_1 <- countdata[, c(2:10)]
dim(countdata_1) #62700     9

# Keep rows where the sum of values in each row is greater than or equal to 10
# and all values in each row are greater than or equal to 1
countdata.filter<-countdata_1[rowSums(countdata_1)>=10&apply(countdata_1,1,function(x){all(x>=1)}),]
dim(countdata.filter) 
# 16522      9 
# 16522/62700 = 26.35% 

colnames(countdata.filter) <- c("F1","F2","F3", "L1", "L2", "L3", "V1","V2", "V3")

## remove the decimal part (.\\d+$) from the end of each row name 
## remove version numbers to keep only the core Ensembl gene IDs.
rownames(countdata.filter) <- lapply(rownames(countdata.filter), sub, pattern = "\\.\\d+$", replacement = "")
head(countdata.filter)

library(AnnotationDbi)
library(org.Hs.eg.db)
countdata.filter$symbol <- mapIds(org.Hs.eg.db, keys = rownames(countdata.filter), keytype = "ENSEMBL", column = "SYMBOL")


sum(is.na(countdata.filter$symbol)) #2304
countdata.filter <- countdata.filter[complete.cases(countdata.filter$symbol), ]
dim(countdata.filter)  ## 14218    10

## use Ensembles ID_gene symbols as the new row name
countdata.filter <- countdata.filter[,c(10,1:9)]
new_row_names <- paste(rownames(countdata.filter), countdata.filter[, 1], sep = "_")
rownames(countdata.filter) <- new_row_names
countdata.filter <- countdata.filter[, 2:10]


#Install ComBat-seq from GitHub
install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
library(sva)

batch <- c(rep(1, 6), rep(2, 3))
countdata.filter_com <- ComBat_seq(countdata.filter, batch=batch, group=NULL)

##prepares a DESeqDataSet object (dds) for further analysis
library(DESeq2)
condition_com <- factor(c(rep("F",3),rep("L",3),rep("V",3))) 
coldata_com <- data.frame(row.names=colnames(countdata.filter_com), condition_com)
dds_com <- DESeqDataSetFromMatrix(countData=countdata.filter_com, colData=coldata_com, design=~condition_com) 

#plot a PCA
library(ggplot2)
rld_com <- rlogTransformation(dds_com)
head(assay(rld_com))
g <- DESeq2::plotPCA(rld_com, intgroup="condition_com")
g + coord_fixed(ratio = 1) + theme_bw()

#plot heatmap
library(RColorBrewer)
library(gplots)

class(condition_com) #factor
condition_com<- data.frame(condition_com)

(mycols <- brewer.pal(8, "Accent")[1:length(unique(condition_com$condition))])
sampleDists <- as.matrix(dist(t(assay(rld_com))))

pdf("heatmap_plot_com.pdf", width = 8, height = 8)
heatmap.2(as.matrix(sampleDists), key=T, trace="none",
          col=colorpanel(100, "#4575B4", "white"),
          ColSideColors=mycols[condition_com$condition],
          RowSideColors=mycols[condition_com$condition],
          main = "Sample Distance Matrix",)
dev.off()








