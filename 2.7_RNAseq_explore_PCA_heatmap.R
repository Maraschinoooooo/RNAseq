
##check featurecounts_summary
file_path <- "FeatureCounts_counts.txt.summary"
featurecounts_summary <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
str(featurecounts_summary)

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

### Optional: manual normalization
### Normalize the counts by dividing each element by the corresponding column sum
### Multiply the result by a scaling factor (20,000,000 in this case)
#column_sums <- colSums(countdata.filter, na.rm = TRUE)
#countdata.filter_1 <- countdata.filter / rep(column_sums, each = nrow(countdata.filter)) * 20000000
#head(countdata.filter_1)
#column_sums_2 <- colSums(countdata.filter_1, na.rm = TRUE)

# Convert all numeric values to integers
#countdata.filter_2 <- as.data.frame(lapply(countdata.filter_1, function(x) as.integer(as.numeric(x))))
#rownames(countdata.filter_2) <- rownames(countdata.filter_1)
#head(countdata.filter_2)
#colSums(countdata.filter_2, na.rm = TRUE)
##  F1       F2       F3       L1       L2       L3       V1       V2       V3 
##19991653 19991662 19991237 19991670 19991204 19989465 19991239 19991976 19991798 


##prepares a DESeqDataSet object (dds) for further analysis
library(DESeq2)
library(ggplot2)
condition <- factor(c(rep("F",3),rep("L",3), rep("V", 3))) 
coldata <- data.frame(row.names=colnames(countdata.filter), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata.filter, colData=coldata, design=~condition) 

#plot a PCA
rld <- rlogTransformation(dds)
head(assay(rld))
pca <- DESeq2::plotPCA(rld, intgroup="condition")
pca + coord_fixed(ratio = 5) + theme_bw()

#sample distance heatmap
library(RColorBrewer)
library(gplots)

class(condition) #factor
condition_df <- data.frame(condition)

(mycols <- brewer.pal(8, "Accent")[1:length(unique(condition_df$condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))

pdf("heatmap_plot.pdf", width = 8, height = 8)
heatmap.2(as.matrix(sampleDists), key=T, trace="none",
          col=colorpanel(100, "#4575B4", "white"),
          ColSideColors=mycols[condition_df$condition],
          RowSideColors=mycols[condition_df$condition],
          main = "Sample Distance Matrix",)
dev.off()



