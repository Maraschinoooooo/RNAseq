##check featurecounts_summary##

setwd("/Users/shaoboqin_home/Desktop/20231214_forR")
file_path <- "FeatureCounts_counts.txt.summary"
featurecounts_summary <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
str(featurecounts_summary)


########### Data Input ##########

countdata <-
  read.table("fcount_all.matrix.txt",
             header = TRUE,
             row.names = 1)
head(countdata)
dim(countdata) #62700    10


############  TPM (transcripts per million)  #######

# Divide numeric columns by the values in the Length column and times 2000
countdata[, 2:10] <- countdata[, 2:10] / countdata$Length * 2000
countdata <- countdata[, 2:10]

### Normalize the counts by dividing each element by the corresponding column sum
### Multiply the result by a scaling factor (1,000,000 in this case)
column_sums <- colSums(countdata, na.rm = TRUE)
countdata_1 <- countdata / rep(column_sums, each = nrow(countdata)) * 1000000
head(countdata_1)
column_sums_2 <- colSums(countdata_1, na.rm = TRUE)

# Convert all numeric values to integers
countdata_2 <- as.data.frame(lapply(countdata_1, function(x) as.integer(as.numeric(x))))
rownames(countdata_2) <- rownames(countdata_1)

colSums(countdata_2, na.rm = TRUE)
 
### F1_sort.bam          F3_sort.bam    F_backup_sort.bam 
### 989636               989188               988523 
### L1_sort.bam          L2_sort.bam    L_backup_sort.bam 
### 988862               988433               987687 
### SRR20810434_sort.bam SRR20810435_sort.bam SRR20810436_sort.bam 
### 988237               988272               988448 

###################filter data #####################

# Keep rows where the sum of values in each row is greater than or equal to 10
# and all values in each row are greater than or equal to 1
countdata.filter<-countdata_2[rowSums(countdata_2)>=10&apply(countdata_2,1,function(x){all(x>=1)}),]
dim(countdata.filter)  ##12073    9   
## 12073/62700 = 19.3%

colnames(countdata.filter) <- c("F1","F2","F3", "L1", "L2", "L3", "V1","V2", "V3")
head(countdata.filter)


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
pca + coord_fixed(ratio = 6) + theme_bw()

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







