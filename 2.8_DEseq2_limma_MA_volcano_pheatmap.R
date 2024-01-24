########### Data Input ##########

# Read count data
countdata <- read.table("fcount_all.matrix.txt", header = TRUE, row.names = 1)
head(countdata)
dim(countdata) #62700 10

# Select columns 2 to 10
countdata_1 <- countdata[, c(2:10)]
dim(countdata_1)

# Filter rows based on sum and values
countdata.filter <- countdata_1[rowSums(countdata_1) >= 10 & apply(countdata_1, 1, function(x) all(x >= 1)), ]
dim(countdata.filter)  ##16522 9

# Rename columns
colnames(countdata.filter) <- c("F1", "F2", "F3", "L1", "L2", "L3", "V1", "V2", "V3")
head(countdata.filter)
dim(countdata.filter) ##16522 9

# Process with only F and L samples
countdata.filter_FL <- countdata.filter[, c(1:6)]
head(countdata.filter_FL)

# Remove decimal part from row names
rownames(countdata.filter_FL) <- lapply(rownames(countdata.filter_FL), sub, pattern = "\\.\\d+$", replacement = "")
head(countdata.filter_FL)

# Add gene symbols with Ensembl gene IDs
library(AnnotationDbi)
library(org.Hs.eg.db)

countdata.filter_FL$symbol <- mapIds(org.Hs.eg.db, keys = rownames(countdata.filter_FL), keytype = "ENSEMBL", column = "SYMBOL")

# Remove NAs in gene symbol column
countdata.filter_FL <- countdata.filter_FL[complete.cases(countdata.filter_FL$symbol), ]
dim(countdata.filter_FL)  ## 14218 7
countdata.filter_FL <- countdata.filter_FL[, c(7, 1:6)]

# Use Ensembl IDs_gene symbols as new row names
new_row_names <- paste(rownames(countdata.filter_FL), countdata.filter_FL[, 1], sep = "_")
rownames(countdata.filter_FL) <- new_row_names
countdata.filter_FL <- countdata.filter_FL[, 2:7]
head(countdata.filter_FL)

## Prepares a DESeqDataSet object (dds) for further analysis
library(DESeq2)
condition_FL <- factor(c(rep("F", 3), rep("L", 3)))
coldata_FL <- data.frame(row.names = colnames(countdata.filter_FL), condition_FL)
dds_FL <- DESeqDataSetFromMatrix(countData = countdata.filter_FL, colData = coldata_FL, design = ~condition_FL)

# Plot a PCA
rld_FL <- rlogTransformation(dds_FL)
head(assay(rld_FL))
g <- DESeq2::plotPCA(rld_FL, intgroup = "condition_FL")
g + coord_fixed(ratio = 2) + theme_bw()

# Plot sample distance heatmap
library(RColorBrewer)
library(gplots)

condition_FL_df <- data.frame(condition_FL)
(mycols <- brewer.pal(8, "Accent")[1:length(unique(condition_FL_df$condition))])
sampleDists <- as.matrix(dist(t(assay(rld_FL))))

pdf("heatmap_plot_FL.pdf", width = 8, height = 8)
heatmap.2(
  as.matrix(sampleDists),
  key = T,
  trace = "none",
  col = colorpanel(100, "#4575B4", "white"),
  ColSideColors = mycols[condition_FL_df$condition],
  RowSideColors = mycols[condition_FL_df$condition],
  main = "Sample Distance Matrix",
)
dev.off()

### DESeq2 for differential expression analysis
dds_FL <- DESeq(dds_FL)
resdata_FL <- results(dds_FL, contrast = c("condition_FL", "F", "L"))
table(resdata_FL$padj < 0.05)  # FALSE 8407 TRUE 5811
res_padj_FL <- resdata_FL[order(resdata_FL$padj), ]

resultsNames(dds_FL) ## "Intercept" "condition_FL_L_vs_F"
write.table(res_padj_FL, "FL_diffexpr_padj_results.txt", quote = F, sep = '\t')

### Normalized counts table
normalized_counts_FL <- as.data.frame(counts(dds_FL, normalized = TRUE))

## Round them to the nearest integer
normalized_counts_FL_1 <- as.data.frame(lapply(normalized_counts_FL, round))
rownames(normalized_counts_FL_1) <- rownames(normalized_counts_FL)
write.csv(normalized_counts_FL_1, file = "FL_normalized.csv")

## Select up/downregulated genes
dim(resdata_FL) # 14218 6
subset(resdata_FL, pvalue < 0.05) -> diff 
subset(diff, log2FoldChange < -0.585) -> down 
subset(diff, log2FoldChange > 0.585) -> up  
dim(down) # 1071 6
dim(up) # 1251 6

## Extract and write the names of up/downregulated genes 
up_names <- rownames(up)
write.table(up_names, 'up_gene.txt', quote = F, sep = '\t', row.names = F)
down_names <- rownames(down)
write.table(down_names, 'down_gene.txt', quote = F, sep = '\t', row.names = F)

# Examine plot of p-values, the MA plot and the Volcano Plot:
hist(resdata_FL$padj, breaks = 300, col = "#4575B4", border = NA)
DESeq2::plotMA(dds_FL, ylim = c(-1, 1))

# Volcano plot
with(resdata_FL,
     plot(
       log2FoldChange,
       -log10(pvalue),
       pch = 20,
       main = "Volcano plot",
       xlim = c(-2.5, 2)
     ))
with(subset(resdata_FL, padj < .0005),
     points(
       log2FoldChange,
       -log10(pvalue),
       pch = 20,
       col = "red"
     ))

####### Volcano plot for publish #####
resdata <- read.table("FL_diffexpr_padj_results.txt", header = T, sep = '\t', row.names = 1)
new_labels <- sapply(strsplit(rownames(resdata)[1:20], "_"), function(x) x[2])
resdata$label <- c(new_labels, rep(NA, (nrow(resdata) - 20)))

library(ggplot2)

pdf("volcano_plot.pdf", width = 9, height = 10)

ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(
    aes(size = -log10(padj), color = -log10(padj)),
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 1),
    alpha = 0.6
  ) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_size_continuous(range = c(0, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.08, 0.9),
    legend.justification = c(0, 1)
  ) +
  guides(col = guide_colourbar(title = "-Log10_q-value"), size = "none") +
  geom_text(
    aes(label = c(label), color = -log10(padj)),
    size = 3,
    vjust = 1.5,
    hjust = 1,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.3)
  ) +
  xlab("Log2FC") + ylab("-Log10(FDR q-value)")

dev.off()

####### Heatmap ######
resdata <- read.table("FL_diffexpr_padj_results.txt", header = T, sep = '\t', row.names = 1)

library(pheatmap) 
choose_gene <- head(rownames(res_padj_FL), 50)  
choose_matrix <- countdata.filter[choose_gene, ]  
choose_matrix <- t(scale(t(choose_matrix)))  

choose_gene <- head(rownames(res_padj_FL), 100)  
choose_gene_1 <- sapply(strsplit(choose_gene, "_"), function(x) x[2])
choose_matrix <- countdata.filter_FL[choose_gene_1, ]  
choose_matrix <- t(scale(t(choose_matrix)))  

png(filename = "DEG_pheatmap.png", width = 600, height = 1000)
pheatmap(choose_matrix)
dev.off()

# Make a mastertable of all count and deseq data
mastertable <- merge(counts_table_short, result_ordered_try, by = 0)
nrow(mastertable) #14425
nrow(counts_table_short) #15502
nrow(result_ordered_try) #14426
mastertable <- mastertable[, c(1, 8, 2:7, 9:ncol(mastertable))]
head(mastertable)
write.csv(as.data.frame(mastertable), file = "RNAseq_mastertable.csv")

# Pheatmap
library(pheatmap)

condition_Low_vs_Full_0.001 <- subset(result_ordered_try, padj < 0.001)
nrow(condition_Low_vs_Full_0.001) #2320
head(condition_Low_vs_Full_0.001)

allsig <- merge(counts_table_short, condition_Low_vs_Full_0.001, by = 0)
sigcount <- allsig[, 2:8]
head(sigcount)
sigcount <- na.omit(sigcount)
row.names(sigcount) <- sigcount$symbol
ncol(sigcount) #7
sigcount <- sigcount[, 1:6]

pheatmap(
  log2(sigcount + 1),
  scale = "row",
  show_rownames = F,
  treeheight_row = 0
)
