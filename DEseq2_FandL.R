########### Data Input ##########

countdata <-
  read.table("fcount_all.matrix.txt",
             header = TRUE,
             row.names = 1)
dim(countdata) #62700    10
countdata_1 <- countdata[, c(2:10)]
dim(countdata_1)

## Keep rows where the sum of values in each row is greater than or equal to 10
## and all values in each row are greater than or equal to 1
countdata.filter <-
  countdata_1[rowSums(countdata_1) >= 10 &
                apply(countdata_1, 1, function(x) {
                  all(x >= 1)
                }), ]
dim(countdata.filter)  ##16522      9
## 12073/62700 = 26.4%

colnames(countdata.filter) <-
  c("F1", "F2", "F3", "L1", "L2", "L3", "V1", "V2", "V3")
dim(countdata.filter) ##16522     9


## process with only F and L samples
countdata.filter_FL <- countdata.filter[, c(1:6)]

## add gene symbols with Ensembl gene IDs
## first, remove the decimal part (.\\d+$) from the end of each row name,
## removing those version numbers to keep only the core Ensembl gene IDs.
rownames(countdata.filter_FL) <-
  lapply(rownames(countdata.filter_FL),
         sub,
         pattern = "\\.\\d+$",
         replacement = "")

library(AnnotationDbi)
library(org.Hs.eg.db)

countdata.filter_FL$symbol <-
  mapIds(
    org.Hs.eg.db,
    keys = rownames(countdata.filter_FL),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  )

sum(is.na(countdata.filter_FL$symbol)) #2304 NAs in gene symbol column
countdata.filter_FL <-
  countdata.filter_FL[complete.cases(countdata.filter_FL$symbol), ]
dim(countdata.filter_FL)  ## 14218     7
countdata.filter_FL <- countdata.filter_FL[, c(7, 1:6)]

## use Ensembles ID_gene symbols as the new row name
new_row_names <-
  paste(rownames(countdata.filter_FL), countdata.filter_FL[, 1], sep = "_")
rownames(countdata.filter_FL) <- new_row_names
countdata.filter_FL <- countdata.filter_FL[, 2:7]
head(countdata.filter_FL)


## prepares a DESeqDataSet object (dds) for further analysis
library(DESeq2)
condition_FL <- factor(c(rep("F", 3), rep("L", 3)))
coldata_FL <-
  data.frame(row.names = colnames(countdata.filter_FL), condition_FL)
dds_FL <-
  DESeqDataSetFromMatrix(
    countData = countdata.filter_FL,
    colData = coldata_FL,
    design =  ~ condition_FL
  )

## plot a PCA
library(ggplot2)
rld_FL <- rlogTransformation(dds_FL)
head(assay(rld_FL))
g <- DESeq2::plotPCA(rld_FL, intgroup = "condition_FL")
g + coord_fixed(ratio = 2) + theme_bw()

## plot sample distance heatmap
library(RColorBrewer)
library(gplots)

class(condition_FL) #factor
condition_FL_df <- data.frame(condition_FL)

(mycols <-
    brewer.pal(8, "Accent")[1:length(unique(condition_FL_df$condition))])
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

###DESeq2 for difference expression analysis
dds_FL <- DESeq(dds_FL)
resdata_FL <-
  results(dds_FL, contrast = c("condition_FL", "F", "L"))
table(resdata_FL$padj < 0.05)  #FALSE 8407  TRUE  5811
res_padj_FL <- resdata_FL[order(resdata_FL$padj), ]

resultsNames(dds_FL) ##"Intercept"  "condition_FL_L_vs_F"
write.table(res_padj_FL,
            "FL_diffexpr_padj_results.txt",
            quote = F,
            sep = '\t')

### normalized counts table
normalized_counts_FL <-
  as.data.frame(counts(dds_FL, normalized = TRUE))

##round them to the nearest integer
normalized_counts_FL_1 <-
  as.data.frame(lapply(normalized_counts_FL, round))
rownames(normalized_counts_FL_1) <- rownames(normalized_counts_FL)
write.csv(normalized_counts_FL_1, file = "FL_normalized.csv")

## select up/down regulated genes
dim(resdata_FL) # 14218 6
subset(resdata_FL, pvalue < 0.05) -> diff
subset(diff, log2FoldChange < -0.585) -> down
subset(diff, log2FoldChange > 0.585) -> up
dim(down) # 1071    6
dim(up) # 1251    6

## Extract and write the names of up/down regulated genes
up_names <- rownames(up)
write.table(
  up_names,
  'up_gene.txt',
  quote = F,
  sep = '\t',
  row.names = F
)
down_names <- rownames(down)
write.table(
  down_names,
  'down_gene.txt',
  quote = F,
  sep = '\t',
  row.names = F
)

## Examine plot of p-values, the MA plot and the Volcano Plot:
hist(
  resdata_FL$padj,
  breaks = 300,
  col = "#4575B4",
  border = NA
)
DESeq2::plotMA(dds_FL, ylim = c(-2, 2))

## MA plot for publish
library(ggpubr)
rownames_FL <- as.character(rownames(resdata_FL))
second_part <-
  sapply(strsplit(rownames_FL, "_"), function(x)
    tail(x, 1))

sum(resdata_FL$padj < 0.05, na.rm = TRUE) #5811
sum(resdata_FL$log2FoldChange > 1, na.rm = TRUE) #517
sum(resdata_FL$log2FoldChange < -1, na.rm = TRUE) #684
sum(resdata_FL$padj < 0.05 &
      resdata_FL$log2FoldChange > 1, na.rm = TRUE) #410
sum(resdata_FL$padj < 0.05 &
      resdata_FL$log2FoldChange < -1, na.rm = TRUE) #410

rownames_FL <- as.character(rownames(resdata_FL))
second_part <-
  sapply(strsplit(rownames_FL, "_"), function(x)
    tail(x, 1))

ggmaplot(
  resdata_FL,
  main = expression("Group 1" %->% "Group 2"),
  fdr = 0.06,
  fc = 2,
  size = 1.5,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(second_part),
  legend = "top",
  top = 20,
  font.label = c("bold", 14),
  label.rectangle = TRUE,
  font.legend = "bold",
  font.main = "bold",
  ggtheme = ggplot2::theme_minimal() +
    
    theme(
      text = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
)

## Volcano plot for fun
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

## volcano plot for publish
resdata <-
  read.table(
    "FL_diffexpr_padj_results.txt",
    header = T,
    sep = '\t',
    row.names = 1
  )
new_labels <-
  sapply(strsplit(rownames(resdata)[1:20], "_"), function(x)
    x[2])
resdata$label <- c(new_labels, rep(NA, (nrow(resdata) - 20)))

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

EnhancedVolcano(resdata_FL,
                lab = second_part,
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(
  resdata_FL,
  lab = second_part,
  x = 'log2FoldChange',
  y = 'pvalue',
  selectLab = c(
    'TFRC',
    'ARPIN',
    'ARRDC4',
    'TMEM79',
    'AGK',
    'ATPAF1',
    'JOSD2',
    'MAFG',
    'GPRASP2',
    'TFT46',
    'SNHG7',
    'ZNF143',
    'ZNF22'
  ),
  xlab = bquote( ~ Log[2] ~ 'fold change'),
  pCutoff = 10e-14,
  FCcutoff = 1,
  pointSize = 4.0,
  labSize = 6.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,
  colAlpha = 4 / 5,
  legendPosition = 'right',
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black'
)

## pheatmap for publish
library(pheatmap)
library(RColorBrewer)
library(viridisLite)
library(viridis)

resdata <- resdata_FL
choose_gene = head(rownames(res_padj_FL), 30)
choose_matrix = countdata.filter_FL[choose_gene, ]
choose_gene_1 <-
  sapply(strsplit(choose_gene, "_"), function(x)
    x[2])
rownames(choose_matrix) <- choose_gene_1

col_groups <- substr(colnames(choose_matrix), 1, 1)
table(col_groups)
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(choose_matrix)

mat_colors <- list(group = c("#BEAED4", "#7FC97F"))
names(mat_colors$group) <- unique(col_groups)

png(filename = "DEG_pheatmap.png",
    width = 3000,
    height = 4000,
    res = 500)

pheatmap(choose_matrix, 
         scale ="row", 
         border_color = NA,
         drop_levels = TRUE,
         annotation_col = mat_col,
         annotation_colors = mat_colors,
         color = colorRampPalette(c("#4575B4", "white", "#d4aeb5"))(50),
         fontsize = 15)
dev.off() 





## make a mastertable of all count and deseq data
mastertable <- merge(counts_table_short, result_ordered_try, by = 0)
nrow(mastertable) #14425
nrow(counts_table_short) #15502
nrow(result_ordered_try) #14426
mastertable <- mastertable[, c(1, 8, 2:7, 9:ncol(mastertable))]
head(mastertable)
write.csv(as.data.frame(mastertable),
          file = "RNAseq_mastertable.csv")

