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


############# limma ###################
countdata.filter_FL 
install.packages("limma")
library(limma)

# Create a Design Matrix:
design <- model.matrix(~0 + factor(c(rep("F", 3), rep("L", 3))))
colnames(design) <- levels(factor(c(rep("F", 3), rep("L", 3))))

# Fit a Linear Model:
fit <- lmFit(countdata.filter_FL , design)

# Contrast Analysis:
contrast_matrix <- makeContrasts(FvsL = "F - L", levels = design)

# Empirical Bayes Moderation:
fit <- eBayes(contrasts.fit(fit, contrast_matrix))

# Extract Results:
results <- topTable(fit, coef = "FvsL", number = Inf)
print(results)

## select up/down regulated genes
dim(results) # 16522 6
subset(results, P.Value < 0.05) -> diff
dim(diff) # 6515 6
subset(diff, logFC < -0.585) -> down
subset(diff, logFC > 0.585) -> up
dim(down) # 1253    6
dim(up) # 5262    6

## Examine plot of p-values, the MA plot and the Volcano Plot:
hist(
  results$adj.P.Val,
  breaks = 300,
  col = "#4575B4",
  border = NA
)


## MA plot for publish
library(ggplot2)
library(ggpubr)
rownames_FL <- as.character(rownames(results))

sum(resdata_FL$padj < 0.05, na.rm = TRUE) #5811
sum(resdata_FL$log2FoldChange > 1, na.rm = TRUE) #517
sum(resdata_FL$log2FoldChange < -1, na.rm = TRUE) #684
sum(resdata_FL$padj < 0.05 &
      resdata_FL$log2FoldChange > 1, na.rm = TRUE) #410
sum(resdata_FL$padj < 0.05 &
      resdata_FL$log2FoldChange < -1, na.rm = TRUE) #410

rownames_FL <- as.character(rownames(results))

ggmaplot(
  results,
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







