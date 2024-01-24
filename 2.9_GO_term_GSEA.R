# Load required libraries
library(tidyr)
library(ggpubr) 

# Section 1: Prepare files for GO term analysis and GSEA analysis
data <- read.table("updown_genes_up.txt", header = TRUE, sep = "\t")

# Separate ENSG and Gene columns
data_split <- separate(data, x, into = c("ENSG", "Gene"), sep = "_")
ENSG <- data_split[,1]
Gene <- data_split[,2]

# Write ENSG and Gene files
write.table(ENSG, file = 'GSEA_up_ENSG.txt', sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Gene, file = 'GSEA_up_Gene.txt', sep = "\t", quote = FALSE, row.names = FALSE)

# Section 2: Preprocess normalized counts data
dim(normalized_counts_FL_1) #14218 6
row_names_first_part <- sapply(strsplit(row.names(normalized_counts_FL_1), "_"), function(x) x[1])
rownames(normalized_counts_FL_1) <- row_names_first_part

# Write preprocessed data
write.table(normalized_counts_FL_1, file = 'normalized_counts_FL_1.txt', sep = "\t", quote = FALSE, row.names = TRUE)

# Section 3: Run GSEA using external software and generate bar plot
GSEA_bar <- data.frame(
  geneset = c("KRAS_SIGNALING_UP", "EPITHELIAL_MESENCHYMAL_TRANSITION", "ADIPOGENESIS", "UV_RESPONSE_DN",
              "IL2_STAT5_SIGNALING", "E2F_TARGETS", "ANGIOGENESIS", "BILE_ACID_METABOLISM", "PANCREAS_BETA_CELLS",
              "DNA_REPAIR", "COAGULATION", "GLYCOLYSIS", "TNFA_SIGNALING_VIA_NFKB", "P53_PATHWAY",
              "CHOLESTEROL_HOMEOSTASIS", "UV_RESPONSE_UP", "MYC_TARGETS_V2", "INFLAMMATORY_RESPONSE",
              "WNT_BETA_CATENIN_SIGNALING", "MTORC1_SIGNALING", "HYPOXIA", "NOTCH_SIGNALING",
              "UNFOLDED_PROTEIN_RESPONSE", "COMPLEMENT"),
  NES = c(-1.3855982, -1.3686199, -1.3094803, -1.2114127, -1.2032882, -1.1646345, -1.1610544, -1.1572849,
          -1.1567018, -1.1348317, -1.1245372, -1.1187463, 1.6034037, 1.551731, 1.3539077, 1.2314007, 1.212223,
          1.1846715, 1.1681427, 1.1423787, 1.1201321, 1.0874139, 1.0867114, 1.080289)
)

# Add a group column based on NES values
GSEA_bar$grp <- factor(ifelse(GSEA_bar$NES < 0, "enriched in minimal medium", "depleted in minimal medium"),
                       levels = c("enriched in minimal medium", "depleted in minimal medium"))

# Save bar plot to PDF
pdf("GSEA_bar.pdf", width = 8, height = 8)
ggbarplot(GSEA_bar, x = "geneset", y = "NES",
          fill = "grp",               # change fill color by grp
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palette. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x-axis texts
          ylab = "normalized enrichment score (NES)",
          xlab = FALSE,
          legend.title = "geneset MsigDB Hallmarks",
          rotate = TRUE,
          ggtheme = theme_minimal()
)
dev.off()
