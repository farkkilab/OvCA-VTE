# RNA-seq Analysis ovca-vte project
# Author(s): Sara Palomino

#### Load packages ##### 
library(DESeq2)
library(GSVA)
library(tximport)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(AnnotationDbi)
library(Biobase)
library(org.Hs.eg.db)
library(enrichplot)
library(ggnewscale)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(openxlsx)


#### PREPROCESSING #### 
#### 01.Load data and metadata
setwd("~/00.Inputs/")
dir_path("~/00.Inputs/RNA_samples")

metadata <- readRDS("metadata.rds")
sample_ids <- metadata$RNA_sample

files_RNA <- file.path(dir_path, sample_ids , "quants.genes.sf")

# Import salmon quantification file
txi_gene <- tximport(files_RNA, type="salmon", txOut=TRUE)

# store the count table
count_table  <- as.data.frame(txi_gene$counts)
names(count_table) <- sample_ids


# transform ENSEMBL into SYMBOL and ENTREZID
geneid <- rownames(count_table)

genes <- select(org.Hs.eg.db, keys=geneid, columns=c("SYMBOL", "ENTREZID"), keytype="ENSEMBL")

genes <- genes[!is.na(genes$SYMBOL)] # keep only those genes with symbol

raw_counts_symbol <- merge(genes, count_table, by.x="ENSEMBL", by.y="row.names")

rownames(raw_counts_symbol) <- raw_counts_symbol$ENSEMBL

#### 02.DESEQ WORKFLOW
#### 02.1 >>> Build  and explore data object
# table of sample information
sampleTable <- data.frame(condition=as.factor(metadata$Group), batch=as.factor(metadata$batch), anticoag=as.factor(metadata$preop_anticoag))
rownames(sampleTable) <- metadata$RNA_sample

# Build the DESeq object
dds <- DESeqDataSetFromTximport(txi_gene, sampleTable, ~batch+ anticoag+ condition ) # we want to measure the effect of the condition, controlling for batch differences
txi_counts <- counts(dds)


# Here we will store the raw counts
txi_counts_filtered <- counts(dds)

# How many reads were sequenced for each sample ( = library sizes)?
colSums(counts(dds)) %>% barplot

dim(dds)

#### >>> Prefiltering
# Remove genes with no reads
keep_genes <- rowSums(counts(dds)) > 0
dds <- dds[ keep_genes, ]
dim(dds)

# Remove genes that were not expressed
smallestGroupSize <- 7
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # 18975


#### 02.3 >>> Normalizing for sequencing depth and RNA composition differences
## define a function to calculate the geometric mean
gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }

## calculate the geometric mean for each gene using that function
## note the use of apply(), which we instruct to apply the gm_mean()
## function per row (this is what the second parameter, 1, indicates)
pseudo_refs <- counts(dds) %>% apply(., 1, gm_mean)

## divide each value by its corresponding pseudo-reference value
pseudo_ref_ratios <- counts(dds) %>% apply(., 2, function(cts){ cts/pseudo_refs})

## if you want to see what that means at the single-gene level,
## compare the result of this:
counts(dds)[1,]/pseudo_refs[1]

## with
pseudo_ref_ratios[1,]

## determine the median value per sample to get the size factor
sf <- apply(pseudo_ref_ratios , 2, median)

# Assign the manually calculated size factors to the DESeqDataSet object
sizeFactors(dds) <- sf # add SFs  to object

sizeFactors(dds) # Check if the size factors were correctly added

## extracting normalized counts
counts.sf_normalized <- counts(dds, normalized=TRUE)

pdf("Count_normalization.pdf", height = 5)
par(mfrow=c(1,2))
## bp of non-normalized
boxplot(counts(dds), notch=TRUE, main = "Raw counts", cex.main = 1,
        ylab="Read counts", cex = .6,
        xaxt="n")

## bp of size-factor normalized values
boxplot(log2(counts(dds, normalize= TRUE) +1), notch=TRUE,
        main = "Size-factor-normalized counts", cex.main = 1,
        ylab="log2(read counts)", cex = .6,
        xaxt="n")
dev.off()


#### UNTARGETED ANALYSIS #### 
    #### >>> 01. Run DESeq2 analysis (DEA) ####
dds <- DESeq(dds)

# Get results from DEA
res <- results(dds, contrast = c("condition", "1", "0"))

resultsNames(dds)

res_df <- as.data.frame(res)

dim(na.omit(res_df[res_df$padj<0.05,]))

res_df$diffexpressed <- "Not significative"
res_df$diffexpressed[res_df$padj<0.05] <- "Adj.P-value"
res_df$diffexpressed[res_df$padj<0.05 & res_df$log2FoldChange>0.58] <- "Adj.P-value & Log2FC +"
res_df$diffexpressed[res_df$padj<0.05 & res_df$log2FoldChange<0.58] <- "Adj.P-value & Log2FC -"


table(res_df$diffexpressed)

res_df <- merge(res_df, raw_counts_symbol[,c(1:3)], by.x="row.names", by.y="ENSEMBL")
res_df_sig <- res_df[!res_df$diffexpressed=="Not.sig",]


###     rlog (Regularized Log Transformation):
rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)


# PCA 
mm <- model.matrix(~condition, colData(rld))
mat <- limma::removeBatchEffect(mat, batch=rld$batch, design=mm)

assay(rld) <- mat

### PCA with DEGs
groups <- as.factor(sampleTable$condition) 
batch <- as.factor(sampleTable$batch) 
# Extract the rlog-transformed counts matrix
rlog_counts <- as.data.frame(rld@assays@data)[,-c(1:2)]
metadata$Group2 <- as.character(metadata$Group)

metadata$Group2[metadata$Group=="1"] <- "VTE"
metadata$Group2[metadata$Group=="0"] <- "Non VTE"

select <- na.omit(res_df[res_df$padj<0.05,]) # significative genes
names(select)[1] <- "ENSEMBL"

write.xlsx(select, 'untargeted_DEA.xlsx')


pca <- prcomp(t(rlog_counts[select$Row.names,]), center = TRUE, scale. = TRUE) # Computes PCA on the dataset

var_prcomp <- pca$sdev^2 # Variance of each principal component (squared standard deviation)

pca_variance_explained <- round(100 * var_prcomp / sum(var_prcomp), 1)  # Convert to %

# Create a data frame showing the proportion of variance explained by each PC
pcvar <- data.frame(
  var = var_prcomp / sum(var_prcomp),  # Normalize variance to get proportion of variance explained
  pc = c(1:length(var_prcomp))  # Create a numeric sequence for principal components
)

# Create a data frame for PCA plot, using the first two principal components (PC1 and PC2)
mds <- data.frame(
  PC1 = pca$x[,1],  # First principal component scores
  PC2 = pca$x[,2],  # Second principal component scores
  #Batch = metadata$batch, 
  Condition = as.factor(metadata$Group2)   # Add grouping information for color mapping
)

# Generate PCA scatter plot using ggplot2 
pca_plot_DEGS <- ggplot(mds, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size=2, alpha=1) +
  scale_color_manual(values = c("orange",  "blue"))+
  labs(
    x = paste0("PC1 (", pca_variance_explained[1], "%)"),
    y = paste0("PC2 (", pca_variance_explained[2], "%)")) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"), legend.text = element_text(size = 9),  # Adjusts the size of the text in the legend
    legend.key.size = unit(0.6, "cm"), legend.title = element_blank() ,legend.position = "bottom" ,legend.margin = margin(t = -10, b = 0)) +
  guides(color = guide_legend(
    override.aes = list(size = 3)  # Adjusts the size of the points in the legend
  )) 


### Heatmap
Group <- c(rep("VTE",7),rep( "Non_VTE",13))
Group_colors <- c(VTE="blue", Non_VTE="orange")

sig_genes <- rlog_counts[select$Row.names,]
sig_genes_scaled <- as.data.frame(scale(sig_genes))
rownames(sig_genes_scaled) <- select$SYMBOL

column_ha <- HeatmapAnnotation(Condition=as.factor(Group), col= list(Condition=Group_colors),show_annotation_name = F, simple_anno_size = unit(0.3, "cm"), show_legend = F)

ht_untargeted <- Heatmap(sig_genes_scaled, name = "Z-score", cluster_columns = F, row_names_gp = gpar(fontsize =0),         
              column_names_gp =  gpar(fontsize = 0),  
              top_annotation = column_ha ,
              show_heatmap_legend =T,
              show_row_dend = F,
              heatmap_legend_param  = 
                
                list( 
                  grid_width=unit(0.3, "cm"),
                  legend_height = unit(3 ,"cm"),
                  labels_gp = gpar(fontsize = 8),
                  title_gp = gpar(fontsize = 8, font=2),         
                  title_position = "lefttop-rot"

                ))
              

grob_untargeted = grid.grabExpr(draw( ht_untargeted, heatmap_legend_side = "left"))  
ggarrange(grob_untargeted)


###  Volcano plot
# Filter top DEGs (e.g., by top 10 genes with the smallest adj. p-value)
top_deg <- res_df %>% 
  arrange(padj) %>%
  head(20)  # Change 10 to the number of genes you want to label
top_deg <- top_deg[-3,]

volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed)) +
  geom_point(size=0.8, alpha=0.5) + 
  guides(color = guide_legend(override.aes = list(size = 2)))+
  geom_vline(xintercept = c(-0.58, 0.58), linetype="dotted", color = "black", size=0.5) +
  geom_hline(yintercept = c(log10(0.05)*(-1)), linetype="dotted", color = "black", size=0.5) +
  theme_minimal() +
  scale_color_manual(values=c(`Not sig`="darkgrey",`Adj.P-value`="chartreuse4",`Log2FC`="darkorange",`Adj.P-value & Log2FC -` ="darkblue",`Adj.P-value & Log2FC +` ="red2")) +
  
  xlim(c(-10, 8)) + 
  ylim (c(0,4))+
  labs(x = "Log2 Fold Change", y = "-log10 Adjusted P-value", color = "DEGs") + 
  
  theme(
    legend.text =   element_text(size=8),
    legend.key.size = unit(0.3, "cm"), 
    legend.title = element_blank()  ,
    axis.text.x = element_text(hjust = 2, size=8),
    axis.text.y = element_text(size=8),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8), legend.position = "bottom", legend.margin = margin(t = -10, b = 0)) +
  geom_text(data = top_deg, aes(label = top_deg$SYMBOL), size =2.3 , vjust = -0.5) +
  geom_text(data = data.frame(x = c(-9, 7), y = c(4, 4), label = c("Non VTE", "VTE")),
            aes(x = x, y = y, label = label), size = 3, color = "black", fontface="bold")

volcano

    #### >>> GSEA ####
set.seed(1234)

# Rank genes by log2FoldChange
res_df_sig <- res_df[res_df$padj<0.05,]

ranked_genes <- res_df$log2FoldChange
names(ranked_genes) <- (res_df$ENTREZID)
ranked_genes <- sort(ranked_genes, decreasing = T)

head(ranked_genes)

# GO Gene Set Enrichment Analysis
go_results <-  gseGO(geneList     = (ranked_genes),
                     OrgDb        = "org.Hs.eg.db",
                     ont          = "BP",
                     minGSSize    = 2,
                     pvalueCutoff = 1, maxGSSize = 500)


go_results_df <- as.data.frame(go_results@result)
rownames(go_results_df) <- 1:dim(go_results_df)[1]

# From the list of the top 50 significative pathways, we select some of them to be visualized
categories <- c(go_results@result$Description[1],go_results@result$Description[2],
                go_results@result$Description[3],go_results@result$Description[6],
                go_results@result$Description[7],go_results@result$Description[8],
                go_results@result$Description[9], go_results@result$Description[11],
                go_results@result$Description[18], go_results@result$Description[20],
                go_results@result$Description[21],go_results@result$Description[25],
                go_results@result$Description[41], go_results@result$Description[42],
                go_results@result$Description[45])

# reorder the categories to group them based on the tree plot
categories <- categories[c(12,2,14,10,11,7,15,5,9,13,6,4,8,3,1)]

# We generate a data.frame containing those selected pathways
rownames(go_results_df) <- go_results_df$Description
go_results_df_plot <- go_results_df[categories,]

go_results_df_plot$Description <- factor(categories, levels=rev(custom_order))

go_plot <- ggplot(go_results_df_plot, aes(NES, Description)) + 
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = setSize)) +
  scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE), name = "Adj.P-value") +
  scale_size_continuous(range = c(2,8), name = "Gene Count", breaks = c(50,100,150), labels= c("50","100", "150")) +
  theme_minimal() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL) + 
  ggtitle("   ")+
  scale_y_discrete(position = "right") +
  theme(
    legend.text = element_text(size = 8),
    legend.title  = element_text(size = 8, face="bold"),
    legend.key.height = unit(0.3, 'cm'), #change legend key height
    legend.key.width = unit(0.5, 'cm'), #change legend key
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(hjust = 2, size=8),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    legend.position = "bottom",  # position to the right
    legend.justification = c(0, 0.5)
  )

# We plot the treeplot to visualize the categories in which the pathways are grouped.
edox2 <- pairwise_termsim(go_results)
treeplot(edox2, showCategory=categories)

# We also create a barplot, based on the categories derived from the treeplot, to add it to the go_plot
cluster <- c(c(rep("cell-cell adhesion", 4), c(rep("antigen receptor signaling pathway",3)),
               c(rep("blood regulation | transmembrane transport",5)), 
               c(rep("immune response", 5))))

bar_plot <- data.frame(Value=as.numeric(as.factor(cluster)), Category=cluster)
bar_plot$count <- ave(bar_plot$Value, bar_plot$Category, FUN = length)

# Count the occurrences of each category
category_count <- table(bar_plot$Category)

# Create a vector with the count for each category in 'bar_plot'
bar_plot$count <- category_count[bar_plot$Category]

bar_plot$Category <- factor(bar_plot$Category, 
                            levels = c("cell-cell adhesion", 
                                       "antigen receptor signaling pathway", 
                                       "immune response", 
                                       "blood regulation | transmembrane transport"))


# Create the plot
bar <- ggplot(bar_plot, aes(x = 1, y = count, fill = Category)) + 
  geom_bar(stat = "identity", width = 0.1) + 
  theme_void() + 
  theme(
    legend.text = element_text(size = 8),
    legend.title  = element_text(size = 8, face="bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "bottom",   # Remove legend
    legend.direction = "vertical",
    axis.text.x = element_blank(),
    axis.title = element_blank(),
  ) +
  scale_fill_manual(values = c("cell-cell adhesion" = "lightpink", 
                               "antigen receptor signaling pathway" = "lightgreen",
                               "immune response" = "darkred",
                               "blood regulation | transmembrane transport" = "skyblue")) 

# Combine the two plots.
gsea_plot <- ggarrange(go_plot, bar, widths = c(1,0.025))

#### TARGETED ANALYSIS #### 
genes <- c("FGB",	"FGG",	"F2",	"F3",	"F5",	"F7",	"F8",	"F9",	"F10",	"F11",	"F12",	"F13A1",	"F13B",	"VWF",	"PROC",	"PROS",	"SERPINC1",
           "KLKB1",	"TFPI",	"SERPIND1",	"KNG1", "THBD",	"PROCR",	"ADAMTS13",	"TFPI2", "PLG",	"CPB2",	"PLAT",	"ANXA2",	"PLAU", "PLAUR",
           "SERPINE1",	"SERPINB2",	"SERPINA5",	"SERPINF2",	"A2M",	"SERPINE2",	"SERPINI1",	"SERPINA1",	"SERPING1",	"SELP",	"PF4",	"TBXA2R",
           "CLEC1B",	"CD40LG",	"PECAM1",	"ITGA2B",	"ITGB3",	"F2R",	"F2RL3",	"FN1",	"VTN",	"PTAFR")

normalized_counts_symbol <- merge(normalized_counts, raw_counts_symbol[,c(1:3)], by.x="row.names", by.y="ENSEMBL")

norm_counts <- normalized_counts_symbol[normalized_counts_symbol$SYMBOL %in% genes, ]
targeted_genes <- norm_counts$SYMBOL

group1_samples <- (metadata$RNA_sample[metadata$Group == "0" ])  # Modify "Group1"
group2_samples <- (metadata$RNA_sample[metadata$Group == "1" ])  # Modify "Group2"

# Ensure sample names match column names in `norm_counts`
group1_samples <- intersect(group1_samples, colnames(norm_counts))
group2_samples <- intersect(group2_samples, colnames(norm_counts))

# Apply Wilcoxon test to each gene
wilcox_results <- apply(norm_counts[,-c(1,22,23)], 1, function(gene_expr) {
  group1_values <- gene_expr[group1_samples]
  group2_values <- gene_expr[group2_samples]
  
  # Run Wilcoxon test only if there are enough non-zero values
  if (length(unique(group1_values)) > 1 & length(unique(group2_values)) > 1) {
    return(wilcox.test(group1_values, group2_values)$p.value)
  } else {
    return(NA)  # Return NA if all values are the same
  }
})


# Convert to data frame
wilcox_df <- data.frame(Gene = targeted_genes, P_value = wilcox_results)

# Adjust for multiple testing (FDR correction)
wilcox_df$adj_P_value <- p.adjust(wilcox_df$P_value, method = "BH")

# View results
wilcox_df[order(wilcox_df$adj_P_value), ]


# Compute logFC based on median expression
logFC_values <- apply(norm_counts[,-c(1,22,23)], 1, function(gene_expr) {
  # Get expression values for each group
  group1_values <- gene_expr[group1_samples]
  group2_values <- gene_expr[group2_samples]
  
  # Compute median expression for each group
  median_group1 <- mean(group1_values, na.rm=TRUE)
  median_group2 <- mean(group2_values, na.rm=TRUE)
  
  # Compute log2 fold change (avoid log(0) issue by adding small value)
  log2FC <- log2((median_group2) / (median_group1))
  
  return(log2FC)
})

# Compute FC based on median expression
FC_values <- apply(norm_counts[,-c(1,22,23)], 1, function(gene_expr) {
  # Get expression values for each group
  group1_values <- gene_expr[group1_samples]
  group2_values <- gene_expr[group2_samples]
  
  # Compute median expression for each group
  median_group1 <- mean(group1_values, na.rm=TRUE)
  median_group2 <- mean(group2_values, na.rm=TRUE)
  
  # Compute log2 fold change (avoid log(0) issue by adding small value)
  FC <- exp(log2((median_group2) / (median_group1)))
  
  return(FC)
})

# Compute % change based on median expression
change_values <- apply(norm_counts[,-c(1,22,23)], 1, function(gene_expr) {
  # Get expression values for each group
  group1_values <- gene_expr[group1_samples]
  group2_values <- gene_expr[group2_samples]
  
  # Compute median expression for each group
  median_group1 <- mean(group1_values, na.rm=TRUE)
  median_group2 <- mean(group2_values, na.rm=TRUE)
  
  # Compute log2 fold change (avoid log(0) issue by adding small value)
  change <- ((median_group2) / (median_group1))-1
  
  return(change)
})

# Convert results to DataFrame
wilcox_df$log2FC <- logFC_values
wilcox_df$FC <- FC_values
wilcox_df$change <- change_values

# View results sorted by significance
wilcox_df <- wilcox_df[order(wilcox_df$adj_P_value), ]
print(wilcox_df)
rownames(wilcox_df) <- wilcox_df$Gene

write.xlsx(wilcox_df, 'targeted_DEA.xlsx')

### PLOT a HEATMAP WITH THE TENDENCIES.
sig_mat <- wilcox_df %>%
  mutate(sig = case_when(
    adj_P_value <= 0.05 ~ 1,
    TRUE ~ 0  # TRUE is the default case when no other conditions are met
  ))

sig_mat <- data.frame(sig=sig_mat[,7], row.names = wilcox_df$Gene)
mat <- data.frame(`VTE vs no VTE`=wilcox_df[,6], row.names = wilcox_df$Gene)

legend <- Legend(
  labels = c("* = Adj. p-value <0.05"),
  labels_gp = gpar(fontsize = 9,font=2))

ht_FC <- Heatmap(mat,    name = "Expression Variation",
                 row_names_gp = gpar(fontsize =6.5),         
                 col = colorRamp2(c(-0.2, 0, 0.2), c("green", "white", "darkorange")),
                 column_names_gp =  gpar(fontsize = 0),  
                 show_row_dend=F,  
                 cell_fun = function(j, i, x, y, width, height, fill) {          
                   if (sig_mat[i, j] == 1) {            
                     grid.points( x, y, pch = 8, size = unit(1.5, "mm")) }},
                 heatmap_legend_param = list( 
                   grid_width=unit(0.3, "cm"),
                   legend_height = unit(3 ,"cm"),
                   title_position = "lefttop-rot",
                   labels_gp = gpar( fontsize = 8),title_gp = gpar(fontsize = 8, font=2))
)

ht_FC

grob_FC = grid.grabExpr(draw( ht_FC,  heatmap_legend_side = "left"))  
ggarrange(grob_FC)


#### GENE SIGNATURES #### 
# Group A : Prothrombotic Coagulation Factors 
group_A <- c("FGA", "FGB", "FGG", "F2", "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F13A1", "F13B", "VWF", "KLKB1", "KNG1", "PTAFR")
group_A <- group_A[group_A %in% rownames(rlog_counts_symbol)]

signature_list_A <- list(group_A)
names(signature_list_A) <- "GeneSetA"

# Group B : Endogenous anticoagulant / antithrombotic factors 
group_B <- c("PROC", "PROS", "SERPINC1", "TFP1", "SERPIND1", "THBD", "PROCR", "TFPI2", "HSPG", "PTGIS")
group_B <- group_B[group_B %in% rownames(rlog_counts_symbol)]

signature_list_B <- list(group_B)
names(signature_list_B) <- "GeneSetB"

#Group C : Profibrinolytic factors 
group_C <- c("PLG", "PLAT", "ANXA2", "PLAU", "PLAUR")
group_C <- group_C[group_C %in% rownames(rlog_counts_symbol)]

signature_list_C <- list(group_C)
names(signature_list_C) <- "GeneSetC"

#Group D : Antifibrinolytic factors 
group_D <- c("CPB2", "SERPINE1", "SERPINB2", "SERPINA5", "SERPINF2", "A2M", "SERPINE2", "SERPINI1", "SERPINA1", "SERPING1")
group_D <- group_D[group_D %in% rownames(rlog_counts_symbol)]

signature_list_D <- list(group_D)
names(signature_list_D) <- "GeneSetD"

#Group E : Platelet activation factors 
group_E <- c("SELP", "PF4", "TBXA2R", "CLEC1B", "CD40LG", "PECAM1", "ITGA2B", "ITGB3", "F2R", "F2RL3")
group_E <- group_E[group_E %in% rownames(rlog_counts_symbol)]

signature_list_E <- list(group_E)
names(signature_list_E) <- "GeneSetE"

# Group F : Balancing endothelial factors  
group_F <- c("FN1", "VTN")
group_F <- group_F[group_F %in% rownames(rlog_counts_symbol)]

signature_list_F <- list(group_F)
names(signature_list_F) <- "GeneSetF"


# Compute signature
gene_list <- list(group_A, group_B, group_C, group_D, group_E, group_F)
names(gene_list) <- c("group_A", "group_B", "group_C", "group_D", "group_E", "group_F")

params <- gsvaParam(as.matrix(rlog_counts_symbol),gene_list)  
gsva_res <- gsva(params)


### PREPARE COLSIDE COLORS
Condition <- c(rep("VTE",7),rep( "Non_VTE",13))
Condition_colors <- c(VTE="blue", Non_VTE="orange")

metadata[,3:14] <- lapply(metadata[, 3:14], as.factor)
str(metadata)

# Rename the factor levels (for example: converting 1 to "Low", 2 to "Medium", and 3 to "High")
levels(metadata[, 5]) <- c("zero", "one") 
levels(metadata[, 6]) <- c("zero", "one")  
levels(metadata[, 7]) <- c("four", "six") 
levels(metadata[, 8]) <- c("zero", "one", "two")  
levels(metadata[, 9]) <- c("one", "three")  
levels(metadata[, 10]) <- c("zero", "one")  
levels(metadata[, 11]) <-c("zero", "one") 
levels(metadata[, 12]) <- c("zero", "one")   

str(metadata)

batch_colors <- c("batch_1"="lightyellow", "batch_2"="lightgreen", "batch_3"="lightblue", "batch_4"="orchid", "batch_5"="lightpink")
HRD_colors <- c("zero" = "darkgoldenrod1", "one" = "cornflowerblue")
FIGO_colors <- c("four" = "coral", "six" = "darkgreen")
preop_anticoag_colors <- c("zero" = "aquamarine", "one" = "darkblue")
preop_statin_colors <- c("zero" = "bisque1", "one" = "cyan4")
intra_bleeding_colors <- c("zero" = "blueviolet", "one" = "chocolate")

rownames(gsva_res) <- c("Prothrombotic Coagulation Factors", "Antithrombotic Factors",
                        "Profibrinolytic Factors", "Antifibrinolytic Factors", 
                        "Platelet Activation Factors", "Balancing Endothelial Factors"
)

column_ha <- HeatmapAnnotation(
  
  Batch = as.factor(metadata$batch),
  preop_statin = metadata$preop_statin,
  preop_anticoag = metadata$preop_anticoag,
  HRD = metadata$HRD,
  intra_bleeding = metadata$intra_bleeding,
  Condition = Condition,
  col = list(
    Batch = batch_colors,
    preop_anticoag = preop_anticoag_colors,
    preop_statin = preop_statin_colors,
    HRD = HRD_colors,
    intra_bleeding = intra_bleeding_colors,
    Condition = Group_colors
  ),
  show_annotation_name = T,
  simple_anno_size = unit(0.3, "cm"),
  annotation_label = c("Batch", "Preop. Statin", "Preop. Anticoagulant", "HRD", "Intra Bleeding", "Condition"),
  annotation_name_gp =gpar(fontsize = 7, font=2) ,
  annotation_names_col=NULL,
  annotation_legend_param = list(
    preop_anticoag = list(
      title = "Preop. Anticoagulant",
      title_gp = gpar(fontsize = 8, font=2),
      labels_gp = gpar(fontsize = 8),
      labels = c("Yes", "No")
    ),
    
    preop_statin= list(
      title = "Preop. Statin",
      title_gp = gpar(fontsize = 8, font=2),
      labels_gp = gpar(fontsize = 8),
      labels = c("Yes", "No")
    ),
    Batch = list(
      title = "Batch",
      title_gp = gpar(fontsize = 8, font=2),
      labels_gp = gpar(fontsize = 8),
      labels = c("1", "2", "3", "4","5")
    ),
    HRD = list(
      title = "HRD",
      title_gp = gpar(fontsize = 8, font=2),
      labels_gp = gpar(fontsize = 8),
      labels = c("Yes", "No")
    ),
    intra_bleeding = list(
      title = "Intra Bleeding",
      title_gp = gpar(fontsize = 8, font=2),
      labels_gp = gpar(fontsize = 8),
      labels = c("Yes", "No")
    ),
      
      Condition = list(
        title = "Condition",
        title_gp = gpar(fontsize = 8, font=2),
        labels_gp = gpar(fontsize = 8),
        labels = c("Non VTE", "VTE")
      )
    )
  )


ht_gsva <- Heatmap(gsva_res, name = "Z-score", cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize =8, font=1 ), 
        column_title=" ",
        top_annotation = column_ha,
        column_split  =c( rep("Non VTE",7),rep("VTE",13)),
        show_column_names = F,
        column_names_gp = gpar(fontsize =0),  # To make column names invisible
        heatmap_legend_param = list(
          legend_height = unit(3 ,"cm"),
          grid_width=unit(0.3, "cm"),
          title_position = "lefttop-rot",
          labels_gp = gpar( fontsize = 8),title_gp = gpar(fontsize = 8, font=2)))
        
grob_gsva = grid.grabExpr(draw( ht_gsva, annotation_legend_side = "bottom", heatmap_legend_side = "left" ))  
ggarrange(grob_gsva)

## Compute if there are differences in the scores among groups
p_results <- apply(gsva_res, 1, function(row)) {
  wilcox.test(row ~metadata$Group)$p.value
}


p_results


#### BUILD PANEL FIGURE #### 
# Figure is built by sections.
TOP_L <- ggarrange(grob_FC, grob_gsva, labels = c("A", "B"), widths = c(0.3,1))
TOP_L

TOP_R <- ggarrange(pca_plot_DEGS, grob_untargeted, labels = c("C", "D"), heights  = c(0.6,1), ncol = 1)
TOP_R

TOP <- ggarrange(TOP_L,TOP_R, widths = c(1,0.5))
TOP

DOWN <- ggarrange(volcano, gsea_plot, labels = c("E", "F"), widths = c(0.6,1) )
DOWN

ALL <- ggarrange(TOP, DOWN, nrow = 2, heights = c(1,1))
ALL

## FINAL FIGURE:
pdf("Panel.pdf", width =12 , height = 10)
ALL
dev.off()
