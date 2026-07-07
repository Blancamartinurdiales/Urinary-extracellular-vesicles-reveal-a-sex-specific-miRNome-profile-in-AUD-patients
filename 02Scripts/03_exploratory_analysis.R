###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 03_exploratory_analysis.R
# Author: Blanca Martin Urdiales
# Description: Exploratory data analysis including 
#              normalization, PCA, and heatmaps
###############################################

# --- Load required packages ---
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(plotly)
library(pheatmap)
library(dplyr)
# Load custom PCA function
source("01Functions/function_pcaGenes_2.r")

# --- Load filtered DESeq2 object ---
dds1_file <- "Results/dds1_filtered.rds"
dds1 <- readRDS(dds1_file)
coldata <- readRDS("metadata/coldata.rds")

# --- Quick check of counts ---
head(counts(dds1, normalized=FALSE))

# --- Estimate size factors ---
dds1 <- estimateSizeFactors(dds1, type = "ratio")
normalized_counts <- counts(dds1, normalized = TRUE)
assay(dds1, "counts.norm") <- normalized_counts

# --- Variance Stabilizing Transformation (VST) ---
assay(dds1, "counts.norm.VST") <- as.data.frame(
  assay(varianceStabilizingTransformation(dds1, blind = TRUE)),
  check.names = FALSE
)

# --- Boxplots of raw, median ratio normalized, and VST counts ---
par(mfrow = c(1,3))
boxplot(log10(as.matrix(assay(dds1, "counts")+1)),
        ylab = expression('Log'[10]~'Read counts'),
        las = 2,
        main = "Raw counts filtered")
boxplot(log10(as.matrix(assay(dds1, "counts.norm")+1)),
        ylab = expression('Log'[10]~'Read counts'),
        las = 2,
        main = "Median ratio")
boxplot(log10(as.matrix(assay(dds1, "counts.norm.VST")+1)),
        ylab = expression('Log'[10]~'Read counts'),
        las = 2,
        main = "VST")

# --- PCA Analysis (using VST data) ---

palette4 <- c("#dd5129", "#0f7ba2", "#43b284", "#fab255")  # MetBrewer palette example

# Run PCA (custom function pcaGenes assumed loaded from functions)
mi.pca <- pcaGenes(assay(dds1, "counts.norm.VST"))
mi.pca.df <- as.data.frame(mi.pca$scores)

mi.pca.df$grupo <- coldata$Condition
mi.pca.df$var.exp <- round(mi.pca$var.exp * 100, 2)
rownames(mi.pca.df) <- coldata$Sample

# --- Interactive 3D PCA plot ---
fig_pca <- plot_ly(mi.pca.df, x = ~V1, y = ~V2, z = ~V3,
                   text = rownames(mi.pca.df),
                   color = ~grupo, colors = palette4) %>%
  add_markers(marker = list(size = 4)) %>%
  layout(title = "<b>PCA</b>",
         legend = list(title = list(text = "<b> Groups </b>")),
         scene = list(
           xaxis = list(title = paste0("PC1: ", mi.pca.df$var.exp[1], "% variance explained")),
           yaxis = list(title = paste0("PC2: ", mi.pca.df$var.exp[2], "% variance explained")),
           zaxis = list(title = paste0("PC3: ", mi.pca.df$var.exp[3], "% variance explained"))
         ))

# --- 2D PCA plots for publication ---
PC1_PC2 <- ggplot(mi.pca.df, aes(x = V1, y = V2, colour = grupo, shape = grupo)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8, size = 2.5) +
  stat_ellipse(aes(group = grupo, colour = grupo), level = 0.95, linetype = 2, linewidth = 1) +
  scale_color_manual(values = palette4) +
  scale_shape_manual(values = c(15,16,17,18)) +
  theme_minimal() +
  xlab(paste("PC1: ", round(mi.pca.df$var.exp[1]), "% explained variance", sep="")) +
  ylab(paste("PC2: ", round(mi.pca.df$var.exp[2]), "% explained variance", sep=""))

PC3_PC2 <- ggplot(mi.pca.df, aes(x = V3, y = V2, colour = grupo, shape = grupo)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8, size = 2.5) +
  stat_ellipse(aes(group = grupo, colour = grupo), level = 0.95, linetype = 2, linewidth = 1) +
  scale_color_manual(values = palette4) +
  scale_shape_manual(values = c(15,16,17,18)) +
  theme_minimal() +
  xlab(paste("PC3: ", round(mi.pca.df$var.exp[3]), "% explained variance", sep="")) +
  ylab(paste("PC2: ", round(mi.pca.df$var.exp[2]), "% explained variance", sep=""))

# Arrange both PCA plots side by side
ggpubr::ggarrange(PC1_PC2, PC3_PC2, ncol=2, nrow=1, common.legend = TRUE)

# --- Correlation heatmap ---
cor_matrix <- cor(assay(dds1, "counts.norm.VST"), method = "spearman")
heatmap <- pheatmap(cor_matrix,
                    border_color = NA,
                    annotation_col = coldata["Condition", drop=F],
                    annotation_row = coldata["Condition", drop=F],
                    annotation_legend = TRUE)



