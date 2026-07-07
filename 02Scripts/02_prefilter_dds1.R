###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 02_prefilter_dds1.R
# Author: Blanca Martin Urdiales
# Description: Load raw DESeq2 object, perform CPM-based pre-filtering,
#              and generate a QC plot of retained vs removed genes
###############################################

# --- Load required packages ---
library(DESeq2)
library(edgeR)   # Used for CPM calculation -> cpm()
library(dplyr)
library(ggplot2)

# --- Load DESeq2 object ---
dds1_file <- "Results/dds1_raw_counts.rds"
dds1 <- readRDS(dds1_file)

# --- Quick check ---
dim(dds1)  # Number of genes x samples

# --- Pre-filtering using CPM ---
# Keep only genes with at least 1 CPM in 4 or more samples
# This removes genes with very low expression (noise)
# Note: CPM is computed from raw counts using edgeR

keep <- rowSums(edgeR::cpm(dds1) >= 1) >= 4
table(keep)  # Count of retained vs removed genes

# --- Save number of genes before and after filtering for QC ---
qc_df <- data.frame(
  Status = c("Removed", "Retained"),
  Count  = c(sum(!keep), sum(keep))
)

# --- Plot QC ---
qc_plot <- ggplot(qc_df, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Pre-filtering QC",
       subtitle = "Genes retained vs removed based on CPM >= 1 in >=4 samples",
       y = "Number of genes",
       x = "") +
  scale_fill_manual(values = c("Removed" = "#D55E00", "Retained" = "#009E73")) +
  theme(legend.position = "none")

# Show plot
print(qc_plot)

# --- Apply filtering ---
dds1 <- dds1[keep, ]

# --- Re-estimate size factors after filtering ---
dds1 <- estimateSizeFactors(dds1, type = "ratio")

# --- Save filtered DESeq2 object ---
saveRDS(dds1, file = "Results/dds1_filtered.rds")
