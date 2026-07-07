###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 04_DEA.R
# Author: Blanca Martin Urdiales
# Description: Differential Expression Analysis (DEA) for Female, Male and Sex interaction
#              Generates volcano plots and a final merged results table
###############################################
getwd()
setwd("C:/Users/blanc/OneDrive/Escritorio/miRNAuEVsAUD/miRNAuEVsAUD-Zenodo")
# --- Load required packages ---
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Load filtered DESeq2 object ---
dds1_file <- "Results/dds1_filtered.rds"
dds1 <- readRDS(dds1_file)
coldata <- readRDS("metadata/coldata.rds")
# --- Run DESeq2 ---
dds1 <- DESeq(dds1)

# ==================== FEMALE (IF) ==========================
resFemale <- results(dds1, contrast = c("Condition", "AUD_Female", "Control_Female"),
                     alpha = 0.05, pAdjustMethod = "BH")
df_resFemale <- as.data.frame(resFemale)
df_resFemale$mirna <- rownames(df_resFemale)
df_resFemale$padj[is.na(df_resFemale$padj)] <- 1

# Classify significance
res_F_1 <- df_resFemale %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Positive LFC",
    padj < 0.05 & log2FoldChange < 0 ~ "Negative LFC",
    TRUE ~ "Not Significant"
  ))

# Volcano plot
volcano_IF <- ggplot(res_F_1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = Significance), color = "#575757ff", shape = 21, size = 2) +
  scale_fill_manual(values = c("Positive LFC"="#542B51", "Negative LFC"="#D3AAD1", "Not Significant"="#808080")) +
  geom_text_repel(data = subset(res_F_1, padj < 0.05), aes(label = mirna),
                  size = 3, max.overlaps = 10) +
  xlab("LFC") + ylab("-log10(p.adjusted)") + theme_bw() + ggtitle("IF")
volcano_IF


# ==================== MALE (IM) ==========================
resMale <- results(dds1, contrast = c("Condition", "AUD_Male", "Control_Male"),
                   alpha = 0.05, pAdjustMethod = "BH")
df_resMale <- as.data.frame(resMale)
df_resMale$mirna <- rownames(df_resMale)
df_resMale$padj[is.na(df_resMale$padj)] <- 1

res_M_2 <- df_resMale %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Positive LFC",
    padj < 0.05 & log2FoldChange < 0 ~ "Negative LFC",
    TRUE ~ "Not Significant"
  ))

volcano_IM <- ggplot(res_M_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = Significance), color = "#575757ff", shape = 21, size = 2) +
  scale_fill_manual(values = c("Positive LFC"="#1C2E4C", "Negative LFC"="#9BADCC", "Not Significant"="#808080")) +
  geom_text_repel(data = subset(res_M_2, padj < 0.05), aes(label = mirna),
                  size = 3, max.overlaps = 10) +
  xlab("LFC") + ylab("-log10(p.adjusted)") + theme_bw() + ggtitle("IM")
volcano_IM


# ==================== SEX INTERACTION (IS) ==========================
resSex <- results(dds1,
                  contrast = list(c("ConditionAUD_Female","ConditionControl_Male"),
                                  c("ConditionAUD_Male","ConditionControl_Female")),
                  alpha = 0.05, pAdjustMethod = "BH")
df_resSex <- as.data.frame(resSex)
df_resSex$mirna <- rownames(df_resSex)
df_resSex$padj[is.na(df_resSex$padj)] <- 1

res_S_3 <- df_resSex %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Positive LFC",
    padj < 0.05 & log2FoldChange < 0 ~ "Negative LFC",
    TRUE ~ "Not Significant"
  ))

volcano_IS <- ggplot(res_S_3, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = Significance), color = "#575757ff", shape = 21, size = 2) +
  scale_fill_manual(values = c("Positive LFC"="#315B2E", "Negative LFC"="#B1DAAE", "Not Significant"="#A19D9F")) +
  geom_text_repel(data = subset(res_S_3, padj < 0.05), aes(label = mirna),
                  size = 3, max.overlaps = 10) +
  xlab("LFC") + ylab("-log10(p.adjusted)") + theme_bw() + ggtitle("IS")
volcano_IS


# ==================== Generate merged results table ==========================
df_IF <- df_resFemale %>%
  select(mirna, log2FoldChange_IF = log2FoldChange, padj_IF = padj)
df_IM <- df_resMale %>%
  select(mirna, log2FoldChange_IM = log2FoldChange, padj_IM = padj)
df_IS <- df_resSex %>%
  select(mirna, log2FoldChange_IS = log2FoldChange, padj_IS = padj)


final_df_orina <- df_IF %>%
  full_join(df_IM, by = "mirna") %>%
  full_join(df_IS, by = "mirna") 


# Save as RDS
saveRDS(final_df_orina, file = "Results/final_df_orina.rds")
saveRDS(df_IF, file = "Results/df_IF.rds")
saveRDS(df_IM, file = "Results/df_IM.rds")
saveRDS(df_IS, file = "Results/df_IS.rds")

