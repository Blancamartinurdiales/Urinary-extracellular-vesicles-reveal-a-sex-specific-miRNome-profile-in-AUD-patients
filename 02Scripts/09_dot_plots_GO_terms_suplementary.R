###############################################
# Project: miRNome Profiling of Urinary EVs in AUD Patients
# Script: 09_dot_plots_GO_terms_suplementary
# Author: Blanca Martin Urdiales
###############################################
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Load data
archivos <- list.files(path = "Results", pattern = "org_db.*\\.RData$", full.names = TRUE)

for (archivo in archivos) {
  load(archivo)
}

listas <- list(
  gsea_org_db_Females = gsea_org_db_Females,
  gsea_org_db_Males = gsea_org_db_Males,
  gsea_org_db_Sex = gsea_org_db_Sex,
  gsea_org_db_Case = gsea_org_db_Case
)


standardize_columns <- function(df) {
  if ("path_name" %in% colnames(df)) {
    colnames(df)[colnames(df) == "path_name"] <- "description"
  }
  if ("summary" %in% colnames(df)) {
    colnames(df)[colnames(df) == "summary"] <- "description"
  }
  return(df)
}

for (i in seq_along(listas)) {
  listas[[i]] <- lapply(listas[[i]], standardize_columns)
}

IF_BP <- gsea_org_db_Females$gsea_bp
IF_BP <- IF_BP %>%
  arrange(padj)

IM_BP <- gsea_org_db_Males$gsea_bp
IM_BP <- IM_BP %>%
  arrange(padj)

IF_CC <- gsea_org_db_Females$gsea_cc
IF_CC <- IF_CC %>%
  arrange(padj)

IM_CC <- gsea_org_db_Males$gsea_cc
IM_CC <- IM_CC %>%
  arrange(padj)

IF_MF <- gsea_org_db_Females$gsea_mf
IF_MF <- IF_MF %>%
  arrange(padj)

IM_MF <- gsea_org_db_Males$gsea_mf
IM_MF <- IM_MF %>%
  arrange(padj)

IF_BP_top20 <- IF_BP[1:20,]
IF_BP_top20$contrast <- "IF"
IM_BP_top20 <- IM_BP[1:20,]
IM_BP_top20$contrast <- "IM"
IF_CC_top20 <- IF_CC[1:20,]
IF_CC_top20$contrast <- "IF"
IM_CC_top20 <- IM_CC[1:20,]
IM_CC_top20$contrast <- "IM"
IF_MF_top20 <- IF_MF[1:20,]
IF_MF_top20$contrast <- "IF"
IM_MF_top20 <- IM_MF[1:20,]
IM_MF_top20$contrast <- "IM"

IF_BP_top20_terms <- IF_BP_top20$description
IM_BP_top20_terms <- IM_BP_top20$description
IF_CC_top20_terms <- IF_CC_top20$description
IM_CC_top20_terms <- IM_CC_top20$description
IF_MF_top20_terms <- IF_MF_top20$description
IM_MF_top20_terms <- IM_MF_top20$description

df_BP <- rbind(IF_BP_top20, IM_BP_top20)
df_BP <- df_BP %>%
  mutate(sig = if_else(padj < 0.05, "sig", "no_sig"))
df_BP <- df_BP[-40,]
df_BP <- df_BP[-37,]

df_CC <- rbind(IF_CC_top20,IM_CC_top20)
df_CC <- df_CC %>%
  mutate(sig = if_else(padj < 0.05, "sig", "no_sig"))

df_MF <- rbind(IF_MF_top20,IM_MF_top20)
df_MF <- df_MF %>%
  mutate(sig = if_else(padj < 0.05, "sig", "no_sig"))


ggplot(df_BP, aes(x = contrast, y = description)) +
  geom_point(aes(fill = lor, shape = sig), color = "grey", size = 3, stroke = 0.8) +
  scale_shape_manual(values = c("sig" = 21, "no_sig" = NA)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "LOR"
  ) +
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "BP pathway",
    shape = "Significant"
  ) +
  ggtitle("Top 20 BP pathways - Comparación IF vs IM") +
  theme(axis.text.y = element_text(size = 8))

ggplot(df_CC, aes(x = contrast, y = description)) +
  geom_point(aes(fill = lor, shape = sig), color = "grey", size = 3, stroke = 0.8) +
  scale_shape_manual(values = c("sig" = 21, "no_sig" = 24)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "LOR"
  ) +
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "CC pathway",
    shape = "Significant"
  ) +
  ggtitle("Top 20 CC pathways - Comparación IF vs IM") +
  theme(axis.text.y = element_text(size = 8))


ggplot(df_MF, aes(x = contrast, y = description)) +
  geom_point(aes(fill = lor, shape = sig), color = "grey", size = 3, stroke = 0.8) +
  scale_shape_manual(values = c("sig" = 21, "no_sig" = NA)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "LOR"
  ) +
  theme_minimal() +
  labs(
    x = "Comparison",
    y = "MF pathway",
    shape = "Significant"
  ) +
  ggtitle("Top 20 MF pathways - Comparación IF vs IM") +
  theme(axis.text.y = element_text(size = 8))

