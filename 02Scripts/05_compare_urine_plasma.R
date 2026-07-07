###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 05_compare_plasma.R
# Author: Blanca Martin Urdiales
# Description: Comparison of significant miRNAs in urine vs plasma, generating Venn diagrams
###############################################

# --- Load required packages ---
library(dplyr)
library(ggvenn)
library(ggplot2)
library(readr)

# --- Load data ---
final_df_plasma <- readRDS("Data/final_df_plasma.rds")
final_df_orina <- readRDS("Results/final_df_orina.rds")


# --- Extract miRNA names ---
mirnas_plasma <- final_df_plasma$mirna
mirnas_orina  <- final_df_orina$mirna

# --- Significant miRNAs ---
sig_IF_plasma <- final_df_plasma %>% filter(padj_IF < 0.05) %>% pull(mirna)
sig_IM_plasma <- final_df_plasma %>% filter(padj_IM < 0.05) %>% pull(mirna)
sig_IS_plasma <- final_df_plasma %>% filter(padj_IS < 0.05) %>% pull(mirna)

sig_IF_orina <- final_df_orina %>% filter(padj_IF < 0.05) %>% pull(mirna)
sig_IM_orina <- final_df_orina %>% filter(padj_IM < 0.05) %>% pull(mirna)
sig_IS_orina <- final_df_orina %>% filter(padj_IS < 0.05) %>% pull(mirna)

# --- Intersections and differences ---
comunes <- intersect(mirnas_plasma, mirnas_orina)
solo_plasma <- setdiff(mirnas_plasma, mirnas_orina)
solo_orina <- setdiff(mirnas_orina, mirnas_plasma)

sig_IF_orina_en_plasma <- intersect(mirnas_plasma, sig_IF_orina)
sig_IM_orina_en_plasma <- intersect(mirnas_plasma, sig_IM_orina)
sig_IS_orina_en_plasma <- intersect(mirnas_plasma, sig_IS_orina)

sig_IF_orina_NO_en_plasma <- setdiff(sig_IF_orina, mirnas_plasma)
sig_IM_orina_NO_en_plasma <- setdiff(sig_IM_orina, mirnas_plasma)
sig_IS_orina_NO_en_plasma <- setdiff(sig_IS_orina, mirnas_plasma)

# --- Create list of sets for Venn diagrams ---
venn_list <- list(
  IF = list(
    Plasma = mirnas_plasma,
    Orina_sig_IF = sig_IF_orina
  ),
  IM = list(
    Plasma = mirnas_plasma,
    Orina_sig_IM = sig_IM_orina
  ),
  IS = list(
    Plasma = mirnas_plasma,
    Orina_sig_IS = sig_IS_orina
  ),
  IF_sig = list(
    Plasma_sig_IF = sig_IF_plasma,
    Orina_sig_IF  = sig_IF_orina
  ),
  IM_sig = list(
    Plasma_sig_IM = sig_IM_plasma,
    Orina_sig_IM  = sig_IM_orina
  ),
  IS_sig = list(
    Plasma_sig_IS = sig_IS_plasma,
    Orina_sig_IS  = sig_IS_orina
  )
)

# --- Function to generate ggvenn plots ---
plot_venn_sig <- function(set_list, title){
  ggvenn(
    set_list,
    fill_color = c("#E69F00", "#56B4E9"),
    stroke_color = NA,
    stroke_size = 0,
    set_name_size = 5,
    text_size = 10,
    show_percentage = FALSE
  ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
}

# --- Generate Venn diagrams ---
venn_IF <- plot_venn_sig(venn_list$IF, "Significant IF miRNAs: Urine vs Plasma")
venn_IM <- plot_venn_sig(venn_list$IM, "Significant IM miRNAs: Urine vs Plasma")
venn_IS <- plot_venn_sig(venn_list$IS, "Significant IS miRNAs: Urine vs Plasma")

pdf("venn_IF.pdf", width = 14, height = 10) 
venn_IF
dev.off()
pdf("venn_IM.pdf", width = 14, height = 10) 
venn_IM
dev.off()

venn_IF_sig <- plot_venn_sig(venn_list$IF_sig, "Significant IF miRNAs: Urine vs Plasma Significant")
venn_IM_sig <- plot_venn_sig(venn_list$IM_sig, "Significant IM miRNAs: Urine vs Plasma Significant")
venn_IS_sig <- plot_venn_sig(venn_list$IS_sig, "Significant IS miRNAs: Urine vs Plasma Significant")

# --- Print summary of intersections ---
cat("Significant IF miRNAs in Urine present in Plasma:\n")
print(sig_IF_orina_en_plasma)
cat("Significant IF miRNAs in Urine NOT in Plasma:\n")
print(sig_IF_orina_NO_en_plasma)

cat("Significant IM miRNAs in Urine present in Plasma:\n")
print(sig_IM_orina_en_plasma)
cat("Significant IM miRNAs in Urine NOT in Plasma:\n")
print(sig_IM_orina_NO_en_plasma)

cat("Significant IS miRNAs in Urine present in Plasma:\n")
print(sig_IS_orina_en_plasma)
cat("Significant IS miRNAs in Urine NOT in Plasma:\n")
print(sig_IS_orina_NO_en_plasma)

