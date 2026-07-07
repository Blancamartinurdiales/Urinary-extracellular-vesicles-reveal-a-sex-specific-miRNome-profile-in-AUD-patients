###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 10_simplify_enrichment.R
# Author: Blanca Martin Urdiales
# Description: Visualization of GO terms using simplifyEnrichment .
###############################################

# --- Load required packages ---
library(dplyr)
library(GO.db)
library(grid)

# Load data
archivos <- list.files(path = "Results", pattern = "org_db.*\\.RData$", full.names = TRUE)
for (archivo in archivos) {
  load(archivo)
}

#===============================================================================================================
# simplify enrichment ####
#===============================================================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("simplifyEnrichment")

# Load the package
library(simplifyEnrichment)

#===================================================================================================================

## all IF ####

go_all_females <- list(bp_down = rownames(gsea_org_db_Females$gsea_bp[gsea_org_db_Females$gsea_bp$padj <= 0.05 &
                                                                        gsea_org_db_Females$gsea_bp$lor <= 0,]),
                       bp_up = rownames(gsea_org_db_Females$gsea_bp[gsea_org_db_Females$gsea_bp$padj <= 0.05 &
                                                                      gsea_org_db_Females$gsea_bp$lor >= 0,]),
                       bp = rownames(gsea_org_db_Females$gsea_bp[gsea_org_db_Females$gsea_bp$padj <= 0.05,]),
                       
                       cc_down = rownames(gsea_org_db_Females$gsea_cc[gsea_org_db_Females$gsea_cc$padj <= 0.05 &
                                                                        gsea_org_db_Females$gsea_cc$lor <= 0,]),
                       cc_up = rownames(gsea_org_db_Females$gsea_cc[gsea_org_db_Females$gsea_cc$padj <= 0.05 &
                                                                      gsea_org_db_Females$gsea_cc$lor >= 0,]),
                       cc = rownames(gsea_org_db_Females$gsea_cc[gsea_org_db_Females$gsea_cc$padj <= 0.05,])
)

set.seed(123)

IF_bp_up = GO_similarity(go_all_females$bp_up, ont = "BP")
pdf("IF_bp_up.pdf", width = 10, height = 10)
simplifyGO(IF_bp_up, method = "binary_cut")
dev.off()


IF_bp_down = GO_similarity(go_all_females$bp_down, ont = "BP")
pdf("IF_bp_down.pdf", width = 10, height = 10)
simplifyGO(IF_bp_down, method = "binary_cut")
dev.off()


IF_cc_up = GO_similarity(go_all_females$cc_up, ont = "CC")
pdf("IF_cc_up.pdf", width = 10, height = 10)
simplifyGO(IF_cc_up, method = "binary_cut")
dev.off()

IF_cc_down = GO_similarity(go_all_females$cc_down, ont = "CC")
pdf("IF_cc_down.pdf", width = 10, height = 10)
simplifyGO(IF_cc_down, method = "binary_cut")
dev.off()


## all IM ####
go_all_males <- list(bp_down = rownames(gsea_org_db_Males$gsea_bp[gsea_org_db_Males$gsea_bp$padj <= 0.05 &
                                                                    gsea_org_db_Males$gsea_bp$lor <= 0,]),
                     bp_up = rownames(gsea_org_db_Males$gsea_bp[gsea_org_db_Males$gsea_bp$padj <= 0.05 &
                                                                  gsea_org_db_Males$gsea_bp$lor >= 0,]),
                     bp = rownames(gsea_org_db_Males$gsea_bp[gsea_org_db_Males$gsea_bp$padj <= 0.05,]),
                     
                     cc_down = rownames(gsea_org_db_Males$gsea_cc[gsea_org_db_Males$gsea_cc$padj <= 0.05 &
                                                                    gsea_org_db_Males$gsea_cc$lor <= 0,]),
                     cc_up = rownames(gsea_org_db_Males$gsea_cc[gsea_org_db_Males$gsea_cc$padj <= 0.05 &
                                                                  gsea_org_db_Males$gsea_cc$lor >= 0,]),
                     cc = rownames(gsea_org_db_Males$gsea_cc[gsea_org_db_Males$gsea_cc$padj <= 0.05,])
)


IM_bp_up = GO_similarity(go_all_males$bp_up, ont = "BP")
pdf("IM_bp_up.pdf", width = 10, height = 10)
simplifyGO(IM_bp_up, method = "binary_cut")
dev.off()

IM_bp_down = GO_similarity(go_all_males$bp_down, ont = "BP")
pdf("IM_bp_down.pdf", width = 10, height = 10)
simplifyGO(IM_bp_down, method = "binary_cut")
dev.off()


IM_cc_up = GO_similarity(go_all_males$cc_up, ont = "CC")
pdf("IM_cc_up.pdf", width = 10, height = 10)
simplifyGO(IM_cc_up, method = "binary_cut")
dev.off()


IM_cc_down = GO_similarity(go_all_males$cc_down, ont = "CC")
pdf("IM_cc_down.pdf", width = 10, height = 10)
simplifyGO(IM_cc_down, method = "binary_cut")
dev.off()
