###############################################
# Project: miRNome Profiling of Urinary EVs in AUD Patients
# Script: 07_funtional analysis_kegg_reactome_description
# Author: Blanca Martin Urdiales
###############################################

# Required packages
library(clusterProfiler)
library(UpSetR)
library(ComplexUpset)
library(enrichplot)
library(RColorBrewer)
library(MetBrewer)
library(ggpubr)
library(ggplot2)
library(dplyr)

# Load data
load("Results/kegg_db.RData")
load("Results/reactome_pathways.RData")


load("Results/gsea_org_db_Males.RData")
gsea_bp_Males <- gsea_org_db_Males$gsea_bp
gsea_mf_Males <- gsea_org_db_Males$gsea_mf
gsea_cc_Males <- gsea_org_db_Males$gsea_cc
gsea_reactome_Males <- gsea_org_db_Males$gsea_reactome
gsea_kegg_Males <- gsea_org_db_Males$gsea_kegg

load("Results/gsea_org_db_Females.RData")
gsea_bp_Females <- gsea_org_db_Females$gsea_bp
gsea_mf_Females <- gsea_org_db_Females$gsea_mf
gsea_cc_Females <- gsea_org_db_Females$gsea_cc
gsea_reactome_Females <- gsea_org_db_Females$gsea_reactome
gsea_kegg_Females <- gsea_org_db_Females$gsea_kegg

load("Results/gsea_org_db_Sex.RData")
gsea_bp_Sex <- gsea_org_db_Sex$gsea_bp
gsea_mf_Sex <- gsea_org_db_Sex$gsea_mf
gsea_cc_Sex <- gsea_org_db_Sex$gsea_cc
gsea_reactome_Sex <- gsea_org_db_Sex$gsea_reactome
gsea_kegg_Sex <- gsea_org_db_Sex$gsea_kegg

load("Results/gsea_org_db_Case.RData")
gsea_bp_Case <- gsea_org_db_Case$gsea_bp
gsea_mf_Case <- gsea_org_db_Case$gsea_mf
gsea_cc_Case <- gsea_org_db_Case$gsea_cc
gsea_reactome_Case <- gsea_org_db_Case$gsea_reactome
gsea_kegg_Case <- gsea_org_db_Case$gsea_kegg


# Add description column to GSEA KEGG results
# Convert rownames (hsaXXXX) into a column named "codigo"
gsea_kegg_Males$codigo <- rownames(gsea_kegg_Males)
gsea_kegg_Females$codigo <- rownames(gsea_kegg_Females)
gsea_kegg_Sex$codigo <- rownames(gsea_kegg_Sex)
gsea_kegg_Case$codigo <- rownames(gsea_kegg_Case)

# Merge with KEGG pathway annotation table
gsea_kegg_Males <- merge(gsea_kegg_Males, kegg_db,
                         by = "codigo", all.x = TRUE)
gsea_kegg_Females <- merge(gsea_kegg_Females, kegg_db,
                           by = "codigo", all.x = TRUE)
gsea_kegg_Sex <- merge(gsea_kegg_Sex, kegg_db,
                       by = "codigo", all.x = TRUE)
gsea_kegg_Case <- merge(gsea_kegg_Case, kegg_db,
                        by = "codigo", all.x = TRUE)
# Rename "summary" column to "description"
colnames(gsea_kegg_Males)[colnames(gsea_kegg_Males) == "summary"] <- "description"
colnames(gsea_kegg_Females)[colnames(gsea_kegg_Females) == "summary"] <- "description"
colnames(gsea_kegg_Sex)[colnames(gsea_kegg_Sex) == "summary"] <- "description"
colnames(gsea_kegg_Case)[colnames(gsea_kegg_Case) == "summary"] <- "description"

# Clean KEGG pathway names by removing species suffix
for (obj in c("gsea_kegg_Males", "gsea_kegg_Females", "gsea_kegg_Sex", "gsea_kegg_Case")) {
  tmp <- get(obj)
  tmp$description <- gsub(" - Homo sapiens \\(human\\)", "", tmp$description)
  assign(obj, tmp)
}


# Add description column for Reactome pathways
# Add pathway ID column
gsea_reactome_Males$codigo <- rownames(gsea_reactome_Males)
gsea_reactome_Females$codigo <- rownames(gsea_reactome_Females)
gsea_reactome_Sex$codigo <- rownames(gsea_reactome_Sex)
gsea_reactome_Case$codigo <- rownames(gsea_reactome_Case)

# Merge with Reactome pathway annotation table
gsea_reactome_Males <- merge(gsea_reactome_Males, reactome_pathways,
                             by = "codigo", all.x = TRUE)
gsea_reactome_Females <- merge(gsea_reactome_Females, reactome_pathways,
                               by = "codigo", all.x = TRUE)
gsea_reactome_Sex <- merge(gsea_reactome_Sex, reactome_pathways,
                           by = "codigo", all.x = TRUE)
gsea_reactome_Case <- merge(gsea_reactome_Case, reactome_pathways,
                            by = "codigo", all.x = TRUE)


for (obj in c("gsea_reactome_Males", "gsea_reactome_Females", "gsea_reactome_Sex", "gsea_reactome_Case")) {
  tmp <- get(obj)
  tmp$description <- gsub("^Homo sapiens: ?", "", tmp$description)
  assign(obj, tmp)
}

gsea_kegg_Males <- gsea_kegg_Males[,-7]

# Example for Males
gsea_org_db_Males$gsea_kegg      <- gsea_kegg_Males
gsea_org_db_Males$gsea_reactome  <- gsea_reactome_Males

# Females
gsea_org_db_Females$gsea_kegg      <- gsea_kegg_Females
gsea_org_db_Females$gsea_reactome  <- gsea_reactome_Females

# Sex
gsea_org_db_Sex$gsea_kegg      <- gsea_kegg_Sex
gsea_org_db_Sex$gsea_reactome  <- gsea_reactome_Sex

# Case
gsea_org_db_Case$gsea_kegg      <- gsea_kegg_Case
gsea_org_db_Case$gsea_reactome  <- gsea_reactome_Case

# Save updated GSEA objects
save(gsea_org_db_Males, file = "Results/gsea_org_db_Males.RData")
save(gsea_org_db_Females, file = "Results/gsea_org_db_Females.RData")
save(gsea_org_db_Sex, file = "Results/gsea_org_db_Sex.RData")
save(gsea_org_db_Case, file = "Results/gsea_org_db_Case.RData")
