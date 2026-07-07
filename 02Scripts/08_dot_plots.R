###############################################
# Project: miRNome Profiling of Urinary EVs in AUD Patients
# Script: 08_dot_plots
# Author: Blanca Martin Urdiales
###############################################
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Load data
# Search for files within the Results folder
archivos <- list.files(path = "Results", pattern = "org_db.*\\.RData$", full.names = TRUE)

# Load each .RData file
for (archivo in archivos) {
  load(archivo)
}

# Create a list with the loaded objects
listas <- list(
  gsea_org_db_Females = gsea_org_db_Females,
  gsea_org_db_Males = gsea_org_db_Males,
  gsea_org_db_Sex = gsea_org_db_Sex,
  gsea_org_db_Case = gsea_org_db_Case
)


# Function to standardize column names
standardize_columns <- function(df) {
  if ("path_name" %in% colnames(df)) {
    colnames(df)[colnames(df) == "path_name"] <- "description"
  }
  if ("summary" %in% colnames(df)) {
    colnames(df)[colnames(df) == "summary"] <- "description"
  }
  return(df)
}

# Apply function to all lists and categories
for (i in seq_along(listas)) {
  listas[[i]] <- lapply(listas[[i]], standardize_columns)
}

IF_kegg <- gsea_org_db_Females$gsea_kegg
IF_kegg <- IF_kegg %>%
  arrange(padj)

IM_kegg <- gsea_org_db_Males$gsea_kegg
IM_kegg <- IM_kegg %>%
  arrange(padj)

IF_reactome <- gsea_org_db_Females$gsea_reactome
IF_reactome <- IF_reactome %>%
  arrange(padj)

IM_reactome <- gsea_org_db_Males$gsea_reactome
IM_reactome <- IM_reactome %>%
  arrange(padj)

IF_kegg_top20 <- IF_kegg[1:20,]
IF_kegg_top20$contrast <- "IF"
IM_kegg_top20 <- IM_kegg[1:20,]
IM_kegg_top20$contrast <- "IM"
IF_reactome_top20 <- IF_reactome[1:20,]
IF_reactome_top20$contrast <- "IF"
IM_reactome_top20 <- IM_reactome[1:20,]
IM_reactome_top20$contrast <- "IM"

IF_kegg_top20_terms <- IF_kegg_top20$description
IM_kegg_top20_terms <- IM_kegg_top20$description
IF_reactome_top20_terms <- IF_reactome_top20$description
IM_reactome_top20_terms <- IM_reactome_top20$description

df_kegg <- rbind(IF_kegg_top20, IM_kegg_top20)
df_kegg <- df_kegg %>%
  mutate(sig = if_else(padj < 0.05, "sig", "no_sig"))
df_reactome <- rbind(IF_reactome_top20,IM_reactome_top20)
df_reactome <- df_reactome %>%
  mutate(sig = if_else(padj < 0.05, "sig", "no_sig"))


ggplot(df_kegg, aes(x = contrast, y = description)) +
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
    y = "KEGG pathway",
    shape = "Significant"
  ) +
  ggtitle("Top 20 KEGG pathways - Comparación IF vs IM") +
  theme(axis.text.y = element_text(size = 8))

ggplot(df_reactome, aes(x = contrast, y = description)) +
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
    y = "Reactome pathway",
    shape = "Significant"
  ) +
  ggtitle("Top 20 Reactome pathways - Comparación IF vs IM") +
  theme(axis.text.y = element_text(size = 8))


kegg_top_terms <- df_kegg$description
intersect(IF_kegg_top20_terms,IM_kegg_top20_terms) # 8 
length(kegg_top_terms) # 40
kegg_top_terms <- unique(kegg_top_terms) 
length(kegg_top_terms) # 32 (40-8)

reactome_top_terms <- df_reactome$description 
intersect(IF_reactome_top20_terms,IM_reactome_top20_terms) # 1
length(reactome_top_terms) # 40
reactome_top_terms <- unique(reactome_top_terms)
length(reactome_top_terms) # 39 (40-1) 

kegg_top_terms
reactome_top_terms



# Categories manually curated based on biological interpretation
# -----------------------------
# KEGG: 
# -----------------------------
df_kegg <- df_kegg %>%
  mutate(category = case_when(
    description %in% c("Pathways in cancer", "Colorectal cancer", "Prostate cancer", 
                       "Renal cell carcinoma", "Endometrial cancer", "Chronic myeloid leukemia",
                       "Pancreatic cancer", "Cell cycle", "Oocyte meiosis", "Spliceosome") ~ "Cancer / Cell proliferation",
    
    description %in% c("Apoptosis") ~ "Apoptosis / Cell death",
    
    description %in% c("Neurotrophin signaling pathway", "Axon guidance",
                       "Long-term potentiation", "Long-term depression",
                       "Olfactory transduction") ~ "Neuronal / Synaptic function",
    
    description %in% c("ErbB signaling pathway", "Insulin signaling pathway",
                       "p53 signaling pathway", "Phosphatidylinositol signaling system",
                       "Gastric acid secretion") ~ "Cell signaling / Growth regulation",
    
    description %in% c("Systemic lupus erythematosus", "Viral myocarditis") ~ "Immune response / Inflammation",
    
    description %in% c("Endocytosis", "SNARE interactions in vesicular transport",
                       "Nucleocytoplasmic transport", "Adherens junction", "Focal adhesion",
                       "Gap junction") ~ "Vesicle / Intracellular transport",
    
    description %in% c("Protein processing in endoplasmic reticulum",
                       "Ubiquitin mediated proteolysis",
                       "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis") ~ "Protein processing / Metabolism",
    
    TRUE ~ "Others"
  ))

# --------------------------------
# Reactome: 
# --------------------------------
df_reactome <- df_reactome %>%
  mutate(category = case_when(
    description %in% c("Apoptosis", "Intrinsic Pathway for Apoptosis",
                       "Activation of BAD and translocation to mitochondria",
                       "Formation of apoptosome", "Cytochrome c-mediated apoptotic response",
                       "Apoptotic cleavage of cellular proteins", "Apoptotic factor-mediated response") ~ "Apoptosis / Cell death",
    
    description %in% c("Translesion synthesis by REV1",
                       "Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template",
                       "Recognition of DNA damage by PCNA-containing replication complex",
                       "Translesion Synthesis by POLH",
                       "Recognition and association of DNA glycosylase with site containing an affected pyrimidine",
                       "Cleavage of the damaged pyrimidine",
                       "Recognition and association of DNA glycosylase with site containing an affected purine",
                       "Cleavage of the damaged purine",
                       "Resolution of AP sites via the multiple-nucleotide patch replacement pathway") ~ "DNA repair / Gene maintenance",
    
    description %in% c("SLBP independent Processing of Histone Pre-mRNAs",
                       "Transcriptional Regulation by TP53",
                       "TP53 Regulates Metabolic Genes",
                       "Regulation of MECP2 expression and activity",
                       "Regulation of TP53 Activity through Phosphorylation",
                       "RUNX1 regulates transcription of genes involved in differentiation of HSCs",
                       "DNA methylation") ~ "Gene regulation / DNA damage response",
    
    description %in% c("PI3K Cascade", "MAPK3 (ERK1) activation",
                       "Chk1/Chk2(Cds1) mediated inactivation of Cyclin B:Cdk1 complex",
                       "Mitotic Prometaphase", "Resolution of Sister Chromatid Cohesion",
                       "EML4 and NUDC in mitotic spindle formation", "RHO GTPases Activate Formins",
                       "Signaling by VEGF", "Post NMDA receptor activation events") ~ "Cell signaling / Growth regulation",
    
    description %in% c("Expression and translocation of olfactory receptors",
                       "Olfactory Signaling Pathway",
                       "Post NMDA receptor activation events") ~ "Neuronal / Sensory signaling",
    
    description %in% c("Interleukin-6 signaling",
                       "SARS-CoV-2 targets host intracellular signalling and regulatory pathways",
                       "SARS-CoV Infections", "SARS-CoV-2 Infection",
                       "Defective pyroptosis") ~ "Immune response / Inflammation",
    
    TRUE ~ "Others"
  ))


# Ensure pathways are ordered by category and significance
df_kegg <- df_kegg %>%
  arrange(category, padj) %>% 
  mutate(description = factor(description, levels = unique(description)))

df_reactome <- df_reactome %>%
  arrange(category, padj) %>% 
  mutate(description = factor(description, levels = unique(description)))


# Create color palette for categories
n_cat_kegg <- length(unique(df_kegg$category))
category_colors_kegg <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_cat_kegg)
names(category_colors_kegg) <- unique(df_kegg$category)

n_cat_reactome <- length(unique(df_reactome$category))
category_colors_reactome <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_cat_reactome)
names(category_colors_reactome) <- unique(df_reactome$category)

# Ensure contrast factor levels
df_kegg$contrast <- factor(df_kegg$contrast, levels = c("IF", "IM"))
df_reactome$contrast <- factor(df_reactome$contrast, levels = c("IF", "IM"))

png("top_terms_KEGG.png", width = 3000, height = 3000, res = 300)

ggplot(df_kegg, aes(y = description)) +
  # -------------------------------
geom_tile(aes(x = "Category", fill = category), width = 0.9, height = 0.9, show.legend = TRUE) +
  scale_fill_manual(values = category_colors_kegg, name = "Category") +
  
  # -------------------------------
ggnewscale::new_scale_fill() +
  geom_point(aes(x = contrast, fill = lor, shape = sig), color = "grey", size = 3, stroke = 0.8) +
  scale_shape_manual(values = c("sig" = 21, "no_sig" = 24)) +
  scale_fill_gradient2(
    low = "deepskyblue3",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "LOR"
  ) +
  
  # -------------------------------
theme_minimal() +
  labs(
    x = "Comparison",
    y = "Pathway",
    shape = "Significant"
  ) +
  ggtitle("Top KEGG pathways") +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
dev.off()

png("top_terms_Reactome.png", width = 3000, height = 3000, res = 300)
# Reactome 
ggplot(df_reactome, aes(y = description)) +
  # -------------------------------

geom_tile(aes(x = "Category", fill = category), width = 0.9, height = 0.9, show.legend = TRUE) +
  scale_fill_manual(values = category_colors_reactome, name = "Category") +
  
  # -------------------------------
ggnewscale::new_scale_fill() +

  geom_point(aes(x = contrast, fill = lor, shape = sig), color = "grey", size = 3, stroke = 0.8) +
  scale_shape_manual(values = c("sig" = 21, "no_sig" = 24)) +
  scale_fill_gradient2(
    low = "deepskyblue3",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "LOR"
  ) +
  
  # -------------------------------
theme_minimal() +
  labs(
    x = "Comparison",
    y = "Pathway",
    shape = "Significant"
  ) +
  ggtitle("Top Reactome pathways") +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
dev.off()
