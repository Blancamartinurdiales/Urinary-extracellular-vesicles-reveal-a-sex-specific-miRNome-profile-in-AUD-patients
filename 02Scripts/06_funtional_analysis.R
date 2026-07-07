###############################################
# Project: miRNome Profiling of Urinary EVs in AUD Patients
# Script: 06_Functional_Analysis_miRNA.R
# Author: Blanca Martin Urdiales
# Description: Functional enrichment analysis (GSEA)
#              using miRNA-target interactions filtered by confidence
###############################################

### Objective:
# The goal of this script is to reduce the number of significant miRNAs
# by restricting the number of targets per miRNA through increasing confidence.
# Only miRNA–gene interactions supported by more than one database
# will be retained as high-confidence targets.
# Based on this filtered dataset, GSEA analyses will be performed.


#### DATA PREPARATION ####
pacman::p_load(multiMiR, mdgsa, readxl, org.Hs.eg.db, GO.db, reactome.db, KEGGREST, KEGG.db, dplyr)

# Connect to GO.db database
con <- GO.db::GO_dbconn()
go_term <- tbl(con, "go_term")
go_db = as.data.frame(go_term)

# Load TarBase database and remove entries with missing values
tarbase_database <- read.table(
  file = "Data/Homo_sapiens_TarBase-v9.tsv.gz",
  header = TRUE,
  sep = "\t"
)

tarbase_database <- tarbase_database[!is.na(tarbase_database$gene_name),] # 4674235
# tarbase_database = tarbase_database[!tarbase_database$gene_name == "",]

tarbase_database <- tarbase_database[!is.na(tarbase_database$mirna_name),] # 4674235
# tarbase_database = tarbase_database[!tarbase_database$mirna_name == "",]

tarbase_database <- tarbase_database[!is.na(tarbase_database$species),] # 4674235
# tarbase_database = tarbase_database[!tarbase_database$species == "",]

# Create a combined miRNA–gene identifier (miRNA_gene)
tarbase_database$combi = paste0(tarbase_database$mirna_name,"_",tarbase_database$gene_name)

# Load miRTarBase and create the same combined identifier
mirtarbase_database <- read.csv("data/hsa_MTI.csv", header = TRUE)
colnames(mirtarbase_database)
mirtarbase_database$combi = paste0(mirtarbase_database$miRNA,"_",mirtarbase_database$`Target.Gene`)

# Extract unique miRNAs from both databases
tarbase_micro = unique(tarbase_database$mirna_name) # 1230 micros
mirtarbase_micro = unique(mirtarbase_database$miRNA) # 3094 micros

# ========================
# GLOBAL PARAMETERS
# ========================

combinar = FALSE    # combinar:
# FALSE → Use union of miRNA–gene interactions (more inclusive)
# TRUE  → Keep only interactions present in both databases (more stringent)

analisis = "GSEA" # analisis:
# "GSEA" → ranking-based enrichment
# "ORA"  → over-representation analysis

universo = TRUE   # si utilizo ORA debería ser TRUE pero en GSEA me da igual  
# Solo tiene efecto en ORA. En GSEA puedo dejarlo TRUE sin problema.

estudio = "sex"     
# estudio:
# Defines which comparison to analyze (sex, male, female, etc.)

ontologia = "BP"    
# ontologia:
# GO ontology:
#   BP → Biological Process
#   MF → Molecular Function
#   CC → Cellular Component

corte_logfc = 0     
# corte_logfc:
# LogFC threshold used only in ORA

database = "hs"     
# database:
# "hs" → org.Hs.eg.db
# "biomart" → biomaRt

# Load differential expression results depending on selected study

#estudio = "sex" 
#estudio = "male" 
estudio = "female" 
#estudio = "AUD" 

if (estudio == "sex") {
  df_res <- readRDS("Results/df_IS.rds")
} else if (estudio == "male") {
  df_res <- readRDS("Results/df_IM.rds")
} else if (estudio == "female") {
  df_res <- readRDS("Results/df_IF.rds")
} else if (estudio == "AUD") {
  df_res <- readRDS("Results/df_AUD.rds")
}

colnames(df_res) <- c("mirna", "logFC", "padj")


mirna <- df_res$mirna 

# Automatically select significant miRNAs based on logFC and adjusted p-value
mirna_sig <- df_res %>%
  dplyr::filter(abs(logFC) > 0.5, padj < 0.05) %>%
  dplyr::pull(mirna)

length(mirna_sig)
mirna_sig

# Strategy:
# Either:
# - Keep only interactions validated in both databases (high confidence)
# - Or combine all validated interactions (broader coverage)

# Quality control:
# Inspect gene symbols for potential problematic entries
# (e.g., pseudogenes or malformed identifiers)

# GSEA section:
# Transform miRNA-level statistics into gene-level rankings
# based on miRNA-target relationships

# Functional enrichment:
# Map genes to GO, Reactome, and KEGG pathways

if (combinar == TRUE){
  # # combine only the common terms
  union = merge(tarbase_database, mirtarbase_database, by = "combi") 
  colnames(union)[3] = "mature_mirna_id"
  colnames(union)[25] = "mature_mirna_id"
  colnames(union)[5] = "target_symbol"
  colnames(union)[27] = "target_symbol"
  union = union[!union$target_symbol == "",] 
  union = union[union$mature_mirna_id %in% mirna,] 
  mirna.union = unique(union$combi)
} else {
  # combine all
  a = tarbase_database[,c(2,4,23)]
  colnames(a) = c("mature_mirna_id", "target_symbol", "combi")
  b = mirtarbase_database[,c(2,4,10)]
  colnames(b) = c("mature_mirna_id", "target_symbol", "combi")
  union = rbind(a,b) 
  union = union[!union$target_symbol == "",] 
  union = union[union$mature_mirna_id %in% mirna,] 
  mirna.union = unique(union$combi)
}

# verify data quality
x = unique(union$target_symbol)
save(x, file = "Results/x_IF.RData")

save(x, file = "Results/x_IM.RData")

patron <- "\\b\\w*P(?=\\d)"

y = x[grep(patron, x, perl=TRUE, value=FALSE)]
y

length(unique(mirna.union))
# 1415899 interactions miRNA–gen
# The object "union" contents all interactions needed for enrichtment analysis.
unique.mirna.union <- unique(union)

### GSEA ###
analisis = "GSEA"

# Extract the estatics values used by mdgsa for the ranking
pvalue <- as.numeric(df_res$padj) # pvalue: vector with adjusted p-values of each miRNA.
names(pvalue) = df_res$mirna 
statistic <- as.numeric(df_res$logFC) # vector with log2FoldChange (change direction).
names(statistic) = df_res$mirna

# generating a new index using pvalue and statistic from miRNA differential expression
rindex0 <- rindexT <- rindex <- list () # miRNAs ranking 
rindex0[[estudio]] <- pval2index(pval = pvalue, sign = statistic) 
rindex0[[2]] <- NULL
rindexT[[2]] <- NULL


# funcTion that extract for each miRNA the associated genes.
tomate = function(mirna){
  a = unique(subset(union, union$mature_mirna_id == mirna, select = target_symbol, drop = TRUE))
  return(a)
} #  target genes of each miRNA.

# Generated a list: for each miRNA we have the associated target genes
tomatito <- lapply(names(rindex0[[estudio]]), tomate)
names(tomatito) = names(rindex0[[estudio]])


# miRNA ranking to genes ranking 

rindexT[[estudio]] <- transferIndex(index = rindex0[[estudio]], # miRNA Ranking
                                    targets = tomatito, # miRNA and targets association
                                    method = "average") #  sum (more terms) or average (less terms)


# Normalization of data distribution
rindex[[estudio]] <- indexTransform(index = rindexT[[estudio]], 
                                    method = "normalize")

genes = unlist(names(rindex[[estudio]]))

error.genes = rindex[[estudio]][grep("^[a-zA-Z0-9.-]+$",genes)] 
diferencia = setdiff(genes, names(error.genes))

for (i in diferencia){
  rindex[[estudio]] = rindex[[estudio]][-(which(names(rindex[[estudio]]) == i))]
}



#### FUNCTIONAL ENRICHMENT ####


## OPTION A: BIOMART ##

if (database == "biomart"){
  library(biomaRt)
  # Pasamos a cargar BIOMART
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # Conectamos con biomart
  
  go_terms <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","go_id", "name_1006", "namespace_1003"),
                    filter = 'hgnc_symbol',
                    values =  names(rindex[[estudio]]), 
                    mart = mart)
  
  go_terms <- go_terms[!(go_terms$go_id == ""),  ]
  

  malitos <- go_terms[which(go_terms$namespace_1003 == ""), 3]
  
  library(dplyr)
  con <- GO.db::GO_dbconn()
  go_term <- tbl(con, "go_term")
  buenitos = as.data.frame(go_term)
  
  arregladitos = buenitos[buenitos$go_id %in% malitos,c(2,3,4)]
  arregladitos$ontology = gsub("MF", "molecular_function", arregladitos$ontology)
  arregladitos$ontology = gsub("BP", "biological_process", arregladitos$ontology)
  arregladitos$ontology = gsub("CC", "cellular_component", arregladitos$ontology)
  colnames(arregladitos) = c("go_id", "name_1006", "namespace_1003")
  
  pruebita = merge(go_terms,arregladitos, by = "go_id")
  pruebita = pruebita[,-c(4,5)]
  pruebita = pruebita[,c(2,3,1,4,5)]
  colnames(pruebita) = colnames(go_terms)
  
  go_terms <- go_terms[-which(go_terms$namespace_1003 == ""), ]
  
  go_terms = rbind(go_terms, pruebita)
  
  
  # Select Go terms type
  go_bp = go_terms[go_terms$namespace_1003 == "biological_process",]
  go_bp = go_bp[,c(1:4)]
  gentogo = split(go_bp$go_id, go_bp$hgnc_symbol)
  gotogen = split( go_bp$hgnc_symbol, go_bp$go_id)
  go_cc = go_terms[go_terms$namespace_1003 == "cellular_component",]
  go_cc = go_cc[,c(1:4)]
  gentoset_cc = split(go_cc$go_id, go_cc$hgnc_symbol)
  settogen_cc = split( go_cc$hgnc_symbol, go_cc$go_id)
  go_mf = go_terms[go_terms$namespace_1003 == "molecular_function",]
  go_mf = go_mf[,c(1:4)]
  gentoset_mf = split(go_mf$go_id, go_mf$hgnc_symbol)
  settogen_mf = split( go_mf$hgnc_symbol, go_mf$go_id)
} else { 
  ## OPTION B: org.Hs.eg.db 
  
  # ################################ GO  #########################################
  go_genes <- AnnotationDbi::select(org.Hs.eg.db, keys=names(rindex[[estudio]]),
                                    columns=c("SYMBOL", "GO"), 
                                    keytype="SYMBOL")
  
  go_term = merge(go_db, go_genes, by.x = "go_id", by.y = "GO")
  go_term_ont = go_term[go_term$ontology == ontologia,]
  
  gentogo = split(go_term_ont$go_id, go_term_ont$SYMBOL)
  gotogen = split( go_term_ont$SYMBOL, go_term_ont$go_id)
  
  # Select the GO term type
  # --- Biological Process (BP) ---
  go_bp = go_term[go_term$ontology == "BP", ]
  go_bp = go_bp[, c("go_id", "term", "definition", "SYMBOL")]
  go_bp = na.omit(go_bp)  # eliminar filas con NA
  gentogo = split(go_bp$go_id, go_bp$SYMBOL)
  gotogen = split(go_bp$SYMBOL, go_bp$go_id)
  
  # --- Cellular Component (CC) ---
  go_cc = go_term[go_term$ontology == "CC", ]
  go_cc = go_cc[, c("go_id", "term", "definition", "SYMBOL")]
  go_cc = na.omit(go_cc)
  gentoset_cc = split(go_cc$go_id, go_cc$SYMBOL)
  settogen_cc = split(go_cc$SYMBOL, go_cc$go_id)
  
  # --- Molecular Function (MF) ---
  go_mf = go_term[go_term$ontology == "MF", ]
  go_mf = go_mf[, c("go_id", "term", "definition", "SYMBOL")]
  go_mf = na.omit(go_mf)
  gentoset_mf = split(go_mf$go_id, go_mf$SYMBOL)
  settogen_mf = split(go_mf$SYMBOL, go_mf$go_id)
}

# ############################## REACTOME ########################################

genes_entrez <- na.omit(AnnotationDbi::select(org.Hs.eg.db, keys=genes,
                                              columns=c("SYMBOL", "ENTREZID"), 
                                              keytype="SYMBOL"))

reactome_db <- AnnotationDbi::select(reactome.db, keys=genes_entrez$ENTREZID,
                                     columns=c("ENTREZID", "PATHID"),
                                     keytype="ENTREZID")

reactome_db = na.omit(merge(reactome_db, genes_entrez, by.x = "ENTREZID", by.y = "ENTREZID"))
genetoreactome = split(reactome_db$PATHID, reactome_db$SYMBOL)
reactometogen = split(reactome_db$SYMBOL, reactome_db$PATHID)

#Load the package
library(reactome.db)

# Extract ID table ↔ Rute names
reactome_pathways <- as.data.frame(reactomePATHID2NAME)

colnames(reactome_pathways) <- c("codigo", "description")

# 🔍 Filter only human rutes (R-HSA-)
reactome_pathways <- reactome_pathways[grep("^R-HSA-", reactome_pathways$codigo), ]


head(reactome_pathways)


# ################################ KEGG ##########################################
pathways.list <- keggList("pathway", "hsa")
kegg_db = data.frame(codigo = names(pathways.list),
                     summary = pathways.list)

kegg_genes <- AnnotationDbi::select(org.Hs.eg.db, keys=genes,
                                    columns=c("SYMBOL", "PATH"),
                                    keytype="SYMBOL")
kegg_genes = na.omit(kegg_genes)
kegg_genes$PATH = paste0("hsa",kegg_genes$PATH)

kegg = merge(kegg_db, kegg_genes, by.x = "codigo", by.y = "PATH")

gentokegg = split(kegg$codigo, kegg$SYMBOL)
keggtogen = split(kegg$SYMBOL, kegg$codigo)

save(kegg_db, file = "Results/kegg_db.RData")
save(reactome_pathways, file = "Results/reactome_pathways.RData")

rm(kegg_genes)
rm(kegg_db)


lista_universo = list(gotogen, settogen_cc, settogen_mf, reactometogen, keggtogen)
list_of_datasets = list()
for (i in 1:length(lista_universo)){
  annot <- annotFilter(lista_universo[[i]], rindex[[estudio]])
  res.uv <- uvGsa(rindex[[estudio]], annot, p.adjust.method = "fdr")
  if (i==1 | i==2 | i==3){
    res.uv$description = getGOnames(res.uv, verbose = TRUE)
  }
  list_of_datasets[[i]] = res.uv[order(res.uv$padj),]
  resultados = rownames(uvSignif(list_of_datasets[[i]]))
  # if (length(resultados) > 0){
  #   dir.create(paste0("/plots","_",lista_universo[[1]]), showWarnings=F)
  #   ruta_4 = paste0(getwd(), "/plots")
  #   for (i in 1:length(resultados)){
  #     plotEnrichment(annot[[resultados[i]]],
  #                    rindex[[contraste]]) + labs(title=resultados[i])
  #     ggsave(paste0(ruta_4,"/", resultados[i],"_rank.png"), width=15, height=10, bg = "white")
  #   }
  # }
}

# Obtain GO terms names
list_of_datasets[[1]]$description = getGOnames(rownames(list_of_datasets[[1]]))
names(list_of_datasets) = c("BP","CC","MF","REACTOME", "KEGG")

gsea_bp = list_of_datasets[[1]]
gsea_cc= list_of_datasets[[2]]
gsea_mf = list_of_datasets[[3]]
gsea_reactome = list_of_datasets[[4]]
gsea_kegg = list_of_datasets[[5]]


# Select one: 

# gsea_org_db_Females <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
# save(file = "Results/gsea_org_db_Females.RData", gsea_org_db_Females)
# 
#gsea_org_db_Males <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
#save(file = "Results/gsea_org_db_Males.RData", gsea_org_db_Males)
# 
gsea_org_db_Sex <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
 save(file = "Results/gsea_org_db_Sex.RData", gsea_org_db_Sex)
 # 
 #gsea_org_db_Case <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
 #save(file = "Results/gsea_org_db_Case.RData", gsea_org_db_Case)
# 
# 
# gsea_biomart_Females <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
# save(file = "functional_results/gsea_biomart_Females.RData", gsea_biomart_Females)
# 
# gsea_biomart_Males <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
# save(file = "functional_results/gsea_biomart_Males.RData", gsea_biomart_Males)
# 
# gsea_biomart_Sex <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
# save(file = "functional_results/gsea_biomart_Sex.RData", gsea_biomart_Sex)
# 
# gsea_biomart_Case <- list("gsea_bp" = gsea_bp, "gsea_mf"= gsea_mf, "gsea_cc" = gsea_cc, "gsea_reactome" = gsea_reactome, "gsea_kegg" = gsea_kegg)
# save(file = "functional_results/gsea_biomart_Case.RData", gsea_biomart_Case)

