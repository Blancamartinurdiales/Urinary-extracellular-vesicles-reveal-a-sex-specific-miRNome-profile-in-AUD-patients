###############################################
# Project: miRNome profiling of urinary EVs in AUD patients
# Script: 01_preprocessing_counts.R
# Author: Blanca Martin Urdiales
# Description: Reads and preprocesses raw count data and metadata
#              for downstream differential expression analysis.
###############################################

# --- Load required packages ---
library(edgeR)
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)

# --- Define directories ---
metadata_path <- "metadata/Metadata_miRNA.csv"
data_dir <- "Data/unique_counts_1"

# --- Load and merge raw count data ---
files <- list.files(data_dir, pattern = "unique_counts.txt$", full.names = TRUE)

# Read and merge count tables into a DGEList object
y <- edgeR::readDGE(files = files)  # Read and merge count tables

# Extract sample names from file names
sample_names <- gsub("_.*", "", basename(files))

# Assign sample names to counts and samples
rownames(y$samples) <- sample_names
colnames(y$counts) <- sample_names

# --- Load sample metadata ---
metadata <- read.csv(metadata_path)
metadata <- metadata[order(metadata$ID), ]

# Match metadata with count data
rownames(y$samples) <- metadata$Sample
colnames(y$counts) <- metadata$Sample
y$counts <- y$counts[, order(colnames(y$counts))]

# --- Prepare count matrix and metadata ---
cts <- as.data.frame(y$counts)

# Remove unwanted sample (example: CTROL6_M -> QC failure)
cts$CTROL6_M <- NULL
metadata2 <- metadata[metadata$Sample != "CTROL6_M", ]

# Prepare DESeq2 input
coldata <- metadata2
rownames(coldata) <- coldata$Sample
coldata <- coldata[order(coldata$Sample), ]

all(colnames(cts) == rownames(coldata))  # Should return TRUE

# Add Condition2 (without sex)
coldata$Condition2 <- ifelse(grepl("AUD", coldata$Condition), "AUD", "Control")
coldata$Condition2 <- factor(coldata$Condition2, levels=c("Control","AUD"))

# Save coldata as RDS
saveRDS(coldata, file = "metadata/coldata.rds")

# --- Create DESeq2 object ---
dds1 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 0 + Condition)

# --- Save DESeq2 object (optional) ---
# Global comparison: AUD vs Control
saveRDS(dds1, file = "Results/dds1_raw_counts.rds")


