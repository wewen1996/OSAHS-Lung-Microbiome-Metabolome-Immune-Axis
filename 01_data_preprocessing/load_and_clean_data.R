# load_and_clean_data.R
# Load and clean raw microbiome and metabolome data for OSAHS analysis

# Load required libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(magrittr)
library(phyloseq)
library(microbiome)
library(vegan)
library(metagenomeSeq)
library(impute)
library(preprocessCore)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load metadata
# ------------------------------------------------------------------------------
metadata <- read.csv("data/metadata/metadata.csv", row.names = 1, stringsAsFactors = FALSE)
print("Metadata loaded successfully:")
print(dim(metadata))
print(head(metadata))

# Check for missing values in metadata
metadata_missing <- colSums(is.na(metadata))
print("Missing values in metadata:")
print(metadata_missing[metadata_missing > 0])

# ------------------------------------------------------------------------------
# 2. Load and process microbiome data
# ------------------------------------------------------------------------------
# Load OTU table (rows = OTUs, columns = samples)
microbiome_otu <- read.csv("data/raw/microbiome_otu_table.csv", row.names = 1, check.names = FALSE)
print("Microbiome OTU table loaded successfully:")
print(dim(microbiome_otu))
print(head(microbiome_otu[, 1:5]))

# Load taxonomy table
microbiome_taxonomy <- read.csv("data/raw/microbiome_taxonomy.csv", row.names = 1, stringsAsFactors = FALSE)
print("Microbiome taxonomy table loaded successfully:")
print(dim(microbiome_taxonomy))
print(head(microbiome_taxonomy))

# Load phylogenetic tree (if available)
if (file.exists("data/raw/microbiome_tree.tree")) {
  microbiome_tree <- read_tree("data/raw/microbiome_tree.tree")
  print("Phylogenetic tree loaded successfully")
} else {
  microbiome_tree <- NULL
  print("Phylogenetic tree not found")
}

# Clean microbiome data
# Remove OTUs with zero abundance across all samples
microbiome_otu <- microbiome_otu[rowSums(microbiome_otu) > 0, ]
print(paste("OTUs after removing zero abundance:", nrow(microbiome_otu)))

# Remove samples with low sequencing depth (less than 1000 reads)
sample_depth <- colSums(microbiome_otu)
print("Sample sequencing depth:")
print(summary(sample_depth))

low_depth_samples <- names(sample_depth[sample_depth < 1000])
if (length(low_depth_samples) > 0) {
  print(paste("Removing", length(low_depth_samples), "samples with low sequencing depth:"))
  print(low_depth_samples)
  microbiome_otu <- microbiome_otu[, !colnames(microbiome_otu) %in% low_depth_samples]
}

print(paste("Samples after removing low depth samples:", ncol(microbiome_otu)))

# Ensure taxonomy table matches OTU table
microbiome_taxonomy <- microbiome_taxonomy[rownames(microbiome_otu), ]

# ------------------------------------------------------------------------------
# 3. Load and process metabolome data
# ------------------------------------------------------------------------------
# Load metabolome data (rows = metabolites, columns = samples)
metabolome_data <- read.csv("data/raw/metabolome_data.csv", row.names = 1, check.names = FALSE)
print("Metabolome data loaded successfully:")
print(dim(metabolome_data))
print(head(metabolome_data[, 1:5]))

# Load metabolite annotations
metabolite_annotations <- read.csv("data/raw/metabolite_annotations.csv", row.names = 1, stringsAsFactors = FALSE)
print("Metabolite annotations loaded successfully:")
print(dim(metabolite_annotations))
print(head(metabolite_annotations))

# Clean metabolome data
# Check for missing values
metabolome_missing <- sum(is.na(metabolome_data))
print(paste("Total missing values in metabolome data:", metabolome_missing))
print(paste("Percentage missing values:", round(metabolome_missing / (nrow(metabolome_data) * ncol(metabolome_data)) * 100, 2), "%"))

# Remove metabolites with more than 50% missing values
metabolite_missing_percentage <- rowMeans(is.na(metabolome_data))
high_missing_metabolites <- names(metabolite_missing_percentage[metabolite_missing_percentage > 0.5])
if (length(high_missing_metabolites) > 0) {
  print(paste("Removing", length(high_missing_metabolites), "metabolites with >50% missing values"))
  metabolome_data <- metabolome_data[!rownames(metabolome_data) %in% high_missing_metabolites, ]
  metabolite_annotations <- metabolite_annotations[!rownames(metabolite_annotations) %in% high_missing_metabolites, ]
}

print(paste("Metabolites after removing high missing value metabolites:", nrow(metabolome_data)))

# Remove samples with more than 30% missing values
sample_missing_percentage <- colMeans(is.na(metabolome_data))
high_missing_samples <- names(sample_missing_percentage[sample_missing_percentage > 0.3])
if (length(high_missing_samples) > 0) {
  print(paste("Removing", length(high_missing_samples), "samples with >30% missing values"))
  metabolome_data <- metabolome_data[, !colnames(metabolome_data) %in% high_missing_samples]
}

print(paste("Samples after removing high missing value samples:", ncol(metabolome_data)))

# ------------------------------------------------------------------------------
# 4. Align samples across datasets
# ------------------------------------------------------------------------------
# Find common samples across all datasets
common_samples <- Reduce(intersect, list(
  rownames(metadata),
  colnames(microbiome_otu),
  colnames(metabolome_data)
))

print(paste("Common samples across all datasets:", length(common_samples)))

# Subset all datasets to common samples
metadata <- metadata[common_samples, ]
microbiome_otu <- microbiome_otu[, common_samples]
metabolome_data <- metabolome_data[, common_samples]

print("Datasets after aligning samples:")
print(paste("Metadata samples:", nrow(metadata)))
print(paste("Microbiome samples:", ncol(microbiome_otu)))
print(paste("Metabolome samples:", ncol(metabolome_data)))

# ------------------------------------------------------------------------------
# 5. Save processed data
# ------------------------------------------------------------------------------
# Save processed metadata
write.csv(metadata, "data/processed/metadata_processed.csv", row.names = TRUE)

# Save processed microbiome data
write.csv(microbiome_otu, "data/processed/microbiome_otu_processed.csv", row.names = TRUE)
write.csv(microbiome_taxonomy, "data/processed/microbiome_taxonomy_processed.csv", row.names = TRUE)
if (!is.null(microbiome_tree)) {
  write.tree(microbiome_tree, "data/processed/microbiome_tree_processed.tree")
}

# Save processed metabolome data
write.csv(metabolome_data, "data/processed/metabolome_data_processed.csv", row.names = TRUE)
write.csv(metabolite_annotations, "data/processed/metabolite_annotations_processed.csv", row.names = TRUE)

print("All processed data saved successfully!")

# ------------------------------------------------------------------------------
# 6. Summary statistics
# ------------------------------------------------------------------------------
cat("\n=== Summary Statistics ===\n")
cat("Metadata:\n")
cat(paste("  Samples:", nrow(metadata), "\n"))
cat(paste("  Variables:", ncol(metadata), "\n"))

cat("\nMicrobiome Data:\n")
cat(paste("  OTUs:", nrow(microbiome_otu), "\n"))
cat(paste("  Samples:", ncol(microbiome_otu), "\n"))
cat(paste("  Total reads:", sum(microbiome_otu), "\n"))
cat(paste("  Mean reads per sample:", round(mean(colSums(microbiome_otu))), "\n"))

cat("\nMetabolome Data:\n")
cat(paste("  Metabolites:", nrow(metabolome_data), "\n"))
cat(paste("  Samples:", ncol(metabolome_data), "\n"))
cat(paste("  Missing values:", sum(is.na(metabolome_data)), "\n"))
cat(paste("  Percentage missing values:", round(sum(is.na(metabolome_data)) / (nrow(metabolome_data) * ncol(metabolome_data)) * 100, 2), "%\n"))
