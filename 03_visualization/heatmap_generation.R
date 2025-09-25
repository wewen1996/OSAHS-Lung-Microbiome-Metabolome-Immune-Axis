# heatmap_generation.R
# Generate various types of heatmaps for microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/figures/heatmaps", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load normalized data
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load differential analysis results
microbiome_significant <- read.csv("results/tables/differential_analysis/microbiome_t_test_significant.csv", stringsAsFactors = FALSE)
metabolome_significant <- read.csv("results/tables/differential_analysis/metabolomics_oplsda_significant.csv", stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Generate microbiome heatmaps at different taxonomic levels
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome heatmaps ===\n")

# Function to generate heatmap at different taxonomic levels
generate_taxonomic_heatmap <- function(phyloseq_obj, tax_level, top_n = 30, title = "") {
  # Agglomerate taxa at the specified level
  phyloseq_agg <- tax_glom(phyloseq_obj, taxrank = tax_level)
  
  # Get abundance matrix
  otu_table <- otu_table(phyloseq_agg) %>% as.matrix()
  
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_agg) %>% as.data.frame()
  
  # Calculate mean abundance by group
  group_means <- aggregate(t(otu_table), by = list(Group = sample_data(phyloseq_agg)$Group), FUN = mean)
  rownames(group_means) <- group_means$Group
  group_means <- group_means[, -1]
  
  # Select top N taxa by variance
  taxa_variance <- apply(group_means, 2, var)
  top_taxa <- names(sort(taxa_variance, decreasing = TRUE))[1:min(top_n, length(taxa_variance))]
  
  # Subset data to top taxa
  top_data <- group_means[, top_taxa]
  
  # Create taxonomic labels
  tax_labels <- apply(taxonomy[top_taxa, ], 1, function(x) {
    paste(x[!is.na(x)], collapse = "|")
  })
  
  # Set row names to taxonomic labels
  colnames(top_data) <- tax_labels
  
  # Create heatmap
  pheatmap(
    t(top_data),
    scale = "row",
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 20,
    treeheight_col = 10,
    fontsize_row = 8,
    main = title,
    color = viridis(100),
    filename = paste0("results/figures/heatmaps/microbiome_", tolower(tax_level), "_heatmap.png"),
    width = 10,
    height = 12,
    dpi = 300
  )
  
  return(top_data)
}

# Generate heatmaps at different taxonomic levels
if (ntaxa(phyloseq_css) > 0) {
  # Phylum level
  generate_taxonomic_heatmap(
    phyloseq_css,
    "Phylum",
    top_n = 15,
    title = "Microbiome Composition at Phylum Level"
  )
  
  # Genus level
  generate_taxonomic_heatmap(
    phyloseq_css,
    "Genus",
    top_n = 30,
    title = "Microbiome Composition at Genus Level"
  )
  
  print("Microbiome taxonomic heatmaps generated successfully!")
} else {
  print("No microbiome data available for heatmap generation.")
}

# ------------------------------------------------------------------------------
# 3. Generate significant features heatmap
# ------------------------------------------------------------------------------
cat("\n=== Generating significant features heatmaps ===\n")

# Function to generate significant features heatmap
generate_significant_heatmap <- function(data, metadata, significant_ids, title, filename, feature_annotations = NULL) {
  # Subset data to significant features
  significant_data <- data[significant_ids, ]
  
  # If there are too many features, select top N by variance
  if (nrow(significant_data) > 50) {
    feature_variance <- apply(significant_data, 1, var)
    top_features <- names(sort(feature_variance, decreasing = TRUE))[1:50]
    significant_data <- significant_data[top_features, ]
    significant_ids <- top_features
  }
  
  # Create feature labels
  if (!is.null(feature_annotations)) {
    if ("Metabolite.Name" %in% colnames(feature_annotations)) {
      feature_labels <- apply(feature_annotations[significant_ids, c("Metabolite.Name", rownames(feature_annotations))], 1, function(x) {
        if (!is.na(x[1]) && x[1] != "") x[1] else x[2]
      })
    } else if ("Genus" %in% colnames(feature_annotations)) {
      feature_labels <- apply(feature_annotations[significant_ids, c("Phylum", "Genus")], 1, function(x) {
        paste(x[!is.na(x)], collapse = "|")
      })
    } else {
      feature_labels <- significant_ids
    }
    rownames(significant_data) <- feature_labels
  }
  
  # Create sample annotations
  sample_annot <- data.frame(Group = metadata$Group)
  rownames(sample_annot) <- colnames(significant_data)
  
  # Create heatmap
  pheatmap(
    significant_data,
    scale = "row",
    annotation_col = sample_annot,
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = title,
    color = viridis(100),
    filename = filename,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  return(significant_data)
}

# Generate significant microbiome features heatmap
if (nrow(microbiome_significant) > 0) {
  # Get microbiome data
  microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()
  
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  # Generate heatmap
  generate_significant_heatmap(
    microbiome_data,
    metadata,
    microbiome_significant$Taxon,
    title = "Significant Microbiome Features (FDR < 0.05)",
    filename = "results/figures/heatmaps/microbiome_significant_heatmap.png",
    feature_annotations = taxonomy
  )
  
  print("Significant microbiome features heatmap generated successfully!")
} else {
  print("No significant microbiome features found. Skipping heatmap generation.")
}

# Generate significant metabolome features heatmap
if (nrow(metabolome_significant) > 0) {
  # Generate heatmap
  generate_significant_heatmap(
    metabolome_data,
    metadata,
    metabolome_significant$Metabolite,
    title = "Significant Metabolite Features (VIP > 1, FDR < 0.05)",
    filename = "results/figures/heatmaps/metabolome_significant_heatmap.png",
    feature_annotations = metabolite_annotations
  )
  
  print("Significant metabolome features heatmap generated successfully!")
} else {
  print("No significant metabolome features found. Skipping heatmap generation.")
}

# ------------------------------------------------------------------------------
# 4. Generate correlation heatmap between microbiome and metabolome
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome-metabolome correlation heatmap ===\n")

# Function to generate correlation heatmap
generate_correlation_heatmap <- function(microbiome_data, metabolome_data, metadata, top_microbiome = 20, top_metabolite = 20, title = "") {
  # Ensure samples are in the same order
  common_samples <- intersect(colnames(microbiome_data), colnames(metabolome_data))
  
  if (length(common_samples) < 3) {
    print("Not enough common samples for correlation analysis.")
    return(NULL)
  }
  
  # Subset data to common samples
  microbiome_data_sub <- microbiome_data[, common_samples]
  metabolome_data_sub <- metabolome_data[, common_samples]
  
  # Select top microbiome features by variance
  microbe_variance <- apply(microbiome_data_sub, 1, var)
  top_microbes <- names(sort(microbe_variance, decreasing = TRUE))[1:min(top_microbiome, length(microbe_variance))]
  microbiome_data_sub <- microbiome_data_sub[top_microbes, ]
  
  # Select top metabolite features by variance
  metabolite_variance <- apply(metabolome_data_sub, 1, var)
  top_metabolites <- names(sort(metabolite_variance, decreasing = TRUE))[1:min(top_metabolite, length(metabolite_variance))]
  metabolome_data_sub <- metabolome_data_sub[top_metabolites, ]
  
  # Calculate correlation matrix
  correlation_matrix <- cor(t(microbiome_data_sub), t(metabolome_data_sub), method = "spearman")
  
  # Create heatmap using ComplexHeatmap
  # Create row annotations (microbiome phylum)
  if (!is.null(taxonomy)) {
    microbe_annot <- taxonomy[top_microbes, "Phylum"]
    names(microbe_annot) <- top_microbes
    
    # Create column annotations (metabolite class if available)
    if (!is.null(metabolite_annotations) && "Class" %in% colnames(metabolite_annotations)) {
      metabolite_annot <- metabolite_annotations[top_metabolites, "Class"]
      names(metabolite_annot) <- top_metabolites
    } else {
      metabolite_annot <- NULL
    }
    
    # Create annotation colors
    if (!is.null(microbe_annot)) {
      phylum_colors <- setNames(brewer.pal(min(12, length(unique(microbe_annot))), "Set3"), unique(microbe_annot))
      row_annot <- rowAnnotation(Phylum = microbe_annot, col = list(Phylum = phylum_colors))
    } else {
      row_annot <- NULL
    }
    
    if (!is.null(metabolite_annot)) {
      class_colors <- setNames(brewer.pal(min(12, length(unique(metabolite_annot))), "Set2"), unique(metabolite_annot))
      col_annot <- columnAnnotation(Class = metabolite_annot, col = list(Class = class_colors))
    } else {
      col_annot <- NULL
    }
  } else {
    row_annot <- NULL
    col_annot <- NULL
  }
  
  # Create heatmap
  png(paste0("results/figures/heatmaps/microbiome_metabolome_correlation.png"), width = 1200, height = 1000, res = 300)
  Heatmap(
    correlation_matrix,
    name = "Spearman Correlation",
    row_title = "Microbiome Features",
    column_title = "Metabolite Features",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    top_annotation = col_annot,
    left_annotation = row_annot,
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Correlation Coefficient",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  )
  dev.off()
  
  return(correlation_matrix)
}

# Generate microbiome-metabolome correlation heatmap
if (ntaxa(phyloseq_css) > 0 && nrow(metabolome_data) > 0) {
  # Get microbiome data
  microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()
  
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  # Generate correlation heatmap
  generate_correlation_heatmap(
    microbiome_data,
    metabolome_data,
    metadata,
    top_microbiome = 30,
    top_metabolite = 30,
    title = "Microbiome-Metabolome Correlation Heatmap"
  )
  
  print("Microbiome-metabolome correlation heatmap generated successfully!")
} else {
  print("Not enough data for microbiome-metabolome correlation analysis.")
}

# ------------------------------------------------------------------------------
# 5. Generate clinical metadata correlation heatmap
# ------------------------------------------------------------------------------
cat("\n=== Generating clinical metadata correlation heatmap ===\n")

# Function to generate clinical metadata correlation heatmap
generate_metadata_correlation_heatmap <- function(metadata, title = "") {
  # Select numeric variables
  numeric_metadata <- metadata %>%
    select_if(is.numeric) %>%
    select(-any_of(c("SampleID", "PatientID")))  # Exclude ID columns
  
  if (ncol(numeric_metadata) < 2) {
    print("Not enough numeric metadata variables for correlation analysis.")
    return(NULL)
  }
  
  # Calculate correlation matrix
  correlation_matrix <- cor(numeric_metadata, method = "pearson")
  
  # Create heatmap
  pheatmap(
    correlation_matrix,
    scale = "none",
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 10,
    fontsize_col = 10,
    main = title,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    filename = "results/figures/heatmaps/metadata_correlation.png",
    width = 10,
    height = 10,
    dpi = 300
  )
  
  return(correlation_matrix)
}

# Generate metadata correlation heatmap
if (ncol(metadata) > 1) {
  generate_metadata_correlation_heatmap(
    metadata,
    title = "Clinical Metadata Correlation Heatmap"
  )
  
  print("Clinical metadata correlation heatmap generated successfully!")
} else {
  print("Not enough metadata variables for correlation analysis.")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== heatmap_generation.R Summary ===\n")
cat("\nGenerated heatmaps:\n")
cat("1. Microbiome taxonomic composition heatmaps (Phylum and Genus levels)\n")

if (nrow(microbiome_significant) > 0) {
  cat("2. Significant microbiome features heatmap\n")
}

if (nrow(metabolome_significant) > 0) {
  cat("3. Significant metabolome features heatmap\n")
}

if (ntaxa(phyloseq_css) > 0 && nrow(metabolome_data) > 0) {
  cat("4. Microbiome-metabolome correlation heatmap\n")
}

if (ncol(metadata) > 1) {
  cat("5. Clinical metadata correlation heatmap\n")
}

cat("\nAll heatmaps saved to: results/figures/heatmaps/\n")
