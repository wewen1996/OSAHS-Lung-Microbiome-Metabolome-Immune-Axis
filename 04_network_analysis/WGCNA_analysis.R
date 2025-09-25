# WGCNA_analysis.R
# Perform WGCNA analysis on microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(WGCNA)
library(flashClust)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(reshape2)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/networks/WGCNA", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/networks/WGCNA", showWarnings = FALSE, recursive = TRUE)

# Enable multi-threading (optional)
enableWGCNAThreads()

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load normalizedized data
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Prepare data for WGCNA
# ------------------------------------------------------------------------------
# Function to prepare data for WGCNA
prepare_WGCNA_data <- function(data, metadata, min_variance = 0.01, min_prevalence = 0.1) {
  
  # Ensure samples are in the same order
  common_samples <- intersect(colnames(data), rownames(metadata))
  data <- data[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Filter low variance features
  feature_variance <- apply(data, 1, var)
  features_to_keep <- feature_variance >= quantile(feature_variance, min_variance)
  
  # Filter low prevalence features
  if (min_prevalence > 0) {
    feature_prevalence <- rowSums(data > 0) / ncol(data)
    features_to_keep <- features_to_keep & (feature_prevalence >= min_prevalence)
  }
  
  data_filtered <- data[features_to_keep, ]
  
  print(paste("Features after filtering:", sum(features_to_keep)))
  print(paste("Original features:", nrow(data)))
  
  # Transpose data (samples as rows, features as columns)
  data_t <- t(data_filtered)
  
  # Convert to numeric matrix
  data_matrix <- as.matrix(data_t)
  
  # Check for missing values
  if (any(is.na(data_matrix))) {
    print("Missing values detected. Imputing using k-nearest neighbors...")
    data_matrix <- impute.knn(data_matrix)$data
  }
  
  return(list(
    data = data_matrix,
    metadata = metadata
  ))
}

# ------------------------------------------------------------------------------
# 3. Perform WGCNA analysis for microbiome data
# ------------------------------------------------------------------------------
cat("\n=== Performing WGCNA analysis for microbiome data ===\n")

# Get microbiome data
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Prepare data for WGCNA
microbiome_WGCNA_data <- prepare_WGCNA_data(
  data = microbiome_data,
  metadata = metadata,
  min_variance = 0.05,
  min_prevalence = 0.2
)

# Check if we have enough features for WGCNA
if (ncol(microbiome_WGCNA_data$data) >= 20) {
  
  # Step 1: Choose soft-thresholding power
  cat("\nStep 1: Choosing soft-thresholding power...\n")
  
  # Set power range
  power_range <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  
  # Calculate scale-free topology fit indices
  sft <- pickSoftThreshold(
    microbiome_WGCNA_data$data,
    powerVector = power_range,
    networkType = "unsigned",
    verbose = 5
  )
  
  # Plot scale-free topology fit indices
  p_sft <- ggplot(data = data.frame(
    Power = sft$fitIndices$Power,
    SFT.R.sq = sft$fitIndices$SFT.R.sq,
    Mean.K = sft$fitIndices$mean.k.,
    Median.K = sft$fitIndices$median.k.
  )) +
    geom_point(aes(x = Power, y = SFT.R.sq), color = "red") +
    geom_line(aes(x = Power, y = SFT.R.sq), color = "red") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
    labs(
      title = "Scale-free Topology Fit Index vs Soft-thresholding Power",
      x = "Soft-thresholding Power",
      y = "Scale-free Topology Fit Index (R^2)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  ggsave("results/figures/networks/WGCNA/microbiome_sft_plot.png", p_sft,
         width = 10, height = 8, dpi = 300)
  
  # Choose soft-thresholding power
  soft_power <- sft$powerEstimate
  
  if (is.na(soft_power)) {
    soft_power <- 6  # Default value if no power is estimated
    print(paste("No soft-thresholding power estimated. Using default value:", soft_power))
  } else {
    print(paste("Chosen soft-thresholding power:", soft_power))
  }
  
  # Step 2: Construct co-expression network
  cat("\nStep 2: Constructing co-expression network...\n")
  
  # Calculate adjacency matrix
  adjacency <- adjacency(
    microbiome_WGCNA_data$data,
    power = soft_power,
    type = "unsigned"
  )
  
  # Convert adjacency to topological overlap matrix (TOM)
  tom <- TOMsimilarity(adjacency)
  diss_tom <- 1 - tom
  
  # Step 3: Identify modules
  cat("\nStep 3: Identifying modules...\n")
  
  # Perform hierarchical clustering
  gene_tree <- flashClust(as.dist(diss_tom), method = "average")
  
  # Plot dendrogram
  png("results/figures/networks/WGCNA/microbiome_gene_dendrogram.png", width = 1200, height = 800, res = 300)
  plot(gene_tree, xlab = "", sub = "", main = "Microbiome Gene Dendrogram", labels = FALSE, hang = 0.04)
  dev.off()
  
  # Identify modules using dynamic tree cut
  min_module_size <- 20
  dynamic_mods <- cutreeDynamic(
    dendro = gene_tree,
    distM = diss_tom,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
  )
  
  # Convert numeric labels to colors
  dynamic_colors <- labels2colors(dynamic_mods)
  
  # Plot dendrogram with module colors
  png("results/figures/networks/WGCNA/microbiome_module_dendrogram.png", width = 1200, height = 800, res = 300)
  plotDendroAndColors(
    gene_tree,
    dynamic_colors,
    "Module Colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  # Step 4: Merge similar modules
  cat("\nStep 4: Merging similar modules...\n")
  
  # Calculate module eigengenes
  me_list <- moduleEigengenes(microbiome_WGCNA_data$data, colors = dynamic_colors)
  mes <- me_list$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  me_diss <- 1 - cor(mes)
  
  # Cluster module eigengenes
  me_tree <- flashClust(as.dist(me_diss), method = "average")
  
  # Plot module eigengene dendrogram
  png("results/figures/networks/WGCNA/microbiome_me_dendrogram.png", width = 1000, height = 600, res = 300)
  plot(me_tree, main = "Microbiome Module Eigengene Dendrogram", xlab = "", sub = "")
  abline(h = 0.25, col = "red")  # Threshold for merging modules
  dev.off()
  
  # Merge modules
  merge_threshold <- 0.25
  merged <- mergeCloseModules(
    microbiome_WGCNA_data$data,
    dynamic_colors,
    cutHeight = merge_threshold,
    verbose = 3
  )
  
  # Get merged module colors and eigengenes
  merged_colors <- merged$colors
  merged_mes <- merged$newMEs
  
  # Step 5: Relate modules to traits
  cat("\nStep 5: Relating modules to traits...\n")
  
  # Create trait matrix
  if ("Group" %in% colnames(microbiome_WGCNA_data$metadata)) {
    # Convert group to numeric (0 for control, 1 for case)
    trait <- ifelse(microbiome_WGCNA_data$metadata$Group == "Control", 0, 1)
    names(trait) <- rownames(microbiome_WGCNA_data$metadata)
    
    # Create trait matrix
    trait_matrix <- matrix(trait, ncol = 1)
    rownames(trait_matrix) <- rownames(microbiome_WGCNA_data$metadata)
    colnames(trait_matrix) <- "Disease"
    
    # Calculate correlation between module eigengenes and traits
    module_trait_cor <- cor(merged_mes, trait_matrix, use = "p")
    module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(microbiome_WGCNA_data$data))
    
    # Plot module-trait relationships
    png("results/figures/networks/WGCNA/microbiome_module_trait_relationship.png", width = 800, height = 1000, res = 300)
    labeledHeatmap(
      Matrix = module_trait_cor,
      xLabels = colnames(trait_matrix),
      yLabels = names(merged_mes),
      ySymbols = names(merged_mes),
      colorLabels = FALSE,
      colors = greenWhiteRed(50),
      textMatrix = paste(signif(module_trait_cor, 2), "\n(",
                         signif(module_trait_pvalue, 1), ")", sep = ""),
      setStdMargins = FALSE,
      cex.text = 0.8,
      zlim = c(-1, 1),
      main = paste("Microbiome Module-Trait Relationships")
    )
    dev.off()
  }
  
  # Step 6: Save results
  cat("\nStep 6: Saving results...\n")
  
  # Create results list
  microbiome_WGCNA_results <- list(
    soft_power = soft_power,
    adjacency = adjacency,
    tom = tom,
    gene_tree = gene_tree,
    dynamic_colors = dynamic_colors,
    merged_colors = merged_colors,
    me_list = me_list,
    merged_mes = merged_mes,
    module_trait_cor = if (exists("module_trait_cor")) module_trait_cor else NULL,
    module_trait_pvalue = if (exists("module_trait_pvalue")) module_trait_pvalue else NULL,
    trait_matrix = if (exists("trait_matrix")) trait_matrix else NULL
  )
  
  # Save results
  saveRDS(microbiome_WGCNA_results, "results/networks/WGCNA/microbiome_WGCNA_results.rds")
  
  # Save module assignments
  module_assignments <- data.frame(
    Feature = colnames(microbiome_WGCNA_data$data),
    Dynamic_Module = dynamic_colors,
    Merged_Module = merged_colors,
    stringsAsFactors = FALSE
  )
  
  # Add taxonomy information
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  module_assignments <- module_assignments %>%
    left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature")
  
  write.csv(module_assignments, "results/networks/WGCNA/microbiome_module_assignments.csv", row.names = FALSE)
  
  print("Microbiome WGCNA analysis completed successfully!")
} else {
  print("Not enough features after filtering. Skipping microbiome WGCNA analysis.")
  microbiome_WGCNA_results <- NULL
}

# ------------------------------------------------------------------------------
# 4. Perform WGCNA analysis for metabolome data
# ------------------------------------------------------------------------------
cat("\n=== Performing WGCNA analysis for metabolome data ===\n")

# Prepare data for WGCNA
metabolome_WGCNA_data <- prepare_WGCNA_data(
  data = metabolome_data,
  metadata = metadata,
  min_variance = 0.05,
  min_prevalence = 0.1
)

# Check if we have enough features for WGCNA
if (ncol(metabolome_WGCNA_data$data) >= 20) {
  
  # Step 1: Choose soft-thresholding power
  cat("\nStep 1: Choosing soft-thresholding power...\n")
  
  # Set power range
  power_range <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  
  # Calculate scale-free topology fit indices
  sft <- pickSoftThreshold(
    metabolome_WGCNA_data$data,
    powerVector = power_range,
    networkType = "unsigned",
    verbose = 5
  )
  
  # Plot scale-free topology fit indices
  p_sft <- ggplot(data = data.frame(
    Power = sft$fitIndices$Power,
    SFT.R.sq = sft$fitIndices$SFT.R.sq,
    Mean.K = sft$fitIndices$mean.k.,
    Median.K = sft$fitIndices$median.k.
  )) +
    geom_point(aes(x = Power, y = SFT.R.sq), color = "blue") +
    geom_line(aes(x = Power, y = SFT.R.sq), color = "blue") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
    labs(
      title = "Scale-free Topology Fit Index vs Soft-thresholding Power",
      x = "Soft-thresholding Power",
      y = "Scale-free Topology Fit Index (R^2)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save plot
  ggsave("results/figures/networks/WGCNA/metabolome_sft_plot.png", p_sft,
         width = 10, height = 8, dpi = 300)
  
  # Choose soft-thresholding power
  soft_power <- sft$powerEstimate
  
  if (is.na(soft_power)) {
    soft_power <- 6  # Default value if no power is estimated
    print(paste("No soft-thresholding power estimated. Using default value:", soft_power))
  } else {
    print(paste("Chosen soft-thresholding power:", soft_power))
  }
  
  # Step 2: Construct co-expression network
  cat("\nStep 2: Constructing co-expression network...\n")
  
  # Calculate adjacency matrix
  adjacency <- adjacency(
    metabolome_WGCNA_data$data,
    power = soft_power,
    type = "unsigned"
  )
  
  # Convert adjacency to topological overlap matrix (TOM)
  tom <- TOMsimilarity(adjacency)
  diss_tom <- 1 - tom
  
  # Step 3: Identify modules
  cat("\nStep 3: Identifying modules...\n")
  
  # Perform hierarchical clustering
  gene_tree <- flashClust(as.dist(diss_tom), method = "average")
  
  # Plot dendrogram
  png("results/figures/networks/WGCNA/metabolome_gene_dendrogram.png", width = 1200, height = 800, res = 300)
  plot(gene_tree, xlab = "", sub = "", main = "Metabolome Gene Dendrogram", labels = FALSE, hang = 0.04)
  dev.off()
  
  # Identify modules using dynamic tree cut
  min_module_size <- 20
  dynamic_mods <- cutreeDynamic(
    dendro = gene_tree,
    distM = diss_tom,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
  )
  
  # Convert numeric labels to colors
  dynamic_colors <- labels2colors(dynamic_mods)
  
  # Plot dendrogram with module colors
  png("results/figures/networks/WGCNA/metabolome_module_dendrogram.png", width = 1200, height = 800, res = 300)
  plotDendroAndColors(
    gene_tree,
    dynamic_colors,
    "Module Colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  # Step 4: Merge similar modules
  cat("\nStep 4: Merging similar modules...\n")
  
  # Calculate module eigengenes
  me_list <- moduleEigengenes(metabolome_WGCNA_data$data, colors = dynamic_colors)
  mes <- me_list$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  me_diss <- 1 - cor(mes)
  
  # Cluster module eigengenes
  me_tree <- flashClust(as.dist(me_diss), method = "average")
  
  # Plot module eigengene dendrogram
  png("results/figures/networks/WGCNA/metabolome_me_dendrogram.png", width = 1000, height = 600, res = 300)
  plot(me_tree, main = "Metabolome Module Eigengene Dendrogram", xlab = "", sub = "")
  abline(h = 0.25, col = "red")  # Threshold for merging modules
  dev.off()
  
  # Merge modules
  merge_threshold <- 0.25
  merged <- mergeCloseModules(
    metabolome_WGCNA_data$data,
    dynamic_colors,
    cutHeight = merge_threshold,
    verbose = 3
  )
  
  # Get merged module colors and eigengenes
  merged_colors <- merged$colors
  merged_mes <- merged$newMEs
  
  # Step 5: Relate modules to traits
  cat("\nStep 5: Relating modules to traits...\n")
  
  # Create trait matrix
  if ("Group" %in% colnames(metabolome_WGCNA_data$metadata)) {
    # Convert group to numeric (0 for control, 1 for case)
    trait <- ifelse(metabolome_WGCNA_data$metadata$Group == "Control", 0, 1)
    names(trait) <- rownames(metabolome_WGCNA_data$metadata)
    
    # Create trait matrix
    trait_matrix <- matrix(trait, ncol = 1)
    rownames(trait_matrix) <- rownames(metabolome_WGCNA_data$metadata)
    colnames(trait_matrix) <- "Disease"
    
    # Calculate correlation between module eigengenes and traits
    module_trait_cor <- cor(merged_mes, trait_matrix, use = "p")
    module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(metabolome_WGCNA_data$data))
    
    # Plot module-trait relationships
    png("results/figures/networks/WGCNA/metabolome_module_trait_relationship.png", width = 800, height = 1000, res = 300)
    labeledHeatmap(
      Matrix = module_trait_cor,
      xLabels = colnames(trait_matrix),
      yLabels = names(merged_mes),
      ySymbols = names(merged_mes),
      colorLabels = FALSE,
      colors = greenWhiteRed(50),
      textMatrix = paste(signif(module_trait_cor, 2), "\n(",
                         signif(module_trait_pvalue, 1), ")", sep = ""),
      setStdMargins = FALSE,
      cex.text = 0.8,
      zlim = c(-1, 1),
      main = paste("Metabolome Module-Trait Relationships")
    )
    dev.off()
  }
  
  # Step 6: Save results
  cat("\nStep 6: Saving results...\n")
  
  # Create results list
  metabolome_WGCNA_results <- list(
    soft_power = soft_power,
    adjacency = adjacency,
    tom = tom,
    gene_tree = gene_tree,
    dynamic_colors = dynamic_colors,
    merged_colors = merged_colors,
    me_list = me_list,
    merged_mes = merged_mes,
    module_trait_cor = if (exists("module_trait_cor")) module_trait_cor else NULL,
    module_trait_pvalue = if (exists("module_trait_pvalue")) module_trait_pvalue else NULL,
    trait_matrix = if (exists("trait_matrix")) trait_matrix else NULL
  )
  
  # Save results
  saveRDS(metabolome_WGCNA_results, "results/networks/WGCNA/metabolome_WGCNA_results.rds")
  
  # Save module assignments
  module_assignments <- data.frame(
    Feature = colnames(metabolome_WGCNA_data$data),
    Dynamic_Module = dynamic_colors,
    Merged_Module = merged_colors,
    stringsAsFactors = FALSE
  )
  
  # Add metabolite annotations
  if (!is.null(metabolite_annotations)) {
    module_assignments <- module_assignments %>%
      left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature")
  }
  
  write.csv(module_assignments, "results/networks/WGCNA/metabolome_module_assignments.csv", row.names = FALSE)
  
  print("Metabolome WGCNA analysis completed successfully!")
} else {
  print("Not enough features after filtering. Skipping metabolome WGCNA analysis.")
  metabolome_WGCNA_results <- NULL
}

# ------------------------------------------------------------------------------
# 5. Perform WGCNA analysis for combined microbiome-metabolome data
# ------------------------------------------------------------------------------
cat("\n=== Performing WGCNA analysis for combined data ===\n")

# Check if we have enough features from both datasets
if (!is.null(microbiome_WGCNA_results) && !is.null(metabolome_WGCNA_results) &&
    ncol(microbiome_WGCNA_data$data) >= 10 && ncol(metabolome_WGCNA_data$data) >= 10) {
  
  # Combine data
  combined_data <- cbind(
    microbiome_WGCNA_data$data,
    metabolome_WGCNA_data$data
  )
  
  # Check for sample overlap
  common_samples <- intersect(
    rownames(microbiome_WGCNA_data$data),
    rownames(metabolome_WGCNA_data$data)
  )
  
  if (length(common_samples) < 3) {
    print("Not enough common samples for combined WGCNA analysis.")
  } else {
    # Subset to common samples
    combined_data <- combined_data[common_samples, ]
    
    # Prepare metadata
    combined_metadata <- metadata[common_samples, ]
    
    # Step 1: Choose soft-thresholding power
    cat("\nStep 1: Choosing soft-thresholding power...\n")
    
    # Set power range
    power_range <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
    
    # Calculate scale-free topology fit indices
    sft <- pickSoftThreshold(
      combined_data,
      powerVector = power_range,
      networkType = "unsigned",
      verbose = 5
    )
    
    # Choose soft-thresholding power
    soft_power <- sft$powerEstimate
    
    if (is.na(soft_power)) {
      soft_power <- 6  # Default value if no power is estimated
      print(paste("No soft-thresholding power estimated. Using default value:", soft_power))
    } else {
      print(paste("Chosen soft-thresholding power:", soft_power))
    }
    
    # Step 2: Construct co-expression network
    cat("\nStep 2: Constructing co-expression network...\n")
    
    # Calculate adjacency matrix
    adjacency <- adjacency(
      combined_data,
      power = soft_power,
      type = "unsigned"
    )
    
    # Convert adjacency to topological overlap matrix (TOM)
    tom <- TOMsimilarity(adjacency)
    diss_tom <- 1 - tom
    
    # Step 3: Identify modules
    cat("\nStep 3: Identifying modules...\n")
    
    # Perform hierarchical clustering
    gene_tree <- flashClust(as.dist(diss_tom), method = "average")
    
    # Identify modules using dynamic tree cut
    min_module_size <- 30
    dynamic_mods <- cutreeDynamic(
      dendro = gene_tree,
      distM = diss_tom,
      deepSplit = 2,
      pamRespectsDendro = FALSE,
      minClusterSize = min_module_size
    )
    
    # Convert numeric labels to colors
    dynamic_colors <- labels2colors(dynamic_mods)
    
    # Step 4: Merge similar modules
    cat("\nStep 4: Merging similar modules...\n")
    
    # Calculate module eigengenes
    me_list <- moduleEigengenes(combined_data, colors = dynamic_colors)
    mes <- me_list$eigengenes
    
    # Calculate dissimilarity of module eigengenes
    me_diss <- 1 - cor(mes)
    
    # Cluster module eigengenes
    me_tree <- flashClust(as.dist(me_diss), method = "average")
    
    # Merge modules
    merge_threshold <- 0.25
    merged <- mergeCloseModules(
      combined_data,
      dynamic_colors,
      cutHeight = merge_threshold,
      verbose = 3
    )
    
    # Get merged module colors and eigengenes
    merged_colors <- merged$colors
    merged_mes <- merged$newMEs
    
    # Step 5: Relate modules to traits
    cat("\nStep 5: Relating modules to traits...\n")
    
    # Create trait matrix
    if ("Group" %in% colnames(combined_metadata)) {
      # Convert group to numeric (0 for control, 1 for case)
      trait <- ifelse(combined_metadata$Group == "Control", 0, 1)
      names(trait) <- rownames(combined_metadata)
      
      # Create trait matrix
      trait_matrix <- matrix(trait, ncol = 1)
      rownames(trait_matrix) <- rownames(combined_metadata)
      colnames(trait_matrix) <- "Disease"
      
      # Calculate correlation between module eigengenes and traits
      module_trait_cor <- cor(merged_mes, trait_matrix, use = "p")
      module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(combined_data))
    }
    
    # Step 6: Analyze module composition
    cat("\nStep 6: Analyzing module composition...\n")
    
    # Create module assignments
    module_assignments <- data.frame(
      Feature = colnames(combined_data),
      Type = ifelse(colnames(combined_data) %in% colnames(microbiome_WGCNA_data$data), "Microbiome", "Metabolome"),
      Dynamic_Module = dynamic_colors,
      Merged_Module = merged_colors,
      stringsAsFactors = FALSE
    )
    
    # Add annotations
    taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
    module_assignments <- module_assignments %>%
      left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature")
    
    if (!is.null(metabolite_annotations)) {
      module_assignments <- module_assignments %>%
        left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature")
    }
    
    # Calculate module composition
    module_composition <- module_assignments %>%
      group_by(Merged_Module, Type) %>%
      summarise(Count = n()) %>%
      ungroup() %>%
      pivot_wider(names_from = Type, values_from = Count, values_fill = 0)
    
    print("Module composition:")
    print(module_composition)
    
    # Step 7: Save results
    cat("\nStep 7: Saving results...\n")
    
    # Create results list
    combined_WGCNA_results <- list(
      soft_power = soft_power,
      adjacency = adjacency,
      tom = tom,
      gene_tree = gene_tree,
      dynamic_colors = dynamic_colors,
      merged_colors = merged_colors,
      me_list = me_list,
      merged_mes = merged_mes,
      module_trait_cor = if (exists("module_trait_cor")) module_trait_cor else NULL,
      module_trait_pvalue = if (exists("module_trait_pvalue")) module_trait_pvalue else NULL,
      trait_matrix = if (exists("trait_matrix")) trait_matrix else NULL,
      module_composition = module_composition
    )
    
    # Save results
    saveRDS(combined_WGCNA_results, "results/networks/WGCNA/combined_WGCNA_results.rds")
    
    # Save module assignments
    write.csv(module_assignments, "results/networks/WGCNA/combined_module_assignments.csv", row.names = FALSE)
    
    # Save module composition
    write.csv(module_composition, "results/networks/WGCNA/combined_module_composition.csv", row.names = FALSE)
    
    print("Combined WGCNA analysis completed successfully!")
  }
} else {
  print("Not enough data for combined WGCNA analysis. Skipping.")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== WGCNA_analysis.R Summary ===\n")
cat("\nGenerated WGCNA results:\n")

if (!is.null(microbiome_WGCNA_results)) {
  cat("1. Microbiome WGCNA analysis\n")
}

if (!is.null(metabolome_WGCNA_results)) {
  cat("2. Metabolome WGCNA analysis\n")
}

if (exists("combined_WGCNA_results")) {
  cat("3. Combined microbiome-metabolome WGCNA analysis\n")
}

cat("\nWGCNA results saved to: results/networks/WGCNA/\n")
cat("WGCNA figures saved to: results/figures/networks/WGCNA/\n")
