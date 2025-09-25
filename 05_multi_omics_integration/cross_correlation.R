# cross_correlation.R
# Perform cross-correlation analysis between microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(corrplot)
library(Hmisc)
library(reshape2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/multi_omics_integration/cross_correlation", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/multi_omics_integration/cross_correlation", showWarnings = FALSE, recursive = TRUE)

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
# 2. Prepare data for cross-correlation analysis
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(sample_names(phyloseq_css), colnames(metabolome_data))

if (length(common_samples) < 3) {
  stop("Not enough common samples for cross-correlation analysis.")
}

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples:", length(common_samples)))

# Get microbiome data matrix
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Filter low abundance/ variance features to reduce dimensionality
# Microbiome filtering
min_abundance <- 0.001  # 0.1% relative abundance
min_prevalence <- 0.2   # Present in at least 20% of samples

# Calculate relative abundance
microbiome_rel <- t(t(microbiome_data) / colSums(microbiome_data))

# Filter taxa
taxa_to_keep <- rowSums(microbiome_rel > min_abundance) >= (min_prevalence * ncol(microbiome_rel))
microbiome_filtered <- microbiome_data[taxa_to_keep, ]

print(paste("Microbiome taxa after filtering:", sum(taxa_to_keep)))
print(paste("Original microbiome taxa:", nrow(microbiome_data)))

# Metabolome filtering
min_variance <- quantile(apply(metabolome_data, 1, var), 0.25)  # Keep top 75% most variable metabolites
metabolite_variance <- apply(metabolome_data, 1, var)
metabolites_to_keep <- metabolite_variance >= min_variance
metabolome_filtered <- metabolome_data[metabolites_to_keep, ]

print(paste("Metabolites after filtering:", sum(metabolites_to_keep)))
print(paste("Original metabolites:", nrow(metabolome_data)))

# Transpose data so that samples are rows and features are columns
microbiome_filtered_t <- t(microbiome_filtered)
metabolome_filtered_t <- t(metabolome_filtered)

print("Data prepared for cross-correlation analysis!")

# ------------------------------------------------------------------------------
# 3. Perform cross-correlation analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing cross-correlation analysis ===\n")

# Combine data matrices
combined_data <- cbind(microbiome_filtered_t, metabolome_filtered_t)

# Calculate correlation matrix using Spearman correlation
# This can be computationally intensive if there are many features
if (ncol(combined_data) <= 500) {
  # Calculate correlation matrix
  correlation_matrix <- cor(combined_data, method = "spearman")
  
  # Calculate p-values for correlation coefficients
  p_matrix <- matrix(NA, nrow = ncol(combined_data), ncol = ncol(combined_data))
  rownames(p_matrix) <- colnames(combined_data)
  colnames(p_matrix) <- colnames(combined_data)
  
  # Use Hmisc to calculate p-values
  rcorr_result <- rcorr(as.matrix(combined_data), type = "spearman")
  correlation_matrix <- rcorr_result$r
  p_matrix <- rcorr_result$P
  
  # Apply FDR correction
  p_values <- p_matrix[lower.tri(p_matrix)]
  p_values_adjusted <- p.adjust(p_values, method = "fdr")
  p_matrix[lower.tri(p_matrix)] <- p_values_adjusted
  p_matrix[upper.tri(p_matrix)] <- t(p_matrix)[upper.tri(p_matrix)]  # Make matrix symmetric
  
  # Save full correlation matrix and p-values
  write.csv(correlation_matrix, "results/multi_omics_integration/cross_correlation/full_correlation_matrix.csv", row.names = TRUE)
  write.csv(p_matrix, "results/multi_omics_integration/cross_correlation/full_p_value_matrix.csv", row.names = TRUE)
  
  print("Full correlation matrix calculated successfully!")
} else {
  print("Too many features for full correlation matrix. Using a different approach...")
  correlation_matrix <- NULL
  p_matrix <- NULL
}

# ------------------------------------------------------------------------------
# 4. Perform pairwise cross-correlation analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing pairwise cross-correlation analysis ===\n")

# Function to calculate pairwise correlations between two datasets
calculate_pairwise_correlations <- function(data1, data2, method = "spearman", p_adjust_method = "fdr") {
  
  # Initialize results data frame
  results <- data.frame(
    Feature1 = character(),
    Feature2 = character(),
    Correlation = numeric(),
    P_Value = numeric(),
    Adjusted_P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate correlations between each feature in data1 and each feature in data2
  for (i in 1:ncol(data1)) {
    for (j in 1:ncol(data2)) {
      # Get feature names
      feature1 <- colnames(data1)[i]
      feature2 <- colnames(data2)[j]
      
      # Calculate correlation
      correlation_test <- cor.test(data1[, i], data2[, j], method = method)
      
      # Add to results
      results <- rbind(results, data.frame(
        Feature1 = feature1,
        Feature2 = feature2,
        Correlation = correlation_test$estimate,
        P_Value = correlation_test$p.value,
        Adjusted_P_Value = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply multiple testing correction
  results$Adjusted_P_Value <- p.adjust(results$P_Value, method = p_adjust_method)
  
  # Order by absolute correlation
  results <- results %>% arrange(desc(abs(Correlation)))
  
  return(results)
}

# Calculate pairwise correlations between microbiome and metabolome
pairwise_correlations <- calculate_pairwise_correlations(
  data1 = microbiome_filtered_t,
  data2 = metabolome_filtered_t,
  method = "spearman",
  p_adjust_method = "fdr"
)

# Save pairwise correlations
write.csv(pairwise_correlations, "results/multi_omics_integration/cross_correlation/pairwise_correlations.csv", row.names = FALSE)

print("Pairwise cross-correlation analysis completed successfully!")
print(paste("Total pairwise correlations calculated:", nrow(pairwise_correlations)))

# ------------------------------------------------------------------------------
# 5. Identify significant correlations
# ------------------------------------------------------------------------------
cat("\n=== Identifying significant correlations ===\n")

# Set significance thresholds
correlation_threshold <- 0.6
p_value_threshold <- 0.05

# Identify significant correlations
significant_correlations <- pairwise_correlations %>%
  filter(abs(Correlation) >= correlation_threshold & Adjusted_P_Value < p_value_threshold)

print(paste("Significant correlations (|r| >=", correlation_threshold, "and FDR <", p_value_threshold, "):", nrow(significant_correlations)))

if (nrow(significant_correlations) > 0) {
  # Save significant correlation
  write.csv(significant_correlations, "results/multi_omics_integration/cross_correlation/significant_correlations.csv", row.names = FALSE)
  
  print("Top 10 significant correlation:")
  print(head(significant_correlations, 10))
} else {
  print("No significant correlation found with current thresholds.")
}

# ------------------------------------------------------------------------------
# 6. Visualize cross-correlation results
# ------------------------------------------------------------------------------
cat("\n=== Visualizing cross-correlation results ===\n")

# 6.1 Heatmap of significant correlation
if (nrow(significant_correlations) > 0) {
  # Get unique features
  unique_microbiome_features <- unique(significant_correlations$Feature1)
  unique_metabolite_features <- unique(significant_correlations$Feature2)
  
  # Create correlation matrix
  correlation_matrix_subset <- matrix(0, 
                                     nrow = length(unique_microbiome_features), 
                                     ncol = length(unique_metabolite_features),
                                     dimnames = list(unique_microbiome_features, unique_metabolite_features))
  
  # Fill correlation matrix
  for (i in 1:nrow(significant_correlations)) {
    correlation_matrix_subset[significant_correlations$Feature1[i], significant_correlations$Feature2[i]] <- 
      significant_correlations$Correlation[i]
  }
  
  # Add taxonomic information to row names
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  row_labels <- apply(taxonomy[rownames(correlation_matrix_subset), c("Phylum", "Genus")], 1, function(x) {
    paste(x[!is.na(x)], collapse = "|")
  })
  
  rownames(correlation_matrix_subset) <- row_labels
  
  # Add metabolite names to column names if available
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    col_labels <- sapply(colnames(correlation_matrix_subset), function(x) {
      if (!is.na(metabolite_annotations[x, "Metabolite.Name"]) && metabolite_annotations[x, "Metabolite.Name"] != "") {
        metabolite_annotations[x, "Metabolite.Name"]
      } else {
        x
      }
    })
    
    colnames(correlation_matrix_subset) <- col_labels
  }
  
  # Create heatmap
  png("results/figures/multi_omics_integration/cross_correlation/significant_correlations_heatmap.png", width = 1200, height = 1000, res = 300)
  pheatmap(
    correlation_matrix_subset,
    scale = "none",
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Significant Microbiome-Metabolome Correlations",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-1, 1, length.out = 101)
  )
  dev.off()
  
  print("Significant correlation heatmap generated successfully!")
}

# 6.2 Volcano plot of correlation
if (nrow(pairwise_correlations) > 0) {
  # Create volcano plot data
  volcano_data <- pairwise_correlations %>%
    mutate(
      Neg_Log10_P = -log10(Adjusted_P_Value),
      Significant = (abs(Correlation) >= correlation_threshold) & (Adjusted_P_Value < p_value_threshold)
    )
  
  # Create volcano plot
  p_volcano <- ggplot(volcano_data, aes(x = Correlation, y = Neg_Log10_P, color = Significant)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c("black", "red")) +
    geom_vline(xintercept = c(-correlation_threshold, correlation_threshold), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "gray50") +
    labs(
      title = "Microbiome-Metabolome Correlation Volcano Plot",
      x = "Spearman Correlation Coefficient",
      y = "-Log10(Adjusted P-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add count annotations
  significant_count <- sum(volcano_data$Significant, na.rm = TRUE)
  positive_count <- sum(volcano_data$Correlation > correlation_threshold & volcano_data$Adjusted_P_Value < p_value_threshold, na.rm = TRUE)
  negative_count <- sum(volcano_data$Correlation < -correlation_threshold & volcano_data$Adjusted_P_Value < p_value_threshold, na.rm = TRUE)
  
  p_volcano <- p_volcano +
    annotate("text", x = max(volcano_data$Correlation, na.rm = TRUE) * 0.9, y = max(volcano_data$Neg_Log10_P, na.rm = TRUE) * 0.9,
             label = paste("Total Significant:", significant_count, "\nPositive:", positive_count, "\nNegative:", negative_count),
             hjust = 1, vjust = 1, size = 4,
             bbox = list(boxstyle = "round,pad=0.3", fill = "white", alpha = 0.8))
  
  # Save plot
  ggsave("results/figures/multi_omics_integration/cross_correlation/correlation_volcano.png", p_volcano,
         width = 12, height = 10, dpi = 300)
  
  print("Correlation volcano plot generated successfully!")
}

# 6.3 Network visualization of significant correlation
if (nrow(significant_correlations) > 0) {
  # Create network data
  network_data <- significant_correlations %>%
    select(Feature1, Feature2, Correlation) %>%
    rename(Source = Feature1, Target = Feature2, Weight = Correlation)
  
  # Create graph object
  graph <- graph_from_data_frame(network_data, directed = FALSE)
  
  # Set node attributes
  V(graph)$type <- ifelse(V(graph)$name %in% colnames(microbiome_filtered_t), "Microbiome", "Metabolite")
  
  # Add taxonomic information to microbiome nodes
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  for (node in V(graph)$name) {
    if (V(graph)[node]$type == "Microbiome" && node %in% rownames(taxonomy)) {
      tax_info <- taxonomy[node, ]
      if (!is.na(tax_info$Genus)) {
        V(graph)[node]$label <- paste(tax_info$Phylum, tax_info$Genus, sep = "|")
      } else if (!is.na(tax_info$Family)) {
        V(graph)[node]$label <- paste(tax_info$Phylum, tax_info$Family, sep = "|")
      } else {
        V(graph)[node]$label <- tax_info$Phylum
      }
    } else if (V(graph)[node]$type == "Metabolite" && !is.null(metabolite_annotations) && node %in% rownames(metabolite_annotations)) {
      if (!is.na(metabolite_annotations[node, "Metabolite.Name"]) && metabolite_annotations[node, "Metabolite.Name"] != "") {
        V(graph)[node]$label <- metabolite_annotations[node, "Metabolite.Name"]
      } else {
        V(graph)[node]$label <- node
      }
    } else {
      V(graph)[node]$label <- node
    }
  }
  
  # Set edge attributes
  E(graph)$color <- ifelse(E(graph)$Weight > 0, "green", "red")
  E(graph)$width <- abs(E(graph)$Weight) * 3
  
  # Create network visualization
  p_network <- ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = factor(sign(Weight))), width = abs(E(graph)$Weight) * 2) +
    geom_node_point(aes(color = type, shape = type), size = 3) +
    geom_node_text(aes(label = label), size = 2, repel = TRUE) +
    labs(
      title = "Microbiome-Metabolome Correlation Network",
      edge_color = "Correlation",
      color = "Node Type",
      shape = "Node Type"
    ) +
    scale_edge_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
    scale_color_manual(values = c("blue", "red")) +
    scale_shape_manual(values = c(21, 22)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Save plot
  ggsave("results/figures/multi_omics_integration/cross_correlation/correlation_network.png", p_network,
         width = 14, height = 12, dpi = 300)
  
  print("Correlation network visualization generated successfully!")
}

# 6.4 Lollipop plot of top correlation
if (nrow(significant_correlations) > 0) {
  # Select top N correlation
  top_n <- min(30, nrow(significant_correlations))
  top_correlations <- significant_correlations %>% head(top_n)
  
  # Create labels for pairs
  top_correlations <- top_correlations %>%
    left_join(taxonomy %>% rownames_to_column("Feature1"), by = "Feature1") %>%
    mutate(
      Microbiome_Label = case_when(
        !is.na(Genus) ~ paste(Phylum, Genus, sep = "|"),
        !is.na(Family) ~ paste(Phylum, Family, sep = "|"),
        TRUE ~ Phylum
      )
    ) %>%
    select(-Phylum, -Class, -Order, -Family, -Genus, -Species)
  
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    top_correlations <- top_correlations %>%
      left_join(metabolite_annotations %>% select(Metabolite.Name) %>% rownames_to_column("Feature2"), by = "Feature2") %>%
      mutate(
        Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Feature2)
      ) %>%
      select(-Metabolite.Name)
  } else {
    top_correlations <- top_correlations %>%
      mutate(Metabolite_Label = Feature2)
  }
  
  # Create pair label
  top_correlations <- top_correlations %>%
    mutate(Pair_Label = paste(Microbiome_Label, "vs", Metabolite_Label, sep = "\n")) %>%
    arrange(desc(abs(Correlation)))
  
  # Create lollipop plot
  p_lollipop <- ggplot(top_correlations, aes(y = reorder(Pair_Label, abs(Correlation)), x = Correlation, color = Correlation > 0)) +
    geom_segment(aes(x = 0, xend = Correlation, y = Pair_Label, yend = Pair_Label), color = "gray50") +
    geom_point(size = 3) +
    scale_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
    labs(
      title = "Top Microbiome-Metabolite Correlations",
      x = "Spearman Correlation Coefficient",
      y = "Feature Pair"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Save plot
  ggsave("results/figures/multi_omics_integration/cross_correlation/top_correlations_lollipop.png", p_lollipop,
         width = 14, height = 12, dpi = 300)
  
  print("Top correlation lollipop plot generated successfully!")
}

# ------------------------------------------------------------------------------
# 7. Functional enrichment analysis of correlated features
# ------------------------------------------------------------------------------
cat("\n=== Performing functional enrichment analysis ===\n")

# 7.1 Microbiome functional enrichment (if annotations are available)
if (nrow(significant_correlations) > 0 && "Phylum" %in% colnames(taxonomy)) {
  # Get significant microbiome features
  significant_microbiome_features <- unique(significant_correlations$Feature1)
  
  # Get their taxonomy
  microbiome_taxonomy <- taxonomy[significant_microbiome_features, ]
  
  # Count phylum distribution
  phylum_counts <- table(microbiome_taxonomy$Phylum)
  
  # Create bar plot
  p_phylum <- ggplot(data.frame(Phylum = names(phylum_counts), Count = as.vector(phylum_counts)),
                     aes(x = reorder(Phylum, Count), y = Count, fill = Phylum)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Phylum Distribution of Significant Microbiome Features",
      x = "Phylum",
      y = "Number of Features"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  ggsave("results/figures/multi_omics_integration/cross_correlation/microbiome_phylum_distribution.png", p_phylum,
         width = 12, height = 8, dpi = 300)
  
  print("Microbiome phylum distribution plot generated successfully!")
}

# 7.2 Metabolite pathway enrichment (if annotations are available)
if (nrow(significant_correlations) > 0 && !is.null(metabolite_annotations) && "Pathway" %in% colnames(metabolite_annotations)) {
  # Get significant metabolite features
  significant_metabolite_features <- unique(significant_correlations$Feature2)
  
  # Get their pathway annotations
  metabolite_pathways <- metabolite_annotations[significant_metabolite_features, "Pathway"]
  
  # Remove NA values
  metabolite_pathways <- metabolite_pathways[!is.na(metabolite_pathways)]
  
  if (length(metabolite_pathways) > 0) {
    # Count pathway distribution
    pathway_counts <- table(metabolite_pathways)
    
    # Keep only pathways with at least 2 metabolites
    pathway_counts <- pathway_counts[pathway_counts >= 2]
    
    if (length(pathway_counts) > 0) {
      # Create bar plot
      p_pathway <- ggplot(data.frame(Pathway = names(pathway_counts), Count = as.vector(pathway_counts)),
                         aes(x = reorder(Pathway, Count), y = Count, fill = Pathway)) +
        geom_bar(stat = "identity") +
        labs(
          title = "Pathway Distribution of Significant Metabolite Features",
          x = "Pathway",
          y = "Number of Features"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      
      ggsave("results/figures/multi_omics_integration/cross_correlation/metabolite_pathway_distribution.png", p_pathway,
             width = 12, height = 8, dpi = 300)
      
      print("Metabolite pathway distribution plot generated successfully!")
    }
  }
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== cross_correlation.R Summary ===\n")
cat("\nGenerated results:\n")
cat("1. Pairwise cross-correlation analysis between microbiome and metabolome\n")
if (ncol(combined_data) <= 500) {
  cat("2. Full correlation matrix and p-value matrix\n")
}
cat("3. Identification of significant correlation\n")
cat("4. Visualization of significant correlation (heatmap)\n")
cat("5. Correlation volcano plot\n")
if (nrow(significant_correlations) > 0) {
  cat("6. Correlation network visualization\n")
  cat("7. Top correlation lollipop plot\n")
  if ("Phylum" %in% colnames(taxonomy)) {
    cat("8. Microbiome phylum distribution plot\n")
  }
  if (!is.null(metabolite_annotations) && "Pathway" %in% colnames(metabolite_annotations)) {
    cat("9. Metabolite pathway distribution plot\n")
  }
}

cat("\nResults saved to: results/multi_omics_integration/cross_correlation/\n")
cat("Figures saved to: results/figures/multi_omics_integration/cross_correlation/\n")
