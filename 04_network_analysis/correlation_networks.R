# correlation_networks.R
# Construct correlation networks between microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(corrplot)
library(igraph)
library(qgraph)
library(WGCNA)
library(ggraph)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/networks", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/networks", showWarnings = FALSE, recursive = TRUE)

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
# 2. Prepare data for correlation analysis
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(sample_names(phyloseq_css), colnames(metabolome_data))

if (length(common_samples) < 3) {
  stop("Not enough common samples for correlation analysis.")
}

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples:", length(common_samples)))

# ------------------------------------------------------------------------------
# 3. Function to create correlation network
# ------------------------------------------------------------------------------
create_correlation_network <- function(data1, data2 = NULL, method = "spearman", threshold = 0.6, p_value = 0.05) {
  
  # Combine data if two datasets are provided
  if (!is.null(data2)) {
    # Ensure samples are in the same order
    common_samples <- intersect(colnames(data1), colnames(data2))
    data1 <- data1[, common_samples]
    data2 <- data2[, common_samples]
    
    # Combine data
    combined_data <- rbind(data1, data2)
    rownames(combined_data) <- c(rownames(data1), rownames(data2))
  } else {
    combined_data <- data1
  }
  
  # Calculate correlation matrix
  correlation_matrix <- cor(t(combined_data), method = method)
  
  # Calculate p-values for correlations coefficients
  p_matrix <- matrix(NA, nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
  rownames(p_matrix) <- rownames(correlation_matrix)
  colnames(p_matrix) <- colnames(correlation_matrix)
  
  for (i in 1:nrow(combined_data)) {
    for (j in i:ncol(combined_data)) {
      if (i == j) {
        p_matrix[i, j] <- 0
      } else {
        test_result <- cor.test(combined_data[i, ], combined_data[j, ], method = method)
        p_matrix[i, j] <- test_result$p.value
        p_matrix[j, i] <- test_result$p.value
      }
    }
  }
  
  # Apply p-value correction
  p_matrix[lower.tri(p_matrix)] <- NA
  p_values <- p_matrix[!is.na(p_matrix)]
  p_values_adjusted <- p.adjust(p_values, method = "fdr")
  p_matrix[!is.na(p_matrix)] <- p_values_adjusted
  
  # Apply threshold to correlation matrix
  correlation_matrix[abs(correlation_matrix) < threshold | p_matrix > p_value] <- 0
  
  # Remove rows and columns with no connections
  correlation_matrix <- correlation_matrix[rowSums(abs(correlation_matrix)) > 0, colSums(abs(correlation_matrix)) > 0]
  
  return(correlation_matrix)
}

# ------------------------------------------------------------------------------
# 4. Create microbiome correlation network
# ------------------------------------------------------------------------------
cat("\n=== Creating microbiome correlation network ===\n")

# Get microbiome data
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Filter low abundance taxa (optional but recommended for network analysis)
min_abundance <- 0.001  # 0.1% relative abundance
min_prevalence <- 0.2   # Present in at least 20% of samples

# Calculate relative abundance
microbiome_rel <- t(t(microbiome_data) / colSums(microbiome_data))

# Filter taxa
taxa_to_keep <- rowSums(microbiome_rel > min_abundance) >= (min_prevalence * ncol(microbiome_rel))
microbiome_filtered <- microbiome_data[taxa_to_keep, ]

print(paste("Microbiome taxa after filtering:", sum(taxa_to_keep)))
print(paste("Original microbiome taxa:", nrow(microbiome_data)))

# Create microbiome correlation network
if (sum(taxa_to_keep) >= 2) {
  microbiome_correlation <- create_correlation_network(
    data1 = microbiome_filtered,
    method = "spearman",
    threshold = 0.6,
    p_value = 0.05
  )
  
  print(paste("Microbiome correlation matrix dimensions:", nrow(microbiome_correlation), "x", ncol(microbiome_correlation)))
  
  # Save correlation matrix
  write.csv(microbiome_correlation, "results/networks/microbiome_correlation_matrix.csv", row.names = TRUE)
  
  print("Microbiome correlation network created successfully!")
} else {
  print("Not enough microbiome taxa after filtering. Skipping microbiome correlation network creation.")
  microbiome_correlation <- NULL
}

# ------------------------------------------------------------------------------
# 5. Create metabolome correlation network
# ------------------------------------------------------------------------------
cat("\n=== Creating metabolome correlation network ===\n")

# Filter low variance metabolites (optional but recommended for network analysis)
min_variance <- quantile(apply(metabolome_data, 1, var), 0.25)  # Keep top 75% most variable metabolites
metabolite_variance <- apply(metabolome_data, 1, var)
metabolites_to_keep <- metabolite_variance >= min_variance
metabolome_filtered <- metabolome_data[metabolites_to_keep, ]

print(paste("Metabolites after filtering:", sum(metabolites_to_keep)))
print(paste("Original metabolites:", nrow(metabolome_data)))

# Create metabolome correlation network
if (sum(metabolites_to_keep) >= 2) {
  metabolome_correlation <- create_correlation_network(
    data1 = metabolome_filtered,
    method = "spearman",
    threshold = 0.6,
    p_value = 0.05
  )
  
  print(paste("Metabolome correlation matrix dimensions:", nrow(metabolome_correlation), "x", ncol(metabolome_correlation)))
  
  # Save correlation matrix
  write.csv(metabolome_correlation, "results/networks/metabolome_correlation_matrix.csv", row.names = TRUE)
  
  print("Metabolome correlation network created successfully!")
} else {
  print("Not enough metabolites after filtering. Skipping metabolome correlation network creation.")
  metabolome_correlation <- NULL
}

# ------------------------------------------------------------------------------
# 6. Create microbiome-metabolome correlation network
# ------------------------------------------------------------------------------
cat("\n=== Creating microbiome-metabolome correlation network ===\n")

# Use filtered data from previous steps
if (sum(taxa_to_keep) >= 2 && sum(metabolites_to_keep) >= 2) {
  # Create combined data matrix
  combined_data <- rbind(microbiome_filtered, metabolome_filtered)
  
  # Create correlation network
  microbe_metabolite_correlation <- create_correlation_network(
    data1 = microbiome_filtered,
    data2 = metabolome_filtered,
    method = "spearman",
    threshold = 0.5,
    p_value = 0.05
  )
  
  print(paste("Microbiome-metabolome correlation matrix dimensions:", nrow(microbe_metabolite_correlation), "x", ncol(microbe_metabolite_correlation)))
  
  # Save correlation matrix
  write.csv(microbe_metabolite_correlation, "results/networks/microbe_metabolite_correlation_matrix.csv", row.names = TRUE)
  
  print("Microbiome-metabolome correlation network created successfully!")
} else {
  print("Not enough data after filtering. Skipping microbiome-metabolome correlation network creation.")
  microbe_metabolite_correlation <- NULL
}

# ------------------------------------------------------------------------------
# 7. Create network visualization
# ------------------------------------------------------------------------------
cat("\n=== Creating network visualization ===\n")

# Function to create network visualization
create_network_visualization <- function(correlation_matrix, 
                                         node_type = NULL, 
                                         title = "Correlation Network",
                                         layout = "fr",
                                         edge_threshold = 0.1) {
  
  # Create adjacency matrix (remove diagonal)
  adjacency_matrix <- correlation_matrix
  diag(adjacency_matrix) <- 0
  
  # Remove edges below threshold
  adjacency_matrix[abs(adjacency_matrix) < edge_threshold] <- 0
  
  # Create graph object
  graph <- graph.adjacency(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Remove isolated nodes
  graph <- delete.vertices(graph, degree(graph) == 0)
  
  if (vcount(graph) == 0) {
    print("No edges left after applying threshold. Skipping network visualization.")
    return(NULL)
  }
  
  # Create node attributes
  V(graph)$name <- V(graph)$name
  
  # Set node color based on type if provided
  if (!is.null(node_type)) {
    node_colors <- ifelse(V(graph)$name %in% node_type$microbe, "red", "blue")
    V(graph)$color <- node_colors
  } else {
    V(graph)$color <- "lightblue"
  }
  
  # Set edge color based on correlation sign
  E(graph)$color <- ifelse(E(graph)$weight > 0, "green", "red")
  
  # Set edge width based on correlation strength
  E(graph)$width <- abs(E(graph)$weight) * 3
  
  # Create network visualization
  p <- ggraph(graph, layout = layout) +
    geom_edge_link(aes(color = factor(sign(weight))), width = abs(E(graph)$weight) * 2) +
    geom_node_point(aes(color = factor(V(graph)$color)), size = 3) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE) +
    labs(
      title = title,
      edge_color = "Correlation"
    ) +
    scale_edge_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
    scale_color_manual(values = if (!is.null(node_type)) c("blue", "red") else "lightblue", 
                       labels = if (!is.null(node_type)) c("Metabolite", "Microbe") else NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    guides(
      color = guide_legend(title = "Node Type"),
      edge_color = guide_legend(title = "Correlation")
    )
  
  return(p)
}

# 7.1 Microbiome network visualization
if (!is.null(microbiome_correlation) && nrow(microbiome_correlation) >= 2) {
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  # Create network visualization
  p_microbiome_network <- create_network_visualization(
    correlation_matrix = microbiome_correlation,
    title = "Microbiome Correlation Network",
    layout = "fr",
    edge_threshold = 0.1
  )
  
  if (!is.null(p_microbiome_network)) {
    # Save plot
    ggsave("results/figures/networks/microbiome_correlation_network.png", p_microbiome_network,
           width = 14, height = 12, dpi = 300)
    
    print("Microbiome network visualization generated successfully!")
  }
}

# 7.2 Metabolome network visualization
if (!is.null(metabolome_correlation) && nrow(metabolome_correlation) >= 2) {
  # Create network visualization
  p_metabolome_network <- create_network_visualization(
    correlation_matrix = metabolome_correlation,
    title = "Metabolome Correlation Network",
    layout = "fr",
    edge_threshold = 0.1
  )
  
  if (!is.null(p_metabolome_network)) {
    # Save plot
    ggsave("results/figures/networks/metabolome_correlation_network.png", p_metabolome_network,
           width = 14, height = 12, dpi = 300)
    
    print("Metabolome network visualization generated successfully!")
  }
}

# 7.3 Microbiome-metabolome network visualization
if (!is.null(microbe_metabolite_correlation) && nrow(microbe_metabolite_correlation) >= 2) {
  # Create node type information
  node_type <- list(
    microbe = rownames(microbiome_filtered),
    metabolite = rownames(metabolome_filtered)
  )
  
  # Create network visualization
  p_microbe_metabolite_network <- create_network_visualization(
    correlation_matrix = microbe_metabolite_correlation,
    node_type = node_type,
    title = "Microbiome-Metabolome Correlation Network",
    layout = "fr",
    edge_threshold = 0.1
  )
  
  if (!is.null(p_microbe_metabolite_network)) {
    # Save plot
    ggsave("results/figures/networks/microbe_metabolite_correlation_network.png", p_microbe_metabolite_network,
           width = 14, height = 12, dpi = 300)
    
    print("Microbiome-metabolome network visualization generated successfully!")
  }
}

# ------------------------------------------------------------------------------
# 8. Network analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing network analysis ===\n")

# Function to perform network analysis
perform_network_analysis <- function(correlation_matrix, title = "Network Analysis") {
  
  # Create adjacency matrix (remove diagonal)
  adjacency_matrix <- correlation_matrix
  diag(adjacency_matrix) <- 0
  
  # Create graph object
  graph <- graph.adjacency(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Remove isolated nodes
  graph <- delete.vertices(graph, degree(graph) == 0)
  
  if (vcount(graph) == 0) {
    print("No edges left. Skipping network analysis.")
    return(NULL)
  }
  
  # Calculate network metrics
  network_metrics <- list(
    Number_of_Nodes = vcount(graph),
    Number_of_Edges = ecount(graph),
    Density = edge_density(graph),
    Average_Degree = mean(degree(graph)),
    Average_Clustering_Coefficient = transitivity(graph, type = "global"),
    Average_Path_Length = average.path.length(graph),
    Modularity = modularity(cluster_louvain(graph))
  )
  
  # Print network metrics
  cat(paste("\n=== ", title, " ===\n", sep = ""))
  for (metric in names(network_metrics)) {
    cat(paste(metric, ": ", network_metrics[[metric]], "\n", sep = ""))
  }
  
  # Identify hub nodes (top 5 nodes by degree)
  node_degree <- degree(graph)
  hub_nodes <- names(sort(node_degree, decreasing = TRUE))[1:min(5, length(node_degree))]
  
  cat("\nTop 5 hub nodes:\n")
  for (node in hub_nodes) {
    cat(paste(node, ": Degree = ", node_degree[node], "\n", sep = ""))
  }
  
  return(list(
    graph = graph,
    metrics = network_metrics,
    hub_nodes = hub_nodes
  ))
}

# 8.1 Microbiome network analysis
if (!is.null(microbiome_correlation) && nrow(microbiome_correlation) >= 2) {
  microbiome_network_analysis <- perform_network_analysis(
    correlation_matrix = microbiome_correlation,
    title = "Microbiome Network Analysis"
  )
  
  if (!is.null(microbiome_network_analysis)) {
    # Save network analysis results
    saveRDS(microbiome_network_analysis, "results/networks/microbiome_network_analysis.rds")
  }
}

# 8.2 Metabolome network analysis
if (!is.null(metabolome_correlation) && nrow(metabolome_correlation) >= 2) {
  metabolome_network_analysis <- perform_network_analysis(
    correlation_matrix = metabolome_correlation,
    title = "Metabolome Network Analysis"
  )
  
  if (!is.null(metabolome_network_analysis)) {
    # Save network analysis results
    saveRDS(metabolome_network_analysis, "results/networks/metabolome_network_analysis.rds")
  }
}

# 8.3 Microbiome-metabolome network analysis
if (!is.null(microbe_metabolite_correlation) && nrow(microbe_metabolite_correlation) >= 2) {
  microbe_metabolite_network_analysis <- perform_network_analysis(
    correlation_matrix = microbe_metabolite_correlation,
    title = "Microbiome-Metabolome Network Analysis"
  )
  
  if (!is.null(microbe_metabolite_network_analysis)) {
    # Save network analysis results
    saveRDS(microbe_metabolite_network_analysis, "results/networks/microbe_metabolite_network_analysis.rds")
  }
}

# ------------------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------------------
cat("\n=== correlation_networks.R Summary ===\n")
cat("\nGenerated networks:\n")

if (!is.null(microbiome_correlation)) {
  cat("1. Microbiome correlation network\n")
}

if (!is.null(metabolome_correlation)) {
  cat("2. Metabolome correlation network\n")
}

if (!is.null(microbe_metabolite_correlation)) {
  cat("3. Microbiome-metabolome correlation network\n")
}

cat("\nNetwork visualizations saved to: results/figures/networks/\n")
cat("Network data saved to: results/networks/\n")
