# graph_visualization.R
# Visualize network graphs structures from correlation and WGCNA analysis

# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tibble)
library(igraph)
library(ggraph)
library(networkD3)
library(corrplot)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)
library(WGCNA)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/figures/networks/graphs", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load correlation matrices
microbiome_correlation <- read.csv("results/networks/microbiome_correlation_matrix.csv", row.names = 1, check.names = FALSE)
metabolome_correlation <- read.csv("results/networks/metabolome_correlation_matrix.csv", row.names = 1, check.names = FALSE)
microbe_metabolite_correlation <- read.csv("results/networks/microbe_metabolite_correlation_matrix.csv", row.names = 1, check.names = FALSE)

# Load WGCNA results
microbiome_WGCNA_results <- tryCatch({
  readRDS("results/networks/WGCNA/microbiome_WGCNA_results.rds")
}, error = function(e) {
  NULL
})

metabolome_WGCNA_results <- tryCatch({
  readRDS("results/networks/WGCNA/metabolome_WGCNA_results.rds")
}, error = function(e) {
  NULL
})

combined_WGCNA_results <- tryCatch({
  readRDS("results/networks/WGCNA/combined_WGCNA_results.rds")
}, error = function(e) {
  NULL
})

# Load module assignments
microbiome_module_assignments <- read.csv("results/networks/WGCNA/microbiome_module_assignments.csv", stringsAsFactors = FALSE)
metabolome_module_assignments <- read.csv("results/networks/WGCNA/metabolome_module_assignments.csv", stringsAsFactors = FALSE)
combined_module_assignments <- read.csv("results/networks/WGCNA/combined_module_assignments.csv", stringsAsFactors = FALSE)

# Load taxonomy and annotations
taxonomy <- tax_table(readRDS("data/processed/phyloseq_object.rds")) %>% as.data.frame()
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Function to create graph from correlation matrix
# ------------------------------------------------------------------------------
create_graph_from_correlation <- function(correlation_matrix, threshold = 0.1) {
  
  # Create adjacency matrix (remove diagonal)
  adjacency_matrix <- as.matrix(correlation_matrix)
  diag(adjacency_matrix) <- 0
  
  # Remove edges below threshold
  adjacency_matrix[abs(adjacency_matrix) < threshold] <- 0
  
  # Create graph object
  graph <- graph.adjacency(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Remove isolated nodes
  graph <- delete.vertices(graph, degree(graph) == 0)
  
  if (vcount(graph) == 0) {
    print("No edges left after applying threshold.")
    return(NULL)
  }
  
  return(graph)
}

# ------------------------------------------------------------------------------
# 3. Visualize correlation networks graphs
# ------------------------------------------------------------------------------
cat("\n=== Visualizing correlation networks ===\n")

# 3.1 Microbiome correlation network
if (nrow(microbiome_correlation) >= 2) {
  # Create graph
  microbiome_graph <- create_graph_from_correlation(microbiome_correlation, threshold = 0.1)
  
  if (!is.null(microbiome_graph)) {
    # Get module information if available
    if (!is.null(microbiome_module_assignments) && "Merged_Module" %in% colnames(microbiome_module_assignments)) {
      # Create node attributes for module
      node_modules <- microbiome_module_assignments$Merged_Module
      names(node_modules) <- microbiome_module_assignments$Feature
      
      # Subset to nodes present in the graph
      node_modules <- node_modules[names(node_modules) %in% V(microbiome_graph)$name]
      
      # Add module information to graph
      V(microbiome_graph)$module <- node_modules[V(microbiome_graph)$name]
      
      # Create color palette for modules
      module_colors <- unique(node_modules)
      names(module_colors) <- module_colors
    } else {
      V(microbiome_graph)$module <- "NA"
      module_colors <- c("NA" = "lightblue")
    }
    
    # Create network visualization
    p_microbiome_network <- ggraph(microbiome_graph, layout = "fr") +
      geom_edge_link(aes(color = factor(sign(weight))), width = abs(E(microbiome_graph)$weight) * 2) +
      geom_node_point(aes(color = module), size = 3) +
      geom_node_text(aes(label = name), size = 2, repel = TRUE) +
      labs(
        title = "Microbiome Correlation Network",
        edge_color = "Correlation"
      ) +
      scale_edge_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
      scale_color_manual(values = module_colors) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      ) +
      guides(
        color = guide_legend(title = "Module"),
        edge_color = guide_legend(title = "Correlation")
      )
    
    # Save plot
    ggsave("results/figures/networks/graphs/microbiome_correlation_network.png", p_microbiome_network,
           width = 14, height = 12, dpi = 300)
    
    print("Microbiome correlation network visualization generated successfully!")
  }
}

# 3.2 Metabolome correlation network
if (nrow(metabolome_correlation) >= 2) {
  # Create graph
  metabolome_graph <- create_graph_from_correlation(metabolome_correlation, threshold = 0.1)
  
  if (!is.null(metabolome_graph)) {
    # Get module information if available
    if (!is.null(metabolome_module_assignments) && "Merged_Module" %in% colnames(metabolome_module_assignments)) {
      # Create node attribute for module
      node_modules <- metabolome_module_assignments$Merged_Module
      names(node_modules) <- metabolome_module_assignments$Feature
      
      # Subset to nodes_modules present in the graph
      node_modules <- node_modules[names(node_modules) %in% V(metabolome_graph)$name]
      
      # Add module information to graph
      V(metabolome_graph)$module <- node_modules[V(metabolome_graph)$name]
      
      # Create color palette for modules
      module_colors <- unique(node_modules)
      names(module_colors) <- module_colors
    } else {
      V(metabolome_graph)$module <- "NA"
      module_colors <- c("NA" = "lightblue")
    }
    
    # Create network visualization
    p_metabolome_network <- ggraph(metabolome_graph, layout = "fr") +
      geom_edge_link(aes(color = factor(sign(weight))), width = abs(E(metabolome_graph)$weight) * 2) +
      geom_node_point(aes(color = module), size = 3) +
      geom_node_text(aes(label = name), size = 2, repel = TRUE) +
      labs(
        title = "Metabolome Correlation Network",
        edge_color = "Correlation"
      ) +
      scale_edge_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
      scale_color_manual(values = module_colors) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      ) +
      guides(
        color = guide_legend(title = "Module"),
        edge_color = guide_legend(title = "Correlation")
      )
    
    # Save plot
    ggsave("results/figures/networks/graphs/metabolome_correlation_network.png", p_metabolome_network,
           width = 14, height = 12, dpi = 300)
    
    print("Metabolome correlation network visualization generated successfully!")
  }
}

# 3.3 Microbiome-metabolome correlation network
if (nrow(microbe_metabolite_correlation) >= 2) {
  # Create graph
  microbe_metabolite_graph <- create_graph_from_correlation(microbe_metabolite_correlation, threshold = 0.1)
  
  if (!is.null(microbe_metabolite_graph)) {
    # Determine node type (microbe or metabolite)
    microbe_nodes <- rownames(microbiome_correlation)
    metabolite_nodes <- rownames(metabolome_correlation)
    
    V(microbe_metabolite_graph)$type <- ifelse(
      V(microbe_metabolite_graph)$name %in% microbe_nodes,
      "Microbe",
      "Metabolite"
    )
    
    # Get module information if available
    if (!is.null(combined_module_assignments) && "Merged_Module" %in% colnames(combined_module_assignments)) {
      # Create node attribute for module
      node_modules <- combined_module_assignments$Merged_Module
      names(node_modules) <- combined_module_assignments$Feature
      
      # Subset to nodes present in the graph
      node_modules <- node_modules[names(node_modules) %in% V(microbe_metabolite_graph)$name]
      
      # Add module information to graph
      V(microbe_metabolite_graph)$module <- node_modules[V(microbe_metabolite_graph)$name]
      
      # Create color palette for modules
      module_colors <- unique(node_modules)
      names(module_colors) <- module_colors
    } else {
      V(microbe_metabolite_graph)$module <- "NA"
      module_colors <- c("NA" = "lightblue")
    }
    
    # Create network visualization
    p_microbe_metabolite_network <- ggraph(microbe_metabolite_graph, layout = "fr") +
      geom_edge_link(aes(color = factor(sign(weight))), width = abs(E(microbe_metabolite_graph)$weight) * 2) +
      geom_node_point(aes(color = module, shape = type), size = 3) +
      geom_node_text(aes(label = name), size = 2, repel = TRUE) +
      labs(
        title = "Microbiome-Metabolome Correlation Network",
        edge_color = "Correlation",
        shape = "Node Type"
      ) +
      scale_edge_color_manual(values = c("red", "green"), labels = c("Negative", "Positive")) +
      scale_color_manual(values = module_colors) +
      scale_shape_manual(values = c(21, 22)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      ) +
      guides(
        color = guide_legend(title = "Module"),
        edge_color = guide_legend(title = "Correlation"),
        shape = guide_legend(title = "Node Type")
      )
    
    # Save plot
    ggsave("results/figures/networks/graphs/microbe_metabolite_correlation_network.png", p_microbe_metabolite_network,
           width = 14, height = 12, dpi = 300)
    
    print("Microbiome-metabolite correlation network visualization generated successfully!")
  }
}

# ------------------------------------------------------------------------------
# 4. Visualize WGCNA network modules
# ------------------------------------------------------------------------------
cat("\n=== Visualizing WGCNA network modules ===\n")

# Function to visualize WGCNA module
visualize_WGCNA_module <- function(WGCNA_results, module_assignments, title = "WGCNA Module Visualization") {
  
  if (is.null(WGCNA_results) || is.null(module_assignments)) {
    return(NULL)
  }
  
  # Get TOM matrix
  tom <- WGCNA_results$tom
  
  if (is.null(tom)) {
    print("TOM matrix not found in WGCNA results.")
    return(NULL)
  }
  
  # Get module colors
  module_colors <- unique(module_assignments$Merged_Module)
  
  # Create plots for each module
  for (module in module_colors) {
    if (module == "grey") {
      next  # Skip grey module (unassigned features)
    }
    
    # Get features in this module
    module_features <- module_assignments$Feature[module_assignments$Merged_Module == module]
    
    # Subset TOM matrix to this module
    if (length(module_features) >= 2) {
      # Check if features are present in TOM matrix
      module_features <- intersect(module_features, rownames(tom))
      
      if (length(module_features) >= 2) {
        module_tom <- tom[module_features, module_features]
        
        # Create graph from TOM matrix
        module_graph <- graph.adjacency(module_tom, mode = "undirected", weighted = TRUE, diag = FALSE)
        
        # Remove isolated nodes
        module_graph <- delete.vertices(module_graph, degree(module_graph) == 0)
        
        if (vcount(module_graph) >= 2) {
          # Create network visualization
          p_module <- ggraph(module_graph, layout = "fr") +
            geom_edge_link(aes(width = weight), color = "grey") +
            geom_node_point(color = module, size = 3) +
            geom_node_text(aes(label = name), size = 2, repel = TRUE) +
            labs(
              title = paste(title, "-", module, "Module"),
              edge_width = "TOM Similarity"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom"
            )
          
          # Save plot
          ggsave(paste0("results/figures/networks/graphs/", tolower(title) %>% gsub(" ", "_", .), "_", module, "_module.png"),
                 p_module, width = 12, height = 10, dpi = 300)
          
          print(paste(module, "module visualization generated successfully!"))
        }
      }
    }
  }
  
  return(TRUE)
}

# 4.1 Visualize microbiome WGCNA modules
if (!is.null(microbiome_WGCNA_results) && nrow(microbiome_module_assignments) > 0) {
  visualize_WGCNA_module(
    WGCNA_results = microbiome_WGCNA_results,
    module_assignments = microbiome_module_assignments,
    title = "Microbiome WGCNA Module"
  )
}

# 4.2 Visualize metabolome WGCNA module
if (!is.null(metabolome_WGCNA_results) && nrow(metabolome_module_assignments) > 0) {
  visualize_WGCNA_module(
    WGCNA_results = metabolome_WGCNA_results,
    module_assignments = metabolome_module_assignments,
    title = "Metabolome WGCNA Module"
  )
}

# 4.3 Visualize combined WGCNA module
if (!is.null(combined_WGCNA_results) && nrow(combined_module_assignments) > 0) {
  # Create visualization for combined modules
  visualize_WGCNA_module(
    WGCNA_results = combined_WGCNA_results,
    module_assignments = combined_module_assignments,
    title = "Combined WGCNA Module"
  )
  
  # Create module composition visualization
  if ("module_composition" %in% names(combined_WGCNA_results)) {
    module_composition <- combined_WGCNA_results$module_composition
    
    if (!is.null(module_composition) && nrow(module_composition) > 0) {
      # Reshape data for plotting
      plot_data <- module_composition %>%
        pivot_longer(cols = -Merged_Module, names_to = "Type", values_to = "Count")
      
      # Create stacked barplot
      p_module_composition <- ggplot(plot_data, aes(x = Merged_Module, y = Count, fill = Type)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(
          title = "Combined Network Module Composition",
          x = "Module",
          y = "Number of Features",
          fill = "Feature Type"
        ) +
        scale_fill_manual(values = c("Microbiome" = "red", "Metabolome" = "blue")) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
        )
      
      # Save plot
      ggsave("results/figures/networks/graphs/combined_module_composition.png", p_module_composition,
             width = 12, height = 8, dpi = 300)
      
      print("Combined module composition visualization generated successfully!")
    }
  }
}

# ------------------------------------------------------------------------------
# 5. Create interactive network visualization
# ------------------------------------------------------------------------------
cat("\n=== Creating interactive network visualization ===\n")

# Function to create interactive network visualization
create_interactive_network <- function(graph, title = "Interactive Network", filename = "interactive_network.html") {
  
  if (is.null(graph) || vcount(graph) == 0) {
    return(NULL)
  }
  
  # Convert igraph to networkD3 format
  network_d3 <- igraph_to_networkD3(graph)
  
  # Add node attributes
  if ("type" %in% vertex_attr_names(graph)) {
    network_d3$nodes$type <- V(graph)$type
  }
  
  if ("module" %in% vertex_attr_names(graph)) {
    network_d3$nodes$module <- V(graph)$module
  }
  
  # Create interactive network
  p <- forceNetwork(
    Links = network_d3$links,
    Nodes = network_d3$nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    Group = if ("module" %in% colnames(network_d3$nodes)) "module" else "group",
    opacity = 0.7,
    zoom = TRUE,
    legend = TRUE,
    bounded = TRUE
  )
  
  # Save interactive network
  htmlwidgets::saveWidget(p, file = paste0("results/figures/networks/graphs/", filename), selfcontained = TRUE)
  
  return(p)
}

# 5.1 Create interactive microbiome network
if (nrow(microbiome_correlation) >= 2) {
  # Create graph
  microbiome_graph <- create_graph_from_correlation(microbiome_correlation, threshold = 0.1)
  
  if (!is.null(microbiome_graph)) {
    create_interactive_network(
      graph = microbiome_graph,
      title = "Interactive Microbiome Correlation Network",
      filename = "interactive_microbiome_network.html"
    )
    
    print("Interactive microbiome network visualization generated successfully!")
  }
}

# 5.2 Create interactive metabolome network
if (nrow(metabolome_correlation) >= 2) {
  # Create graph
  metabolome_graph <- create_graph_from_correlation(metabolome_correlation, threshold = 0.1)
  
  if (!is.null(metabolome_graph)) {
    create_interactive_network(
      graph = metabolome_graph,
      title = "Interactive Metabolome Correlation Network",
      filename = "interactive_metabolome_network.html"
    )
    
    print("Interactive metabolome network visualization generated successfully!")
  }
}

# 5.3 Create interactive microbiome-metabolome network
if (nrow(microbe_metabolite_correlation) >= 2) {
  # Create graph
  microbe_metabolite_graph <- create_graph_from_correlation(microbe_metabolite_correlation, threshold = 0.1)
  
  if (!is.null(microbe_metabolite_graph)) {
    # Determine node type (microbe or metabolite)
    microbe_nodes <- rownames(microbiome_correlation)
    metabolite_nodes <- rownames(metabolome_correlation)
    
    V(microbe_metabolite_graph)$type <- ifelse(
      V(microbe_metabolite_graph)$name %in% microbe_nodes,
      "Microbe",
      "Metabolite"
    )
    
    create_interactive_network(
      graph = microbe_metabolite_graph,
      title = "Interactive Microbiome-Metabolome Correlation Network",
      filename = "interactive_microbe_metabolite_network.html"
    )
    
    print("Interactive microbiome-metabolite network visualization generated successfully!")
  }
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== graph_visualization.R Summary ===\n")
cat("\nGenerated network visualizations:\n")

if (nrow(microbiome_correlation) >= 2) {
  cat("1. Microbiome correlation network visualization\n")
  cat("2. Interactive microbiome network visualization\n")
}

if (nrow(metabolome_correlation) >= 2) {
  cat("3. Metabolome correlation network visualization\n")
  cat("4. Interactive metabolome network visualization\n")
}

if (nrow(microbe_metabolite_correlation) >= 2) {
  cat("5. Microbiome-metabolite correlation network visualization\n")
  cat("6. Interactive microbiome-metabolite network visualization\n")
}

if (!is.null(microbiome_WGCNA_results) && nrow(microbiome_module_assignments) > 0) {
  cat("7. Microbiome WGCNA module visualizations\n")
}

if (!is.null(metabolome_WGCNA_results) && nrow(metabolome_module_assignments) > 0) {
  cat("8. Metabolome WGCNA module visualizations\n")
}

if (!is.null(combined_WGCNA_results) && nrow(combined_module_assignments) > 0) {
  cat("9. Combined WGCNA module visualizations\n")
  cat("10. Combined module composition visualization\n")
}

cat("\nNetwork visualizations saved to: results/figures/networks/graphs/\n")
