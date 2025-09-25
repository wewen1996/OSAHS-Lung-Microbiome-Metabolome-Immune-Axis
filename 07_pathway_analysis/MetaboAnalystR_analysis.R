# MetaboAnalystR_analysis.R
# Perform pathway analysis using MetaboAnalystR package

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(MetaboAnalystR)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/pathway_analysis/MetaboAnalystR", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/pathway_analysis/MetaboAnalystR", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load metabolome data
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Prepare data for MetaboAnalystR
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(colnames(metabolome_data), rownames(metadata))

if (length(common_samples) < 3) {
  stop("Not enough common samples for pathway analysis.")
}

# Subset data to common samples
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples IDs:", length(common_samples)))

# Transpose data so that samples are rows and metabolites are columns
metabolome_data_t <- t(metabolome_data) %>% as.data.frame()

# Create sample metadata
sample_metadata <- data.frame(
  Sample = rownames(metabolome_data_t),
  Group = metadata$Group
)

# Create metabolite annotations table
if (!is.null(metabolite_annotations) && "HMDB.ID" %in% colnames(metabolite_annotations)) {
  metabolite_info <- metabolite_annotations %>%
    select(HMDB.ID, Metabolite.Name, Class, Formula) %>%
    rownames_to_column("Metabolite") %>%
    rename(
      hmdb = HMDB.ID,
      name = Metabolite.Name,
      super_class = Class,
      formula = Formula
    )
  
  # Remove rows with missing HMDB IDs
  metabolite_info <- metabolite_info %>% filter(!is.na(hmdb) & hmdb != "")
  
  print(paste("Metabolites with HMDB IDs:", nrow(metabolite_info)))
} else {
  stop("HMDB IDs not found in metabolite annotations. Cannot perform pathway analysis.")
}

print("Data prepared for MetaboAnalystR!")

# ------------------------------------------------------------------------------
# 3. Perform pathway analysis using MetaboAnalystR
# ------------------------------------------------------------------------------
cat("\n=== Performing pathway analysis using MetaboAnalystR ===\n")

# Create MetaboAnalystR object
mSet <- InitDataObjects("conc", "mset", FALSE)

# Set data
mSet <- SetupData(mSet, metabolome_data_t, sample_metadata, "group", FALSE)

# Set metabolite annotations
mSet <- Setup.MapData(mSet, metabolite_info, "hmdb", "name")

# Perform pathway analysis
mSet <- PerformPathwayAnalysis(mSet, "pathora", "hsa", "fdr", 0.05)

# Save MetaboAnalystR results
saveRDS(mSet, "results/pathway_analysis/MetaboAnalystR/metaboanalyst_results.rds")

print("Pathway analysis completed successfully!")

# ------------------------------------------------------------------------------
# 4. Extract and visualize pathway analysis results
# ------------------------------------------------------------------------------
cat("\n=== Extracting and visualizing pathway analysis results ===\n")

# Extract pathway analysis results
pathway_results <- mSet$analSet$pathora$pathway.res %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  arrange(desc(Impact))

# Save pathway results
write.csv(pathway_results, "results/pathway_analysis/MetaboAnalystR/pathway_analysis_results.csv", row.names = FALSE)

print("Pathway analysis results extracted successfully!")
print(paste("Total pathways analyzed:", nrow(pathway_results)))

# Identify significant pathways
significant_pathways <- pathway_results %>%
  filter(Adjust.P.Value < 0.05) %>%
  arrange(desc(Impact))

print(paste("Significant pathways (FDR < 0.05):", nrow(significant_pathways)))

if (nrow(significant_pathways) > 0) {
  # Save significant pathways
  write.csv(significant_pathways, "results/pathway_analysis/MetaboAnalystR/significant_pathways.csv", row.names = FALSE)
  
  print("Top 10 significant pathways:")
  print(head(significant_pathways[, c("Pathway", "P.Value", "Adjust.P.Value", "Impact")], 10))
}

# ------------------------------------------------------------------------------
# 5. Visualize pathway analysis results
# ------------------------------------------------------------------------------
cat("\n=== Visualizing pathway analysis results ===\n")

# 5.1 Pathway impact plot
if (nrow(pathway_results) > 0) {
  # Create bubble plot
  p_bubble <- ggplot(pathway_results, aes(x = log10(P.Value), y = Impact, size = log10(Count), color = Adjust.P.Value)) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = log10(0.05), linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "Pathway Impact Plot",
      x = "-Log10(P-value)",
      y = "Pathway Impact",
      size = "-Log10(Metabolite Count)",
      color = "Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Add labels for significant pathways
  if (nrow(significant_pathways) > 0) {
    top_pathways <- significant_pathways %>% head(10)
    
    p_bubble <- p_bubble +
      geom_text_repel(
        data = top_pathways,
        aes(label = Pathway),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.5,
        segment.color = "gray50"
      )
  }
  
  # Save plot
  ggsave("results/figures/pathway_analysis/MetaboAnalystR/pathway_impact_plot.png", p_bubble,
         width = 12, height = 10, dpi = 300)
  
  print("Pathway impact plot generated successfully!")
}

# 5.2 Pathway enrichment bar plot
if (nrow(significant_pathways) > 0) {
  # Select top N pathways
  top_n <- min(20, nrow(significant_pathways))
  top_pathways <- significant_pathways %>% head(top_n)
  
  # Create bar plot
  p_bar <- ggplot(top_pathways, aes(x = reorder(Pathway, Impact), y = Impact, fill = Adjust.P.Value)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(direction = -1) +
    coord_flip() +
    labs(
      title = "Top Pathway Impact",
      x = "Pathway",
      y = "Impact",
      fill = "Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Save plot
  ggsave("results/figures/pathway_analysis/MetaboAnalystR/pathway_enrichment_barplot.png", p_bar,
         width = 12, height = 10, dpi = 300)
  
  print("Pathway enrichment bar plot generated successfully!")
}

# 5.3 Pathway heatmap
if (nrow(significant_pathways) > 0 && nrow(metabolite_data_t) > 0) {
  # Get significant pathway names
  significant_pathway_names <- significant_pathways$Pathway
  
  # Get pathway-metabolite association
  pathway_metabolite <- mSet$analSet$pathora$pathway.mat %>%
    as.data.frame() %>%
    rownames_to_column("Metabolite") %>%
    filter(Pathway %in% significant_pathway_names)
  
  if (nrow(pathway_metabolite) > 0) {
    # Create heatmap data
    heatmap_data <- pathway_metabolite %>%
      pivot_wider(names_from = Pathway, values_from = Value, values_fill = 0) %>%
      column_to_rownames("Metabolite")
    
    # Add metabolite names if available
    if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
      rownames(heatmap_data) <- sapply(rownames(heatmap_data), function(x) {
        if (!is.na(metabolite_annotations[x, "Metabolite.Name"]) && metabolite_annotations[x, "Metabolite.Name"] != "") {
          metabolite_annotations[x, "Metabolite.Name"]
        } else {
          x
        }
      })
    }
    
    # Create heatmap
    p_heatmap <- pheatmap(
      heatmap_data,
      scale = "row",
      show_rownames = TRUE,
      show_colnames = TRUE,
      treeheight_row = 20,
      treeheight_col = 20,
      fontsize_row = 8,
      fontsize_col = 8,
      main = "Significant Pathway Metabolite Heatmap",
      color = viridis(100),
      filename = "results/figures/pathway_analysis/MetaboAnalystR/pathway_metabolite_heatmap.png",
      width = 14,
      height = 12,
      dpi = 300
    )
    
    print("Pathway metabolite heatmap generated successfully!")
  }
}

# 5.4 Pathway network
if (nrow(significant_pathways) > 0) {
  # Try to generate pathway network if available
  tryCatch({
    # This function is specific to MetaboAnalystR and may require additional setup
    png("results/figures/pathway_analysis/MetaboAnalystR/pathway_network.png", width = 1200, height = 1000, res = 300)
    PlotPathwayNetwork(mSet, "pathora", 0.05)
    dev.off()
    
    print("Pathway network plot generated successfully!")
  }, error = function(e) {
    print("Could not generate pathway network plot. This may require additional setup or data.")
    print(e$message)
  })
}

# ------------------------------------------------------------------------------
# 6. Pathway analysis by group
# ------------------------------------------------------------------------------
cat("\n=== Performing pathway analysis by group ===\n")

# Check if there are at least two groups
if (length(unique(metadata$Group)) >= 2) {
  # For each group, perform pathway analysis
  for (group in unique(metadata$Group)) {
    cat(paste("\nAnalyzing pathway for", group, "group...\n"))
    
    # Subset data to current group
    group_samples <- metadata$Group == group
    metabolome_group <- metabolome_data_t[group_samples, ]
    
    print(paste("Sample in", group, "group:", sum(group_samples)))
    
    if (sum(group_samples) >= 3) {
      # Create MetaboAnalystR object for group
      mSet_group <- InitDataObjects("conc", "mset", FALSE)
      
      # Set data - we need to create a dummy group variable for the single group
      sample_metadata_group <- data.frame(
        Sample = rownames(metabolome_group),
        Group = rep("Group1", sum(group_samples))
      )
      
      mSet_group <- SetupData(mSet_group, metabolome_group, sample_metadata_group, "group", FALSE)
      
      # Set metabolite annotations
      mSet_group <- Setup.MapData(mSet_group, metabolite_info, "hmdb", "name")
      
      # Perform pathway analysis
      mSet_group <- PerformPathwayAnalysis(mSet_group, "pathora", "hsa", "fdr", 0.05)
      
      # Extract pathway analysis results
      pathway_results_group <- mSet_group$analSet$pathora$pathway.res %>%
        as.data.frame() %>%
        rownames_to_column("Pathway") %>%
        arrange(desc(Impact)) %>%
        mutate(Group = group)
      
      # Save pathway results for group
      write.csv(pathway_results_group, paste0("results/pathway_analysis/MetaboAnalystR/pathway_analysis_results_", group, ".csv"), row.names = FALSE)
      
      # Identify significant pathways for group
      significant_pathways_group <- pathway_results_group %>%
        filter(Adjust.P.Value < 0.05) %>%
        arrange(desc(Impact))
      
      print(paste("Significant pathways for", group, "group (FDR < 0.05):", nrow(significant_pathways_group)))
      
      if (nrow(significant_pathways_group) > 0) {
        # Save significant pathways for group
        write.csv(significant_pathways_group, paste0("results/pathway_analysis/MetaboAnalystR/significant_pathways_", group, ".csv"), row.names = FALSE)
      }
    } else {
      print(paste("Not enough samples in", group, "group for pathway analysis."))
    }
  }
  
  # Compare pathway analysis results between groups
  if (length(unique(metadata$Group)) == 2) {
    group1 <- unique(metadata$Group)[1]
    group2 <- unique(metadata$Group)[2]
    
    # Load pathway results for both groups
    pathway_results_group1 <- read.csv(paste0("results/pathway_analysis/MetaboAnalystR/pathway_analysis_results_", group1, ".csv"), stringsAsFactors = FALSE)
    pathway_results_group2 <- read.csv(paste0("results/pathway_analysis/MetaboAnalystR/pathway_analysis_results_", group2, ".csv"), stringsAsFactors = FALSE)
    
    # Merge results
    pathway_comparison <- pathway_results_group1 %>%
      full_join(pathway_results_group2, by = "Pathway", suffix = c(paste0("_", group1), paste0("_", group2)))
    
    # Save comparison results
    write.csv(pathway_comparison, "results/pathway_analysis/MetaboAnalystR/pathway_comparison_between_groups.csv", row.names = FALSE)
    
    print("Pathway comparison between groups completed successfully!")
  }
} else {
  print("Only one group found. Skipping pathway analysis by group.")
}

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------
cat("\n=== MetaboAnalystR_analysis.R Summary ===\n")
cat("\nGenerated results:\n")
cat("1. Pathway analysis results using MetaboAnalystR\n")
cat("2. Significant pathway identification\n")
cat("3. Pathway impact plot\n")
if (nrow(significant_pathways) > 0) {
  cat("4. Pathway enrichment bar plot\n")
  if (exists("p_heatmap")) {
    cat("5. Pathway metabolite heatmap\n")
  }
  if (file.exists("results/figures/pathway_analysis/MetaboAnalystR/pathway_network.png")) {
    cat("6. Pathway network plot\n")
  }
}
if (length(unique(metadata$Group)) >= 2) {
  cat("7. Pathway analysis by group\n")
  if (length(unique(metadata$Group)) == 2) {
    cat("8. Pathway comparison between groups\n")
  }
}

cat("\nResults saved to: results/pathway_analysis/MetaboAnalystR/\n")
cat("Figures saved to: results/figures/pathway_analysis/MetaboAnalystR/\n")
