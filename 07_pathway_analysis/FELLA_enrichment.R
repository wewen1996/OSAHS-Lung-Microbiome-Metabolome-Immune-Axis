# FELLA_enrichment.R
# Perform functional enrichment analysis using FELLA package

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(FELLA)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/pathway_analysis/FELLA", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/pathway_analysis/FELLA", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load differential abundance data
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load differential differential analysis results
microBIOME_t_test_results <- read.csv("results/tables/differential_analysis/microbiome_t_test_results.csv", stringsAsFactors = FALSE)
ancombc_results <- read.csv("results/tables/differential_analysis/ancombc_results.csv", stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Prepare data for FELLA
# ------------------------------------------------------------------------------
# Get taxonomy table
taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()

# Function to map taxonomy to KEGG IDs
# This requires a mapping file between taxonomy and KEGG IDs
map_taxonomy_to_kegg <- function(taxonomy_df, mapping_file = "data/metadata/taxonomy_to_kegg_mapping.csv") {
  
  if (!file.exists(mapping_file)) {
    warning("Taxonomy to KEGG mapping file not found. Please create this file with columns: Taxonomy, KEGG_ID")
    return(NULL)
  }
  
  # Load mapping file
  kegg_mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
  
  # Create taxonomy string
  taxonomy_df <- taxonomy_df %>%
    rownames_to_column("OTU") %>%
    mutate(
      Taxonomy_String = case_when(
        !is.na(Species) ~ paste(Phylum, Class, Order, Family, Genus, Species, sep = "|"),
        !is.na(Genus) ~ paste(Phylum, Class, Order, Family, Genus, sep = "|"),
        !is.na(Family) ~ paste(Phylum, Class, Order, Family, sep = "|"),
        !is.na(Order) ~ paste(Phylum, Class, Order, sep = "|"),
        !is.na(Class) ~ paste(Phylum, Class, sep = "|"),
        TRUE ~ Phylum
      )
    )
  
  # Map to KEGG IDs
  mapped_taxonomy <- taxonomy_df %>%
    left_join(kegg_mapping, by = c("Taxonomy_String" = "Taxonomy"))
  
  return(mapped_taxonomy)
}

# Map taxonomy to KEGG IDs
mapped_taxonomy <- map_taxonomy_to_kegg(taxonomy)

if (!is.null(mapped_taxonomy)) {
  # Remove rows with missing KEGG IDs
  mapped_taxonomy <- mapped_taxonomy %>% filter(!is.na(KEgg_ID))
  
  print(paste("Taxa with KEGG IDs:", nrow(mapped_taxonomy)))
  
  if (nrow(mapped_taxonomy) == 0) {
    stop("No taxa could be mapped to KEGG IDs. Cannot perform FELLA enrichment analysis.")
  }
} else {
  stop("Could not map taxonomy to KEGG IDs. Please check your mapping file.")
}

# ------------------------------------------------------------------------------
# 3. Prepare differential abundance data for enrichment analysis
# ------------------------------------------------------------------------------
# Get significantly differential abundant taxa from t-test
significant_taxa_t_test <- microbiome_t_test_results %>%
  filter(Adjusted_P_Value < 0.05) %>%
  pull(Taxon)

# Get significantly differential abundant taxa from ANCOM-BC
significant_taxa_ancombc <- ancombc_results %>%
  filter(q_val < 0.05) %>%
  pull(OTU)

# Combine significant taxa from different methods
significant_taxa <- unique(c(significant_taxa_t_test, significant_taxa_ancombc))

print(paste("Significant taxa from differential analysis:", length(significant_taxa)))

# Map significant taxa to KEGG IDs
significant_kegg_ids <- mapped_taxonomy %>%
  filter(OTU %in% significant_taxa & !is.na(KEGG_ID)) %>%
  pull(KEGG_ID) %>%
  unique()

print(paste("Significant taxa with KEGG IDs:", length(significant_kegg_ids)))

if (length(significant_kegg_ids) == 0) {
  stop("No significant taxa could be mapped to KEGG IDs. Cannot perform FELLA enrichment analysis.")
}

# Get all KEGG IDs for background
background_kegg_ids <- mapped_taxonomy %>%
  filter(!is.na(KEGG_ID)) %>%
  pull(KEGG_ID) %>%
  unique()

print(paste("Background KEGG IDs:", length(background_kegg_ids)))

# ------------------------------------------------------------------------------
# 4. Perform FELLA enrichment analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing FELLA enrichment analysis ===\n")

# Load KEGG database for FELLA
# This requires downloading the database first using FELLA::downloadKEGG()
if (!file.exists("data/databases/FELLA/hsa")) {
  cat("Downloading KEGG database for FELLA...\n")
  FELLA::downloadKEGG(
    organism = "hsa",
    databaseDir = "data/databases/FELLA",
    internalDir = TRUE
  )
}

# Load FELLA database
fella_db <- FELLA::loadKEGG(
  organism = "hsa",
  databaseDir = "data/databases/FELLA",
  internalDir = TRUE
)

# Perform enrichment analysis
set.seed(123)
fella_results <- FELLA::runEnrichment(
  ids = significant_kegg_ids,
  data = fella_db,
  method = "hypergeom",
  background = background_kegg_ids
)

# Save FELLA results
saveRDS(fella_results, "results/pathway_analysis/FELLA/fella_enrichment_results.rds")

print("FELLA enrichment analysis completed successfully!")

# ------------------------------------------------------------------------------
# 5. Extract and visualize FELLA results
# ------------------------------------------------------------------------------
cat("\n=== Extracting and visualizing FELLA results ===\n")

# Extract enrichment results
enrichment_results <- FELLA::getResults(
  fella_results,
  level = "pathway",
  format = "data.frame"
)

if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
  # Add pathway names
  pathway_names <- FELLA::getNames(fella_db, ids = enrichment_results$id, level = "pathway")
  enrichment_results$Pathway_Name <- pathway_names[match(enrichment_results$id, names(pathway_names))]
  
  # Reorder columns
  enrichment_results <- enrichment_results %>%
    select(id, Pathway_Name, p.value, padj, overlap, size, everything()) %>%
    arrange(padj)
  
  # Save enrichment results
  write.csv(enrichment_results, "results/pathway_analysis/FELLA/fella_enrichment_results.csv", row.names = FALSE)
  
  print("FELLA enrichment results extracted successfully!")
  print(paste("Total pathways tested:", nrow(enrichment_results)))
  
  # Identify significant pathway
  significant_pathways <- enrichment_results %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  print(paste("Significant pathway (FDR < 0.05):", nrow(significant_pathways)))
  
  if (nrow(significant_pathways) > 0) {
    # Save significant pathway
    write.csv(significant_pathways, "results/pathway_analysis/FELLA/fella_significant_pathways.csv", row.names = FALSE)
    
    print("Top 10 significant pathway:")
    print(head(significant_pathways[, c("id", "Pathway_Name", "p.value", "padj", "overlap", "size")], 10))
  }
} else {
  print("No enrichment results found.")
  significant_pathways <- data.frame()
}

# ------------------------------------------------------------------------------
# 6. Visualize FELLA results
# ------------------------------------------------------------------------------
cat("\n=== Visualizing FELLA results ===\n")

# 6.1 Enrichment bar plot
if (nrow(significant_pathways) > 0) {
  # Select top N pathway
  top_n <- min(20, nrow(significant_pathways))
  top_pathways <- significant_pathways %>% head(top_n)
  
  # Create bar plot
  p_bar <- ggplot(top_pathways, aes(x = reorder(Pathway_Name, -log10(padj)), y = -log10(padj), fill = overlap/size)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(name = "Overlap Ratio                        (Hit/Size)") +
    coord_flip() +
    labs(
      title = "FELLA Enrichment Results",
      x = "Pathway",
      y = "-Log10(Adjusted P-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
  
  # Save plot
  ggsave("results/figures/pathway_analysis/FELLA/fella_enrichment_barplot.png", p_bar,
         width = 14, height = 12, dpi = 300)
  
  print("FELLA enrichment bar plot generated successfully!")
}

# 6.2 Enrichment bubble plot
if (nrow(enrichment_results) > 0) {
  # Create bubble plot
  p_bubble <- ggplot(enrichment_results, aes(x = -log10(p.value), y = overlap/size, size = size, color = padj)) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_viridis_c(direction = -1, name = "Adjusted P-value") +
    scale_size(name = "Pathway Size") +
    labs(
      title = "FELLA Enrichment Bubble Plot",
      x = "-Log10(P-value)",
      y = "Overlap Ratio (Hit/Size)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Add labels for significant pathway
  if (nrow(significant_pathways) > 0) {
    top_pathways <- significant_pathways %>% head(10)
    
    p_bubble <- p_bubble +
      geom_text_repel(
        data = top_pathways,
        aes(label = Pathway_Name),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.5,
        segment.color = "gray50"
      )
  }
  
  # Save plot
  ggsave("results/figures/pathway_analysis/FELLA/fella_enrichment_bubble.png", p_bubble,
         width = 12, height = 10, dpi = 300)
  
  print("FELLA enrichment bubble plot generated successfully!")
}

# 6.3 Pathway network
if (nrow(significant_pathways) > 0) {
  # Try to generate pathway network if available
  tryCatch({
    # Create network plot
    png("results/figures/pathway_analysis/FELLA/fella_pathway_network.png", width = 1200, height = 1000, res = 300)
    FELLA::plotNetwork(
      fella_results,
      data = fella_db,
      level = "pathway",
      nlimit = 20,
      layout = "force",
      nodeColor = "padj",
      nodeSize = "size",
      showLabels = TRUE,
      labelSize = 3
    )
    dev.off()
    
    print("FELLA pathway network plot generated successfully!")
  }, error = function(e) {
    print("Could not generate FELLA pathway network plot.")
    print(e$message)
  })
}

# 6.4 Heatmap of pathway coverage
if (nrow(significant_pathways) > 0) {
  # Get pathway-gene association
  pathway_gene <- FELLA::getNeighbors(fella_db, ids = significant_pathways$id, level = "pathway", method = "children")
  
  # Create heatmap data
  if (!is.null(pathway_gene) && length(pathway_gene) > 0) {
    # Convert to data frame
    heatmap_data <- map_dfr(names(pathway_gene), function(pathway) {
      data.frame(
        Pathway = pathway,
        Gene = pathway_gene[[pathway]],
        stringsAsFactors = FALSE
      )
    }) %>%
      mutate(Value = 1) %>%
      pivot_wider(names_from = Pathway, values_from = Value, values_fill = 0)
    
    # Get pathway names
    pathway_names <- FELLA::getNames(fella_db, ids = colnames(heatmap_data)[-1], level = "pathway")
    
    # Rename columns with pathway names
    colnames(heatmap_data)[-1] <- pathway_names[match(colnames(heatmap_data)[-1], names(pathway_names))]
    
    # Create heatmap
    p_heatmap <- pheatmap(
      as.matrix(heatmap_data[, -1]),
      scale = "none",
      show_rownames = TRUE,
      show_colnames = TRUE,
      treeheight_row = 20,
      treeheight_col = 20,
      fontsize_row = 8,
      fontsize_col = 8,
      main = "Pathway Gene Coverage Heatmap",
      color = c("white", "red"),
      filename = "results/figures/pathway_analysis/FELLA/fella_pathway_gene_heatmap.png",
      width = 14,
      height = 12,
      dpi = 300
    )
    
    print("FELLA pathway gene coverage heatmap generated successfully!")
  }
}

# ------------------------------------------------------------------------------
# 7. Compare enrichment results between groups (if possible)
# ------------------------------------------------------------------------------
cat("\n=== Comparing enrichment results between groups ===\n")

# Check if there are multiple groups and if we have results for each group
if (length(unique(metadata$Group)) >= 2) {
  # This would require running FELLA separately for each group
  # and then comparing the results
  
  print("Group comparison not implemented in this version.")
  print("To compare between groups, run FELLA separately for each group and then compare the results.")
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== FELLA_enrichment.R Summary ===\n")
cat("\nGenerated results:\n")
cat("1. FELLA enrichment analysis results\n")
cat("2. Significant pathway identification\n")
if (nrow(significant_pathways) > 0) {
  cat("3. FELLA enrichment bar plot\n")
  cat("4. FELLA enrichment bubble plot\n")
  if (file.exists("results/figures/pathway_analysis/FELLA/fella_pathway_network.png")) {
    cat("5. FELLA pathway network plot\n")
  }
  if (exists("p_heatmap")) {
    cat("6. FELLA pathway gene coverage heatmap\n")
  }
}

cat("\nResults saved to: results/pathway_analysis/FELLA/\n")
cat("Figures saved to: results/figures/pathway_analysis/FELLA/\n")
