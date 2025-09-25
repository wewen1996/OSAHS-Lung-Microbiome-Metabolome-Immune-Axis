# phyloseq_object_creation.R
# Create phyloseq objects for microbiome data analysis

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ape)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# ------------------------------------------------------------------------------
# 1. Load processed data
# ------------------------------------------------------------------------------
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)
microbiome_otu <- read.csv("data/processed/microbiome_otu_processed.csv", row.names = 1, check.names = FALSE)
microbiome_taxonomy <- read.csv("data/processed/microbiome_taxonomy_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load phylogenetic tree if available
if (file.exists("data/processed/microbiome_tree_processed.tree")) {
  microbiome_tree <- read_tree("data/processed/microbiome_tree_processed.tree")
  print("Phylogenetic tree loaded successfully")
} else {
  microbiome_tree <- NULL
  print("Phylogenetic tree not found")
}

# ------------------------------------------------------------------------------
# 2. Create phyloseq components
# ------------------------------------------------------------------------------
# Create OTU matrix (rows = OTUs, columns = samples IDs)
otu_matrix <- as.matrix(microbiome_otu)
otu_table <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Create taxonomy table
tax_matrix <- as.matrix(microbiome_taxonomy)
tax_table <- tax_table(tax_matrix)

# Create sample data
sample_data <- sample_data(metadata)

# Create phylogenetic tree if available
phy_tree <- if (!is.null(microbiome_tree)) {
  phy_tree(microbiome_tree)
} else {
  NULL
}

# ------------------------------------------------------------------------------
# 3. Create phyloseq object
# ------------------------------------------------------------------------------
# Create phyloseq object
if (!is.null(phy_tree)) {
  phyloseq_obj <- phyloseq(otu_table, tax_table, sample_data, phy_tree)
} else {
  phyloseq_obj <- phyloseq(otu_table, tax_table, sample_data)
}

print("Phyloseq object created successfully:")
print(phyloseq_obj)

# ------------------------------------------------------------------------------
# 4. Basic exploration of the phyloseq object
# ------------------------------------------------------------------------------
# Summary of the phyloseq object
cat("\n=== Phyloseq Object Summary ===\n")
cat(paste("Number of OTUs:", ntaxa(phyloseq_obj), "\n"))
cat(paste("Number of samples:", nsamples(phyloseq_obj), "\n"))
cat(paste("Number of sample variables:", ncol(sample_data(phyloseq_obj)), "\n"))
cat(paste("Number of taxonomic ranks:", rank_names(phyloseq_obj), "\n"))

# Check sample metadata
cat("\nSample metadata variables:\n")
print(colnames(sample_data(phyloseq_obj)))

# Check taxonomic ranks
cat("\nTaxonomic ranks:\n")
print(rank_names(phyloseq_obj))

# ------------------------------------------------------------------------------
# 5. Prune the phyloseq object
# ------------------------------------------------------------------------------
# Remove OTUs with zero abundance
phyloseq_obj <- prune_taxa(taxa_sums(phyloseq_obj) > 0, phyloseq_obj)
print(paste("OTUs after removing zero abundance:", ntaxa(phyloseq_obj)))

# Remove samples with zero reads
phyloseq_obj <- prune_samples(sample_sums(phyloseq_obj) > 0, phyloseq_obj)
print(paste("Samples after removing zero read samples:", nsamples(phyloseq_obj)))

# ------------------------------------------------------------------------------
# 6. Agglomerate taxa at different taxonomic levels
# ------------------------------------------------------------------------------
# Agglomerate at phylum level
phyloseq_phylum <- tax_glom(phyloseq_obj, taxrank = "Phylum")
print(paste("Phyla:", ntaxa(phyloseq_phylum)))

# Agglomerate at genus level
phyloseq_genus <- tax_glom(phyloseq_obj, taxrank = "Genus")
print(paste("Genera:", ntaxa(phyloseq_genus)))

# ------------------------------------------------------------------------------
# 7. Save phyloseq objects
# ------------------------------------------------------------------------------
# Save the main phyloseq object
saveRDS(phyloseq_obj, "data/processed/phyloseq_object.rds")

# Save agglomerated objects
saveRDS(phyloseq_phylum, "data/processed/phyloseq_phylum.rds")
saveRDS(phyloseq_genus, "data/processed/phyloseq_genus.rds")

print("Phyloseq objects saved successfully!")

# ------------------------------------------------------------------------------
# 8. Generate basic plots for quality control
# ------------------------------------------------------------------------------
# Create output directory for plots
dir.create("results/figures/quality_control", showWarnings = FALSE, recursive = TRUE)

# Plot sequencing depth
p_depth <- ggplot(data.frame(Depth = sample_sums(phyloseq_obj), 
                             Group = sample_data(phyloseq_obj)$Group),
                  aes(x = Group, y = Depth, fill = Group)) +
  geom_boxplot() +
  labs(title = "Sequencing Depth Distribution",
       x = "Group", y = "Sequencing Depth") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("results/figures/quality_control/sequencing_depth.png", p_depth, 
       width = 10, height = 6, dpi = 300)

# Plot taxonomic distribution at phylum level
p_phylum <- plot_bar(phyloseq_phylum, x = "Group", fill = "Phylum") +
  labs(title = "Taxonomic Distribution at Phylum Level",
       x = "Group", y = "Relative Abundance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

ggsave("results/figures/quality_control/phylum_distribution.png", p_phylum, 
       width = 12, height = 8, dpi = 300)

# Plot alpha diversity (Shannon index)
p_alpha <- plot_richness(phyloseq_obj, x = "Group", measures = c("Shannon", "Simpson")) +
  labs(title = "Alpha Diversity Indices",
       x = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("results/figures/quality_control/alpha_diversity.png", p_alpha, 
       width = 12, height = 6, dpi = 300)

print("Quality control plots generated successfully!")

# ------------------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------------------
cat("\n=== Phyloseq Object Creation Summary ===\n")
cat("Main phyloseq object created with:\n")
cat(paste("  -", ntaxa(phyloseq_obj), "OTUs\n"))
cat(paste("  -", nsamples(phyloseq_obj), "samples\n"))
cat("\nAgglomerated objects:\n")
cat(paste("  - Phylum level:", ntaxa(phyloseq_phylum), "taxa\n"))
cat(paste("  - Genus level:", ntaxa(phyloseq_genus), "taxa\n"))
cat("\nQuality control plots saved to: results/figures/quality_control/\n")
