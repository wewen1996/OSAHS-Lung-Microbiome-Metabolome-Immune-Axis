# microbiome_t_test.R
# Perform t-test analysis on microbiome data to identify differentially abundant taxa

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(ggpubr)
library(permute)
library(vegan)
library(multcomp)
library(broom)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/tables/differential_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/differential_analysis", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load normalized phyloseq objects
phyloseq_tss <- readRDS("data/processed/normalized/phyloseq_tss.rds")
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")

# Load metadata
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Check the grouping variable
print("Grouping variable levels:")
print(table(sample_data(phyloseq_tss)$Group))

# ------------------------------------------------------------------------------
# 2. Prepare data for t-test
# ------------------------------------------------------------------------------
# Use CSS normalized data for differential analysis
phyloseq_obj <- phyloseq_css

# Extract OTU table and metadata
otu_table <- otu_table(phyloseq_obj) %>% as.matrix()
group <- sample_data(phyloseq_obj)$Group

# Ensure samples are in the same order
stopifnot(all(colnames(otu_table) == rownames(metadata)))

# Split samples into groups
groups <- unique(group)
if (length(groups) != 2) {
  stop("This script is designed for two-group comparison only.")
}

group1_samples <- names(group[group == groups[1]])
group2_samples <- names(group[group == groups[2]])

print(paste("Comparing", groups[1], "vs", groups[2]))
print(paste("Group 1 samples:", length(group1_samples)))
print(paste("Group 2 samples:", length(group2_samples)))

# ------------------------------------------------------------------------------
# 3. Perform t-test on each taxon
# ------------------------------------------------------------------------------
cat("\nPerforming t-test on each taxon...\n")

# Function to perform t-test for each taxon
perform_t_test <- function(otu_table, group1_samples, group2_samples) {
  results <- data.frame(
    Taxon = rownames(otu_table),
    Mean_Group1 = NA,
    Mean_Group2 = NA,
    Fold_Change = NA,
    Log2_Fold_Change = NA,
    T_Statistic = NA,
    P_Value = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(otu_table)) {
    taxon <- otu_table[i, ]
    group1_values <- taxon[group1_samples]
    group2_values <- taxon[group2_samples]
    
    # Perform t-test
    t_test_result <- t.test(group1_values, group2_values)
    
    # Calculate means
    mean_group1 <- mean(group1_values)
    mean_group2 <- mean(group2_values)
    
    # Calculate fold change (group1 / group2)
    fold_change <- mean_group1 / mean_group2
    log2_fold_change <- log2(fold_change)
    
    # Store results
    results$Mean_Group1[i] <- mean_group1
    results$Mean_Group2[i] <- mean_group2
    results$Fold_Change[i] <- fold_change
    results$Log2_Fold_Change[i] <- log2_fold_change
    results$T_Statistic[i] <- t_test_result$statistic
    results$P_Value[i] <- t_test_result$p.value
  }
  
  # Add adjusted p-values (FDR)
  results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
  
  # Add significance level
  results$Significance <- case_when(
    results$Adjusted_P_Value < 0.001 ~ "***",
    results$Adjusted_P_Value < 0.01 ~ "**",
    results$Adjusted_P_Value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # Order by adjusted p-value
  results <- results %>% arrange(Adjusted_P_Value)
  
  return(results)
}

# Perform t-test
t_test_results <- perform_t_test(otu_table, group1_samples, group2_samples)

# Add taxonomy information
taxonomy <- tax_table(phyloseq_obj) %>% as.data.frame()
t_test_results <- cbind(t_test_results, taxonomy[match(t_test_results$Taxon, rownames(taxonomy)), ])

# Print summary of results
cat("\n=== T-test Results Summary ===\n")
print(paste("Total taxa tested:", nrow(t_test_results)))
print(paste("Significant taxa (FDR < 0.05):", sum(t_test_results$Adjusted_P_Value < 0.05)))
print(paste("Highly significant taxa (FDR < 0.01):", sum(t_test_results$Adjusted_P_Value < 0.01)))

# ------------------------------------------------------------------------------
# 4. Save results
# ------------------------------------------------------------------------------
# Save full results
write.csv(t_test_results, "results/tables/differential_analysis/microbiome_t_test_results.csv", row.names = FALSE)

# Save significant results (FDR < 0.05)
significant_results <- t_test_results %>% filter(Adjusted_P_Value < 0.05)
write.csv(significant_results, "results/tables/differential_analysis/microbiome_t_test_significant.csv", row.names = FALSE)

print(paste("Significant results saved. Found", nrow(significant_results), "significant taxa."))

# ------------------------------------------------------------------------------
# 5. Generate volcano plot
# ------------------------------------------------------------------------------
cat("\nGenerating volcano plot...\n")

# Prepare data for volcano plot
volcano_data <- t_test_results %>%
  mutate(
    Neg_Log10_P_Adj = -log10(Adjusted_P_Value),
    Significance = case_when(
      Adjusted_P_Value < 0.05 & abs(Log2_Fold_Change) > 1 ~ "Significant (FDR < 0.05, |log2FC| > 1)",
      Adjusted_P_Value < 0.05 ~ "Significant (FDR < 0.05)",
      TRUE ~ "Not significant"
    )
  )

# Create volcano plot
p_volcano <- ggplot(volcano_data, aes(x = Log2_Fold_Change, y = Neg_Log10_P_Adj, color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Volcano Plot of Differential Abundance Analysis",
       x = paste("Log2 Fold Change (", groups[1], " vs ", groups[2], ")", sep = ""),
       y = "-Log10(Adjusted P-value)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("black", "orange", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Add labels for top significant taxa
top_taxa <- volcano_data %>%
  filter(Adjusted_P_Value < 0.01 & abs(Log2_Fold_Change) > 1.5) %>%
  arrange(Adjusted_P_Value) %>%
  head(10)

if (nrow(top_taxa) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = top_taxa,
      aes(label = Taxon),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.5,
      segment.color = "gray50"
    )
}

# Save volcano plot
ggsave("results/figures/differential_analysis/microbiome_volcano_plot.png", p_volcano,
       width = 12, height = 8, dpi = 300)

print("Volcano plot generated successfully!")

# ------------------------------------------------------------------------------
# 6. Generate heatmap of significant taxa
# ------------------------------------------------------------------------------
cat("\nGenerating heatmap of significant taxa...\n")

# If there are significant taxa, create heatmap
if (nrow(significant_results) > 0) {
  # Get top N significant taxa (up to 50)
  top_n <- min(50, nrow(significant_results))
  top_taxa <- significant_results$Taxon[1:top_n]
  
  # Subset OTU table to top taxa
  top_otu_table <- otu_table[top_taxa, ]
  
  # Add taxonomy information to row names
  tax_labels <- apply(taxonomy[top_taxa, c("Phylum", "Genus")], 1, function(x) {
    paste(x[!is.na(x)], collapse = "|")
  })
  rownames(top_otu_table) <- tax_labels
  
  # Create annotation for samples groups
  sample_annot <- data.frame(Group = group)
  rownames(sample_annot) <- colnames(top_otu_table)
  
  # Create heatmap
  pheatmap(
    top_otu_table,
    scale = "row",
    annotation_col = sample_annot,
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = paste("Top", top_n, "Significant Taxa (FDR < 0.05)"),
    filename = "results/figures/differential_analysis/microbiome_heatmap_significant.png",
    width = 10,
    height = 10,
    dpi = 300
  )
  
  print(paste("Heatmap generated with top", top_n, "significant taxa."))
} else {
  print("No significant taxa found. Skipping heatmap generation.")
}

# ------------------------------------------------------------------------------
# 7. Generate boxplots for top significant taxa
# ------------------------------------------------------------------------------
cat("\nGenerating boxplots for top significant taxa...\n")

# Get top 5 significant taxa
top_5_taxa <- significant_results$Taxon[1:min(5, nrow(significant_results))]

if (length(top_5_taxa) > 0) {
  # Prepare data for boxplots
  boxplot_data <- as.data.frame(t(otu_table[top_5_taxa, ])) %>%
    rownames_to_column("Sample") %>%
    mutate(Group = group[match(Sample, names(group))]) %>%
    pivot_longer(cols = -c(Sample, Group), names_to = "Taxon", values_to = "Abundance")
  
  # Add taxonomy information
  tax_labels <- apply(taxonomy[top_5_taxa, c("Phylum", "Genus")], 1, function(x) {
    paste(x[!is.na(x)], collapse = "|")
  })
  names(tax_labels) <- top_5_taxa
  boxplot_data$Taxon_Label <- tax_labels[boxplot_data$Taxon]
  
  # Create boxplots
  p_boxplot <- ggplot(boxplot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~Taxon_Label, scales = "free_y") +
    labs(title = "Top Significant Taxa Abundance by Group",
         x = "Group",
         y = "Relative Abundance (CSS Normalized)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Save boxplots
  ggsave("results/figures/differential_analysis/microbiome_boxplots_top_significant.png", p_boxplot,
         width = 15, height = 10, dpi = 300)
  
  print("Boxplots generated for top significant taxa.")
} else {
  print("No significant taxa found. Skipping boxplot generation.")
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== Microbiome_t_test.R Summary ===\n")
print(paste("Comparison:", groups[1], "vs", groups[2]))
print(paste("Total taxa tested:", nrow(t_test_results)))
print(paste("Significant taxa (FDR < 0.05):", sum(t_test_results$Adjusted_P_Value < 0.05)))
print(paste("Highly significant taxa (FDR < 0.01):", sum(t_test_results$Adjusted_P_Value < 0.01)))

if (sum(t_test_results$Adjusted_P_Value < 0.05) > 0) {
  print("\nTop 5 significant taxa:")
  print(significant_results %>%
          select(Taxon, Phylum, Genus, Log2_Fold_Change, Adjusted_P_Value) %>%
          head(5))
}

print("\nFiles generated:")
print("1. results/tables/differential_analysis/microbiome_t_test_results.csv")
print("2. results/tables/differential_analysis/microbiome_t_test_significant.csv")
print("3. results/figures/differential_analysis/microbiome_volcano_plot.png")
if (nrow(significant_results) > 0) {
  print("4. results/figures/differential_analysis/microbiome_heatmap_significant.png")
  if (length(top_5_taxa) > 0) {
    print("5. results/figures/differential_analysis/microbiome_boxplots_top_significant.png")
  }
}
