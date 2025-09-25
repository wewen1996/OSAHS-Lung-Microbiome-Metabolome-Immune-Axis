# corncob_maaslin_ancombc.R
# Compare multiple differential abundance analysis methods for microbiome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(corncob)
library(Maaslin2)
library(ANCOMBC)
library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(pheatmap)
library(ggpubr)
library(VennDiagram)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/tables/differential_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/differential_analysis", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load raw phyloseq object (not normalized for these methods)
phyloseq_obj <- readRDS("data/processed/phyloseq_object.rds")

# Load metadata
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")
print(paste("Microbiome samples:", nsamples(phyloseq_obj)))
print(paste("Microbiome features:", ntaxa(phyloseq_obj)))

# Check grouping variable
print("Grouping variable levels:")
print(table(sample_data(phyloseq_obj)$Group))

# ------------------------------------------------------------------------------
# 2. Prepare data for differential abundance analysis
# ------------------------------------------------------------------------------
# Extract OTU table and metadata
otu_table <- otu_table(phyloseq_obj) %>% as.matrix()
taxonomy <- tax_table(phyloseq_obj) %>% as.data.frame()
sample_data <- sample_data(phyloseq_obj) %>% as.data.frame()

# Ensure samples IDs match
stopifnot(all(colnames(otu_table) == rownames(sample_data)))

# Filter low abundance OTUs (optional but recommended for some methods)
# Keep OTUs with at least 10 reads counts in at least 10% of samples
min_reads <- 10
min_samples <- 0.1 * ncol(otu_table)
otu_table_filtered <- otu_table[rowSums(otu_table >= min_reads) >= min_samples, ]

print(paste("OTUs after filtering:", nrow(otu_table_filtered)))
print(paste("Original OTUs:", nrow(otu_table)))

# Update taxonomy table
taxonomy_filtered <- taxonomy[rownames(otu_table_filtered), ]

# ------------------------------------------------------------------------------
# 3. Differential abundance analysis using ANCOM-BC
# ------------------------------------------------------------------------------
cat("\n=== Performing ANCOM-BC analysis ===\n")

# Prepare data for ANCOM-BC
# Create a SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = otu_table_filtered),
  colData = DataFrame(sample_data)
)

# Run ANCOM-BC
set.seed(123)
ancombc_result <- ancombc(
  data = se,
  tax_level = "OTU",
  formula = "Group",
  p_adj_method = "fdr",
  zero_cut = 0.90,
  lib_cut = 1000,
  group = "Group",
  out_cut = 0.05,
  verbose = FALSE
)

# Extract results
ancombc_results <- ancombc_result$res

# Add taxonomy information
ancombc_results <- ancombc_results %>%
  rownames_to_column("OTU") %>%
  left_join(taxonomy_filtered %>% rownames_to_column("OTU"), by = "OTU") %>%
  arrange(q_val)

# Identify significant taxa (FDR < 0.05)
ancombc_significant <- ancombc_results %>%
  filter(q_val < 0.05) %>%
  arrange(q_val)

print(paste("ANCOM-BC significant taxa (FDR < 0.05):", nrow(ancombc_significant)))

# ------------------------------------------------------------------------------
# 4. Differential abundance analysis using Maaslin2
# ------------------------------------------------------------------------------
cat("\n=== Performing Maaslin2 analysis ===\n")

# Prepare data for Maaslin2
# Transpose OTU table (samples as rows, OTUs as columns)
maaslin_input <- data.frame(t(otu_table_filtered), check.names = FALSE)

# Add metadata
maaslin_input <- cbind(maaslin_input, sample_data)

# Run Maaslin2
set.seed(123)
maaslin_result <- Maaslin2(
  input_data = maaslin_input,
  output = "results/tables/differential_analysis/maaslin2_output",
  fixed_effects = "Group",
  reference = "Control",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "fdr",
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

# Extract results
maaslin_results <- maaslin_result$results %>%
  rownames_to_column("OTU") %>%
  left_join(taxonomy_filtered %>% rownames_to_column("OTU"), by = "OTU") %>%
  arrange(qval)

# Identify significant taxa (FDR < 0.05)
maaslin_significant <- maaslin_results %>%
  filter(qval < 0.05) %>%
  arrange(qval)

print(paste("Maaslin2 significant taxa (FDR < 0.05):", nrow(maaslin_significant)))

# ------------------------------------------------------------------------------
# 5. Differential abundance analysis using corncob
# ------------------------------------------------------------------------------
cat("\n=== Performing corncob analysis ===\n")

# Prepare data for corncob
# Convert to DGEList for normalization
dge <- DGEList(counts = otu_table_filtered)

# Apply TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Create phyloseq object with normalized data
phyloseq_norm <- phyloseq_obj
otu_table(phyloseq_norm) <- otu_table(cpm(dge, log = TRUE), taxa_are_rows = TRUE)

# Run corncob
set.seed(123)
corncob_result <- differentialTest(
  formula = ~ Group,
  phi.formula = ~ Group,
  formula_null = ~ 1,
  phi.formula_null = ~ Group,
  test = "Wald",
  boot = FALSE,
  data = phyloseq_norm,
  fdr_cutoff = 0.05
)

# Extract results
corncob_results <- corncob_result$pvalues %>%
  rownames_to_column("OTU") %>%
  left_join(taxonomy_filtered %>% rownames_to_column("OTU"), by = "OTU") %>%
  arrange(p_fdr)

# Identify significant taxa (FDR < 0.05)
corncob_significant <- corncob_results %>%
  filter(p_fdr < 0.05) %>%
  arrange(p_fdr)

print(paste("corncob significant taxa (FDR < 0.05):", nrow(corncob_significant)))

# ------------------------------------------------------------------------------
# 6. Compare results across methods
# ------------------------------------------------------------------------------
cat("\n=== Comparing results across methods ===\n")

# Get significant OTUs from each method
ancombc_otus <- ancombc_significant$OTU
maaslin_otus <- maaslin_significant$OTU
corncob_otus <- corncob_significant$OTU

# Create Venn diagram data
venn_data <- list(
  ANCOM_BC = ancombc_otus,
  Maaslin2 = maaslin_otus,
  corncob = corncob_otus
)

# Print overlap statistics
cat("\nMethod comparison statistics:\n")
print(paste("ANCOM-BC significant OTUs:", length(ancombc_otus)))
print(paste("Maaslin2 significant OTUs:", length(maaslin_otus)))
print(paste("corncob significant OTUs:", length(corncob_otus)))
print(paste("Common to all methods:", length(Reduce(intersect, venn_data))))

# Calculate pairwise overlaps
pairwise_overlaps <- combn(names(venn_data), 2, function(x) {
  overlap <- length(intersect(venn_data[[x[1]]], venn_data[[x[2]]]))
  return(paste(x[1], "vs", x[2], ":", overlap))
})

print("Pairwise overlaps:")
print(pairwise_overlaps)

# ------------------------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------------------------
# Save full results
write.csv(ancombc_results, "results/tables/differential_analysis/ancombc_results.csv", row.names = FALSE)
write.csv(maaslin_results, "results/tables/differential_analysis/maaslin2_results.csv", row.names = FALSE)
write.csv(corncob_results, "results/tables/differential_analysis/corncob_results.csv", row.names = FALSE)

# Save significant results
write.csv(ancombc_significant, "results/tables/differential_analysis/ancombc_significant.csv", row.names = FALSE)
write.csv(maaslin_significant, "results/tables/differential_analysis/maaslin2_significant.csv", row.names = FALSE)
write.csv(corncob_significant, "results/tables/differential_analysis/corncob_significant.csv", row.names = FALSE)

# Save comparison results
comparison_results <- data.frame(
  Method = c("ANCOM-BC", "Maaslin2", "corncob"),
  Significant_Taxa = c(length(ancombc_otus), length(maaslin_otus), length(corncob_otus)),
  stringsAsFactors = FALSE
)

write.csv(comparison_results, "results/tables/differential_analysis/method_comparison.csv", row.names = FALSE)

print("\nResults saved successfully!")

# ------------------------------------------------------------------------------
# 8. Generate visualization
# ------------------------------------------------------------------------------
# 8.1 Venn diagram of significant OTUs
cat("\nGenerating Venn diagram...\n")

# Create Venn diagram
venn_plot <- venn.diagram(
  x = venn_data,
  filename = NULL,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cat.col = c("red", "blue", "green"),
  cat.cex = 1.2,
  main = "Comparison of Significant OTUs Across Methods",
  main.cex = 1.5
)

# Save Venn diagram
png("results/figures/differential_analysis/method_comparison_venn.png", width = 800, height = 800, res = 300)
grid.draw(venn_plot)
dev.off()

print("Venn diagram generated successfully!")

# 8.2 Bar plot of method comparison
cat("\nGenerating method comparison bar plot...\n")

# Create bar plot
p_method_comparison <- ggplot(comparison_results, aes(x = Method, y = Significant_Taxa, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Significant_Taxa), vjust = -0.3, size = 5) +
  labs(title = "Number of Significant Taxa Identified by Each Method",
       x = "Method",
       y = "Number of Significant Taxa (FDR < 0.05)") +
  scale_fill_manual(values = c("red", "blue", "green")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

# Save bar plot
ggsave("results/figures/differential_analysis/method_comparison_bar.png", p_method_comparison,
       width = 10, height = 8, dpi = 300)

print("Method comparison bar plot generated successfully!")

# 8.3 Heatmap of common significant OTUs
cat("\nGenerating heatmap of common significant OTUs...\n")

# Get common significant OTUs across all methods
common_otus <- Reduce(intersect, venn_data)

if (length(common_otus) > 0) {
  # Get normalized abundance data (CSS normalized)
  phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
  otu_css <- otu_table(phyloseq_css) %>% as.matrix()
  
  # Subset to common significant OTUs
  common_otu_data <- otu_css[common_otus, ]
  
  # Add taxonomy information to row names
  tax_labels <- apply(taxonomy[common_otus, c("Phylum", "Genus")], 1, function(x) {
    paste(x[!is.na(x)], collapse = "|")
  })
  rownames(common_otu_data) <- tax_labels
  
  # Create annotation for samples groups
  sample_annot <- data.frame(Group = sample_data$Group)
  rownames(sample_annot) <- colnames(common_otu_data)
  
  # Create heatmap
  pheatmap(
    common_otu_data,
    scale = "row",
    annotation_col = sample_annot,
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = "Common Significant OTUs Across All Methods",
    filename = "results/figures/differential_analysis/common_significant_otus_heatmap.png",
    width = 10,
    height = 8,
    dpi = 300
  )
  
  print(paste("Heatmap generated for", length(common_otus), "common significant OTUs."))
} else {
  print("No common significant OTUs found across all methods. Skipping heatmap generation.")
}

# 8.4 Volcano plots for each method
cat("\nGenerating volcano plots for each method...\n")

# Function to create volcano plot
create_volcano_plot <- function(results, method_name, x_col, p_col, title) {
  volcano_data <- results %>%
    mutate(
      Neg_Log10_P = -log10(!!sym(p_col)),
      Significant = !!sym(p_col) < 0.05
    ) %>%
    mutate(
      Phylum_Genus = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), Phylum)
    )
  
  p <- ggplot(volcano_data, aes(x = !!sym(x_col), y = Neg_Log10_P, color = Significant)) +
    geom_point(alpha = 0.7, size = 2) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10(Adjusted P-value)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add labels for top significant taxa
  top_taxa <- volcano_data %>%
    filter(Significant) %>%
    arrange(desc(Neg_Log10_P)) %>%
    head(5)
  
  if (nrow(top_taxa) > 0) {
    p <- p +
      geom_text_repel(
        data = top_taxa,
        aes(label = Phylum_Genus),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.5,
        segment.color = "gray50"
      )
  }
  
  return(p)
}

# Create volcano plots
if (nrow(ancombc_results) > 0) {
  p_ancombc_volcano <- create_volcano_plot(
    ancombc_results,
    "ANCOM-BC",
    "lfc",
    "q_val",
    "ANCOM-BC Differential Abundance Analysis"
  )
  
  ggsave("results/figures/differential_analysis/ancombc_volcano.png", p_ancombc_volcano,
         width = 12, height = 8, dpi = 300)
}

if (nrow(maaslin_results) > 0) {
  p_maaslin_volcano <- create_volcano_plot(
    maaslin_results,
    "Maaslin2",
    "coef",
    "qval",
    "Maaslin2 Differential Abundance Analysis"
  )
  
  ggsave("results/figures/differential_analysis/maaslin2_volcano.png", p_maaslin_volcano,
         width = 12, height = 8, dpi = 300)
}

if (nrow(corncob_results) > 0) {
  # For corncob, we need to calculate log2 fold change
  corncob_results <- corncob_results %>%
    mutate(
      log2FoldChange = log2(estimate / estimate_null)
    )
  
  p_corncob_volcano <- create_volcano_plot(
    corncob_results,
    "corncob",
    "log2FoldChange",
    "p_fdr",
    "corncob Differential Abundance Analysis"
  )
  
  ggsave("results/figures/differential_analysis/corncob_volcano.png", p_corncob_volcano,
         width = 12, height = 8, dpi = 300)
}

print("Volcano plots generated successfully!")

# ------------------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------------------
cat("\n=== corncob_maaslin_ancombc.R Summary ===\n")
cat("\nMethod comparison results:\n")
print(comparison_results)

cat("\nOverlap statistics:\n")
print(paste("Common to all methods:", length(Reduce(intersect, venn_data))))
print("Pairwise overlaps:")
print(pairwise_overlaps)

if (length(common_otus) > 0) {
  print(paste("\nTop common significant OTUs (", length(common_otus), " total):", sep = ""))
  
  # Get common OTUs results from ANCOM-BC (as an example)
  common_results <- ancombc_significant %>%
    filter(OTU %in% common_otus) %>%
    select(OTU, Phylum, Genus, lfc, q_val) %>%
    head(10)
  
  print(common_results)
}

print("\nFiles generated:")
print("1. results/tables/differential_analysis/ancombc_results.csv")
print("2. results/tables/differential_analysis/maaslin2_results.csv")
print("3. results/tables/differential_analysis/corncob_results.csv")
print("4. results/tables/differential_analysis/ancombc_significant.csv")
print("5. results/tables/differential_analysis/maaslin2_significant.csv")
print("6. results/tables/differential_analysis/corncob_significant.csv")
print("7. results/tables/differential_analysis/method_comparison.csv")
print("8. results/figures/differential_analysis/method_comparison_venn.png")
print("9. results/figures/differential_analysis/method_comparison_bar.png")
if (length(common_otus) > 0) {
  print("10. results/figures/differential_analysis/common_significant_otus_heatmap.png")
}
print("11. results/figures/differential_analysis/ancombc_volcano.png")
print("12. results/figures/differential_analysis/maaslin2_volcano.png")
print("13. results/figures/differential_analysis/corncob_volcano.png")
