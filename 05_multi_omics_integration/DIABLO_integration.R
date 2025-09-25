# DIABLO_integration.R
# Perform multi-omics integration using DIABLO method from mixOmics

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(mixOmics)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/multi_omics_integration/DIABLO", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/multi_omics_integration/DIABLO", showWarnings = FALSE, recursive = TRUE)

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
# 2. Prepare data for DIablo analysis
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(sample_names(phyloseq_css), colnames(metabolome_data))

if (length(common_samples) < 3) {
  stop("Not enough common samples for multi-omics integration.")
}

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples IDs:", length(common_samples)))

# Get microbiome data matrix
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Transpose data so that samples are rows and features are columns
microbiome_data_t <- t(microbiome_data)
metabolome_data_t <- t(metabolome_data)

# Ensure sample order matches
stopifnot(all(rownames(microbiome_data_t) == rownames(metabolome_data_t)))
stopifnot(all(rownames(microbiome_data_t) == rownames(metadata)))

# Create list of data matrices
data_list <- list(
  microbiome = microbiome_data_t,
  metabolome = metabolome_data_t
)

# Create design matrix for multi-omics integration
# For diagonal design (each dataset is connected to every other dataset)
design <- matrix(1, ncol = length(data_list), nrow = length(data_list),
                 dimnames = list(names(data_list), names(data_list)))
diag(design) <- 0

print("Data prepared for DIablo analysis!")
print(paste("Microbiome features:", ncol(microbiome_data_t)))
print(paste("Metabolite features:", ncol(metabolome_data_t)))
print(paste("Sample:", nrow(microbiome_data_t)))

# ------------------------------------------------------------------------------
# 3. Perform diablo analysis
# ------------------------------------------------------------------------------
cat("\n=== performing diablo analysis ===\n")

# Set number of component
ncomp <- 3

# Set number of variables to keep for each block on each component
# These values may need to be adjusted based on your data
keep <- list(
  microbiome = c(50, 30, 20),
  metabolome = c(50, 30, 20)
)

# perform diablo
set.seed(123)
diablo_result <- block.splsda(
  X = data_list,
  Y = metadata$Group,
  ncomp = ncomp,
  keepX = keep,
  design = design
)

# Save diablo results
saveRDS(diablo_result, "results/multi_omics_integration/DIABLO/diablo_result.rds")

print("diablo analysis completed successfully!")

# ------------------------------------------------------------------------------
# 4. Visualize diablo results
# ------------------------------------------------------------------------------
cat("\n=== Visualizing diablo results ===\n")

# 4.1 Sample plot
png("results/figures/multi_omics_integration/DIABLO/diablo_sample_plot.png", width = 1000, height = 800, res = 300)
plotIndiv(
  diablo_result,
  ind.names = FALSE,
  group = metadata$Group,
  legend = TRUE,
  title = "diablo Sample Plot",
  ellipse = TRUE
)
dev.off()

# 4.2 Correlation circle plot
png("results/figures/multi_omics_integration/DIABLO/diablo_correlation_circle.png", width = 1000, height = 800, res = 300)
plotVar(
  diablo_result,
  col = c("blue", "red"),
  legend = TRUE,
  title = "diablo Correlation Circle"
)
dev.off()

# 4.3 Network plot
png("results/figures/multi_omics_integration/DIABLO/diablo_network.png", width = 1200, height = 1000, res = 300)
plotNetwork(
  diablo_result,
  blocks = c(1, 2),
  comp = 1,
  cutoff = 0.7,
  color.blocks = c("blue", "red"),
  shape.blocks = c("circle", "square"),
  title = "diablo Network Plot"
)
dev.off()

# 4.4 Heatmap of block weights
png("results/figures/multi_omics_integration/DIABLO/diablo_block_weights.png", width = 1000, height = 800, res = 300)
plotLoadings(
  diablo_result,
  block = "microbiome",
  comp = 1,
  method = "mean",
  main = "Microbiome Block Weights - Component 1"
)
dev.off()

png("results/figures/multi_omics_integration/DIABLO/metabolome_block_weights.png", width = 1000, height = 800, res = 300)
plotLoadings(
  diablo_result,
  block = "metabolome",
  comp = 1,
  method = "mean",
  main = "Metabolome Block Weights - Component 1"
)
dev.off()

# 4.5 Circos plot
png("results/figures/multi_omics_integration/DIABLO/diablo_circos.png", width = 1200, height = 1200, res = 300)
circosPlot(
  diablo_result,
  cutoff = 0.7,
  color.blocks = c("blue", "red"),
  title = "diablo Circos Plot"
)
dev.off()

print("diablo visualization generated successfully!")

# ------------------------------------------------------------------------------
# 5. Perform permutation test for diablo
# ------------------------------------------------------------------------------
cat("\n=== performing permutation test for diablo ===\n")

set.seed(123)
diablo_perm <- perf(
  diablo_result,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)

# Save permutation test results
saveRDS(diablo_perm, "results/multi_omics_integration/DIABLO/diablo_permutation_test.rds")

# Plot permutation test results
png("results/figures/multi_omics_integration/DIABLO/diablo_permutation_test.png", width = 1000, height = 800, res = 300)
plot(diablo_perm, col = color.mixo(1:2), sd = TRUE, legend.position = "horizontal")
dev.off()

print("diablo permutation test completed successfully!")

# ------------------------------------------------------------------------------
# 6. Extract and analyze significant variables
# ------------------------------------------------------------------------------
cat("\n=== Extracting significant variables ===\n")

# Function to extract significant variables from diablo results
extract_significant_variables <- function(diablo_result, block_name, comp = 1, cutoff = 0.7) {
  
  # Get variable weights
  weights <- diablo_result$loadings[[block_name]][, comp]
  
  # Select variables above cutoff
  significant_vars <- names(weights)[abs(weights) > cutoff]
  
  # Create data frame of significant variables
  significant_df <- data.frame(
    Feature = significant_vars,
    Weight = weights[significant_vars],
    Block = block_name,
    Component = comp,
    stringsAsFactors = FALSE
  )
  
  return(significant_df)
}

# Extract significant variables for each block
microbiome_significant <- extract_significant_variables(
  diablo_result = diablo_result,
  block_name = "microbiome",
  comp = 1,
  cutoff = 0.7
)

metabolome_significant <- extract_significant_variables(
  diablo_result = diablo_result,
  block_name = "metabolome",
  comp = 1,
  cutoff = 0.7
)

# Combine significant variables
significant_variables <- bind_rows(microbiome_significant, metabolome_significant)

# Add annotations
taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
significant_variables <- significant_variables %>%
  left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature")

if (!is.null(metabolite_annotations)) {
  significant_variables <- significant_variables %>%
    left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature")
}

# Save significant variables
write.csv(significant_variables, "results/multi_omics_integration/DIABLO/diablo_significant_variables.csv", row.names = FALSE)

print("Significant variables extracted successfully!")
print(paste("Significant microbiome features:", nrow(microbiome_significant)))
print(paste("Significant metabolite features:", nrow(metabolome_significant)))

# ------------------------------------------------------------------------------
# 7. Visualize significant variables
# ------------------------------------------------------------------------------
cat("\n=== Visualizing significant variables ===\n")

# 7.1 Microbiome significant variables
if (nrow(microbiome_significant) > 0) {
  # Create taxonomic labels
  microbiome_significant <- microbiome_significant %>%
    left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature") %>%
    mutate(
      Taxonomic_Label = case_when(
        !is.na(Genus) ~ paste(Phylum, Genus, sep = "|"),
        !is.na(Family) ~ paste(Phylum, Family, sep = "|"),
        !is.na(Order) ~ paste(Phylum, Order, sep = "|"),
        !is.na(Class) ~ paste(Phylum, Class, sep = "|"),
        TRUE ~ Phylum
      )
    )
  
  # Order by absolute weight
  microbiome_significant <- microbiome_significant %>%
    arrange(desc(abs(Weight)))
  
  # Create lollipop plot
  p_microbiome_significant <- ggplot(microbiome_significant, aes(y = reorder(Taxonomic_Label, abs(Weight)), x = Weight)) +
    geom_segment(aes(x = 0, xend = Weight, y = Taxonomic_Label, yend = Taxonomic_Label), color = "gray50") +
    geom_point(color = "blue", size = 3) +
    labs(
      title = "Significant Microbiome Features (diablo)",
      x = "Weight",
      y = "Feature"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave("results/figures/multi_omics_integration/DIABLO/microbiome_significant_variables.png", p_microbiome_significant,
         width = 12, height = 10, dpi = 300)
}

# 7.2 Metabolome significant variables
if (nrow(metabolome_significant) > 0) {
  # Add metabolite names
  metabolome_significant <- metabolome_significant %>%
    left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature") %>%
    mutate(
      Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Feature)
    )
  
  # Order by absolute weight
  metabolome_significant <- metabolome_significant %>%
    arrange(desc(abs(Weight)))
  
  # Create lollipop plot
  p_metabolome_significant <- ggplot(metabolome_significant, aes(y = reorder(Metabolite_Label, abs(Weight)), x = Weight)) +
    geom_segment(aes(x = 0, xend = Weight, y = Metabolite_Label, yend = Metabolite_Label), color = "gray50") +
    geom_point(color = "red", size = 3) +
    labs(
      title = "Significant Metabolite Features (diablo)",
      x = "Weight",
      y = "Feature"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave("results/figures/multi_omics_integration/DIABLO/metabolome_significant_variables.png", p_metabolome_significant,
         width = 12, height = 10, dpi = 300)
}

# 7.3 Combined significant variables heatmap
if (nrow(microbiome_significant) > 0 && nrow(metabolome_significant) > 0) {
  # Get significant features
  significant_microbiome_features <- microbiome_significant$Feature
  significant_metabolome_features <- metabolome_significant$Feature
  
  # Subset data
  significant_microbiome_data <- microbiome_data_t[, significant_microbiome_features]
  significant_metabolome_data <- metabolome_data_t[, significant_metabolome_features]
  
  # Combine data
  significant_data <- cbind(significant_microbiome_data, significant_metabolome_data)
  
  # Add group information
  significant_data <- cbind(significant_data, Group = metadata$Group)
  
  # Create heatmap
  p_heatmap <- pheatmap(
    t(significant_data[, -ncol(significant_data)]),
    scale = "row",
    annotation_col = significant_data[, "Group", drop = FALSE],
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = "Significant Features from diablo Integration",
    color = viridis(100),
    filename = "results/figures/multi_omics_integration/DIABLO/significant_features_heatmap.png",
    width = 12,
    height = 10,
    dpi = 300
  )
}

print("Significant variables visualization generated successfully!")

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== DIABLO_integration.R Summary ===\n")
cat("\nGenerated results:\n")
cat("1. diablo multi-omics integration analysis\n")
cat("2. diablo sample plot\n")
cat("3. diablo correlation circle plot\n")
cat("4. diablo network plot\n")
cat("5. diablo block weights heatmap\n")
cat("6. diablo circos plot\n")
cat("7. diablo permutation test\n")
cat("8. Significant variable extraction and analysis\n")
cat("9. Significant variable visualization\n")

cat("\nResults saved to: results/multi_omics_integration/DIABLO/\n")
cat("Figures saved to: results/figures/multi_omics_integration/DIABLO/\n")
