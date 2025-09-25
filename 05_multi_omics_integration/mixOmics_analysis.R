# mixOmics_analysis.R
# Perform multi-omics integration using mixOmics package

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
dir.create("results/multi_omics_integration", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/multi_omics_integration", showWarnings = FALSE, recursive = TRUE)

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
# 2. Prepare data for mixOmics analysis
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

print("Data prepared for mixOmics analysis!")
print(paste("Microbiome features:", ncol(microbiome_data_t)))
print(paste("Metabolite features:", ncol(metabolome_data_t)))
print(paste("Samples:", nrow(microbiome_data_t)))

# ------------------------------------------------------------------------------
# 3. Perform PCA on individual datasets
# ------------------------------------------------------------------------------
cat("\n=== Performing PCA on individual datasets ===\n")

# 3.1 Microbiome PCA
microbiome_pca <- pca(
  X = microbiome_data_t,
  ncomp = 5,
  center = TRUE,
  scale = FALSE
)

# Plot microbiome PCA
png("results/figures/multi_omics_integration/microbiome_pca.png", width = 1000, height = 800, res = 300)
plotIndiv(
  microbiome_pca,
  group = metadata$Group,
  ind.names = FALSE,
  legend = TRUE,
  title = "Microbiome PCA",
  ellipse = TRUE
)
dev.off()

# 3.2 Metabolome PCA
metabolome_pca <- pca(
  X = metabolome_data_t,
  ncomp = 5,
  center = TRUE,
  scale = FALSE
)

# Plot metabolome PCA
png("results/figures/multi_omics_integration/metabolome_pca.png", width = 1000, height = 800, res = 300)
plotIndiv(
  metabolome_pca,
  group = metadata$Group,
  ind.names = FALSE,
  legend = TRUE,
  title = "Metabolome PCA",
  ellipse = TRUE
)
dev.off()

print("PCA analysis completed successfully!")

# ------------------------------------------------------------------------------
# 4. Perform sPLS analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing sPLS analysis ===\n")

# Perform sPLS between microbiome and metabolome
spls_result <- spls(
  X = microbiome_data_t,
  Y = metabolome_data_t,
  ncomp = 3,
  keepX = c(50, 30, 20),  # Number of variables to keep for X (microbiome) on each component
  keepY = c(50, 30, 20),  # Number of variables to keep for Y (metabolome) on each component
  mode = "regression"
)

# Save sPLS results
saveRDS(spls_result, "results/multi_omics_integration/spls_result.rds")

# Plot sPLS sample plot
png("results/figures/multi_omics_integration/spls_sample_plot.png", width = 1000, height = 800, res = 300)
plotIndiv(
  spls_result,
  group = metadata$Group,
  ind.names = FALSE,
  legend = TRUE,
  title = "sPLS Sample Plot",
  ellipse = TRUE
)
dev.off()

# Plot sPLS correlation circle
png("results/figures/multi_omics_integration/spls_correlation_circle.png", width = 1000, height = 800, res = 300)
plotVar(
  spls_result,
  col = c("blue", "red"),
  cex = 0.7,
  legend = TRUE,
  title = "sPLS Correlation Circle"
)
dev.off()

# Plot sPLS loadings
png("results/figures/multi_omics_integration/spls_loadings.png", width = 1200, height = 1000, res = 300)
plotLoadings(
  spls_result,
  comp = 1,
  method = "mean",
  main = "sPLS Loadings - Component 1"
)
dev.off()

# Extract significant variables
spls_significant <- selectVar(spls_result, comp = 1)

# Create data frame of significant variables
microbiome_significant <- data.frame(
  Feature = names(spls_significant$X$value),
  Loading = spls_significant$X$value,
  Dataset = "Microbiome",
  stringsAsFactors = FALSE
)

metabolome_significant <- data.frame(
  Feature = names(spls_significant$Y$value),
  Loading = spls_significant$Y$value,
  Dataset = "Metabolome",
  stringsAsFactors = FALSE
)

spls_significant_df <- bind_rows(microbiome_significant, metabolome_significant)

# Add annotations
taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
spls_significant_df <- spls_significant_df %>%
  left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature")

if (!is.null(metabolite_annotations)) {
  spls_significant_df <- spls_significant_df %>%
    left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature")
}

# Save significant variables
write.csv(spls_significant_df, "results/multi_omics_integration/spls_significant_variables.csv", row.names = FALSE)

print("sPLS analysis completed successfully!")
print(paste("Significant microbiome features:", nrow(microbiome_significant)))
print(paste("Significant metabolome features:", nrow(metabolome_significant)))

# ------------------------------------------------------------------------------
# 5. Perform multilevel sPLS analysis (if there are repeated measures)
# ------------------------------------------------------------------------------
if ("PatientID" %in% colnames(metadata)) {
  cat("\n=== Performing multilevel sPLS analysis ===\n")
  
  # Create sample metadata with PatientID as the repeated measure
  sample_metadata <- data.frame(
    Sample = rownames(metadata),
    PatientID = metadata$PatientID,
    Group = metadata$Group
  )
  
  # Perform multilevel sPLS
  ml_spls_result <- spls(
    X = microbiome_data_t,
    Y = metabolome_data_t,
    ncomp = 3,
    keepX = c(50, 30, 20),
    keepY = c(50, 30, 20),
    mode = "regression",
    multilevel = sample_metadata[, "PatientID"]
  )
  
  # Save multilevel sPLS results
  saveRDS(ml_spls_result, "results/multi_omics_integration/ml_spls_result.rds")
  
  # Plot multilevel sPLS sample plot
  png("results/figures/multi_omics_integration/ml_spls_sample_plot.png", width = 1000, height = 800, res = 300)
  plotIndiv(
    ml_spls_result,
    group = metadata$Group,
    ind.names = FALSE,
    legend = TRUE,
    title = "Multilevel sPLS Sample Plot",
    ellipse = TRUE
  )
  dev.off()
  
  print("Multilevel sPLS analysis completed successfully!")
} else {
  print("No PatientID found in metadata. Skipping multilevel sPLS analysis.")
}

# ------------------------------------------------------------------------------
# 6. Perform PLS-DA analysis
# ------------------------------------------------------------------------------
cat("\n=== Performing PLS-DA analysis ===\n")

# Combine data for PLS-DA
combined_data <- cbind(microbiome_data_t, metabolome_data_t)

# Perform PLS-DA
plsda_result <- plsda(
  X = combined_data,
  Y = metadata$Group,
  ncomp = 5
)

# Save PLS-DA results
saveRDS(plsda_result, "results/multi_omics_integration/plsda_result.rds")

# Plot PLS-DA sample plot
png("results/figures/multi_omics_integration/plsda_sample_plot.png", width = 1000, height = 800, res = 300)
plotIndiv(
  plsda_result,
  group = metadata$Group,
  ind.names = FALSE,
  legend = TRUE,
  title = "PLS-DA Sample Plot",
  ellipse = TRUE
)
dev.off()

# Perform permutation test for PLS-DA
set.seed(123)
plsda_perm <- perf(
  plsda_result,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  progressBar = TRUE
)

# Save permutation test results
saveRDS(plsda_perm, "results/multi_omics_integration/plsda_permutation_test.rds")

# Plot permutation test results
png("results/figures/multi_omics_integration/plsda_permutation_test.png", width = 1000, height = 800, res = 300)
plot(plsda_perm, col = color.mixo(1:2), sd = TRUE, legend.position = "horizontal")
dev.off()

# Extract variable importance
vip <- vip(plsda_result, plot = FALSE)

# Create VIP plot
vip_df <- data.frame(
  Feature = rownames(vip),
  VIP = vip[, 1],
  Dataset = ifelse(rownames(vip) %in% colnames(microbiome_data_t), "Microbiome", "Metabolome"),
  stringsAsFactors = FALSE
)

# Add annotations
vip_df <- vip_df %>%
  left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature")

if (!is.null(metabolite_annotations)) {
  vip_df <- vip_df %>%
    left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature")
}

# Order by VIP score
vip_df <- vip_df %>% arrange(desc(VIP))

# Save VIP results
write.csv(vip_df, "results/multi_omics_integration/plsda_vip_scores.csv", row.names = FALSE)

# Plot top VIP scores
top_vip <- vip_df %>% head(30)

p_vip <- ggplot(top_vip, aes(x = VIP, y = reorder(Feature, VIP), fill = Dataset)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Top 30 Variables Importance in Projection (VIP)",
    x = "VIP Score",
    y = "Feature"
  ) +
  scale_fill_manual(values = c("Microbiome" = "blue", "Metabolome" = "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave("results/figures/multi_omics_integration/plsda_vip_plot.png", p_vip,
       width = 12, height = 10, dpi = 300)

print("PLS-DA analysis completed successfully!")

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------
cat("\n=== mixOmics_analysis.R Summary ===\n")
cat("\nGenerated results:\n")
cat("1. PCA analysis for microbiome and metabolome\n")
cat("2. sPLS analysis between microbiome and metabolome\n")
if ("PatientID" %in% colnames(metadata)) {
  cat("3. Multilevel sPLS analysis (accounting for repeated measures)\n")
}
cat("4. PLS-DA analysis on combined data\n")
cat("5. Permutation test for PLS-DA\n")
cat("6. Variable Importance in Projection (VIP) analysis\n")

cat("\nResults saved to: results/multi_omics_integration/\n")
cat("Figures saved to: results/figures/multi_omics_integration/\n")
