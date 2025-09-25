# data_normalization.R
# Normalize microbiome and metabolome data for downstream analysis

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(metagenomeSeq)
library(preprocessCore)
library(impute)
library(DESeq2)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("data/processed/normalized", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load processed data
# ------------------------------------------------------------------------------
# Load phyloseq object
phyloseq_obj <- readRDS("data/processed/phyloseq_object.rds")
print("Phyloseq object loaded successfully:")
print(phyloseq_obj)

# Load metabolome data
metabolome_data <- read.csv("data/processed/metabolome_data_processed.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Metabolome data loaded successfully:")
print(dim(metabolome_data))

# ------------------------------------------------------------------------------
# 2. Normalize microbiome data
# ------------------------------------------------------------------------------
cat("\n=== Normalizing Microbiome Data ===\n")

# Extract OTU table from phyloseq object
otu_table_raw <- otu_table(phyloseq_obj) %>% as.matrix()

# a) Total sum scaling (TSS) normalization
otu_tss <- t(t(otu_table_raw) / colSums(otu_table_raw))
print("TSS normalization completed")

# b) Cumulative sum scaling (CSS) normalization using metagenomeSeq
# Create MRexperiment object
mrexp <- newMRexperiment(otu_table_raw)

# Apply CSS normalization
mrexp_css <- cumNorm(mrexp, p = 0.5)
otu_css <- t(scale(t(MRcounts(mrexp_css, norm = TRUE)), center = FALSE, scale = TRUE))
print("CSS normalization completed")

# c) DESeq2 variance stabilizing transformation (VST)
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = otu_table_raw,
  colData = metadata,
  design = ~ 1
)

# Apply VST normalization
dds_vst <- varianceStabilizingTransformation(dds)
otu_vst <- assay(dds_vst)
print("VST normalization completed")

# d) Relative abundance (percentage)
otu_relative <- otu_tss * 100
print("Relative abundance calculation completed")

# e) Log transformation (after adding a pseudocount)
otu_log <- log10(otu_table_raw + 1)
print("Log transformation completed")

# ------------------------------------------------------------------------------
# 3. Normalize metabolome data
# ------------------------------------------------------------------------------
cat("\n=== Normalizing Metabolome Data ===\n")

# a) Remove missing values (impute if necessary)
# Check for missing values
metabolome_missing <- sum(is.na(metabolome_data))
print(paste("Missing values in metabolome data:", metabolome_missing))

if (metabolome_missing > 0) {
  # Impute missing values using k-nearest neighbors
  metabolome_imputed <- impute.knn(as.matrix(metabolome_data))$data
  print("Missing values imputed using k-nearest neighbors method")
} else {
  metabolome_imputed <- metabolome_data
}

# b) Log transformation
metabolome_log <- log2(metabolome_imputed)
print("Log2 transformation completed")

# c) Pareto scaling (mean centering and scaling by square root of standard deviation)
metabolome_pareto <- preprocessCore::normalize.quantiles.use.target(
  metabolome_log,
  target = apply(metabolome_log, 1, function(x) {
    (x - mean(x)) / sqrt(sd(x))
  })
)
rownames(metabolome_pareto) <- rownames(metabolome_log)
colnames(metabolome_pareto) <- colnames(metabolome_log)
print("Pareto scaling completed")

# d) Autoscaling (mean centering and unit variance scaling)
metabolome_auto <- preprocessCore::normalize.quantiles.use.target(
  metabolome_log,
  target = apply(metabolome_log, 1, function(x) {
    (x - mean(x)) / sd(x)
  })
)
rownames(metabolome_auto) <- rownames(metabolome_log)
colnames(metabolome_auto) <- colnames(metabolome_log)
print("Autoscaling completed")

# e) Quantile normalization
metabolome_quantile <- preprocessCore::normalize.quantiles(as.matrix(metabolome_log))
rownames(metabolome_quantile) <- rownames(metabolome_log)
colnames(metabolome_quantile) <- colnames(metabolome_log)
print("Quantile normalization completed")

# ------------------------------------------------------------------------------
# 4. Save normalized data
# ------------------------------------------------------------------------------
# Save normalized microbiome data
write.csv(otu_tss, "data/processed/normalized/microbiome_tss_normalized.csv", row.names = TRUE)
write.csv(otu_css, "data/processed/normalized/microbiome_css_normalized.csv", row.names = TRUE)
write.csv(otu_vst, "data/processed/normalized/microbiome_vst_normalized.csv", row.names = TRUE)
write.csv(otu_relative, "data/processed/normalized/microbiome_relative_abundance.csv", row.names = TRUE)
write.csv(otu_log, "data/processed/normalized/microbiome_log_transformed.csv", row.names = TRUE)

# Save normalized metabolome data
write.csv(metabolome_imputed, "data/processed/normalized/metabolome_imputed.csv", row.names = TRUE)
write.csv(metabolome_log, "data/processed/normalized/metabolome_log_transformed.csv", row.names = TRUE)
write.csv(metabolome_pareto, "data/processed/normalized/metabolome_pareto_scaled.csv", row.names = TRUE)
write.csv(metabolome_auto, "data/processed/normalized/metabolome_autoscaled.csv", row.names = TRUE)
write.csv(metabolome_quantile, "data/processed/normalized/metabolome_quantile_normalized.csv", row.names = TRUE)

# Save normalized phyloseq objects
# Create phyloseq object with TSS normalized data
phyloseq_tss <- phyloseq_obj
otu_table(phyloseq_tss) <- otu_table(otu_tss, taxa_are_rows = TRUE)
saveRDS(phyloseq_tss, "data/processed/normalized/phyloseq_tss.rds")

# Create phyloseq object with CSS normalized data
phyloseq_css <- phyloseq_obj
otu_table(phyloseq_css) <- otu_table(otu_css, taxa_are_rows = TRUE)
saveRDS(phyloseq_css, "data/processed/normalized/phyloseq_css.rds")

print("All normalized data saved successfully!")

# ------------------------------------------------------------------------------
# 5. Generate normalization quality control plots
# ------------------------------------------------------------------------------
# Create output directory for plots
dir.create("results/figures/normalization", showWarnings = FALSE, recursive = TRUE)

# a) Microbiome normalization comparison (boxplots of sample distributions)
# Prepare data for plotting
microbiome_distributions <- data.frame(
  Raw = as.vector(otu_table_raw),
  TSS = as.vector(otu_tss),
  CSS = as.vector(otu_css),
  VST = as.vector(otu_vst),
  Log = as.vector(otu_log)
)

# Melt data for ggplot2
microbiome_distributions_melt <- reshape2::melt(microbiome_distributions)

# Create boxplot
p_microbiome_norm <- ggplot(microbiome_distributions_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(title = "Microbiome Data Normalization Comparison",
       x = "Normalization Method", y = "Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/normalization/microbiome_normalization_comparison.png", p_microbiome_norm,
       width = 12, height = 8, dpi = 300)

# b) Metabolome normalization comparison
# Prepare data for plotting
metabolome_distributions <- data.frame(
  Raw = as.vector(metabolome_data),
  Imputed = as.vector(metabolome_imputed),
  Log = as.vector(metabolome_log),
  Pareto = as.vector(metabolome_pareto),
  Auto = as.vector(metabolome_auto),
  Quantile = as.vector(metabolome_quantile)
)

# Melt data for ggplot2
metabolome_distributions_melt <- reshape2::melt(metabolome_distributions)

# Create boxplot
p_metabolome_norm <- ggplot(metabolome_distributions_melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(title = "Metabolome Data Normalization Comparison",
       x = "Normalization Method", y = "Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/normalization/metabolome_normalization_comparison.png", p_metabolome_norm,
       width = 14, height = 8, dpi = 300)

print("Normalization quality control plots generated successfully!")

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== Data Normalization Summary ===\n")
cat("Microbiome normalization methods applied:\n")
cat("  - Total Sum Scaling (TSS)\n")
cat("  - Cumulative Sum Scaling (CSS)\n")
cat("  - Variance Stabilizing Transformation (VST)\n")
cat("  - Relative Abundance (percentage)\n")
cat("  - Log transformation (log10(count + 1))\n")

cat("\nMetabolome normalization methods applied:\n")
cat("  - Missing value imputation (k-nearest neighbors)\n")
cat("  - Log2 transformation\n")
cat("  - Pareto scaling\n")
cat("  - Autoscaling (mean centering and unit variance)\n")
cat("  - Quantile normalization\n")

cat("\nNormalized data saved to: data/processed/normalized/\n")
cat("Normalization quality control plotsssaved to: results/figures/normalization/\n")
