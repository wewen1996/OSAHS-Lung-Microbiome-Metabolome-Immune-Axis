# metabolomics_OPLSDA.R
# Perform OPLS-DA analysis on metabolome data to identify differentially abundant metabolites

# Load required libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(ropls)
library(mixOmics)
library(caret)
library(pROC)
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
# Load normalized metabolome data (Pareto scaled is recommended for OPLS-DA)
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Metabolome data loaded successfully:")
print(dim(metabolome_data))
print(head(metabolome_data[, 1:5]))

# Check the grouping variable
print("Grouping variable levels:")
print(table(metadata$Group))

# ------------------------------------------------------------------------------
# 2. Prepare data for OPLS-DA
# ------------------------------------------------------------------------------
# Transpose metabolome data so that samples are rows and metabolites are columns
metabolome_data_t <- t(metabolome_data) %>% as.data.frame()

# Ensure samples are in the same order
stopifnot(all(rownames(metabolome_data_t) == rownames(metadata)))

# Extract grouping variable
group <- metadata$Group

# Check for two groups
groups <- unique(group)
if (length(groups) != 2) {
  stop("This script is designed for two-group comparison only.")
}

print(paste("Comparing", groups[1], "vs", groups[2]))
print(paste("Group 1 samples:", sum(group == groups[1])))
print(paste("Group 2 samples:", sum(group == groups[2])))

# Convert group to factor
group_factor <- factor(group, levels = groups)

# ------------------------------------------------------------------------------
# 3. Perform OPLS-DA analysis
# ------------------------------------------------------------------------------
cat("\nPerforming OPLS-DA analysis...\n")

# Create OPLS-DA model
oplsda_model <- opls(
  x = metabolome_data_t,
  y = group_factor,
  predI = 1,          # Number of predictive components
  orthoI = NA,        # Number of orthogonal components (NA = automatic selection)
  scaleC = "none",    # No scaling (data is already Pareto scaled)
  info.txtC = "none"  # No verbose output
)

# Print model summary
print(oplsda_model)

# Extract model statistics
r2x <- oplsda_model@summaryDF$R2X[1]          # R2X (explained variance in X)
r2y <- oplsda_model@summaryDF$R2Y[1]          # R2Y (explained variance in Y)
q2 <- oplsda_model@summaryDF$Q2[1]            # Q2 (predictive ability)

cat("\n=== OPLS-DA Model Statistics ===\n")
cat(paste("R2X:", round(r2x, 3), "\n"))
cat(paste("R2Y:", round(r2y, 3), "\n"))
cat(paste("Q2:", round(q2, 3), "\n"))

# ------------------------------------------------------------------------------
# 4. Model validation using permutation test
# ------------------------------------------------------------------------------
cat("\nPerforming permutation test for model validation...\n")

# Perform permutation test (100 permutations)
set.seed(123)
perm_results <- validate(oplsda_model, nperm = 100)

# Print permutation test results
print(perm_results)

# Extract p-values from permutation test
p_r2y <- perm_results$pR2Y
p_q2 <- perm_results$pQ2

cat("\n=== Permutation Test Results ===\n")
cat(paste("P-value (R2Y):", round(p_r2y, 4), "\n"))
cat(paste("P-value (Q2):", round(p_q2, 4), "\n"))

# ------------------------------------------------------------------------------
# 5. Identify significant metabolites
# ------------------------------------------------------------------------------
cat("\nIdentifying significant metabolites...\n")

# Extract variable importance in projection (VIP) scores
vip_scores <- oplsda_model@vipVn
vip_df <- data.frame(
  Metabolite = names(vip_scores),
  VIP = vip_scores,
  stringsAsFactors = FALSE
)

# Extract coefficients from the model
coefficients <- oplsda_model@coefMN
coef_df <- data.frame(
  Metabolite = names(coefficients),
  Coefficient = coefficients,
  stringsAsFactors = FALSE
)

# Combine VIP scores and coefficients
metabolite_results <- vip_df %>%
  left_join(coef_df, by = "Metabolite") %>%
  # Add p-values using t-test
  mutate(
    Group1_Mean = apply(metabolome_data[, group == groups[1]], 1, mean)[match(Metabolite, rownames(metabolome_data))],
    Group2_Mean = apply(metabolome_data[, group == groups[2]], 1, mean)[match(Metabolite, rownames(metabolome_data))],
    Fold_Change = Group1_Mean / Group2_Mean,
    Log2_Fold_Change = log2(Fold_Change)
  )

# Add p-values using t-test
for (i in 1:nrow(metabolite_results)) {
  metabolite <- metabolite_results$Metabolite[i]
  group1_values <- metabolome_data[metabolite, group == groups[1]]
  group2_values <- metabolome_data[metabolite, group == groups[2]]
  
  t_test_result <- t.test(group1_values, group2_values)
  metabolite_results$P_Value[i] <- t_test_result$p.value
}

# Add adjusted p-values (FDR)
metabolite_results$Adjusted_P_Value <- p.adjust(metabolite_results$P_Value, method = "fdr")

# Add significance level
metabolite_results$Significance <- case_when(
  metabolite_results$Adjusted_P_Value < 0.001 ~ "***",
  metabolite_results$Adjusted_P_Value < 0.01 ~ "**",
  metabolite_results$Adjusted_P_Value < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Define significant metabolites (VIP > 1 and FDR < 0.05)
metabolite_results$Significant <- (metabolite_results$VIP > 1) & (metabolite_results$Adjusted_P_Value < 0.05)

# Add metabolite annotations
if (!is.null(metabolite_annotations)) {
  metabolite_results <- metabolite_results %>%
    left_join(metabolite_annotations, by = c("Metabolite" = rownames(metabolite_annotations)))
}

# Order results by adjusted p-value
metabolite_results <- metabolite_results %>% arrange(Adjusted_P_Value)

# Print summary of results
cat("\n=== OPLS-DA Results Summary ===\n")
print(paste("Total metabolites tested:", nrow(metabolite_results)))
print(paste("Significant metabolites (VIP > 1, FDR < 0.05):", sum(metabolite_results$Significant)))

# ------------------------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------------------------
# Save full results
write.csv(metabolite_results, "results/tables/differential_analysis/metabolomics_oplsda_results.csv", row.names = FALSE)

# Save significant results
significant_metabolites <- metabolite_results %>% filter(Significant)
write.csv(significant_metabolites, "results/tables/differential_analysis/metabolomics_oplsda_significant.csv", row.names = FALSE)

print(paste("Significant results saved. Found", nrow(significant_metabolites), "significant metabolites."))

# ------------------------------------------------------------------------------
# 7. Generate visualization
# ------------------------------------------------------------------------------
# 7.1 OPLS-DA score plot
cat("\nGenerating OPLS-DA score plot...\n")

# Extract scores values
score_data <- oplsda_model@scoreMN %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(Group = group[match(Sample, names(group))])

# Create score plot
p_score <- ggplot(score_data, aes(x = p1, y = o1, color = Group, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "OPLS-DA Score Plot",
       x = paste("Predictive Component 1 (", round(oplsda_model@summaryDF$R2X[1] * 100, 1), "%)", sep = ""),
       y = paste("Orthogonal Component 1 (", round(oplsda_model@summaryDF$R2X[2] * 100, 1), "%)", sep = "")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  stat_ellipse(type = "t", level = 0.95)

# Add model statistics to plot
p_score <- p_score +
  annotate("text", x = Inf, y = Inf, 
           label = paste("R2X =", round(r2x, 3), "\nR2Y =", round(r2y, 3), "\nQ2 =", round(q2, 3)),
           hjust = 1.1, vjust = 1.1, size = 3.5)

# Save score plot
ggsave("results/figures/differential_analysis/metabolomics_oplsda_score_plot.png", p_score,
       width = 10, height = 8, dpi = 300)

# 7.2 VIP plot
cat("\nGenerating VIP plot...\n")

# Get top 20 metabolites by VIP score
top_vip_metabolites <- metabolite_results %>%
  arrange(desc(VIP)) %>%
  head(20) %>%
  mutate(Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Metabolite)) %>%
  mutate(Metabolite_Label = factor(Metabolite_Label, levels = rev(Metabolite_Label)))

# Create VIP plot
p_vip <- ggplot(top_vip_metabolites, aes(x = VIP, y = Metabolite_Label, fill = Significant)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Top 20 Variables Importance in Projection (VIP) Scores",
       x = "VIP Score",
       y = "Metabolite") +
  scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save VIP plot
ggsave("results/figures/differential_analysis/metabolomics_oplsda_vip_plot.png", p_vip,
       width = 12, height = 10, dpi = 300)

# 7.3 Volcano plot
cat("\nGenerating volcano plot...\n")

# Prepare data for volcano plot
volcano_data <- metabolite_results %>%
  mutate(
    Neg_Log10_P_Adj = -log10(Adjusted_P_Value),
    Significance = case_when(
      VIP > 1 & Adjusted_P_Value < 0.05 ~ "Significant (VIP > 1, FDR < 0.05)",
      VIP > 1 ~ "VIP > 1",
      Adjusted_P_Value < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    )
  )

# Create volcano plot
p_volcano <- ggplot(volcano_data, aes(x = Log2_Fold_Change, y = Neg_Log10_P_Adj, color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Volcano Plot of OPLS-DA Results",
       x = paste("Log2 Fold Change (", groups[1], " vs ", groups[2], ")", sep = ""),
       y = "-Log10(Adjusted P-value)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("black", "blue", "orange", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Add labels for top significant metabolites
top_significant <- volcano_data %>%
  filter(Significance == "Significant (VIP > 1, FDR < 0.05)") %>%
  arrange(desc(Neg_Log10_P_Adj)) %>%
  head(10) %>%
  mutate(Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Metabolite))

if (nrow(top_significant) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = top_significant,
      aes(label = Metabolite_Label),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.5,
      segment.color = "gray50"
    )
}

# Save volcano plot
ggsave("results/figures/differential_analysis/metabolomics_oplsda_volcano_plot.png", p_volcano,
       width = 12, height = 8, dpi = 300)

# 7.4 Heatmap of significant metabolites
cat("\nGenerating heatmap of significant metabolites...\n")

if (nrow(significant_metabolites) > 0) {
  # Get top N significant metabolites (up to 50)
  top_n <- min(50, nrow(significant_metabolites))
  top_metabolites <- significant_metabolites$Metabolite[1:top_n]
  
  # Subset metabolome data to top metabolites
  top_metabolite_data <- metabolome_data[top_metabolites, ]
  
  # Add metabolite names to row names if available
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    row_labels <- apply(metabolite_annotations[top_metabolites, c("Metabolite.Name", "Metabolite")], 1, function(x) {
      if (!is.na(x[1]) && x[1] != "") x[1] else x[2]
    })
    rownames(top_metabolite_data) <- row_labels
  }
  
  # Create annotation for samples groups
  sample_annot <- data.frame(Group = group)
  rownames(sample_annot) <- colnames(top_metabolite_data)
  
  # Create heatmap
  pheatmap(
    top_metabolite_data,
    scale = "row",
    annotation_col = sample_annot,
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = paste("Top", top_n, "Significant Metabolites (VIP > 1, FDR < 0.05)"),
    filename = "results/figures/differential_analysis/metabolomics_heatmap_significant.png",
    width = 10,
    height = 10,
    dpi = 300
  )
  
  print(paste("Heatmap generated with top", top_n, "significant metabolites."))
} else {
  print("No significant metabolites found. Skipping heatmap generation.")
}

# 7.5 Boxplots for top significant metabolites
cat("\nGenerating boxplots for top significant metabolites...\n")

# Get top 5 significant metabolites
top_5_metabolites <- significant_metabolites$Metabolite[1:min(5, nrow(significant_metabolites))]

if (length(top_5_metabolites) > 0) {
  # Prepare data for boxplots
  boxplot_data <- as.data.frame(t(metabolome_data[top_5_metabolites, ])) %>%
    rownames_to_column("Sample") %>%
    mutate(Group = group[match(Sample, names(group))]) %>%
    pivot_longer(cols = -c(Sample, Group), names_to = "Metabolite", values_to = "Abundance")
  
  # Add metabolite names if available
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    boxplot_data$Metabolite_Label <- apply(metabolite_annotations[boxplot_data$Metabolite, c("Metabolite.Name", "Metabolite")], 1, function(x) {
      if (!is.na(x[1]) && x[1] != "") x[1] else x[2]
    })
  } else {
    boxplot_data$Metabolite_Label <- boxplot_data$Metabolite
  }
  
  # Create boxplots
  p_boxplot <- ggplot(boxplot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_wrap(~Metabolite_Label, scales = "free_y") +
    labs(title = "Top Significant Metabolites Abundance by Group",
         x = "Group",
         y = "Normalized Abundance (Pareto Scaled)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # Save boxplots
  ggsave("results/figures/differential_analysis/metabolomics_boxplots_top_significant.png", p_boxplot,
         width = 15, height = 10, dpi = 300)
  
  print("Boxplots generated for top significant metabolites.")
} else {
  print("No significant metabolites found. Skipping boxplot generation.")
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== Metabolomics_OPLSDA.R Summary ===\n")
print(paste("Comparison:", groups[1], "vs", groups[2]))
print(paste("Total metabolites tested:", nrow(metabolite_results)))
print(paste("Significant metabolites (VIP > 1, FDR < 0.05):", sum(metabolite_results$Significant)))

cat("\n=== OPLS-DA Model Performance ===\n")
cat(paste("R2X:", round(r2x, 3), "(Explained variance in X)\n"))
cat(paste("R2Y:", round(r2y, 3), "(Explained variance in Y)\n"))
cat(paste("Q2:", round(q2, 3), "(Predictive ability)\n"))
cat(paste("Permutation test p-value (R2Y):", round(p_r2y, 4), "\n"))
cat(paste("Permutation test p-value (Q2):", round(p_q2, 4), "\n"))

if (sum(metabolite_results$Significant) > 0) {
  print("\nTop 5 significant metabolites:")
  print(significant_metabolites %>%
          select(Metabolite, Metabolite.Name, VIP, Log2_Fold_Change, Adjusted_P_Value) %>%
          head(5))
}

print("\nFiles generated:")
print("1. results/tables/differential_analysis/metabolomics_oplsda_results.csv")
print("2. results/tables/differential_analysis/metabolomics_oplsda_significant.csv")
print("3. results/figures/differential_analysis/metabolomics_oplsda_score_plot.png")
print("4. results/figures/differential_analysis/metabolomics_oplsda_vip_plot.png")
print("5. results/figures/differential_analysis/metabolomics_oplsda_volcano_plot.png")
if (nrow(significant_metabolites) > 0) {
  print("6. results/figures/differential_analysis/metabolomics_heatmap_significant.png")
  if (length(top_5_metabolites) > 0) {
    print("7. results/figures/differential_analysis/metabolomics_boxplots_top_significant.png")
  }
}
