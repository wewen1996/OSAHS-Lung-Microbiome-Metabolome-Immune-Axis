# age_adjusted_analysis.R
# Perform age-adjusted differential analysis on microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)
library(broom)
library(ggpubr)
library(gridExtra)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/tables/differential_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/differential_analysis", showWarnings = FALSE, recursive = TRUE)

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
print(paste("Microbiome samples:", nsamples(phyloseq_css)))
print(paste("Metabolome samples:", ncol(metabolome_data)))
print(paste("Metadata samples:", nrow(metadata)))

# Check if age variable exists
if (!"Age" %in% colnames(metadata)) {
  stop("Age variable not found in metadata. Please ensure metadata contains an 'Age' column.")
}

# Check if group variable exists
if (!"Group" %in% colnames(metadata)) {
  stop("Group variable not found in metadata. Please ensure metadata contains a 'Group' column.")
}

# ------------------------------------------------------------------------------
# 2. Prepare data for age-adjusted analysis
# ------------------------------------------------------------------------------
# Ensure samples IDs match across datasets
common_samples <- Reduce(intersect, list(
  sample_names(phyloseq_css),
  colnames(metabolome_data),
  rownames(metadata)
))

print(paste("Common samples IDs across all datasets:", length(common_samples)))

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

# Check for missing values in age
if (sum(is.na(metadata$Age)) > 0) {
  print(paste("Removing", sum(is.na(metadata$Age)), "samples with missing age values"))
  metadata <- metadata[!is.na(metadata$Age), ]
  phyloseq_css <- prune_samples(rownames(metadata), phyloseq_css)
  metabolome_data <- metabolome_data[, rownames(metadata)]
}

print("Final sample counts after quality control:")
print(table(metadata$Group))

# ------------------------------------------------------------------------------
# 3. Age-adjusted differential analysis for microbiome data
# ------------------------------------------------------------------------------
cat("\n=== Performing age-adjusted differential analysis for microbiome data ===\n")

# Extract OTU table and metadata from phyloseq object
otu_table <- otu_table(phyloseq_css) %>% as.matrix()
taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
group <- sample_data(phyloseq_css)$Group
age <- sample_data(phyloseq_css)$Age

# Ensure samples order matches
stopifnot(all(colnames(otu_table) == names(group)))
stopifnot(all(colnames(otu_table) == names(age)))

# Convert data to long format for mixed effects model
microbiome_long <- as.data.frame(t(otu_table)) %>%
  rownames_to_column("Sample") %>%
  mutate(Group = group[match(Sample, names(group))],
         Age = age[match(Sample, names(age))]) %>%
  pivot_longer(cols = -c(Sample, Group, age), names_to = "OTU", values_to = "Abundance")

# Function to perform age-adjusted differential analysis for each taxon
perform_age_adjusted_microbiome_analysis <- function(otu_table, group, age, taxonomy) {
  results <- data.frame(
    OTU = rownames(otu_table),
    Estimate = NA,
    Std.Error = NA,
    DF = NA,
    t.value = NA,
    P.Value = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(otu_table)) {
    otu_id <- rownames(otu_table)[i]
    otu_data <- otu_table[i, ]
    
    # Create data frame for the model
    model_data <- data.frame(
      Abundance = as.numeric(otu_data),
      Group = factor(group, levels = c("Control", "OSAHS")),
      Age = age
    )
    
    # Fit linear model adjusted for age
    model <- lm(Abundance ~ Group + Age, data = model_data)
    
    # Extract coefficients
    coef_summary <- summary(model)$coefficients["GroupOSAHS", ]
    
    results$Estimate[i] <- coef_summary["Estimate"]
    results$Std.Error[i] <- coef_summary["Std. Error"]
    results$DF[i] <- summary(model)$df[2]
    results$t.value[i] <- coef_summary["t value"]
    results$P.value[i] <- coef_summary["Pr(>|t|)"]
  }
  
  # Add adjusted p-values (FDR)
  results$Adjusted_P_Value <- p.adjust(results$P.value, method = "fdr")
  
  # Add significance level
  results$Significance <- case_when(
    results$Adjusted_P_Value < 0.001 ~ "***",
    results$Adjusted_P_Value < 0.01 ~ "**",
    results$Adjusted_P_Value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # Add taxonomy information
  results <- cbind(results, taxonomy[match(results$OTU, rownames(taxonomy)), ])
  
  # Order by adjusted p-value
  results <- results %>% arrange(Adjusted_P_Value)
  
  return(results)
}

# Perform age-adjusted analysis
microbiome_age_adjusted_results <- perform_age_adjusted_microbiome_analysis(otu_table, group, age, taxonomy)

# Print summary of results
cat("\n=== Microbiome Age-adjusted Results Summary ===\n")
print(paste("Total OTUs tested:", nrow(microbiome_age_adjusted_results)))
print(paste("Significant OTUs (FDR < 0.05):", sum(microbiome_age_adjusted_results$Adjusted_P_Value < 0.05)))

# ------------------------------------------------------------------------------
# 4. age-adjusted differential analysis for metabolome data
# ------------------------------------------------------------------------------
cat("\n=== Performing age-adjusted differential analysis for metabolome data ===\n")

# Transpose metabolome data so samples are rows and metabolites are columns
metabolome_data_t <- t(metabolome_data) %>% as.data.frame()

# Ensure sample order match
stopifnot(all(rownames(metabolome_data_t) == rownames(metadata)))

# Function to perform age-adjusted differential analysis for each metabolite
perform_age_adjusted_metabolome_analysis <- function(metabolome_data, group, age, metabolite_annotations = NULL) {
  results <- data.frame(
    Metabolite = colnames(metabolome_data),
    Estimate = NA,
    Std.Error = NA,
    DF = NA,
    t.value = NA,
    P.value = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in 1:ncol(metabolome_data)) {
    metabolite_id <- colnames(metabolome_data)[i]
    metabolite_data <- metabolome_data[, i]
    
    # Create data frame for the model
    model_data <- data.frame(
      Abundance = as.numeric(metabolite_data),
      Group = factor(group, level = c("Control", "OSAHS")),
      Age = age
    )
    
    # Fit linear model adjusted for age
    model <- lm(Abundance ~ Group + age, data = model_data)
    
    # Extract coefficients
    coef_summary <- summary(model)$coefficients["GroupOSAHS", ]
    
    results$Estimate[i] <- coef_summary["Estimate"]
    results$Std.Error[i] <- coef_summary["Std. Error"]
    results$DF[i] <- summary(model)$df[2]
    results$t.value[i] <- coef_summary["t value"]
    results$P.value[i] <- coef_summary["Pr(>|t|)"]
  }
  
  # Add adjusted p-values (FDR)
  results$Adjusted_P_Value <- p.adjust(results$P.value, method = "fdr")
  
  # Add significance level
  results$Significance <- case_when(
    results$Adjusted_P_Value < 0.001 ~ "***",
    results$Adjusted_P_Value < 0.01 ~ "**",
    results$Adjusted_P_Value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # Add metabolite annotations if available
  if (!is.null(metabolite_annotations)) {
    results <- results %>%
      left_join(metabolite_annotations, by = c("Metabolite" = rownames(metabolite_annotations)))
  }
  
  # Order by adjusted p-value
  results <- results %>% arrange(Adjusted_P_Value)
  
  return(results)
}

# Perform age-adjusted analysis
metabolome_age_adjusted_results <- perform_age_adjusted_metabolome_analysis(
  metabolome_data_t,
  metadata$Group,
  metadata$Age,
  metabolite_annotations
)

# Print summary of results
cat("\n=== Metabolome age-adjusted Results Summary ===\n")
print(paste("Total metabolites tested:", nrow(metabolome_age_adjusted_results)))
print(paste("Significant metabolites (FDR < 0.05):", sum(metabolome_age_adjusted_results$Adjusted_P_Value < 0.05)))

# ------------------------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------------------------
# Save microbiome results
write.csv(microbiome_age_adjusted_results, "results/tables/differential_analysis/microbiome_age_adjusted_results.csv", row.names = FALSE)

# Save significant microbiome results
microbiome_significant <- microbiome_age_adjusted_results %>% filter(Adjusted_P_Value < 0.05)
write.csv(microbiome_significant, "results/tables/differential_analysis/microbiome_age_adjusted_significant.csv", row.names = FALSE)

# Save metabolome results
write.csv(metabolome_age_adjusted_results, "results/tables/differential_analysis/metabolome_age_adjusted_results.csv", row.names = FALSE)

# Save significant metabolome results
metabolome_significant <- metabolome_age_adjusted_results %>% filter(Adjusted_P_Value < 0.05)
write.csv(metabolome_significant, "results/tables/differential_analysis/metabolome_age_adjusted_significant.csv", row.names = FALSE)

print("\nResults saved successfully!")
print(paste("Significant microbiome features (FDR < 0.05):", nrow(microbiome_significant)))
print(paste("Significant metabolome features (FDR < 0.05):", nrow(metabolome_significant)))

# ------------------------------------------------------------------------------
# 6. Generate visualization
# ------------------------------------------------------------------------------
# 6.1 Volcano plot for microbiome data
cat("\nGenerating volcano plots...\n")

# Prepare data for volcano plot
microbiome_volcano <- microbiome_age_adjusted_results %>%
  mutate(
    Neg_Log10_P_Adj = -log10(Adjusted_P_Value),
    Significant = Adjusted_P_Value < 0.05
  ) %>%
  mutate(
    Phylum_Genus = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), Phylum)
  )

# Create volcano plot
p_microbiome_volcano <- ggplot(microbiome_volcano, aes(x = Estimate, y = Neg_Log10_P_Adj, color = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Age-adjusted Differential Abundance Analysis (Microbiome)",
       x = "Estimate (OSAHS vs Control)",
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
top_microbiome <- microbiome_volcano %>%
  filter(Significant) %>%
  arrange(desc(Neg_Log10_P_Adj)) %>%
  head(10)

if (nrow(top_microbiome) > 0) {
  p_microbiome_volcano <- p_microbiome_volcano +
    geom_text_repel(
      data = top_microbiome,
      aes(label = Phylum_Genus),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.5,
      segment.color = "gray50"
    )
}

# Save volcano plot
ggsave("results/figures/differential_analysis/microbiome_age_adjusted_volcano.png", p_microbiome_volcano,
       width = 12, height = 8, dpi = 300)

# 6.2 Volcano plot for metabolome data
# Prepare data for volcano plot
metabolome_volcano <- metabolome_age_adjusted_results %>%
  mutate(
    Neg_Log10_P_Adj = -log10(Adjusted_P_Value),
    Significant = Adjusted_P_Value < 0.05,
    Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Metabolite)
  )

# Create volcano plot
p_metabolome_volcano <- ggplot(metabolome_volcano, aes(x = Estimate, y = Neg_Log10_P_Adj, color = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Age-adjusted Differential Abundance Analysis (Metabolome)",
       x = "Estimate (OSAHS vs Control)",
       y = "-Log10(Adjusted P-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Add labels for top significant metabolites
top_metabolite <- metabolome_volcano %>%
  filter(Significant) %>%
  arrange(desc(Neg_Log10_P_Adj)) %>%
  head(10)

if (nrow(top_metabolite) > 0) {
  p_metabolome_volcano <- p_metabolome_volcano +
    geom_text_repel(
      data = top_metabolite,
      aes(label = Metabolite_Label),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.5,
      segment.color = "gray50"
    )
}

# Save volcano plot
ggsave("results/figures/differential_analysis/metabolome_age_adjusted_volcano.png", p_metabolome_volcano,
       width = 12, height = 8, dpi = 300)

# 6.3 Boxplots for top significant features
cat("\nGenerating boxplots for top significant features...\n")

# Function to create boxplots with age as a covariate
create_age_adjusted_boxplots <- function(data, metadata, feature_ids, feature_type = "microbiome") {
  # Subset data to top features
  if (feature_type == "microbiome") {
    plot_data <- as.data.frame(t(data[feature_ids, ])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Abundance") %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
    
    # Add taxonomy information
    taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
    plot_data <- plot_data %>%
      left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature") %>%
      mutate(Feature_Label = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), Feature))
  } else {
    plot_data <- as.data.frame(t(data[feature_ids, ])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Abundance") %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
    
    # Add metabolite names if available
    if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
      plot_data <- plot_data %>%
        left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature") %>%
        mutate(Feature_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Feature))
    } else {
      plot_data <- plot_data %>%
        mutate(Feature_Label = Feature)
    }
  }
  
  # Create boxplots with age as a point size
  p <- ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot() +
    geom_point(aes(size = Age), alpha = 0.6, position = position_jitter(width = 0.2)) +
    facet_wrap(~Feature_Label, scales = "free_y") +
    labs(title = paste("Top Significant", str_to_title(feature_type), "Features (Age-adjusted)"),
         x = "Group",
         y = "Normalized Abundance",
         size = "Age") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  return(p)
}

# Create microbiome boxplots
if (nrow(microbiome_significant) > 0) {
  top_microbiome_features <- microbiome_significant$OTU[1:min(5, nrow(microbiome_significant))]
  p_microbiome_boxplot <- create_age_adjusted_boxplots(
    otu_table,
    metadata,
    top_microbiome_features,
    feature_type = "microbiome"
  )
  
  ggsave("results/figures/differential_analysis/microbiome_age_adjusted_boxplots.png", p_microbiome_boxplot,
         width = 15, height = 10, dpi = 300)
  
  print("Microbiome boxplots generated successfully!")
} else {
  print("No significant microbiome features found. Skipping boxplot generation.")
}

# Create metabolome boxplots
if (nrow(metabolome_significant) > 0) {
  top_metabolome_features <- metabolome_significant$Metabolite[1:min(5, nrow(metabolome_significant))]
  p_metabolome_boxplot <- create_age_adjusted_boxplots(
    metabolome_data,
    metadata,
    top_metabolome_features,
    feature_type = "metabolome"
  )
  
  ggsave("results/figures/differential_analysis/metabolome_age_adjusted_boxplots.png", p_metabolome_boxplot,
         width = 15, height = 10, dpi = 300)
  
  print("Metabolome boxplots generated successfully!")
} else {
  print("No significant metabolome features found. Skipping boxplot generation.")
}

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------
cat("\n=== Age_adjusted_analysis.R Summary ===\n")
cat("\nMicrobiome Analysis:\n")
print(paste("Total OTUs tested:", nrow(microbiome_age_adjusted_results)))
print(paste("Significant OTUs (FDR < 0.05):", sum(microbiome_age_adjusted_results$Adjusted_P_Value < 0.05)))

cat("\nMetabolome Analysis:\n")
print(paste("Total metabolites tested:", nrow(metabolome_age_adjusted_results)))
print(paste("Significant metabolites (FDR < 0.05):", sum(metabolome_age_adjusted_results$Adjusted_P_Value < 0.05)))

if (sum(microbiome_age_adjusted_results$Adjusted_P_Value < 0.05) > 0) {
  print("\nTop 5 significant microbiome features:")
  print(microbiome_significant %>%
          select(OTU, Phylum, Genus, Estimate, Adjusted_P_Value) %>%
          head(5))
}

if (sum(metabolome_age_adjusted_results$Adjusted_P_Value < 0.05) > 0) {
  print("\nTop 5 significant metabolome features:")
  print(metabolome_significant %>%
          select(Metabolite, Metabolite.Name, Estimate, Adjusted_P_Value) %>%
          head(5))
}

print("\nFiles generated:")
print("1. results/tables/differential_analysis/microbiome_age_adjusted_results.csv")
print("2. results/tables/differential_analysis/microbiome_age_adjusted_significant.csv")
print("3. results/tables/differential_analysis/metabolome_age_adjusted_results.csv")
print("4. results/tables/differential_analysis/metabolome_age_adjusted_significant.csv")
print("5. results/figures/differential_analysis/microbiome_age_adjusted_volcano.png")
print("6. results/figures/differential_analysis/metabolome_age_adjusted_volcano.png")
if (nrow(microbiome_significant) > 0) {
  print("7. results/figures/differential_analysis/microbiome_age_adjusted_boxplots.png")
}
if (nrow(metabolome_significant) > 0) {
  print("8. results/figures/differential_analysis/metabolome_age_adjusted_boxplots.png")
}
