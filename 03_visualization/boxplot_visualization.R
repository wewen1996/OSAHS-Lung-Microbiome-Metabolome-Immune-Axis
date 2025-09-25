# boxplot_visualization.R
# Generate boxplots for microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/figures/boxplots", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load normalized data
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load differential analysis results
microbiome_significant <- read.csv("results/tables/differential_analysis/microbiome_t_test_significant.csv", stringsAsFactors = FALSE)
metabolome_significant <- read.csv("results/tables/differential_analysis/metabolomics_oplsda_significant.csv", stringsAsfactor = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsfactor = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Function to create boxplots
# ------------------------------------------------------------------------------
create_boxplot <- function(data, metadata, feature_ids, feature_type = "microbiome", group_var = "Group", facet_by = NULL) {
  # Transpose data if features are rows
  if (feature_type == "Microbiome" || feature_type == "Metabolome") {
    plot_data <- as.data.frame(t(data[feature_ids, ])) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Abundance") %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
  } else {
    plot_data <- data %>%
      filter(Feature %in% feature_ids)
  }
  
  # Create feature labels
  if (feature_type == "Microbiome") {
    taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
    plot_data <- plot_data %>%
      left_join(taxonomy %>% rownames_to_column("Feature"), by = "Feature") %>%
      mutate(Feature_Label = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), 
                                    ifelse(!is.na(Family), paste(Phylum, Family, sep = "|"), Feature)))
  } else if (feature_type == "Metabolome") {
    if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
      plot_data <- plot_data %>%
        left_join(metabolite_annotations %>% rownames_to_column("Feature"), by = "Feature") %>%
        mutate(Feature_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Feature))
    } else {
      plot_data <- plot_data %>%
        mutate(Feature_Label = Feature)
    }
  } else {
    plot_data <- plot_data %>%
      mutate(Feature_Label = Feature)
  }
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = !!sym(group_var), y = Abundance, fill = !!sym(group_var))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    labs(title = paste("Distribution of", feature_type, "Features by", str_to_title(group_var)),
         x = str_to_title(group_var),
         y = "Normalized Abundance") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    scale_fill_viridis(discrete = TRUE, option = "plasma")
  
  # Add significance labels
  if (length(unique(plot_data[[group_var]])) == 2) {
    p <- p +
      stat_compare_means(method = "t.test", label = "p.signif", label.y = max(plot_data$Abundance) * 1.1)
  } else {
    p <- p +
      stat_compare_means(method = "anova", label = "p.signif", label.y = max(plot_data$Abundance) * 1.1)
  }
  
  # Add facet if specified
  if (!is.null(facet_by)) {
    p <- p +
      facet_wrap(as.formula(paste("~", facet_by)), scales = "free_y")
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. Generate microbiome boxplots
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome boxplots ===\n")

# Function to generate microbiome boxplots at different taxonomic levels
generate_microbiome_boxplots <- function(phyloseq_obj, tax_level, top_n = 10, group_var = "Group") {
  # Agglomerate taxa at the specified level
  phyloseq_agg <- tax_glom(phyloseq_obj, taxrank = tax_level)
  
  # Get abundance matrix
  otu_table <- otu_table(phyloseq_agg) %>% as.matrix()
  
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_agg) %>% as.data.frame()
  
  # Calculate mean abundance by group
  group_means <- aggregate(t(otu_table), by = list(Group = sample_data(phyloseq_agg)$Group), FUN = mean)
  rownames(group_means) <- group_means$Group
  group_means <- group_means[, -1]
  
  # Calculate fold change for two-group comparison
  if (nrow(group_means) == 2) {
    fold_change <- group_means[1, ] / group_means[2, ]
    top_features <- names(sort(abs(log2(fold_change)), decreasing = TRUE))[1:min(top_n, length(fold_change))]
  } else {
    # For more than two groups, use variance to select top features
    feature_variance <- apply(otu_table, 1, var)
    top_features <- names(sort(feature_variance, decreasing = TRUE))[1:min(top_n, length(feature_variance))]
  }
  
  # Create boxplot
  p <- create_boxplot(
    data = otu_table,
    metadata = sample_data(phyloseq_agg) %>% as.data.frame(),
    feature_ids = top_features,
    feature_type = "Microbiome",
    group_var = group_var,
    facet_by = "Feature_Label"
  )
  
  # Save plot
  ggsave(paste0("results/figures/boxplots/microbiome_", tolower(tax_level), "_boxplots.png"), p,
         width = 15, height = 8, dpi = 300)
  
  return(p)
}

# Generate boxplots at different taxonomic level
if (ntaxa(phyloseq_css) > 0) {
  # Phylum level
  generate_microbiome_boxplots(
    phyloseq_css,
    "Phylum",
    top_n = 10,
    group_var = "Group"
  )
  
  # Genus level
  generate_microbiome_boxplots(
    phyloseq_css,
    "Genus",
    top_n = 10,
    group_var = "Group"
  )
  
  print("Microbiome taxonomic boxplots generated successfully!")
} else {
  print("No microbiome data available for boxplot generation.")
}

# ------------------------------------------------------------------------------
# 4. Generate significant features boxplots
# ------------------------------------------------------------------------------
cat("\n=== Generating significant features boxplots ===\n")

# Generate significant microbiome features boxplots
if (nrow(microbiome_significant) > 0) {
  # Get microbiome data
  microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()
  
  # Select top N significant features
  top_n <- min(10, nrow(microbiome_significant))
  top_features <- microbiome_significant$Taxon[1:top_n]
  
  # Create boxplot
  p_microbiome_significant <- create_boxplot(
    data = microbiome_data,
    metadata = metadata,
    feature_ids = top_features,
    feature_type = "Microbiome",
    group_var = "Group",
    facet_by = "Feature_Label"
  )
  
  # Save plot
  ggsave("results/figures/boxplots/microbiome_significant_boxplots.png", p_microbiome_significant,
         width = 15, height = 8, dpi = 300)
  
  print("Significant microbiome features boxplots generated successfully!")
} else {
  print("No significant microbiome features found. Skipping boxplot generation.")
}

# Generate significant metabolome features boxplots
if (nrow(metabolome_significant) > 0) {
  # Select top N significant features
  top_n <- min(10, nrow(metabolome_significant))
  top_features <- metabolome_significant$Metabolite[1:top_n]
  
  # Create boxplot
  p_metabolome_significant <- create_boxplot(
    data = metabolome_data,
    metadata = metadata,
    feature_ids = top_features,
    feature_type = "Metabolome",
    group_var = "Group",
    facet_by = "Feature_Label"
  )
  
  # Save plot
  ggsave("results/figures/boxplots/metabolome_significant_boxplots.png", p_metabolome_significant,
         width = 15, height = 8, dpi = 300)
  
  print("Significant metabolome features boxplots generated successfully!")
} else {
  print("No significant metabolome features found. Skipping boxplot generation.")
}

# ------------------------------------------------------------------------------
# 5. Generate boxplots for alpha diversity indices
# ------------------------------------------------------------------------------
cat("\n=== Generating alpha diversity boxplots ===\n")

# Function to generate alpha diversity boxplots
generate_alpha_diversity_boxplots <- function(phyloseq_obj, measures = c("Shannon", "Simpson", "Chao1"), group_var = "Group") {
  # Calculate alpha diversity
  alpha_diversity <- estimate_richness(phyloseq_obj, measures = measures)
  
  # Add metadata
  alpha_diversity <- alpha_diversity %>%
    rownames_to_column("Sample") %>%
    left_join(sample_data(phyloseq_obj) %>% as.data.frame() %>% rownames_to_column("Sample"), by = "Sample")
  
  # Melt data for plotting
  alpha_diversity_melt <- alpha_diversity %>%
    pivot_longer(cols = all_of(measures), names_to = "Measure", values_to = "Value")
  
  # Create boxplot
  p <- ggplot(alpha_diversity_melt, aes(x = !!sym(group_var), y = Value, fill = !!sym(group_var))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    facet_wrap(~Measure, scales = "free_y") +
    labs(title = "Alpha Diversity Indices by Group",
         x = str_to_title(group_var),
         y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    scale_fill_viridis(discrete = TRUE, option = "plasma") +
    stat_compare_means(method = "t.test", label = "p.signif", label.y = "max")
  
  # Save plot
  ggsave("results/figures/boxplots/alpha_diversity_boxplots.png", p,
         width = 15, height = 8, dpi = 300)
  
  return(p)
}

# Generate alpha diversity boxplots
if (ntaxa(phyloseq_css) > 0) {
  generate_alpha_diversity_boxplots(
    phyloseq_css,
    measures = c("Shannon", "Simpson", "Chao1"),
    group_var = "Group"
  )
  
  print("Alpha diversity boxplots generated successfully!")
} else {
  print("No microbiome data available for alpha diversity analysis.")
}

# ------------------------------------------------------------------------------
# 6. Generate boxplots for clinical metadata
# ------------------------------------------------------------------------------
cat("\n=== Generating clinical metadata boxplots ===\n")

# Function to generate clinical metadata boxplots
generate_metadata_boxplots <- function(metadata, numeric_vars = NULL, group_var = "Group") {
  # Select numeric variables if not specified
  if (is.null(numeric_vars)) {
    numeric_vars <- metadata %>%
      select_if(is.numeric) %>%
      select(-any_of(c("SampleID", "PatientID"))) %>%  # Exclude ID columns
      colnames()
  }
  
  if (length(numeric_vars) == 0) {
    print("No numeric metadata variables found for boxplot generation.")
    return(NULL)
  }
  
  # Select top 6 numeric variables
  numeric_vars <- numeric_vars[1:min(6, length(numeric_vars))]
  
  # Melt data for plotting
  metadata_melt <- metadata %>%
    pivot_longer(cols = all_of(numeric_vars), names_to = "Variable", values_to = "Value")
  
  # Create boxplot
  p <- ggplot(metadata_melt, aes(x = !!sym(group_var), y = Value, fill = !!sym(group_var))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    facet_wrap(~Variable, scales = "free_y") +
    labs(title = "Clinical Metadata by Group",
         x = str_to_title(group_var),
         y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    scale_fill_viridis(discrete = TRUE, option = "plasma") +
    stat_compare_means(method = "t.test", label = "p.signif", label.y = "max")
  
  # Save plot
  ggsave("results/figures/boxplots/metadata_boxplots.png", p,
         width = 15, height = 8, dpi = 300)
  
  return(p)
}

# Generate metadata boxplots
if (ncol(metadata) > 1) {
  generate_metadata_boxplots(
    metadata,
    numeric_vars = NULL,
    group_var = "Group"
  )
  
  print("Clinical metadata boxplots generated successfully!")
} else {
  print("Not enough metadata variables for boxplot generation.")
}

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------
cat("\n=== boxplot_visualization.R Summary ===\n")
cat("\nGenerated boxplots:\n")

if (ntaxa(phyloseq_css) > 0) {
  cat("1. Microbiome taxonomic boxplots (Phylum and Genus levels)\n")
  cat("2. Alpha diversity boxplots\n")
}

if (nrow(microbiome_significant) > 0) {
  cat("3. Significant microbiome features boxplots\n")
}

if (nrow(metabolome_significant) > 0) {
  cat("4. Significant metabolome features boxplots\n")
}

if (ncol(metadata) > 1) {
  cat("5. Clinical metadata boxplots\n")
}

cat("\nAll boxplots saved to: results/figures/boxplots/\n")
