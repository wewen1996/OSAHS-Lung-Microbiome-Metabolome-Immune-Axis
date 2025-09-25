# stacked_barplots.R
# Generate stacked barplots for microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(grid)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/figures/stacked_barplots", showWarnings = FALSE, recursive = TRUE)

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
# 2. Function to create stacked barplot
# ------------------------------------------------------------------------------
create_stacked_barplot <- function(data, 
                                   metadata, 
                                   group_var = "Group", 
                                   tax_level = "Phylum", 
                                   top_n = 10, 
                                   title = "Stacked Barplot", 
                                   xlab = "Sample", 
                                   ylab = "Relative Abundance",
                                   color_palette = "viridis") {
  
  # Agglomerate taxa at the specified level
  if (inherits(data, "phyloseq")) {
    data_agg <- tax_glom(data, taxrank = tax_level)
    
    # Convert to relative abundance
    data_agg_rel <- transform_sample_counts(data_agg, function(x) x / sum(x) * 100)
    
    # Convert to data frame
    plot_data <- psmelt(data_agg_rel)
    
    # Rename columns
    colnames(plot_data)[colnames(plot_data) == tax_level] <- "Taxon"
    
  } else {
    # For metabolome data
    # Transpose data
    data_t <- t(data) %>% as.data.frame()
    
    # Add metadata
    plot_data <- data_t %>%
      rownames_to_column("Sample") %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample") %>%
      pivot_longer(cols = -c(Sample, !!sym(group_var)), names_to = "Metabolite", values_to = "Abundance")
    
    # If metabolite annotations are available, use class information
    if (!is.null(metabolite_annotations) && "Class" %in% colnames(metabolite_annotations)) {
      plot_data <- plot_data %>%
        left_join(metabolite_annotations %>% select(Class) %>% rownames_to_column("Metabolite"), by = "Metabolite") %>%
        rename(Taxon = Class)
    } else {
      plot_data <- plot_data %>%
        rename(Taxon = Metabolite)
    }
  }
  
  # Calculate total abundance for each taxon
  taxon_totals <- plot_data %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(Total))
  
  # Select top N taxa
  top_taxa <- taxon_totals$Taxon[1:min(top_n, nrow(taxon_totals))]
  
  # Create "Other" category for remaining taxa
  plot_data <- plot_data %>%
    mutate(
      Taxon_Category = ifelse(Taxon %in% top_taxa, Taxon, "Other")
    )
  
  # Calculate abundances for categories
  plot_data_summary <- plot_data %>%
    group_by(Sample, !!sym(group_var), Taxon_Category) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
    ungroup()
  
  # Order samples by group
  plot_data_summary <- plot_data_summary %>%
    mutate(Sample = factor(Sample, levels = unique(Sample[order(!!sym(group_var))])))
  
  # Create stacked barplot
  p <- ggplot(plot_data_summary, aes(x = Sample, y = Abundance, fill = Taxon_Category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  # Apply color palette
  if (color_palette == "viridis") {
    p <- p +
      scale_fill_viridis(discrete = TRUE, option = "plasma")
  } else {
    # Use custom color palette
    n_colors <- length(unique(plot_data_summary$Taxon_Category))
    p <- p +
      scale_fill_brewer(palette = color_palette, direction = -1)
  }
  
  # Add group labels
  if (!is.null(group_var)) {
    # Calculate group positions
    group_positions <- plot_data_summary %>%
      group_by(!!sym(group_var)) %>%
      summarise(Position = mean(as.numeric(Sample))) %>%
      ungroup()
    
    # Add group labels
    p <- p +
      geom_vline(data = group_positions, aes(xintercept = Position + 0.5), color = "black", linetype = "dashed") +
      annotate("text", x = group_positions$Position, y = max(plot_data_summary$Abundance) * 1.05,
               label = group_positions[[group_var]], size = 4)
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. Generate microbiome stacked barplots
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome stacked barplots ===\n")

# 3.1 Phylum level stacked barplot
if (ntaxa(phyloseq_css) > 0) {
  p_phylum <- create_stacked_barplot(
    data = phyloseq_css,
    metadata = metadata,
    group_var = "Group",
    tax_level = "Phylum",
    top_n = 10,
    title = "Microbiome Composition at Phylum Level",
    xlab = "Sample",
    ylab = "Relative Abundance (%)",
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/stacked_barplots/microbiome_phylum_stacked_barplot.png", p_phylum,
         width = 15, height = 10, dpi = 300)
  
  print("Microbiome phylum level stacked barplot generated successfully!")
} else {
  print("No microbiome data available for stacked barplot generation.")
}

# 3.2 Genus level stacked barplot
if (ntaxa(phyloseq_css) > 0) {
  p_genus <- create_stacked_barplot(
    data = phyloseq_css,
    metadata = metadata,
    group_var = "Group",
    tax_level = "Genus",
    top_n = 15,
    title = "Microbiome Composition at Genus Level",
    xlab = "Sample",
    ylab = "Relative Abundance (%)",
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/stacked_barplots/microbiome_genus_stacked_barplot.png", p_genus,
         width = 15, height = 10, dpi = 300)
  
  print("Microbiome genus level stacked barplot generated successfully!")
} else {
  print("No microbiome data available for stacked barplot generation.")
}

# 3.3 Group mean stacked barplot
cat("\n=== Generating group mean stacked barplots ===\n")

# Function to create group mean stacked barplot
create_group_mean_stacked_barplot <- function(data, 
                                             metadata, 
                                             group_var = "Group", 
                                             tax_level = "Phylum", 
                                             top_n = 10, 
                                             title = "Group Mean Stacked Barplot", 
                                             xlab = "Group", 
                                             ylab = "Mean Relative Abundance (%)",
                                             color_palette = "viridis") {
  
  # Agglomerate taxa at the specified level
  data_agg <- tax_glom(data, taxrank = tax_level)
  
  # Convert to relative abundance
  data_agg_rel <- transform_sample_counts(data_agg, function(x) x / sum(x) * 100)
  
  # Convert to data frame
  plot_data <- psmelt(data_agg_rel)
  
  # Rename columns
  colnames(plot_data)[colnames(plot_data) == tax_level] <- "Taxon"
  
  # Calculate total abundance for each taxon
  taxon_totals <- plot_data %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(Total))
  
  # Select top N taxa
  top_taxa <- taxon_totals$Taxon[1:min(top_n, nrow(taxon_totals))]
  
  # Create "Other" category for remaining taxa
  plot_data <- plot_data %>%
    mutate(
      Taxon_Category = ifelse(Taxon %in% top_taxa, Taxon, "Other")
    )
  
  # Calculate mean abundances by group
  plot_data_summary <- plot_data %>%
    group_by(!!sym(group_var), Taxon_Category) %>%
    summarise(
      Mean_Abundance = mean(Abundance, na.rm = TRUE),
      SD_Abundance = sd(Abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Create stacked barplot
  p <- ggplot(plot_data_summary, aes(x = !!sym(group_var), y = Mean_Abundance, fill = Taxon_Category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  # Apply color palette
  if (color_palette == "viridis") {
    p <- p +
      scale_fill_viridis(discrete = TRUE, option = "plasma")
  } else {
    # Use custom color palette
    n_colors <- length(unique(plot_data_summary$Taxon_Category))
    p <- p +
      scale_fill_brewer(palette = color_palette, direction = -1)
  }
  
  return(p)
}

# Generate group mean stacked barplot at phylum level
if (ntaxa(phyloseq_css) > 0) {
  p_group_mean_phylum <- create_group_mean_stacked_barplot(
    data = phyloseq_css,
    metadata = metadata,
    group_var = "Group",
    tax_level = "Phylum",
    top_n = 10,
    title = "Mean Microbiome Composition at Phylum Level",
    xlab = "Group",
    ylab = "Mean Relative Abundance (%)",
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/stacked_barplots/microbiome_group_mean_phylum_stacked_barplot.png", p_group_mean_phylum,
         width = 12, height = 8, dpi = 300)
  
  print("Microbiome group mean phylum level stacked barplot generated successfully!")
} else {
  print("No microbiome data available for group mean stacked barplot generation.")
}

# Generate group mean stacked barplot at genus level
if (ntaxa(phyloseq_css) > 0) {
  p_group_mean_genus <- create_group_mean_stacked_barplot(
    data = phyloseq_css,
    metadata = metadata,
    group_var = "Group",
    tax_level = "Genus",
    top_n = 15,
    title = "Mean Microbiome Composition at Genus Level",
    xlab = "Group",
    ylab = "Mean Relative Abundance (%)",
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/stacked_barplots/microbiome_group_mean_genus_stacked_barplot.png", p_group_mean_genus,
         width = 12, height = 8, dpi = 300)
  
  print("Microbiome group mean genus level stacked barplot generated successfully!")
} else {
  print("No microbiome data available for group mean stacked barplot generation.")
}

# ------------------------------------------------------------------------------
# 4. Generate metabolome stacked barplots
# ------------------------------------------------------------------------------
cat("\n=== Generating metabolome stacked barplots ===\n")

# 4.1 Metabolite class stacked barplot (if annotations are available)
if (nrow(metabolome_data) > 0 && !is.null(metabolite_annotations) && "Class" %in% colnames(metabolite_annotations)) {
  p_metabolite_class <- create_stacked_barplot(
    data = metabolome_data,
    metadata = metadata,
    group_var = "Group",
    top_n = 10,
    title = "Metabolite Composition by Class",
    xlab = "Sample",
    ylab = "Normalized Abundance",
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/stacked_barplots/metabolite_class_stacked_barplot.png", p_metabolite_class,
         width = 15, height = 10, dpi = 300)
  
  print("Metabolite class stacked barplot generated successfully!")
} else {
  print("No metabolite class information available. Skipping metabolite class stacked barplot generation.")
}

# 4.2 Top metabolites stacked barplot
if (nrow(metabolome_data) > 0) {
  # Transpose data
  metabolome_data_t <- t(metabolome_data) %>% as.data.frame()
  
  # Add metadata
  plot_data <- metabolome_data_t %>%
    rownames_to_column("Sample") %>%
    left_join(metadata %>% rownames_to_column("Sample"), by = "Sample") %>%
    pivot_longer(cols = -c(Sample, Group), names_to = "Metabolite", values_to = "Abundance")
  
  # Calculate total abundance for each metabolite
  metabolite_totals <- plot_data %>%
    group_by(Metabolite) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(Total))
  
  # Select top N metabolites
  top_n <- 15
  top_metabolites <- metabolite_totals$Metabolite[1:min(top_n, nrow(metabolite_totals))]
  
  # Create "Other" category for remaining metabolites
  plot_data <- plot_data %>%
    mutate(
      Metabolite_Category = ifelse(Metabolite %in% top_metabolites, Metabolite, "Other")
    )
  
  # Add metabolite names if available
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    plot_data <- plot_data %>%
      left_join(
        metabolite_annotations %>%
          select(Metabolite.Name) %>%
          rownames_to_column("Metabolite"),
        by = "Metabolite"
      ) %>%
      mutate(
        Metabolite_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Metabolite)
      )
    
    # Update category labels
    plot_data <- plot_data %>%
      mutate(
        Metabolite_Category = ifelse(Metabolite %in% top_metabolites, Metabolite_Label, "Other")
      )
  }
  
  # Order samples by group
  plot_data <- plot_data %>%
    mutate(Sample = factor(Sample, levels = unique(Sample[order(Group)])))
  
  # Create stacked barplot
  p_top_metabolites <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Metabolite_Category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = "Top Metabolite Composition",
      x = "Sample",
      y = "Normalized Abundance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    ) +
    scale_fill_viridis(discrete = TRUE, option = "plasma")
  
  # Add group labels
  group_positions <- plot_data %>%
    group_by(Group) %>%
    summarise(Position = mean(as.numeric(Sample))) %>%
    ungroup()
  
  p_top_metabolites <- p_top_metabolites +
    geom_vline(data = group_positions, aes(xintercept = Position + 0.5), color = "black", linetype = "dashed") +
    annotate("text", x = group_positions$Position, y = max(plot_data$Abundance) * 1.05,
             label = group_positions$Group, size = 4)
  
  # Save plot
  ggsave("results/figures/stacked_barplots/metabolite_top_stacked_barplot.png", p_top_metabolites,
         width = 15, height = 10, dpi = 300)
  
  print("Top metabolites stacked barplot generated successfully!")
} else {
  print("No metabolome data available for stacked barplot generation.")
}

# ------------------------------------------------------------------------------
# 5. Generate stacked barplots with clinical metadata
# ------------------------------------------------------------------------------
cat("\n=== Generating stacked barplots with clinical metadata ===\n")

# Function to create stacked barplot with clinical metadata
create_clinical_metadata_stacked_barplot <- function(data, 
                                                    metadata, 
                                                    clinical_var,
                                                    group_var = "Group", 
                                                    tax_level = "Phylum", 
                                                    top_n = 10, 
                                                    title = "Stacked Barplot with Clinical Metadata", 
                                                    xlab = "Clinical Variable", 
                                                    ylab = "Mean Relative Abundance (%)",
                                                    color_palette = "viridis") {
  
  # Check if clinical variable exists
  if (!clinical_var %in% colnames(metadata)) {
    print(paste("Clinical variable", clinical_var, "not found in metadata."))
    return(NULL)
  }
  
  # Agglomerate taxa at the specified level
  data_agg <- tax_glom(data, taxrank = tax_level)
  
  # Convert to relative abundance
  data_agg_rel <- transform_sample_counts(data_agg, function(x) x / sum(x) * 100)
  
  # Convert to data frame
  plot_data <- psmelt(data_agg_rel)
  
  # Rename columns
  colnames(plot_data)[colnames(plot_data) == tax_level] <- "Taxon"
  
  # Add clinical metadata
  plot_data <- plot_data %>%
    left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
  
  # Calculate total abundance for each taxon
  taxon_totals <- plot_data %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
    arrange(desc(Total))
  
  # Select top N taxa
  top_taxa <- taxon_totals$Taxon[1:min(top_n, nrow(taxon_totals))]
  
  # Create "Other" category for remaining taxa
  plot_data <- plot_data %>%
    mutate(
      Taxon_Category = ifelse(Taxon %in% top_taxa, Taxon, "Other")
    )
  
  # Calculate mean abundances by clinical variable and group
  plot_data_summary <- plot_data %>%
    group_by(!!sym(clinical_var), !!sym(group_var), Taxon_Category) %>%
    summarise(
      Mean_Abundance = mean(Abundance, na.rm = TRUE),
      SD_Abundance = sd(Abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Create stacked barplot
  p <- ggplot(plot_data_summary, aes(x = !!sym(clinical_var), y = Mean_Abundance, fill = Taxon_Category)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(as.formula(paste("~", group_var))) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  # Apply color palette
  if (color_palette == "viridis") {
    p <- p +
      scale_fill_viridis(discrete = TRUE, option = "plasma")
  } else {
    # Use custom color palette
    n_colors <- length(unique(plot_data_summary$Taxon_Category))
    p <- p +
      scale_fill_brewer(palette = color_palette, direction = -1)
  }
  
  return(p)
}

# Generate stacked barplot with clinical metadata (e.g., AgeGroup)
if (ntaxa(phyloseq_css) > 0 && "Age" %in% colnames(metadata)) {
  # Create age groups
  metadata$AgeGroup <- cut(metadata$Age, breaks = quantile(metadata$Age, seq(0, 1, 0.25)), include.lowest = TRUE)
  
  p_clinical_metadata <- create_clinical_metadata_stacked_barplot(
    data = phyloseq_css,
    metadata = metadata,
    clinical_var = "AgeGroup",
    group_var = "Group",
    tax_level = "Phylum",
    top_n = 10,
    title = "Microbiome Composition by Age Group and Disease Status",
    xlab = "Age Group",
    ylab = "Mean Relative Abundance (%)",
    color_palette = "viridis"
  )
  
  if (!is.null(p_clinical_metadata)) {
    # Save plot
    ggsave("results/figures/stacked_barplots/microbiome_age_group_stacked_barplot.png", p_clinical_metadata,
           width = 15, height = 10, dpi = 300)
    
    print("Microbiome age group stacked barplot generated successfully!")
  }
} else {
  print("No age metadata available. Skipping clinical metadata stacked barplot generation.")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== stacked_barplots.R Summary ===\n")
cat("\nGenerated stacked barplots:\n")

if (ntaxa(phyloseq_css) > 0) {
  cat("1. Microbiome phylum level stacked barplot\n")
  cat("2. microbiome genus level stacked barplot\n")
  cat("3. microbiome group mean phylum level stacked barplot\n")
  cat("4. microbiome group mean genus level stacked barplot\n")
}

if (nrow(metabolome_data) > 0 && !is.null(metabolite_annotations) && "Class" %in% colnames(metabolite_annotations)) {
  cat("5. Metabolite class stacked barplot\n")
}

if (nrow(metabolome_data) > 0) {
  cat("6. Top metabolites stacked barplot\n")
}

if (ntaxa(phyloseq_css) > 0 && "Age" %in% colnames(metadata)) {
  cat("7. microbiome age group stacked barplot\n")
}

cat("\nAll stacked barplots saved to: results/figures/stacked_barplots/\n")
