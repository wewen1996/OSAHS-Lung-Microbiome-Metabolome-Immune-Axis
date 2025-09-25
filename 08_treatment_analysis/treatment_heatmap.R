# treatment_heatmap.R
# Analyze and visualize treatment effects on microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directories if they don't exist
dir.create("results/treatment_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/treatment_analysis", showWarnings = FALSE, recursive = TRUE)

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
# 2. Check for treatment data
# ------------------------------------------------------------------------------
# Check if treatment metadata is available
required_cols <- c("Treatment", "TimePoint")

if (!all(required_cols %in% colnames(metadata))) {
  stop(paste("Treatment metadata not found. Required columns:", paste(required_cols, collapse = ", ")))
}

# Check if there are at least two time points
if (length(unique(metadata$TimePoint)) < 2) {
  stop("At least two time points are required for treatment analysis.")
}

# Check if there are paired samples (same individual at different time points)
if (!"PatientID" %in% colnames(metadata)) {
  warning("PatientID not found in metadata. Cannot perform paired analysis.")
  paired_analysis <- FALSE
} else {
  paired_analysis <- TRUE
  
  # Check if there are patients with samples at multiple time points
  patient_timepoint_count <- metadata %>%
    group_by(PatientID) %>%
    summarise(TimePoint_Count = n_distinct(TimePoint)) %>%
    filter(TimePoint_Count >= 2)
  
  if (nrow(patient_timepoint_count) == 0) {
    warning("No patient with samples at multiple time points. Cannot perform paired analysis.")
    paired_analysis <- FALSE
  } else {
    print(paste("Patient with samples at multiple time points:", nrow(patient_timepoint_count)))
  }
}

print("Treatment data checked successfully!")

# ------------------------------------------------------------------------------
# 3. Prepare data for treatment analysis
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(sample_names(phyloseq_css), colnames(metabolome_data))
common_samples <- intersect(common_samples, rownames(metadata))

if (length(common_samples) < 3) {
  stop("Not enough common samples for treatment analysis.")
}

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples:", length(common_samples)))

# Create combined metadata for samples
sample_metadata <- data.frame(
  Sample = common_samples,
  Treatment = metadata$Treatment[common_samples],
  TimePoint = metadata$TimePoint[common_samples],
  Group = metadata$Group[common_samples]
)

if (paired_analysis) {
  sample_metadata$PatientID <- metadata$PatientID[common_samples]
}

print("Data prepared for treatment analysis!")

# ------------------------------------------------------------------------------
# 4. Function to create treatment heatmap
# ------------------------------------------------------------------------------
create_treatment_heatmap <- function(data, metadata, features, title = "Treatment Heatmap") {
  
  # Subset data to selected features
  data_subset <- data[features, ]
  
  # Create annotation data frame
  annotation_col <- metadata %>%
    select(Treatment, TimePoint, Group) %>%
    distinct()
  
  rownames(annotation_col) <- annotation_col$Sample
  annotation_col$Sample <- NULL
  
  # Create heatmap
  p <- pheatmap(
    data_subset,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    show_colnames = FALSE,
    treeheight_row = 20,
    treeheight_col = 20,
    fontsize_row = 8,
    main = title,
    color = viridis(100)
  )
  
  return(p)
}

# ------------------------------------------------------------------------------
# 5. Analyze treatment effects on microbiome
# ------------------------------------------------------------------------------
cat("\n=== Analyzing treatment effects on microbiome ===\n")

# 5.1 Prepare microbiome data
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Filter low abundance/ variance features
min_abundance <- 0.001  # 0.1% relative abundance
min_prevalence <- 0.2   # Present in at least 20% of samples

# Calculate relative abundance
microbiome_rel <- t(t(microbiome_data) / colSums(microbiome_data))

# Filter taxa
taxa_to_keep <- rowSums(microbiome_rel > min_abundance) >= (min_prevalence * ncol(microbiome_rel))
microbiome_filtered <- microbiome_data[taxa_to_keep, ]

print(paste("Microbiome taxa after filtering:", sum(taxa_to_keep)))

# 5.2 Identify treatment-responsive taxa
if (paired_analysis) {
  # For paired analysis, we can use paired t-test or mixed effects models
  
  # Get patient with samples at multiple time points
  patients_with_multiple_timepoints <- metadata %>%
    group_by(PatientID) %>%
    summarise(TimePoint_Count = n_distinct(TimePoint)) %>%
    filter(TimePoint_Count >= 2) %>%
    pull(PatientID)
  
  # Subset metadata to these patients
  paired_metadata <- metadata %>%
    filter(PatientID %in% patients_with_multiple_timepoints)
  
  # Get paired samples
  paired_samples <- paired_metadata$Sample
  
  # Subset microbiome data
  microbiome_paired <- microbiome_filtered[, paired_samples]
  
  print(paste("Paired samples for microbiome analysis:", length(paired_samples)))
  
  if (length(paired_samples) >= 6) {  # Need at least 3 pairs for meaningful analysis
    # Function to perform paired analysis for each taxon
    perform_paired_analysis <- function(data, metadata, patient_id_col, timepoint_col) {
      
      results <- data.frame(
        Taxon = character(),
        TimePoint1 = character(),
        TimePoint2 = character(),
        Log2_Fold_Change = numeric(),
        P_Value = numeric(),
        Adjusted_P_Value = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Get unique time points
      time_points <- sort(unique(metadata[[timepoint_col]]))
      
      if (length(time_points) >= 2) {
        # Compare first and last time point
        timepoint1 <- time_points[1]
        timepoint2 <- time_points[length(time_points)]
        
        # Get samples for each time point
        samples_timepoint1 <- metadata %>%
          filter(!!sym(timepoint_col) == timepoint1) %>%
          pull(Sample)
        
        samples_timepoint2 <- metadata %>%
          filter(!!sym(timepoint_col) == timepoint2) %>%
          pull(Sample)
        
        # Get patients with samples at both time points
        patients_timepoint1 <- metadata %>%
          filter(Sample %in% samples_timepoint1) %>%
          pull(!!sym(patient_id_col))
        
        patients_timepoint2 <- metadata %>%
          filter(Sample %in% samples_timepoint2) %>%
          pull(!!sym(patient_id_col))
        
        common_patients <- intersect(patients_timepoint1, patients_timepoint2)
        
        if (length(common_patients) >= 3) {
          # Get samples for common patients at both time points
          paired_samples_timepoint1 <- metadata %>%
            filter(!!sym(patient_id_col) %in% common_patients & !!sym(timepoint_col) == timepoint1) %>%
            pull(Sample)
          
          paired_samples_timepoint2 <- metadata %>%
            filter(!!sym(patient_id_col) %in% common_patients & !!sym(timepoint_col) == timepoint2) %>%
            pull(Sample)
          
          # Ensure samples are in the same order
          paired_samples_timepoint1 <- paired_samples_timepoint1[order(metadata$PatientID[match(paired_samples_timepoint1, metadata$Sample)])]
          paired_samples_timepoint2 <- paired_samples_timepoint2[order(metadata$PatientID[match(paired_samples_timepoint2, metadata$Sample)])]
          
          # Perform paired t-test for each taxon
          for (taxon in rownames(data)) {
            # Get data for this taxon
            taxon_data_timepoint1 <- data[taxon, paired_samples_timepoint1]
            taxon_data_timepoint2 <- data[taxon, paired_samples_timepoint2]
            
            # Perform paired t-test
            t_test_result <- t.test(taxon_data_timepoint1, taxon_data_timepoint2, paired = TRUE)
            
            # Calculate log2 fold change (timepoint2 - timepoint1)
            log2_fold_change <- log2(mean(taxon_data_timepoint2) / mean(taxon_data_timepoint1))
            
            # Add to results
            results <- rbind(results, data.frame(
              Taxon = taxon,
              TimePoint1 = timepoint1,
              TimePoint2 = timepoint2,
              Log2_Fold_Change = log2_fold_change,
              P_Value = t_test_result$p.value,
              Adjusted_P_Value = NA,
              stringsAsFactors = FALSE
            ))
          }
          
          # Apply multiple testing correction
          results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
          
          # Order by adjusted p-value
          results <- results %>% arrange(Adjusted_P_Value)
        }
      }
      
      return(results)
    }
    
    # Perform paired analysis
    microbiome_treatment_results <- perform_paired_analysis(
      data = microbiome_paired,
      metadata = paired_metadata,
      patient_id_col = "PatientID",
      timepoint_col = "TimePoint"
    )
    
    if (nrow(microbiome_treatment_results) > 0) {
      # Save results
      write.csv(microbiome_treatment_results, "results/treatment_analysis/microbiome_treatment_results.csv", row.names = FALSE)
      
      # Identify significant taxa
      significant_taxa <- microbiome_treatment_results %>%
        filter(Adjusted_P_Value < 0.05) %>%
        pull(Taxon)
      
      print(paste("Significant treatment-responsive microbiome taxa:", length(significant_taxa)))
      
      if (length(significant_taxa) > 0) {
        # Get taxonomy
        taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
        
        # Add taxonomy information
        microbiome_treatment_results <- microbiome_treatment_results %>%
          left_join(taxonomy %>% rownames_to_column("Taxon"), by = "Taxon")
        
        # Save significant results
        microbiome_treatment_significant <- microbiome_treatment_results %>%
          filter(Adjusted_P_Value < 0.05) %>%
          arrange(desc(abs(Log2_Fold_Change)))
        
        write.csv(microbiome_treatment_significant, "results/treatment_analysis/microbiome_treatment_significant.csv", row.names = FALSE)
        
        print("Top 10 treatment-responsive microbiome taxa:")
        print(head(microbiome_treatment_significant[, c("Taxon", "Phylum", "Genus", "Log2_Fold_Change", "Adjusted_P_Value")], 10))
      }
    }
  } else {
    print("Not enough paired samples for microbiome treatment analysis.")
  }
} else {
  # For unpaired analysis, we can use ANOVA or Kruskal-Wallis test
  
  print("Performing unpaired treatment analysis for microbiome...")
  
  # Function to perform ANOVA for each taxon
  perform_anova_analysis <- function(data, metadata, timepoint_col) {
    
    results <- data.frame(
      Taxon = character(),
      F_Statistic = numeric(),
      P_Value = numeric(),
      Adjusted_P_Value = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Perform ANOVA for each taxon
    for (taxon in rownames(data)) {
      # Create data frame for this taxon
      taxon_data <- data.frame(
        Value = as.numeric(data[taxon, ]),
        TimePoint = metadata[[timepoint_col]]
      )
      
      # Perform ANOVA
      anova_result <- aov(Value ~ TimePoint, data = taxon_data)
      anova_summary <- summary(anova_result)
      
      # Add to results
      results <- rbind(results, data.frame(
        Taxon = taxon,
        F_Statistic = anova_summary[[1]]$F[1],
        P_Value = anova_summary[[1]]$Pr[1],
        Adjusted_P_Value = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Apply multiple testing correction
    results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
    
    # Order by adjusted p-value
    results <- results %>% arrange(Adjusted_P_Value)
    
    return(results)
  }
  
  # Perform ANOVA analysis
  microbiome_treatment_results <- perform_anova_analysis(
    data = microbiome_filtered,
    metadata = metadata,
    timepoint_col = "TimePoint"
  )
  
  if (nrow(microbiome_treatment_results) > 0) {
    # Save results
    write.csv(microbiome_treatment_results, "results/treatment_analysis/microbiome_treatment_results.csv", row.names = FALSE)
    
    # Identify significant taxa
    significant_taxa <- microbiome_treatment_results %>%
      filter(Adjusted_P_Value < 0.05) %>%
      pull(Taxon)
    
    print(paste("Significant treatment-responsive microbiome taxa:", length(significant_taxa)))
    
    if (length(significant_taxa) > 0) {
      # Get taxonomy
      taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
      
      # Add taxonomy information
      microbiome_treatment_results <- microbiome_treatment_results %>%
        left_join(taxonomy %>% rownames_to_column("Taxon"), by = "Taxon")
      
      # Save significant results
      microbiome_treatment_significant <- microbiome_treatment_results %>%
        filter(Adjusted_P_Value < 0.05) %>%
        arrange(desc(F_Statistic))
      
      write.csv(microbiome_treatment_significant, "results/treatment_analysis/microbiome_treatment_significant.csv", row.names = FALSE)
      
      print("Top 10 treatment-responsive microbiome taxa:")
      print(head(microbiome_treatment_significant[, c("Taxon", "Phylum", "Genus", "F_Statistic", "Adjusted_P_Value")], 10))
    }
  }
}

# 5.3 Create treatment heatmap for microbiome
if (exists("significant_taxa") && length(significant_taxa) > 0) {
  # Select top N taxa
  top_n <- min(50, length(significant_taxa))
  top_taxa <- significant_taxa[1:top_n]
  
  # Create heatmap
  p_microbiome_heatmap <- create_treatment_heatmap(
    data = microbiome_filtered[top_taxa, ],
    metadata = sample_metadata,
    features = top_taxa,
    title = "Treatment-Responsive Microbiome Taxa"
  )
  
  # Save heatmap
  png("results/figures/treatment_analysis/microbiome_treatment_heatmap.png", width = 14, height = 12, dpi = 300)
  print(p_microbiome_heatmap)
  dev.off()
  
  print("Microbiome treatment heatmap generated successfully!")
} else if (nrow(microbiome_filtered) > 0) {
  # If no significant taxa, create heatmap of top variable taxa
  top_variable_taxa <- names(sort(apply(microbiome_filtered, 1, var), decreasing = TRUE))[1:min(50, nrow(microbiome_filtered))]
  
  p_microbiome_heatmap <- create_treatment_heatmap(
    data = microbiome_filtered[top_variable_taxa, ],
    metadata = sample_metadata,
    features = top_variable_taxa,
    title = "Top Variable Microbiome Taxa (Treatment Analysis)"
  )
  
  # Save heatmap
  png("results/figures/treatment_analysis/microbiome_top_variable_heatmap.png", width = 14, height = 12, dpi = 300)
  print(p_microbiome_heatmap)
  dev.off()
  
  print("Microbiome top variable heatmap generated successfully!")
}

# ------------------------------------------------------------------------------
# 6. Analyze treatment effects on metabolome
# ------------------------------------------------------------------------------
cat("\n=== Analyzing treatment effects on metabolome ===\n")

# 6.1 Prepare metabolome data
# Filter low variance metabolites
min_variance <- quantile(apply(metabolome_data, 1, var), 0.25)  # Keep top 75% most variable metabolites
metabolite_variance <- apply(metabolome_data, 1, var)
metabolites_to_keep <- metabolite_variance >= min_variance
metabolome_filtered <- metabolome_data[metabolites_to_keep, ]

print(paste("Metabolites after filtering:", sum(metabolites_to_keep)))

# 6.2 Identify treatment-responsive metabolites
if (paired_analysis) {
  # For paired analysis, we can use paired t-test
  
  # Get patient with samples at multiple time points
  patients_with_multiple_timepoints <- metadata %>%
    group_by(PatientID) %>%
    summarise(TimePoint_Count = n_distinct(TimePoint)) %>%
    filter(TimePoint_Count >= 2) %>%
    pull(PatientID)
  
  # Subset metadata to these patients
  paired_metadata <- metadata %>%
    filter(PatientID %in% patients_with_multiple_timepoints)
  
  # Get paired samples
  paired_samples <- paired_metadata$Sample
  
  # Subset metabolome data
  metabolome_paired <- metabolome_filtered[, paired_samples]
  
  print(paste("Paired samples for metabolome analysis:", length(paired_samples)))
  
  if (length(paired_samples) >= 6) {  # Need at least 3 pairs for meaningful analysis
    # Perform paired analysis
    metabolome_treatment_results <- perform_paired_analysis(
      data = metabolome_paired,
      metadata = paired_metadata,
      patient_id_col = "PatientID",
      timepoint_col = "TimePoint"
    )
    
    if (nrow(metabolome_treatment_results) > 0) {
      # Save results
      write.csv(metabolome_treatment_results, "results/treatment_analysis/metabolome_treatment_results.csv", row.names = FALSE)
      
      # Identify significant metabolites
      significant_metabolites <- metabolome_treatment_results %>%
        filter(Adjusted_P_Value < 0.05) %>%
        pull(Taxon)
      
      print(paste("Significant treatment-responsive metabolites:", length(significant_metabolites)))
      
      if (length(significant_metabolites) > 0) {
        # Add metabolite annotations
        metabolome_treatment_results <- metabolome_treatment_results %>%
          rename(Metabolite = Taxon) %>%
          left_join(metabolite_annotations %>% rownames_to_column("Metabolite"), by = "Metabolite")
        
        # Save significant results
        metabolome_treatment_significant <- metabolome_treatment_results %>%
          filter(Adjusted_P_Value < 0.05) %>%
          arrange(desc(abs(Log2_Fold_Change)))
        
        write.csv(metabolome_treatment_significant, "results/treatment_analysis/metabolome_treatment_significant.csv", row.names = FALSE)
        
        print("Top 10 treatment-responsive metabolites:")
        print(head(metabolome_treatment_significant[, c("Metabolite", "Metabolite.Name", "Class", "Log2_Fold_Change", "Adjusted_P_Value")], 10))
      }
    }
  } else {
    print("Not enough paired samples for metabolome treatment analysis.")
  }
} else {
  # For unpaired analysis, we can use ANOVA
  
  print("Performing unpaired treatment analysis for metabolome...")
  
  # Perform ANOVA analysis
  metabolome_treatment_results <- perform_anova_analysis(
    data = metabolome_filtered,
    metadata = metadata,
    timepoint_col = "TimePoint"
  )
  
  if (nrow(metabolome_treatment_results) > 0) {
    # Save results
    write.csv(metabolome_treatment_results, "results/treatment_analysis/metabolome_treatment_results.csv", row.names = FALSE)
    
    # Identify significant metabolites
    significant_metabolites <- metabolome_treatment_results %>%
      filter(Adjusted_P_Value < 0.05) %>%
      pull(Taxon)
    
    print(paste("Significant treatment-responsive metabolites:", length(significant_metabolites)))
    
    if (length(significant_metabolites) > 0) {
      # Add metabolite annotations
      metabolome_treatment_results <- metabolome_treatment_results %>%
        rename(Metabolite = Taxon) %>%
        left_join(metabolite_annotations %>% rownames_to_column("Metabolite"), by = "Metabolite")
      
      # Save significant results
      metabolome_treatment_significant <- metabolome_treatment_results %>%
        filter(Adjusted_P_Value < 0.05) %>%
        arrange(desc(F_Statistic))
      
      write.csv(metabolome_treatment_significant, "results/treatment_analysis/metabolome_treatment_significant.csv", row.names = FALSE)
      
      print("Top 10 treatment-responsive metabolites:")
      print(head(metabolome_treatment_significant[, c("Metabolite", "Metabolite.Name", "Class", "F_Statistic", "Adjusted_P_Value")], 10))
    }
  }
}

# 6.3 Create treatment heatmap for metabolome
if (exists("significant_metabolites") && length(significant_metabolites) > 0) {
  # Select top N metabolites
  top_n <- min(50, length(significant_metabolites))
  top_metabolites <- significant_metabolites[1:top_n]
  
  # Create heatmap
  p_metabolome_heatmap <- create_treatment_heatmap(
    data = metabolome_filtered[top_metabolites, ],
    metadata = sample_metadata,
    features = top_metabolites,
    title = "Treatment-Responsive Metabolites"
  )
  
  # Save heatmap
  png("results/figures/treatment_analysis/metabolome_treatment_heatmap.png", width = 14, height = 12, dpi = 300)
  print(p_metabolome_heatmap)
  dev.off()
  
  print("Metabolome treatment heatmap generated successfully!")
} else if (nrow(metabolome_filtered) > 0) {
  # If no significant metabolites, create heatmap of top variable metabolites
  top_variable_metabolites <- names(sort(apply(metabolome_filtered, 1, var), decreasing = TRUE))[1:min(50, nrow(metabolome_filtered))]
  
  p_metabolome_heatmap <- create_treatment_heatmap(
    data = metabolome_filtered[top_variable_metabolites, ],
    metadata = sample_metadata,
    features = top_variable_metabolites,
    title = "Top Variable Metabolites (Treatment Analysis)"
  )
  
  # Save heatmap
  png("results/figures/treatment_analysis/metabolome_top_variable_heatmap.png", width = 14, height = 12, dpi = 300)
  print(p_metabolome_heatmap)
  dev.off()
  
  print("Metabolome top variable heatmap generated successfully!")
}

# ------------------------------------------------------------------------------
# 7. Visualize treatment effects
# ------------------------------------------------------------------------------
cat("\n=== Visualizing treatment effects ===\n")

# 7.1 Treatment response trajectory plot for top taxa/metabolites
if (exists("significant_taxa") && length(significant_taxa) > 0) {
  # Select top 5 taxa
  top_taxa <- significant_taxa[1:5]
  
  # Get taxonomy
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  # Create trajectory plot
  plot_list <- list()
  
  for (i in 1:length(top_taxa)) {
    taxon <- top_taxa[i]
    
    # Get taxon data
    taxon_data <- data.frame(
      Sample = colnames(microbiome_filtered),
      Abundance = as.numeric(microbiome_filtered[taxon, ])
    ) %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
    
    # Create taxonomic label
    taxonomic_label <- if (!is.na(taxonomy[taxon, "Genus"])) {
      paste(taxonomy[taxon, "Phylum"], taxonomy[taxon, "Genus"], sep = "|")
    } else {
      taxon
    }
    
    # Create plot
    p <- ggplot(taxon_data, aes(x = TimePoint, y = Abundance, color = Group, group = interaction(PatientID, Group))) +
      geom_line(alpha = 0.5) +
      geom_point(size = 2) +
      geom_smooth(aes(group = Group), method = "lm", se = TRUE, alpha = 0.8) +
      labs(
        title = taxonomic_label,
        x = "Time Point",
        y = "Abundance"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      )
    
    plot_list[[i]] <- p
  }
  
  # Combine plots
  if (length(plot_list) > 0) {
    p_combined <- grid.arrange(grobs = plot_list, ncol = 1)
    
    # Save plot
    ggsave("results/figures/treatment_analysis/microbiome_treatment_trajectory.png", p_combined,
           width = 12, height = 15, dpi = 300)
    
    print("Microbiome treatment trajectory plot generated successfully!")
  }
}

if (exists("significant_metabolites") && length(significant_metabolites) > 0) {
  # Select top 5 metabolites
  top_metabolites <- significant_metabolites[1:5]
  
  # Create trajectory plot
  plot_list <- list()
  
  for (i in 1:length(top_metabolites)) {
    metabolite <- top_metabolites[i]
    
    # Get metabolite data
    metabolite_data <- data.frame(
      Sample = colnames(metabolome_filtered),
      Abundance = as.numeric(metabolome_filtered[metabolite, ])
    ) %>%
      left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
    
    # Create metabolite label
    metabolite_label <- if (!is.null(metabolite_annotations) && 
                           metabolite %in% rownames(metabolite_annotations) && 
                           !is.na(metabolite_annotations[metabolite, "Metabolite.Name"])) {
      metabolite_annotations[metabolite, "Metabolite.Name"]
    } else {
      metabolite
    }
    
    # Create plot
    p <- ggplot(metabolite_data, aes(x = TimePoint, y = Abundance, color = Group, group = interaction(PatientID, Group))) +
      geom_line(alpha = 0.5) +
      geom_point(size = 2) +
      geom_smooth(aes(group = Group), method = "lm", se = TRUE, alpha = 0.8) +
      labs(
        title = metabolite_label,
        x = "Time Point",
        y = "Normalized Abundance"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom"
      )
    
    plot_list[[i]] <- p
  }
  
  # Combine plots
  if (length(plot_list) > 0) {
    p_combined <- grid.arrange(grobs = plot_list, ncol = 1)
    
    # Save plot
    ggsave("results/figures/treatment_analysis/metabolome_treatment_trajectory.png", p_combined,
           width = 12, height = 15, dpi = 300)
    
    print("Metabolome treatment trajectory plot generated successfully!")
  }
}

# ------------------------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------------------------
cat("\n=== treatment_heatmap.R Summary ===\n")
cat("\nGenerated results:\n")

if (exists("microbiome_treatment_results")) {
  cat("1. Microbiome treatment response analysis results\n")
  if (exists("significant_taxa") && length(significant_taxa) > 0) {
    cat("2. Significant treatment-responsive microbiome taxa\n")
    cat("3. Microbiome treatment heatmap\n")
    if (exists("p_combined")) {
      cat("4. Microbiome treatment trajectory plot\n")
    }
  } else {
    cat("2. Microbiome top variable heatmap (treatment analysis)\n")
  }
}

if (exists("metabolome_treatment_results")) {
  cat("5. Metabolome treatment response analysis results\n")
  if (exists("significant_metabolites") && length(significant_metabolites) > 0) {
    cat("6. Significant treatment-responsive metabolites\n")
    cat("7. Metabolome treatment heatmap\n")
    if (exists("p_combined")) {
      cat("8. Metabolome treatment trajectory plot\n")
    }
  } else {
    cat("6. Metabolome top variable heatmap (treatment analysis)\n")
  }
}

cat("\nResults saved to: results/treatment_analysis/\n")
cat("Figures saved to: results/figures/treatment_analysis/\n")
