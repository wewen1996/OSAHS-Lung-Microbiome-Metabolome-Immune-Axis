# correlation_tables.R
# Generate correlation tables for microbiome and metabolome data

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(corrplot)
library(Hmisc)
library(reshape2)
library(gridExtra)
library(grid)
library(kableExtra)
library(formattable)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/tables/correlation_tables", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load normalizedized data
phyloseq_css <- readRDS("data/processed/normalized/phyloseq_css.rds")
metabolome_data <- read.csv("data/processed/normalized/metabolome_pareto_scaled.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsfactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Prepare data for correlation analysis
# ------------------------------------------------------------------------------
# Ensure samples are in the same order
common_samples <- intersect(sample_names(phyloseq_css), colnames(metabolome_data))

if (length(common_samples) < 3) {
  stop("Not enough common samples for correlation analysis.")
}

# Subset data to common samples
phyloseq_css <- prune_samples(common_samples, phyloseq_css)
metabolome_data <- metabolome_data[, common_samples]
metadata <- metadata[common_samples, ]

print(paste("Common samples:", length(common_samples)))

# Get microbiome data matrix
microbiome_data <- otu_table(phyloseq_css) %>% as.matrix()

# Filter low abundance/ variance features to reduce dimensionality
# Microbiome filtering
min_abundance <- 0.001  # 0.1% relative abundance
min_prevalence <- 0.2   # Present in at least 20% of samples

# Calculate relative abundance
microbiome_rel <- t(t(microbiome_data) / colSums(microbiome_data))

# Filter taxa
taxa_to_keep <- rowSums(microbiome_rel > min_abundance) >= (min_prevalence * ncol(microbiome_rel))
microbiome_filtered <- microbiome_data[taxa_to_keep, ]

print(paste("Microbiome taxa after filtering:", sum(taxa_to_keep)))
print(paste("Original microbiome taxa:", nrow(microbiome_data)))

# Metabolome filtering
min_variance <- quantile(apply(metabolome_data, 1, var), 0.25)  # Keep top 75% most variable metabolites
metabolite_variance <- apply(metabolome_data, 1, var)
metabolites_to_keep <- metabolite_variance >= min_variance
metabolome_filtered <- metabolome_data[metabolites_to_keep, ]

print(paste("Metabolites after filtering:", sum(metabolites_to_keep)))
print(paste("Original metabolites:", nrow(metabolome_data)))

# Transpose data so that samples are rows and features are columns
microbiome_filtered_t <- t(microbiome_filtered)
metabolome_filtered_t <- t(metabolome_filtered)

print("Data prepared for correlation analysis!")

# ------------------------------------------------------------------------------
# 3. Calculate correlations between microbiome and clinical variables
# ------------------------------------------------------------------------------
cat("\n=== Calculating correlation between microbiome and clinical variables ===\n")

# Select clinical variables for correlation analysis
clinical_vars <- c("Age", "BMI", "AHI", "ODI", "MinSaO2", "MeanSaO2", "ESS")
clinical_vars <- intersect(clinical_vars, colnames(metadata))

if (length(clinical_vars) > 0) {
  # Get clinical data
  clinical_data <- metadata[, clinical_vars, drop = FALSE]
  
  # Remove samples with missing values
  complete_samples <- complete.cases(clinical_data)
  clinical_data_complete <- clinical_data[complete_samples, ]
  microbiome_filtered_t_complete <- microbiome_filtered_t[complete_samples, ]
  
  print(paste("Sample with complete clinical data:", nrow(clinical_data_complete)))
  
  if (nrow(clinical_data_complete) >= 3) {
    # Function to calculate correlation between microbiome and clinical variables
    calculate_microbiome_clinical_correlation <- function(microbiome_data, clinical_data, method = "spearman") {
      
      # Initialize results data frame
      results <- data.frame(
        Taxon = character(),
        Clinical_Variable = character(),
        Correlation = numeric(),
        P_Value = numeric(),
        Adjusted_P_Value = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Calculate correlation between each taxon and each clinical variable
      for (i in 1:ncol(microbiome_data)) {
        for (j in 1:ncol(clinical_data)) {
          # Get feature names
          taxon <- colnames(microbiome_data)[i]
          clinical_var <- colnames(clinical_data)[j]
          
          # Calculate correlation
          correlation_test <- cor.test(microbiome_data[, i], clinical_data[, j], method = method)
          
          # Add to results
          results <- rbind(results, data.frame(
            Taxon = taxon,
            Clinical_Variable = clinical_var,
            Correlation = correlation_test$estimate,
            P_Value = correlation_test$p.value,
            Adjusted_P_Value = NA,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Apply multiple testing correction
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
      
      # Order by absolute correlation
      results <- results %>% arrange(desc(abs(Correlation)))
      
      return(results)
    }
    
    # Calculate correlation
    microbiome_clinical_correlation <- calculate_microbiome_clinical_correlation(
      microbiome_data = microbiome_filtered_t_complete,
      clinical_data = clinical_data_complete,
      method = "spearman"
    )
    
    # Save results
    write.csv(microbiome_clinical_correlation, "results/tables/correlation_tables/microbiome_clinical_correlation.csv", row.names = FALSE)
    
    # Add taxonomy information
    taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
    microbiome_clinical_correlation <- microbiome_clinical_correlation %>%
      left_join(taxonomy %>% rownames_to_column("Taxon"), by = "Taxon")
    
    # Identify significant correlation
    significant_correlation <- microbiome_clinical_correlation %>%
      filter(Adjusted_P_Value < 0.05) %>%
      arrange(desc(abs(Correlation)))
    
    # Save significant correlation
    write.csv(significant_correlation, "results/tables/correlation_tables/microbiome_clinical_significant_correlation.csv", row.names = FALSE)
    
    print(paste("Significant microbiome-clinical correlation found:", nrow(significant_correlation)))
    
    if (nrow(significant_correlation) > 0) {
      # Create formatted table for publication
      formatted_table <- significant_correlation %>%
        select(Taxon, Phylum, Genus, Clinical_Variable, Correlation, Adjusted_P_Value) %>%
        mutate(
          Correlation = round(Correlation, 3),
          Adjusted_P_Value = format.pval(Adjusted_P_Value, digits = 3),
          Taxon_Label = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), Taxon)
        ) %>%
        select(-Taxon, -Phylum, -Genus) %>%
        rename(Taxon = Taxon_Label) %>%
        arrange(desc(abs(Correlation)))
      
      # Save formatted table
      write.csv(formatted_table, "results/tables/correlation_tables/microbiome_clinical_correlation_formatted.csv", row.names = FALSE)
      
      # Create HTML table
      html_table <- kable(formatted_table, format = "html", caption = "Significant Correlation Between Microbiome and Clinical Variables") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(5, color = ifelse(formatted_table$Adjusted_P_Value < 0.01, "red", "black"))
      
      cat(html_table, file = "results/tables/correlation_tables/microbiome_clinical_correlation.html")
      
      print("Microbiome-clinical correlation tables created successfully!")
    }
  } else {
    print("Not enough samples with complete clinical data for correlation analysis.")
  }
} else {
  print("No clinical variables found for correlation analysis.")
}

# ------------------------------------------------------------------------------
# 4. Calculate correlation between metabolome and clinical variables
# ------------------------------------------------------------------------------
cat("\n=== Calculating correlation between metabolome and clinical variables ===\n")

if (length(clinical_vars) > 0) {
  # Get clinical data
  clinical_data <- metadata[, clinical_vars, drop = FALSE]
  
  # Remove samples with missing values
  complete_samples <- complete.cases(clinical_data)
  clinical_data_complete <- clinical_data[complete_samples, ]
  metabolome_filtered_t_complete <- metabolome_filtered_t[complete_samples, ]
  
  print(paste("Samples with complete clinical data:", nrow(clinical_data_complete)))
  
  if (nrow(clinical_data_complete) >= 3) {
    # Function to calculate correlation between metabolome and clinical variables
    calculate_metabolome_clinical_correlation <- function(metabolome_data, clinical_data, method = "spearman") {
      
      # Initialize results data frame
      results <- data.frame(
        Metabolite = character(),
        Clinical_Variable = character(),
        Correlation = numeric(),
        P_Value = numeric(),
        Adjusted_P_Value = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Calculate correlation between each metabolite and each clinical variable
      for (i in 1:ncol(metabolome_data)) {
        for (j in 1:ncol(clinical_data)) {
          # Get feature names
          metabolite <- colnames(metabolome_data)[i]
          clinical_var <- colnames(clinical_data)[j]
          
          # Calculate correlation
          correlation_test <- cor.test(metabolome_data[, i], clinical_data[, j], method = method)
          
          # Add to results
          results <- rbind(results, data.frame(
            Metabolite = metabolite,
            Clinical_Variable = clinical_var,
            Correlation = correlation_test$estimate,
            P_Value = correlation_test$p.value,
            Adjusted_P_Value = NA,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Apply multiple testing correction
      results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")
      
      # Order by absolute correlation
      results <- results %>% arrange(desc(abs(Correlation)))
      
      return(results)
    }
    
    # Calculate correlation
    metabolome_clinical_correlation <- calculate_metabolome_clinical_correlation(
      metabolome_data = metabolome_filtered_t_complete,
      clinical_data = clinical_data_complete,
      method = "spearman"
    )
    
    # Save results
    write.csv(metabolome_clinical_correlation, "results/tables/correlation_tables/metabolome_clinical_correlation.csv", row.names = FALSE)
    
    # Add metabolite annotations
    if (!is.null(metabolite_annotations)) {
      metabolome_clinical_correlation <- metabolome_clinical_correlation %>%
        left_join(metabolite_annotations %>% rownames_to_column("Metabolite"), by = "Metabolite")
    }
    
    # Identify significant correlation
    significant_correlation <- metabolome_clinical_correlation %>%
      filter(Adjusted_P_Value < 0.05) %>%
      arrange(desc(abs(Correlation)))
    
    # Save significant correlation
    write.csv(significant_correlation, "results/tables/correlation_tables/metabolome_clinical_significant_correlation.csv", row.names = FALSE)
    
    print(paste("Significant metabolome-clinical correlation found:", nrow(significant_correlation)))
    
    if (nrow(significant_correlation) > 0) {
      # Create formatted table for publication
      formatted_table <- significant_correlation %>%
        select(Metabolite, Metabolite.Name, Class, Clinical_Variable, Correlation, Adjusted_P_Value) %>%
        mutate(
          Correlation = round(Correlation, 3),
          Adjusted_P_Value = format.pval(Adjusted_P_Value, digits = 3),
          Metabolite_Label = ifelse(!is.na(Metabolite.Name), paste(Metabolite.Name, "(", Metabolite, ")", sep = ""), Metabolite)
        ) %>%
        select(-Metabolite, -Metabolite.Name) %>%
        rename(Metabolite = Metabolite_Label) %>%
        arrange(desc(abs(Correlation)))
      
      # Save formatted table
      write.csv(formatted_table, "results/tables/correlation_tables/metabolome_clinical_correlation_formatted.csv", row.names = FALSE)
      
      # Create HTML table
      html_table <- kable(formatted_table, format = "html", caption = "Significant Correlation Between Metabolome and Clinical Variables") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(5, color = ifelse(formatted_table$Adjusted_P_Value < 0.01, "red", "black"))
      
      cat(html_table, file = "results/tables/correlation_tables/metabolome_clinical_correlation.html")
      
      print("Metabolome-clinical correlation tables created successfully!")
    }
  } else {
    print("Not enough samples with complete clinical data for correlation analysis.")
  }
} else {
  print("No clinical variables found for correlation analysis.")
}

# ------------------------------------------------------------------------------
# 5. Calculate correlation between microbiome and metabolome
# ------------------------------------------------------------------------------
cat("\n=== Calculating correlation between microbiome and metabolome ===\n")

# Function to calculate pairwise correlations between two datasets
calculate_pairwise_correlations <- function(data1, data2, method = "spearman", p_adjust_method = "fdr") {
  
  # Initialize results data frame
  results <- data.frame(
    Feature1 = character(),
    Feature2 = character(),
    Correlation = numeric(),
    P_Value = numeric(),
    Adjusted_P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Calculate correlation between each feature in data1 and each feature in data2
  for (i in 1:ncol(data1)) {
    for (j in 1:ncol(data2)) {
      # Get feature names
      feature1 <- colnames(data1)[i]
      feature2 <- colnames(data2)[j]
      
      # Calculate correlation
      correlation_test <- cor.test(data1[, i], data2[, j], method = method)
      
      # Add to results
      results <- rbind(results, data.frame(
        Feature1 = feature1,
        Feature2 = feature2,
        Correlation = correlation_test$estimate,
        P_Value = correlation_test$p.value,
        Adjusted_P_Value = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply multiple testing correction
  results$Adjusted_P_Value <- p.adjust(results$P_Value, method = p_adjust_method)
  
  # Order by absolute correlation
  results <- results %>% arrange(desc(abs(Correlation)))
  
  return(results)
}

# Calculate pairwise correlations between microbiome and metabolome
pairwise_correlations <- calculate_pairwise_correlations(
  data1 = microbiome_filtered_t,
  data2 = metabolome_filtered_t,
  method = "spearman",
  p_adjust_method = "fdr"
)

# Save pairwise correlations
write.csv(pairwise_correlations, "results/tables/correlation_tables/microbiome_metabolome_pairwise_correlations.csv", row.names = FALSE)

# Identify significant correlations
correlation_threshold <- 0.6
p_value_threshold <- 0.05

significant_correlations <- pairwise_correlations %>%
  filter(abs(Correlation) >= correlation_threshold & Adjusted_P_Value < p_value_threshold)

# Save significant correlations
write.csv(significant_correlations, "results/tables/correlation_tables/microbiome_metabolome_significant_correlations.csv", row.names = FALSE)

print(paste("Significant microbiome-metabolome correlations found (|r| >=", correlation_threshold, "and FDR <", p_value_threshold, "):", nrow(significant_correlations)))

if (nrow(significant_correlations) > 0) {
  # Add taxonomy and annotation information
  taxonomy <- tax_table(phyloseq_css) %>% as.data.frame()
  
  significant_correlations <- significant_correlations %>%
    left_join(taxonomy %>% rownames_to_column("Feature1"), by = "Feature1") %>%
    rename(Taxon = Feature1)
  
  if (!is.null(metabolite_annotations)) {
    significant_correlations <- significant_correlations %>%
      left_join(metabolite_annotations %>% rownames_to_column("Feature2"), by = "Feature2") %>%
      rename(Metabolite = Feature2)
  }
  
  # Create formatted table for publication
  formatted_table <- significant_correlations %>%
    select(Taxon, Phylum, Genus, Metabolite, Metabolite.Name, Class, Correlation, Adjusted_P_Value) %>%
    mutate(
      Correlation = round(Correlation, 3),
      Adjusted_P_Value = format.pval(Adjusted_P_Value, digits = 3),
      Taxon_Label = ifelse(!is.na(Genus), paste(Phylum, Genus, sep = "|"), Taxon),
      Metabolite_Label = ifelse(!is.na(Metabolite.Name), paste(Metabolite.Name, "(", Metabolite, ")", sep = ""), Metabolite)
    ) %>%
    select(-Taxon, -Phylum, -Genus, -Metabolite, -Metabolite.Name) %>%
    rename(Taxon = Taxon_Label, Metabolite = Metabolite_Label) %>%
    arrange(desc(abs(Correlation)))
  
  # Save formatted table
  write.csv(formatted_table, "results/tables/correlation_tables/microbiome_metabolome_correlation_formatted.csv", row.names = FALSE)
  
  # Create HTML table
  html_table <- kable(formatted_table, format = "html", caption = "Significant Correlation Between Microbiome and Metabolome") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    column_spec(7, color = ifelse(formatted_table$Adjusted_P_Value < 0.01, "red", "black"))
  
  cat(html_table, file = "results/tables/correlation_tables/microbiome_metabolome_correlation.html")
  
  print("Microbiome-metabolome correlation tables created successfully!")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== correlation_tables.R Summary ===\n")
cat("\nGenerated tables:\n")

if (length(clinical_vars) > 0 && exists("microbiome_clinical_correlation")) {
  cat("1. Microbiome-clinical variable correlation table\n")
  if (nrow(significant_correlation) > 0) {
    cat("2. Significant microbiome-clinical variable correlation table (formatted for publication)\n")
  }
}

if (length(clinical_vars) > 0 && exists("metabolome_clinical_correlation")) {
  cat("3. Metabolome-clinical variable correlation table\n")
  if (exists("significant_correlation") && nrow(significant_correlation) > 0) {
    cat("4. Significant metabolome-clinical variable correlation table (formatted for publication)\n")
  }
}

cat("5. Microbiome-metabolome pairwise correlation table\n")
if (nrow(significant_correlations) > 0) {
  cat("6. Significant microbiome-metabolome correlation table (formatted for publication)\n")
}

cat("\nAll correlation tables saved to: results/tables/correlation_tables/\n")
