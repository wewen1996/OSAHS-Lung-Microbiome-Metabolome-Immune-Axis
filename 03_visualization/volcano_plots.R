# volcano_plots.R
# Generate volcano plots for differential analysis results

# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(viridis)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/figures/volcano_plots", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load differential analysis results
microBIOME_t_test_results <- read.csv("results/tables/differential_analysis/microbiome_t_test_results.csv", stringsAsFactors = FALSE)
metabolomics_oplsda_results <- read.csv("results/tables/differential_analysis/metabolomics_oplsda_results.csv", stringsAsFactors = FALSE)
microBIOME_age_adjusted_results <- read.csv("results/tables/differential_analysis/microbiome_age_adjusted_results.csv", stringsAsfactor = FALSE)
metabolome_age_adjusted_results <- read.csv("results/tables/differential_analysis/metabolome_age_adjusted_results.csv", stringsAsfactor = FALSE)
ancombc_results <- read.csv("results/tables/differential_analysis/ancombc_results.csv", stringsAsfactor = FALSE)
maaslin2_results <- read.csv("results/tables/differential_analysis/maaslin2_results.csv", stringsAsfactor = FALSE)
corncob_results <- read.csv("results/tables/differential_analysis/corncob_results.csv", stringsAsfactor = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsfactor = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Function to create volcano plot
# ------------------------------------------------------------------------------
create_volcano_plot <- function(results, 
                                log2fc_col, 
                                pval_col, 
                                feature_id_col = "Taxon",
                                feature_name_col = NULL,
                                title = "Volcano Plot",
                                xlab = "Log2 Fold Change",
                                ylab = "-Log10(Adjusted P-value)",
                                log2fc_threshold = 1,
                                pval_threshold = 0.05,
                                top_n_labels = 10) {
  
  # Create volcano data
  volcano_data <- results %>%
    mutate(
      Neg_Log10_P = -log10(!!sym(pval_col)),
      Significant = (abs(!!sym(log2fc_col)) > log2fc_threshold) & (!!sym(pval_col) < pval_threshold),
      Upregulated = (!!sym(log2fc_col) > log2fc_threshold) & (!!sym(pval_col) < pval_threshold),
      Downregulated = (!!sym(log2fc_col) < -log2fc_threshold) & (!!sym(pval_col) < pval_threshold)
    )
  
  # Add feature names if available
  if (!is.null(feature_name_col)) {
    volcano_data <- volcano_data %>%
      mutate(Feature_Label = ifelse(!is.na(!!sym(feature_name_col)), !!sym(feature_name_col), !!sym(feature_id_col)))
  } else {
    volcano_data <- volcano_data %>%
      mutate(Feature_Label = !!sym(feature_id_col))
  }
  
  # For microbiome data, create taxonomic labels
  if ("Phylum" %in% colnames(volcano_data) && "Genus" %in% colnames(volcano_data)) {
    volcano_data <- volcano_data %>%
      mutate(
        Taxonomic_Label = case_when(
          !is.na(Genus) ~ paste(Phylum, Genus, sep = "|"),
          !is.na(Family) ~ paste(Phylum, Family, sep = "|"),
          !is.na(Order) ~ paste(Phylum, Order, sep = "|"),
          !is.na(Class) ~ paste(Phylum, Class, sep = "|"),
          TRUE ~ Phylum
        )
      ) %>%
      mutate(Feature_Label = Taxonomic_Label)
  }
  
  # Create volcano plot
  p <- ggplot(volcano_data, aes(x = !!sym(log2fc_col), y = Neg_Log10_P)) +
    geom_point(aes(color = Significant), alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "TRUE" = "red",
      "false" = "black"
    )) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add labels for top significant features
  if (sum(volcano_data$Significant) > 0) {
    top_features <- volcano_data %>%
      filter(Significant) %>%
      arrange(desc(Neg_Log10_P)) %>%
      head(top_n_labels)
    
    p <- p +
      geom_text_repel(
        data = top_features,
        aes(label = Feature_Label),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.5,
        segment.color = "gray50"
      )
  }
  
  # Add counts
  up_count <- sum(volcano_data$Upregulated, na.rm = TRUE)
  down_count <- sum(volcano_data$Downregulated, na.rm = TRUE)
  
  p <- p +
    annotate("text", x = max(volcano_data[[log2fc_col]], na.rm = TRUE) * 0.9, y = max(volcano_data$Neg_Log10_P, na.rm = TRUE) * 0.9,
             label = paste("Upregulated:", up_count, "\nDownregulated:", down_count),
             hjust = 1, vjust = 1, size = 4,
             bbox = list(boxstyle = "round,pad=0.3", fill = "white", alpha = 0.8))
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. Generate microbiome volcano plots
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome volcano plots ===\n")

# 3.1 Microbiome t-test volcano plot
if (nrow(microbiome_t_test_results) > 0) {
  p_microbiome_t_test <- create_volcano_plot(
    results = microbiome_t_test_results,
    log2fc_col = "Log2_Fold_Change",
    pval_col = "Adjusted_P_Value",
    feature_id_col = "Taxon",
    title = "Microbiome Differential Abundance (t-test)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 1,
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/microbiome_t_test_volcano.png", p_microbiome_t_test,
         width = 12, height = 8, dpi = 300)
  
  print("Microbiome t-test volcano plot generated successfully!")
} else {
  print("No microbiome t-test results available. Skipping volcano plot generation.")
}

# 3.2 Age-adjusted microbiome volcano plot
if (nrow(microbiome_age_adjusted_results) > 0) {
  p_microbiome_age_adjusted <- create_volcano_plot(
    results = microbiome_age_adjusted_results,
    log2fc_col = "Estimate",
    pval_col = "Adjusted_P_Value",
    feature_id_col = "OTU",
    title = "Age-adjusted Microbiome Differential Abundance",
    xlab = "Estimate (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 0.5,  # Using a smaller threshold for estimates
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/microbiome_age_adjusted_volcano.png", p_microbiome_age_adjusted,
         width = 12, height = 8, dpi = 300)
  
  print("Age-adjusted microbiome volcano plot generated successfully!")
} else {
  print("No age-adjusted microbiome results available. Skipping volcano plot generation.")
}

# 3.3 ANCOM-BC volcano plot
if (nrow(ancombc_results) > 0) {
  p_ancombc <- create_volcano_plot(
    results = ancombc_results,
    log2fc_col = "lfc",
    pval_col = "q_val",
    feature_id_col = "OTU",
    title = "Microbiome Differential Abundance (ANCOM-BC)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 1,
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/ancombc_volcano.png", p_ancombc,
         width = 12, height = 8, dpi = 300)
  
  print("ANCOM-BC volcano plot generated successfully!")
} else {
  print("No ANCOM-BC results available. Skipping volcano plot generation.")
}

# 3.4 Maaslin2 volcano plot
if (nrow(maaslin2_results) > 0) {
  p_maaslin2 <- create_volcano_plot(
    results = maaslin2_results,
    log2fc_col = "coef",
    pval_col = "qval",
    feature_id_col = "OTU",
    title = "Microbiome Differential Abundance (Maaslin2)",
    xlab = "Coefficient (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 0.5,  # Using a smaller threshold for coefficients
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/maaslin2_volcano.png", p_maaslin2,
         width = 12, height = 8, dpi = 300)
  
  print("Maaslin2 volcano plot generated successfully!")
} else {
  print("No Maaslin2 results available. Skipping volcano plot generation.")
}

# 3.5 corncob volcano plot
if (nrow(corncob_results) > 0) {
  # Calculate log2 fold change for corncob results
  corncob_results <- corncob_results %>%
    mutate(
      log2FoldChange = log2(estimate / estimate_null)
    )
  
  p_corncob <- create_volcano_plot(
    results = corncob_results,
    log2fc_col = "log2FoldChange",
    pval_col = "p_fdr",
    feature_id_col = "OTU",
    title = "Microbiome Differential Abundance (corncob)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 1,
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/corncob_volcano.png", p_corncob,
         width = 12, height = 8, dpi = 300)
  
  print("corncob volcano plot generated successfully!")
} else {
  print("No corncob results available. Skipping volcano plot generation.")
}

# ------------------------------------------------------------------------------
# 4. Generate metabolome volcano plots
# ------------------------------------------------------------------------------
cat("\n=== Generating metabolome volcano plots ===\n")

# 4.1 Metabolomics OPLS-DA volcano plot
if (nrow(metabolomics_oplsda_results) > 0) {
  p_metabolomics_oplsda <- create_volcano_plot(
    results = metabolomics_oplsda_results,
    log2fc_col = "Log2_Fold_Change",
    pval_col = "Adjusted_P_Value",
    feature_id_col = "Metabolite",
    feature_name_col = "Metabolite.Name",
    title = "Metabolome Differential Abundance (OPLS-DA)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 1,
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/metabolomics_oplsda_volcano.png", p_metabolomics_oplsda,
         width = 12, height = 8, dpi = 300)
  
  print("Metabolomics OPLS-DA volcano plot generated successfully!")
} else {
  print("No metabolomics OPLS-DA results available. Skipping volcano plot generation.")
}

# 4.2 Age-adjusted metabolome volcano plot
if (nrow(metabolome_age_adjusted_results) > 0) {
  p_metabolome_age_adjusted <- create_volcano_plot(
    results = metabolome_age_adjusted_results,
    log2fc_col = "Estimate",
    pval_col = "Adjusted_P_Value",
    feature_id_col = "Metabolite",
    feature_name_col = "Metabolite.Name",
    title = "Age-adjusted Metabolome Differential Abundance",
    xlab = "Estimate (OSAHS vs Control)",
    ylab = "-Log10(FDR)",
    log2fc_threshold = 0.5,  # Using a smaller threshold for estimates
    pval_threshold = 0.05,
    top_n_labels = 10
  )
  
  # Save plot
  ggsave("results/figures/volcano_plots/metabolome_age_adjusted_volcano.png", p_metabolome_age_adjusted,
         width = 12, height = 8, dpi = 300)
  
  print("Age-adjusted metabolome volcano plot generated successfully!")
} else {
  print("No age-adjusted metabolome results available. Skipping volcano plot generation.")
}

# ------------------------------------------------------------------------------
# 5. Generate combined volcano plots
# ------------------------------------------------------------------------------
cat("\n=== Generating combined volcano plots ===\n")

# 5.1 Combined microbiome volcano plots (t-test vs ANCOM-BC)
if (nrow(microbiome_t_test_results) > 0 && nrow(ancombc_results) > 0) {
  # Create combined data
  combined_microbiome <- microbiome_t_test_results %>%
    select(Taxon, Log2_Fold_Change, Adjusted_P_Value) %>%
    rename(
      log2fc_t_test = Log2_Fold_Change,
      pval_t_test = Adjusted_P_Value
    ) %>%
    left_join(
      ancombc_results %>%
        select(OTU, lfc, q_val) %>%
        rename(
          Taxon = OTU,
          log2fc_ancombc = lfc,
          pval_ancombc = q_val
        ),
      by = "Taxon"
    ) %>%
    filter(!is.na(log2fc_t_test) & !is.na(log2fc_ancombc)) %>%
    mutate(
      Significant_t_test = (abs(log2fc_t_test) > 1) & (pval_t_test < 0.05),
      Significant_ancombc = (abs(log2fc_ancombc) > 1) & (pval_ancombc < 0.05),
      Agreement = case_when(
        Significant_t_test & Significant_ancombc ~ "Both Significant",
        Significant_t_test & !Significant_ancombc ~ "Only t-test Significant",
        !Significant_t_test & Significant_ancombc ~ "Only ANCOM-BC Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create scatter plot of log2 fold changes
  p_combined_microbiome <- ggplot(combined_microbiome, aes(x = log2fc_t_test, y = log2fc_ancombc, color = Agreement)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    labs(
      title = "Microbiome Differential Abundance: t-test vs ANCOM-BC",
      x = "Log2 Fold Change (t-test)",
      y = "Log2 Fold Change (ANCOM-BC)"
    ) +
    scale_color_manual(values = c(
      "Both Significant" = "red",
      "Only t-test Significant" = "blue",
      "Only ANCOM-BC Significant" = "green",
      "Not Significant" = "black"
    )) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add counts
  agreement_counts <- table(combined_microbiome$Agreement)
  p_combined_microbiome <- p_combined_microbiome +
    annotate("text", x = max(combined_microbiome$log2fc_t_test) * 0.9, y = min(combined_microbiome$log2fc_ancombc) * 0.9,
             label = paste(paste(names(agreement_counts), agreement_counts, sep = ": "), collapse = "\n"),
             hjust = 1, vjust = 0, size = 3,
             bbox = list(boxstyle = "round,pad=0.3", fill = "white", alpha = 0.8))
  
  # Save plot
  ggsave("results/figures/volcano_plots/combined_microbiome_volcano.png", p_combined_microbiome,
         width = 12, height = 10, dpi = 300)
  
  print("Combined microbiome volcano plot generated successfully!")
} else {
  print("Not enough data to generate combined microbiome volcano plot.")
}

# 5.2 Combined metabolome volcano plots (OPLS-DA vs age-adjusted)
if (nrow(metabolomics_oplsda_results) > 0 && nrow(metabolome_age_adjusted_results) > 0) {
  # Create combined data
  combined_metabolome <- metabolomics_oplsda_results %>%
    select(Metabolite, Log2_Fold_Change, Adjusted_P_Value) %>%
    rename(
      log2fc_oplsda = Log2_Fold_Change,
      pval_oplsda = Adjusted_P_Value
    ) %>%
    left_join(
      metabolome_age_adjusted_results %>%
        select(Metabolite, Estimate, Adjusted_P_Value) %>%
        rename(
          log2fc_age_adjusted = Estimate,
          pval_age_adjusted = Adjusted_P_Value
        ),
      by = "Metabolite"
    ) %>%
    filter(!is.na(log2fc_oplsda) & !is.na(log2fc_age_adjusted)) %>%
    mutate(
      Significant_oplsda = (abs(log2fc_oplsda) > 1) & (pval_oplsda < 0.05),
      Significant_age_adjusted = (abs(log2fc_age_adjusted) > 0.5) & (pval_age_adjusted < 0.05),
      Agreement = case_when(
        Significant_oplsda & Significant_age_adjusted ~ "Both Significant",
        Significant_oplsda & !Significant_age_adjusted ~ "Only OPLS-DA Significant",
        !Significant_oplsda & Significant_age_adjusted ~ "Only Age-adjusted Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # Add metabolite names
  if (!is.null(metabolite_annotations) && "Metabolite.Name" %in% colnames(metabolite_annotations)) {
    combined_metabolome <- combined_metabolome %>%
      left_join(
        metabolite_annotations %>%
          select(Metabolite.Name) %>%
          rownames_to_column("Metabolite"),
        by = "Metabolite"
      ) %>%
      mutate(
        Feature_Label = ifelse(!is.na(Metabolite.Name), Metabolite.Name, Metabolite)
      )
  } else {
    combined_metabolome <- combined_metabolome %>%
      mutate(
        Feature_Label = Metabolite
      )
  }
  
  # Create scatter plot of log2 fold changes
  p_combined_metabolome <- ggplot(combined_metabolome, aes(x = log2fc_oplsda, y = log2fc_age_adjusted, color = Agreement)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
    labs(
      title = "Metabolome Differential Abundance: OPLS-DA vs Age-adjusted",
      x = "Log2 Fold Change (OPLS-DA)",
      y = "Estimate (Age-adjusted)"
    ) +
    scale_color_manual(values = c(
      "Both Significant" = "red",
      "Only OPLS-DA Significant" = "blue",
      "Only Age-adjusted Significant" = "green",
      "Not Significant" = "black"
    )) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Add labels for top significant features
  top_features <- combined_metabolome %>%
    filter(Agreement == "Both Significant") %>%
    arrange(desc(abs(log2fc_oplsda))) %>%
    head(10)
  
  if (nrow(top_features) > 0) {
    p_combined_metabolome <- p_combined_metabolome +
      geom_text_repel(
        data = top_features,
        aes(label = Feature_Label),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.5,
        segment.color = "gray50"
      )
  }
  
  # Add counts
  agreement_counts <- table(combined_metabolome$Agreement)
  p_combined_metabolome <- p_combined_metabolome +
    annotate("text", x = max(combined_metabolome$log2fc_oplsda) * 0.9, y = min(combined_metabolome$log2fc_age_adjusted) * 0.9,
             label = paste(paste(names(agreement_counts), agreement_counts, sep = ": "), collapse = "\n"),
             hjust = 1, vjust = 0, size = 3,
             bbox = list(boxstyle = "round,pad=0.3", fill = "white", alpha = 0.8))
  
  # Save plot
  ggsave("results/figures/volcano_plots/combined_metabolome_volcano.png", p_combined_metabolome,
         width = 14, height = 10, dpi = 300)
  
  print("Combined metabolome volcano plot generated successfully!")
} else {
  print("Not enough data to generate combined metabolome volcano plot.")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== volcano_plots.R Summary ===\n")
cat("\nGenerated volcano plots:\n")

if (nrow(microbiome_t_test_results) > 0) {
  cat("1. Microbiome t-test volcano plot\n")
}

if (nrow(microbiome_age_adjusted_results) > 0) {
  cat("2. Age-adjusted microbiome volcano plot\n")
}

if (nrow(ancombc_results) > 0) {
  cat("3. ANCOM-BC volcano plot\n")
}

if (nrow(maaslin2_results) > 0) {
  cat("4. Maaslin2 volcano plot\n")
}

if (nrow(corncob_results) > 0) {
  cat("5. corncob volcano plot\n")
}

if (nrow(metabolomics_oplsda_results) > 0) {
  cat("6. Metabolomics OPLS-DA volcano plot\n")
}

if (nrow(metabolome_age_adjusted_results) > 0) {
  cat("7. Age-adjusted metabolome volcano plot\n")
}

if (nrow(microbiome_t_test_results) > 0 && nrow(ancombc_results) > 0) {
  cat("8. Combined microbiome volcano plot (t-test vs ANCOM-BC)\n")
}

if (nrow(metabolomics_oplsda_results) > 0 && nrow(metabolome_age_adjusted_results) > 0) {
  cat("9. Combined metabolome volcano plot (OPLS-DA vs age-adjusted)\n")
}

cat("\nAll volcano plots saved to: results/figures/volcano_plots/\n")
