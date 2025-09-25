# lollipop_charts.R
# Generate lollipop charts for visualizing significant features

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
dir.create("results/figures/lollipop_charts", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load differential analysis results
microbiome_t_test_results <- read.csv("results/tables/differential_analysis/microbiome_t_test_results.csv", stringsAsFactors = FALSE)
metabolomics_oplsda_results <- read.csv("results/tables/differential_analysis/metabolomics_oplsda_results.csv", stringsAsFactors = FALSE)
microbiome_age_adjusted_results <- read.csv("results/tables/differential_analysis/microbiome_age_adjusted_results.csv", stringsAsFactors = FALSE)
metabolome_age_adjusted_results <- read.csv("results/tables/differential_analysis/metabolome_age_adjusted_results.csv", stringsAsFactors = FALSE)
ancombc_results <- read.csv("results/tables/differential_analysis/ancombc_results.csv", stringsAsFactors = FALSE)
maaslin2_results <- read.csv("results/tables/differential_analysis/maaslin2_results.csv", stringsAsFactors = FALSE)
corncob_results <- read.csv("results/tables/differential_analysis/corncob_results.csv", stringsAsFactors = FALSE)

# Load metabolite annotations
metabolite_annotations <- read.csv("data/processed/metabolite_annotations_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")

# ------------------------------------------------------------------------------
# 2. Function to create lollipop chart
# ------------------------------------------------------------------------------
create_lollipop_chart <- function(results, 
                                  x_col, 
                                  y_col, 
                                  feature_id_col = "Taxon",
                                  feature_name_col = NULL,
                                  color_col = NULL,
                                  title = "Lollipop Chart",
                                  xlab = "Feature",
                                  ylab = "Value",
                                  top_n = 20,
                                  horizontal = TRUE,
                                  color_palette = "viridis") {
  
  # Select top N features
  results_sorted <- results %>%
    arrange(desc(abs(!!sym(x_col)))) %>%
    head(top_n)
  
  # Add feature names if available
  if (!is.null(feature_name_col)) {
    results_sorted <- results_sorted %>%
      mutate(Feature_Label = ifelse(!is.na(!!sym(feature_name_col)), !!sym(feature_name_col), !!sym(feature_id_col)))
  } else {
    results_sorted <- results_sorted %>%
      mutate(Feature_Label = !!sym(feature_id_col))
  }
  
  # For microbiome data, create taxonomic labels
  if ("Phylum" %in% colnames(results_sorted) && "Genus" %in% colnames(results_sorted)) {
    results_sorted <- results_sorted %>%
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
  
  # Create factor for ordering
  results_sorted <- results_sorted %>%
    mutate(Feature_Label = factor(Feature_Label, levels = rev(unique(Feature_Label))))
  
  # Create lollipop chart
  if (horizontal) {
    p <- ggplot(results_sorted, aes(y = Feature_Label, x = !!sym(x_col)))
  } else {
    p <- ggplot(results_sorted, aes(x = Feature_Label, y = !!sym(x_col)))
  }
  
  # Add lines
  p <- p +
    geom_segment(
      aes(
        x = ifelse(horizontal, 0, as.numeric(Feature_Label)), 
        y = ifelse(horizontal, Feature_Label, 0),
        xend = !!sym(x_col), 
        yend = ifelse(horizontal, Feature_Label, !!sym(x_col))
      ),
      color = "gray50"
    )
  
  # Add points
  if (!is.null(color_col)) {
    if (color_palette == "viridis") {
      p <- p +
        geom_point(aes(color = !!sym(color_col)), size = 4) +
        scale_color_viridis_c()
    } else {
      p <- p +
        geom_point(aes(color = !!sym(color_col)), size = 4) +
        scale_color_gradient(low = "blue", high = "red")
    }
  } else {
    p <- p +
      geom_point(color = "red", size = 4)
  }
  
  # Add labels
  if (horizontal) {
    p <- p +
      geom_text(aes(label = round(!!sym(x_col), 2)), hjust = ifelse(results_sorted[[x_col]] > 0, 0.1, -0.1), size = 3)
  } else {
    p <- p +
      geom_text(aes(label = round(!!sym(x_col), 2)), vjust = ifelse(results_sorted[[x_col]] > 0, -0.5, 1.5), size = 3)
  }
  
  # Add theme and labels
  p <- p +
    labs(
      title = title,
      x = ifelse(horizontal, xlab, ""),
      y = ifelse(horizontal, "", ylab)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = ifelse(horizontal, element_text(), element_text(angle = 45, hjust = 1)),
      axis.text.y = ifelse(horizontal, element_text(size = 8), element_text()),
      legend.position = if (!is.null(color_col)) "bottom" else "none",
      legend.title = element_blank()
    )
  
  # Add zero line
  if (horizontal) {
    p <- p +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  } else {
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. Generate microbiome lollipop charts
# ------------------------------------------------------------------------------
cat("\n=== Generating microbiome lollipop charts ===\n")

# 3.1 Microbiome t-test lollipop chart
if (nrow(microbiome_t_test_results) > 0) {
  p_microbiome_t_test <- create_lollipop_chart(
    results = microbiome_t_test_results,
    x_col = "Log2_Fold_Change",
    y_col = "Adjusted_P_Value",
    feature_id_col = "Taxon",
    color_col = "Adjusted_P_Value",
    title = "Top Significant Microbiome Features (t-test)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "Adjusted P-value",
    top_n = 20,
    horizontal = TRUE,
    color_palette = "gradient"
  )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/microbiome_t_test_lollipop.png", p_microbiome_t_test,
         width = 12, height = 10, dpi = 300)
  
  print("Microbiome t-test lollipop chart generated successfully!")
} else {
  print("No microbiome t-test results available. Skipping lollipop chart generation.")
}

# 3.2 Age-adjusted microbiome lollipop chart
if (nrow(microbiome_age_adjusted_results) > 0) {
  p_microbiome_age_adjusted <- create_lollipop_chart(
    results = microbiome_age_adjusted_results,
    x_col = "Estimate",
    y_col = "Adjusted_P_Value",
    feature_id_col = "OTU",
    color_col = "Adjusted_P_Value",
    title = "Top Age-adjusted Microbiome Features",
    xlab = "Estimate (OSAHS vs Control)",
    ylab = "Adjusted P-value",
    top_n = 20,
    horizontal = TRUE,
    color_palette = "gradient"
  )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/microbiome_age_adjusted_lollipop.png", p_microbiome_age_adjusted,
         width = 12, height = 10, dpi = 300)
  
  print("Age-adjusted microbiome lollipop chart generated successfully!")
} else {
  print("No age-adjusted microbiome results available. Skipping lollipop chart generation.")
}

# 3.3 ANCOM-BC lollipop chart
if (nrow(ancombc_results) > 0) {
  p_ancombc <- create_lollipop_chart(
    results = ancombc_results,
    x_col = "lfc",
    y_col = "q_val",
    feature_id_col = "OTU",
    color_col = "q_val",
    title = "Top microbiome Features (ANCOM-BC)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "Adjusted P-value",
    top_n = 20,
    horizontal = TRUE,
    color_palette = "gradient"
  )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/ancombc_lollipop.png", p_ancombc,
         width = 12, height = 10, dpi = 300)
  
  print("ANCOM-BC lollipop chart generated successfully!")
} else {
  print("No ANCOM-BC results available. Skipping lollipop chart generation.")
}

# ------------------------------------------------------------------------------
# 4. Generate metabolome lollipop charts
# ------------------------------------------------------------------------------
cat("\n=== Generating metabolome lollipop charts ===\n")

# 4.1 Metabolomics OPLS-DA lollipop chart
if (nrow(metabolomics_oplsda_results) > 0) {
  p_metabolomics_oplsda <- create_lollipop_chart(
    results = metabolomics_oplsda_results,
    x_col = "Log2_Fold_Change",
    y_col = "VIP",
    feature_id_col = "Metabolite",
    feature_name_col = "Metabolite.Name",
    color_col = "VIP",
    title = "Top Metabolite Features (OPLS-DA)",
    xlab = "Log2 Fold Change (OSAHS vs Control)",
    ylab = "VIP Score",
    top_n = 20,
    horizontal = TRUE,
    color_palette = "viridis"
  )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/metabolomics_oplsda_lollipop.png", p_metabolomics_oplsda,
         width = 12, height = 10, dpi = 300)
  
  print("Metabolomics OPLS-DA lollipop chart generated successfully!")
} else {
  print("No metabolomics OPLS-DA results available. Skipping lollipop chart generation.")
}

# 4.2 Age-adjusted metabolome lollipop chart
if (nrow(metabolome_age_adjusted_results) > 0) {
  p_metabolome_age_adjusted <- create_lollipop_chart(
    results = metabolome_age_adjusted_results,
    x_col = "Estimate",
    y_col = "Adjusted_P_Value",
    feature_id_col = "Metabolite",
    feature_name_col = "Metabolite.Name",
    color_col = "Adjusted_P_Value",
    title = "Top Age-adjusted Metabolite Features",
    xlab = "Estimate (OSAHS vs Control)",
    ylab = "Adjusted P-value",
    top_n = 20,
    horizontal = TRUE,
    color_palette = "gradient"
  )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/metabolome_age_adjusted_lollipop.png", p_metabolome_age_adjusted,
         width = 12, height = 10, dpi = 300)
  
  print("Age-adjusted metabolome lollipop chart generated successfully!")
} else {
  print("No age-adjusted metabolome results available. Skipping lollipop chart generation.")
}

# ------------------------------------------------------------------------------
# 5. Generate combined lollipop charts
# ------------------------------------------------------------------------------
cat("\n=== Generating combined lollipop charts ===\n")

# 5.1 Combined microbiome lollipop chart (top features from all methods)
if (nrow(microbiome_t_test_results) > 0 && nrow(ancombc_results) > 0) {
  # Get top features from each method
  top_t_test <- microbiome_t_test_results %>%
    arrange(Adjusted_P_Value) %>%
    head(10) %>%
    select(Taxon, Log2_Fold_Change, Adjusted_P_Value) %>%
    mutate(Method = "t-test")
  
  top_ancombc <- ancombc_results %>%
    arrange(q_val) %>%
    head(10) %>%
    select(OTU, lfc, q_val) %>%
    rename(
      Taxon = OTU,
      Log2_Fold_Change = lfc,
      Adjusted_P_Value = q_val
    ) %>%
    mutate(Method = "ANCOM-BC")
  
  # Combine results
  combined_microbiome <- bind_rows(top_t_test, top_ancombc)
  
  # Add taxonomy information
  taxonomy <- tax_table(readRDS("data/processed/phyloseq_object.rds")) %>% as.data.frame()
  combined_microbiome <- combined_microbiome %>%
    left_join(taxonomy %>% rownames_to_column("Taxon"), by = "Taxon") %>%
    mutate(
      Taxonomic_Label = case_when(
        !is.na(Genus) ~ paste(Phylum, Genus, sep = "|"),
        !is.na(Family) ~ paste(Phylum, Family, sep = "|"),
        !is.na(Order) ~ paste(Phylum, Order, sep = "|"),
        !is.na(Class) ~ paste(Phylum, Class, sep = "|"),
        TRUE ~ Phylum
      )
    )
  
  # Create lollipop chart
  p_combined_microbiome <- ggplot(combined_microbiome, aes(y = Taxonomic_Label, x = Log2_Fold_Change, color = Method)) +
    geom_segment(
      aes(
        x = 0, 
        y = Taxonomic_Label,
        xend = Log2_Fold_Change, 
        yend = Taxonomic_Label
      ),
      color = "gray50"
    ) +
    geom_point(size = 4) +
    geom_text(aes(label = round(Log2_Fold_Change, 2)), hjust = ifelse(combined_microbiome$Log2_Fold_Change > 0, 0.1, -0.1), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "Top microbiome Features from Different Methods",
      x = "Log2 Fold Change (OSAHS vs Control)",
      y = ""
    ) +
    scale_color_manual(values = c("t-test" = "blue", "ANCOM-BC" = "red")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/combined_microbiome_lollipop.png", p_combined_microbiome,
         width = 14, height = 12, dpi = 300)
  
  print("Combined microbiome lollipop chart generated successfully!")
} else {
  print("Not enough data to generate combined microbiome lollipop chart.")
}

# 5.2 Combined metabolome lollipop chart (top features from all methods)
if (nrow(metabolomics_oplsda_results) > 0 && nrow(metabolome_age_adjusted_results) > 0) {
  # Get top features from each method
  top_oplsda <- metabolomics_oplsda_results %>%
    arrange(Adjusted_P_Value) %>%
    head(10) %>%
    select(Metabolite, Log2_Fold_Change, Adjusted_P_Value, VIP) %>%
    mutate(Method = "OPLS-DA")
  
  top_age_adjusted <- metabolome_age_adjusted_results %>%
    arrange(Adjusted_P_Value) %>%
    head(10) %>%
    select(Metabolite, Estimate, Adjusted_P_Value) %>%
    rename(Log2_Fold_Change = Estimate) %>%
    mutate(Method = "Age-adjusted", VIP = NA)
  
  # Combine results
  combined_metabolome <- bind_rows(top_oplsda, top_age_adjusted)
  
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
  
  # Create lollipop chart
  p_combined_metabolome <- ggplot(combined_metabolome, aes(y = Feature_Label, x = Log2_Fold_Change, color = Method)) +
    geom_segment(
      aes(
        x = 0, 
        y = Feature_Label,
        xend = Log2_Fold_Change, 
        yend = Feature_Label
      ),
      color = "gray50"
    ) +
    geom_point(size = 4) +
    geom_text(aes(label = round(Log2_Fold_Change, 2)), hjust = ifelse(combined_metabolome$Log2_Fold_Change > 0, 0.1, -0.1), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "Top Metabolite Features from Different Methods",
      x = "Log2 Fold Change / Estimate (OSAHS vs Control)",
      y = ""
    ) +
    scale_color_manual(values = c("OPLS-DA" = "green", "Age-adjusted" = "purple")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  # Save plot
  ggsave("results/figures/lollipop_charts/combined_metabolome_lollipop.png", p_combined_metabolome,
         width = 14, height = 12, dpi = 300)
  
  print("Combined metabolome lollipop chart generated successfully!")
} else {
  print("Not enough data to generate combined metabolome lollipop chart.")
}

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== lollipop_charts.R Summary ===\n")
cat("\nGenerated lollipop charts:\n")

if (nrow(microbiome_t_test_results) > 0) {
  cat("1. Microbiome t-test lollipop chart\n")
}

if (nrow(microbiome_age_adjusted_results) > 0) {
  cat("2. Age-adjusted microbiome lollipop chart\n")
}

if (nrow(ancombc_results) > 0) {
  cat("3. ANCOM-BC lollipop chart\n")
}

if (nrow(metabolomics_oplsda_results) > 0) {
  cat("4. Metabolomics OPLS-DA lollipop chart\n")
}

if (nrow(metabolome_age_adjusted_results) > 0) {
  cat("5. Age-adjusted metabolome lollipop chart\n")
}

if (nrow(microbiome_t_test_results) > 0 && nrow(ancombc_results) > 0) {
  cat("6. Combined microbiome lollipop chart\n")
}

if (nrow(metabolomics_oplsda_results) > 0 && nrow(metabolome_age_adjusted_results) > 0) {
  cat("7. Combined metabolome lollipop chart\n")
}

cat("\nAll lollipop charts saved to: results/figures/lollipop_charts/\n")
