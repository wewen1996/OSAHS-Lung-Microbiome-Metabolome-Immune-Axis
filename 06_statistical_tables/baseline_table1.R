# baseline_table1.R
# Generate baseline characteristics table (Table 1) for study participants

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(tibble)
library(tableone)
library(gridExtra)
library(grid)
library(kableExtra)
library(formattable)

# Set working directory
setwd("/home/user/vibecoding/workspace/OSAHS-Microbiome-Metabolome-Analysis")

# Create output directory if it doesn't exist
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
# Load metadata
metadata <- read.csv("data/processed/metadata_processed.csv", row.names = 1, stringsAsFactors = FALSE)

print("Data loaded successfully!")
print(paste("Total samples:", nrow(metadata)))

# ------------------------------------------------------------------------------
# 2. Prepare data for Table 1
# ------------------------------------------------------------------------------
# Select variables for Table 1
# This should be customized based on your specific metadata
demographic_vars <- c("Age", "Gender", "BMI")
clinical_vars <- c("AHI", "ODI", "MinSaO2", "MeanSaO2", "ESS")
comorbidity_vars <- c("Hypertension", "Diabetes", "Dyslipidemia", "CHD", "Stroke")

# Combine variables
table1_vars <- c(demographic_vars, clinical_vars, comorbidity_vars)

# Check which variables are present in metadata
available_vars <- intersect(table1_vars, colnames(metadata))
missing_vars <- setdiff(table1_vars, colnames(metadata))

if (length(missing_vars) > 0) {
  print(paste("Warning: The following variables are not present in metadata:", paste(missing_vars, collapse = ", ")))
}

# Subset metadata to available variables
table1_data <- metadata[, available_vars, drop = FALSE]

# Add Group variable
if ("Group" %in% colnames(metadata)) {
  table1_data$Group <- metadata$Group
} else {
  stop("Group variable not found in metadata.")
}

# Convert character variables to factors
for (var in colnames(table1_data)) {
  if (is.character(table1_data[[var]]) && var != "Group") {
    table1_data[[var]] <- as.factor(table1_data[[var]])
  }
}

print("Data prepared for Table 1!")
print(paste("Variables included:", paste(available_vars, collapse = ", ")))

# ------------------------------------------------------------------------------
# 3. Create Table 1
# ------------------------------------------------------------------------------
cat("\n=== Creating Table 1 ===")

# Define variable labels (customize based on your data)
var_labels <- list(
  Age = "Age, years",
  Gender = "Gender",
  BMI = "Body Mass Index, kg/m²",
  AHI = "Apnea-Hypopnea Index, events/h",
  ODI = "Oxygen Desaturation Index, events/h",
  MinSaO2 = "Minimum Oxygen Saturation, %",
  MeanSaO2 = "Mean Oxygen Saturation, %",
  ESS = "Epworth Sleepiness Scale, score",
  Hypertension = "Hypertension",
  Diabetes = "Diabetes",
  Dyslipidemia = "Dyslipidemia",
  CHD = "Coronary Heart Disease",
  Stroke = "Stroke"
)

# Select only available labels
available_labels <- var_labels[available_vars]

# Create Table 1 using tableone package
table1 <- CreateTableOne(
  vars = available_vars,
  strata = "Group",
  data = table1_data,
  factorVars = setdiff(available_vars, c(demographic_vars[1:3], clinical_vars))
)

# Save Table 1 as text
write.csv(print(table1), "results/tables/table1_raw.csv", row.names = TRUE)

print("Table 1 created successfully!")

# ------------------------------------------------------------------------------
# 4. Format Table 1 for publication
# ------------------------------------------------------------------------------
cat("\n=== Formatting Table 1 for publication ===")

# Convert Table 1 to data frame
table1_df <- as.data.frame(table1)

# Clean up row names
rownames(table1_df) <- gsub("^\\d+\\.", "", rownames(table1_df))

# Add variable labels
table1_df <- table1_df %>%
  rownames_to_column("Variable") %>%
  mutate(Variable = ifelse(Variable %in% names(available_labels), available_labels[Variable], Variable)) %>%
  column_to_rownames("Variable")

# Save formatted Table 1
write.csv(table1_df, "results/tables/table1_formatted.csv", row.names = TRUE)

# Create LaTeX version of Table 1
table1_latex <- print(table1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# Save LaTeX table
write.table(table1_latex, "results/tables/table1.tex", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Create HTML version of Table 1 using kableExtra
table1_html <- kable(table1_df, format = "html", caption = "Baseline Characteristics of Study Participants") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  add_header_above(c(" " = 1, "Group" = (ncol(table1_df) - 1))) %>%
  footnote(general = "Continuous variables are presented as mean (standard deviation) or median [interquartile range] as appropriate. Categorical variables are presented as count (percentage).",
           general_title = "Note: ")

# Save HTML table
cat(table1_html, file = "results/tables/table1.html")

print("Formatted Table 1 saved in multiple formats!")

# ------------------------------------------------------------------------------
# 5. Create supplementary tables
# ------------------------------------------------------------------------------
cat("\n=== Creating supplementary tables ===")

# 5.1 Sample size by group
sample_size <- table(metadata$Group) %>%
  as.data.frame() %>%
  rename(Group = Var1, Count = Freq) %>%
  mutate(Percentage = Count / sum(Count) * 100)

write.csv(sample_size, "results/tables/sample_size.csv", row.names = FALSE)

# 5.2 Missing data summary
missing_data <- metadata %>%
  summarise_all(~sum(is.na(.))) %>%
  gather(Variable, Missing_Count) %>%
  mutate(
    Missing_Percentage = Missing_Count / nrow(metadata) * 100,
    Total_Samples = nrow(metadata)
  ) %>%
  arrange(desc(Missing_Count))

write.csv(missing_data, "results/tables/missing_data_summary.csv", row.names = FALSE)

# 5.3 Microbiome sequencing summary
if (file.exists("data/processed/phyloseq_object.rds")) {
  phyloseq_obj <- readRDS("data/processed/phyloseq_object.rds")
  
  # Sample sequencing depth
  seq_depth <- sample_sums(phyloseq_obj) %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    rename(Sequencing_Depth = ".") %>%
    left_join(metadata %>% rownames_to_column("Sample"), by = "Sample")
  
  # Sequencing depth summary by group
  seq_depth_summary <- seq_depth %>%
    group_by(Group) %>%
    summarise(
      Mean_Depth = mean(Sequencing_Depth),
      Median_Depth = median(Sequencing_Depth),
      SD_Depth = sd(Sequencing_Depth),
      Min_Depth = min(Sequencing_Depth),
      Max_Depth = max(Sequencing_Depth)
    ) %>%
    ungroup()
  
  write.csv(seq_depth_summary, "results/tables/sequencing_depth_summary.csv", row.names = FALSE)
  
  # Number of taxa at different levels
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxa_counts <- data.frame(
    Taxonomic_Level = tax_levels,
    Number_of_Taxa = sapply(tax_levels, function(level) {
      if (level %in% rank_names(phyloseq_obj)) {
        length(get_taxa_unique(phyloseq_obj, level))
      } else {
        NA
      }
    })
  ) %>%
    filter(!is.na(Number_of_Taxa))
  
  write.csv(taxa_counts, "results/tables/taxa_counts.csv", row.names = FALSE)
  
  print("Microbiome sequencing summary tables created successfully!")
}

print("Supplementary tables created successfully!")

# ------------------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------------------
cat("\n=== baseline_table1.R Summary ===\n")
cat("\nGenerated tables:\n")
cat("1. Table 1 (baseline characteristics) in raw, formatted, LaTeX, and HTML formats\n")
cat("2. Sample size by group\n")
cat("3. Missing data summary\n")
if (file.exists("data/processed/phyloseq_object.rds")) {
  cat("4. Sequencing depth summary by group\n")
  cat("5. Number of taxa at different taxonomic levels\n")
}

cat("\nAll tables saved to: results/tables/\n")
