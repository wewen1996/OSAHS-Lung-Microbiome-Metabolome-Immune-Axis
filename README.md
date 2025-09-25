# OSAHS-Microbiome-Metabolome-Analysis

This repository contains a comprehensive analysis pipeline for integrating microbiome and metabolome data in Obstructive Sleep Apnea-Hypopnea Syndrome (OSAHS) research. The pipeline includes data preprocessing, differential analysis, visualization, network analysis, multi-omics integration, statistical tables, pathway analysis, and treatment analysis.

## Project Structure

```
OSAHS-Microbiome-Metabolome-Analysis/
│
├── README.md                           # Project documentation
├── requirements.R                      # R package dependencies list
│
├── 01_data_preprocessing/              # Data preprocessing scripts
│   ├── load_and_clean_data.R
│   ├── phyloseq_object_creation.R
│   └── data_normalization.R
│
├── 02_differential_analysis/           # Differential analysis scripts
│   ├── microbiome_t_test.R
│   ├── metabolomics_OPLSDA.R
│   ├── age_adjusted_analysis.R
│   └── corncob_maaslin_ancombc.R
│
├── 03_visualization/                   # Visualization scripts
│   ├── heatmap_generation.R
│   ├── boxplot_visualization.R
│   ├── volcano_plots.R
│   ├── lollipop_charts.R
│   └── stacked_barplots.R
│
├── 04_network_analysis/                # Network analysis script
│   ├── correlation_networks.R
│   ├── WGCNA_analysis.R
│   └── graph_visualization.R
│
├── 05_multi_omics_integration/         # Multi-omics integration script
│   ├── mixOmics_analysis.R
│   ├── DIABLO_integration.R
│   └── cross_correlation.R
│
├── 06_statistical_tables/              # Statistical tables script
│   ├── baseline_table1.R
│   └── correlation_tables.R
│
├── 07_pathway_analysis/                # Pathway analysis script
│   ├── MetaboAnalystR_analysis.R
│   └── FELLA_enrichment.R
│
├── 08_treatment_analysis/              # Treatment analysis script
│   └── treatment_heatmap.R
│
├── data/                               # Data directory
│   ├── raw/                           # Raw data
│   ├── processed/                     # Processed data
│   └── metadata/                      # Metadata
│
└── results/                           # Results output directory
    ├── figures/                       # Generated figures
    ├── tables/                        # Statistical tables
    └── networks/                      # Network files
```

## Prerequisites

Before running the pipeline, ensure you have the following installed:

- R (version 4.0.0 or higher)
- Required R packages (see `requirements.R`)

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/OSAHS-Microbiome-Metabolome-Analysis.git
   cd OSAHS-Microbiome-Metabolome-Analysis
   ```

2. Install required R packages:
   ```
   Rscript requirements.R
   ```

## Data Preparation

1. Place your raw sequencing data in `data/raw/`
2. Place your metadata in `data/metadata/`
3. Ensure your data follows the required format:
   - Microbiome data: Fastq files or OTU table
   - Metabolome data: CSV file with metabolites as rows and samples as columns
   - Metadata: CSV file with samples as rows and variables as columns

## Running the Pipeline

The pipeline is designed to be run in order. Each script generates outputs that are used by subsequent scripts.

### 1. Data Preprocessing

```
cd 01_data_preprocessing
Rscript load_and_clean_data.R
Rscript phyloseq_object_creation.R
Rscript data_normalization.R
cd ..
```

### 2. Differential Analysis

```
cd 02_differential_analysis
Rscript microbiome_t_test.R
Rscript metabolomics_OPLSDA.R
Rscript age_adjusted_analysis.R
Rscript corncob_maaslin_ancombc.R
cd ..
```

### 3. Visualization

```
cd 03_visualization
Rscript heatmap_generation.R
Rscript boxplot_visualization.R
Rscript volcano_plots.R
Rscript lollipop_charts.R
Rscript stacked_barplots.R
cd ..
```

### 4. Network Analysis

```
cd 04_network_analysis
Rscript correlation_networks.R
Rscript WGCNA_analysis.R
Rscript graph_visualization.R
cd ..
```

### 5. Multi-omics Integration

```
cd 05_multi_omics_integration
Rscript mixOmics_analysis.R
Rscript DIABLO_integration.R
Rscript cross_correlation.R
cd ..
```

### 6. Statistical Tables

```
cd 06_statistical_tables
Rscript baseline_table1.R
Rscript correlation_tables.R
cd ..
```

### 7. Pathway Analysis

```
cd 07_pathway_analysis
Rscript MetaboAnalystR_analysis.R
Rscript FELLA_enrichment.R
cd ..
```

### 8. Treatment Analysis

```
cd 08_treatment_analysis
Rscript treatment_heatmap.R
cd ..
```

## Output

All results are saved in the `results/` directory:
- `results/figures/`: Generated figures in PNG format
- `results/tables/`: Statistical tables in CSV, HTML, and LaTeX formats
- `results/networks/`: Network analysis results

## Customization

This pipeline is designed for OSAHS microbiome-metabolome integration analysis but can be adapted for other studies:

1. Modify metadata columns in `06_statistical_tables/baseline_table1.R`
2. Adjust filtering parameters in data preprocessing scripts
3. Modify statistical methods in differential analysis script
4. Update pathway databases in pathway analysis script

## Troubleshooting

- If you encounter package installation issues, try installing them directly from CRAN or Bioconductor
- For memory issues with large datasets, consider increasing memory allocation in R
- If you get errors related to missing files, ensure you have run the previous script in the pipeline

## Contributing

Contributions are welcome! Please fork this repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- This pipeline uses several R packages from CRAN and Bioconductor
- Special thanks to the developers of phyloseq, microbiome, mixOmics, and other packages used in this pipeline
- This work was supported by [your funding source]
