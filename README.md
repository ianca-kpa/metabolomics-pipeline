# Metabolomics Pipeline

Modular **R pipeline** for analysis of **untargeted LC–MS metabolomics data** exported from **Compound Discoverer**.

The pipeline performs preprocessing, normalization, filtering, statistical analysis and visualization, generating reproducible outputs such as PCA plots, volcano plots, heatmaps and statistical tables.

---

## Main Features

- Modular R architecture
- Support for **Compound Discoverer exports**
- Weight normalization
- Missing value filtering
- QC-based filtering
- PCA visualization
- Statistical testing
- Volcano plots
- Heatmaps
- Export of results for **MetaboAnalyst**
- Organized output structure

---

## Project Structure

```text
metabolomics-pipeline/
│
├── config/
│   └── settings.example.R     # Example configuration file
│
├── data/                      # Input data (not tracked in Git)
│
├── R/                         # Pipeline modules
│   ├── 00_packages.R
│   ├── 01_validation.R
│   ├── 02_comparisons.R
│   ├── 03_helpers_io_log.R
│   ├── 04_metadata.R
│   ├── 05_features_assay.R
│   ├── 06_normalization_filters.R
│   ├── 07_duplicates.R
│   ├── 08_exports.R
│   ├── 09_pca.R
│   ├── 10_stats_volcano.R
│   ├── 11_heatmaps.R
│   └── 12_main_pipeline.R
│
├── run_pipeline.R             # Main script to run the pipeline
│
├── .gitignore
├── LICENSE
└── README.md
