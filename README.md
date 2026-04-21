# Metabolomics Pipeline

Modular **R pipeline** for analysis of **untargeted LC–MS/MS metabolomics data** exported from **Compound Discoverer**.

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
├── images/
│   ├── pca_example.png
│   ├── volcano_example.png
│   └── heatmap_example.png
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
```
------------------------------------------------------------------------

## Requirements

R version: **4.5.2** => https://cran.r-project.org/bin/windows/base/old/4.5.2/R-4.5.2-win.exe

Required R packages:

tidyverse\
readxl\
openxlsx\
pheatmap\
ggrepel\
stringi\
RColorBrewer

------------------------------------------------------------------------

## Configuration

Copy the example configuration file:
```
config/settings.example.R
```

and rename it to:
```
config/settings.R
```

Then edit the file and set:

```
cd_file_path <- "path/to/compound_discoverer_export.xlsx"
metadata_path <- "path/to/metadata.xlsx"
output_dir <- "path/to/output_folder"
```

------------------------------------------------------------------------

## Running the Pipeline

```
source("run_pipeline.R")
```

------------------------------------------------------------------------

## Output

The pipeline generates:

-   PCA plots
-   Volcano plots
-   Heatmaps
-   Statistical tables
-   MetaboAnalyst-ready datasets
-   Log files

Outputs are saved in the directory defined in `settings.R`.

------------------------------------------------------------------------

## Example Outputs

### PCA

![PCA](images/pca_example.png)

### Volcano Plot

![Volcano](images/volcano_example.png)

### Heatmap

![Heatmap](images/heatmap_example.png)

------------------------------------------------------------------------

## License

MIT License
