# ==========================================================
# SETTINGS FILE
# Copy this file to:
# config/settings.R
#
# Then edit the paths below before running the pipeline.
# ==========================================================

# Input files
cd_file_path <- "data/Compounds_POS_example.xlsx"
cd_sheet     <- 1  # options: sheet index or sheet name

metadata_path  <- "data/metadata_example.xlsx"
metadata_sheet <- 1  # options: sheet index or sheet name

# Output directory (root)
output_dir <- "output"

# Weight normalization
use_weight_normalization <- TRUE   # options: TRUE/FALSE
stop_on_invalid_weight   <- TRUE   # options: TRUE/FALSE
invalid_weight_to_NA     <- TRUE   # options: TRUE/FALSE

# Filters
missing_exclusion_max_fraction <- 0.50  # options: 0..1 (set >=1 to disable)
presence_filter_min_fraction   <- 0.00  # options: 0..1
impute_half_min <- TRUE                 # options: TRUE/FALSE

# QC RSD thresholds (variant creation)
rsd_thresholds <- c(20)          # options: c(10,15,20,30,...) etc.
active_variant <- "QC_RSD20"     # options: "BASE" or one of paste0("QC_RSD", rsd_thresholds)

# Low-variance filter
low_variance_filter_method   <- "none"  # options: "none" or "iqr"
low_variance_filter_fraction <- 0.20     # options: 0..1

# Transformation
log2_offset <- 1  # options: 0, 1, 0.5 ... (avoid 0 if you may have zeros)

# Statistical thresholds
alpha_sig <- 0.05
fc_cutoff_log2 <- 1

# Known-only filter
use_only_known <- TRUE  # options: TRUE/FALSE

# Duplicate metabolite handling
duplicate_name_strategy <- "collapse_best_qc_rsd"
# options: "keep_separate"; "collapse_mean"; "collapse_sum"; "collapse_best_qc_rsd"

# Duplicate rounding
dup_mz_digits <- 4
dup_rt_digits <- 2

# Metrics to run
run_metrics <- c("FDR", "p_value")
heatmap_rank_metrics <- c("FDR", "p_value")

# Exports
export_metaboanalyst_ready <- TRUE
save_stats_excel_per_model <- TRUE
make_volcano_plots <- TRUE

# PCA / Heatmaps
pca_scaling <- "pareto"  # options: "none","pareto","autoscale"

make_heatmap_by_model <- TRUE
make_heatmap_by_model_sex <- TRUE
heatmap_top_n <- 80

heatmap_scale_method <- "zscore"  # options: "none","zscore","pareto"
heatmap_order_samples_by_group <- TRUE
heatmap_cluster_distance <- "euclidean"
heatmap_cluster_method   <- "ward.D2"

heatmap_palette_n <- 101
heatmap_breaks_symmetric <- TRUE
heatmap_breaks_limit <- 5

# Significant heatmaps
make_sig_heatmap_by_model <- TRUE
make_sig_heatmap_by_model_sex <- TRUE
make_sig_heatmap_FvsM_within_group <- TRUE
sig_heatmap_max_features <- 70
sig_heatmap_require_fc_cutoff <- TRUE

# Text sanitation
sanitize_names_for_exports <- TRUE
sanitize_mode <- "greek_latin_ascii"  # options: "greek_latin_ascii" or "ascii_translit"

# Output control
minimal_outputs <- FALSE
