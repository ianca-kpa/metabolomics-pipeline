# ==========================================================
# SETTINGS FILE
# Copy this file to:
# config/settings.R
#
# Then edit the paths below before running the pipeline.
# ==========================================================

# Input files
cd_file_path <- "data/Compounds_POS_example.xlsx"
cd_sheet <- 1 # options: sheet index or sheet name

metadata_path <- "data/metadata_example.xlsx"
metadata_sheet <- 1 # options: sheet index or sheet name

# Output directory (root)
output_dir <- "output"

# Weight normalization
use_weight_normalization <- TRUE # options: TRUE/FALSE
stop_on_invalid_weight <- TRUE # options: TRUE/FALSE
invalid_weight_to_NA <- TRUE # options: TRUE/FALSE

# Filters
missing_exclusion_max_fraction <- 0.50 # options: 0..1 (set >=1 to disable)
presence_filter_min_fraction <- 0.00 # options: 0..1
impute_half_min <- TRUE # options: TRUE/FALSE

# QC RSD thresholds (variant creation)
rsd_thresholds <- c(20) # options: c(10,15,20,30,...) etc.
active_variant <- "QC_RSD20" # options: "BASE" or one of paste0("QC_RSD", rsd_thresholds)

# Low-variance filter
low_variance_filter_method <- "none" # options: "none" or "iqr"
low_variance_filter_fraction <- 0.20 # options: 0..1

# Transformation
log2_offset <- 1 # options: 0, 1, 0.5 ... (avoid 0 if you may have zeros)

# Statistical thresholds
alpha_sig <- 0.05
fc_cutoff_log2 <- 0

# Known-only filter
use_only_known <- TRUE # options: TRUE/FALSE

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
pca_scaling <- "pareto" # options: "none","pareto","autoscale"

make_heatmap_by_model <- TRUE
make_heatmap_by_model_sex <- TRUE
heatmap_top_n <- 80

heatmap_scale_method <- "zscore" # options: "none","zscore","pareto"
heatmap_order_samples_by_group <- TRUE
heatmap_cluster_distance <- "euclidean"
heatmap_cluster_method <- "ward.D2"

heatmap_palette_n <- 101
heatmap_breaks_symmetric <- TRUE
heatmap_breaks_limit <- 5

# Significant heatmaps
make_sig_heatmap_by_model <- TRUE
make_sig_heatmap_by_model_sex <- TRUE
make_sig_heatmap_FvsM_within_group <- TRUE
sig_heatmap_max_features <- 70
sig_heatmap_require_fc_cutoff <- TRUE

# =============================================================================
# Volcano plot settings
# =============================================================================

# Main volcano style:
# "classic" = categorical publication-like volcano
# "gradual" = continuous gradient volcano
# "both"    = export both styles
volcano_style <- "both" # options: "classic", "gradual", "both"

# Automatic axis scaling per plot
volcano_auto_axis <- TRUE # options: TRUE/FALSE

# Axis expansion factor when volcano_auto_axis = TRUE
volcano_axis_expand_mult <- 0.08 # options: numeric > 0

# Labels
volcano_add_labels <- TRUE # options: TRUE/FALSE
volcano_label_number <- 10 # options: integer >= 0

# Use NULL for automatic labels.
# Example:
# volcano_custom_labels <- c("L-Tryptophan", "Corticosterone")
volcano_custom_labels <- NULL # options: NULL or character vector

# Threshold lines
volcano_add_cutoff_lines <- TRUE # options: TRUE/FALSE

# -------------------------
# Classic style settings
# -------------------------
volcano_classic_point_size <- 2.5 # options: numeric > 0
volcano_classic_point_shape <- 21 # options: 21, 16, etc.

# IMPORTANT ORDER:
# c("Down", "Normal", "Up")
volcano_classic_fills <- c("#3B82F6", "#BDBDBD", "#EF4444")
volcano_classic_colors <- c("#1D4ED8", "#7A7A7A", "#B91C1C")
volcano_classic_legend_title <- "Regulation"

# -------------------------
# Gradual style settings
# -------------------------
volcano_gradual_point_shape <- 21
volcano_gradual_point_size_range <- c(1.5, 6)

# Default internal palettes
volcano_gradual_fills <- c("#39489F", "#39BBEC", "#F9ED36", "#F38466", "#B81F25")
volcano_gradual_colors <- c("#17194E", "#68BFE7", "#F9ED36", "#A22F27", "#211F1F")

# If TRUE, use RColorBrewer instead of the vectors above
volcano_gradual_use_RColorBrewer <- FALSE
volcano_gradual_brewer_palette <- "RdYlBu"
volcano_gradual_brewer_n <- 5
volcano_gradual_reverse_brewer <- TRUE

volcano_gradual_legend_title <- "Significance"
volcano_gradual_legend_breaks <- c(1, 2, 3, 4, 5)
volcano_gradual_legend_limits <- c(0, 5)

# Text sanitation
sanitize_names_for_exports <- TRUE
sanitize_mode <- "greek_latin_ascii" # options: "greek_latin_ascii" or "ascii_translit"

# Output control
minimal_outputs <- FALSE
