# =============================================================================
# 01_validation.R
# Settings validation
# =============================================================================

validate_settings <- function() {
  stopifnot(pca_scaling %in% c("none","pareto","autoscale"))
  stopifnot(low_variance_filter_method %in% c("none","iqr"))
  stopifnot(duplicate_name_strategy %in% c("keep_separate","collapse_mean","collapse_sum","collapse_best_qc_rsd"))
  stopifnot(sanitize_mode %in% c("greek_latin_ascii","ascii_translit"))
  stopifnot(all(run_metrics %in% c("FDR","p_value")))
  stopifnot(all(heatmap_rank_metrics %in% c("FDR","p_value")))
  stopifnot(heatmap_scale_method %in% c("none","zscore","pareto"))
}