# =============================================================================
# 01_validation.R
# Settings validation
# =============================================================================

validate_settings <- function() {
  # ---------------------------------------------------------------------------
  # PCA
  # ---------------------------------------------------------------------------
  stopifnot(pca_scaling %in% c("none", "pareto", "autoscale"))

  # ---------------------------------------------------------------------------
  # Filters
  # ---------------------------------------------------------------------------
  stopifnot(low_variance_filter_method %in% c("none", "iqr"))

  # ---------------------------------------------------------------------------
  # Duplicate handling
  # ---------------------------------------------------------------------------
  stopifnot(
    duplicate_name_strategy %in% c(
      "keep_separate",
      "collapse_mean",
      "collapse_sum",
      "collapse_best_qc_rsd"
    )
  )

  # ---------------------------------------------------------------------------
  # Name sanitation
  # ---------------------------------------------------------------------------
  stopifnot(
    sanitize_mode %in% c(
      "greek_latin_ascii",
      "ascii_translit"
    )
  )

  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------
  stopifnot(all(run_metrics %in% c("FDR", "p_value")))
  stopifnot(all(heatmap_rank_metrics %in% c("FDR", "p_value")))

  # ---------------------------------------------------------------------------
  # Heatmap
  # ---------------------------------------------------------------------------
  stopifnot(heatmap_scale_method %in% c("none", "zscore", "pareto"))

  # ---------------------------------------------------------------------------
  # Volcano main style
  # ---------------------------------------------------------------------------
  stopifnot(
    volcano_style %in% c(
      "classic",
      "gradual",
      "both"
    )
  )

  # ---------------------------------------------------------------------------
  # Volcano classic legend
  # ---------------------------------------------------------------------------
  stopifnot(is.character(volcano_classic_legend_title))
  stopifnot(length(volcano_classic_legend_title) == 1)

  # classic colors must be exactly 3
  stopifnot(length(volcano_classic_fills) == 3)
  stopifnot(length(volcano_classic_colors) == 3)

  # ---------------------------------------------------------------------------
  # Volcano gradual legend
  # ---------------------------------------------------------------------------
  stopifnot(is.character(volcano_gradual_legend_title))
  stopifnot(length(volcano_gradual_legend_title) == 1)

  stopifnot(is.numeric(volcano_gradual_legend_breaks))
  stopifnot(length(volcano_gradual_legend_breaks) >= 2)

  stopifnot(is.numeric(volcano_gradual_legend_limits))
  stopifnot(length(volcano_gradual_legend_limits) == 2)

  # ---------------------------------------------------------------------------
  # Volcano gradual palette
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_gradual_use_RColorBrewer))

  if (!volcano_gradual_use_RColorBrewer) {
    stopifnot(length(volcano_gradual_fills) >= 3)
    stopifnot(length(volcano_gradual_colors) >= 3)
  }

  # ---------------------------------------------------------------------------
  # Volcano labels
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_add_labels))
  stopifnot(is.numeric(volcano_label_number))
  stopifnot(volcano_label_number >= 0)

  # ---------------------------------------------------------------------------
  # Axis scaling
  # ---------------------------------------------------------------------------
  stopifnot(is.logical(volcano_auto_axis))
  stopifnot(is.numeric(volcano_axis_expand_mult))
  stopifnot(volcano_axis_expand_mult > 0)
}
