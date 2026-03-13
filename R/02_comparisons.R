# =============================================================================
# 02_comparisons.R
# Comparison definitions
# =============================================================================

COMPARISON_CONFIGS <- list(
  ALL_TGvsWT = list(
    meta_filter = function(m) dplyr::filter(m, group %in% c("TG", "WT")),
    stats_compare_var = "group",
    stats_den = "WT",
    stats_num = "TG",
    pca_color_var = "group",
    pca_shape_var = "sex",
    prefix = "ALL",
    label = "TG vs WT | sex=ALL"
  ),
  
  F_TGvsWT = list(
    meta_filter = function(m) dplyr::filter(m, sex == "F", group %in% c("TG", "WT")),
    stats_compare_var = "group",
    stats_den = "WT",
    stats_num = "TG",
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "sex",
    label = "TG vs WT | sex=F"
  ),
  
  M_TGvsWT = list(
    meta_filter = function(m) dplyr::filter(m, sex == "M", group %in% c("TG", "WT")),
    stats_compare_var = "group",
    stats_den = "WT",
    stats_num = "TG",
    pca_color_var = "group",
    pca_shape_var = NULL,
    prefix = "sex",
    label = "TG vs WT | sex=M"
  ),
  
  TG_FvsM = list(
    meta_filter = function(m) dplyr::filter(m, group == "TG", sex %in% c("F", "M")),
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "sex_FvsM",
    label = "F vs M within TG"
  ),
  
  WT_FvsM = list(
    meta_filter = function(m) dplyr::filter(m, group == "WT", sex %in% c("F", "M")),
    stats_compare_var = "sex",
    stats_den = "M",
    stats_num = "F",
    pca_color_var = "sex",
    pca_shape_var = NULL,
    prefix = "sex_FvsM",
    label = "F vs M within WT"
  )
)

COMPARISON_NAMES <- names(COMPARISON_CONFIGS)