# =============================================================================
# 07_duplicates.R
# Duplicate handling
# =============================================================================

collapse_duplicate_names <- function(mat, feature_tbl,
                                     strategy = c("keep_separate","collapse_mean","collapse_sum","collapse_best_qc_rsd"),
                                     qc_rsd = NULL,
                                     audit_path = NULL) {
  strategy <- match.arg(strategy)
  
  ft <- feature_tbl %>%
    mutate(
      is_named = !is.na(Name_canon) & Name_canon != "",
      dup_key = if_else(is_named, Name_canon, featureID)
    )
  
  audit <- ft %>%
    group_by(dup_key) %>%
    mutate(n_in_group = n()) %>%
    ungroup() %>%
    mutate(will_collapse = (is_named & n_in_group > 1))
  
  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)
  
  if (strategy == "keep_separate") {
    ft2 <- ft %>% mutate(display_name = make.unique(display_name, sep = "__dup"))
    return(list(mat = mat, feature = ft2, audit = audit))
  }
  
  named_groups <- split(audit$featureID[audit$is_named], audit$dup_key[audit$is_named])
  unnamed_ids <- audit$featureID[!audit$is_named]
  
  keep_named_ids <- character(0)
  collapsed_vecs <- list()
  collapsed_feats <- list()
  
  if (strategy == "collapse_best_qc_rsd") {
    if (is.null(qc_rsd)) stop("collapse_best_qc_rsd requires qc_rsd (named by featureID).")
    
    for (k in names(named_groups)) {
      ids <- named_groups[[k]]
      
      if (length(ids) == 1) {
        keep_named_ids <- c(keep_named_ids, ids)
      } else {
        rsd <- unname(qc_rsd[ids])
        rsd_sort <- ifelse(is.finite(rsd), rsd, Inf)
        chosen <- ids[order(rsd_sort, ids)][1]
        keep_named_ids <- c(keep_named_ids, chosen)
      }
    }
    
    keep_ids <- c(keep_named_ids, unnamed_ids)
    mat2 <- mat[, keep_ids, drop = FALSE]
    
    ft2 <- feature_tbl %>%
      filter(featureID %in% keep_ids) %>%
      mutate(
        display_name = if_else(!is.na(Name_canon) & Name_canon != "", Name_canon, display_name),
        display_name = make.unique(display_name, sep = "__dup")
      )
    
    return(list(mat = mat2, feature = ft2, audit = audit))
  }
  
  for (k in names(named_groups)) {
    ids <- named_groups[[k]]
    
    if (length(ids) == 1) {
      keep_named_ids <- c(keep_named_ids, ids)
    } else {
      sub <- mat[, ids, drop = FALSE]
      agg <- if (strategy == "collapse_sum") rowSums(sub, na.rm = TRUE) else rowMeans(sub, na.rm = TRUE)
      
      new_id <- make.unique(paste0("COLLAPSED_", make.names(k)), sep = "_r")
      collapsed_vecs[[new_id]] <- agg
      
      base_row <- feature_tbl %>% filter(featureID %in% ids) %>% slice(1)
      base_row$featureID <- new_id
      base_row$Name <- k
      base_row$Name_canon <- k
      base_row$display_name <- k
      collapsed_feats[[new_id]] <- base_row
    }
  }
  
  mat_keep <- mat[, c(keep_named_ids, unnamed_ids), drop = FALSE]
  
  if (length(collapsed_vecs) > 0) {
    mat_coll <- do.call(cbind, collapsed_vecs)
    mat_out <- cbind(mat_keep, mat_coll)
    ft_keep <- feature_tbl %>% filter(featureID %in% colnames(mat_keep))
    ft_coll <- bind_rows(collapsed_feats)
    ft_out <- bind_rows(ft_keep, ft_coll)
  } else {
    mat_out <- mat_keep
    ft_out <- feature_tbl %>% filter(featureID %in% colnames(mat_out))
  }
  
  ft_out <- ft_out %>%
    mutate(
      display_name = if_else(!is.na(Name_canon) & Name_canon != "", Name_canon, display_name),
      display_name = make.unique(display_name, sep = "__dup")
    )
  
  list(mat = mat_out, feature = ft_out, audit = audit)
}