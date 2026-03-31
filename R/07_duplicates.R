# =============================================================================
# 07_duplicates.R
# Duplicate handling
# =============================================================================

# Functions for handling duplicate features in the feature table, based on the
# settings specified in config/settings.R.
#
# The main function is collapse_duplicate_names(), which takes the assay matrix
# and feature table, identifies duplicates based on the "Name_canon" column,
# and handles them according to the specified strategy:
# - keep_separate
# - collapse_mean
# - collapse_sum
# - collapse_best_qc_rsd
#
# An audit table is also generated to show how duplicates were handled.

collapse_duplicate_names <- function(mat, feature_tbl,
                                     strategy = c("keep_separate", "collapse_mean", "collapse_sum", "collapse_best_qc_rsd"),
                                     qc_rsd = NULL,
                                     audit_path = NULL) {
  strategy <- match.arg(strategy)
  orig_feature_cols <- colnames(feature_tbl)

  # ---------------------------------------------------------------------------
  # Basic validation
  # ---------------------------------------------------------------------------
  required_cols <- c("featureID", "Name_canon", "display_name")
  missing_cols <- setdiff(required_cols, names(feature_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "collapse_duplicate_names(): feature_tbl is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (is.null(colnames(mat))) {
    stop("collapse_duplicate_names(): mat must have column names matching featureID.")
  }

  if (anyDuplicated(feature_tbl$featureID) > 0) {
    dup_ids <- unique(feature_tbl$featureID[duplicated(feature_tbl$featureID)])
    stop(
      "collapse_duplicate_names(): feature_tbl has duplicated featureID values: ",
      paste(utils::head(dup_ids, 10), collapse = ", "),
      if (length(dup_ids) > 10) " ..." else ""
    )
  }

  if (anyDuplicated(colnames(mat)) > 0) {
    dup_cols <- unique(colnames(mat)[duplicated(colnames(mat))])
    stop(
      "collapse_duplicate_names(): mat has duplicated column names: ",
      paste(utils::head(dup_cols, 10), collapse = ", "),
      if (length(dup_cols) > 10) " ..." else ""
    )
  }

  missing_in_mat <- setdiff(feature_tbl$featureID, colnames(mat))
  if (length(missing_in_mat) > 0) {
    stop(
      "collapse_duplicate_names(): some featureIDs from feature_tbl are missing in mat: ",
      paste(utils::head(missing_in_mat, 10), collapse = ", "),
      if (length(missing_in_mat) > 10) " ..." else ""
    )
  }

  extra_in_mat <- setdiff(colnames(mat), feature_tbl$featureID)
  if (length(extra_in_mat) > 0) {
    stop(
      "collapse_duplicate_names(): mat has columns not present in feature_tbl$featureID: ",
      paste(utils::head(extra_in_mat, 10), collapse = ", "),
      if (length(extra_in_mat) > 10) " ..." else ""
    )
  }

  make_unique_id <- function(base_id, used_ids) {
    candidate <- base_id
    i <- 1
    while (candidate %in% used_ids) {
      candidate <- paste0(base_id, "_r", i)
      i <- i + 1
    }
    candidate
  }

  ft <- feature_tbl %>%
    dplyr::mutate(
      is_named = !is.na(Name_canon) & Name_canon != "",
      dup_key = dplyr::if_else(is_named, Name_canon, featureID)
    )

  audit <- ft %>%
    dplyr::group_by(dup_key) %>%
    dplyr::mutate(n_in_group = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      will_collapse = (is_named & n_in_group > 1),
      duplicate_action = dplyr::case_when(
        strategy == "keep_separate" ~ "kept_separate",
        !is_named ~ "kept_unnamed",
        is_named & n_in_group == 1 ~ "kept_single_named",
        strategy == "collapse_best_qc_rsd" & is_named & n_in_group > 1 ~ "selected_best_qc_rsd",
        strategy %in% c("collapse_mean", "collapse_sum") & is_named & n_in_group > 1 ~ "collapsed",
        TRUE ~ "unknown"
      )
    )

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  # ---------------------------------------------------------------------------
  # Keep all features separate
  # ---------------------------------------------------------------------------
  if (strategy == "keep_separate") {
    mat2 <- mat[, feature_tbl$featureID, drop = FALSE]

    ft2 <- ft %>%
      dplyr::mutate(display_name = make.unique(display_name, sep = "__dup")) %>%
      dplyr::select(dplyr::all_of(orig_feature_cols))

    return(list(mat = mat2, feature = ft2, audit = audit))
  }

  named_groups <- split(audit$featureID[audit$is_named], audit$dup_key[audit$is_named])
  unnamed_ids <- audit$featureID[!audit$is_named]

  keep_named_ids <- character(0)
  collapsed_vecs <- list()
  collapsed_feats <- list()
  used_ids <- feature_tbl$featureID

  # ---------------------------------------------------------------------------
  # Strategy: choose the best representative by QC RSD
  # ---------------------------------------------------------------------------
  if (strategy == "collapse_best_qc_rsd") {
    if (is.null(qc_rsd)) {
      stop("collapse_best_qc_rsd requires qc_rsd (named by featureID).")
    }

    if (is.null(names(qc_rsd))) {
      stop("collapse_best_qc_rsd requires qc_rsd to be a named vector by featureID.")
    }

    named_ids <- unique(audit$featureID[audit$is_named])
    missing_qc <- setdiff(named_ids, names(qc_rsd))

    if (length(named_ids) > 0 && length(missing_qc) == length(named_ids)) {
      stop("collapse_best_qc_rsd found no qc_rsd values for named features.")
    }

    if (length(missing_qc) > 0) {
      warning(
        paste0(
          "collapse_best_qc_rsd: ",
          length(missing_qc),
          " named featureIDs are missing qc_rsd values. ",
          "Missing/non-finite values will be treated as Inf, ",
          "and ties will be resolved by featureID order."
        ),
        call. = FALSE
      )
    }

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
      dplyr::filter(featureID %in% keep_ids) %>%
      dplyr::slice(match(keep_ids, featureID)) %>%
      dplyr::mutate(
        display_name = dplyr::if_else(
          !is.na(Name_canon) & Name_canon != "",
          Name_canon,
          display_name
        ),
        display_name = make.unique(display_name, sep = "__dup")
      ) %>%
      dplyr::select(dplyr::all_of(orig_feature_cols))

    return(list(mat = mat2, feature = ft2, audit = audit))
  }

  # ---------------------------------------------------------------------------
  # Strategies: collapse_mean / collapse_sum
  # ---------------------------------------------------------------------------
  for (k in names(named_groups)) {
    ids <- named_groups[[k]]

    if (length(ids) == 1) {
      keep_named_ids <- c(keep_named_ids, ids)
    } else {
      sub <- mat[, ids, drop = FALSE]
      agg <- if (strategy == "collapse_sum") {
        rowSums(sub, na.rm = TRUE)
      } else {
        rowMeans(sub, na.rm = TRUE)
      }

      # Preserve NA when all duplicate values are NA for a sample row
      all_na_rows <- apply(sub, 1, function(x) all(is.na(x)))
      agg[all_na_rows] <- NA

      base_id <- paste0("COLLAPSED_", make.names(k))
      new_id <- make_unique_id(base_id, used_ids)
      used_ids <- c(used_ids, new_id)
      collapsed_vecs[[new_id]] <- agg

      base_row <- feature_tbl %>%
        dplyr::filter(featureID %in% ids) %>%
        dplyr::slice(1)

      base_row$featureID <- new_id
      if ("Name" %in% names(base_row)) {
        base_row$Name <- k
      }
      base_row$Name_canon <- k
      base_row$display_name <- k

      collapsed_feats[[new_id]] <- base_row
    }
  }

  mat_keep <- mat[, c(keep_named_ids, unnamed_ids), drop = FALSE]

  if (length(collapsed_vecs) > 0) {
    mat_coll <- do.call(cbind, collapsed_vecs)
    mat_out <- cbind(mat_keep, mat_coll)

    ft_keep <- feature_tbl %>%
      dplyr::filter(featureID %in% colnames(mat_keep)) %>%
      dplyr::slice(match(colnames(mat_keep), featureID))

    ft_coll <- dplyr::bind_rows(collapsed_feats)
    ft_out <- dplyr::bind_rows(ft_keep, ft_coll)
  } else {
    mat_out <- mat_keep
    ft_out <- feature_tbl %>%
      dplyr::filter(featureID %in% colnames(mat_out)) %>%
      dplyr::slice(match(colnames(mat_out), featureID))
  }

  ft_out <- ft_out %>%
    dplyr::filter(featureID %in% colnames(mat_out)) %>%
    dplyr::slice(match(colnames(mat_out), featureID)) %>%
    dplyr::mutate(
      display_name = dplyr::if_else(
        !is.na(Name_canon) & Name_canon != "",
        Name_canon,
        display_name
      ),
      display_name = make.unique(display_name, sep = "__dup")
    ) %>%
    dplyr::select(dplyr::all_of(orig_feature_cols))

  list(mat = mat_out, feature = ft_out, audit = audit)
}