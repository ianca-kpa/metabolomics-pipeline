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
# - reference_or_best_qc_rsd
# - keep_separate
# - collapse_mean
# - collapse_sum
# - collapse_best_qc_rsd
#
# An audit table is also generated to show how duplicates were handled.

collapse_duplicate_names <- function(mat, feature_tbl,
                                     strategy = c("reference_or_best_qc_rsd", "keep_separate", "collapse_mean", "collapse_sum", "collapse_best_qc_rsd"),
                                     reference_tbl = NULL,
                                     sanitize_mode = "greek_latin_ascii",
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
        strategy == "reference_or_best_qc_rsd" & is_named & n_in_group > 1 ~ "selected_reference_or_best_qc_rsd",
        strategy == "collapse_best_qc_rsd" & is_named & n_in_group > 1 ~ "selected_best_qc_rsd",
        strategy %in% c("collapse_mean", "collapse_sum") & is_named & n_in_group > 1 ~ "collapsed",
        TRUE ~ "unknown"
      )
    )

  # ---------------------------------------------------------------------------
  # Keep all features separate
  # ---------------------------------------------------------------------------
  if (strategy == "keep_separate") {
    mat2 <- mat[, feature_tbl$featureID, drop = FALSE]

    ft2 <- ft %>%
      dplyr::mutate(display_name = make.unique(display_name, sep = "__dup")) %>%
      dplyr::select(dplyr::all_of(orig_feature_cols))

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    return(list(mat = mat2, feature = ft2, audit = audit))
  }

  named_groups <- split(audit$featureID[audit$is_named], audit$dup_key[audit$is_named])
  unnamed_ids <- audit$featureID[!audit$is_named]

  keep_named_ids <- character(0)
  collapsed_vecs <- list()
  collapsed_feats <- list()
  used_ids <- feature_tbl$featureID

  # ---------------------------------------------------------------------------
  # Strategy: use comparison file or best QC RSD to choose representative
  # ---------------------------------------------------------------------------
  if (strategy == "reference_or_best_qc_rsd") {
    normalize_text <- function(x) {
      if (exists("sanitize_text_for_exports", mode = "function")) {
        x <- sanitize_text_for_exports(x, mode = sanitize_mode)
      }
      x <- as.character(x)
      x <- trimws(tolower(x))
      x[x %in% c("", "na", "nan", "null")] <- NA_character_
      x
    }

    parse_num_local <- function(x) {
      x <- as.character(x)
      x <- gsub(",", ".", x, fixed = TRUE)
      suppressWarnings(as.numeric(x))
    }

    find_col_by_pattern <- function(df, patterns, fallback_first = FALSE) {
      nms <- names(df)
      nms_norm <- tolower(gsub("[^a-z0-9]+", "", nms))
      for (p in patterns) {
        idx <- which(grepl(p, nms_norm, perl = TRUE))
        if (length(idx) > 0) return(nms[idx[1]])
      }
      if (fallback_first && length(nms) > 0) return(nms[1])
      NA_character_
    }

    choose_best_qc <- function(ids, qc_rsd_vec) {
      if (!is.null(qc_rsd_vec) && !is.null(names(qc_rsd_vec))) {
        rsd <- unname(qc_rsd_vec[ids])
        rsd_sort <- ifelse(is.finite(rsd), rsd, Inf)
        return(ids[order(rsd_sort, ids)][1])
      }
      ids[order(ids)][1]
    }

    mz_digits_local <- if (exists("dup_mz_digits", inherits = TRUE)) {
      as.integer(get("dup_mz_digits", inherits = TRUE))
    } else {
      4L
    }

    rt_digits_local <- if (exists("dup_rt_digits", inherits = TRUE)) {
      as.integer(get("dup_rt_digits", inherits = TRUE))
    } else {
      2L
    }

    if (is.null(reference_tbl)) {
      comparison_path_local <- if (exists("comparison_path", inherits = TRUE)) {
        get("comparison_path", inherits = TRUE)
      } else {
        NULL
      }

      comparison_sheet_local <- if (exists("comparison_sheet", inherits = TRUE)) {
        get("comparison_sheet", inherits = TRUE)
      } else {
        1
      }

      if (!is.null(comparison_path_local) &&
          nzchar(as.character(comparison_path_local)) &&
          exists("read_any_table", mode = "function")) {
        reference_tbl <- tryCatch(
          read_any_table(comparison_path_local, comparison_sheet_local),
          error = function(e) NULL
        )
      }
    }

    if (is.null(reference_tbl) || nrow(reference_tbl) == 0) {
      warning(
        "reference_or_best_qc_rsd: no reference table found. Falling back to best QC RSD for all duplicate groups.",
        call. = FALSE
      )
    }

    if (!is.null(qc_rsd) && is.null(names(qc_rsd))) {
      stop("reference_or_best_qc_rsd requires qc_rsd to be a named vector by featureID when provided.")
    }

    ft_local <- feature_tbl
    ft_local$.name_norm <- normalize_text(ft_local$Name_canon)
    ft_local$.mz_round <- round(ft_local$mz, digits = mz_digits_local)
    ft_local$.rt_round <- round(ft_local$RT, digits = rt_digits_local)

    ref_proc <- NULL
    if (!is.null(reference_tbl) && nrow(reference_tbl) > 0) {
      metab_col <- find_col_by_pattern(reference_tbl, c("metabolite", "compound", "name"), fallback_first = TRUE)
      refion_col <- find_col_by_pattern(reference_tbl, c("referenceion", "refion", "^refion$", "ion"))
      mz_col <- find_col_by_pattern(reference_tbl, c("mz", "masstocharge", "moverz", "masscharge"))
      rt_col <- find_col_by_pattern(reference_tbl, c("^rt$", "^rtmin$", "retentiontime", "retention"))

      if (is.na(metab_col) || is.na(refion_col) || is.na(mz_col) || is.na(rt_col)) {
        warning(
          "reference_or_best_qc_rsd: could not detect required columns in reference table (metabolite, Ref ion, m/z, RT). Falling back to best QC RSD.",
          call. = FALSE
        )
      } else {
        ref_proc <- reference_tbl %>%
          dplyr::mutate(
            .name_norm = normalize_text(.data[[metab_col]]),
            .refion_norm = normalize_text(.data[[refion_col]]),
            .mz_round = round(parse_num_local(.data[[mz_col]]), digits = mz_digits_local),
            .rt_round = round(parse_num_local(.data[[rt_col]]), digits = rt_digits_local)
          ) %>%
          dplyr::filter(!is.na(.name_norm), !is.na(.refion_norm), is.finite(.mz_round), is.finite(.rt_round))
      }
    }

    if ("Ref ion" %in% names(ft_local)) {
      ft_local$.refion_norm <- normalize_text(ft_local[["Ref ion"]])
    } else if ("Reference ion" %in% names(ft_local)) {
      ft_local$.refion_norm <- normalize_text(ft_local[["Reference ion"]])
    } else if ("ref_ion" %in% names(ft_local)) {
      ft_local$.refion_norm <- normalize_text(ft_local[["ref_ion"]])
    } else {
      ft_local$.refion_norm <- NA_character_
    }

    group_keys <- names(named_groups)
    ref_option_count <- setNames(integer(length(group_keys)), group_keys)
    full_match_option_count <- setNames(integer(length(group_keys)), group_keys)
    selection_source <- setNames(rep("", length(group_keys)), group_keys)

    for (k in names(named_groups)) {
      ids <- named_groups[[k]]

      refs_k <- NULL
      if (!is.null(ref_proc) && nrow(ref_proc) > 0) {
        refs_k <- ref_proc %>% dplyr::filter(.name_norm == normalize_text(k))
      }

      n_ref <- if (is.null(refs_k)) 0L else nrow(refs_k)
      ref_option_count[k] <- n_ref

      candidates <- ft_local %>% dplyr::filter(featureID %in% ids)
      matched_ids <- character(0)

      if (n_ref > 0) {
        has_match <- vapply(seq_len(nrow(candidates)), function(i) {
          rr <- candidates[i, , drop = FALSE]
          if (is.na(rr$.refion_norm) || !is.finite(rr$.mz_round) || !is.finite(rr$.rt_round)) {
            return(FALSE)
          }
          any(
            refs_k$.refion_norm == rr$.refion_norm &
              refs_k$.mz_round == rr$.mz_round &
              refs_k$.rt_round == rr$.rt_round,
            na.rm = TRUE
          )
        }, logical(1))
        matched_ids <- candidates$featureID[has_match]
      }

      full_match_option_count[k] <- length(matched_ids)

      chosen <- NA_character_
      if (length(ids) == 1) {
        chosen <- ids[1]
        selection_source[k] <- "kept_single_named"
      } else if (length(matched_ids) == 1) {
        chosen <- matched_ids[1]
        selection_source[k] <- "reference_full_match"
      } else if (length(matched_ids) > 1) {
        chosen <- choose_best_qc(matched_ids, qc_rsd)
        selection_source[k] <- "reference_tie_best_qc_rsd"
      } else {
        chosen <- choose_best_qc(ids, qc_rsd)
        selection_source[k] <- "best_qc_rsd_fallback"
      }

      keep_named_ids <- c(keep_named_ids, chosen)
    }

    group_stats <- tibble::tibble(
      dup_key = names(ref_option_count),
      reference_has_metabolite = as.logical(ref_option_count > 0),
      reference_option_count = as.integer(ref_option_count),
      full_match_option_count = as.integer(full_match_option_count),
      selection_source = as.character(selection_source)
    )

    audit <- audit %>%
      dplyr::left_join(group_stats, by = "dup_key")

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

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

    return(list(mat = mat2, feature = ft2, audit = audit))
  }
  
  
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

    if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

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

  if (!is.null(audit_path)) write_csv_safe(audit, audit_path)

  list(mat = mat_out, feature = ft_out, audit = audit)
}