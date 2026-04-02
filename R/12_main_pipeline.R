# =============================================================================
# 12_main_pipeline.R
# Main pipeline runner
# =============================================================================

run_untargeted_pipeline <- function() {
  pipeline_t0 <- Sys.time()

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(output_dir, "PIPELINE_LOG.txt")

  append_log_line(log_path, "============================================================")
  append_log_line(log_path, "PIPELINE LOG — Untargeted Metabolomics - Data Processing")
  append_log_line(log_path, paste0("Date/time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  append_log_line(log_path, "============================================================")
  append_log_line(log_path, "")
  append_log_line(log_path, "[INPUTS]")
  append_log_line(log_path, paste0("cd_file_path: ", cd_file_path))
  append_log_line(log_path, paste0("cd_sheet: ", cd_sheet))
  append_log_line(log_path, paste0("metadata_path: ", metadata_path))
  append_log_line(log_path, paste0("metadata_sheet: ", metadata_sheet))
  append_log_line(log_path, paste0("output_dir: ", output_dir))
  append_log_line(log_path, paste0("use_only_known: ", use_only_known))
  append_log_line(log_path, paste0("duplicate_name_strategy: ", duplicate_name_strategy))
  append_log_line(log_path, paste0("active_variant: ", active_variant))
  append_log_line(log_path, "")
  append_log_line(log_path, "[CONSOLE OUTPUT]")

  with_console_capture_to_file(log_path, {
    # -------------------------------------------------------------------------
    # STEP 02 — Reading input files
    # -------------------------------------------------------------------------
    t_step <- step_start(2, 15, "Reading input files")

    # message("Step 3: Reading input files...")
    cd_raw <- read_any_table(cd_file_path, cd_sheet)
    metadata <- read_any_table(metadata_path, metadata_sheet)
    reference_tbl <- NULL

    if (exists("comparison_path", inherits = TRUE) &&
        !is.null(comparison_path) &&
        nzchar(as.character(comparison_path))) {
      comparison_sheet_local <- if (exists("comparison_sheet", inherits = TRUE)) comparison_sheet else 1
      reference_tbl <- read_any_table(comparison_path, comparison_sheet_local)
      step_info("Reference table: rows=", nrow(reference_tbl), " cols=", ncol(reference_tbl))
    } else {
      step_info("Reference table: not provided (comparison_path empty)")
    }

    step_info("Compound Discoverer table: rows=", nrow(cd_raw), " cols=", ncol(cd_raw))
    step_info("Metadata table: rows=", nrow(metadata), " cols=", ncol(metadata))
    step_ok("Reading input files", t_step)

    # -------------------------------------------------------------------------
    # STEP 03 — Standardizing metadata
    # -------------------------------------------------------------------------
    t_step <- step_start(3, 15, "Standardizing metadata")

    metadata <- clean_metadata(metadata)
    metadata_aligned <- metadata

    step_info("Metadata after cleaning: rows=", nrow(metadata), " cols=", ncol(metadata))
    step_ok("Standardizing metadata", t_step)

    # -------------------------------------------------------------------------
    # STEP 04 — Creating per-model output directory structure
    # -------------------------------------------------------------------------
    t_step <- step_start(4, 15, "Creating per-model output directory structure")

    models_detected <- get_models_from_metadata(metadata_aligned)
    step_info("Models detected: ", paste(models_detected, collapse = ", "))

    paths <- setup_output_dirs(output_dir, model_names = models_detected)
    filter_summary <- tibble()
    summary_csv <- file.path(paths$global$audits, "filter_summary.csv")

    step_ok("Creating per-model output directory structure", t_step)

    # -------------------------------------------------------------------------
    # STEP 05 — Building feature table
    # -------------------------------------------------------------------------
    t_step <- step_start(5, 15, "Building feature table")

    feature_tbl <- build_feature_table(
      cd_raw,
      sanitize_names_for_exports = sanitize_names_for_exports,
      sanitize_mode = sanitize_mode,
      mz_digits = dup_mz_digits,
      rt_digits = dup_rt_digits
    )

    out_feat_map <- file.path(paths$global$exports, "02_featureID_to_display_name_map.csv")
    write_csv_safe(
      feature_tbl %>% select(featureID, display_name, Name, Name_canon, dplyr::any_of("Ref ion"), Formula, mz, RT),
      out_feat_map
    )
    log_written_object(log_path, out_feat_map, feature_tbl, note = "feature table")

    step_info("Feature table created: rows=", nrow(feature_tbl), " cols=", ncol(feature_tbl))
    step_ok("Building feature table", t_step)

    # -------------------------------------------------------------------------
    # STEP 06 — Building raw assay
    # -------------------------------------------------------------------------
    t_step <- step_start(6, 15, "Building raw assay")

    assay_bundle <- build_assay_from_cd(cd_raw, feature_tbl, metadata_aligned, paths)

    step_info("DEBUG — assay_obj class: ", paste(class(assay_bundle), collapse = ", "))
    step_info("DEBUG — assay_obj names: ", paste(names(assay_bundle), collapse = ", "))
    assay_num_raw <- assay_bundle$assay_num_raw
    metadata_aligned <- assay_bundle$metadata_aligned
    qc_idx <- assay_bundle$qc_idx
    sample_idx <- assay_bundle$sample_idx

    step_info("Raw assay created: rows=", nrow(assay_num_raw), " cols=", ncol(assay_num_raw))
    step_info("Sample index length: ", length(sample_idx))
    step_info("QC index length: ", length(qc_idx))
    step_ok("Building raw assay", t_step)

    # -------------------------------------------------------------------------
    # STEP 07 — Applying weight normalization
    # -------------------------------------------------------------------------
    t_step <- step_start(7, 15, "Applying weight normalization")

    assay_num_weight <- normalize_by_weight(
      assay_num_raw,
      metadata_aligned,
      sample_idx,
      use_weight_normalization = use_weight_normalization,
      stop_on_invalid_weight   = stop_on_invalid_weight,
      invalid_weight_to_NA     = invalid_weight_to_NA
    )


    step_info("Weight normalization applied. NAs introduced for invalid weights: ", sum(is.na(assay_num_weight)))
    step_ok("Applying weight normalization", t_step)

    # -------------------------------------------------------------------------
    # STEP 08 — Applying PQN normalization
    # -------------------------------------------------------------------------
    t_step <- step_start(8, 15, "Applying PQN normalization")

    pqn_bundle <- normalize_pqn_qc_ref(assay_num_weight, qc_idx)
    assay_num_pqn <- pqn_bundle$assay_num_pqn

    if (!minimal_outputs) {
      out_pqn <- file.path(paths$global$exports, "06_pqn_factors_weight_then_PQN.csv")
      df_pqn <- bind_cols(metadata_aligned %>% select(sample, type), pqn_bundle$pqn_tbl)
      write_csv_safe(df_pqn, out_pqn)
      log_written_object(log_path, out_pqn, df_pqn, note = "PQN factors table")
    }

    step_info("PQN normalization applied.")
    step_ok("Applying PQN normalization", t_step)

    # -------------------------------------------------------------------------
    # Initialize filter summary table
    # -------------------------------------------------------------------------
    filter_summary <- tibble::tibble(
      step = character(),
      n_features_before = integer(),
      n_features_after = integer(),
      n_removed = integer(),
      pct_removed = numeric()
    )

    # -------------------------------------------------------------------------
    # STEP 09 — Applying missing-value exclusion
    # -------------------------------------------------------------------------
    t_step <- step_start(9, 15, "Applying missing-value exclusion")

    n0 <- ncol(assay_num_pqn)

    miss_bundle <- filter_missing_exclusion(
      assay_num_pqn,
      feature_tbl,
      sample_idx,
      max_missing_fraction = missing_exclusion_max_fraction,
      audit_path = file.path(paths$global$audits, "missing_exclusion_audit.csv")
    )

    assay_work <- miss_bundle$assay
    feature_tbl_work <- miss_bundle$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "missing_exclusion",
      n0,
      ncol(assay_work),
      out_csv = summary_csv
    )

    step_info("Features before filter: ", n0)
    step_info("Features after filter: ", ncol(assay_work))
    step_ok("Applying missing-value exclusion", t_step)

    # -------------------------------------------------------------------------
    # STEP 10 — Applying presence filter and imputation
    # -------------------------------------------------------------------------
    t_step <- step_start(10, 15, "Applying presence filter and imputation")

    n1 <- ncol(assay_work)

    miss2 <- presence_filter_and_impute(
      assay_work,
      feature_tbl_work,
      sample_idx,
      min_fraction = presence_filter_min_fraction,
      impute_half_min = impute_half_min,
      audit_path = file.path(paths$global$audits, "presence_filter_audit.csv")
    )

    assay_work <- miss2$assay
    feature_tbl_work <- miss2$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "presence_filter",
      n1,
      ncol(assay_work),
      out_csv = summary_csv
    )

    step_info("Presence filter min fraction: ", presence_filter_min_fraction)
    step_info("Features before filter: ", n1)
    step_info("Features after filter: ", ncol(assay_work))
    step_ok("Applying presence filter and imputation", t_step)

    # -------------------------------------------------------------------------
    # STEP 11 — Applying feature-level filters and duplicate handling
    # -------------------------------------------------------------------------
    t_step <- step_start(11, 15, "Applying feature-level filters and duplicate handling")

    n2 <- ncol(assay_work)

    known <- filter_known(
      assay_work,
      feature_tbl_work,
      use_only_known,
      audit_path = file.path(paths$global$audits, "known_filter_audit.csv")
    )

    assay_work <- known$assay
    feature_tbl_work <- known$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "known_filter",
      n2,
      ncol(assay_work),
      out_csv = summary_csv
    )

    qc_rsd_all <- calc_qc_rsd(assay_work, qc_idx)

    df_qc_rsd_pre <- tibble(
      featureID = names(qc_rsd_all),
      qc_rsd = as.numeric(qc_rsd_all)
    ) %>% arrange(qc_rsd)


    out_qc_rsd_pre <- file.path(paths$global$audits, "qc_rsd_values_pre_variants.csv")
    write_csv_safe(df_qc_rsd_pre, out_qc_rsd_pre)
    log_written_object(log_path, out_qc_rsd_pre, df_qc_rsd_pre, note = "QC RSD values (pre-variants)")

    variants <- list()
    variants$BASE <- list(mat = assay_work, feature = feature_tbl_work)

    for (thr in rsd_thresholds) {
      keep <- names(qc_rsd_all)[!is.na(qc_rsd_all) & qc_rsd_all <= thr]
      variants[[paste0("QC_RSD", thr)]] <- list(
        mat = assay_work[, keep, drop = FALSE],
        feature = feature_tbl_work %>% filter(featureID %in% keep)
      )
    }

    step_info("Variants created based on RSD thresholds: ", paste(names(variants), collapse = ", "))

    if (!active_variant %in% names(variants)) {
      stop(
        "active_variant not found: ", active_variant,
        "\nAvailable: ", paste(names(variants), collapse = ", ")
      )
    }

    variants$ACTIVE <- variants[[active_variant]]
    step_info("Active variant selected: ", active_variant, " | n_features=", ncol(variants$ACTIVE$mat))

    if (!minimal_outputs) {
      out_post_rsd_mat <- file.path(
        paths$global$exports,
        paste0("10_MATRIX_post_", active_variant, "_postRSD_preLowVar_preDup_ALL.csv")
      )
      dir.create(dirname(out_post_rsd_mat), recursive = TRUE, showWarnings = FALSE)
      write.csv(variants$ACTIVE$mat, file = out_post_rsd_mat, row.names = TRUE)
      log_written_object(
        log_path,
        out_post_rsd_mat,
        variants$ACTIVE$mat,
        note = "post-RSD, pre-lowVar, pre-dup"
      )

      df_post_rsd_named <- as.data.frame(variants$ACTIVE$mat, check.names = FALSE) %>%
        tibble::rownames_to_column("sample")

      df_post_rsd_named <- rename_feature_cols(
        df_post_rsd_named,
        variants$ACTIVE$feature,
        sample_col = "sample"
      )

      out_post_rsd_named <- file.path(
        paths$global$exports,
        paste0("10_TABLE_post_", active_variant, "_postRSD_preLowVar_preDup_ALL_NAMED.csv")
      )
      write_csv_safe(df_post_rsd_named, out_post_rsd_named)
      log_written_object(
        log_path,
        out_post_rsd_named,
        df_post_rsd_named,
        note = "post-RSD named table"
      )
    }

    step_info("Applying feature-level filters and duplicate handling", t_step)

    n3 <- ncol(variants$ACTIVE$mat)

    lv <- filter_low_variance_deterministic(
      variants$ACTIVE$mat,
      variants$ACTIVE$feature,
      method = low_variance_filter_method,
      frac = low_variance_filter_fraction,
      sample_idx = sample_idx,
      audit_path = file.path(paths$global$audits, "low_variance_iqr_audit_ACTIVE.csv")
    )

    variants$ACTIVE$mat <- lv$assay
    variants$ACTIVE$feature <- lv$feature

    filter_summary <- append_filter_summary(
      filter_summary,
      "low_variance_iqr_ACTIVE",
      n3,
      ncol(variants$ACTIVE$mat),
      out_csv = summary_csv
    )

    step_info("Applying duplicate name handling on ACTIVE (collapse ONLY named by Name_canon)...")

    qc_rsd_active <- calc_qc_rsd(variants$ACTIVE$mat, qc_idx)

    dup_out <- collapse_duplicate_names(
      mat = variants$ACTIVE$mat,
      feature_tbl = variants$ACTIVE$feature,
      strategy = duplicate_name_strategy,
      reference_tbl = reference_tbl,
      sanitize_mode = sanitize_mode,
      qc_rsd = qc_rsd_active,
      audit_path = file.path(
        paths$global$audits,
        paste0("duplicate_name_audit_", duplicate_name_strategy, ".csv")
      )
    )

    variants$ACTIVE$mat <- dup_out$mat
    variants$ACTIVE$feature <- dup_out$feature

    if (duplicate_name_strategy == "reference_or_best_qc_rsd" &&
        all(c("reference_has_metabolite", "reference_option_count", "full_match_option_count", "selection_source") %in% names(dup_out$audit))) {
      dup_ref_summary <- dup_out$audit %>%
        dplyr::filter(is_named, n_in_group > 1) %>%
        dplyr::distinct(
          dup_key,
          reference_has_metabolite,
          reference_option_count,
          full_match_option_count,
          selection_source
        )

      out_dup_ref_summary <- file.path(
        paths$global$audits,
        "duplicate_name_reference_summary_reference_or_best_qc_rsd.csv"
      )
      write_csv_safe(dup_ref_summary, out_dup_ref_summary)
      log_written_object(
        log_path,
        out_dup_ref_summary,
        dup_ref_summary,
        note = "duplicate reference summary"
      )

      n_dup_groups <- nrow(dup_ref_summary)
      n_in_reference <- sum(dup_ref_summary$reference_has_metabolite, na.rm = TRUE)
      n_not_in_reference <- n_dup_groups - n_in_reference

      step_info(
        "Reference check (duplicate metabolites): groups=",
        n_dup_groups,
        " | in list=",
        n_in_reference,
        " | not in list=",
        n_not_in_reference
      )

      if (n_dup_groups > 0) {
        ref_opts <- dup_ref_summary$reference_option_count
        match_opts <- dup_ref_summary$full_match_option_count

        step_info(
          "Reference options per metabolite (min/median/max): ",
          paste0(min(ref_opts), "/", round(stats::median(ref_opts), 1), "/", max(ref_opts))
        )

        step_info(
          "Full-match options per metabolite (min/median/max): ",
          paste0(min(match_opts), "/", round(stats::median(match_opts), 1), "/", max(match_opts))
        )

        src_counts <- table(dup_ref_summary$selection_source, useNA = "ifany")
        step_info(
          "Selection source counts: ",
          paste0(names(src_counts), "=", as.integer(src_counts), collapse = ", ")
        )
      }
    }

    step_info("ACTIVE variant selected: ", active_variant)
    step_info("Final ACTIVE matrix: ", fmt_dims(variants$ACTIVE$mat))
    step_info("Duplicate strategy: ", duplicate_name_strategy)
    step_ok("Applying feature-level filters and duplicate handling", t_step)

    # -------------------------------------------------------------------------
    # STEP 12 — Applying log2 transform and MetaboAnalyst export
    # -------------------------------------------------------------------------
    t_step <- step_start(12, 15, "Applying log2 transform and MetaboAnalyst export")

    mat_prelog_base <- variants$ACTIVE$mat
    mat_log2_base <- log2_transform(mat_prelog_base, log2_offset)
    rownames(mat_log2_base) <- rownames(mat_prelog_base)

    base_log2_df <- as.data.frame(mat_log2_base) %>%
      tibble::rownames_to_column("sample")

    if (isTRUE(export_metaboanalyst_ready)) {
      models_for_export <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))

      for (m in models_for_export) {
        mp <- get_model_paths(paths, m)

        export_metaboanalyst_one_model(
          log2_df = base_log2_df,
          metadata_aligned = metadata_aligned,
          feature_tbl = variants$ACTIVE$feature,
          model_name = m,
          export_dir = mp$exports$metaboanalyst,
          log_path = log_path
        )
      }
    }

    step_info("MetaboAnalyst export enabled: ", export_metaboanalyst_ready)
    step_ok("Applying log2 transform and MetaboAnalyst export", t_step)

    # -------------------------------------------------------------------------
    # STEP 13 — Generating PCA plots
    # -------------------------------------------------------------------------
    t_step <- step_start(13, 15, "Generating PCA plots")

    suppressWarnings({
      plot_pca_per_model(
        mat_log2_base,
        metadata_aligned,
        paths,
        pca_scaling = pca_scaling
      )

      step_info("PCA scaling: ", pca_scaling)
      step_ok("Generating PCA plots", t_step)

      # -------------------------------------------------------------------------
      # STEP 14 — Preparing statistics-driven heatmaps
      # -------------------------------------------------------------------------
      t_step <- step_start(14, 15, "Preparing statistics-driven heatmaps")

      heatmap_rank_metrics_cfg <- if (exists("heatmap_rank_metrics", inherits = TRUE)) {
        heatmap_rank_metrics
      } else {
        run_metrics
      }

      step_info("Top/significant heatmaps will be generated after statistics.")
      step_info("Top heatmap metrics configured: ", paste(heatmap_rank_metrics_cfg, collapse = ", "))
      step_info("Significant heatmap metrics configured: ", paste(run_metrics, collapse = ", "))
      step_ok("Preparing statistics-driven heatmaps", t_step)

      # -------------------------------------------------------------------------
      # STEP 15 — Generating statistics, volcano plots, stats Excel, and heatmaps
      # -------------------------------------------------------------------------
      t_step <- step_start(15, 15, "Generating statistics, volcano plots, stats Excel, and heatmaps")
      metric_mode_resolver <- if (exists("normalize_metric_modes", mode = "function")) {
        normalize_metric_modes
      } else {
        function(metrics) {
          valid <- c("FDR", "p_value")
          if (is.null(metrics) || length(metrics) == 0) return(valid)
          expanded <- unlist(lapply(as.character(metrics), function(m) {
            if (identical(m, "FDR_and_p_value")) c("FDR", "p_value") else m
          }), use.names = FALSE)
          expanded <- unique(expanded)
          expanded <- expanded[expanded %in% valid]
          if (length(expanded) == 0) valid else expanded
        }
      }
      run_metrics_expanded <- metric_mode_resolver(run_metrics)
      heatmap_rank_metrics_expanded <- metric_mode_resolver(heatmap_rank_metrics_cfg)

      # -----------------------------------------------------------------------
      # Top heatmaps by model
      # -----------------------------------------------------------------------
      if (isTRUE(make_heatmap_by_model)) {
        step_info("TOP heatmaps per model (loop rank metrics)...")

        suppressWarnings({
          for (rk in heatmap_rank_metrics_expanded) {
            plot_heatmap_top_ttest_per_model(
              mat_log2_base,
              metadata_aligned,
              variants$ACTIVE$feature,
              paths,
              top_n = heatmap_top_n,
              rank_by = rk,
              split_by_sex = FALSE,
              order_samples_by_group = heatmap_order_samples_by_group,
              scale_method = heatmap_scale_method
            )
          }
        })
      }

      # -----------------------------------------------------------------------
      # Top heatmaps by model and sex
      # -----------------------------------------------------------------------
      if (isTRUE(make_heatmap_by_model_sex)) {
        step_info("TOP heatmaps per model split by sex (loop rank metrics)...")

        suppressWarnings({
          for (rk in heatmap_rank_metrics_expanded) {
            plot_heatmap_top_ttest_per_model(
              mat_log2_base,
              metadata_aligned,
              variants$ACTIVE$feature,
              paths,
              top_n = heatmap_top_n,
              rank_by = rk,
              split_by_sex = TRUE,
              order_samples_by_group = heatmap_order_samples_by_group,
              scale_method = heatmap_scale_method
            )
          }
        })
      }

      step_info("Stats (5 comparisons) + Volcano plots for BOTH metrics...")

      stats_5sets_by_model <- run_all_stats_5sets_per_model(
        mat_log2 = mat_log2_base,
        mat_prelog = mat_prelog_base,
        metadata_aligned = metadata_aligned,
        feat_info = variants$ACTIVE$feature,
        paths = paths,
        alpha_sig = alpha_sig,
        fc_cutoff_log2 = fc_cutoff_log2,
        run_metrics = run_metrics,
        make_volcano_plots = make_volcano_plots,
        comparison_configs = COMPARISON_CONFIGS
      )

      if (isTRUE(save_stats_excel_per_model)) {
        step_info("Exporting stats to Excel (README + 5 tabs per model)...")

        export_stats_excel_by_model(
          stats_5sets_by_model,
          paths = paths,
          alpha_sig = alpha_sig,
          fc_cutoff_log2 = fc_cutoff_log2,
          active_variant = active_variant,
          log_path = log_path
        )
      }

      models <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))
      sexes <- c("F", "M")

      # -----------------------------------------------------------------------
      # Significant heatmaps by model (ALL sex TG vs WT)
      # -----------------------------------------------------------------------
      if (isTRUE(make_sig_heatmap_by_model)) {
        step_info("Significant heatmaps (ALL sex TG/WT) for BOTH metrics...")

        for (met in run_metrics_expanded) {
          for (m in models) {
            mp <- get_model_paths(paths, m)

            meta_m <- metadata_aligned %>%
              dplyr::filter(type == "Sample", model == m, group %in% c("TG", "WT"))

            st <- stats_5sets_by_model[[m]]$ALL_TGvsWT
            if (is.null(st)) next

            out_png <- file.path(
              mp$plots$heatmap_SIG_ALL,
              paste0(
                "HEATMAP_SIG_ACTIVE_ALL_TGvsWT_", met, "_lt_", alpha_sig,
                "_model_", m, "_scale_", heatmap_scale_method, ".png"
              )
            )

            plot_sig_heatmap_from_stats(
              mat_log2_base,
              meta_m,
              variants$ACTIVE$feature,
              st,
              sig_metric = met,
              alpha_sig = alpha_sig,
              fc_cutoff_log2 = fc_cutoff_log2,
              require_fc_cutoff = sig_heatmap_require_fc_cutoff,
              sig_max = sig_heatmap_max_features,
              scale_method = heatmap_scale_method,
              order_samples_by_group = heatmap_order_samples_by_group,
              out_png = out_png,
              title_main = paste0(
                "SIG (TG/WT) | model=", m,
                " | sex=ALL | ", met, "<", alpha_sig,
                if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                " | scale=", heatmap_scale_method
              )
            )
          }
        }
      }

      # -----------------------------------------------------------------------
      # Significant heatmaps by model and sex (TG vs WT within F and within M)
      # -----------------------------------------------------------------------
      if (isTRUE(make_sig_heatmap_by_model_sex)) {
        step_info("Significant heatmaps (TG/WT) BY sex for BOTH metrics...")

        for (met in run_metrics_expanded) {
          for (m in models) {
            mp <- get_model_paths(paths, m)

            for (sx in sexes) {
              meta_m <- metadata_aligned %>%
                dplyr::filter(type == "Sample", model == m, sex == sx, group %in% c("TG", "WT"))

              st <- if (sx == "F") {
                stats_5sets_by_model[[m]]$F_TGvsWT
              } else {
                stats_5sets_by_model[[m]]$M_TGvsWT
              }

              if (is.null(st)) next

              out_png <- file.path(
                mp$plots$heatmap_SIG_sex,
                paste0(
                  "HEATMAP_SIG_ACTIVE_", sx, "_TGvsWT_", met, "_lt_", alpha_sig,
                  "_model_", m, "_scale_", heatmap_scale_method, ".png"
                )
              )

              plot_sig_heatmap_from_stats(
                mat_log2_base,
                meta_m,
                variants$ACTIVE$feature,
                st,
                sig_metric = met,
                alpha_sig = alpha_sig,
                fc_cutoff_log2 = fc_cutoff_log2,
                require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                sig_max = sig_heatmap_max_features,
                scale_method = heatmap_scale_method,
                order_samples_by_group = heatmap_order_samples_by_group,
                out_png = out_png,
                title_main = paste0(
                  "SIG (TG/WT) | model=", m,
                  " | sex=", sx,
                  " | ", met, "<", alpha_sig,
                  if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                  " | scale=", heatmap_scale_method
                )
              )
            }
          }
        }
      }

      # -----------------------------------------------------------------------
      # Significant heatmaps by model and sex (F vs M within TG and within WT)
      # -----------------------------------------------------------------------
      if (isTRUE(make_sig_heatmap_FvsM_within_group)) {
        step_info("Significant heatmaps (F vs M within TG and WT) for BOTH metrics...")

        for (met in run_metrics_expanded) {
          for (m in models) {
            mp <- get_model_paths(paths, m)

            meta_tg <- metadata_aligned %>%
              dplyr::filter(type == "Sample", model == m, group == "TG", sex %in% c("F", "M"))
            st_tg <- stats_5sets_by_model[[m]]$TG_FvsM

            if (!is.null(st_tg)) {
              out_png_tg <- file.path(
                mp$plots$heatmap_SIG_sex_FvsM,
                paste0(
                  "HEATMAP_SIG_ACTIVE_TG_FvsM_", met, "_lt_", alpha_sig,
                  "_model_", m, "_scale_", heatmap_scale_method, ".png"
                )
              )

              plot_sig_heatmap_from_stats(
                mat_log2_base,
                meta_tg,
                variants$ACTIVE$feature,
                st_tg,
                sig_metric = met,
                alpha_sig = alpha_sig,
                fc_cutoff_log2 = fc_cutoff_log2,
                require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                sig_max = sig_heatmap_max_features,
                scale_method = heatmap_scale_method,
                order_samples_by_group = FALSE,
                out_png = out_png_tg,
                title_main = paste0(
                  "SIG (F/M within TG) | model=", m,
                  " | ", met, "<", alpha_sig,
                  if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                  " | log2FC=log2(F/M) | scale=", heatmap_scale_method
                )
              )
            }

            meta_wt <- metadata_aligned %>%
              dplyr::filter(type == "Sample", model == m, group == "WT", sex %in% c("F", "M"))
            st_wt <- stats_5sets_by_model[[m]]$WT_FvsM

            if (!is.null(st_wt)) {
              out_png_wt <- file.path(
                mp$plots$heatmap_SIG_sex_FvsM,
                paste0(
                  "HEATMAP_SIG_ACTIVE_WT_FvsM_", met, "_lt_", alpha_sig,
                  "_model_", m, "_scale_", heatmap_scale_method, ".png"
                )
              )

              plot_sig_heatmap_from_stats(
                mat_log2_base,
                meta_wt,
                variants$ACTIVE$feature,
                st_wt,
                sig_metric = met,
                alpha_sig = alpha_sig,
                fc_cutoff_log2 = fc_cutoff_log2,
                require_fc_cutoff = sig_heatmap_require_fc_cutoff,
                sig_max = sig_heatmap_max_features,
                scale_method = heatmap_scale_method,
                order_samples_by_group = FALSE,
                out_png = out_png_wt,
                title_main = paste0(
                  "SIG (F/M within WT) | model=", m,
                  " | ", met, "<", alpha_sig,
                  if (sig_heatmap_require_fc_cutoff) paste0(" & |log2FC|>=", fc_cutoff_log2) else "",
                  " | log2FC=log2(F/M) | scale=", heatmap_scale_method
                )
              )
            }
          }
        }
      }

      step_info("Volcano plots enabled: ", make_volcano_plots)
      step_info("Volcano metrics: ", paste(run_metrics_expanded, collapse = ", "))
      step_info("Top heatmaps metrics: ", paste(heatmap_rank_metrics_expanded, collapse = ", "))
      step_info("Significant heatmaps metrics: ", paste(run_metrics_expanded, collapse = ", "))
      step_info("Stats Excel export enabled: ", save_stats_excel_per_model)
      step_info("Top heatmaps by model enabled: ", make_heatmap_by_model)
      step_info("Top heatmaps by model and sex enabled: ", make_heatmap_by_model_sex)
      step_info("Top heatmaps F vs M")
      step_info("Significant heatmaps by model enabled: ", make_sig_heatmap_by_model)
      step_info("Significant heatmaps by model and sex enabled: ", make_sig_heatmap_by_model_sex)
      step_info("Significant heatmaps F vs M enabled: ", make_sig_heatmap_FvsM_within_group)
      step_ok("Generating statistics, volcano plots, stats Excel, and heatmaps", t_step)
      # -------------------------------------------------------------------------
      # Finalize filter summary table and save
      # -------------------------------------------------------------------------

      if (file.exists(summary_csv)) {
        df_sum <- readr::read_csv(summary_csv, show_col_types = FALSE)
        log_written_object(log_path, summary_csv, df_sum, note = "Filter summary table")
      }

      remove_empty_directories(output_dir)

      message("\n", console_rule())
      message("[PIPELINE DONE] Total runtime: ", fmt_time_sec(pipeline_t0), " sec")
      message("[OUTPUT DIR] ", output_dir)
      message("[LOG FILE] ", log_path)
      message(console_rule())

      return(list(
        cd_raw = cd_raw,
        metadata = metadata,
        metadata_aligned = metadata_aligned,
        models_detected = models_detected,
        paths = paths,
        feature_tbl = feature_tbl,
        assay_bundle = assay_bundle,
        assay_num_raw = assay_num_raw,
        assay_num_weight = assay_num_weight,
        assay_num_pqn = assay_num_pqn,
        variants = variants,
        filter_summary = filter_summary,
        stats_5sets_by_model = stats_5sets_by_model,
        log_path = log_path
      ))
    })
  })
}
