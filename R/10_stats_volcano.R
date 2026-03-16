# =============================================================================
# 10_stats_volcano.R
# Statistics + volcano + Excel export
# =============================================================================

safe_var_equal_test <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 2 || length(y) < 2) {
    return(list(p = NA_real_, var_equal = NA))
  }

  vt <- tryCatch(var.test(x, y), error = function(e) NULL)
  if (is.null(vt)) {
    return(list(p = NA_real_, var_equal = NA))
  }

  p <- as.numeric(vt$p.value)
  list(p = p, var_equal = is.finite(p) && (p > 0.05))
}

FLOOR_P <- 1e-300

compute_ttest_stats_general <- function(mat_log2, mat_prelog, meta_sub, feat_info,
                                        compare_var = c("group", "sex"),
                                        num_level, den_level) {
  compare_var <- match.arg(compare_var)

  meta_sub <- meta_sub %>%
    filter(sample %in% rownames(mat_log2), sample %in% rownames(mat_prelog)) %>%
    mutate(.ord = match(sample, rownames(mat_log2))) %>%
    arrange(.ord) %>%
    select(-.ord)

  if (nrow(meta_sub) < 4) {
    return(NULL)
  }

  v <- meta_sub[[compare_var]]
  if (sum(v == den_level, na.rm = TRUE) < 2 || sum(v == num_level, na.rm = TRUE) < 2) {
    return(NULL)
  }

  s <- meta_sub$sample
  sub_log2 <- mat_log2[s, , drop = FALSE]
  sub_pre <- mat_prelog[s, , drop = FALSE]
  v2 <- meta_sub[[compare_var]]

  mean_den_pre <- colMeans(sub_pre[v2 == den_level, , drop = FALSE], na.rm = TRUE)
  mean_num_pre <- colMeans(sub_pre[v2 == num_level, , drop = FALSE], na.rm = TRUE)

  fc <- rep(NA_real_, length(mean_den_pre))
  ok <- is.finite(mean_den_pre) & is.finite(mean_num_pre) & mean_den_pre > 0 & mean_num_pre > 0
  fc[ok] <- mean_num_pre[ok] / mean_den_pre[ok]
  log2FC <- log2(fc)

  pvals <- rep(NA_real_, ncol(sub_log2))
  var_equal_p <- rep(NA_real_, ncol(sub_log2))
  var_equal <- rep(NA, ncol(sub_log2))

  for (j in seq_len(ncol(sub_log2))) {
    x <- sub_log2[v2 == den_level, j]
    y <- sub_log2[v2 == num_level, j]

    if (all(is.na(x)) || all(is.na(y))) next
    if (length(unique(na.omit(c(x, y)))) < 2) next

    vt <- safe_var_equal_test(x, y)
    var_equal_p[j] <- vt$p
    var_equal[j] <- vt$var_equal

    ve <- if (isTRUE(vt$var_equal)) TRUE else FALSE
    pvals[j] <- tryCatch(
      t.test(y, x, var.equal = ve)$p.value,
      error = function(e) NA_real_
    )
  }

  fdr <- p.adjust(pvals, method = "BH")

  tibble(
    featureID = colnames(sub_log2),
    FC_num_over_den = as.numeric(fc),
    log2FC_num_over_den = as.numeric(log2FC),
    p_value = pvals,
    FDR = fdr,
    var_equal_p = as.numeric(var_equal_p),
    var_equal = as.logical(var_equal)
  ) %>%
    left_join(
      feat_info %>%
        select(
          featureID,
          any_of(c(
            "display_name", "mz", "RT", "Name", "Name_canon",
            "Metabolika_pathways", "Formula"
          ))
        ),
      by = "featureID"
    )
}

plot_volcano_metric <- function(stats_df, title, out_path,
                                metric = c("FDR", "p_value"),
                                alpha = alpha_sig,
                                fc_cutoff_log2 = fc_cutoff_log2,
                                xlab = "log2FC") {
  metric <- match.arg(metric)

  if (is.null(stats_df) || nrow(stats_df) == 0) {
    stop("stats_df is NULL or empty.")
  }

  df <- stats_df %>%
    dplyr::mutate(
      metric_val = .data[[metric]],
      metric_plot = dplyr::if_else(
        is.na(metric_val),
        NA_real_,
        pmax(metric_val, FLOOR_P)
      ),
      minus_log10_metric = dplyr::if_else(
        is.na(metric_plot),
        NA_real_,
        -log10(metric_plot)
      ),
      label = dplyr::case_when(
        !is.na(display_name) & trimws(display_name) != "" ~ display_name,
        !is.na(Name) & trimws(Name) != "" ~ Name,
        TRUE ~ featureID
      ),

      # ---------------------------------------------------------
      # STRICT LOGIC:
      # Up/Down must respect BOTH significance and FC cutoff
      # ---------------------------------------------------------
      regulation = dplyr::case_when(
        !is.na(metric_val) &
          metric_val < alpha &
          !is.na(log2FC_num_over_den) &
          log2FC_num_over_den >= fc_cutoff_log2 ~ "Up",
        !is.na(metric_val) &
          metric_val < alpha &
          !is.na(log2FC_num_over_den) &
          log2FC_num_over_den <= -fc_cutoff_log2 ~ "Down",
        TRUE ~ "Normal"
      ),
      regulation = factor(regulation, levels = c("Down", "Normal", "Up"))
    )

  # -------------------------------------------------------------
  # Labels only for points that pass BOTH significance and FC cutoff
  # -------------------------------------------------------------
  df_labels <- df %>%
    dplyr::filter(
      !is.na(metric_val),
      !is.na(log2FC_num_over_den),
      metric_val < alpha,
      abs(log2FC_num_over_den) >= fc_cutoff_log2
    ) %>%
    dplyr::arrange(metric_val, dplyr::desc(abs(log2FC_num_over_den)))

  if (nrow(df_labels) > 15) {
    df_labels <- df_labels %>% dplyr::slice_head(n = 15)
  }

  # -------------------------------------------------------------
  # Axis limits
  # -------------------------------------------------------------
  x_vals <- df$log2FC_num_over_den[is.finite(df$log2FC_num_over_den)]
  y_vals <- df$minus_log10_metric[is.finite(df$minus_log10_metric)]

  max_abs_x <- if (length(x_vals) > 0) max(abs(x_vals), na.rm = TRUE) else 1
  max_abs_x <- max(max_abs_x, fc_cutoff_log2, 0.5)
  max_abs_x <- max_abs_x * 1.15

  sig_y <- -log10(alpha)
  max_y <- if (length(y_vals) > 0) max(y_vals, na.rm = TRUE) else 1
  max_y <- max(max_y, sig_y, 1) * 1.10

  # -------------------------------------------------------------
  # Dummy points only for legend
  # This forces Down / Normal / Up to always appear in legend
  # -------------------------------------------------------------
  legend_dummy <- tibble::tibble(
    log2FC_num_over_den = c(-999, -999, -999),
    minus_log10_metric = c(-999, -999, -999),
    regulation = factor(c("Down", "Normal", "Up"),
      levels = c("Down", "Normal", "Up")
    )
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log2FC_num_over_den,
      y = minus_log10_metric
    )
  ) +
    ggplot2::geom_hline(
      yintercept = sig_y,
      linetype = "dashed",
      linewidth = 0.5,
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff_log2, fc_cutoff_log2),
      linetype = "dashed",
      linewidth = 0.5,
      color = "black"
    ) +

    # real points
    ggplot2::geom_point(
      ggplot2::aes(fill = regulation),
      shape = 21,
      size = 3.0,
      color = "black",
      stroke = 0.25,
      alpha = 0.90,
      na.rm = TRUE
    ) +

    # dummy legend points
    ggplot2::geom_point(
      data = legend_dummy,
      ggplot2::aes(
        x = log2FC_num_over_den,
        y = minus_log10_metric,
        fill = regulation
      ),
      shape = 21,
      size = 3.0,
      color = "black",
      stroke = 0.25,
      alpha = 0.90,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = df_labels,
      ggplot2::aes(label = label),
      size = 3.8,
      box.padding = 0.35,
      point.padding = 0.20,
      segment.color = "grey40",
      segment.size = 0.30,
      max.overlaps = 100,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        Down = "#3B82F6",
        Normal = "grey75",
        Up = "#EF4444"
      ),
      drop = FALSE,
      name = "Regulation"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-max_abs_x, max_abs_x),
      ylim = c(0, max_y),
      expand = TRUE
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = paste0("-log10(", metric, ")")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(shape = 21, size = 3.5, color = "black", alpha = 1)
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      axis.title = ggplot2::element_text(face = "bold", size = 13),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = out_path,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

save_placeholder_volcano <- function(out_path, title, reason) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  p <- ggplot() +
    theme_void() +
    labs(title = title, subtitle = reason) +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = "black"),
      plot.subtitle = element_text(size = 12, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  ggsave(
    filename = out_path,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

run_all_stats_5sets_per_model <- function(mat_log2, mat_prelog, metadata_aligned, feat_info, paths,
                                          alpha_sig = 0.05, fc_cutoff_log2 = 1,
                                          run_metrics = c("FDR", "p_value"),
                                          make_volcano_plots = TRUE,
                                          comparison_configs = COMPARISON_CONFIGS) {
  models <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))
  out <- list()

  get_volcano_dir <- function(prefix, mp) {
    if (prefix == "ALL") {
      return(mp$plots$volcano_ALL)
    }
    if (prefix == "sex") {
      return(mp$plots$volcano_sex)
    }
    if (prefix == "sex_FvsM") {
      return(mp$plots$volcano_sex_FvsM)
    }
    stop("Unknown comparison prefix: ", prefix)
  }

  for (m in models) {
    mp <- get_model_paths(paths, m)
    meta_m <- metadata_aligned %>% filter(type == "Sample", model == m)
    out[[m]] <- list()

    for (comp_name in names(comparison_configs)) {
      cfg <- comparison_configs[[comp_name]]
      meta_sub <- cfg$meta_filter(meta_m)

      st <- compute_ttest_stats_general(
        mat_log2 = mat_log2,
        mat_prelog = mat_prelog,
        meta_sub = meta_sub,
        feat_info = feat_info,
        compare_var = cfg$stats_compare_var,
        num_level = cfg$stats_num,
        den_level = cfg$stats_den
      )

      out[[m]][[comp_name]] <- st

      if (!isTRUE(make_volcano_plots)) next

      for (met in run_metrics) {
        out_dir <- get_volcano_dir(cfg$prefix, mp)

        out_path <- file.path(
          out_dir,
          paste0("volcano_ACTIVE_model_", m, "_", comp_name, "_metric_", met, ".png")
        )

        title <- paste0(
          "Volcano | model=", m,
          " | ", cfg$label,
          " | metric=", met,
          " | log2FC=log2(", cfg$stats_num, "/", cfg$stats_den, ")"
        )

        if (is.null(st) || nrow(st) == 0) {
          save_placeholder_volcano(
            out_path = out_path,
            title = title,
            reason = "Not enough samples for this comparison (need >=2 per group and >=4 total), or no valid tests after filtering."
          )
          next
        }

        tryCatch(
          {
            plot_volcano_metric(
              stats_df = st,
              title = title,
              out_path = out_path,
              metric = met,
              alpha = alpha_sig,
              fc_cutoff = fc_cutoff_log2
            )
          },
          error = function(e) {
            message("  - Volcano failed: ", out_path)
            message("    Error: ", conditionMessage(e))

            save_placeholder_volcano(
              out_path = out_path,
              title = title,
              reason = paste0("Volcano failed with error: ", conditionMessage(e))
            )
          }
        )
      }
    }
  }

  out
}

export_stats_excel_by_model <- function(stats_5sets_by_model, paths, alpha_sig, fc_cutoff_log2,
                                        active_variant, log_path = NULL) {
  comparisons <- COMPARISON_NAMES

  col_idx <- function(df, colname) {
    idx <- match(colname, names(df))
    if (is.na(idx)) {
      return(NA_integer_)
    }
    idx
  }

  for (m in names(stats_5sets_by_model)) {
    mp <- get_model_paths(paths, m)
    out_dir <- mp$exports$stats
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    wb <- openxlsx::createWorkbook()

    readme <- tibble(
      field = c(
        "model", "ACTIVE_variant", "alpha_sig", "fc_cutoff_log2",
        "log2FC_definition", "stats_definition", "significance_logic"
      ),
      value = c(
        m, active_variant, alpha_sig, fc_cutoff_log2,
        "log2FC = log2(mean(num_prelog)/mean(den_prelog)); for FvsM: log2(F/M)",
        "t-test on log2 matrix; var.equal if var.test suggests equal variances else Welch; BH correction per comparison",
        "yellow row = (p_value<alpha OR FDR<alpha); green row = yellow row AND |log2FC|>=cutoff; red cell = p_value or FDR < alpha"
      )
    )

    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)

    tabs <- stats_5sets_by_model[[m]]

    for (nm in comparisons) {
      df <- tabs[[nm]]
      openxlsx::addWorksheet(wb, nm)

      if (is.null(df) || nrow(df) == 0) {
        tmp <- tibble(note = "Not enough samples or no valid tests")
        openxlsx::writeData(wb, nm, tmp)
        next
      }

      df_clean <- df %>%
        rename(
          FC = FC_num_over_den,
          log2FC = log2FC_num_over_den
        ) %>%
        mutate(
          row_sig_p = as.integer(!is.na(p_value) & p_value < alpha_sig),
          row_sig_fdr = as.integer(!is.na(FDR) & FDR < alpha_sig),
          row_sig_any = as.integer(row_sig_p == 1 | row_sig_fdr == 1),
          row_sig_and_fc = as.integer(
            row_sig_any == 1 &
              !is.na(log2FC) &
              abs(log2FC) >= fc_cutoff_log2
          )
        )

      openxlsx::writeData(wb, nm, df_clean)
      openxlsx::freezePane(wb, nm, firstRow = TRUE)
      openxlsx::addFilter(wb, nm, rows = 1, cols = 1:ncol(df_clean))

      style_red <- openxlsx::createStyle(fgFill = "#FF0000", fontColour = "#FFFFFF")
      style_yellow <- openxlsx::createStyle(fgFill = "#FFD966")
      style_green <- openxlsx::createStyle(fgFill = "#00B050", fontColour = "#FFFFFF")

      cols_full_row <- 1:ncol(df_clean)

      green_rows <- which(df_clean$row_sig_and_fc == 1) + 1
      yellow_rows <- which(df_clean$row_sig_any == 1 & df_clean$row_sig_and_fc == 0) + 1

      if (length(yellow_rows) > 0) {
        openxlsx::addStyle(
          wb, nm,
          style = style_yellow,
          rows = yellow_rows,
          cols = cols_full_row,
          gridExpand = TRUE,
          stack = FALSE
        )
      }

      if (length(green_rows) > 0) {
        openxlsx::addStyle(
          wb, nm,
          style = style_green,
          rows = green_rows,
          cols = cols_full_row,
          gridExpand = TRUE,
          stack = FALSE
        )
      }

      p_col <- col_idx(df_clean, "p_value")
      if (!is.na(p_col)) {
        p_rows <- which(!is.na(df_clean$p_value) & df_clean$p_value < alpha_sig) + 1
        if (length(p_rows) > 0) {
          openxlsx::addStyle(
            wb, nm,
            style = style_red,
            rows = p_rows,
            cols = p_col,
            gridExpand = TRUE,
            stack = FALSE
          )
        }
      }

      fdr_col <- col_idx(df_clean, "FDR")
      if (!is.na(fdr_col)) {
        fdr_rows <- which(!is.na(df_clean$FDR) & df_clean$FDR < alpha_sig) + 1
        if (length(fdr_rows) > 0) {
          openxlsx::addStyle(
            wb, nm,
            style = style_red,
            rows = fdr_rows,
            cols = fdr_col,
            gridExpand = TRUE,
            stack = FALSE
          )
        }
      }

      helper_cols <- c(
        col_idx(df_clean, "row_sig_p"),
        col_idx(df_clean, "row_sig_fdr"),
        col_idx(df_clean, "row_sig_any"),
        col_idx(df_clean, "row_sig_and_fc")
      )
      helper_cols <- helper_cols[!is.na(helper_cols)]

      visible_cols <- setdiff(seq_len(ncol(df_clean)), helper_cols)

      if (length(visible_cols) > 0) {
        openxlsx::setColWidths(wb, nm, cols = visible_cols, widths = "auto")
      }

      if (length(helper_cols) > 0) {
        openxlsx::setColWidths(wb, nm, cols = helper_cols, widths = 0)
      }
    }

    openxlsx::setColWidths(wb, "README", cols = 1:2, widths = "auto")

    out_xlsx <- file.path(out_dir, paste0("STATS_model_", m, ".xlsx"))
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

    if (!is.null(log_path)) {
      append_log_line(
        log_path,
        paste0("- ", basename(out_xlsx), " -> (Excel workbook) | path: ", out_xlsx)
      )
    }
  }

  message("  ✓ Stats Excel per model saved in: ", out_dir)
}
