# =============================================================================
# 09_pca.R
# PCA functions
# =============================================================================

scale_uv <- function(mat) {
  out <- scale(mat, center = TRUE, scale = TRUE)
  out <- as.matrix(out)
  out[is.na(out)] <- 0
  out
}

scale_pareto <- function(mat) {
  mu <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, stats::sd, na.rm = TRUE)
  denom <- sqrt(sdv)
  denom[is.na(denom) | denom == 0] <- 1

  out <- sweep(sweep(mat, 2, mu, "-"), 2, denom, "/")
  out <- as.matrix(out)
  out[is.na(out)] <- 0
  out
}

apply_pca_scaling <- function(mat, method = c("pareto", "autoscale", "none")) {
  method <- match.arg(method)

  if (method == "none") {
    out <- as.matrix(mat)
    out[is.na(out)] <- 0
    return(list(mat = out, label = "none"))
  }

  if (method == "pareto") {
    return(list(mat = scale_pareto(mat), label = "pareto"))
  }

  list(mat = scale_uv(mat), label = "autoscale")
}

# -----------------------------------------------------------------------------
# Remove invalid / constant columns before PCA
# -----------------------------------------------------------------------------
prepare_matrix_for_pca <- function(mat) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) {
    return(NULL)
  }

  keep_cols <- apply(mat, 2, function(v) {
    vv <- stats::na.omit(v)
    length(vv) > 0 && length(unique(vv)) > 1
  })

  mat2 <- mat[, keep_cols, drop = FALSE]

  if (ncol(mat2) < 2) {
    return(NULL)
  }

  mat2 <- as.matrix(mat2)
  mat2[is.na(mat2)] <- 0
  mat2
}

# -----------------------------------------------------------------------------
# Draw and save one PCA plot
# -----------------------------------------------------------------------------
plot_one_pca_subset <- function(mat_log2,
                                meta,
                                out_png,
                                title_main,
                                pca_scaling = "pareto",
                                color_var = "group",
                                shape_var = "sex",
                                log_path = NULL) {
  meta <- meta %>%
    dplyr::filter(sample %in% rownames(mat_log2)) %>%
    dplyr::distinct(sample, .keep_all = TRUE)

  if (nrow(meta) < 3) {
    message("  - PCA skipped: fewer than 3 samples for ", basename(out_png))
    return(FALSE)
  }

  mat_sub <- mat_log2[meta$sample, , drop = FALSE]
  mat_sub <- prepare_matrix_for_pca(mat_sub)

  if (is.null(mat_sub)) {
    message("  - PCA skipped: matrix too small or without variable features for ", basename(out_png))
    return(FALSE)
  }

  sc <- apply_pca_scaling(mat_sub, method = pca_scaling)

  pca <- tryCatch(
    stats::prcomp(sc$mat, center = FALSE, scale. = FALSE),
    error = function(e) NULL
  )

  if (is.null(pca)) {
    message("  - PCA failed for ", basename(out_png))
    return(FALSE)
  }

  if (is.null(pca$x) || ncol(pca$x) < 2) {
    message("  - PCA skipped: PC1/PC2 not available for ", basename(out_png))
    return(FALSE)
  }

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1 <- round(100 * var_exp[1], 1)
  pc2 <- round(100 * var_exp[2], 1)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta, by = "sample")

  if (!(color_var %in% names(scores))) {
    scores[[color_var]] <- "ALL"
  }
  if (!(shape_var %in% names(scores))) {
    scores[[shape_var]] <- "ALL"
  }

  scores[[color_var]] <- as.factor(scores[[color_var]])
  scores[[shape_var]] <- as.factor(scores[[shape_var]])

  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = .data[[color_var]],
      shape = .data[[shape_var]]
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = sample),
      size = 3,
      show.legend = FALSE,
      max.overlaps = 20
    ) +
    ggplot2::labs(
      title = title_main,
      x = paste0("PC1 (", pc1, "%)"),
      y = paste0("PC2 (", pc2, "%)"),
      color = color_var,
      shape = shape_var
    ) +
    ggplot2::theme_minimal()

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )

  if (!is.null(log_path)) {
    log_written_object(
      log_path,
      out_png,
      scores,
      note = paste0(
        "PCA plot | title=", title_main,
        " | scaling=", sc$label,
        " | color=", color_var,
        " | shape=", shape_var
      )
    )
  }

  message("  ✓ PCA saved: ", out_png)
  TRUE
}

# -----------------------------------------------------------------------------
# PCA driver
# - one PCA per model
# - one PCA per model restricted to WT
# - one PCA per model restricted to TG
# - one PCA per model restricted to F
# - one PCA per model restricted to M
# -----------------------------------------------------------------------------
plot_pca_per_model <- function(mat_log2,
                               metadata_aligned,
                               paths,
                               pca_scaling = "pareto",
                               log_path = NULL) {
  models <- sort(unique(metadata_aligned$model[metadata_aligned$type == "Sample"]))
  n_done <- 0

  for (m in models) {
    mp <- get_model_paths(paths, m)
    out_dir <- mp$plots$pca_ALL
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    meta_model <- metadata_aligned %>%
      dplyr::filter(type == "Sample", model == m)

    # -------------------------------------------------------------------------
    # 1) PCA with all samples from the model
    # color = group | shape = sex
    # -------------------------------------------------------------------------
    ok <- plot_one_pca_subset(
      mat_log2 = mat_log2,
      meta = meta_model,
      out_png = file.path(
        out_dir,
        paste0("PCA_ACTIVE_model_", m, "_scaling_", pca_scaling, ".png")
      ),
      title_main = paste0("PCA - model=", m, " (", pca_scaling, ")"),
      pca_scaling = pca_scaling,
      color_var = "group",
      shape_var = "sex",
      log_path = log_path
    )
    if (isTRUE(ok)) n_done <- n_done + 1

    # -------------------------------------------------------------------------
    # 2) PCA within WT and within TG
    # color = sex | shape = sex
    # -------------------------------------------------------------------------
    for (grp in c("WT", "TG")) {
      meta_grp <- meta_model %>%
        dplyr::filter(group == grp)

      ok <- plot_one_pca_subset(
        mat_log2 = mat_log2,
        meta = meta_grp,
        out_png = file.path(
          out_dir,
          paste0("PCA_ACTIVE_model_", m, "_group_", grp, "_scaling_", pca_scaling, ".png")
        ),
        title_main = paste0("PCA - model=", m, " | group=", grp, " (", pca_scaling, ")"),
        pca_scaling = pca_scaling,
        color_var = "sex",
        shape_var = "sex",
        log_path = log_path
      )
      if (isTRUE(ok)) n_done <- n_done + 1
    }

    # -------------------------------------------------------------------------
    # 3) PCA within F and within M
    # color = group | shape = group
    # -------------------------------------------------------------------------
    for (sx in c("F", "M")) {
      meta_sex <- meta_model %>%
        dplyr::filter(sex == sx)

      ok <- plot_one_pca_subset(
        mat_log2 = mat_log2,
        meta = meta_sex,
        out_png = file.path(
          out_dir,
          paste0("PCA_ACTIVE_model_", m, "_sex_", sx, "_scaling_", pca_scaling, ".png")
        ),
        title_main = paste0("PCA - model=", m, " | sex=", sx, " (", pca_scaling, ")"),
        pca_scaling = pca_scaling,
        color_var = "group",
        shape_var = "group",
        log_path = log_path
      )
      if (isTRUE(ok)) n_done <- n_done + 1
    }
  }

  message("  ✓ PCA plots created: ", n_done)
}