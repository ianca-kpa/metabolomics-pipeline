# =============================================================================
# run_pipeline.R
# Entry point for the untargeted metabolomics pipeline
# =============================================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("
==================================================
")
cat("Untargeted metabolomics pipeline
")
cat("Runner started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
cat("==================================================
")

tryCatch({
  this_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "run_pipeline.R"
  project_dir <- normalizePath(dirname(this_file), winslash = "/", mustWork = FALSE)
  if (!dir.exists(project_dir)) project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  setwd(project_dir)

  cat("
Current working directory:
")
  cat(getwd(), "
")

  if (!file.exists("config/settings.R")) {
    stop("config/settings.R not found. Copy config/settings.example.R to config/settings.R and edit it before running the pipeline.")
  }

  cat("
[1/4] Loading settings.R ...
")
  source("config/settings.R", local = .GlobalEnv)
  cat("settings.R loaded successfully.
")

  cat("
[2/4] Loading pipeline modules ...
")

  module_files <- c(
    "R/00_packages.R",
    "R/01_validation.R",
    "R/02_comparisons.R",
    "R/03_helpers_io_log.R",
    "R/04_metadata.R",
    "R/05_features_assay.R",
    "R/06_normalization_filters.R",
    "R/07_duplicates.R",
    "R/08_exports.R",
    "R/09_pca.R",
    "R/10_stats_volcano.R",
    "R/11_heatmaps.R",
    "R/12_main_pipeline.R"
  )

  for (f in module_files) {
    cat("Loading:", f, "
")
    source(f, local = .GlobalEnv)
  }

  cat("
All modules loaded successfully.
")

  cat("
[3/4] Validating configuration ...
")
  validate_settings()
  cat("Configuration OK.
")

  cat("
[4/4] Running pipeline ...
")
  pipeline_result <- run_untargeted_pipeline()

  cat("
==================================================
")
  cat("Pipeline finished successfully.
")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
  cat("==================================================
")

  if (is.list(pipeline_result)) {
    cat("
Returned objects:
")
    cat(paste(names(pipeline_result), collapse = "
"), "
")
  }

  if (exists("output_dir")) {
    cat("
Output directory:
")
    cat(output_dir, "
")
  }

  if (exists("pipeline_result") && is.list(pipeline_result) && "log_path" %in% names(pipeline_result)) {
    cat("
Log file:
")
    cat(pipeline_result$log_path, "
")
  }

}, error = function(e) {
  cat("
==================================================
")
  cat("PIPELINE FAILED
")
  cat("Message:", conditionMessage(e), "
")
  cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
")
  cat("==================================================
")
  stop(e)
})
