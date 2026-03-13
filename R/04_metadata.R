# =============================================================================
# 04_metadata.R
# Metadata helpers
# =============================================================================

normalize_sex_value <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s[s %in% c("", "na", "n/a", "nan", "unknown", "unk", "-", "null")] <- NA_character_
  
  dplyr::case_when(
    is.na(s) ~ NA_character_,
    s %in% c("f","fem","female") ~ "F",
    s %in% c("m","mal","male") ~ "M",
    TRUE ~ NA_character_
  )
}

clean_metadata <- function(metadata) {
  needed <- c("sample","weight_mg","group","sex","model")
  miss <- setdiff(needed, names(metadata))
  
  if (length(miss) > 0) {
    stop("Metadata missing required columns: ", paste(miss, collapse = ", "))
  }
  
  metadata %>%
    mutate(
      sample = str_trim(as.character(sample)),
      group  = str_trim(as.character(group)),
      model  = str_trim(as.character(model)),
      sex    = normalize_sex_value(sex),
      weight_mg = suppressWarnings(as.numeric(weight_mg)),
      type = if (!("type" %in% names(metadata))) {
        if_else(str_detect(sample, "^QC"), "QC", "Sample")
      } else {
        str_trim(as.character(type))
      }
    )
}