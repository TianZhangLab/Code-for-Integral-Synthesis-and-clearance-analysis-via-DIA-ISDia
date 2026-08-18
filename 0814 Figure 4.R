
library(readxl)
library(RColorBrewer)
library(colorspace)
library(plotly)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(minpack.lm)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(rlang)
library(purrr)

# -----------------------------------------------------------------------------
# 1. File paths
# -----------------------------------------------------------------------------
path_growth_xlsx      <- "D:/personal/UVA/SILAC/second rebuttal/Cell growth under ER stress.xlsx"
path_pooled_factors   <- "D:/personal/UVA/SILAC/second rebuttal/pooled_factors.csv"

path_df3_report       <- "D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20250402_121853_0327_A012_IT3.5_combine_Report.tsv"
path_df4_report       <- "D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20250420_162456_H730_SILAC_IT3.5_combined_Report.tsv"

path_df3_cv           <- "D:/personal/UVA/SILAC/rebuttal/df3_cv.csv"
path_df4_cv           <- "D:/personal/UVA/SILAC/rebuttal/df4_cv.csv"
path_df3_cv_corrected <- "D:/personal/UVA/SILAC/rebuttal/df3_cv_with_correct_light.csv"
path_df4_cv_corrected <- "D:/personal/UVA/SILAC/rebuttal/df4_cv_with_correct_light.csv"

# -----------------------------------------------------------------------------
# 2. Generic numeric helper functions
# -----------------------------------------------------------------------------
safe_mean3 <- function(x) if (sum(!is.na(x)) == 3) mean(x, na.rm = TRUE) else NA_real_
safe_se3   <- function(x) if (sum(!is.na(x)) == 3) sd(x, na.rm = TRUE) / sqrt(3) else NA_real_
safe_sum2  <- function(x, y) dplyr::if_else(!is.na(x) & !is.na(y), x + y, NA_real_)
safe_ratio <- function(num, den) dplyr::if_else(!is.na(num) & !is.na(den), num / den, NA_real_)
safe_cv    <- function(vec) {
  if (sum(!is.na(vec)) >= 3) sd(vec, na.rm = TRUE) / mean(vec, na.rm = TRUE) * 100 else NA_real_
}
safe_ttest_3 <- function(x, y) {
  if (sum(!is.na(x)) >= 3 && sum(!is.na(y)) >= 3) {
    tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

# -----------------------------------------------------------------------------
# 3. Per-protein CV / log2FC / t-test table (Light, Heavy, Total, H/L)
#    calculate_cv() is defined but NOT called by default -- df3_cv/df4_cv are
#    read from the cached csv files instead (see section 3b below)
# -----------------------------------------------------------------------------
calculate_cv <- function(df) {
  df %>%
    dplyr::mutate(
      con1_total = safe_sum2(con1_Channel1, con1_Channel2),
      con2_total = safe_sum2(con2_Channel1, con2_Channel2),
      con3_total = safe_sum2(con3_Channel1, con3_Channel2),
      tun1_total = safe_sum2(tun1_Channel1, tun1_Channel2),
      tun2_total = safe_sum2(tun2_Channel1, tun2_Channel2),
      tun3_total = safe_sum2(tun3_Channel1, tun3_Channel2),
      thg1_total = safe_sum2(thg1_Channel1, thg1_Channel2),
      thg2_total = safe_sum2(thg2_Channel1, thg2_Channel2),
      thg3_total = safe_sum2(thg3_Channel1, thg3_Channel2),
      
      con1_HL = safe_ratio(con1_Channel2, con1_Channel1),
      con2_HL = safe_ratio(con2_Channel2, con2_Channel1),
      con3_HL = safe_ratio(con3_Channel2, con3_Channel1),
      tun1_HL = safe_ratio(tun1_Channel2, tun1_Channel1),
      tun2_HL = safe_ratio(tun2_Channel2, tun2_Channel1),
      tun3_HL = safe_ratio(tun3_Channel2, tun3_Channel1),
      thg1_HL = safe_ratio(thg1_Channel2, thg1_Channel1),
      thg2_HL = safe_ratio(thg2_Channel2, thg2_Channel1),
      thg3_HL = safe_ratio(thg3_Channel2, thg3_Channel1)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      con_light_mean = safe_mean3(dplyr::c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_mean = safe_mean3(dplyr::c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_mean = safe_mean3(dplyr::c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      con_light_SE   = safe_se3(dplyr::c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_SE   = safe_se3(dplyr::c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_SE   = safe_se3(dplyr::c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      
      con_heavy_mean = safe_mean3(dplyr::c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_mean = safe_mean3(dplyr::c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_mean = safe_mean3(dplyr::c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      con_heavy_SE   = safe_se3(dplyr::c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_SE   = safe_se3(dplyr::c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_SE   = safe_se3(dplyr::c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      
      con_total_mean = safe_mean3(dplyr::c_across(c(con1_total, con2_total, con3_total))),
      tun_total_mean = safe_mean3(dplyr::c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_mean = safe_mean3(dplyr::c_across(c(thg1_total, thg2_total, thg3_total))),
      con_total_SE   = safe_se3(dplyr::c_across(c(con1_total, con2_total, con3_total))),
      tun_total_SE   = safe_se3(dplyr::c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_SE   = safe_se3(dplyr::c_across(c(thg1_total, thg2_total, thg3_total))),
      
      con_HL_mean = safe_mean3(dplyr::c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_mean = safe_mean3(dplyr::c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_mean = safe_mean3(dplyr::c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      con_HL_SE   = safe_se3(dplyr::c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_SE   = safe_se3(dplyr::c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_SE   = safe_se3(dplyr::c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      
      con_light_CV = safe_cv(dplyr::c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_CV = safe_cv(dplyr::c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_CV = safe_cv(dplyr::c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      
      con_heavy_CV = safe_cv(dplyr::c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_CV = safe_cv(dplyr::c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_CV = safe_cv(dplyr::c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      
      con_total_CV = safe_cv(dplyr::c_across(c(con1_total, con2_total, con3_total))),
      tun_total_CV = safe_cv(dplyr::c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_CV = safe_cv(dplyr::c_across(c(thg1_total, thg2_total, thg3_total))),
      
      con_HL_CV = safe_cv(dplyr::c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_CV = safe_cv(dplyr::c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_CV = safe_cv(dplyr::c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      
      light_tun_vs_con = log2(tun_light_mean / con_light_mean),
      light_thg_vs_con = log2(thg_light_mean / con_light_mean),
      heavy_tun_vs_con = log2(tun_heavy_mean / con_heavy_mean),
      heavy_thg_vs_con = log2(thg_heavy_mean / con_heavy_mean),
      total_tun_vs_con = log2(tun_total_mean / con_total_mean),
      total_thg_vs_con = log2(thg_total_mean / con_total_mean),
      HL_tun_vs_con    = log2(tun_HL_mean / con_HL_mean),
      HL_thg_vs_con    = log2(thg_HL_mean / con_HL_mean),
      
      light_tun_vs_con_p = safe_ttest_3(dplyr::c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1)),
                                        dplyr::c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      light_thg_vs_con_p = safe_ttest_3(dplyr::c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1)),
                                        dplyr::c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      heavy_tun_vs_con_p = safe_ttest_3(dplyr::c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2)),
                                        dplyr::c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      heavy_thg_vs_con_p = safe_ttest_3(dplyr::c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2)),
                                        dplyr::c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      total_tun_vs_con_p = safe_ttest_3(dplyr::c_across(c(tun1_total, tun2_total, tun3_total)),
                                        dplyr::c_across(c(con1_total, con2_total, con3_total))),
      total_thg_vs_con_p = safe_ttest_3(dplyr::c_across(c(thg1_total, thg2_total, thg3_total)),
                                        dplyr::c_across(c(con1_total, con2_total, con3_total))),
      HL_tun_vs_con_p = safe_ttest_3(dplyr::c_across(c(tun1_HL, tun2_HL, tun3_HL)),
                                     dplyr::c_across(c(con1_HL, con2_HL, con3_HL))),
      HL_thg_vs_con_p = safe_ttest_3(dplyr::c_across(c(thg1_HL, thg2_HL, thg3_HL)),
                                     dplyr::c_across(c(con1_HL, con2_HL, con3_HL)))
    ) %>%
    dplyr::ungroup()
}

# -----------------------------------------------------------------------------
# 3b. Read df3_cv / df4_cv directly from the cached csv files (no recompute)
#     If you need to recompute from the raw DIA reports, uncomment the block
#     below labeled "Recompute from raw reports"
# -----------------------------------------------------------------------------
df3_cv <- read.csv(path_df3_cv, check.names = FALSE)   # HeLa WT
df4_cv <- read.csv(path_df4_cv, check.names = FALSE)   # HeLa PERK KO

# --- Recompute from raw reports (uncomment if needed) ---
# read_and_rename_report <- function(path) {
#   read.delim(path, header = TRUE, sep = "\t") %>%
#     dplyr::rename_with(
#       ~ gsub("X\\.[0-9]+\\.\\.[0-9]+_(con|tun|thg)([0-9]).*MS2Channel([0-9])",
#              "\\1\\2_Channel\\3", .x),
#       dplyr::starts_with("X.")
#     )
# }
# df3 <- read_and_rename_report(path_df3_report)
# df4 <- read_and_rename_report(path_df4_report)
# df3_cv <- calculate_cv(df3)
# df4_cv <- calculate_cv(df4)

# =============================================================================
# 4. Cell growth curve analysis -> Light-channel "growth dilution" correction
#    factors
# =============================================================================
growth_raw <- readxl::read_excel(path_growth_xlsx, sheet = 1)

growth_dat <- growth_raw %>%
  dplyr::transmute(
    elapsed_h = as.numeric(Elapsed),
    WT_DMSO = as.numeric(`A012...3`),
    WT_Tun  = as.numeric(`A012...7`),
    WT_Thg  = as.numeric(`A012...11`),
    PERK_KO_DMSO = as.numeric(`H730...4`),
    PERK_KO_Tun  = as.numeric(`H730...8`),
    PERK_KO_Thg  = as.numeric(`H730...12`)
  )

growth_long <- growth_dat %>%
  tidyr::pivot_longer(cols = -elapsed_h, names_to = "group", values_to = "confluence") %>%
  dplyr::mutate(
    cell_line = dplyr::case_when(
      stringr::str_detect(group, "^WT_") ~ "WT",
      stringr::str_detect(group, "^PERK_KO_") ~ "PERK_KO"
    ),
    treatment = dplyr::case_when(
      stringr::str_detect(group, "DMSO$") ~ "DMSO",
      stringr::str_detect(group, "Tun$")  ~ "Tun",
      stringr::str_detect(group, "Thg$")  ~ "Thg"
    )
  ) %>%
  dplyr::filter(!is.na(confluence)) %>%
  dplyr::arrange(cell_line, treatment, elapsed_h)

growth_norm <- growth_long %>%
  dplyr::group_by(cell_line, treatment) %>%
  dplyr::arrange(elapsed_h, .by_group = TRUE) %>%
  dplyr::mutate(
    baseline_time_h = dplyr::first(elapsed_h),
    baseline_confluence = dplyr::first(confluence),
    growth_factor = confluence / baseline_confluence,
    # Expected fraction of the original Light pool remaining after
    # growth-associated dilution only
    light_dilution_factor = 1 / growth_factor
  ) %>%
  dplyr::ungroup()

growth_norm %>% dplyr::distinct(cell_line, treatment, baseline_time_h) %>% print()

target_times <- c(6, 12, 24, 48)

growth_target <- growth_norm %>%
  dplyr::group_by(cell_line, treatment) %>%
  dplyr::group_modify(~ {
    dplyr::bind_rows(lapply(target_times, function(tt) {
      .x %>%
        dplyr::slice_min(order_by = abs(elapsed_h - tt), n = 1, with_ties = FALSE) %>%
        dplyr::mutate(target_time_h = tt)
    }))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::select(cell_line, treatment, target_time_h,
                matched_elapsed_h = elapsed_h,
                confluence, growth_factor, light_dilution_factor)

# Time-point-specific factors, kept for inspection; the pooled version below
# is what actually gets applied downstream
timepoint_factors <- growth_target %>%
  tidyr::pivot_wider(names_from = treatment, values_from = c(growth_factor, light_dilution_factor), names_sep = "_") %>%
  dplyr::mutate(
    Tun_growth_only_light_retention_factor = light_dilution_factor_Tun / light_dilution_factor_DMSO,
    Thg_growth_only_light_retention_factor = light_dilution_factor_Thg / light_dilution_factor_DMSO,
    Tun_growth_only_light_offset_log2 = log2(Tun_growth_only_light_retention_factor),
    Thg_growth_only_light_offset_log2 = log2(Thg_growth_only_light_retention_factor),
    Tun_light_correction_multiplier = 1 / Tun_growth_only_light_retention_factor,
    Thg_light_correction_multiplier = 1 / Thg_growth_only_light_retention_factor,
    Tun_light_correction_log2 = log2(Tun_light_correction_multiplier),
    Thg_light_correction_log2 = log2(Thg_light_correction_multiplier)
  )

# Average across the 4 time points (6/12/24/48 h) -> factor actually applied
pooled_dilution <- growth_target %>%
  dplyr::group_by(cell_line, treatment) %>%
  dplyr::summarise(pooled_light_dilution_factor = mean(light_dilution_factor, na.rm = TRUE), .groups = "drop")

pooled_factors <- pooled_dilution %>%
  tidyr::pivot_wider(names_from = treatment, values_from = pooled_light_dilution_factor,
                     names_prefix = "pooled_dilution_") %>%
  dplyr::mutate(
    Tun_growth_only_light_retention_factor = pooled_dilution_Tun / pooled_dilution_DMSO,
    Thg_growth_only_light_retention_factor = pooled_dilution_Thg / pooled_dilution_DMSO,
    Tun_growth_only_light_offset_log2 = log2(Tun_growth_only_light_retention_factor),
    Thg_growth_only_light_offset_log2 = log2(Thg_growth_only_light_retention_factor),
    Tun_light_correction_multiplier = pooled_dilution_DMSO / pooled_dilution_Tun,
    Thg_light_correction_multiplier = pooled_dilution_DMSO / pooled_dilution_Thg,
    Tun_light_correction_log2 = log2(Tun_light_correction_multiplier),
    Thg_light_correction_log2 = log2(Thg_light_correction_multiplier)
  )

write.csv(pooled_factors, path_pooled_factors, row.names = FALSE)

# =============================================================================
# 5. Apply the growth-based correction to the Light channel
#    (used later by the "corrected Pink" GO groups)
# =============================================================================
select_pooled_factor_row <- function(pooled_factors, dataset_col, dataset_value) {
  if (!dataset_col %in% names(pooled_factors)) {
    stop("Column '", dataset_col, "' was not found in pooled_factors.")
  }
  row <- pooled_factors %>% dplyr::filter(as.character(.data[[dataset_col]]) == dataset_value)
  if (nrow(row) != 1) {
    stop("Expected exactly one pooled_factors row for ", dataset_col, " = ", dataset_value,
         ", but found ", nrow(row), ".")
  }
  row
}

pooled_factors_df3 <- select_pooled_factor_row(pooled_factors, "cell_line", "WT")
pooled_factors_df4 <- select_pooled_factor_row(pooled_factors, "cell_line", "PERK_KO")

get_light_factors <- function(pooled_factor_row, treatment = c("tun", "thg")) {
  treatment <- match.arg(treatment)
  if (nrow(pooled_factor_row) != 1) stop("pooled_factor_row must contain exactly one row.")
  
  prefix <- ifelse(treatment == "tun", "Tun", "Thg")
  correction_col <- paste0(prefix, "_light_correction_log2")
  offset_col     <- paste0(prefix, "_growth_only_light_offset_log2")
  
  correction_log2 <- as.numeric(pooled_factor_row[[correction_col]][1])
  growth_offset   <- as.numeric(pooled_factor_row[[offset_col]][1])
  
  if (!is.finite(correction_log2)) correction_log2 <- -growth_offset
  if (!is.finite(growth_offset))   growth_offset   <- -correction_log2
  
  if (!is.finite(correction_log2) && !is.finite(growth_offset)) {
    stop("No finite correction factor was found for ", treatment, ".")
  }
  if (abs(growth_offset + correction_log2) > 1e-8) {
    warning(prefix, " growth offset and correction log2 are not exact opposites; ",
            "using correction_log2 to derive the multiplier.")
    growth_offset <- -correction_log2
  }
  
  tibble::tibble(
    Treatment = treatment,
    correction_log2 = correction_log2,
    correction_multiplier = 2^correction_log2,
    growth_offset_log2 = growth_offset
  )
}

# Returns just the growth offset (log2), independent of any corrected table
get_growth_offset <- function(pooled_factor_row, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  get_light_factors(pooled_factor_row, treatment)$growth_offset_log2[1]
}

add_corrected_light_results <- function(raw_df, cv_df, pooled_factor_row) {
  
  required_raw_cols <- c("PG.ProteinAccessions", "PG.Genes",
                         paste0("con", 1:3, "_Channel1"),
                         paste0("tun", 1:3, "_Channel1"),
                         paste0("thg", 1:3, "_Channel1"))
  missing_raw <- setdiff(required_raw_cols, names(raw_df))
  if (length(missing_raw) > 0) {
    stop("The following raw Light columns are missing from raw_df: ", paste(missing_raw, collapse = ", "),
         "\n(check whether df3_cv.csv / df4_cv.csv still retain the raw Channel1 intensity columns)")
  }
  
  if (nrow(pooled_factor_row) != 1) stop("pooled_factor_row must contain exactly one row.")
  if (nrow(raw_df) != nrow(cv_df)) stop("raw_df and cv_df have different numbers of rows.")
  if (!identical(as.character(raw_df$PG.ProteinAccessions), as.character(cv_df$PG.ProteinAccessions))) {
    stop("Protein-accession row order differs between raw_df and cv_df.")
  }
  
  tun_factors <- get_light_factors(pooled_factor_row, "tun")
  thg_factors <- get_light_factors(pooled_factor_row, "thg")
  tun_mult <- tun_factors$correction_multiplier[1]
  thg_mult <- thg_factors$correction_multiplier[1]
  
  corrected <- raw_df %>%
    dplyr::transmute(
      con1_correct_light = con1_Channel1,
      con2_correct_light = con2_Channel1,
      con3_correct_light = con3_Channel1,
      tun1_correct_light = tun1_Channel1 * tun_mult,
      tun2_correct_light = tun2_Channel1 * tun_mult,
      tun3_correct_light = tun3_Channel1 * tun_mult,
      thg1_correct_light = thg1_Channel1 * thg_mult,
      thg2_correct_light = thg2_Channel1 * thg_mult,
      thg3_correct_light = thg3_Channel1 * thg_mult
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      con_correct_light_mean = safe_mean3(dplyr::c_across(c(con1_correct_light, con2_correct_light, con3_correct_light))),
      tun_correct_light_mean = safe_mean3(dplyr::c_across(c(tun1_correct_light, tun2_correct_light, tun3_correct_light))),
      thg_correct_light_mean = safe_mean3(dplyr::c_across(c(thg1_correct_light, thg2_correct_light, thg3_correct_light))),
      
      con_correct_light_SE = safe_se3(dplyr::c_across(c(con1_correct_light, con2_correct_light, con3_correct_light))),
      tun_correct_light_SE = safe_se3(dplyr::c_across(c(tun1_correct_light, tun2_correct_light, tun3_correct_light))),
      thg_correct_light_SE = safe_se3(dplyr::c_across(c(thg1_correct_light, thg2_correct_light, thg3_correct_light))),
      
      con_correct_light_CV = safe_cv(dplyr::c_across(c(con1_correct_light, con2_correct_light, con3_correct_light))),
      tun_correct_light_CV = safe_cv(dplyr::c_across(c(tun1_correct_light, tun2_correct_light, tun3_correct_light))),
      thg_correct_light_CV = safe_cv(dplyr::c_across(c(thg1_correct_light, thg2_correct_light, thg3_correct_light))),
      
      correct_light_tun_vs_con = log2(tun_correct_light_mean / con_correct_light_mean),
      correct_light_thg_vs_con = log2(thg_correct_light_mean / con_correct_light_mean),
      
      correct_light_tun_vs_con_p = safe_ttest_3(
        dplyr::c_across(c(tun1_correct_light, tun2_correct_light, tun3_correct_light)),
        dplyr::c_across(c(con1_correct_light, con2_correct_light, con3_correct_light))
      ),
      correct_light_thg_vs_con_p = safe_ttest_3(
        dplyr::c_across(c(thg1_correct_light, thg2_correct_light, thg3_correct_light)),
        dplyr::c_across(c(con1_correct_light, con2_correct_light, con3_correct_light))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Tun_light_correction_log2_applied = tun_factors$correction_log2[1],
      Thg_light_correction_log2_applied = thg_factors$correction_log2[1],
      Tun_light_correction_multiplier_applied = tun_mult,
      Thg_light_correction_multiplier_applied = thg_mult,
      Tun_growth_only_light_offset_log2_applied = tun_factors$growth_offset_log2[1],
      Thg_growth_only_light_offset_log2_applied = thg_factors$growth_offset_log2[1]
    )
  
  cv_df_clean <- cv_df %>% dplyr::select(-dplyr::any_of(names(corrected)))
  dplyr::bind_cols(cv_df_clean, corrected)
}

df3_cv_correct_light <- add_corrected_light_results(df3_cv, df3_cv, pooled_factors_df3)
df4_cv_correct_light <- add_corrected_light_results(df4_cv, df4_cv, pooled_factors_df4)

write.csv(df3_cv_correct_light, path_df3_cv_corrected, row.names = FALSE)
write.csv(df4_cv_correct_light, path_df4_cv_corrected, row.names = FALSE)

# =============================================================================
# 6. ColorGroup classification -- always uses raw (uncorrected) Light
# =============================================================================
get_fc_cols <- function(treatment) {
  list(
    light_fc = paste0("light_", treatment, "_vs_con"),
    light_p  = paste0("light_", treatment, "_vs_con_p"),
    heavy_fc = paste0("heavy_", treatment, "_vs_con"),
    heavy_p  = paste0("heavy_", treatment, "_vs_con_p"),
    total_fc = paste0("total_", treatment, "_vs_con"),
    total_p  = paste0("total_", treatment, "_vs_con_p"),
    hl       = paste0("HL_", treatment, "_vs_con")
  )
}

select_rename_fc <- function(df, cols) {
  df %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes,
                  !!rlang::sym(cols$light_fc), !!rlang::sym(cols$light_p),
                  !!rlang::sym(cols$heavy_fc), !!rlang::sym(cols$heavy_p),
                  !!rlang::sym(cols$total_fc), !!rlang::sym(cols$total_p),
                  !!rlang::sym(cols$hl)) %>%
    dplyr::rename(
      Light_FC = !!rlang::sym(cols$light_fc), Light_p = !!rlang::sym(cols$light_p),
      Heavy_FC = !!rlang::sym(cols$heavy_fc), Heavy_p = !!rlang::sym(cols$heavy_p),
      Total_FC = !!rlang::sym(cols$total_fc), Total_p = !!rlang::sym(cols$total_p),
      HL_FC    = !!rlang::sym(cols$hl)
    )
}

# Used for counting: full 8 groups (A-H) + NS
classify_full <- function(df, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  cols <- get_fc_cols(treatment)
  select_rename_fc(df, cols) %>%
    dplyr::mutate(
      ColorGroup = dplyr::case_when(
        Light_p > 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "D",
        Light_p > 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "F",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "G",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "B",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "C",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "H",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "A",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "E",
        TRUE ~ "NS"
      )
    ) %>%
    dplyr::filter(!is.na(Light_FC), !is.na(Heavy_FC), !is.na(HL_FC), !is.na(Total_FC))
}

# Used for plotting/heatmap: 6 groups (A-F) + NS -- G/H are not shown, as
# requested (any protein matching the G/H conditions falls into NS here)
classify_plot <- function(df, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  cols <- get_fc_cols(treatment)
  select_rename_fc(df, cols) %>%
    dplyr::mutate(
      ColorGroup = dplyr::case_when(
        Light_p > 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "D",
        Light_p > 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "F",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "B",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "C",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "A",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "E",
        TRUE ~ "NS"
      )
    ) %>%
    dplyr::filter(!is.na(Light_FC), !is.na(Heavy_FC), !is.na(HL_FC), !is.na(Total_FC))
}

count_HL_total_groups <- function(df, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  classify_full(df, treatment) %>%
    dplyr::count(ColorGroup, name = "Count") %>%
    dplyr::mutate(ColorGroup = factor(ColorGroup, levels = c("A","B","C","D","E","F","G","H","NS"))) %>%
    dplyr::arrange(ColorGroup)
}

cell_lines <- list(WT = df3_cv, PERK_KO = df4_cv)
treatments <- c("tun", "thg")

# Counts for all 4 combinations (WT/PERK_KO x tun/thg), raw Light
group_counts_all <- tidyr::expand_grid(cell_line = names(cell_lines), treatment = treatments) %>%
  purrr::pmap(function(cell_line, treatment) {
    count_HL_total_groups(cell_lines[[cell_line]], treatment) %>%
      dplyr::mutate(cell_line = cell_line, treatment = treatment, .before = 1)
  }) %>%
  dplyr::bind_rows()

group_counts_all


# =============================================================================
# 7. Scatter plot: Total FC vs H/L FC (raw Light grouping)
# =============================================================================
color_map_6 <- c(A = "#FF1F5B", B = "#00CD6C", C = "#F28522",
                 D = "#AF58BA", E = "#FFC61E", F = "#009ADE", NS = "grey80")

df_name_map <- list("df3_cv" = "WT", "df4_cv" = expression("PERK"^"-/-"))

make_title <- function(df_name, treatment) {
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")
  if (df_name %in% names(df_name_map)) {
    substitute(group ~ x * "/DMSO", list(group = df_name_map[[df_name]], x = treat_label))
  } else {
    paste0(df_name, " ", treat_label, "/DMSO")
  }
}

plot_HL_vs_total_by_group <- function(df, df_name, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  df_plot <- classify_plot(df, treatment)
  
  cor_res <- cor.test(df_plot$Total_FC, df_plot$HL_FC, use = "complete.obs")
  r_val <- round(cor_res$estimate, 2)
  p_val_display <- ifelse(cor_res$p.value < 1e-16, "< 2.2e-16", formatC(cor_res$p.value, format = "e", digits = 2))
  
  ggplot(df_plot, aes(x = Total_FC, y = HL_FC, color = ColorGroup)) +
    geom_point(data = df_plot %>% dplyr::filter(ColorGroup == "NS"), color = "grey80", alpha = 0.5, size = 1.8) +
    geom_point(data = df_plot %>% dplyr::filter(ColorGroup != "NS"), aes(color = ColorGroup), alpha = 0.8, size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50", linewidth = 0.8) +
    annotate("text", x = -3.5, y = 2, label = paste0("R = ", r_val, "\nP = ", p_val_display),
             hjust = 0, size = 4, fontface = "bold") +
    scale_color_manual(values = color_map_6) +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(-4, 4)) +
    annotate("text", x = -2.8, y = 3.8, label = paste0("n = ", nrow(df_plot)), hjust = 1, vjust = 1, size = 4) +
    ggtitle(make_title(df_name, treatment)) +
    labs(x = expression("Total log"[2] * "(FC)"), y = expression("H/L ratio log"[2] * "(FC)"), color = NULL) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = "black"), axis.text = element_text(size = 12),
          axis.title = element_text(size = 14), legend.position = "top", legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
}

plots_HL_vs_total <- list(
  WT_tun      = plot_HL_vs_total_by_group(df3_cv, "df3_cv", "tun"),
  WT_thg      = plot_HL_vs_total_by_group(df3_cv, "df3_cv", "thg"),
  PERK_KO_tun = plot_HL_vs_total_by_group(df4_cv, "df4_cv", "tun"),
  PERK_KO_thg = plot_HL_vs_total_by_group(df4_cv, "df4_cv", "thg")
)

# =============================================================================
# 8. Heavy vs Light scatter plot
#    Coloring: raw Light (classify_plot)
#    Red dashed line: cell-division (growth-rate) offset threshold, drawn for
#    WT/PERK_KO x Tun/Thg
# =============================================================================
# Given a PG.Genes string (which can contain multiple ";"-separated gene
# symbols per protein group) and a vector of genes of interest, return the
# first matching gene symbol for each row (or NA if none match)
extract_matched_gene <- function(pg_genes, genes_of_interest) {
  purrr::map_chr(strsplit(pg_genes, ";"), function(g) {
    hit <- intersect(trimws(g), genes_of_interest)
    if (length(hit) > 0) hit[1] else NA_character_
  })
}

# Given a PG.ProteinAccessions string (";"-separated) and a NAMED vector
# mapping label -> accession (e.g. c(HYOU1 = "Q9Y4L1")), return the label for
# each row whose accession matches (or NA if none match). Use this when you
# want to pin a label to one specific protein group instead of matching by
# gene symbol (which can be ambiguous if several protein groups share a gene
# name, or if PG.Genes doesn't exactly match the symbol you expect).
extract_matched_gene_by_accession <- function(pg_accessions, accession_map) {
  purrr::map_chr(strsplit(pg_accessions, ";"), function(a) {
    a <- trimws(a)
    hit_idx <- which(accession_map %in% a)
    if (length(hit_idx) > 0) names(accession_map)[hit_idx[1]] else NA_character_
  })
}

plot_light_heavy <- function(df, df_name, treatment = c("thg", "tun"), pooled_factor_row,
                             highlight_genes = NULL, highlight_accessions = NULL) {
  treatment <- match.arg(treatment)
  df_plot <- classify_plot(df, treatment)
  
  cor_res <- cor.test(df_plot$Heavy_FC, df_plot$Light_FC, use = "complete.obs")
  r_val <- round(cor_res$estimate, 2)
  p_val_display <- ifelse(cor_res$p.value < 1e-16, "< 2.2e-16", formatC(cor_res$p.value, format = "e", digits = 2))
  
  growth_offset <- get_growth_offset(pooled_factor_row, treatment)
  adjusted_lines <- c(-0.58 + growth_offset, 0.58 + growth_offset)
  
  p <- ggplot(df_plot, aes(x = Heavy_FC, y = Light_FC, color = ColorGroup)) +
    geom_point(data = df_plot %>% dplyr::filter(ColorGroup == "NS"), color = "grey80", alpha = 0.5, size = 1.8) +
    geom_point(data = df_plot %>% dplyr::filter(ColorGroup != "NS"), aes(color = ColorGroup), alpha = 0.8, size = 1.8) +
    # Original +-0.58 threshold lines (grey dashed)
    geom_hline(yintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50", linewidth = 0.8) +
    # Growth-rate (cell division) corrected adjusted threshold lines (red)
    geom_hline(yintercept = adjusted_lines, linetype = "longdash", color = "red", linewidth = 0.9) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey50", linewidth = 0.8) +
    annotate("text", x = -3.5, y = 2, label = paste0("R = ", r_val, "\nP = ", p_val_display),
             hjust = 0, size = 4, fontface = "bold") +
    annotate("text", x = 3.8, y = -3.8, label = paste0("Cell division offset = ", round(growth_offset, 3)),
             hjust = 1, vjust = 0, size = 3.5, color = "red") +
    scale_color_manual(values = color_map_6) +
    scale_x_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4, 2), oob = scales::squish) +
    scale_y_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4, 2), oob = scales::squish) +
    annotate("text", x = -2.8, y = 3.8, label = paste0("n = ", nrow(df_plot)), hjust = 1, vjust = 1, size = 4) +
    ggtitle(make_title(df_name, treatment)) +
    labs(x = expression("Heavy log"[2] * "(FC)"), y = expression("Light log"[2] * "(FC)"), color = NULL) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = "black"), axis.text = element_text(size = 12),
          axis.title = element_text(size = 14), legend.position = "top", legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  
  # Highlight & label specific genes/proteins of interest, if requested
  if (!is.null(highlight_genes) || !is.null(highlight_accessions)) {
    highlight_by_gene <- NULL
    if (!is.null(highlight_genes)) {
      highlight_by_gene <- df_plot %>%
        dplyr::mutate(highlight_label = extract_matched_gene(PG.Genes, highlight_genes)) %>%
        dplyr::filter(!is.na(highlight_label))
    }
    
    highlight_by_accession <- NULL
    if (!is.null(highlight_accessions)) {
      highlight_by_accession <- df_plot %>%
        dplyr::mutate(highlight_label = extract_matched_gene_by_accession(PG.ProteinAccessions, highlight_accessions)) %>%
        dplyr::filter(!is.na(highlight_label))
    }
    
    highlight_df <- dplyr::bind_rows(highlight_by_gene, highlight_by_accession) %>%
      dplyr::distinct(PG.ProteinAccessions, .keep_all = TRUE)
    
    p <- p +
      geom_point(data = highlight_df, shape = 21, color = "black", fill = NA, size = 3.2, stroke = 1.1) +
      ggrepel::geom_text_repel(data = highlight_df, aes(label = highlight_label),
                               color = "black", size = 3.8, fontface = "bold",
                               box.padding = 0.5, point.padding = 0.3,
                               max.overlaps = Inf, segment.color = "black", min.segment.length = 0)
  }
  
  p
}

# Genes to highlight on the WT Thg (Th) Heavy-vs-Light plot
# case the gene-symbol match was ambiguous or picked the wrong protein group
highlight_genes_th      <- c("THBS1", "HSPA5", "CTH", "NIBAN1")
highlight_accessions_th <- c(HYOU1 = "Q9Y4L1")

# All 4 combinations: WT/PERK_KO x tun/thg
plots_light_heavy <- list(
  WT_tun      = plot_light_heavy(df3_cv, "df3_cv", "tun", pooled_factors_df3),
  WT_thg      = plot_light_heavy(df3_cv, "df3_cv", "thg", pooled_factors_df3,
                                 highlight_genes = highlight_genes_th,
                                 highlight_accessions = highlight_accessions_th),
  PERK_KO_tun = plot_light_heavy(df4_cv, "df4_cv", "tun", pooled_factors_df4),
  PERK_KO_thg = plot_light_heavy(df4_cv, "df4_cv", "thg", pooled_factors_df4)
)

# -----------------------------------------------------------------------------
# 8b. Combined 2x2 figures, one for Th (Thg) and one for Tm (Tun):
#     top row    = WT     Total-vs-H/L scatter | WT     Heavy-vs-Light scatter
#     bottom row = PERK_KO Total-vs-H/L scatter | PERK_KO Heavy-vs-Light scatter
# -----------------------------------------------------------------------------
combined_thg <- (plots_HL_vs_total$WT_thg      | plots_light_heavy$WT_thg) /
  (plots_HL_vs_total$PERK_KO_thg | plots_light_heavy$PERK_KO_thg) +
  patchwork::plot_annotation(title = "Th (Thapsigargin) / DMSO")

combined_tun <- (plots_HL_vs_total$WT_tun      | plots_light_heavy$WT_tun) /
  (plots_HL_vs_total$PERK_KO_tun | plots_light_heavy$PERK_KO_tun) +
  patchwork::plot_annotation(title = "Tm (Tunicamycin) / DMSO")

combined_thg
combined_tun


# =============================================================================
# 9. Count of orange (C group: LightUp + HeavyDown + TotalNS) proteins that
#    fall below the growth-corrected adjusted line (raw Light for grouping,
#    growth offset for the threshold, same convention as the plot above)
# =============================================================================
count_orange_below_adjusted_line <- function(df, df_name, treatment = c("thg", "tun"), pooled_factor_row) {
  treatment <- match.arg(treatment)
  growth_offset <- get_growth_offset(pooled_factor_row, treatment)
  adjusted_line <- 0.58 + growth_offset
  
  orange_df <- classify_plot(df, treatment) %>% dplyr::filter(ColorGroup == "C")
  orange_total <- nrow(orange_df)
  orange_below <- sum(orange_df$Light_FC < adjusted_line, na.rm = TRUE)
  
  tibble::tibble(
    Dataset = df_name, Treatment = treatment,
    Original_Light_cutoff = 0.58, Growth_offset = growth_offset,
    Adjusted_Light_line = adjusted_line,
    Orange_total = orange_total, Orange_below_adjusted_line = orange_below,
    Orange_at_or_above_adjusted_line = orange_total - orange_below,
    Percent_below_adjusted_line = dplyr::if_else(orange_total > 0, orange_below / orange_total * 100, NA_real_)
  )
}

orange_group_counts <- dplyr::bind_rows(
  count_orange_below_adjusted_line(df3_cv, "WT", "tun", pooled_factors_df3),
  count_orange_below_adjusted_line(df3_cv, "WT", "thg", pooled_factors_df3),
  count_orange_below_adjusted_line(df4_cv, "PERK_KO", "tun", pooled_factors_df4),
  count_orange_below_adjusted_line(df4_cv, "PERK_KO", "thg", pooled_factors_df4)
)
orange_group_counts

# =============================================================================
# 10. Group-median heatmap (raw Light)
# =============================================================================
plot_group_median_heatmap <- function(df, treatment = c("tun", "thg"), dataset_label) {
  treatment <- match.arg(treatment)
  df_plot <- classify_plot(df, treatment)
  
  med <- df_plot %>%
    dplyr::select(PG.ProteinAccessions, ColorGroup, Light_FC, Heavy_FC, Total_FC, HL_FC) %>%
    dplyr::mutate(ColorGroup = dplyr::if_else(is.na(ColorGroup), "NS", ColorGroup)) %>%
    dplyr::group_by(ColorGroup) %>%
    dplyr::summarise(Light_FC_med = median(Light_FC, na.rm = TRUE),
                     Heavy_FC_med = median(Heavy_FC, na.rm = TRUE),
                     Total_FC_med = median(Total_FC, na.rm = TRUE),
                     HL_FC_med = median(HL_FC, na.rm = TRUE),
                     n = dplyr::n(), .groups = "drop")
  
  group_order <- c("A", "B", "C", "D", "E", "F", "NS")
  med <- med %>%
    dplyr::mutate(ColorGroup = factor(ColorGroup, levels = group_order)) %>%
    dplyr::arrange(ColorGroup) %>%
    dplyr::filter(!is.na(ColorGroup))
  
  med_mat <- med %>%
    dplyr::select(ColorGroup, Light_FC_med, Heavy_FC_med, Total_FC_med, HL_FC_med) %>%
    tibble::column_to_rownames("ColorGroup") %>%
    as.matrix()
  
  cap <- 2.5
  med_mat_cap <- pmax(pmin(med_mat, cap), -cap)
  bk <- seq(-cap, cap, length.out = 101)
  hm_cols <- colorRampPalette(c("#4A90E2", "#FFFFFF", "#E84A5F"))(length(bk) - 1)
  
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")
  pheatmap(
    med_mat_cap, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
    color = hm_cols, breaks = bk, cellwidth = 22, cellheight = 22,
    main = paste0("Median log2FC per ColorGroup (", dataset_label, " ", treat_label, ")"),
    display_numbers = TRUE, number_format = "%.2f", number_color = "black"
  )
  invisible(med)
}

heatmap_medians <- list(
  WT_tun      = plot_group_median_heatmap(df3_cv, "tun", "WT"),
  WT_thg      = plot_group_median_heatmap(df3_cv, "thg", "WT"),
  PERK_KO_tun = plot_group_median_heatmap(df4_cv, "tun", "PERK_KO"),
  PERK_KO_thg = plot_group_median_heatmap(df4_cv, "thg", "PERK_KO")
)

# =============================================================================
# 11. GO enrichment analysis
#     Blue / Purple                 -> original uncorrected code (raw Light)
#     Pink (H/L > 0) / Pink (H/L < 0) -> run once with raw Light and once with
#     growth-corrected Light
# =============================================================================
run_go_analysis <- function(df, group_name) {
  prot_ids <- df$PG.ProteinAccessions
  prot_ids_clean <- sub(";.*", "", prot_ids)
  prot_ids_clean <- unique(prot_ids_clean[!is.na(prot_ids_clean) & prot_ids_clean != ""])
  
  if (length(prot_ids_clean) == 0) { message(group_name, ": no valid UniProt IDs."); return(NULL) }
  
  id_map <- clusterProfiler::bitr(prot_ids_clean, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (is.null(id_map) || nrow(id_map) == 0) { message(group_name, ": no IDs mapped to ENTREZID."); return(NULL) }
  
  go_result <- clusterProfiler::enrichGO(gene = unique(id_map$ENTREZID), OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  if (is.null(go_result) || nrow(go_result@result) == 0) { message(group_name, ": no GO terms."); return(NULL) }
  
  go_filtered <- go_result@result %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    dplyr::mutate(Group = group_name) %>%
    dplyr::arrange(p.adjust)
  
  if (nrow(go_filtered) == 0) { message(group_name, ": no GO terms with p.adjust < 0.05."); return(NULL) }
  go_filtered
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) != 2) return(NA_real_)
    num <- as.numeric(z[1]); den <- as.numeric(z[2])
    if (!is.finite(num) || !is.finite(den) || den == 0) return(NA_real_)
    num / den
  }, numeric(1))
}

add_fold_enrichment <- function(go_df) {
  if (is.null(go_df) || nrow(go_df) == 0) return(NULL)
  go_df %>% dplyr::mutate(GeneRatio_numeric = ratio_to_numeric(GeneRatio),
                          BgRatio_numeric = ratio_to_numeric(BgRatio),
                          FoldEnrichment = GeneRatio_numeric / BgRatio_numeric)
}

# Blue / Purple / Pink(raw): original uncorrected Light
get_go_groups_raw <- function(df, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  light_col   <- paste0("light_", treatment, "_vs_con")
  light_p_col <- paste0("light_", treatment, "_vs_con_p")
  heavy_col   <- paste0("heavy_", treatment, "_vs_con")
  heavy_p_col <- paste0("heavy_", treatment, "_vs_con_p")
  total_col   <- paste0("total_", treatment, "_vs_con")
  total_p_col <- paste0("total_", treatment, "_vs_con_p")
  hl_col      <- paste0("HL_", treatment, "_vs_con")
  
  blue <- df %>% dplyr::filter(.data[[light_p_col]] > 0.05,
                               .data[[heavy_col]] > 0.58, .data[[heavy_p_col]] < 0.05,
                               .data[[total_col]] > 0.58, .data[[total_p_col]] < 0.05)
  purple <- df %>% dplyr::filter(.data[[light_p_col]] > 0.05,
                                 .data[[heavy_col]] < -0.58, .data[[heavy_p_col]] < 0.05,
                                 .data[[total_col]] < -0.58, .data[[total_p_col]] < 0.05)
  pink_hl_pos <- df %>% dplyr::filter(.data[[light_col]] > 0.58, .data[[light_p_col]] < 0.05,
                                      .data[[heavy_col]] > 0.58, .data[[heavy_p_col]] < 0.05,
                                      .data[[total_col]] > 0.58, .data[[total_p_col]] < 0.05,
                                      .data[[hl_col]] > 0)
  pink_hl_neg <- df %>% dplyr::filter(.data[[light_col]] > 0.58, .data[[light_p_col]] < 0.05,
                                      .data[[heavy_col]] > 0.58, .data[[heavy_p_col]] < 0.05,
                                      .data[[total_col]] > 0.58, .data[[total_p_col]] < 0.05,
                                      .data[[hl_col]] < 0)
  
  list(Blue = blue, Purple = purple,
       `Pink raw (H/L > 0)` = pink_hl_pos, `Pink raw (H/L < 0)` = pink_hl_neg)
}

# Pink(corrected): growth-corrected Light
get_go_groups_corrected <- function(df, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  cl_col   <- paste0("correct_light_", treatment, "_vs_con")
  cl_p_col <- paste0("correct_light_", treatment, "_vs_con_p")
  heavy_col   <- paste0("heavy_", treatment, "_vs_con")
  heavy_p_col <- paste0("heavy_", treatment, "_vs_con_p")
  total_col   <- paste0("total_", treatment, "_vs_con")
  total_p_col <- paste0("total_", treatment, "_vs_con_p")
  hl_col      <- paste0("HL_", treatment, "_vs_con")
  
  pink_hl_pos <- df %>% dplyr::filter(.data[[cl_col]] > 0.58, .data[[cl_p_col]] < 0.05,
                                      .data[[heavy_col]] > 0.58, .data[[heavy_p_col]] < 0.05,
                                      .data[[total_col]] > 0.58, .data[[total_p_col]] < 0.05,
                                      .data[[hl_col]] > 0)
  pink_hl_neg <- df %>% dplyr::filter(.data[[cl_col]] > 0.58, .data[[cl_p_col]] < 0.05,
                                      .data[[heavy_col]] > 0.58, .data[[heavy_p_col]] < 0.05,
                                      .data[[total_col]] > 0.58, .data[[total_p_col]] < 0.05,
                                      .data[[hl_col]] < 0)
  
  list(`Pink corrected (H/L > 0)` = pink_hl_pos, `Pink corrected (H/L < 0)` = pink_hl_neg)
}

# Run for WT only, both treatments: (Blue/Purple/Pink-raw) x tun/thg
# (PERK_KO is skipped per request -- only WT Tun/Thg plots are needed)
go_cell_lines_to_run <- c("WT")

go_results_raw <- list()
for (cl in go_cell_lines_to_run) {
  for (tr in treatments) {
    groups <- get_go_groups_raw(cell_lines[[cl]], tr)
    for (g in names(groups)) {
      label <- paste(cl, tr, g, sep = " | ")
      go_results_raw[[label]] <- run_go_analysis(groups[[g]], label)
    }
  }
}

# Run for WT only, both treatments: Pink-corrected x tun/thg
cell_lines_corrected <- list(WT = df3_cv_correct_light)  # PERK_KO not needed; only WT is run below
go_results_corrected <- list()
for (cl in go_cell_lines_to_run) {
  for (tr in treatments) {
    groups <- get_go_groups_corrected(cell_lines_corrected[[cl]], tr)
    for (g in names(groups)) {
      label <- paste(cl, tr, g, sep = " | ")
      go_results_corrected[[label]] <- run_go_analysis(groups[[g]], label)
    }
  }
}

# Base colors reused from the scatter-plot ColorGroup palette (color_map_6):
#   Blue GO group   <-> ColorGroup F pattern (Light NS, Heavy up, Total up)   -> "#009ADE"
#   Purple GO group <-> ColorGroup D pattern (Light NS, Heavy down, Total down) -> "#AF58BA"
#   Pink GO group   <-> ColorGroup A pattern (Light up, Heavy up, Total up)   -> "#FF1F5B"
go_group_base_colors <- c(Blue = "#009ADE", Purple = "#AF58BA", Pink = "#FF1F5B")

# Pick the right base color from a go_results_* list name such as
# "WT | thg | Blue" or "WT | thg | Pink corrected (H/L > 0)"
get_go_group_color <- function(label) {
  if (grepl("Blue", label))   return(go_group_base_colors[["Blue"]])
  if (grepl("Purple", label)) return(go_group_base_colors[["Purple"]])
  if (grepl("Pink", label))   return(go_group_base_colors[["Pink"]])
  go_group_base_colors[["Pink"]]  # fallback
}

plot_go_bubble <- function(go_df, top_n = 10, title_text = "GO Enrichment Bubble Plot",
                           base_color = "#FF1F5B") {
  if (is.null(go_df) || nrow(go_df) == 0) { message("No GO terms to plot."); return(NULL) }
  
  go_df <- add_fold_enrichment(go_df)
  go_top <- go_df %>%
    dplyr::filter(!is.na(p.adjust), is.finite(FoldEnrichment)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(log_p = -log10(pmax(p.adjust, .Machine$double.xmin)),
                  Description = factor(Description, levels = rev(Description)))
  
  if (nrow(go_top) == 0) { message("No GO terms remain after filtering."); return(NULL) }
  
  # Gradient from a light tint of base_color up to the full base_color itself
  low_color <- colorspace::lighten(base_color, amount = 0.75)
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +
    scale_fill_gradient(low = low_color, high = base_color, name = "-log10(p.adjust)") +
    scale_size_continuous(name = "Gene Count") +
    labs(title = title_text, x = "Fold Enrichment", y = "GO Term") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      
      plot.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}

# -----------------------------------------------------------------------------
# 11c. Actually generate & display the GO bubble plots
#      (go_results_raw / go_results_corrected only hold the enrichGO tables;
#      plot_go_bubble() has to be called explicitly to draw anything)
# -----------------------------------------------------------------------------
go_bubble_plots <- list()

for (label in names(go_results_raw)) {
  p <- plot_go_bubble(go_results_raw[[label]], top_n = 10, title_text = label,
                      base_color = get_go_group_color(label))
  if (!is.null(p)) go_bubble_plots[[label]] <- p
}

for (label in names(go_results_corrected)) {
  p <- plot_go_bubble(go_results_corrected[[label]], top_n = 10, title_text = label,
                      base_color = get_go_group_color(label))
  if (!is.null(p)) go_bubble_plots[[label]] <- p
}

# go_bubble_plots now holds one ggplot per group that had >=1 significant GO
# term (p.adjust < 0.05). Print them all (each print() opens a new plotting
# device page / RStudio Plots pane entry):
for (label in names(go_bubble_plots)) {
  print(go_bubble_plots[[label]])
}

# To view just one specific combination instead of looping through all of
# them, index directly, e.g.:
# go_bubble_plots[["WT | thg | Blue"]]
# go_bubble_plots[["WT | thg | Pink corrected (H/L > 0)"]]

# -----------------------------------------------------------------------------
# 11d. Summary table: every significant GO term (p.adjust < 0.05) for WT,
#      with the actual enriched gene symbols (geneID column, "/" separated)
# -----------------------------------------------------------------------------
path_go_summary_wt <- "D:/personal/UVA/SILAC/second rebuttal/GO_summary_WT.csv"
 
go_summary_wt <- dplyr::bind_rows(c(go_results_raw, go_results_corrected)) %>%
  add_fold_enrichment() %>%
  dplyr::select(Group, ID, Description, GeneRatio, BgRatio, FoldEnrichment,
                pvalue, p.adjust, qvalue, Count, geneID) %>%
  dplyr::arrange(Group, p.adjust)
 
write.csv(go_summary_wt, path_go_summary_wt, row.names = FALSE)
View(go_summary_wt)
 
# Split into one gene per row instead of a single "/"-separated string, in
# case you want to filter/count by individual gene:
go_summary_wt_long <- go_summary_wt %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(Gene = geneID)

# =============================================================================
# 12. light_up / heavy_down / total_no_change plot (kept as originally
#     provided, raw Light only)
# =============================================================================
plot_lightup_heavydown_totalns <- function(df3_cv, df4_cv, channel = c("light", "heavy", "total"), treatment = c("tun", "thg")) {
  channel <- match.arg(channel)
  treatment <- match.arg(treatment)
  
  ch_treat_vs_con     <- paste0(channel, "_", treatment, "_vs_con")
  ch_treat_vs_con_p   <- paste0(channel, "_", treatment, "_vs_con_p")
  light_col           <- paste0("light_", treatment, "_vs_con")
  light_p_col         <- paste0("light_", treatment, "_vs_con_p")
  heavy_col           <- paste0("heavy_", treatment, "_vs_con")
  heavy_p_col         <- paste0("heavy_", treatment, "_vs_con_p")
  total_col           <- paste0("total_", treatment, "_vs_con")
  total_p_col         <- paste0("total_", treatment, "_vs_con_p")
  
  selected_proteins <- df3_cv %>%
    dplyr::filter(
      .data[[light_col]] > 0.58, .data[[light_p_col]] < 0.05,
      .data[[heavy_col]] < -0.58, .data[[heavy_p_col]] < 0.05,
      .data[[total_p_col]] > 0.05
    ) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes)
  
  # df3 subset
  df3_selected <- df3_cv %>%
    dplyr::filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes,
                  !!light_col, !!heavy_col, !!total_col,
                  !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # df4 subset
  df4_selected <- df4_cv %>%
    dplyr::filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes,
                  !!light_col, !!heavy_col, !!total_col,
                  !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # rename for merging
  df3_selected <- df3_selected %>%
    dplyr::rename(
      light_df3 = !!light_col, heavy_df3 = !!heavy_col, total_df3 = !!total_col,
      light_p_df3 = !!light_p_col, heavy_p_df3 = !!heavy_p_col, total_p_df3 = !!total_p_col
    )
  
  df4_selected <- df4_selected %>%
    dplyr::rename(
      light_df4 = !!light_col, heavy_df4 = !!heavy_col, total_df4 = !!total_col,
      light_p_df4 = !!light_p_col, heavy_p_df4 = !!heavy_p_col, total_p_df4 = !!total_p_col
    )
  
  # merge
  merged_selected <- dplyr::inner_join(df3_selected, df4_selected, by = c("PG.ProteinAccessions", "PG.Genes"))
  
  # Color group
  merged_selected <- merged_selected %>%
    dplyr::mutate(ColorGroup = dplyr::case_when(
      .data[[paste0(channel, "_p_df3")]] < 0.05 &
        (is.na(.data[[paste0(channel, "_p_df4")]]) | .data[[paste0(channel, "_p_df4")]] >= 0.05)~ "sig_df3",
      .data[[paste0(channel, "_p_df3")]] >= 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "sig_df4",
      .data[[paste0(channel, "_p_df3")]] < 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "both",
      TRUE ~ "ns"
    ))
  
  # Color map
  color_map <- c("sig_df3" = "#009ade", "sig_df4" = "#ff1f5b", "both" = "#ffc61e", "ns" = "grey80")
  print(any(merged_selected$ColorGroup == "ns"))
  # plotting columns
  fc_df3_col <- paste0(channel, "_df3")
  fc_df4_col <- paste0(channel, "_df4")
  
  p <- ggplot(merged_selected, aes(x = .data[[fc_df3_col]], y = .data[[fc_df4_col]], color = ColorGroup)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = color_map) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    labs(
      title = paste0("Scatter of ", str_to_title(channel), " FC: df3 vs df4 (", treatment, ") LightUp_HeavyDown_TotalNS"),
      x = paste0("df3 ", channel, "_", treatment, "_vs_con"),
      y = paste0("df4 ", channel, "_", treatment, "_vs_con"),
      color = "Significance"
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}


plot_lightup_heavydown_totalns <- function(df3_cv, df4_cv, channel = c("light", "heavy", "total"), treatment = c("tun", "thg")) {
  channel <- match.arg(channel)
  treatment <- match.arg(treatment)
  
  ch_treat_vs_con     <- paste0(channel, "_", treatment, "_vs_con")
  ch_treat_vs_con_p   <- paste0(channel, "_", treatment, "_vs_con_p")
  light_col           <- paste0("light_", treatment, "_vs_con")
  light_p_col         <- paste0("light_", treatment, "_vs_con_p")
  heavy_col           <- paste0("heavy_", treatment, "_vs_con")
  heavy_p_col         <- paste0("heavy_", treatment, "_vs_con_p")
  total_col           <- paste0("total_", treatment, "_vs_con")
  total_p_col         <- paste0("total_", treatment, "_vs_con_p")
  
  selected_proteins <- df3_cv %>%
    dplyr::filter(
      .data[[light_col]] > 0.58, .data[[light_p_col]] < 0.05,
      .data[[heavy_col]] < -0.58, .data[[heavy_p_col]] < 0.05,
      .data[[total_p_col]] > 0.05
    ) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes)
  
  # df3 subset
  df3_selected <- df3_cv %>%
    dplyr::filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes,
                  !!light_col, !!heavy_col, !!total_col,
                  !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # df4 subset
  df4_selected <- df4_cv %>%
    dplyr::filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes,
                  !!light_col, !!heavy_col, !!total_col,
                  !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # rename for merging
  df3_selected <- df3_selected %>%
    dplyr::rename(
      light_df3 = !!light_col, heavy_df3 = !!heavy_col, total_df3 = !!total_col,
      light_p_df3 = !!light_p_col, heavy_p_df3 = !!heavy_p_col, total_p_df3 = !!total_p_col
    )
  
  df4_selected <- df4_selected %>%
    dplyr::rename(
      light_df4 = !!light_col, heavy_df4 = !!heavy_col, total_df4 = !!total_col,
      light_p_df4 = !!light_p_col, heavy_p_df4 = !!heavy_p_col, total_p_df4 = !!total_p_col
    )
  
  # merge
  merged_selected <- dplyr::inner_join(df3_selected, df4_selected, by = c("PG.ProteinAccessions", "PG.Genes"))
  
  # Color group
  merged_selected <- merged_selected %>%
    dplyr::mutate(ColorGroup = dplyr::case_when(
      .data[[paste0(channel, "_p_df3")]] < 0.05 &
        (is.na(.data[[paste0(channel, "_p_df4")]]) | .data[[paste0(channel, "_p_df4")]] >= 0.05)~ "sig_df3",
      .data[[paste0(channel, "_p_df3")]] >= 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "sig_df4",
      .data[[paste0(channel, "_p_df3")]] < 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "both",
      TRUE ~ "ns"
    ))
  
  # Color map + display labels (PERK -/- rendered with a superscript via plotmath)
  color_map <- c("sig_df3" = "#009ade", "sig_df4" = "#ff1f5b", "both" = "#ffc61e", "ns" = "grey80")
  color_breaks <- c("sig_df3", "sig_df4", "both", "ns")
  color_labels_txt <- c("WT", "PERK^{\"-/-\"}", "Both", "NS")
  print(any(merged_selected$ColorGroup == "ns"))
  # plotting columns
  fc_df3_col <- paste0(channel, "_df3")
  fc_df4_col <- paste0(channel, "_df4")
  
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")
  panel_title <- switch(channel, total = "H+L", light = "Light", heavy = "Heavy")
  x_lab <- substitute(WT ~ log[2] * "(" * tl * "/DMSO)", list(tl = treat_label))
  y_lab <- substitute(PERK^{"-/-"} ~ log[2] * "(" * tl * "/DMSO)", list(tl = treat_label))
  
  p <- ggplot(merged_selected, aes(x = .data[[fc_df3_col]], y = .data[[fc_df4_col]], color = ColorGroup)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = color_map, breaks = color_breaks, labels = parse(text = color_labels_txt)) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    labs(
      title = panel_title,
      x = x_lab,
      y = y_lab,
      color = NULL
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "top",
      legend.text = element_text(size = 8)
    )
  
  return(p)
}
p1<-plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "light", treatment = "tun")
p2<-plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "heavy", treatment = "tun")
p3<-plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "total", treatment = "tun")
combined_plot <- p3 + p1+ p2 +plot_layout(ncol = 3)
combined_plot

